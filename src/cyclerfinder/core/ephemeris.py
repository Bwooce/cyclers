"""Heliocentric inertial planet-state provider.

The :class:`Ephemeris` class supplies ``(r_km, v_km_s)`` heliocentric inertial
state vectors for the planets defined in :mod:`cyclerfinder.core.constants`.
Two backends:

- ``"circular"`` (M1) — every planet on its mean-motion circle in the
  ecliptic plane. At ``t_sec = 0`` every planet sits at ``theta = 0``.
- ``"astropy"`` (M6-slice, 2026-06-01) — real heliocentric position/velocity
  from astropy's bundled JPL DE440 ephemeris, computed at the absolute epoch
  ``2000-01-01T12:00:00 TDB + t_sec``. ``+x`` toward the J2000 vernal equinox,
  ``+z`` along the ecliptic north pole.

The two backends share the ``Ephemeris.state(body, t_sec)`` signature. Callers
that don't care about the epoch convention can swap one for the other freely.
Callers that DO care (M6b phase-matching, real launch-window finders) must use
``"astropy"`` because only it returns true heliocentric positions on real
dates; ``"circular"`` is a per-cycler idealisation with no absolute calendar
meaning.

References
----------
Spec §4, §6, §12.1 (phase-match). Plans:
``docs/phases/m1-core-mechanics/plan.md`` §3.1 (circular backend),
``docs/phases/m4-enumeration-scoring/plan.md`` (M4 consumer notes),
and the launch-windows slice (2026-06-01) that brought the astropy backend
forward from full M6.
"""

from __future__ import annotations

from datetime import UTC, datetime
from math import cos, pi, sin
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, PLANETS, SECONDS_PER_DAY

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

# astropy epoch convention for the astropy backend: t_sec is seconds since
# J2000 (2000-01-01T12:00:00 TDB). Match the most common heliocentric
# convention in the cycler literature.
_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

# Body-code → astropy body name map. The circular backend uses the same V/E/M
# codes; the astropy backend translates to astropy's own naming.
_ASTROPY_BODY_NAMES: dict[str, str] = {
    "V": "venus",
    "E": "earth",
    "M": "mars",
}


class _Backend(Protocol):
    """Strategy interface — both backends conform to a single state() method."""

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:  # pragma: no cover
        ...


class _CircularBackend:
    """Circular, coplanar, prograde planet motion in the J2000 ecliptic plane."""

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        # KeyError is the right error if the caller passes an unknown body.
        planet = PLANETS[body]
        a_km = planet.sma_au * AU_KM
        # Convert tabulated mean motion (deg/day) to rad/s.
        n_rad_s = planet.mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        theta = n_rad_s * t_sec
        cos_t = cos(theta)
        sin_t = sin(theta)
        # Circular orbital speed = a * n.
        speed = a_km * n_rad_s
        r = np.array([a_km * cos_t, a_km * sin_t, 0.0], dtype=np.float64)
        v = np.array([-speed * sin_t, speed * cos_t, 0.0], dtype=np.float64)
        return r, v


class _AstropyBackend:
    """Real heliocentric states from astropy's bundled JPL DE440 ephemeris.

    Returns ICRS heliocentric ``(r_km, v_km_s)`` at the astropy ``Time``
    corresponding to ``J2000 + t_sec`` (TDB). ICRS is barycentric-aligned with
    J2000 ecliptic; for cycler work (V/E/M heliocentric positions on
    interplanetary scales) the J2000 ecliptic alignment is what callers want.
    """

    def __init__(self) -> None:
        # Import lazily so packages that omit the runtime dep can still import
        # this module's circular path.
        from astropy.coordinates import solar_system_ephemeris

        # DE440 is bundled with astropy as of 6.x+; opt in once.
        solar_system_ephemeris.set("de440")

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        from astropy.coordinates import get_body_barycentric_posvel
        from astropy.time import Time

        if body not in _ASTROPY_BODY_NAMES:
            raise KeyError(f"unknown body code {body!r}; astropy backend handles V/E/M")
        astropy_name = _ASTROPY_BODY_NAMES[body]
        # Time arithmetic via timedelta keeps timezone-aware datetimes consistent.
        epoch = _J2000_EPOCH.timestamp() + t_sec
        t = Time(epoch, format="unix", scale="tdb")
        # Heliocentric posvel: subtract Sun's barycentric posvel from the
        # body's. astropy returns CartesianRepresentation in km and km/s when
        # we cast through .xyz.to(km).value.
        body_pos, body_vel = get_body_barycentric_posvel(astropy_name, t)
        sun_pos, sun_vel = get_body_barycentric_posvel("sun", t)
        r_helio = (body_pos - sun_pos).xyz.to("km").value
        v_helio = (body_vel - sun_vel).xyz.to("km/s").value
        return (
            np.asarray(r_helio, dtype=np.float64),
            np.asarray(v_helio, dtype=np.float64),
        )


class Ephemeris:
    """Planet-state provider; selects the backend at construction time.

    Parameters
    ----------
    model:
        Backend name.

        - ``"circular"`` (default) — analytic states for the planets in
          :data:`cyclerfinder.core.constants.PLANETS`; perfect circles in the
          ecliptic, every planet at ``theta = 0`` at ``t_sec = 0``.
        - ``"astropy"`` — real heliocentric states from astropy's bundled
          JPL DE440 ephemeris. ``t_sec = 0`` corresponds to the J2000 epoch
          (2000-01-01T12:00:00 TDB).
    """

    def __init__(self, model: str = "circular") -> None:
        if model == "circular":
            self._backend: _Backend = _CircularBackend()
        elif model == "astropy":
            self._backend = _AstropyBackend()
        else:
            raise ValueError(f"unknown ephemeris model {model!r}; expected 'circular' or 'astropy'")
        self._model: str = model

    @property
    def model(self) -> str:
        """Name of the backend in use."""
        return self._model

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        """Return ``(r_km, v_km_s)`` heliocentric inertial state for ``body``.

        Parameters
        ----------
        body:
            One-letter body code from
            :data:`cyclerfinder.core.constants.PLANETS` (``"V"``, ``"E"``,
            ``"M"``). Unknown codes raise :class:`KeyError`.
        t_sec:
            Time since the M1 reference epoch, seconds. May be negative.

        Returns
        -------
        (r_km, v_km_s):
            Two ``numpy.ndarray`` of shape ``(3,)`` and ``dtype=float64``, in
            the heliocentric ecliptic J2000 frame described in the module
            docstring.
        """
        return self._backend.state(body, t_sec)
