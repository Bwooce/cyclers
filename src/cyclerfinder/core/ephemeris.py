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

from dataclasses import replace
from datetime import UTC, datetime
from math import cos, pi, radians, sin
from typing import Final, Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, PLANETS, SECONDS_PER_DAY, PlanetData

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

# astropy epoch convention for the astropy backend: t_sec is seconds since
# J2000 (2000-01-01T12:00:00 TDB). Match the most common heliocentric
# convention in the cycler literature.
_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

# Mean obliquity of the ecliptic at J2000 (IAU 2006 / IERS 2010). Used by
# the astropy backend to rotate ICRS (BCRS-aligned, equatorial) into the
# J2000 ecliptic frame the rest of cyclerfinder expects (z-axis along the
# ecliptic north pole).
_J2000_OBLIQUITY_RAD: Final[float] = radians(23.4392911)

# Body-code → astropy body name map. Derived from the central PLANETS registry
# so any (sourced) body added there is automatically resolvable here: astropy's
# solar-system body names are the lowercased planet names (``"venus"``,
# ``"earth"``, ``"mars"``, ``"jupiter"`` …). Bodies whose astropy slug differs
# from ``name.lower()`` (e.g. a barycenter or a moon) get an explicit override.
_ASTROPY_NAME_OVERRIDES: Final[dict[str, str]] = {}
_ASTROPY_BODY_NAMES: dict[str, str] = {
    code: _ASTROPY_NAME_OVERRIDES.get(code, data.name.lower()) for code, data in PLANETS.items()
}


# Sourced J2000 inclination / longitude-of-ascending-node, Standish & Williams,
# "Approximate Positions of the Planets", JPL Solar System Dynamics, Table 1
# (valid 1800-2050 AD) — the same source quoted in the PLANETS comments
# (constants.py). Earth defines the ecliptic, so its inc/lan stay 0.0. These
# live here (not in PLANETS) so the live coplanar table is never mutated.
_INCLINED_ELEMENTS_J2000: Final[dict[str, tuple[float, float]]] = {
    "V": (3.39467605, 76.67984255),
    "M": (1.84969142, 49.55953891),
}


def inclined_planets() -> dict[str, PlanetData]:
    """Return a NEW planet table with the sourced J2000 inc/Ω filled in.

    A copy of :data:`~cyclerfinder.core.constants.PLANETS` where each body's
    ``inc_deg``/``lan_deg`` are replaced with the Standish & Williams Table 1
    values (Venus i=3.39467605°, Ω=76.67984255°; Mars i=1.84969142°,
    Ω=49.55953891°). Bodies without a tabulated entry (Earth, the
    ecliptic-defining body) keep their live ``0.0`` defaults.

    The live ``PLANETS`` dict is **never mutated** — this returns a fresh dict of
    fresh :class:`PlanetData` records (via :func:`dataclasses.replace`), so the
    coplanar default stays byte-identical for every existing caller. Intended to
    be injected into :class:`_InclinedCircularBackend` (see
    :meth:`Ephemeris.inclined_circular`).
    """
    out: dict[str, PlanetData] = {}
    for code, data in PLANETS.items():
        if code in _INCLINED_ELEMENTS_J2000:
            inc_deg, lan_deg = _INCLINED_ELEMENTS_J2000[code]
            out[code] = replace(data, inc_deg=inc_deg, lan_deg=lan_deg)
        else:
            out[code] = replace(data)
    return out


class _Backend(Protocol):
    """Strategy interface — both backends conform to a single state() method."""

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:  # pragma: no cover
        ...


def _circular_inplane_state(planet: PlanetData, t_sec: float) -> tuple[Vec3, Vec3]:
    """In-plane circular state in the orbital plane (z == 0, prograde).

    The shared kernel for both the flat and inclined backends. The returned
    state lives in the orbital plane with the reference direction (where the
    body sits at ``t_sec == 0``, i.e. the ascending-node direction for an
    inclined orbit) along ``+x``.
    """
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


class _CircularBackend:
    """Circular, prograde planet motion in the J2000 ecliptic frame.

    Coplanar (z == 0) by default. For a planet whose ``inc_deg != 0.0`` this
    delegates to :class:`_InclinedCircularBackend`, which rotates the in-plane
    circular state into a tilted orbital plane. The ``inc_deg == 0.0`` test is
    an EXACT float comparison so the coplanar path stays byte-identical for any
    body that has not been given a (sourced) non-zero inclination.
    """

    def __init__(self) -> None:
        self._inclined = _InclinedCircularBackend()

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        # KeyError is the right error if the caller passes an unknown body.
        planet = PLANETS[body]
        if planet.inc_deg != 0.0:
            return self._inclined.state(body, t_sec)
        return _circular_inplane_state(planet, t_sec)


class _InclinedCircularBackend:
    """Circular planet motion in an orbital plane inclined to the ecliptic.

    The in-plane circular state (z == 0, body along ``+x`` at ``t_sec == 0``)
    is rotated into the J2000 ecliptic frame by ``R_z(-lan) @ R_x(-inc)``:
    first tilt the plane about the node line (the ``+x`` reference direction)
    by the inclination, then rotate about the ecliptic ``+z`` by the longitude
    of the ascending node. With this convention the body sits on the ascending
    node at ``t_sec == 0`` (z == 0 there), and the orbit normal points along
    ``n_hat = (-sin(lan) sin(inc), cos(lan) sin(inc), cos(inc))``.

    For ``inc_deg == 0.0`` the rotation is the identity, so the result is
    numerically identical to the flat in-plane state.
    """

    def __init__(self, planets: dict[str, PlanetData] | None = None) -> None:
        self._planets: dict[str, PlanetData] = PLANETS if planets is None else planets

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        planet = self._planets[body]
        r_plane, v_plane = _circular_inplane_state(planet, t_sec)
        rot = self._rotation(planet)
        r = rot @ r_plane
        v = rot @ v_plane
        return (
            np.asarray(r, dtype=np.float64),
            np.asarray(v, dtype=np.float64),
        )

    @staticmethod
    def _rotation(planet: PlanetData) -> NDArray[np.float64]:
        """Orbital-plane -> ecliptic rotation ``R_z(-lan) @ R_x(-inc)``."""
        inc = radians(planet.inc_deg)
        lan = radians(planet.lan_deg)
        ci, si = cos(inc), sin(inc)
        cl, sl = cos(lan), sin(lan)
        # R_x(-inc): tilt the orbital plane about the node line (+x) so the
        # orbit normal acquires a +y, +z component before the node rotation.
        rx = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, ci, si],
                [0.0, -si, ci],
            ],
            dtype=np.float64,
        )
        # R_z(+lan): rotate the node line to its ecliptic longitude.
        rz = np.array(
            [
                [cl, -sl, 0.0],
                [sl, cl, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        return rz @ rx


class _AstropyBackend:
    """Real heliocentric states from astropy's bundled JPL DE440 ephemeris.

    Returns heliocentric ``(r_km, v_km_s)`` in the **J2000 ecliptic
    frame** (z-axis along the ecliptic north pole) at the astropy
    ``Time`` corresponding to ``J2000 + t_sec`` (TDB).

    astropy's :func:`get_body_barycentric_posvel` returns BCRS
    coordinates, which are aligned with ICRS — the equatorial mean
    pole and equinox of J2000. The cyclerfinder ``circular`` backend
    and the rest of the package (frames, propagator, dynamic frame)
    work in the **ecliptic** frame, so we rotate the ICRS output by
    the J2000 obliquity (-23.4393 deg about +x) to produce ecliptic-
    aligned coordinates. Without this rotation Earth's heliocentric
    z-component oscillates by ~ ±58 million km over a year (the
    23.4393 deg equator/ecliptic tilt), which silently breaks any
    rotating-frame transform that assumes z is the ecliptic normal —
    exactly the M6a spec §10 risk.
    """

    def __init__(self) -> None:
        # Import lazily so packages that omit the runtime dep can still import
        # this module's circular path.
        from astropy.coordinates import solar_system_ephemeris

        # DE440 is bundled with astropy as of 6.x+; opt in once.
        solar_system_ephemeris.set("de440")
        # Pre-compute the ICRS->ecliptic rotation about +x by -obliquity.
        eps = _J2000_OBLIQUITY_RAD
        self._r_icrs_to_ecl: NDArray[np.float64] = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, cos(eps), sin(eps)],
                [0.0, -sin(eps), cos(eps)],
            ],
            dtype=np.float64,
        )

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        from astropy.coordinates import get_body_barycentric_posvel
        from astropy.time import Time

        if body not in _ASTROPY_BODY_NAMES:
            raise KeyError(
                f"unknown body code {body!r}; astropy backend handles {tuple(_ASTROPY_BODY_NAMES)}"
            )
        astropy_name = _ASTROPY_BODY_NAMES[body]
        # Time arithmetic via timedelta keeps timezone-aware datetimes consistent.
        epoch = _J2000_EPOCH.timestamp() + t_sec
        t = Time(epoch, format="unix", scale="tdb")
        # Heliocentric posvel: subtract Sun's barycentric posvel from the
        # body's. astropy returns CartesianRepresentation in km and km/s when
        # we cast through .xyz.to(km).value.
        body_pos, body_vel = get_body_barycentric_posvel(astropy_name, t)
        sun_pos, sun_vel = get_body_barycentric_posvel("sun", t)
        r_helio_icrs = (body_pos - sun_pos).xyz.to("km").value
        v_helio_icrs = (body_vel - sun_vel).xyz.to("km/s").value
        # Rotate ICRS (equatorial) to J2000 ecliptic so the rest of
        # cyclerfinder gets z-along-ecliptic-pole positions matching
        # the circular backend's convention.
        r_helio = self._r_icrs_to_ecl @ np.asarray(r_helio_icrs, dtype=np.float64)
        v_helio = self._r_icrs_to_ecl @ np.asarray(v_helio_icrs, dtype=np.float64)
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
