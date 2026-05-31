"""Heliocentric inertial planet-state provider.

The :class:`Ephemeris` class supplies ``(r_km, v_km_s)`` heliocentric inertial
state vectors for the three planets defined in
:mod:`cyclerfinder.core.constants`. M1 only implements the ``"circular"``
backend (every planet on its mean-motion circle in the ecliptic plane). The
``"astropy"`` backend is an M6 deliverable; constructing it raises
:class:`NotImplementedError`.

Reference frame
---------------
A right-handed heliocentric frame with ``+x`` toward the J2000 vernal equinox
and ``+z`` along the ecliptic north pole. The circular backend places every
planet in the ``z = 0`` plane, so the returned ``r[2]`` and ``v[2]`` components
are exactly zero. Epoch convention: at ``t_sec = 0`` every planet sits at
``theta = 0`` (i.e. ``+x`` direction). M1 does not need a real J2000 phase
calibration; later phases (M3 Aldrin reproduction, M6 JPL backend) choose
their own epoch.

References
----------
Spec §4, §6. Plan: ``docs/phases/m1-core-mechanics/plan.md`` §3.1.
"""

from __future__ import annotations

from math import cos, pi, sin
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, PLANETS, SECONDS_PER_DAY

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64


class _Backend(Protocol):
    """Strategy interface so the M6 astropy backend can drop in unchanged."""

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


class Ephemeris:
    """Planet-state provider; selects the backend at construction time.

    Parameters
    ----------
    model:
        Backend name. ``"circular"`` (default) returns analytic states for the
        planets in :data:`cyclerfinder.core.constants.PLANETS`. ``"astropy"``
        is reserved for the M6 JPL backend and currently raises
        :class:`NotImplementedError`.
    """

    def __init__(self, model: str = "circular") -> None:
        if model == "circular":
            self._backend: _Backend = _CircularBackend()
        elif model == "astropy":
            raise NotImplementedError("astropy backend lands in M6; use model='circular' for now")
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
