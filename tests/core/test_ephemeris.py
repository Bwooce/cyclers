"""Tests for :mod:`cyclerfinder.core.ephemeris` (circular backend)."""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris


def test_circular_earth_period() -> None:
    """Propagating one Julian year via mean motion returns Earth close to start.

    Earth's mean motion in :data:`PLANETS` is derived from its semi-major axis,
    so the orbital period is not exactly 365.25 d — it's the Keplerian period
    for a = 1.00000261 AU around the Sun, ~365.26 d. After exactly one such
    period, the position closes to within < 1 km.
    """
    eph = Ephemeris(model="circular")
    earth = PLANETS["E"]
    a_km = earth.sma_au * AU_KM
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)

    r0, v0 = eph.state("E", 0.0)
    r1, v1 = eph.state("E", period_s)

    assert float(np.linalg.norm(r1 - r0)) < 1.0
    assert float(np.linalg.norm(v1 - v0)) < 1.0e-6


def test_circular_speeds() -> None:
    """Each planet's heliocentric speed equals sqrt(mu / r) for a circular orbit."""
    eph = Ephemeris(model="circular")
    for body, data in PLANETS.items():
        r, v = eph.state(body, 0.0)
        r_n = float(np.linalg.norm(r))
        v_n = float(np.linalg.norm(v))
        # Expected circular speed.
        v_circ = sqrt(MU_SUN_KM3_S2 / r_n)
        # Tolerance: 1 m/s.
        assert abs(v_n - v_circ) < 1.0e-3, f"{data.name}: |v|={v_n:.6f}, expected {v_circ:.6f} km/s"


def test_planar() -> None:
    """Circular backend keeps every planet exactly in the ecliptic plane."""
    eph = Ephemeris(model="circular")
    for body in PLANETS:
        for t_days in (0.0, 73.0, 365.0, 1000.0):
            r, v = eph.state(body, t_days * SECONDS_PER_DAY)
            assert r[2] == 0.0
            assert v[2] == 0.0


def test_astropy_backend_constructs() -> None:
    """The astropy backend was implemented in the 2026-06-01 launch-windows slice.

    Was previously expected to raise NotImplementedError (M6 deferral); the
    backend now instantiates cleanly. Detailed behaviour tests live in
    ``tests/core/test_ephemeris_astropy.py``.
    """
    ephem = Ephemeris(model="astropy")
    assert ephem.model == "astropy"


def test_unknown_body() -> None:
    """Unknown body codes propagate the dict's :class:`KeyError`."""
    eph = Ephemeris(model="circular")
    with pytest.raises(KeyError):
        eph.state("Pluto", 0.0)


def test_unknown_model() -> None:
    """Constructing with an unknown model name is a :class:`ValueError`."""
    with pytest.raises(ValueError):
        Ephemeris(model="spice")


def test_t_zero_at_x_axis() -> None:
    """Epoch convention: at t=0 every planet sits on the +x axis."""
    eph = Ephemeris(model="circular")
    for body, data in PLANETS.items():
        r, _ = eph.state(body, 0.0)
        assert r[1] == 0.0
        assert r[0] == pytest.approx(data.sma_au * AU_KM, rel=1e-15)
