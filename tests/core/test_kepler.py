"""Tests for :mod:`cyclerfinder.core.kepler` (universal-variable propagator)."""

from __future__ import annotations

from math import pi, sqrt

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate


def test_kepler_zero_dt_identity() -> None:
    """``dt == 0`` returns the original state exactly (as a copy)."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    r, v = propagate(r0, v0, 0.0)
    assert np.array_equal(r, r0)
    assert np.array_equal(v, v0)
    # And it's a copy, not a view.
    r[0] = 12345.0
    assert r0[0] != 12345.0


def test_kepler_circular_period() -> None:
    """Earth's circular state propagates one Keplerian period back to itself."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    a_km = float(np.linalg.norm(r0))
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
    r1, v1 = propagate(r0, v0, period_s)
    # 1 km positional tolerance per plan §4.6.
    assert float(np.linalg.norm(r1 - r0)) < 1.0
    assert float(np.linalg.norm(v1 - v0)) < 1.0e-6


def test_kepler_reversibility() -> None:
    """``propagate(propagate(s, +dt), -dt) == s`` within numerical noise."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("M", 0.0)
    dt = 137.0 * SECONDS_PER_DAY
    r1, v1 = propagate(r0, v0, dt)
    r2, v2 = propagate(r1, v1, -dt)
    assert float(np.linalg.norm(r2 - r0)) < 1.0
    assert float(np.linalg.norm(v2 - v0)) < 1.0e-6


def test_kepler_energy_conservation() -> None:
    """Specific orbital energy is conserved to high relative precision."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    energy0 = float(np.dot(v0, v0)) / 2.0 - MU_SUN_KM3_S2 / float(np.linalg.norm(r0))
    for t_days in (10.0, 100.0, 365.0, 700.0):
        r, v = propagate(r0, v0, t_days * SECONDS_PER_DAY)
        energy = float(np.dot(v, v)) / 2.0 - MU_SUN_KM3_S2 / float(np.linalg.norm(r))
        # Relative tolerance 1e-8 per plan §4.6.
        assert abs(energy - energy0) / abs(energy0) < 1.0e-8


def test_kepler_hyperbolic() -> None:
    """A fabricated hyperbolic state round-trips identity under +dt then -dt."""
    # Start at 1 AU with velocity 1.4x escape (well into the hyperbolic regime).
    r0 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    v_esc = sqrt(2.0 * MU_SUN_KM3_S2 / AU_KM)
    v0 = np.array([0.0, 1.4 * v_esc, 0.0], dtype=np.float64)
    dt = 30.0 * SECONDS_PER_DAY
    r1, v1 = propagate(r0, v0, dt)
    r2, v2 = propagate(r1, v1, -dt)
    assert float(np.linalg.norm(r2 - r0)) < 1.0
    assert float(np.linalg.norm(v2 - v0)) < 1.0e-6


def test_kepler_negative_dt_directly() -> None:
    """Backward propagation from a future state lands at the past state."""
    eph = Ephemeris(model="circular")
    r_at_zero, v_at_zero = eph.state("E", 0.0)
    dt = 200.0 * SECONDS_PER_DAY
    r_at_dt, v_at_dt = eph.state("E", dt)
    # Propagate the future state backward by dt; should land at t=0.
    r_back, v_back = propagate(r_at_dt, v_at_dt, -dt)
    assert float(np.linalg.norm(r_back - r_at_zero)) < 1.0
    assert float(np.linalg.norm(v_back - v_at_zero)) < 1.0e-6
