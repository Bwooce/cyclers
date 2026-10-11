"""Tier-1 Phase 2: planet-centred circular moon ephemeris (plan Phase 2 Task 2.0)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import SATELLITES


def test_europa_state_is_about_jupiter_at_its_sma() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    r, v = ephem.state("Europa", 0.0)
    # |r| ~ Europa's SMA about Jupiter (km), NOT a heliocentric AU-scaled radius.
    assert np.linalg.norm(r) == pytest.approx(SATELLITES["Europa"].sma_km, rel=1e-9)
    # circular speed |v| = sqrt(mu_jup / a); ~13.7 km/s.
    assert 13.0 < np.linalg.norm(v) < 14.5


def test_heliocentric_default_unchanged() -> None:
    # No center -> the existing heliocentric circular backend, byte-identical.
    helio = Ephemeris(model="circular")
    r, _ = helio.state("E", 0.0)
    assert 1.40e8 < np.linalg.norm(r) < 1.55e8  # ~1 AU in km


def test_pluto_small_moon_rides_the_sourced_mean_motion_circle() -> None:
    """#1022d: Nix's circle has the Brozovic-Jacobson 2024 Table 11 mean motion (14.484 deg/day)
    and its Kepler radius, so the speed is sqrt(mu / r)."""
    ephem = Ephemeris(model="circular", center="Pluto")
    mu_system = 975.5
    n_rad_s = np.radians(14.484) / 86400.0  # literal from the source, not from the registry
    a_expected = (mu_system / n_rad_s**2) ** (1.0 / 3.0)
    r, v = ephem.state("Nix", 0.0)
    assert np.linalg.norm(r) == pytest.approx(a_expected, rel=1e-9)
    assert np.linalg.norm(v) == pytest.approx(np.sqrt(mu_system / a_expected), rel=1e-9)
    # a quarter of a period later it is at 90 degrees
    t_quarter = (np.pi / 2.0) / n_rad_s
    r2, _ = ephem.state("Nix", t_quarter)
    assert r2[0] == pytest.approx(0.0, abs=1e-6 * a_expected)
    assert r2[1] == pytest.approx(a_expected, rel=1e-9)


def test_moon_without_a_sourced_rate_is_unchanged() -> None:
    """Europa keeps radius sma_km and the Kepler rate of sma_km (bit-identical to before #1022d)."""
    ephem = Ephemeris(model="circular", center="Jupiter")
    r, v = ephem.state("Europa", 12345.0)
    sat = SATELLITES["Europa"]
    n = sat.mean_motion_deg_day * (np.pi / 180.0) / 86400.0
    assert np.linalg.norm(r) == pytest.approx(sat.sma_km, rel=1e-15)
    assert np.linalg.norm(v) == pytest.approx(sat.sma_km * n, rel=1e-15)
