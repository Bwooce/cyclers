"""Tier-1 Phase 1: SatelliteData shape + Kepler-III mean motion (plan Phase 1)."""

from __future__ import annotations

import math

from cyclerfinder.core.satellites import SatelliteData, mean_motion_deg_day_about


def test_satellitedata_fields_present() -> None:
    s = SatelliteData(
        name="Europa",
        code="Europa",
        primary="Jupiter",
        mu_km3_s2=3203.0,
        radius_eq_km=1560.8,
        sma_km=671100.0,
        mean_motion_deg_day=mean_motion_deg_day_about(671100.0, mu_primary=1.26686534e8),
        safe_alt_km=100.0,
    )
    assert s.code == "Europa"  # full-name scheme (data/README.md:65-67)
    assert s.primary == "Jupiter"
    assert s.sma_km == 671100.0


def test_mean_motion_matches_kepler_third_law() -> None:
    # n = 360 / (2 pi sqrt(a^3/mu) / 86400) deg/day; Europa about Jupiter.
    mu_jup = 1.26686534e8
    a = 671100.0
    period_s = 2.0 * math.pi * math.sqrt(a**3 / mu_jup)
    expected = 360.0 / (period_s / 86400.0)
    assert mean_motion_deg_day_about(a, mu_primary=mu_jup) == expected
