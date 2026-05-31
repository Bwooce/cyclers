"""Smoke tests for the astropy ephemeris backend.

The 2026-06-01 launch-windows slice brought astropy forward from full M6.
These tests verify the backend returns sensible heliocentric positions; they
do NOT cross-validate against an independent source — that's M6 work.
"""

from __future__ import annotations

from datetime import UTC, datetime

import numpy as np

from cyclerfinder.core.constants import AU_KM, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris

J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def _seconds_since_j2000(dt: datetime) -> float:
    return (dt - J2000).total_seconds()


def test_astropy_backend_instantiates() -> None:
    ephem = Ephemeris(model="astropy")
    assert ephem.model == "astropy"


def test_astropy_earth_position_is_near_1_au() -> None:
    """Earth's heliocentric distance must be ~1 AU on any date (eccentricity ~0.017)."""
    ephem = Ephemeris(model="astropy")
    # Pick three dates across the year to exercise eccentricity.
    for dt in (
        datetime(2026, 1, 1, tzinfo=UTC),
        datetime(2026, 4, 1, tzinfo=UTC),
        datetime(2026, 7, 1, tzinfo=UTC),
    ):
        r, _ = ephem.state("E", _seconds_since_j2000(dt))
        distance_au = float(np.linalg.norm(r)) / AU_KM
        # Earth perihelion ~0.983 AU, aphelion ~1.017 AU.
        assert 0.97 < distance_au < 1.03, f"Earth distance {distance_au:.4f} AU on {dt!r}"


def test_astropy_mars_position_is_in_mars_orbit() -> None:
    """Mars's heliocentric distance must lie in [1.38, 1.67] AU (perihelion to aphelion)."""
    ephem = Ephemeris(model="astropy")
    for dt in (
        datetime(2026, 1, 1, tzinfo=UTC),
        datetime(2027, 1, 1, tzinfo=UTC),
        datetime(2028, 1, 1, tzinfo=UTC),
    ):
        r, _ = ephem.state("M", _seconds_since_j2000(dt))
        distance_au = float(np.linalg.norm(r)) / AU_KM
        # Mars perihelion 1.381 AU, aphelion 1.666 AU.
        assert 1.35 < distance_au < 1.70, f"Mars distance {distance_au:.4f} AU on {dt!r}"


def test_astropy_venus_position_is_in_venus_orbit() -> None:
    """Venus heliocentric distance must lie near 0.72 AU."""
    ephem = Ephemeris(model="astropy")
    for dt in (
        datetime(2026, 1, 1, tzinfo=UTC),
        datetime(2026, 6, 1, tzinfo=UTC),
    ):
        r, _ = ephem.state("V", _seconds_since_j2000(dt))
        distance_au = float(np.linalg.norm(r)) / AU_KM
        # Venus's orbit is nearly circular; allow generous tolerance.
        assert 0.71 < distance_au < 0.73, f"Venus distance {distance_au:.4f} AU on {dt!r}"


def test_astropy_earth_orbital_speed_matches_kepler() -> None:
    """|v_Earth| ≈ sqrt(mu_Sun / a_Earth) ≈ 29.78 km/s."""
    ephem = Ephemeris(model="astropy")
    _r, v = ephem.state("E", _seconds_since_j2000(datetime(2026, 3, 21, tzinfo=UTC)))
    speed = float(np.linalg.norm(v))
    # Mean orbital speed of Earth; tolerance allows for eccentricity-driven variation.
    assert 29.0 < speed < 30.5, f"Earth speed {speed:.3f} km/s out of range"


def test_circular_backend_still_works() -> None:
    """Adding the astropy backend must not break the M1 circular backend."""
    ephem = Ephemeris(model="circular")
    r, _v = ephem.state("E", 0.0)
    a_au = float(np.linalg.norm(r)) / AU_KM
    assert abs(a_au - PLANETS["E"].sma_au) < 1e-9


def test_unknown_model_raises() -> None:
    import pytest

    with pytest.raises(ValueError, match="unknown ephemeris model"):
        Ephemeris(model="not-a-real-backend")
