"""Smoke tests for cyclerfinder.core.constants.

These guard against typos and gross corruption of the physical constants.
They are deliberately loose; downstream physics tests (M1+) provide the real
numerical anchors against published values.
"""

from __future__ import annotations

from cyclerfinder.core import constants


def test_au_matches_iau() -> None:
    assert constants.AU_KM == 1.49597870700e8  # exact IAU 2012


def test_sun_gm_in_range() -> None:
    # IAU 2015 nominal, allow +/-0.1% drift if a later source is adopted
    assert 1.326e11 < constants.MU_SUN_KM3_S2 < 1.328e11


def test_mars_sma_in_au_range() -> None:
    assert 1.50 < constants.PLANETS["M"].sma_au < 1.55


def test_planet_codes_unique_and_match_names() -> None:
    codes = [p.code for p in constants.PLANETS.values()]
    assert len(codes) == len(set(codes))
    for code, p in constants.PLANETS.items():
        assert code == p.code


def test_all_planet_gms_positive() -> None:
    assert all(p.mu_km3_s2 > 0 for p in constants.PLANETS.values())


def test_safe_perihelion_above_planet_radius() -> None:
    for code, p in constants.PLANETS.items():
        assert constants.SAFE_PERIHELION_KM[code] > p.radius_eq_km
