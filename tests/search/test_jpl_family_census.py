"""Tests for the JPL periodic-orbit DISCOVERY-INPUT mining module (#667).

No network: :func:`fetch_family_window` is exercised with an injected fake
``query_fn`` (same convention ``tests/search/test_jpl_family_check.py`` already
uses); :func:`propagate_min_distances_km`/:func:`classify_secondary_approach`
are exercised against small SYNTHETIC states with hand-computable expected
distances (arithmetic/scaling unit tests of this module's OWN new code, not a
claim about any real orbit's physical properties -- the golden-tests-sourced-
only discipline applies to claims about real published orbits, not to unit
tests of a helper function's internal distance/unit-conversion arithmetic).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

from cyclerfinder.genome.hill_screen import PASS_HILL_FRACTION, earth_sun_hill_radius_km
from cyclerfinder.search.jpl_family_census import (
    CLOSE_APPROACH_HILL_FRACTION,
    classify_secondary_approach,
    fetch_family_window,
    hill_radius_km,
    propagate_min_distances_km,
)
from cyclerfinder.verify.jpl_periodic_orbits import JplOrbit, JplSystemConstants

# ---------------------------------------------------------------------------
# fetch_family_window (no network -- injected fake query_fn)
# ---------------------------------------------------------------------------


def _fake_query_fn(captured: dict[str, Any]) -> Any:
    def _fn(system: str, family: str, **kwargs: Any) -> tuple[Any, Any]:
        captured["system"] = system
        captured["family"] = family
        captured.update(kwargs)
        constants = JplSystemConstants(
            name="fake",
            mu=0.01,
            lunit_km=1000.0,
            tunit_s=1000.0,
            radius_secondary_km=10.0,
            libration_points={},
        )
        return constants, []

    return _fn


def test_fetch_family_window_rejects_unsupported_system() -> None:
    with pytest.raises(ValueError, match="does not catalog the 'sun-jupiter' system"):
        fetch_family_window("sun-jupiter", "dro", query_fn=_fake_query_fn({}))


def test_fetch_family_window_rejects_unsupported_family() -> None:
    with pytest.raises(ValueError, match="does not catalog the 'not-a-family' family"):
        fetch_family_window("mars-phobos", "not-a-family", query_fn=_fake_query_fn({}))


def test_fetch_family_window_requires_libr_for_halo() -> None:
    with pytest.raises(ValueError, match="requires a libr"):
        fetch_family_window("mars-phobos", "halo", branch="N", query_fn=_fake_query_fn({}))


def test_fetch_family_window_requires_branch_for_halo() -> None:
    with pytest.raises(ValueError, match="requires a branch"):
        fetch_family_window("mars-phobos", "halo", libr=1, query_fn=_fake_query_fn({}))


def test_fetch_family_window_passes_bulk_range_params_through() -> None:
    captured: dict[str, Any] = {}
    fetch_family_window(
        "Mars-Phobos",
        "DRO",
        jacobi_min=2.5,
        jacobi_max=3.2,
        period_min=0.1,
        period_max=10.0,
        query_fn=_fake_query_fn(captured),
    )
    assert captured["system"] == "mars-phobos"
    assert captured["family"] == "dro"
    assert captured["jacobimin"] == 2.5
    assert captured["jacobimax"] == 3.2
    assert captured["periodmin"] == 0.1
    assert captured["periodmax"] == 10.0
    assert captured["periodunits"] == "TU"


def test_fetch_family_window_forwards_libr_and_branch() -> None:
    captured: dict[str, Any] = {}
    fetch_family_window(
        "mars-phobos", "halo", libr=1, branch="N", query_fn=_fake_query_fn(captured)
    )
    assert captured["libr"] == 1
    assert captured["branch"] == "N"


# ---------------------------------------------------------------------------
# hill_radius_km
# ---------------------------------------------------------------------------


def test_hill_radius_km_matches_hill_screen_earth_sun_case() -> None:
    """Cross-check against the existing, already-vetted
    ``hill_screen.earth_sun_hill_radius_km`` -- same formula
    (``a*(GM2/(3*GM1))**(1/3) == a*(mu/(3*(1-mu)))**(1/3)``), generalized
    from hardcoded Earth-Sun constants to an arbitrary ``(mu, l_km)`` pair."""
    from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
    from cyclerfinder.core.satellites import PRIMARIES

    a_es_km = float(PLANETS["E"].sma_au) * AU_KM
    gm_emb = float(PRIMARIES["Earth"])
    mu_emb_about_sun = gm_emb / (gm_emb + MU_SUN_KM3_S2)

    expected = earth_sun_hill_radius_km()
    got = hill_radius_km(mu_emb_about_sun, a_es_km)
    assert got == pytest.approx(expected, rel=1e-9)


def test_hill_radius_km_rejects_invalid_mu() -> None:
    with pytest.raises(ValueError, match="mu must be"):
        hill_radius_km(0.0, 1000.0)
    with pytest.raises(ValueError, match="mu must be"):
        hill_radius_km(1.5, 1000.0)


def test_hill_radius_km_rejects_nonpositive_l_km() -> None:
    with pytest.raises(ValueError, match="l_km must be positive"):
        hill_radius_km(0.01, 0.0)


# ---------------------------------------------------------------------------
# propagate_min_distances_km -- synthetic, hand-computable expected distances
# ---------------------------------------------------------------------------


def _tiny_period_orbit(
    state0: np.ndarray, mu: float, lunit_km: float
) -> tuple[JplOrbit, JplSystemConstants]:
    """A JplOrbit with a vanishingly short 'period' so DOP853 barely moves the
    state -- isolates the distance formula + km scaling from any real
    dynamics, since over dt=1e-9 (nondim CR3BP time units) the velocity
    change is O(1e-9), far below any tolerance used below."""
    orbit = JplOrbit(
        state0=np.asarray(state0, dtype=np.float64), jacobi=3.0, period=1e-9, stability=1.0
    )
    constants = JplSystemConstants(
        name="synthetic",
        mu=mu,
        lunit_km=lunit_km,
        tunit_s=1000.0,
        radius_secondary_km=None,
        libration_points={},
    )
    return orbit, constants


def test_propagate_min_distances_km_secondary_offset() -> None:
    mu = 0.1
    lunit_km = 1000.0
    dx, dy, dz = 0.01, 0.02, 0.03
    state0 = np.array([1.0 - mu + dx, dy, dz, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)

    _d_primary_km, d_secondary_km = propagate_min_distances_km(orbit, constants)
    expected_km = math.sqrt(dx**2 + dy**2 + dz**2) * lunit_km
    assert d_secondary_km == pytest.approx(expected_km, rel=1e-3)


def test_propagate_min_distances_km_primary_offset() -> None:
    mu = 0.1
    lunit_km = 1000.0
    dx, dy, dz = 0.02, -0.01, 0.0
    state0 = np.array([-mu + dx, dy, dz, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)

    d_primary_km, _d_secondary_km = propagate_min_distances_km(orbit, constants)
    expected_km = math.sqrt(dx**2 + dy**2 + dz**2) * lunit_km
    assert d_primary_km == pytest.approx(expected_km, rel=1e-3)


# ---------------------------------------------------------------------------
# classify_secondary_approach
# ---------------------------------------------------------------------------


def test_classify_secondary_approach_close_and_valid() -> None:
    mu = 0.01
    lunit_km = 1000.0
    r_hill = hill_radius_km(mu, lunit_km)
    # Sit well inside the close-approach ceiling (10% of the Hill radius, and
    # well clear of a tiny 1 km secondary radius).
    offset_km = 0.1 * r_hill
    dx = offset_km / lunit_km
    state0 = np.array([1.0 - mu + dx, 0.0, 0.0, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)

    verdict = classify_secondary_approach(
        orbit, constants, system="mars-phobos", family="dro", radius_secondary_km_override=1.0
    )
    assert verdict.physically_valid is True
    assert verdict.is_close_approach is True
    assert verdict.hill_fraction == pytest.approx(0.1, rel=1e-2)
    assert verdict.system == "mars-phobos"
    assert verdict.family == "dro"


def test_classify_secondary_approach_remote_fails_close_ceiling() -> None:
    mu = 0.01
    lunit_km = 1000.0
    r_hill = hill_radius_km(mu, lunit_km)
    # Sit well BEYOND the close-approach ceiling (5x the Hill radius).
    offset_km = 5.0 * r_hill
    dx = offset_km / lunit_km
    state0 = np.array([1.0 - mu + dx, 0.0, 0.0, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)

    verdict = classify_secondary_approach(
        orbit, constants, system="mars-phobos", family="halo", radius_secondary_km_override=1.0
    )
    assert verdict.physically_valid is True
    assert verdict.is_close_approach is False
    assert verdict.hill_fraction > CLOSE_APPROACH_HILL_FRACTION


def test_classify_secondary_approach_inside_body_fails_physically() -> None:
    mu = 0.01
    lunit_km = 1000.0
    # 1 km from secondary centre, but the (overridden) body radius is 50 km.
    dx = 1.0 / lunit_km
    state0 = np.array([1.0 - mu + dx, 0.0, 0.0, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)

    verdict = classify_secondary_approach(
        orbit, constants, system="mars-phobos", family="dro", radius_secondary_km_override=50.0
    )
    assert verdict.physically_valid is False
    assert verdict.is_close_approach is False


def test_classify_secondary_approach_unsourced_radius_is_honest_nonclaim() -> None:
    mu = 0.01
    lunit_km = 1000.0
    state0 = np.array([1.0 - mu + 0.001, 0.0, 0.0, 0.0, 0.0, 0.0])
    orbit, constants = _tiny_period_orbit(state0, mu, lunit_km)
    assert constants.radius_secondary_km is None

    verdict = classify_secondary_approach(orbit, constants, system="mars-phobos", family="dro")
    assert verdict.radius_secondary_km is None
    assert verdict.physically_valid is True
    assert "unevaluated" in verdict.notes.lower() or "not evaluated" in verdict.notes.lower()


def test_close_approach_hill_fraction_matches_hill_screen_constant() -> None:
    """Documents the deliberate reuse-by-analogy this module's docstring
    describes -- if hill_screen's own PASS_HILL_FRACTION ever changes, this
    module's default should be a conscious decision, not a silent drift."""
    assert CLOSE_APPROACH_HILL_FRACTION == PASS_HILL_FRACTION
