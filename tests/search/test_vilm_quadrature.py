"""Tier-1 Phase 5 GOLDEN: VILM ΔV-min quadrature vs published Part-1 Table 1
(no-GA, mining note A2, lines 358-381) + Table 2 (with-GA, A3) + Europa A6 scalar.

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Table 1/2 ΔV_min (km/s) and the A6
worked scalar. Model caveat (mining note 491-493): Part-1 ΔV are linked-conic,
upper-bound-ish vs CR3BP — the band is a tolerance, never equality, and a
CR3BP-disagreeing value <=10% is NOT a rejection."""

from __future__ import annotations

import pytest

from cyclerfinder.search.vilm import europa_endgame_dv, vilm_dv_floor, vilm_dv_min

# (moon_a, moon_b, ΔV_min km/s) — Part-1 Table 1 (note lines 362-378).
_DV_MIN = [
    ("Callisto", "Ganymede", 1.81),
    ("Ganymede", "Europa", 1.71),
    ("Europa", "Io", 1.76),
    ("Titan", "Rhea", 1.15),
    ("Rhea", "Dione", 0.52),
    ("Tethys", "Enceladus", 0.34),
]


@pytest.mark.parametrize("a,b,dv", _DV_MIN)
def test_quadrature_dv_min_matches_published_table1(a: str, b: str, dv: float) -> None:
    # 10% band per the linked-conic vs CR3BP model caveat (mining note 491-493).
    assert vilm_dv_min(a, b) == pytest.approx(dv, rel=0.10)


def test_ga_routed_dv_min_matches_table2() -> None:
    # Callisto-G-Europa ΔV_min 1.61 (Part-1 Table 2, note line 389).
    assert vilm_dv_min("Callisto", "Europa", via=("Ganymede",)) == pytest.approx(1.61, rel=0.10)


def test_europa_3vilm_endgame_scalar() -> None:
    # A6 (note 436-438): Europa endgame V_inf_bar 1.8 -> 0.77 km/s. The published
    # DISCRETE 3-VILM design costs 154 m/s over ~46 days; the CR3BP re-optimised
    # long-transfer is 147 m/s (note line 439). europa_endgame_dv() computes the
    # CONTINUOUS-VILM Eq.(13) THEORETICAL MINIMUM over the same V_inf bounds — a
    # valid LOWER BOUND on any finite-VILM design (a discrete sequence cannot beat
    # the infinite-VILM floor). So the computed floor must be <= the published
    # discrete scalar; it sits in the physically-expected band below 154.
    dv_ms, days = europa_endgame_dv()
    # Continuous floor is below the discrete 154 m/s design (discrete > floor)...
    assert dv_ms <= 154.0
    # ...and within the linked-conic band of it (not arbitrarily low) — the
    # 3-VILM design overshoots the floor by the finite-VILM penalty (~17%).
    assert dv_ms == pytest.approx(154.0, rel=0.20)
    # Duration is the published phasing scalar (the phase-free quadrature does not
    # predict ToF), carried verbatim from A6.
    assert days == pytest.approx(46.0, abs=5.0)


def test_dv_floor_is_admissible_lower_bound() -> None:
    # DEVIATION (task #76): the admissible floor is escape+capture, which is <=
    # EVERY routing. A gravity assist REDUCES ΔV (Table 2 < Table 1), so with-GA
    # < no-GA — the plan's "floor <= with-GA" only holds because the true floor
    # (escape+capture) is below both, not because no-GA is a lower bound. See
    # vilm_dv_floor docstring.
    floor = vilm_dv_floor("Callisto", "Europa")
    no_ga = vilm_dv_min("Callisto", "Europa")
    with_ga = vilm_dv_min("Callisto", "Europa", via=("Ganymede",))
    assert floor <= with_ga <= no_ga  # GA reduces ΔV; floor below both
    assert floor <= no_ga
