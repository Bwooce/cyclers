"""Two-arc free-return CHAIN — mechanics gates + the sourced-anchor probe (#163).

The chain (:mod:`cyclerfinder.search.free_return_chain`) builds Russell's actual
generic-return-arc construction: two distinct Earth-to-Earth free-return arcs
(``g`` + ``G``) patched at an intermediate Earth flyby, each crossing Mars's radius.
Its objective is ANCHOR-RESPECTING (the residual IS ``emerged V_inf - sourced
anchor`` plus the per-arc descriptor-ToF match), reusing the
:mod:`cyclerfinder.search.free_return` primitive that CLOSED the symmetric rows —
deliberately NOT the dV-budget objective of :mod:`dsm_leg` that floored at ~9 km/s.

GOLDEN DISCIPLINE: the SOURCED anchor is EXPECTED; the emerged V_inf is EVIDENCE,
never imposed. The mechanics tests assert OUR construction's internal consistency
(degenerate reduction, symmetric reproduction) — non-golden, pinned to guard against
regression. The probe compares emerged-vs-sourced and applies the three-way gate; a
clean EMPTY-SET / AMBIGUOUS is a SUCCESS, not a softening target.

See ``docs/notes/2026-06-08-free-return-chain-results.md``.
"""

from __future__ import annotations

import pytest

from cyclerfinder.search.free_return import free_return_geometry
from cyclerfinder.search.free_return_chain import (
    _arc_ee_time_years,
    _best_n_rev,
    free_return_chain_correct,
    single_arc_degenerate,
)

# Sourced anchors (EXPECTED side; never imposed in the solve).
# russell-ch4-6.44Gg3 (Russell 2004 Table 4.13): aphel 1.54, g(2.087) + G(4.3191),
# v_inf E = 6.44, M = 3.74.
_644_APHELION = 1.54
_644_G_TOF = 2.087
_644_BIGG_TOF = 4.3191
_644_VINF_E = 6.44
_644_VINF_M = 3.74
# russell-ch4-4.991gG2 / S1L1 (Russell 2004 Table 4.9): aphel 1.64,
# g(1.4612) + G(2.8096), v_inf E = 4.99, M = 5.10.
_4991_APHELION = 1.64
_4991_G_TOF = 1.4612
_4991_BIGG_TOF = 2.8096
_4991_VINF_E = 4.99
_4991_VINF_M = 5.10


# ---------------------------------------------------------------------------
# Mechanics gate 1 — multi-rev arc-time arithmetic.
# ---------------------------------------------------------------------------


def test_arc_ee_time_increases_by_one_period_per_rev() -> None:
    """Each added revolution adds exactly one orbital period to the arc ToF."""
    a, e = 1.223, 0.2602
    t0 = _arc_ee_time_years(a, e, 0)
    t1 = _arc_ee_time_years(a, e, 1)
    g = free_return_geometry(a, e)
    period_years = g.period_days / 365.25
    assert t1 - t0 == pytest.approx(period_years, rel=1e-9)


def test_best_n_rev_selects_closest() -> None:
    """The discrete n_rev DOF picks the revolution count nearest the target ToF."""
    a, e = 1.223, 0.2602
    # n=1 -> ~2.43 yr; targeting 2.4 yr must pick n=1, targeting 3.8 yr -> n=2.
    assert _best_n_rev(a, e, 2.4) == 1
    assert _best_n_rev(a, e, 3.8) == 2


# ---------------------------------------------------------------------------
# Mechanics gate 2 — degenerate single-arc reduction == free_return.
# ---------------------------------------------------------------------------


def test_degenerate_reduces_to_single_free_return_ellipse() -> None:
    """With arc1 == arc2 the chain collapses to ONE free-return ellipse: the two
    arcs converge to the same (a, e), and the intermediate-flyby continuity and
    turn are exactly zero (no patch needed)."""
    d = single_arc_degenerate(_644_APHELION, _644_G_TOF, _644_VINF_E, _644_VINF_M)
    a1, e1 = d.arcs[0].a_au, d.arcs[0].e
    a2, e2 = d.arcs[1].a_au, d.arcs[1].e
    assert a1 == pytest.approx(a2, abs=1e-6)
    assert e1 == pytest.approx(e2, abs=1e-6)
    # Continuity and turn vanish for coincident arcs (the degenerate patch).
    assert d.vinf_continuity_kms < 1e-6
    assert d.intermediate_turn_deg < 1e-6
    assert d.intermediate_flyby_feasible
    # And the emerged geometry matches a direct single free_return_geometry call.
    g = free_return_geometry(a1, e1)
    assert d.arcs[0].vinf_e == pytest.approx(g.vinf["E"], rel=1e-9)
    assert d.arcs[0].vinf_m == pytest.approx(g.vinf["M"], rel=1e-9)


# ---------------------------------------------------------------------------
# Mechanics gate 3 — distinct-arc binding (no single-ellipse collapse).
# ---------------------------------------------------------------------------


def test_two_arc_descriptors_drive_distinct_revolution_counts() -> None:
    """The g- and G-arc descriptor ToFs (different by construction) force DIFFERENT
    revolution counts, so the chain is a genuine two-arc object, not a single
    ellipse wearing two hats. (Pins the ToF-binding that forbids collapse.)"""
    r = free_return_chain_correct(
        _644_APHELION, _644_G_TOF, _644_BIGG_TOF, _644_VINF_E, _644_VINF_M
    )
    assert r.arcs[0].n_rev != r.arcs[1].n_rev
    # The emerged arc ToFs straddle their distinct descriptor targets.
    assert r.arcs[0].arc_tof_years < r.arcs[1].arc_tof_years


# ---------------------------------------------------------------------------
# THE PROBE — sourced-anchor gate, three-way verdict (@slow).
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_644gg3_two_arc_chain_sourced_anchor_probe() -> None:
    """6.44Gg3 two-arc chain at the SOURCED anchors (E 6.44 / M 3.74).

    The decisive result (see the results note): the V_inf anchors ARE reachable
    (vinf_res << 0.5 km/s at every encounter) and the intermediate Earth turn is
    trivially bend-feasible, but the per-arc DESCRIPTOR ToFs are NOT simultaneously
    reachable — the arc ToF is quantised by the discrete revolution count and the
    V_inf-fixed ellipse period leaves the g-arc (2.087 yr) and G-arc (4.3191 yr)
    targets in the gaps. This is the AMBIGUOUS regime: far below the dsm_leg ~9 and
    the single-arc 3.01/3.06 floors on the V_inf axis, but not a full structural
    closure. EMPTY-SET-flavoured on the ToF axis. Pinned as evidence, not golden."""
    r = free_return_chain_correct(
        _644_APHELION, _644_G_TOF, _644_BIGG_TOF, _644_VINF_E, _644_VINF_M
    )
    # EVIDENCE: the emerged V_inf is at the sourced anchors (the basin exists).
    assert r.vinf_residual_kms < 0.5
    assert abs(r.arcs[0].vinf_m - _644_VINF_M) < 0.5
    assert abs(r.arcs[1].vinf_m - _644_VINF_M) < 0.5
    assert abs(r.arcs[0].vinf_e - _644_VINF_E) < 0.5
    assert abs(r.arcs[1].vinf_e - _644_VINF_E) < 0.5
    # Beat the prior floors decisively on the V_inf axis (single-arc 3.01/3.06,
    # dsm_leg ~9) — emerged Mars V_inf is at 3.74, not 3.06 or 8-16.
    assert r.vinf_residual_kms < 3.0
    # The intermediate Earth turn fits the ballistic cone (bend is not the wall).
    assert r.intermediate_turn_deg < r.intermediate_max_turn_deg
    # But the STRUCTURAL descriptor ToF does NOT close: the g-arc target 2.087 yr is
    # unreachable below the n_rev=1 minimum, leaving a multi-tenth-year residual.
    # This is the EMPTY-SET-flavoured half — assert the regime, not a manufactured 0.
    assert r.tof_residual_years > 0.2
    # NO catalogue writeback; this test asserts evidence only.


@pytest.mark.slow
def test_4991gg2_two_arc_chain_sourced_anchor_probe() -> None:
    """S1L1 / 4.991gG2 two-arc chain at its OWN sourced anchors (E 4.99 / M 5.10).

    This row's anchors AND descriptor ToFs are both nearly reachable: emerged V_inf
    matches to <0.1 km/s and the G-arc ToF (2.8096 yr) lands essentially exactly at
    n_rev=1; only the g-arc (1.4612 yr) carries a ~0.14 yr residual (its n_rev=0
    arc time is ~1.32 yr). The closest the two-arc chain comes to a full structural
    closure of a multi-arc row — the continuation-seed candidate."""
    r = free_return_chain_correct(
        _4991_APHELION, _4991_G_TOF, _4991_BIGG_TOF, _4991_VINF_E, _4991_VINF_M
    )
    assert r.vinf_residual_kms < 0.5
    assert abs(r.arcs[1].vinf_m - _4991_VINF_M) < 0.5
    assert abs(r.arcs[1].vinf_e - _4991_VINF_E) < 0.5
    # The G-arc ToF lands near-exactly on its descriptor (the strong half).
    assert abs(r.arcs[1].arc_tof_years - _4991_BIGG_TOF) < 0.1
    # The intermediate Earth turn fits the cone.
    assert r.intermediate_turn_deg < r.intermediate_max_turn_deg
