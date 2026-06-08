"""Two-arc free-return CHAIN continuation — mechanics gates + the decisive S1L1 probe (#164).

The continuation chain (:mod:`cyclerfinder.search.continuation_chain`) walks the
#163 two-arc free-return seed (:mod:`cyclerfinder.search.free_return_chain`) from the
circular-coplanar planet model out to the real DE440-consistent J2000 eccentric /
inclined model, reusing the #158 ramp machinery
(:mod:`cyclerfinder.search.continuation`). It answers the S1L1/#94 frontier
question: does real planet eccentricity close the g-arc descriptor-ToF gap (~0.14 yr
in circular) while holding both arcs' V_inf at the SOURCED anchors?

GOLDEN DISCIPLINE: the row's OWN SOURCED anchor is EXPECTED; emerged V_inf AND ToF
are evidence, compared non-circularly. A CLOSE must satisfy BOTH halves (V_inf
within ~0.5 of the anchors AND both descriptor ToFs reached) — V_inf alone is the
#163 spurious-collapse trap. The mechanics tests assert OUR construction's internal
consistency (bit-identical lam=0 reduction, continuous first ramp step) — non-golden,
pinned against regression. The probe applies the three-way gate.

See ``docs/notes/2026-06-08-continuation-chain-s1l1-results.md``.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.search.continuation_chain import (
    _ramped_arc_ee_time_years,
    _ramped_arc_geometry,
    _solve_chain_step,
    continuation_chain_correct,
)
from cyclerfinder.search.free_return import free_return_geometry
from cyclerfinder.search.free_return_chain import (
    SECONDS_PER_YEAR,
    _arc_ee_time_years,
    free_return_chain_correct,
)

# S1L1 / 4.991gG2 sourced anchors (Russell 2004 Table 4.9, the row's OWN anchors —
# NOT the CPOM 5.65/3.05 framing): aphel 1.64, g(1.4612) + G(2.8096), E 4.99 / M 5.10.
_4991_APHELION = 1.64
_4991_G_TOF = 1.4612
_4991_BIGG_TOF = 2.8096
_4991_VINF_E = 4.99
_4991_VINF_M = 5.10
# 6.44Gg3 (Russell 2004 Table 4.13): aphel 1.54, g(2.087) + G(4.3191), E 6.44 / M 3.74.
_644_APHELION = 1.54
_644_G_TOF = 2.087
_644_BIGG_TOF = 4.3191
_644_VINF_E = 6.44
_644_VINF_M = 3.74


def _s1l1_seed() -> tuple[float, float, float, float]:
    """The #163 circular two-arc seed (a1, e1, a2, e2) for S1L1."""
    r = free_return_chain_correct(
        _4991_APHELION, _4991_G_TOF, _4991_BIGG_TOF, _4991_VINF_E, _4991_VINF_M
    )
    return r.arcs[0].a_au, r.arcs[0].e, r.arcs[1].a_au, r.arcs[1].e


# ---------------------------------------------------------------------------
# Mechanics gate 1 — lam=0 reproduces the #163 circular seed bit-identically.
# ---------------------------------------------------------------------------


def test_lambda0_ramped_geometry_is_bit_identical_to_circular_chain() -> None:
    """At lam_e = lam_i = 0 the ramped two-arc geometry equals the circular
    free_return_geometry bit-for-bit (V_inf, crossing nu, leg ToF) and the ramped
    arc-time equals the circular arc-time exactly — the bit-identical-seed gate."""
    a1, e1, a2, e2 = _s1l1_seed()
    for a, e in ((a1, e1), (a2, e2)):
        g_circ = free_return_geometry(a, e)
        g_ramp = _ramped_arc_geometry(a, e, 0.0, 0.0, 0.0)
        assert g_ramp.vinf["E"] == g_circ.vinf["E"]
        assert g_ramp.vinf["M"] == g_circ.vinf["M"]
        assert g_ramp.nu["E"] == g_circ.nu["E"]
        assert g_ramp.nu["M"] == g_circ.nu["M"]
        assert g_ramp.tof_em_days == g_circ.tof_em_days
        for n_rev in (0, 1, 2):
            assert _ramped_arc_ee_time_years(a, e, 0.0, 0.0, 0.0, n_rev) == _arc_ee_time_years(
                a, e, n_rev
            )


def test_lambda0_seed_solve_reproduces_163_close_leaning() -> None:
    """The lam=0 seed solve reproduces the #163 CLOSE-LEANING numbers: V_inf at the
    anchors (<0.1 km/s) and the g-arc carrying the ~0.14 yr ToF residual."""
    from cyclerfinder.core.constants import MU_SUN_KM3_S2

    a1, e1, a2, e2 = _s1l1_seed()
    x0 = np.array([a1, e1, a2, e2, 0.0])
    s = _solve_chain_step(
        x0,
        0.0,
        0.0,
        vinf_e_anchor=_4991_VINF_E,
        vinf_m_anchor=_4991_VINF_M,
        arc1_tof_years=_4991_G_TOF,
        arc2_tof_years=_4991_BIGG_TOF,
        mu=MU_SUN_KM3_S2,
        tol_kms=0.5,
    )
    assert s is not None
    assert s.vinf_residual_kms < 0.1
    # g-arc circular ToF ~1.325 yr vs 1.4612 target -> ~0.14 yr residual.
    assert abs(s.arc1.arc_tof_years - _4991_G_TOF) > 0.1
    # G-arc circular ToF near-exact at n_rev=1.
    assert abs(s.arc2.arc_tof_years - _4991_BIGG_TOF) < 0.05


# ---------------------------------------------------------------------------
# Mechanics gate 2 — one e-ramp step moves the solution continuously (no jump).
# ---------------------------------------------------------------------------


def test_one_e_ramp_step_moves_continuously_no_collapse() -> None:
    """A single small e-ramp step (lam_e = 0.1) moves the solution continuously: the
    V_inf residual stays bounded (no off-family jump to ~tens of km/s) and the per-arc
    n_rev ToF-binding term keeps the two arcs DISTINCT (no single-ellipse collapse —
    the #163 spurious-close trap). Guards the frame-consistent-ramp fix."""
    from cyclerfinder.core.constants import MU_SUN_KM3_S2

    a1, e1, a2, e2 = _s1l1_seed()
    x0 = np.array([a1, e1, a2, e2, 0.0])
    s0 = _solve_chain_step(
        x0,
        0.0,
        0.0,
        vinf_e_anchor=_4991_VINF_E,
        vinf_m_anchor=_4991_VINF_M,
        arc1_tof_years=_4991_G_TOF,
        arc2_tof_years=_4991_BIGG_TOF,
        mu=MU_SUN_KM3_S2,
        tol_kms=0.5,
    )
    assert s0 is not None
    x1 = np.array([s0.a1, s0.e1, s0.a2, s0.e2, s0.t0_sec / SECONDS_PER_YEAR])
    s1 = _solve_chain_step(
        x1,
        0.1,
        0.0,
        vinf_e_anchor=_4991_VINF_E,
        vinf_m_anchor=_4991_VINF_M,
        arc1_tof_years=_4991_G_TOF,
        arc2_tof_years=_4991_BIGG_TOF,
        mu=MU_SUN_KM3_S2,
        tol_kms=0.5,
    )
    assert s1 is not None
    # Continuous, NOT a collapse: V_inf residual stays small (the #163 spurious
    # off-family jump would push this to tens of km/s).
    assert s1.vinf_residual_kms < 0.5
    # Arcs stay distinct (the two shapes don't fuse to one ellipse).
    assert s1.arc1.n_rev != s1.arc2.n_rev or abs(s1.a1 - s1.a2) > 1e-4


# ---------------------------------------------------------------------------
# THE DECISIVE PROBE — S1L1 two-arc continuation to DE440, three-way gate (@slow).
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_s1l1_two_arc_continuation_to_ephemeris_closes() -> None:
    """S1L1 / 4.991gG2 two-arc continuation circular -> real ephemeris (#94 probe).

    THE DECISIVE RESULT (see the results note): walking the #163 two-arc seed from
    the circular-coplanar model to the real J2000 eccentric/inclined model CLOSES
    the row — both arcs' V_inf land at the SOURCED anchors (E 4.99 / M 5.10) to
    <0.1 km/s AND both descriptor ToFs are reached (the g-arc 0.14 yr circular gap
    CLOSES to ~0.03 yr, tof_residual < 0.1 yr) AND the intermediate Earth flyby is
    bend-feasible. Real Mars eccentricity at the converged encounter epoch breaks
    the circular integer-rev ToF quantization — exactly the blocker's prediction.

    BOTH halves of the gate must pass (V_inf AND ToF); V_inf alone is the #163
    spurious-collapse trap. Pinned as EVIDENCE (the row's own sourced anchor is the
    EXPECTED side). NO catalogue writeback — this is the continuation seed for a V3
    evidence chain; any promotion is the main session's call after review."""
    a1, e1, a2, e2 = _s1l1_seed()
    res = continuation_chain_correct(
        a1,
        e1,
        a2,
        e2,
        0.0,
        _4991_G_TOF,
        _4991_BIGG_TOF,
        _4991_VINF_E,
        _4991_VINF_M,
        ladder=(1, 3, 9),
    )
    bf = res.best_final
    assert bf is not None
    # V_inf half: both arcs at the sourced anchors.
    assert res.vinf_close
    assert bf.vinf_residual_kms < 0.1
    assert abs(bf.arc1.vinf_e - _4991_VINF_E) < 0.5
    assert abs(bf.arc1.vinf_m - _4991_VINF_M) < 0.5
    assert abs(bf.arc2.vinf_e - _4991_VINF_E) < 0.5
    assert abs(bf.arc2.vinf_m - _4991_VINF_M) < 0.5
    # ToF half: the g-arc gap CLOSES (circular ~0.14 yr -> ephemeris < 0.1 yr).
    assert res.tof_close
    assert bf.tof_residual_years < 0.1
    assert abs(bf.arc1.arc_tof_years - _4991_G_TOF) < 0.1
    assert abs(bf.arc2.arc_tof_years - _4991_BIGG_TOF) < 0.1
    # Intermediate Earth flyby bend-feasible.
    assert bf.intermediate_flyby_feasible
    assert bf.intermediate_turn_deg < bf.intermediate_max_turn_deg
    # FULL three-way gate: CLOSE (both halves + bend).
    assert res.closed


@pytest.mark.slow
def test_644gg3_two_arc_continuation_secondary() -> None:
    """6.44Gg3 secondary: same machinery on the bigger-gap row (circular g-arc
    0.50 yr off). Asserts the EVIDENCE regime — V_inf stays reachable on the real
    ephemeris and the ToF residual is recorded (whether or not the bigger gap fully
    closes). Pinned as evidence, not a manufactured close."""
    r = free_return_chain_correct(
        _644_APHELION, _644_G_TOF, _644_BIGG_TOF, _644_VINF_E, _644_VINF_M
    )
    a1, e1, a2, e2 = r.arcs[0].a_au, r.arcs[0].e, r.arcs[1].a_au, r.arcs[1].e
    res = continuation_chain_correct(
        a1,
        e1,
        a2,
        e2,
        0.0,
        _644_G_TOF,
        _644_BIGG_TOF,
        _644_VINF_E,
        _644_VINF_M,
        ladder=(1, 3, 9),
    )
    bf = res.best_final
    assert bf is not None
    # V_inf stays reachable on the real ephemeris (the basin survives the ramp).
    assert bf.vinf_residual_kms < 1.0
    # Record the ToF residual regime (evidence, not a manufactured 0).
    assert bf.tof_residual_years >= 0.0
