"""Tests for the #754 Jupiter-Europa 3:4-LO homoclinic connection + Anderson & Lo
2011 Table 2 gate (docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md
Sec. 5, docs/notes/2026-07-28-754-jupiter-europa-3-4-lo-homoclinic-table2-gate.md).

Sourced-golden discipline: :data:`jrc.TABLE2_STATE` traces verbatim to Anderson & Lo
2011 Table 2 (p.190); the gate tolerance (:data:`jrc.TABLE2_GATE_ABS_TOL`) is
justified in the module docstring (the paper's own state is interpolation-limited,
not Newton-corrected). The Table-2 gate test is deliberately NOT marked
``@pytest.mark.slow`` (``feedback_delegation_fresh_agent_not_fork``: a discovery-
verdict-bearing evidence test must not be silently skipped by CI) -- it uses the
ONE (branch_u, branch_s, k_u, k_s) combination this task's own coarse scan found to
converge (documented in the results note), not a re-run of the full scan, to keep
runtime bounded while still exercising the real ``correct_connection``/
``crosscheck_cycle``/``gate_table2`` code paths end-to-end.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
from cyclerfinder.genome.heteroclinic_cycle import (
    _planar_floquet_pair,
    assemble_cycle,
    correct_connection,
    crosscheck_cycle,
)
from cyclerfinder.search.jovian_resonant_families import ResonantFamilyCandidate


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return jrc.jupiter_europa_system()


@pytest.fixture(scope="module")
def node_34lo(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    _sys, node = jrc.build_3_4_lo_node(system)
    return node


# ---------------------------------------------------------------------------
# (1) Sourced-constant reproduce-before-trust gate.
# ---------------------------------------------------------------------------


def test_table2_state_matches_paper_p190() -> None:
    """Anderson & Lo 2011 Table 2, p.190, verbatim."""
    assert jrc.TABLE2_STATE == {
        "x": -1.28427733,
        "xdot": 0.00000009,
        "ydot": 0.46372205,
    }


def test_section_convention_matches_paper_p170_171() -> None:
    """Eq. 7 / p.170-171: one-sided {y=0}, x<0, ydot>0 -- both Table 2/3 rows agree."""
    assert jrc.SECTION_X_SIGN == -1
    assert jrc.SECTION_YDOT_SIGN == +1
    assert jrc.TABLE2_STATE["x"] < 0
    assert jrc.TABLE2_STATE["ydot"] > 0


def test_ydot_from_section_eq7_matches_jacobi_constant(system: cr3bp.CR3BPSystem) -> None:
    """Eq. 7's ydot recovery must agree with cr3bp.jacobi_constant's own formula.

    Round-trip: take a state's own (x, xdot, ydot), recompute ydot via Eq. 7 from
    (x, xdot, jacobi_constant(state)) -- must reproduce the original ydot (the two
    formulas are algebraically the same equation solved for the same variable).
    """
    state = np.array([-1.3, 0.0, 0.0, 0.01, 0.5, 0.0])
    c = cr3bp.jacobi_constant(state, system.mu)
    ydot_rec = jrc.ydot_from_section_eq7(system, -1.3, 0.01, c)
    assert abs(ydot_rec - 0.5) < 1e-12


def test_ydot_from_section_eq7_raises_on_unreachable_point(system: cr3bp.CR3BPSystem) -> None:
    """A (x, xdot, C) triple with negative radicand is not a physical section point."""
    with pytest.raises(ValueError, match="negative radicand"):
        jrc.ydot_from_section_eq7(system, -1.3, 0.01, jacobi=100.0)


# ---------------------------------------------------------------------------
# (2) ResonantNode adapter: eigenvector recomputation matches direct Floquet analysis.
# ---------------------------------------------------------------------------


def test_resonant_node_matches_direct_floquet_pair(system: cr3bp.CR3BPSystem) -> None:
    """ResonantNode.from_candidate's saddle pair equals a direct
    _planar_floquet_pair call on the same (state0, period) -- no hidden
    transformation, just the existing machinery applied to the candidate's stored
    IC, exactly as the #757 scoping note specifies.
    """
    from cyclerfinder.search.jovian_resonant_families import recover_table1_candidate

    cand = recover_table1_candidate("3:4-LO", system)
    node = jrc.ResonantNode.from_candidate(system, cand)

    state0_direct = np.array([cand.x0, 0.0, 0.0, 0.0, cand.ydot0, 0.0])
    lam_u, v_u, _lam_s, v_s = _planar_floquet_pair(system, state0_direct, cand.period)

    assert node.label == "3:4-LO"
    assert np.array_equal(node.state0, state0_direct)
    assert node.period == cand.period
    assert node.jacobi == cand.jacobi
    assert np.array_equal(node.unstable_eigvec, v_u)
    assert np.array_equal(node.stable_eigvec, v_s)
    assert node.converged
    # And the recomputed unstable eigenvalue agrees with the candidate's own
    # Barden-classified max_eigenvalue (the #755-confirmed 1036.116... family).
    rel_err = abs(lam_u - cand.max_eigenvalue) / cand.max_eigenvalue
    assert rel_err < 1e-6, f"rel_err={rel_err:.3e}"


def test_resonant_node_rejects_non_real_unstable_candidate(system: cr3bp.CR3BPSystem) -> None:
    """A candidate whose nontrivial eigenpair is complex (marginally stable, e.g. a
    5:6-LI-like near-unit-modulus pair) is NOT a valid saddle for the
    _planar_floquet_pair magnitude-only convention -- must raise, not silently
    build a bogus node.
    """
    fake = ResonantFamilyCandidate(
        label="fake-complex",
        x0=-1.4304078294961569,
        ydot0=0.6727850993335511,
        period=25.672528919046933,
        jacobi=2.99163956830415,
        ydot0_sign=1.0,
        half_crossings=6,
        crossing_residual=1e-13,
        barden_eigenvalue=complex(0.963, 0.270),
        max_eigenvalue=1.0,
        is_real_unstable=False,
        planar_floquet_eigenvalue=1.0,
        period_over_2pi=4.0859,
    )
    with pytest.raises(ValueError, match="is_real_unstable"):
        jrc.ResonantNode.from_candidate(system, fake)


def test_resonant_node_rejects_stale_eigenvalue_mismatch(system: cr3bp.CR3BPSystem) -> None:
    """A candidate whose stored max_eigenvalue disagrees with the state/period's OWN
    recomputed Floquet eigenvalue is stale/mismatched data -- must raise loudly.
    """
    fake = ResonantFamilyCandidate(
        label="fake-mismatched",
        x0=-1.4304078294961569,  # genuine 3:4-LO IC ...
        ydot0=0.6727850993335511,
        period=25.672528919046933,
        jacobi=2.99163956830415,
        ydot0_sign=1.0,
        half_crossings=6,
        crossing_residual=1e-13,
        barden_eigenvalue=complex(5.0, 0.0),
        max_eigenvalue=5.0,  # ... but a WRONG max_eigenvalue (true is ~1036.116)
        is_real_unstable=True,
        planar_floquet_eigenvalue=5.0,
        period_over_2pi=4.0859,
    )
    with pytest.raises(ValueError, match="disagrees"):
        jrc.ResonantNode.from_candidate(system, fake)


# ---------------------------------------------------------------------------
# (3) Homoclinic trivial-solution ghost guard.
# ---------------------------------------------------------------------------


def test_ghost_distance_rejects_own_point_accepts_genuine_intersection() -> None:
    """_ghost_distance correctly separates a trivial self-shadow from a genuine,
    well-separated homoclinic point (using the module's own GHOST_GUARD_DELTA).
    """
    own_pts = [np.array([-1.4304078294961569, 0.0])]
    trivial = np.array([-1.4304078294961569 + 1e-7, 5e-8])  # numerically == own point
    genuine = np.array([jrc.TABLE2_STATE["x"], jrc.TABLE2_STATE["xdot"]])  # real Table 2 point

    assert jrc._ghost_distance(trivial, own_pts) < jrc.GHOST_GUARD_DELTA
    assert jrc._ghost_distance(genuine, own_pts) > jrc.GHOST_GUARD_DELTA


def test_own_section_points_is_the_orbits_own_ic(node_34lo: jrc.ResonantNode) -> None:
    """3:4-LO's own qualifying section point is its own IC (x0=-1.4304...), and the
    Table-2 homoclinic point sits well outside the ghost-guard radius of it -- the
    ~0.146 margin the #757 scoping note documents.
    """
    system = jrc.jupiter_europa_system()
    pts = jrc.own_section_points(system, node_34lo)
    assert len(pts) == 1
    assert abs(float(pts[0][0]) - node_34lo.state0[0]) < 1e-9
    assert abs(float(pts[0][1])) < 1e-9  # xdot=0 at this perpendicular IC

    genuine = np.array([jrc.TABLE2_STATE["x"], jrc.TABLE2_STATE["xdot"]])
    d = jrc._ghost_distance(genuine, pts)
    assert d > 0.14, f"expected the ~0.146 margin documented in #757; got {d:.4f}"
    assert d > jrc.GHOST_GUARD_DELTA


# ---------------------------------------------------------------------------
# (4) The Table-2 gate itself (NOT @pytest.mark.slow -- must run in CI).
# ---------------------------------------------------------------------------

# The one (branch_u, branch_s, k_u, k_s) combination found by this task's own coarse
# scan over branch_u, branch_s in {+1,-1} and k_u, k_s in 1..6 to converge to a
# genuine (non-ghost) homoclinic self-intersection -- see the results note for the
# full scan log (every other combination tried either failed to converge or landed
# on the orbit's own trivial section point).
_KNOWN_HIT = {"branch_u": +1, "branch_s": -1, "k_u": 3, "k_s": 3}


def test_table2_gate_honest_result(system: cr3bp.CR3BPSystem, node_34lo: jrc.ResonantNode) -> None:
    """End-to-end Table-2 gate on the one converged homoclinic self-connection.

    Reports the HONEST result: Newton-converged, ghost-guard-passed, independently
    Radau-cross-checked -- but the (x, xdot) match against Anderson & Lo's own
    Table 2 state MISSES the 1e-4 gate on x (by ~1.1e-3, an order of magnitude),
    while xdot matches to <1e-4. This is a genuine, reported FAIL -- not fudged or
    loosened (see results note for the full discussion of why this is still a
    meaningful positive signal for the #755 3:4-LO identification).
    """
    conn = correct_connection(
        system,
        node_34lo,
        node_34lo,
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=_KNOWN_HIT["branch_u"],
        branch_s=_KNOWN_HIT["branch_s"],
        k_u=_KNOWN_HIT["k_u"],
        k_s=_KNOWN_HIT["k_s"],
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
        max_time_factor=3.0,
        scan_n=8,
        tol=1e-7,
        max_iter=20,
    )
    assert conn.converged, f"expected the known-hit combo to converge; notes={conn.notes}"
    assert conn.residual < 1e-7

    # Ghost guard: this must NOT be the orbit trivially shadowing itself.
    own_pts = jrc.own_section_points(system, node_34lo)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    assert d_ghost > jrc.GHOST_GUARD_DELTA, (
        f"converged crossing at {conn.crossing_xv} is within the ghost-guard radius "
        f"of the orbit's own section point(s) {own_pts} (d={d_ghost:.3e}) -- "
        "this would be a trivial non-connection, not a genuine homoclinic point"
    )

    # Independent Radau re-derivation (mandatory before trusting the crossing).
    cycle = assemble_cycle(
        system,
        [node_34lo],
        tol=1e-7,
        connection_kwargs={
            "epsilon": jrc.ANDERSON_LO_EPSILON,
            "branch_u": _KNOWN_HIT["branch_u"],
            "branch_s": _KNOWN_HIT["branch_s"],
            "k_u": _KNOWN_HIT["k_u"],
            "k_s": _KNOWN_HIT["k_s"],
            "ydot_sign_u": jrc.SECTION_YDOT_SIGN,
            "ydot_sign_s": jrc.SECTION_YDOT_SIGN,
            "x_sign_u": jrc.SECTION_X_SIGN,
            "x_sign_s": jrc.SECTION_X_SIGN,
            "max_time_factor": 3.0,
            "scan_n": 8,
            "max_iter": 20,
        },
    )
    assert cycle.closed
    checked = crosscheck_cycle(
        system, [node_34lo], cycle, method="Radau", rtol=1e-11, atol=1e-11, max_time_factor=3.0
    )
    assert checked.independent_residual < 1e-6, (
        f"DOP853 vs Radau disagreement {checked.independent_residual:.3e} exceeds 1e-6"
    )

    candidate = jrc.HomoclinicCandidate(
        connection=conn,
        ghost_distance=d_ghost,
        dist_to_table2=float(
            np.linalg.norm(
                conn.crossing_xv - np.array([jrc.TABLE2_STATE["x"], jrc.TABLE2_STATE["xdot"]])
            )
        ),
    )
    gate = jrc.gate_table2(system, [candidate])

    # Honest report: x misses tolerance, xdot passes, overall gate FAILS.
    assert gate.candidate is not None
    assert gate.x_err < 2e-3, f"x_err={gate.x_err:.3e} unexpectedly far from the known result"
    assert gate.xdot_passed, f"xdot_err={gate.xdot_err:.3e} should be well within tolerance"
    assert not gate.x_passed, (
        f"x_err={gate.x_err:.3e} unexpectedly passed 1e-4 -- if this regresses to a real "
        "pass, update the results note (this would be genuinely good news, not a bug)"
    )
    assert not gate.passed


def test_gate_table2_reports_clean_fail_on_empty_candidates(
    system: cr3bp.CR3BPSystem,
) -> None:
    """No surviving candidates -> an honest, explicit FAIL, never a fabricated pass."""
    gate = jrc.gate_table2(system, [])
    assert not gate.passed
    assert gate.candidate is None
    assert gate.notes
