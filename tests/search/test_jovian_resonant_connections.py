"""Tests for the #754 Jupiter-Europa 3:4-LO homoclinic connection + Anderson & Lo
2011 Table 2 gate (docs/notes/2026-07-28-757-task-b-rescoping-confirmed-families.md
Sec. 5, docs/notes/2026-07-28-754-jupiter-europa-3-4-lo-homoclinic-table2-gate.md),
the #759 Table-3 heteroclinic gate (Wu(3:4-LO) intersect Ws(5:6-LO),
docs/notes/2026-07-28-759-jupiter-europa-3-4-lo-5-6-lo-heteroclinic-table3-gate.md),
and the #766 homoclinic self-connection AT Kumar et al. 2021's own seed energy
(C=3.0041, docs/notes/2026-07-29-766-torus-seed-homoclinic-connection-c30041.md).

Sourced-golden discipline: :data:`jrc.TABLE2_STATE`/:data:`jrc.TABLE3_STATE` trace
verbatim to Anderson & Lo 2011 Tables 2/3 (p.190-191); the gate tolerances
(:data:`jrc.TABLE2_GATE_ABS_TOL`/:data:`jrc.TABLE3_GATE_ABS_TOL`) are justified in
the module docstring (the paper's own states are interpolation-limited, not
Newton-corrected). The gate tests are deliberately NOT marked ``@pytest.mark.slow``
(``feedback_delegation_fresh_agent_not_fork``: a discovery-verdict-bearing evidence
test must not be silently skipped by CI) -- they use the ONE (branch_u, branch_s,
k_u, k_s) combination each task's own coarse scan found (documented in the results
notes), not a re-run of the full scan, to keep runtime bounded while still
exercising the real ``correct_connection``/``crosscheck_cycle``/``gate_table2``/
``gate_table3`` code paths end-to-end.

The #766 tests are NOT gated against a published state (Kumar 2021 reports no
homoclinic connection state for its own seed orbit) -- they instead assert the
SELF-CONSISTENCY evidence directly (Newton residual, ghost-guard margin,
independent Radau cross-check, forward/backward re-approach), per that task's
own honest framing.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    HeteroclinicCycle,
    _planar_floquet_pair,
    _section_crossing,
    _seed_on_manifold,
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


@pytest.fixture(scope="module")
def node_56lo(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    _sys, node = jrc.build_5_6_lo_node(system)
    return node


@pytest.fixture(scope="module")
def node_kumar_c(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    _sys, node, _endpoint = jrc.build_34lo_kumar_c_node(system)
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


# ---------------------------------------------------------------------------
# (5) #759: the Table-3 HETEROCLINIC connection Wu(3:4-LO) -> Ws(5:6-LO).
# ---------------------------------------------------------------------------


def test_table3_state_matches_paper_p191() -> None:
    """Anderson & Lo 2011 Table 3, p.191, verbatim."""
    assert jrc.TABLE3_STATE == {
        "x": -1.43029175,
        "xdot": 0.00018678,
        "ydot": 0.67262261,
    }
    assert jrc.TABLE3_STATE["x"] < 0
    assert jrc.TABLE3_STATE["ydot"] > 0


def test_5_6_lo_node_builds_from_758_candidate(
    system: cr3bp.CR3BPSystem, node_56lo: jrc.ResonantNode
) -> None:
    """build_5_6_lo_node uses #758's reviewer-confirmed candidate, and shares the
    SAME Jacobi constant as 3:4-LO (both evaluated at ANDERSON_LO_C_FLYBY) -- the
    precondition ``correct_connection`` requires for any connection between them.
    """
    assert node_56lo.label == "758-table2-seeded-5:6-LO"
    assert node_56lo.converged
    assert abs(node_56lo.jacobi - jrc.ANDERSON_LO_C_FLYBY) < 1e-9


def test_own_section_points_5_6_lo_is_its_own_ic(node_56lo: jrc.ResonantNode) -> None:
    system = jrc.jupiter_europa_system()
    pts = jrc.own_section_points(system, node_56lo)
    assert len(pts) == 1
    assert abs(float(pts[0][0]) - node_56lo.state0[0]) < 1e-9
    assert abs(float(pts[0][1])) < 1e-9


def test_hetero_ghost_guard_rejects_near_degenerate_accepts_genuine(
    node_34lo: jrc.ResonantNode, node_56lo: jrc.ResonantNode
) -> None:
    """HETERO_GHOST_GUARD_DELTA correctly separates a near-zero-propagation
    degenerate crossing (numerically indistinguishable from an orbit's own IC --
    e.g. the (branch_u=-1, k_u=1) combination this task's own scan found landing
    at x=-1.430408617, only ~6e-7 from 3:4-LO's own IC) from a genuine, converged
    heteroclinic candidate (the (branch=+1,+1, k=4,4) hit this task's scan
    certified, ~0.037-0.127 away from either orbit's own point).
    """
    system = jrc.jupiter_europa_system()
    own34 = jrc.own_section_points(system, node_34lo)
    own56 = jrc.own_section_points(system, node_56lo)

    degenerate = np.array([-1.430408617159821, -1.2679107888463603e-06])
    genuine = np.array([-1.3068706902782006, 0.02869075783271665])

    d_degenerate = min(
        jrc._ghost_distance(degenerate, own34), jrc._ghost_distance(degenerate, own56)
    )
    d_genuine = min(jrc._ghost_distance(genuine, own34), jrc._ghost_distance(genuine, own56))

    assert d_degenerate < jrc.HETERO_GHOST_GUARD_DELTA
    assert d_genuine > jrc.HETERO_GHOST_GUARD_DELTA
    # And landing CLOSE to (but numerically distinct from) 3:4-LO's own point is NOT,
    # by itself, degenerate -- the paper's own text (p.191) says the genuine Table-3
    # intersection sits "near the 3:4 orbit" (#757's own ~1.16e-4 margin finding).
    near_but_genuine = np.array([jrc.TABLE3_STATE["x"], jrc.TABLE3_STATE["xdot"]])
    assert jrc._ghost_distance(near_but_genuine, own34) > jrc.HETERO_GHOST_GUARD_DELTA


# The one (branch_u, branch_s, k_u, k_s) combination this task's own systematic
# 144-combination coarse scan (branch_u, branch_s in {+1,-1}, k_u, k_s in 1..6)
# found to converge Newton-cleanly AND land closest to Table 3's own state -- see
# the results note for the full scan log (5 combinations converged in total; every
# other combination either failed to reach the section, or converged to a point
# considerably farther from Table 3's own state).
_KNOWN_HIT_T3 = {"branch_u": +1, "branch_s": +1, "k_u": 4, "k_s": 4}


def test_find_heteroclinic_returns_known_combo(
    system: cr3bp.CR3BPSystem, node_34lo: jrc.ResonantNode, node_56lo: jrc.ResonantNode
) -> None:
    """find_heteroclinic's own scan/ghost-guard/ranking plumbing, restricted to the
    ONE known-converging combination (for runtime) -- verifies the function returns
    exactly this candidate, ghost-guard-passed, ranked/reported correctly.
    """
    candidates = jrc.find_heteroclinic(
        system,
        node_34lo,
        node_56lo,
        branches=(_KNOWN_HIT_T3["branch_u"],),
        k_range=range(_KNOWN_HIT_T3["k_u"], _KNOWN_HIT_T3["k_u"] + 1),
        max_time_factor=3.0,
        scan_n=8,
        tol=1e-7,
        max_iter=20,
    )
    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.connection.converged
    assert cand.connection.residual < 1e-7
    assert cand.ghost_distance > jrc.HETERO_GHOST_GUARD_DELTA
    assert abs(float(cand.connection.crossing_xv[0]) - (-1.3068706902782006)) < 1e-6
    assert abs(float(cand.connection.crossing_xv[1]) - 0.02869075783271665) < 1e-6


def test_table3_gate_honest_result(
    system: cr3bp.CR3BPSystem, node_34lo: jrc.ResonantNode, node_56lo: jrc.ResonantNode
) -> None:
    """End-to-end Table-3 gate on the best converged heteroclinic candidate this
    task's own systematic scan found.

    Reports the HONEST result: Newton-converged (residual ~1.2e-8), ghost-guard-
    passed (well clear of either orbit's own section point), independently
    Radau-cross-checked -- but its (x, xdot) misses Anderson & Lo's own Table 3
    state by roughly two orders of magnitude MORE than Table 2's own honest miss
    (#754: 1.13e-3 in x). This is reported plainly, not fudged -- see the results
    note for the striking, separately-documented corroborating evidence that
    Wu(3:4-LO) ALONE (without a certified Ws(5:6-LO) match) reproduces Table 3's
    own state to ~4.5e-7, far inside tolerance -- just not as part of a certified
    two-manifold connection.
    """
    conn = correct_connection(
        system,
        node_34lo,
        node_56lo,
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=_KNOWN_HIT_T3["branch_u"],
        branch_s=_KNOWN_HIT_T3["branch_s"],
        k_u=_KNOWN_HIT_T3["k_u"],
        k_s=_KNOWN_HIT_T3["k_s"],
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

    own34 = jrc.own_section_points(system, node_34lo)
    own56 = jrc.own_section_points(system, node_56lo)
    d_ghost = min(
        jrc._ghost_distance(conn.crossing_xv, own34), jrc._ghost_distance(conn.crossing_xv, own56)
    )
    assert d_ghost > jrc.HETERO_GHOST_GUARD_DELTA

    # Independent Radau re-derivation -- a manual one-leg (non-wraparound) "cycle",
    # since Table 3 is a one-way connection, not a closed cycle (assemble_cycle's
    # own wraparound would incorrectly also try to certify the UNSCANNED reverse
    # leg Ws(3:4-LO)<-Wu(5:6-LO), which this task never scanned for).
    cycle = HeteroclinicCycle(
        orbits=[node_34lo.label, node_56lo.label],
        connections=[conn],
        jacobi=node_34lo.jacobi,
        closed=conn.converged,
        max_leg_residual=conn.residual,
        independent_residual=float("nan"),
        symbol_sequence=[node_34lo.label, node_56lo.label],
        notes="",
    )
    checked = crosscheck_cycle(
        system,
        [node_34lo, node_56lo],
        cycle,
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
        max_time_factor=3.0,
    )
    assert checked.independent_residual < 1e-6, (
        f"DOP853 vs Radau disagreement {checked.independent_residual:.3e} exceeds 1e-6"
    )

    candidate = jrc.HeteroclinicCandidate(
        connection=conn,
        ghost_distance=d_ghost,
        dist_to_table3=float(
            np.linalg.norm(
                conn.crossing_xv - np.array([jrc.TABLE3_STATE["x"], jrc.TABLE3_STATE["xdot"]])
            )
        ),
    )
    gate = jrc.gate_table3(system, [candidate])

    # Honest report: both x and xdot miss tolerance badly -- overall gate FAILS.
    assert gate.candidate is not None
    assert gate.x_err > 0.1, f"x_err={gate.x_err:.3e} unexpectedly small for the known result"
    assert not gate.x_passed
    assert not gate.xdot_passed
    assert not gate.passed


def test_gate_table3_reports_clean_fail_on_empty_candidates(
    system: cr3bp.CR3BPSystem,
) -> None:
    """No surviving candidates -> an honest, explicit FAIL, never a fabricated pass."""
    gate = jrc.gate_table3(system, [])
    assert not gate.passed
    assert gate.candidate is None
    assert gate.notes


def test_wu_3_4_lo_alone_reproduces_table3_state_almost_exactly(
    system: cr3bp.CR3BPSystem, node_34lo: jrc.ResonantNode
) -> None:
    """Striking, honestly-caveated corroborating finding (NOT a certified two-
    manifold connection -- see the results note): tracing Wu(3:4-LO)'s OWN section
    curve directly (branch=+1, k=1) -- the SAME interpolation method the paper
    itself used (p.190-191) -- at phase tau=12.74694 lands within ~4.5e-7 of
    Anderson & Lo's own published Table-3 state, an order of magnitude tighter
    than even this module's own 1e-4 gate tolerance. This does NOT by itself
    certify a heteroclinic connection (Ws(5:6-LO) must ALSO be shown to pass
    through this point, which this task's own extensive search could not certify
    to Newton-convergence -- see :func:`test_table3_gate_honest_result` and the
    results note) -- but it is strong, quantitative, independently-computed
    evidence that the genuine intersection sits almost exactly where the paper
    says, and a fourth independent corroborating axis (beyond #755's eigenvalue
    match, #755's shape/close-approach match, and #754's homoclinic Table-2 match)
    for the #755-confirmed 3:4-LO identification.
    """
    tau = 12.74694
    max_time = 3.0 * max(node_34lo.period, 38.76527763438627)
    seed = _seed_on_manifold(
        system, node_34lo, tau=tau, direction="unstable", branch=+1, epsilon=jrc.ANDERSON_LO_EPSILON
    )
    p = _section_crossing(
        system,
        seed,
        direction="unstable",
        k=1,
        max_time=max_time,
        ydot_sign=jrc.SECTION_YDOT_SIGN,
        x_sign=jrc.SECTION_X_SIGN,
    )
    assert p is not None
    target = np.array([jrc.TABLE3_STATE["x"], jrc.TABLE3_STATE["xdot"]])
    x_err = abs(float(p[0]) - target[0])
    xdot_err = abs(float(p[1]) - target[1])
    assert x_err < 1e-5, f"x_err={x_err:.3e}"
    assert xdot_err < 1e-5, f"xdot_err={xdot_err:.3e}"
    # far inside the module's own formal Table-3 gate tolerance
    assert x_err < jrc.TABLE3_GATE_ABS_TOL
    assert xdot_err < jrc.TABLE3_GATE_ABS_TOL


# ---------------------------------------------------------------------------
# (6) #766: homoclinic self-connection AT Kumar et al. 2021's own seed energy
# (C=3.0041) -- self-consistency evidence only, NO published state to gate
# against (see docs/notes/2026-07-29-766-torus-seed-homoclinic-connection-c30041.md).
# ---------------------------------------------------------------------------


def test_build_34lo_kumar_c_node_matches_761_finding(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """The continued 3:4-LO member AT KUMAR_2021_C reproduces #761's own
    recorded numbers (docs/notes/2026-07-29-761-torus-seed-continuation-tractability.md):
    x0=-1.3852484456241585, period=25.31211964876615, |lambda|~54.5898, a
    genuine real saddle -- NOT a re-derivation of a hardcoded seed, a fresh
    continuation gauntlet re-run every call (see build_34lo_kumar_c_node's
    own docstring).
    """
    import cyclerfinder.search.jovian_resonant_families as jrf

    _sys, node, endpoint = jrc.build_34lo_kumar_c_node(system)
    assert abs(node.jacobi - jrf.KUMAR_2021_C) < 1e-9
    assert abs(float(node.state0[0]) - (-1.3852484456241585)) < 1e-6
    assert abs(node.period - 25.31211964876615) < 1e-5
    assert endpoint.is_real_unstable
    assert abs(endpoint.max_eigenvalue - 54.589750588974) < 1e-3
    # And node_kumar_c (the module-scoped fixture) matches, confirming
    # reproducibility across independent calls.
    assert abs(node.jacobi - node_kumar_c.jacobi) < 1e-12
    assert abs(float(node.state0[0]) - float(node_kumar_c.state0[0])) < 1e-9


def test_own_section_points_kumar_c_is_its_own_ic(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """Same structural check as the C_flyby node: the C=3.0041 member's own
    qualifying section point is its own IC."""
    pts = jrc.own_section_points(system, node_kumar_c)
    assert len(pts) == 1
    assert abs(float(pts[0][0]) - node_kumar_c.state0[0]) < 1e-9
    assert abs(float(pts[0][1])) < 1e-9


# The one (branch_u, branch_s, k_u, k_s) combination this task's own coarse scan
# (branch_u, branch_s in {+1,-1}, k_u, k_s in 1..8, max_time_factor=8.0 -- needed
# because the ~19x weaker instability at this energy, |lambda|~54.6 vs the
# C_flyby family's own 1036, makes the manifold grow away from the orbit far
# more slowly than #754's own scan needed) found to converge Newton-cleanly to
# a genuine (ghost-guard-passed) homoclinic self-intersection -- see the
# results note for the full scan log. Structurally analogous to Anderson & Lo's
# own Table-2 self-connection (#754): xdot lands at ~0, the on-symmetry-axis
# case (branch_u == branch_s, k_u == k_s).
_KNOWN_HIT_766 = {"branch_u": +1, "branch_s": +1, "k_u": 6, "k_s": 6}
_KNOWN_HIT_766_TAU = {"tau_u": 10.72913431175392, "tau_s": 14.582984371438714}

# A mirror pair the same scan also found (k_u, k_s swapped, off-symmetry-axis
# xdot != 0) -- additional independent corroboration that genuine transversal
# homoclinic self-intersections exist at this energy (not just an isolated
# fluke at one specific (k_u, k_s)).
_MIRROR_HIT_766_A = {"k_u": 5, "k_s": 6, "tau_u": 6.991901899057191, "tau_s": 19.244249476652467}
_MIRROR_HIT_766_B = {"k_u": 6, "k_s": 5, "tau_u": 6.067870972987947, "tau_s": 18.32021607044972}


def _correct_known_hit_766(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    k_u: int,
    k_s: int,
    tau_u: float,
    tau_s: float,
) -> HeteroclinicConnection:
    return correct_connection(
        system,
        node,
        node,
        k_u=k_u,
        k_s=k_s,
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=+1,
        branch_s=+1,
        tau_u0=tau_u,
        tau_s0=tau_s,
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
        max_time_factor=8.0,
        tol=1e-9,
        max_iter=40,
    )


def test_kumar_c_homoclinic_known_hit_converges_with_real_ghost_margin(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """The (branch_u=+1, branch_s=+1, k_u=6, k_s=6) combination converges to a
    Newton residual far tighter than this chain's own demonstrated corrector
    precision elsewhere (<=1e-9), and the crossing sits well clear of the
    orbit's own trivial section point -- a REAL margin (>10x the
    GHOST_GUARD_DELTA=1e-3 threshold), not a razor-thin one.
    """
    conn = _correct_known_hit_766(
        system, node_kumar_c, _KNOWN_HIT_766["k_u"], _KNOWN_HIT_766["k_s"], **_KNOWN_HIT_766_TAU
    )
    assert conn.converged, f"expected the known-hit combo to converge; notes={conn.notes}"
    assert conn.residual < 1e-9, f"residual={conn.residual:.3e}"

    own_pts = jrc.own_section_points(system, node_kumar_c)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    assert d_ghost > jrc.GHOST_GUARD_DELTA, (
        f"converged crossing at {conn.crossing_xv} is within the ghost-guard radius "
        f"of the orbit's own section point(s) {own_pts} (d={d_ghost:.3e})"
    )
    # A REAL margin -- an order of magnitude past the guard threshold, not a
    # borderline pass (feedback_verify_automated_ghost_guard_booleans).
    assert d_ghost > 10 * jrc.GHOST_GUARD_DELTA, f"d_ghost={d_ghost:.3e} -- margin too thin"
    # Structurally the same on-symmetry-axis type as #754's own Table-2 point:
    # xdot lands at ~0.
    assert abs(float(conn.crossing_xv[1])) < 1e-6


def test_kumar_c_homoclinic_independent_radau_crosscheck(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """Independent Radau re-derivation of the known #766 hit (mandatory before
    trusting the crossing, per this chain's own orbit-closure discipline).
    """
    conn = _correct_known_hit_766(
        system, node_kumar_c, _KNOWN_HIT_766["k_u"], _KNOWN_HIT_766["k_s"], **_KNOWN_HIT_766_TAU
    )
    assert conn.converged

    cycle = assemble_cycle(
        system,
        [node_kumar_c],
        tol=1e-9,
        connection_kwargs={
            "epsilon": jrc.ANDERSON_LO_EPSILON,
            "branch_u": _KNOWN_HIT_766["branch_u"],
            "branch_s": _KNOWN_HIT_766["branch_s"],
            "k_u": _KNOWN_HIT_766["k_u"],
            "k_s": _KNOWN_HIT_766["k_s"],
            "tau_u0": conn.tau_u,
            "tau_s0": conn.tau_s,
            "ydot_sign_u": jrc.SECTION_YDOT_SIGN,
            "ydot_sign_s": jrc.SECTION_YDOT_SIGN,
            "x_sign_u": jrc.SECTION_X_SIGN,
            "x_sign_s": jrc.SECTION_X_SIGN,
            "max_time_factor": 8.0,
            "max_iter": 40,
        },
    )
    assert cycle.closed
    checked = crosscheck_cycle(
        system, [node_kumar_c], cycle, method="Radau", rtol=1e-11, atol=1e-11, max_time_factor=8.0
    )
    assert checked.independent_residual < 1e-6, (
        f"DOP853 vs Radau disagreement {checked.independent_residual:.3e} exceeds 1e-6"
    )


def test_kumar_c_homoclinic_forward_backward_reapproach(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """Forward/backward re-approach self-consistency (#766's own honest gate,
    since no published state exists to reproduce here): the found intersection
    state, propagated BACKWARD by the unstable leg's own elapsed transit time,
    reproduces the epsilon-scale unstable-manifold seed it was found from --
    genuine homoclinic behaviour, not a numerical artifact. Symmetric forward
    check on the stable leg.
    """
    conn = _correct_known_hit_766(
        system, node_kumar_c, _KNOWN_HIT_766["k_u"], _KNOWN_HIT_766["k_s"], **_KNOWN_HIT_766_TAU
    )
    assert conn.converged
    own_pts = jrc.own_section_points(system, node_kumar_c)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    candidate = jrc.HomoclinicCandidate(
        connection=conn, ghost_distance=d_ghost, dist_to_table2=float("nan")
    )
    reap = jrc.homoclinic_reapproach_check(system, node_kumar_c, candidate, max_time_factor=8.0)
    # Both legs took ~5.3 orbital periods to develop (weak instability at this
    # energy) -- reported plainly, not assumed to be a fixed "3 periods" the
    # way #754's own (much more strongly unstable) case used.
    assert reap.t_u > 5.0 * node_kumar_c.period
    assert abs(reap.t_s) > 5.0 * node_kumar_c.period
    assert reap.backward_distance < 1e-5, f"backward_distance={reap.backward_distance:.3e}"
    assert reap.forward_distance < 1e-3, f"forward_distance={reap.forward_distance:.3e}"


def test_kumar_c_homoclinic_mirror_pair_also_converges(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """A second, independent (k_u, k_s)-swapped mirror pair ((5,6) and (6,5))
    ALSO converges to a genuine, ghost-guard-passed homoclinic point --
    corroborating that C=3.0041 genuinely supports transversal homoclinic
    self-intersections, not just a single isolated fluke at (6,6).
    """
    own_pts = jrc.own_section_points(system, node_kumar_c)
    conns = []
    for hit in (_MIRROR_HIT_766_A, _MIRROR_HIT_766_B):
        conn = _correct_known_hit_766(
            system,
            node_kumar_c,
            int(hit["k_u"]),
            int(hit["k_s"]),
            float(hit["tau_u"]),
            float(hit["tau_s"]),
        )
        assert conn.converged, f"expected {hit} to converge; notes={conn.notes}"
        assert conn.residual < 1e-8
        d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
        assert d_ghost > 10 * jrc.GHOST_GUARD_DELTA
        conns.append(conn)
    # The two mirror hits land at (x, +-xdot) of each other (reflection
    # symmetry), not the same point twice.
    conn_a, conn_b = conns
    assert abs(float(conn_a.crossing_xv[0]) - float(conn_b.crossing_xv[0])) < 1e-6
    assert abs(float(conn_a.crossing_xv[1]) + float(conn_b.crossing_xv[1])) < 1e-6


def test_find_homoclinic_rank_by_residual_reports_nan_distance(
    system: cr3bp.CR3BPSystem, node_kumar_c: jrc.ResonantNode
) -> None:
    """find_homoclinic's #766 extension: rank_by_residual=True (the honest
    ranking when no published target exists at this energy) restricts the
    scan to the known combination (for runtime) and reports dist_to_table2 as
    nan rather than a distance to an irrelevant Table-2 target.
    """
    candidates = jrc.find_homoclinic(
        system,
        node_kumar_c,
        branches=(_KNOWN_HIT_766["branch_u"],),
        k_range=range(_KNOWN_HIT_766["k_u"], _KNOWN_HIT_766["k_u"] + 1),
        max_time_factor=8.0,
        scan_n=8,
        max_iter=20,
        tol=1e-7,
        rank_by_residual=True,
    )
    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.connection.converged
    assert np.isnan(cand.dist_to_table2)
    assert cand.ghost_distance > jrc.GHOST_GUARD_DELTA


def test_find_homoclinic_default_target_unchanged(system: cr3bp.CR3BPSystem) -> None:
    """The #766 target/rank_by_residual extension must not change ANY
    pre-existing caller's behaviour: default (target=None, rank_by_residual=False)
    still ranks by distance to TABLE2_STATE, exactly as before #766.
    """
    _sys, node = jrc.build_3_4_lo_node(system)
    candidates = jrc.find_homoclinic(
        system,
        node,
        branches=(+1, -1),
        k_range=range(3, 4),
        max_time_factor=3.0,
        scan_n=8,
        max_iter=20,
        tol=1e-7,
    )
    assert len(candidates) == 1
    assert not np.isnan(candidates[0].dist_to_table2)
    expected = float(
        np.linalg.norm(
            candidates[0].connection.crossing_xv
            - np.array([jrc.TABLE2_STATE["x"], jrc.TABLE2_STATE["xdot"]])
        )
    )
    assert abs(candidates[0].dist_to_table2 - expected) < 1e-12
