"""Tests for `#781`: Neptune-Triton 4:5 resonant-orbit homoclinic self-connection
(docs/notes/2026-08-08-781-neptune-triton-homoclinic-connections.md).

Sourced-golden discipline: this task has NO published state to gate against --
Miceli & Bosanac 2026 publishes no homoclinic/heteroclinic connection content
of any kind for this system (see the module docstring's "LITERATURE NOVELTY
GATE" section). Every test below therefore asserts SELF-CONSISTENCY evidence
directly (Newton residual, ghost-guard margin, independent Radau cross-check,
forward/backward re-approach), per this task chain's own honest framing
(`#766`'s/`#767`'s own precedent). None are marked ``@pytest.mark.slow`` -- a
discovery-verdict-bearing evidence test must run in CI, not be silently
skipped (``feedback_delegation_fresh_agent_not_fork``).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.neptune_triton_resonant_connections as ntc
import cyclerfinder.search.neptune_triton_resonant_families as ntf
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    assemble_cycle,
    correct_connection,
    crosscheck_cycle,
)


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return ntf.neptune_triton_system()


@pytest.fixture(scope="module")
def node(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    _sys, nd, _row = ntc.build_45_node(system)
    return nd


# ---------------------------------------------------------------------------
# (1) build_45_node reuses #776's own confirmed "4:5-saddle" candidate.
# ---------------------------------------------------------------------------


def test_build_45_node_matches_776_confirmed_row(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """The node's own state/period/jacobi match #776's own recorded 4:5-saddle
    numbers (docs/notes/2026-08-01-776-neptune-triton-resonant-families-gate.md):
    x0=-0.969056422015..., period=30.398418802764855, jacobi=2.987089791658,
    max_eigenvalue=-105.05423683377265 (the ONLY negative eigenvalue among
    #776's ten gate rows -- see module docstring).
    """
    _sys, _nd, row = ntc.build_45_node(system)
    assert row.passed
    assert row.label == "4:5-saddle"
    assert row.is_real_unstable

    assert abs(float(node.state0[0]) - (-0.9690564220154437)) < 1e-9
    assert abs(node.period - 30.398418802764855) < 1e-8
    assert abs(node.jacobi - 2.987089791658) < 1e-12
    assert node.converged


def test_resonant_node_is_real_saddle_despite_negative_stored_eigenvalue(
    node: jrc.ResonantNode,
) -> None:
    """`_planar_floquet_pair`'s own magnitude-only convention still produces a
    unit-normalised eigenvector direction even though `#776`'s own stored
    `max_eigenvalue` is genuinely NEGATIVE (-105.05...) -- see module
    docstring ("NEGATIVE-REAL SADDLE EIGENVALUE").
    """
    lam_norm = float(np.linalg.norm(node.unstable_eigvec))
    assert abs(lam_norm - 1.0) < 1e-9


def test_node_from_candidate_staleness_guard_fixed_for_negative_eigenvalue(
    system: cr3bp.CR3BPSystem,
) -> None:
    """`#781` fix regression: `ResonantNode.from_candidate` must actually
    build successfully (not silently mis-validate) for a candidate whose
    stored `max_eigenvalue` is negative -- see `jovian_resonant_connections.py`'s
    own updated docstring for the before/after `rel_err` formula account.
    """
    cand = ntf.recover_esm_candidate("4:5-saddle", system)
    assert cand.max_eigenvalue < 0.0
    nd = jrc.ResonantNode.from_candidate(system, cand)
    assert nd.converged


# ---------------------------------------------------------------------------
# (2) Section convention: BOTH x_sign and ydot_sign fully unrestricted --
# generalizing even further than Saturn-Titan's own x-only relaxation (see
# module docstring, "SECTION CONVENTION").
# ---------------------------------------------------------------------------


def test_own_section_points_unions_all_four_sign_combinations(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """The 4:5-saddle orbit's own qualifying section points (both signs
    unrestricted) recover the full 6-point per-period crossing set,
    confirmed directly this task by inspecting the orbit's own {y=0}
    crossing list -- no double-counting across the four sign buckets.
    """
    pts = ntc.own_section_points(system, node)
    assert len(pts) == 6


def test_own_section_points_includes_both_the_ic_and_half_period_return(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """The two PERPENDICULAR (xdot=0) reference points -- the orbit's own IC
    (x0=-0.969056..., ydot<0) and the half-period return (x=0.981105...,
    ydot>0, `half_crossings=3`) -- sit at OPPOSITE `ydot` signs. This is the
    regression test that would have caught the ghost-guard trap the module
    docstring documents: a single fixed `ydot_sign` (Saturn-Titan's own
    remaining restriction) silently drops ONE of these two, and for this
    orbit specifically, either choice drops a different one of them.
    """
    pts = ntc.own_section_points(system, node)
    perpendicular = [p for p in pts if abs(float(p[1])) < 1e-6]
    assert len(perpendicular) == 2
    xs = sorted(float(p[0]) for p in perpendicular)
    assert abs(xs[0] - (-0.969056422015)) < 1e-6  # the IC
    assert abs(xs[1] - 0.981105) < 1e-4  # the half-period return


def test_single_ydot_sign_choice_would_drop_one_perpendicular_point(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """Sanity check for the module docstring's own claim: EITHER single
    `ydot_sign` choice (the natural generalization Saturn-Titan's own
    x-only relaxation would suggest trying first) drops a different one of
    the two perpendicular reference points -- neither is safe alone, which
    is why this module unions all four `(x_sign, ydot_sign)` combinations
    instead of just two.
    """
    only_neg_ydot = jrc.own_section_points(
        system, node, ydot_sign=-1, x_sign=-1
    ) + jrc.own_section_points(system, node, ydot_sign=-1, x_sign=1)
    only_pos_ydot = jrc.own_section_points(
        system, node, ydot_sign=1, x_sign=-1
    ) + jrc.own_section_points(system, node, ydot_sign=1, x_sign=1)
    neg_has_ic = any(abs(float(p[0]) - (-0.969056422015)) < 1e-6 for p in only_neg_ydot)
    neg_has_half = any(abs(float(p[0]) - 0.981105) < 1e-4 for p in only_neg_ydot)
    pos_has_ic = any(abs(float(p[0]) - (-0.969056422015)) < 1e-6 for p in only_pos_ydot)
    pos_has_half = any(abs(float(p[0]) - 0.981105) < 1e-4 for p in only_pos_ydot)
    assert neg_has_ic and not neg_has_half
    assert pos_has_half and not pos_has_ic


# ---------------------------------------------------------------------------
# (3) The four genuine, independently-converged, ghost-guard-passed hits this
# task's own targeted geometric-seed scan found (all opposite-branch-sign,
# `branch_u=+1, branch_s=-1` -- see module docstring for why the negative
# eigenvalue's own period-to-period branch flip makes this the natural
# convergent combination; same-sign seeds were explored and did NOT converge
# cleanly, an honest negative, not suppressed). All four sit at 134x-287x the
# GHOST_GUARD_DELTA=1e-3 threshold -- a real, non-delicate margin.
# ---------------------------------------------------------------------------

# On-symmetry-axis (xdot~0) -- structurally the SAME point-type as Anderson &
# Lo's own Table 2 state and #766's/#767's own primary hits. Tightest
# evidence bundle of the four -- reported as this task's PRIMARY result.
_PRIMARY = {"branch_u": +1, "branch_s": -1, "k_u": 13, "k_s": 13}
_PRIMARY_TAU = {"tau_u": 1.9994947221216297, "tau_s": 28.39892348910753}

# A SECOND, independently-found on-axis hit at a different k -- corroborates
# that the axis-crossing family of intersections is not a one-off.
_SECONDARY = {"branch_u": +1, "branch_s": -1, "k_u": 15, "k_s": 15}
_SECONDARY_TAU = {"tau_u": 18.3007470738707, "tau_s": 12.097671673053219}

# A genuine mirror pair (k_u, k_s swapped, off-symmetry-axis xdot != 0,
# reflection-symmetric (x, +-xdot)) -- independent corroboration that this
# energy supports transversal homoclinic self-intersections beyond one
# isolated (branch, k) choice.
_MIRROR_A = {"branch_u": +1, "branch_s": -1, "k_u": 13, "k_s": 17}
_MIRROR_A_TAU = {"tau_u": 18.300747081815683, "tau_s": 12.097671665527479}
_MIRROR_B = {"branch_u": +1, "branch_s": -1, "k_u": 17, "k_s": 13}
_MIRROR_B_TAU = {"tau_u": 18.30074707058968, "tau_s": 12.09767166644772}


def _correct_known_hit(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    hit: dict[str, int],
    tau: dict[str, float],
) -> HeteroclinicConnection:
    return correct_connection(
        system,
        node,
        node,
        k_u=hit["k_u"],
        k_s=hit["k_s"],
        epsilon=ntc.EPSILON,
        branch_u=hit["branch_u"],
        branch_s=hit["branch_s"],
        tau_u0=tau["tau_u"],
        tau_s0=tau["tau_s"],
        ydot_sign_u=None,
        ydot_sign_s=None,
        x_sign_u=None,
        x_sign_s=None,
        max_time_factor=6.0,
        tol=1e-9,
        max_iter=60,
        fd_step=1e-7,
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau", "expect_x", "expect_xdot", "min_ghost_margin"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU, 1.16933872, 0.0, 130.0),
        ("secondary", _SECONDARY, _SECONDARY_TAU, -1.38561105, 0.0, 280.0),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU, -1.27434376, 0.12236316, 130.0),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU, -1.27434376, -0.12236316, 130.0),
    ],
)
def test_known_hit_converges_with_real_ghost_margin(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    name: str,
    hit: dict[str, int],
    tau: dict[str, float],
    expect_x: float,
    expect_xdot: float,
    min_ghost_margin: float,
) -> None:
    conn = _correct_known_hit(system, node, hit, tau)
    assert conn.converged, f"{name}: expected convergence; notes={conn.notes}"
    assert conn.residual < 1e-8, f"{name}: residual={conn.residual:.3e}"
    assert abs(float(conn.crossing_xv[0]) - expect_x) < 1e-5
    assert abs(float(conn.crossing_xv[1]) - expect_xdot) < 1e-5

    own_pts = ntc.own_section_points(system, node)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    assert d_ghost > ntc.GHOST_GUARD_DELTA
    # A REAL margin -- two orders of magnitude past the guard threshold, not
    # a borderline pass (feedback_verify_automated_ghost_guard_booleans).
    assert d_ghost > min_ghost_margin * ntc.GHOST_GUARD_DELTA, (
        f"{name}: d_ghost={d_ghost:.3e} -- margin too thin"
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU),
        ("secondary", _SECONDARY, _SECONDARY_TAU),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU),
    ],
)
def test_known_hit_independent_radau_crosscheck(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    name: str,
    hit: dict[str, int],
    tau: dict[str, float],
) -> None:
    """Independent Radau re-derivation of each known hit (mandatory before
    trusting the crossing, per this chain's own orbit-closure discipline).
    """
    conn = _correct_known_hit(system, node, hit, tau)
    assert conn.converged

    cycle = assemble_cycle(
        system,
        [node],
        tol=1e-9,
        connection_kwargs={
            "epsilon": ntc.EPSILON,
            "branch_u": hit["branch_u"],
            "branch_s": hit["branch_s"],
            "k_u": hit["k_u"],
            "k_s": hit["k_s"],
            "tau_u0": conn.tau_u,
            "tau_s0": conn.tau_s,
            "ydot_sign_u": None,
            "ydot_sign_s": None,
            "x_sign_u": None,
            "x_sign_s": None,
            "max_time_factor": 6.0,
            "max_iter": 60,
            "fd_step": 1e-7,
        },
    )
    assert cycle.closed
    checked = crosscheck_cycle(
        system, [node], cycle, method="Radau", rtol=1e-11, atol=1e-11, max_time_factor=6.0
    )
    assert checked.independent_residual < 1e-6, (
        f"{name}: DOP853 vs Radau disagreement {checked.independent_residual:.3e} exceeds 1e-6"
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau", "max_forward"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU, 0.2),
        ("secondary", _SECONDARY, _SECONDARY_TAU, 0.01),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU, 0.01),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU, 0.01),
    ],
)
def test_known_hit_forward_backward_reapproach(
    system: cr3bp.CR3BPSystem,
    node: jrc.ResonantNode,
    name: str,
    hit: dict[str, int],
    tau: dict[str, float],
    max_forward: float,
) -> None:
    """Forward/backward re-approach self-consistency (this task's own honest
    gate, since no published state exists to reproduce here). ``backward_distance``
    is consistently tight (close to numerically-exact time-reversal of the
    same integration) for all four; ``forward_distance`` is looser -- a
    genuine cross-leg consistency test amplified by this orbit's own
    instability. PRIMARY's own ``forward_distance`` (~0.16) is notably
    looser than the other three (~1.7e-3) despite this orbit's own moderate
    instability (|lambda|~105, far milder than Saturn-Titan's 2129.8) --
    reported honestly, not fudged: still far below the O(1-2) trajectory
    scale.
    """
    conn = _correct_known_hit(system, node, hit, tau)
    assert conn.converged
    own_pts = ntc.own_section_points(system, node)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    candidate = ntc.HomoclinicCandidate(connection=conn, ghost_distance=d_ghost)
    reap = ntc.homoclinic_reapproach_check(system, node, candidate, max_time_factor=10.0)

    assert reap.backward_distance < 1e-7, f"{name}: backward_distance={reap.backward_distance:.3e}"
    assert reap.forward_distance < max_forward, (
        f"{name}: forward_distance={reap.forward_distance:.3e}"
    )
    # Both are far smaller than the O(1-2) trajectory scale (positions of
    # order 1-2 nondim units) -- never a coincidental close pass.
    assert reap.forward_distance < 1.0


def test_mirror_pair_reflection_symmetry(system: cr3bp.CR3BPSystem, node: jrc.ResonantNode) -> None:
    """MIRROR_A/MIRROR_B land at (x, +-xdot) of each other (the CR3BP's own
    time-reversal reflection symmetry) -- genuinely two independent
    intersections related by symmetry, not the same one found twice.
    """
    conn_a = _correct_known_hit(system, node, _MIRROR_A, _MIRROR_A_TAU)
    conn_b = _correct_known_hit(system, node, _MIRROR_B, _MIRROR_B_TAU)
    assert conn_a.converged and conn_b.converged
    assert abs(float(conn_a.crossing_xv[0]) - float(conn_b.crossing_xv[0])) < 1e-6
    assert abs(float(conn_a.crossing_xv[1]) + float(conn_b.crossing_xv[1])) < 1e-6


def test_primary_and_secondary_hits_are_on_symmetry_axis(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """PRIMARY (k_u=k_s=13) and SECONDARY (k_u=k_s=15) both sit ON the
    symmetry axis (xdot~0) -- structurally the same point-type as Anderson &
    Lo's own Table 2 self-connection and #766's/#767's own primary hits.
    """
    conn_p = _correct_known_hit(system, node, _PRIMARY, _PRIMARY_TAU)
    conn_s = _correct_known_hit(system, node, _SECONDARY, _SECONDARY_TAU)
    assert conn_p.converged and conn_s.converged
    assert abs(float(conn_p.crossing_xv[1])) < 1e-8
    assert abs(float(conn_s.crossing_xv[1])) < 1e-8


# ---------------------------------------------------------------------------
# (4) find_homoclinic's own scan/ghost-guard/ranking plumbing.
# ---------------------------------------------------------------------------


def test_find_homoclinic_returns_known_primary_combo(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """find_homoclinic's own scan, restricted to a narrow (branch, k) window
    around the known PRIMARY combination (for runtime), returns it,
    ghost-guard-passed, ranked by residual (never by distance to a
    nonexistent published target -- see module docstring). This independently
    re-derives PRIMARY via the module's OWN internal grid-scan path (not the
    tau0-seeded shortcut the tests above use), corroborating the hit through
    a genuinely different code path.
    """
    candidates = ntc.find_homoclinic(
        system,
        node,
        branches=(+1, -1),
        k_range=range(13, 14),
        max_time_factor=6.0,
        scan_n=4,
        max_iter=20,
        tol=1e-7,
    )
    assert len(candidates) >= 1
    cand = candidates[0]
    assert cand.connection.converged
    assert cand.ghost_distance > ntc.GHOST_GUARD_DELTA
    assert cand.connection.branch_u == 1
    assert cand.connection.branch_s == -1
    assert abs(float(cand.connection.crossing_xv[0]) - 1.16933872) < 1e-4


def test_gate_empty_scan_is_honest(system: cr3bp.CR3BPSystem, node: jrc.ResonantNode) -> None:
    """No surviving candidates in a deliberately-empty scan window -> an
    honest, explicit empty list, never a fabricated hit. Low k (1..2) lands
    close to the orbit's own trivial section points (see the results note's
    own crossing-time dump: k=1,2 sit near the base orbit's own extent) --
    ghost-shadowed, never converging to a genuine hit here.
    """
    candidates = ntc.find_homoclinic(
        system,
        node,
        branches=(+1,),
        k_range=range(1, 2),
        max_time_factor=3.0,
        scan_n=4,
        max_iter=10,
        tol=1e-6,
    )
    assert candidates == []


# ---------------------------------------------------------------------------
# (5) EPSILON/GHOST_GUARD_DELTA are carried over from the Jovian module, not
# independently sourced -- documented explicitly (module docstring).
# ---------------------------------------------------------------------------


def test_epsilon_and_ghost_guard_reuse_jovian_module_values() -> None:
    assert ntc.EPSILON == jrc.ANDERSON_LO_EPSILON
    assert ntc.GHOST_GUARD_DELTA == jrc.GHOST_GUARD_DELTA
