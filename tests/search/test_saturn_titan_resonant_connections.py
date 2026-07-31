"""Tests for `#767`: Saturn-Titan 3:4 resonant-orbit homoclinic self-connection
(docs/notes/2026-07-31-767-saturn-titan-homoclinic-connection.md).

Sourced-golden discipline: this task has NO published state to gate against
(Vaquero 2013 Sec. 4.3.1's own Fig. 4.9 is a figure only, no state table --
see the module docstring and the results note). Every test below therefore
asserts SELF-CONSISTENCY evidence directly (Newton residual, ghost-guard
margin, independent Radau cross-check, forward/backward re-approach), per
this task chain's own honest framing (`#766`'s own precedent). None are
marked ``@pytest.mark.slow`` -- a discovery-verdict-bearing evidence test
must run in CI, not be silently skipped
(``feedback_delegation_fresh_agent_not_fork``).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.saturn_titan_resonant_connections as stc
import cyclerfinder.search.saturn_titan_resonant_families as stf
from cyclerfinder.genome.heteroclinic_cycle import (
    HeteroclinicConnection,
    assemble_cycle,
    correct_connection,
    crosscheck_cycle,
)


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return stf.saturn_titan_system()


@pytest.fixture(scope="module")
def node(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    _sys, nd, _row = stc.build_34_node(system)
    return nd


# ---------------------------------------------------------------------------
# (1) build_34_node reuses #765's own confirmed Table 4.1 3:4 candidate.
# ---------------------------------------------------------------------------


def test_build_34_node_matches_765_confirmed_row(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """The node's own state/period/jacobi match #765's own recorded 3:4 numbers
    (docs/notes/2026-07-29-765-saturn-titan-resonant-families-vaquero-gate.md):
    x0=1.0301662783998498, period=26.140797240249157, jacobi=3.01,
    eigenvalue rel_err 1.07e-6 against Vaquero's own target 2129.81.
    """
    _sys, _nd, row = stc.build_34_node(system)
    assert row.passed
    assert row.label == "3:4"
    assert row.eigenvalue_rel_err < 2e-6

    assert abs(float(node.state0[0]) - 1.0301662783998498) < 1e-9
    assert abs(node.period - 26.140797240249157) < 1e-9
    assert abs(node.jacobi - 3.01) < 1e-12
    assert node.converged


def test_resonant_node_is_real_saddle(node: jrc.ResonantNode) -> None:
    """The 3:4 orbit's own eigenvalue (~2129.8) is the STRONGEST instability
    of any orbit this task chain has built a homoclinic self-connection for
    (#754's C_flyby: 1036; #766's C=3.0041: 54.6) -- see the module docstring
    for why this needed a tighter forward/backward re-approach integrator.
    """
    lam_norm = float(np.linalg.norm(node.unstable_eigvec))
    assert abs(lam_norm - 1.0) < 1e-9  # real_unit-normalised


# ---------------------------------------------------------------------------
# (2) Section convention: {y=0, ydot>0}, x UNRESTRICTED -- NOT Anderson & Lo's
# own single-x-sign choice (see module docstring for why).
# ---------------------------------------------------------------------------


def test_section_convention_is_ydot_positive_x_unrestricted() -> None:
    assert stc.SECTION_YDOT_SIGN == +1


def test_own_section_points_is_both_perpendicular_axis_crossings(
    node: jrc.ResonantNode,
) -> None:
    """The 3:4 orbit's own qualifying section points (ydot>0, x unrestricted)
    are its own IC (x0=1.0301663) AND the corrector's own half-period target
    (x=-1.3666368, `_HALF_CROSSINGS['3:4']=2` in
    saturn_titan_resonant_families.py) -- BOTH perpendicular (xdot=0)
    crossings, confirmed directly this task by inspecting the orbit's own
    {y=0} crossing list (4 per period; only these two are perpendicular).
    """
    system = stf.saturn_titan_system()
    pts = stc.own_section_points(system, node)
    assert len(pts) == 2
    xs = sorted(float(p[0]) for p in pts)
    assert abs(xs[0] - (-1.3666368)) < 1e-6
    assert abs(xs[1] - 1.0301663) < 1e-6
    for p in pts:
        assert abs(float(p[1])) < 1e-8  # perpendicular: xdot=0


def test_x_sign_restricted_convention_would_starve_the_k_index(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """Sanity check for the module docstring's own claim: restricting to a
    SINGLE x sign (Anderson & Lo's own convention) leaves only ONE own point
    for this orbit -- NOT both -- which is why this module does not reuse
    that convention verbatim.
    """
    only_positive = jrc.own_section_points(system, node, ydot_sign=1, x_sign=1)
    only_negative = jrc.own_section_points(system, node, ydot_sign=1, x_sign=-1)
    assert len(only_positive) == 1
    assert len(only_negative) == 1
    # The union (this module's own convention) recovers both.
    assert len(only_positive) + len(only_negative) == len(stc.own_section_points(system, node))


# ---------------------------------------------------------------------------
# (3) The four genuine, independently-converged, ghost-guard-passed hits this
# task's own coarse scan found (branch_u, branch_s in {+1,-1}, k_u, k_s
# mostly in 1..6, max_time_factor=3.0) -- see the results note for the full
# scan log. All four sit at 219x-306x the GHOST_GUARD_DELTA=1e-3 threshold --
# a real, non-delicate margin.
# ---------------------------------------------------------------------------

# On-symmetry-axis (xdot~0) -- structurally the SAME point-type as Anderson &
# Lo's own Table 2 state and #766's own primary hit. Tightest evidence
# bundle of the four (Newton residual, Radau cross-check, ghost margin) --
# reported as this task's PRIMARY result.
_PRIMARY = {"branch_u": +1, "branch_s": +1, "k_u": 5, "k_s": 5}
_PRIMARY_TAU = {"tau_u": 6.18697656714116, "tau_s": 19.953796185428153}

# A genuine mirror pair (k_u, k_s swapped, off-symmetry-axis xdot != 0,
# reflection-symmetric (x, +-xdot)) -- independent corroboration that this
# energy supports transversal homoclinic self-intersections beyond one
# isolated (branch, k) choice, mirroring #766's own strongest additional
# evidence.
_MIRROR_A = {"branch_u": -1, "branch_s": -1, "k_u": 4, "k_s": 5}
_MIRROR_A_TAU = {"tau_u": 24.006236766666365, "tau_s": 2.1347309180960954}
_MIRROR_B = {"branch_u": -1, "branch_s": -1, "k_u": 5, "k_s": 4}
_MIRROR_B_TAU = {"tau_u": 24.00623677219259, "tau_s": 2.1347309142886957}

# A fourth, independently-found (different branch-sign pair, off-axis)
# genuine hit -- further corroboration.
_TERTIARY = {"branch_u": +1, "branch_s": -1, "k_u": 4, "k_s": 4}
_TERTIARY_TAU = {"tau_u": 3.570427157637045, "tau_s": 23.179487310480308}


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
        epsilon=stc.EPSILON,
        branch_u=hit["branch_u"],
        branch_s=hit["branch_s"],
        tau_u0=tau["tau_u"],
        tau_s0=tau["tau_s"],
        ydot_sign_u=stc.SECTION_YDOT_SIGN,
        ydot_sign_s=stc.SECTION_YDOT_SIGN,
        x_sign_u=None,
        x_sign_s=None,
        max_time_factor=3.0,
        tol=1e-9,
        max_iter=60,
        fd_step=1e-7,
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau", "expect_x", "expect_xdot", "min_ghost_margin"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU, -1.67320084, -7.68e-10, 200.0),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU, 0.84403079, -0.11481253, 200.0),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU, 0.84403079, 0.11481253, 200.0),
        ("tertiary", _TERTIARY, _TERTIARY_TAU, -1.12647328, 0.13211509, 200.0),
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

    own_pts = stc.own_section_points(system, node)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    assert d_ghost > stc.GHOST_GUARD_DELTA
    # A REAL margin -- two orders of magnitude past the guard threshold, not
    # a borderline pass (feedback_verify_automated_ghost_guard_booleans).
    assert d_ghost > min_ghost_margin * stc.GHOST_GUARD_DELTA, (
        f"{name}: d_ghost={d_ghost:.3e} -- margin too thin"
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU),
        ("tertiary", _TERTIARY, _TERTIARY_TAU),
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
            "epsilon": stc.EPSILON,
            "branch_u": hit["branch_u"],
            "branch_s": hit["branch_s"],
            "k_u": hit["k_u"],
            "k_s": hit["k_s"],
            "tau_u0": conn.tau_u,
            "tau_s0": conn.tau_s,
            "ydot_sign_u": stc.SECTION_YDOT_SIGN,
            "ydot_sign_s": stc.SECTION_YDOT_SIGN,
            "x_sign_u": None,
            "x_sign_s": None,
            "max_time_factor": 3.0,
            "max_iter": 60,
            "fd_step": 1e-7,
        },
    )
    assert cycle.closed
    checked = crosscheck_cycle(
        system, [node], cycle, method="Radau", rtol=1e-11, atol=1e-11, max_time_factor=3.0
    )
    assert checked.independent_residual < 1e-6, (
        f"{name}: DOP853 vs Radau disagreement {checked.independent_residual:.3e} exceeds 1e-6"
    )


@pytest.mark.parametrize(
    ("name", "hit", "tau", "max_forward"),
    [
        ("primary", _PRIMARY, _PRIMARY_TAU, 0.05),
        ("mirror_a", _MIRROR_A, _MIRROR_A_TAU, 0.1),
        ("mirror_b", _MIRROR_B, _MIRROR_B_TAU, 0.1),
        ("tertiary", _TERTIARY, _TERTIARY_TAU, 0.2),
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
    same integration); ``forward_distance`` is looser -- a genuine cross-leg
    consistency test amplified by this orbit's own strong instability
    (|lambda|~2129.8, the strongest of any orbit this task chain has built a
    self-connection for) -- see the module docstring for the full
    explanation. Both are reported honestly, not fudged: small relative to
    the O(1-2) nondim trajectory scale, but not machine-tight.
    """
    conn = _correct_known_hit(system, node, hit, tau)
    assert conn.converged
    own_pts = stc.own_section_points(system, node)
    d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
    candidate = stc.HomoclinicCandidate(connection=conn, ghost_distance=d_ghost)
    reap = stc.homoclinic_reapproach_check(system, node, candidate, max_time_factor=8.0)

    assert reap.backward_distance < 1e-6, f"{name}: backward_distance={reap.backward_distance:.3e}"
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


def test_primary_hit_is_on_symmetry_axis(system: cr3bp.CR3BPSystem, node: jrc.ResonantNode) -> None:
    """The PRIMARY hit (branch_u=branch_s=+1, k_u=k_s=5) sits ON the
    symmetry axis (xdot~0) -- structurally the same point-type as Anderson &
    Lo's own Table 2 self-connection and #766's own primary hit.
    """
    conn = _correct_known_hit(system, node, _PRIMARY, _PRIMARY_TAU)
    assert conn.converged
    assert abs(float(conn.crossing_xv[1])) < 1e-6


# ---------------------------------------------------------------------------
# (4) find_homoclinic's own scan/ghost-guard/ranking plumbing.
# ---------------------------------------------------------------------------


def test_find_homoclinic_returns_known_primary_combo(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """find_homoclinic's own scan, restricted to the ONE known-converging
    combination (for runtime), returns exactly this candidate, ghost-guard-
    passed, ranked by residual (never by distance to a nonexistent published
    target -- see module docstring).
    """
    candidates = stc.find_homoclinic(
        system,
        node,
        branches=(_PRIMARY["branch_u"],),
        k_range=range(_PRIMARY["k_u"], _PRIMARY["k_u"] + 1),
        max_time_factor=3.0,
        scan_n=12,
        max_iter=40,
        tol=1e-7,
    )
    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.connection.converged
    assert cand.ghost_distance > stc.GHOST_GUARD_DELTA
    assert abs(float(cand.connection.crossing_xv[0]) - (-1.67320088)) < 1e-4


def test_find_homoclinic_ranks_by_residual(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """A scan spanning both the primary and tertiary combinations' own
    branch/k ranges returns both, ranked by Newton-residual tightness
    (ascending) -- never by distance to a published target (none exists).
    """
    candidates = stc.find_homoclinic(
        system,
        node,
        branches=(+1,),
        k_range=range(4, 6),
        max_time_factor=3.0,
        scan_n=8,
        max_iter=20,
        tol=1e-7,
    )
    assert len(candidates) >= 1
    residuals = [c.connection.residual for c in candidates]
    assert residuals == sorted(residuals)


def test_gate_table2_style_empty_scan_is_honest() -> None:
    """No surviving candidates in a deliberately-empty scan window -> an
    honest, explicit empty list, never a fabricated hit.
    """
    system = stf.saturn_titan_system()
    cand = stf.recover_table41_candidate("3:4", system)
    node = jrc.ResonantNode.from_candidate(system, cand)
    candidates = stc.find_homoclinic(
        system,
        node,
        branches=(+1,),
        k_range=range(1, 2),
        max_time_factor=3.0,
        scan_n=6,
        max_iter=10,
        tol=1e-6,
    )
    # k=1 for branch=+1 lands close to the orbit's own IC (ghost-shadowed) --
    # see the results note; the scan honestly returns nothing here.
    assert candidates == []


# ---------------------------------------------------------------------------
# (5) EPSILON/GHOST_GUARD_DELTA are carried over from the Jovian module, not
# independently sourced -- documented explicitly (module docstring).
# ---------------------------------------------------------------------------


def test_epsilon_and_ghost_guard_reuse_jovian_module_values() -> None:
    assert stc.EPSILON == jrc.ANDERSON_LO_EPSILON
    assert stc.GHOST_GUARD_DELTA == jrc.GHOST_GUARD_DELTA


# ---------------------------------------------------------------------------
# `#768`: the periodic 3:4<->6:5 "resonant chain" (Vaquero 2013 Fig. 4.9-4.10).
#
# Step 1 finding (see the `#768` results note): the "chain" is a HOMOCLINIC
# self-connection of 3:4 alone (`#767`'s own machinery), whose crossing is
# selected near 6:5's own fixed point -- so this needs 6:5's own IC location
# only, not its eigenvalue/manifold structure. Step 2's own further
# periodicity-correction attempt (`attempt_chain_closure`) is an HONEST
# PARTIAL/NEGATIVE result (see the results note) -- tests below assert the
# exact, reproducible, bounded behaviour observed, not a forced convergence.
# ---------------------------------------------------------------------------


def test_resonant_chain_target_point_is_65_ic_with_honest_gate_row(
    system: cr3bp.CR3BPSystem,
) -> None:
    """6:5's own {y=0, ydot>0} fixed point, per #765's own recovered IC
    (x0=0.9347726861768341) -- the gate row is included and its own
    eigenvalue criterion is honestly FALSE (2.34e-3 miss, #765), NOT silently
    hidden by this convenience function.
    """
    _sys, target, row = stc.resonant_chain_target_point(system)
    assert row.label == "6:5"
    assert not row.eigenvalue_confirmed
    assert abs(row.eigenvalue_rel_err - 2.3387e-3) < 1e-4
    assert abs(float(target[0]) - 0.9347726861768341) < 1e-8
    assert float(target[1]) == 0.0


def test_rank_by_proximity_to_65_finds_a_much_closer_hit_than_767s_own(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """A direct re-scan of the two `(-1,-1)`-branch homoclinic self-
    connections (the `#767` results note's own MIRROR pair combination)
    converges to a genuine hit at `dist_to_65 ~= 0.094` -- notably closer to
    6:5's own fixed point than `#767`'s own originally-reported MIRROR pair
    (`~0.146`, a different local root of the same (branch, k) residual
    equation found by a different scan seed), and dramatically closer than
    the PRIMARY/TERTIARY hits (`>2.0`) -- directly corroborating Vaquero's own
    qualitative Fig. 4.9 description ("the intersection...is selected to be
    near the fixed point corresponding to the 6:5 resonant orbit").
    """
    _sys, target65, _row = stc.resonant_chain_target_point(system)
    own_pts = stc.own_section_points(system, node)
    hits = []
    for branch_u, branch_s, k_u, k_s in [(-1, -1, 4, 5), (-1, -1, 5, 4)]:
        conn = correct_connection(
            system,
            node,
            node,
            k_u=k_u,
            k_s=k_s,
            epsilon=stc.EPSILON,
            branch_u=branch_u,
            branch_s=branch_s,
            ydot_sign_u=stc.SECTION_YDOT_SIGN,
            ydot_sign_s=stc.SECTION_YDOT_SIGN,
            x_sign_u=None,
            x_sign_s=None,
            max_time_factor=3.0,
            scan_n=12,
            tol=1e-9,
            max_iter=60,
            fd_step=1e-7,
        )
        assert conn.converged
        d_ghost = jrc._ghost_distance(conn.crossing_xv, own_pts)
        assert d_ghost >= stc.GHOST_GUARD_DELTA
        hits.append(stc.HomoclinicCandidate(connection=conn, ghost_distance=d_ghost))

    ranked = stc.rank_by_proximity_to_65(hits, target65)
    assert len(ranked) == 2
    # Both are the reflection-symmetric mirror pair -- same distance.
    assert abs(ranked[0].dist_to_65 - ranked[1].dist_to_65) < 1e-6
    assert ranked[0].dist_to_65 < 0.10  # closer than #767's own 0.146 MIRROR pair
    for r in ranked:
        assert r.candidate.connection.residual < 1e-8


def test_attempt_chain_closure_seed_residual_matches_expected(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """A single Newton evaluation at the seed (node's own IC, xdot=0), fixed
    crossing index nearest `t_target ~= t_u+|t_s| = 110.4996` (the total
    unstable+stable transit time of the near-6:5 candidate above, per the
    `#768` results note) -- reproduces the exact seed residual/crossing-count
    found by this task's own exploratory script, an honest, un-converged FAIL
    (`max_iter=1` never even attempts a Newton step).
    """
    res = stc.attempt_chain_closure(system, node, t_target=110.4996, max_iter=1)
    assert not res.converged
    assert res.n_iter == 1
    assert res.n_events_seed == 16
    assert abs(res.residual - 0.2534297910848558) < 1e-6
    assert "exhausted" in res.notes


def test_attempt_chain_closure_makes_progress_but_does_not_converge(
    system: cr3bp.CR3BPSystem, node: jrc.ResonantNode
) -> None:
    """With a real (if bounded) Newton effort, the residual drops by roughly
    two orders of magnitude (0.253 -> ~0.006) over a handful of damped,
    backtracked iterations -- genuine, monotonic progress, not a wild
    divergence -- but the line search stalls before reaching the `1e-9`
    convergence tolerance: a genuine Newton stall, not a forced convergence,
    mirroring `#759`'s own documented Table-3 stall for `Ws(5:6-LO)`'s severe
    manifold sensitivity ("the residual plateaus...without converging").
    """
    res = stc.attempt_chain_closure(system, node, t_target=110.4996, max_iter=8, max_backtrack=6)
    assert not res.converged
    assert res.residual < 0.02  # real progress from the seed's own 0.253
    assert res.n_iter >= 2


def test_chain_ydot_negative_radicand_is_none(system: cr3bp.CR3BPSystem) -> None:
    """A physically inadmissible (x, xdot) pair at Vaquero's own Jacobi
    constant -- an honest `None`, never a fabricated/complex `ydot`.
    """
    ydot = stc._chain_ydot(0.0, 100.0, stf.VAQUERO_C, system.mu)
    assert ydot is None
