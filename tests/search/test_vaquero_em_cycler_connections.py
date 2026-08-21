"""Tests for `#822`: Vaquero 2:1 <-> 3:1 Earth-Moon free transfer (heteroclinic)
(docs/notes/2026-08-11-822-vaquero-em-free-transfer.md).

Sourced-golden discipline: Vaquero (2013, pp. 171-172) asserts the connection's
EXISTENCE but prints NO transfer state, so there is no published digit to gate
against -- every test below asserts SELF-CONSISTENCY evidence directly (Newton
residual, full-4-state gap incl. the ydot-sign hard gate, ghost-guard margin,
independent Radau cross-check, forward/backward re-approach), per the
`#766`/`#767`/`#781` precedent. The vendored KNOWN-HIT constants below are this
task's own converged values (DERIVE, archived in
``data/found/822_vaquero_em_free_transfer/results.json``) used as regression
targets + Newton starts, never as sourced goldens. None are ``@pytest.mark.slow``
-- verdict-bearing evidence tests must run in CI
(``feedback_delegation_fresh_agent_not_fork``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.vaquero_em_cycler_connections as vcc
import cyclerfinder.search.vaquero_em_cyclers as vem
from cyclerfinder.genome.heteroclinic_cycle import HeteroclinicConnection


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


@pytest.fixture(scope="module")
def node21(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    return vcc.build_vaquero_overlap_node(system, 2, 2.60)


@pytest.fixture(scope="module")
def node31(system: cr3bp.CR3BPSystem) -> jrc.ResonantNode:
    return vcc.build_vaquero_overlap_node(system, 3, 2.60)


# ---------------------------------------------------------------------------
# (1) Sourced/derived constants.
# ---------------------------------------------------------------------------


def test_overlap_band_is_intersection_of_sourced_family_ranges() -> None:
    """OVERLAP_C_RANGE is DERIVED from the two SOURCED printed C ranges (#799's
    own vendored Vaquero p.171 prose values) -- and equals the [2.54, 2.66]
    band Vaquero's own free-transfer prose names (pp. 171-172)."""
    assert vcc.OVERLAP_C_RANGE == (2.54, 2.66)
    assert vcc.OVERLAP_C_RANGE[0] == max(vem.VAQUERO_C_RANGE_21[0], vem.VAQUERO_C_RANGE_31[0])
    assert vcc.OVERLAP_C_RANGE[1] == min(vem.VAQUERO_C_RANGE_21[1], vem.VAQUERO_C_RANGE_31[1])


def test_epsilon_and_ghost_guard_reuse_sibling_module_values() -> None:
    """Neither is Vaquero-sourced; both are the sibling connection modules'
    own generic defaults, reused verbatim (module docstring)."""
    assert vcc.EPSILON == jrc.ANDERSON_LO_EPSILON
    assert vcc.GHOST_GUARD_DELTA == jrc.GHOST_GUARD_DELTA


def test_overlap_grid_covers_both_families_at_every_shared_grid_point() -> None:
    """13 shared 0.01-grid points per family, both endpoints inclusive."""
    cs = sorted({c for (_p, c) in vcc.OVERLAP_GRID_ICS})
    assert cs == [round(2.54 + 0.01 * i, 4) for i in range(13)]
    for c in cs:
        assert (2, c) in vcc.OVERLAP_GRID_ICS
        assert (3, c) in vcc.OVERLAP_GRID_ICS


# ---------------------------------------------------------------------------
# (2) Node re-derivation (never a raw trust of the vendored guesses).
# ---------------------------------------------------------------------------


def test_nodes_rederive_at_equal_jacobi(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> None:
    """Both nodes re-converge at EXACTLY the target C (the corrector holds C
    fixed) -- the equal-energy precondition of a heteroclinic connection, and
    of Vaquero's own free-transfer statement."""
    assert node21.converged and node31.converged
    assert abs(node21.jacobi - 2.60) < 1e-12
    assert abs(node31.jacobi - 2.60) < 1e-12
    # Fresh correction reproduces #799's archived ICs (not merely echoes them):
    x0_21, ydot0_21, t_21, _ = vcc.OVERLAP_GRID_ICS[(2, 2.60)]
    assert abs(float(node21.state0[0]) - x0_21) < 1e-10
    assert abs(float(node21.state0[4]) - ydot0_21) < 1e-10
    assert abs(node21.period - t_21) < 1e-9


def test_node_21_saddle_is_negative_real_and_31_positive_real(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> None:
    """The 2:1 family's saddle eigenvalue is genuinely NEGATIVE real across the
    band (Barden nu < 0 -- `#781`'s 4:5-saddle situation, module docstring);
    the 3:1's is positive real. Checked from the actual full-period planar
    monodromy, not the magnitude-only Floquet helper."""
    for node, lo, hi in ((node21, -5.0, -4.9), (node31, 11.9, 12.1)):
        arc = cr3bp.propagate(
            system, node.state0, node.period, with_stm=True, rtol=1e-13, atol=1e-13
        )
        assert arc.stm is not None
        phi4 = arc.stm[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])]
        eigs = np.linalg.eigvals(phi4)
        lam = eigs[int(np.argmax(np.abs(eigs)))]
        assert abs(float(np.imag(lam))) < 1e-6
        assert lo < float(np.real(lam)) < hi


def test_node_build_rejects_out_of_band_request(system: cr3bp.CR3BPSystem) -> None:
    with pytest.raises(KeyError):
        vcc.build_vaquero_overlap_node(system, 2, 2.30)


# ---------------------------------------------------------------------------
# (3) THE KNOWN HIT (this task's primary verified free transfer, C=2.60):
# Wu(2:1) -> Ws(3:1), branch (+1, +1), k=(45, 30), crossing at
# (x, xdot) = (0.86597176, -0.47992143), ydot < 0 on both legs. DERIVE
# regression targets from this task's own sweep (results.json).
# ---------------------------------------------------------------------------

_PRIMARY_TAU_U = 3.6477845733834244
_PRIMARY_TAU_S = 5.846080701611376
_PRIMARY_X = 0.86597176
_PRIMARY_XDOT = -0.47992143


def _refine_primary(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> HeteroclinicConnection:
    seed = vcc.ConnectionSeed(
        distance=0.0029,
        branch_u=+1,
        branch_s=+1,
        k_u=45,
        k_s=30,
        tau_u=_PRIMARY_TAU_U,
        tau_s=_PRIMARY_TAU_S,
        x=_PRIMARY_X,
        xdot=_PRIMARY_XDOT,
        ydot=-0.607,
    )
    return vcc.refine_connection(system, node21, node31, seed)


def test_primary_hit_converges(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> None:
    conn = _refine_primary(system, node21, node31)
    assert conn.converged, conn.notes
    assert conn.residual < 1e-8
    assert abs(float(conn.crossing_xv[0]) - _PRIMARY_X) < 1e-5
    assert abs(float(conn.crossing_xv[1]) - _PRIMARY_XDOT) < 1e-5


def test_primary_hit_full_verification_battery(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> None:
    """The complete evidence bundle in one pass: full-4-state gap (incl. the
    ydot-sign hard gate the 2-D section residual cannot see), ghost margins
    against BOTH orbits, independent Radau re-derivation, forward/backward
    re-approach, per-leg Jacobi drift. Margins asserted well past the gates
    (feedback_verify_automated_ghost_guard_booleans: never trust a razor-thin
    boolean)."""
    conn = _refine_primary(system, node21, node31)
    assert conn.converged
    ev = vcc.verify_connection(system, node21, node31, conn)
    assert ev.passed, ev.notes

    # Real margins, not borderline passes (observed: 5.7e-6 / 4.2e-6 / 5.2e-8):
    assert ev.full_state_gap < 0.2 * vcc.FULL_STATE_GAP_TOL
    assert ev.radau_gap < 0.2 * vcc.RADAU_GAP_TOL
    assert ev.backward_distance < 0.1 * vcc.BACKWARD_TOL
    assert ev.forward_distance < vcc.FORWARD_TOL
    # ydot-sign gate: both legs arrive with ydot < 0 (a genuine full-state
    # match, not the |ydot| section ambiguity).
    assert ev.ydot_signs_match
    assert ev.ydot_u < 0.0 and ev.ydot_s < 0.0
    # Ghost margins: ~480x-520x the guard radius.
    assert ev.ghost_distance_from > 400.0 * vcc.GHOST_GUARD_DELTA
    assert ev.ghost_distance_to > 400.0 * vcc.GHOST_GUARD_DELTA
    # Energy: both legs conserve their seed's Jacobi to ~1e-11, and the seeds'
    # own O(eps^2) energy-surface offsets are at the eps=0.5e-5 quadratic scale.
    assert ev.jacobi_drift_u < vcc.JACOBI_DRIFT_TOL
    assert ev.jacobi_drift_s < vcc.JACOBI_DRIFT_TOL
    assert ev.seed_jacobi_offset_u < 100.0 * vcc.EPSILON**2
    assert ev.seed_jacobi_offset_s < 100.0 * vcc.EPSILON**2
    # Transit times: unstable leg forward (~43.5 nd), stable leg backward.
    assert ev.t_u > 0.0 > ev.t_s


def test_verify_connection_rejects_unconverged(
    system: cr3bp.CR3BPSystem, node21: jrc.ResonantNode, node31: jrc.ResonantNode
) -> None:
    bogus = HeteroclinicConnection(
        orbit_from=node21.label,
        orbit_to=node31.label,
        jacobi=node21.jacobi,
        tau_u=1.0,
        tau_s=1.0,
        k_u=1,
        k_s=1,
        branch_u=1,
        branch_s=1,
        crossing_xv=np.zeros(2),
        residual=float("inf"),
        converged=False,
        n_iter=0,
    )
    with pytest.raises(ValueError, match="converged"):
        vcc.verify_connection(system, node21, node31, bogus)


# ---------------------------------------------------------------------------
# (4) Seed-pairing plumbing: the ydot-sign gate and ordering.
# ---------------------------------------------------------------------------


def test_find_connection_seeds_drops_opposite_ydot_sign_pairs() -> None:
    """An (x, xdot)-identical pair with OPPOSITE ydot signs is NOT a candidate
    (module docstring, 'THE ydot-SIGN GATE') -- only the same-sign partner
    qualifies, however much farther it sits."""
    cu = [vcc.ManifoldCrossing(tau=0.1, k=3, t=5.0, x=0.5, xdot=0.1, ydot=-0.7)]
    cs_opp = [vcc.ManifoldCrossing(tau=0.2, k=4, t=-5.0, x=0.5, xdot=0.1, ydot=+0.7)]
    cs_same = [vcc.ManifoldCrossing(tau=0.3, k=5, t=-6.0, x=0.501, xdot=0.1, ydot=-0.7)]

    assert vcc.find_connection_seeds(cu, cs_opp, branch_u=1, branch_s=1) == []
    seeds = vcc.find_connection_seeds(cu, cs_opp + cs_same, branch_u=1, branch_s=1)
    assert len(seeds) == 1
    assert seeds[0].k_s == 5  # the same-sign partner, not the closer opposite-sign one
    assert seeds[0].distance == pytest.approx(0.001)


def test_find_connection_seeds_sorted_and_distance_capped() -> None:
    cu = [
        vcc.ManifoldCrossing(tau=0.1, k=1, t=1.0, x=0.5, xdot=0.0, ydot=-1.0),
        vcc.ManifoldCrossing(tau=0.2, k=2, t=2.0, x=1.5, xdot=0.0, ydot=-1.0),
        vcc.ManifoldCrossing(tau=0.3, k=3, t=3.0, x=9.0, xdot=0.0, ydot=-1.0),
    ]
    cs = [
        vcc.ManifoldCrossing(tau=0.4, k=1, t=-1.0, x=0.508, xdot=0.0, ydot=-1.0),
        vcc.ManifoldCrossing(tau=0.5, k=2, t=-2.0, x=1.501, xdot=0.0, ydot=-1.0),
    ]
    seeds = vcc.find_connection_seeds(cu, cs, branch_u=1, branch_s=-1, max_distance=0.02)
    assert len(seeds) == 2  # the x=9.0 crossing has no partner within max_distance
    assert seeds[0].distance <= seeds[1].distance
    assert seeds[0].k_u == 2 and seeds[0].k_s == 2
    assert seeds[0].branch_u == 1 and seeds[0].branch_s == -1


def test_manifold_section_crossings_short_horizon_plumbing(
    system: cr3bp.CR3BPSystem, node31: jrc.ResonantNode
) -> None:
    """Cheap plumbing check on a short horizon: crossings are indexed 1..k in
    time order per trajectory, carry the signed clock (negative for the stable
    direction), and the ~6-crossings-per-period rate of this family's own
    topology shows up."""
    crossings = vcc.manifold_section_crossings(
        system, node31, direction="stable", branch=+1, n_tau=2, n_periods=1.5
    )
    assert crossings
    per_tau: dict[float, list[vcc.ManifoldCrossing]] = {}
    for c in crossings:
        per_tau.setdefault(c.tau, []).append(c)
    assert len(per_tau) == 2
    for pts in per_tau.values():
        assert [c.k for c in pts] == list(range(1, len(pts) + 1))
        assert all(c.t < 0.0 for c in pts)  # stable leg runs backward
        # At least the base orbit's own ~6-crossings-per-period rate over 1.5
        # periods; manifold trajectories can add near-tangential extra
        # crossings beyond it (observed: 9-11), so only loosely bounded above.
        assert 7 <= len(pts) <= 18


# ---------------------------------------------------------------------------
# (5) `#840`: REVERSE direction (Wu(3:1) -> Ws(2:1)) round-trip
# demonstrations at the two catalogued band-edge Jacobi constants, seeded
# from this task's own recorded phase indices -- mirrors the primary-hit
# pattern above (3), never the full blind seed search in CI.
# ---------------------------------------------------------------------------


def test_predict_reverse_crossing_flips_xdot_sign_only() -> None:
    """CR3BP time-reversal symmetry ``(x, y, xdot, ydot, t) -> (x, -y,
    -xdot, ydot, -t)``: on the ``{y=0}`` section this flips ``xdot`` and
    keeps ``x`` fixed."""
    assert vcc.predict_reverse_crossing((0.5, -1.2)) == (0.5, 1.2)
    assert vcc.predict_reverse_crossing((-0.3, 0.9)) == (-0.3, -0.9)


def test_reverse_c254_converges_and_matches_symmetry_prediction(
    system: cr3bp.CR3BPSystem,
) -> None:
    """The C=2.54 reverse connection (`#840`): converges immediately at the
    sibling-module default epsilon on its own recorded phase indices, and
    lands almost exactly on the CR3BP time-reversal prediction of the
    FORWARD entry's own crossing (~6e-9) -- genuinely the mirror-image
    manifold intersection, not a coincidentally nearby different one."""
    node21 = vcc.build_vaquero_overlap_node(system, 2, 2.54)
    node31 = vcc.build_vaquero_overlap_node(system, 3, 2.54)
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=-1,
        branch_s=1,
        k_u=26,
        k_s=38,
        tau_u=4.258037700729926,
        tau_s=0.316633984139843,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=1.0e-4)
    assert conn.converged, conn.notes
    assert conn.residual < 1e-8
    assert abs(float(conn.crossing_xv[0]) - 0.23322508223108054) < 1e-6
    assert abs(float(conn.crossing_xv[1]) - 1.7855964271675904) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=1.0e-4)
    assert ev.passed, ev.notes
    assert ev.ydot_signs_match
    assert ev.full_state_gap < 1e-5
    assert ev.forward_distance < vcc.FORWARD_TOL
    assert ev.jacobi_drift_u < vcc.JACOBI_DRIFT_TOL
    assert ev.jacobi_drift_s < vcc.JACOBI_DRIFT_TOL

    # The forward entry's own recorded crossing (em-vaquero-hetero-wu21c254-
    # ws31c254-2026, data/manifold_connections.yaml) predicts this reverse
    # crossing via time-reversal symmetry to within ~6e-9.
    predicted = vcc.predict_reverse_crossing((0.23322507758762506, -1.785596430399187))
    dist = math.hypot(
        float(conn.crossing_xv[0]) - predicted[0], float(conn.crossing_xv[1]) - predicted[1]
    )
    assert dist < 1e-6


def test_reverse_c266_converges_at_larger_epsilon_different_crossing(
    system: cr3bp.CR3BPSystem,
) -> None:
    """The C=2.66 reverse connection (`#840`): the sibling-module default
    ``epsilon=1e-4`` converges on this candidate but fails
    ``verify_connection``'s forward-reapproach gate (the same documented
    `#822` evidence-quality tradeoff), diagnosed and fixed by raising
    ``epsilon`` to ``2e-4`` -- which converges on a DIFFERENT crossing from
    the symmetry prediction (this family pair's own "rich" intersection
    structure, module docstring), still a genuine verified connection."""
    node21 = vcc.build_vaquero_overlap_node(system, 2, 2.66)
    node31 = vcc.build_vaquero_overlap_node(system, 3, 2.66)
    seed = vcc.ConnectionSeed(
        distance=0.0,
        branch_u=-1,
        branch_s=1,
        k_u=29,
        k_s=33,
        tau_u=4.122557122466014,
        tau_s=2.595563242258196,
        x=0.0,
        xdot=0.0,
        ydot=0.0,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=2.0e-4)
    assert conn.converged, conn.notes
    assert conn.residual < 1e-8
    assert abs(float(conn.crossing_xv[0]) - 0.9590289887986131) < 1e-6
    assert abs(float(conn.crossing_xv[1]) - 0.32719383414376546) < 1e-6

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=2.0e-4)
    assert ev.passed, ev.notes
    assert ev.ydot_signs_match
    assert ev.forward_distance < vcc.FORWARD_TOL
    assert ev.jacobi_drift_u < vcc.JACOBI_DRIFT_TOL
    assert ev.jacobi_drift_s < vcc.JACOBI_DRIFT_TOL

    # The forward entry's own recorded crossing (em-vaquero-hetero-wu21c266-
    # ws31c266-2026) predicts a DIFFERENT crossing than the one actually
    # verified here -- confirmed genuinely different, not merely epsilon-shifted.
    predicted = vcc.predict_reverse_crossing((0.8691637553536629, -0.5754661265104016))
    dist = math.hypot(
        float(conn.crossing_xv[0]) - predicted[0], float(conn.crossing_xv[1]) - predicted[1]
    )
    assert dist > 0.1


def test_reverse_c266_at_default_epsilon_fails_forward_gate_not_other_gates(
    system: cr3bp.CR3BPSystem,
) -> None:
    """Regression pin for the diagnosed C=2.66 pathology this task found:
    the closest-to-predicted candidate at the sibling-module default
    ``epsilon=1e-4`` Newton-converges cleanly and passes every OTHER gate,
    but fails specifically the forward-reapproach evidence-quality gate --
    confirming the fix (raising epsilon) targets the right mechanism rather
    than papering over a genuine non-connection."""
    node21 = vcc.build_vaquero_overlap_node(system, 2, 2.66)
    node31 = vcc.build_vaquero_overlap_node(system, 3, 2.66)
    seed = vcc.ConnectionSeed(
        distance=0.0022718899144641,
        branch_u=-1,
        branch_s=-1,
        k_u=32,
        k_s=37,
        tau_u=2.8801447617277045,
        tau_s=0.07865060534415375,
        x=0.8669622060749299,
        xdot=0.5769024162175868,
        ydot=-0.45520176895324843,
    )
    conn = vcc.refine_connection(system, node31, node21, seed, epsilon=1.0e-4)
    assert conn.converged, conn.notes

    ev = vcc.verify_connection(system, node31, node21, conn, epsilon=1.0e-4)
    assert not ev.passed
    assert "forward" in ev.notes
    # Every OTHER gate passes with good margin -- this is a real,
    # zero-Delta-v heteroclinic geometry, just an evidence-quality-floor
    # issue at this epsilon (module docstring, "Manifold-offset epsilon as
    # an evidence-quality control"), not a spurious Newton false-positive.
    assert ev.ydot_signs_match
    assert ev.full_state_gap < vcc.FULL_STATE_GAP_TOL
    assert ev.radau_gap < vcc.RADAU_GAP_TOL
    assert ev.backward_distance < vcc.BACKWARD_TOL
    assert ev.forward_distance > vcc.FORWARD_TOL  # the one gate that fails
