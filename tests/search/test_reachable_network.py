"""Unit tests for the reduced (x, y, theta) reachable-set accessibility scorer.

Mechanics tests (this file) are fully hand-checkable and source-anchored to the
Braik-Ross 2026 equations (Eq. 8 speed relation, Eq. 14 time-reversal, Eq. 26
turn cost) and to elementary graph-theory facts on tiny hand-built graphs. The
heavy C_J=3.1294 method-validation gate lives in
``test_reachable_network_gate.py`` (marked ``slow``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.reachable_network as rn

MU = 0.01215058439469525


# ---------------------------------------------------------------------------
# Reduced model: speed relation (Eq. 8), turn cost (Eq. 26), time-reversal.
# ---------------------------------------------------------------------------


def test_reduced_speed_matches_jacobi_definition() -> None:
    # v = sqrt(-2 Ubar - C_J) must reproduce the Jacobi speed: pick a state,
    # compute C from the full jacobi_constant, then recover v at that position.
    x, y, vx, vy = 0.5, 0.1, 0.2, -0.3
    state = np.array([x, y, 0.0, vx, vy, 0.0])
    c_j = cr3bp.jacobi_constant(state, MU)
    v_expected = math.hypot(vx, vy)
    v = rn.reduced_speed(x, y, MU, c_j)
    assert v == pytest.approx(v_expected, abs=1e-12)


def test_reduced_speed_forbidden_outside_hill_region() -> None:
    # A very high C_J makes most positions forbidden (radicand negative).
    with pytest.raises(ValueError, match="forbidden"):
        rn.reduced_speed(0.5, 0.1, MU, 10.0)
    assert not rn.is_admissible(0.5, 0.1, MU, 10.0)


def test_dv_turn_eq26() -> None:
    # dV_turn = 2 v sin(|d|/2): a 180-deg reversal at speed v costs 2v (full
    # diameter); a 60-deg turn costs 2 v sin(30deg) = v; 0 costs 0.
    v = 1.3
    assert rn.dv_turn(v, math.pi) == pytest.approx(2.0 * v)
    assert rn.dv_turn(v, math.radians(60.0)) == pytest.approx(v)
    assert rn.dv_turn(v, 0.0) == pytest.approx(0.0)
    # Symmetric in sign of delta (|d|).
    assert rn.dv_turn(v, -0.7) == pytest.approx(rn.dv_turn(v, 0.7))


def test_heading_roundtrip() -> None:
    v, theta = 0.9, 1.1
    vx, vy = rn.velocity_from_heading(v, theta)
    assert math.hypot(vx, vy) == pytest.approx(v)
    assert rn.heading(vx, vy) == pytest.approx(theta)


def test_time_reversal_eq14() -> None:
    # R(x, y, theta) = (x, -y, pi - theta), wrapped to (-pi, pi].
    x, y, th = 0.4, 0.3, math.radians(40.0)
    rx, ry, rt = rn.time_reversal(x, y, th)
    assert rx == pytest.approx(x)
    assert ry == pytest.approx(-y)
    assert rt == pytest.approx(rn.wrap_angle(math.pi - th))
    # Involution: applying R twice returns the original (x, y, theta).
    rx2, ry2, rt2 = rn.time_reversal(rx, ry, rt)
    assert (rx2, ry2) == pytest.approx((x, y))
    assert rn.angular_diff(rt2, th) == pytest.approx(0.0, abs=1e-12)


def test_wrap_angle() -> None:
    assert rn.wrap_angle(math.pi) == pytest.approx(math.pi)
    assert rn.wrap_angle(-math.pi) == pytest.approx(math.pi)
    assert rn.wrap_angle(3.0 * math.pi) == pytest.approx(math.pi)
    assert rn.wrap_angle(0.0) == pytest.approx(0.0)
    assert rn.wrap_angle(1.5 * math.pi) == pytest.approx(-0.5 * math.pi)


# ---------------------------------------------------------------------------
# Voxel grid indexing.
# ---------------------------------------------------------------------------


def test_voxel_index_and_center_consistency() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(5.0))
    # A point maps to a voxel whose center is within half a cell of the point.
    x, y, th = 0.37, -0.22, math.radians(12.3)
    idx = grid.index(x, y, th)
    cx, cy, ct = grid.center(idx)
    assert abs(cx - x) <= 0.5 * grid.dx + 1e-12
    assert abs(cy - y) <= 0.5 * grid.dy + 1e-12
    assert abs(rn.angular_diff(ct, th)) <= 0.5 * grid.dtheta + 1e-12


def test_voxel_theta_wraps_modular() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(5.0))
    # theta = pi and theta = -pi index to the same modular theta cell.
    assert grid.index(0.1, 0.1, math.pi)[2] == grid.index(0.1, 0.1, -math.pi)[2]
    assert grid.n_theta == 72


def test_voxel_index_deterministic() -> None:
    grid = rn.VoxelGrid(dx=0.05, dy=0.05, dtheta=math.radians(10.0))
    assert grid.index(0.0, 0.0, 0.0) == grid.index(0.0, 0.0, 0.0)
    # Distinct points in different cells get different indices.
    assert grid.index(0.0, 0.0, 0.0) != grid.index(0.2, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Pairwise overlap -> proxy.
# ---------------------------------------------------------------------------


def test_pair_proxy_no_overlap_inaccessible() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(10.0))
    fa = rn.ReachableSet()
    fa._log((1, 1, 1), 0.1, 0.5)
    bb = rn.ReachableSet()
    bb._log((9, 9, 9), 0.2, 0.5)
    p = rn.pair_proxy(fa, bb, grid)
    assert not p.accessible
    assert math.isinf(p.proxy_dv)


def test_pair_proxy_overlap_sums_costs_plus_patch() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(10.0))
    fa = rn.ReachableSet()
    fa._log((2, 2, 2), 0.10, 0.3, heading=0.20, speed=0.5)
    fa._log((3, 3, 3), 0.40, 0.1)
    bb = rn.ReachableSet()
    bb._log((2, 2, 2), 0.05, 0.2, heading=0.30, speed=0.5)  # shared voxel
    bb._log((7, 7, 7), 0.01, 0.1)
    p = rn.pair_proxy(fa, bb, grid)
    assert p.accessible
    assert p.voxel == (2, 2, 2)
    # Physical patch: local speed (0.5) x actual heading mismatch (|0.20-0.30|).
    patch = rn.dv_turn(0.5, 0.20 - 0.30)
    assert p.proxy_dv == pytest.approx(0.10 + 0.05 + patch)
    assert p.proxy_time == pytest.approx(0.3 + 0.2)


def test_pair_proxy_picks_min_dv_voxel() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(10.0))
    fa = rn.ReachableSet()
    fa._log((1, 1, 1), 0.50, 0.1)
    fa._log((2, 2, 2), 0.05, 0.1)
    bb = rn.ReachableSet()
    bb._log((1, 1, 1), 0.01, 0.1)
    bb._log((2, 2, 2), 0.05, 0.1)
    p = rn.pair_proxy(fa, bb, grid)
    # voxel (2,2,2): 0.05+0.05 = 0.10; voxel (1,1,1): 0.50+0.01 = 0.51 -> pick (2,2,2)
    assert p.voxel == (2, 2, 2)


def test_reachable_log_keeps_min_cost() -> None:
    rs = rn.ReachableSet()
    rs._log((0, 0, 0), 0.5, 1.0)
    rs._log((0, 0, 0), 0.2, 2.0)  # cheaper cost wins even with later time
    rs._log((0, 0, 0), 0.9, 0.1)  # more expensive -> ignored
    assert rs.voxels[(0, 0, 0)] == pytest.approx(0.2)
    assert rs.times[(0, 0, 0)] == pytest.approx(2.0)


def test_mirror_reachable_set_is_time_reversal() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(10.0))
    fwd = rn.ReachableSet()
    idx = grid.index(0.3, 0.2, math.radians(35.0))
    fwd._log(idx, 0.123, 4.5, heading=math.radians(35.0), speed=0.42)
    back = rn.mirror_reachable_set(fwd, grid)
    cx, cy, ct = grid.center(idx)
    mx, my, mt = rn.time_reversal(cx, cy, ct)
    midx = grid.index(mx, my, mt)
    assert midx in back.voxels
    assert back.voxels[midx] == pytest.approx(0.123)
    assert back.times[midx] == pytest.approx(4.5)
    # Heading time-reverses (pi - theta); speed is preserved.
    assert back.headings[midx] == pytest.approx(rn.wrap_angle(math.pi - math.radians(35.0)))
    assert back.speeds[midx] == pytest.approx(0.42)


# ---------------------------------------------------------------------------
# Proxy matrix assembly.
# ---------------------------------------------------------------------------


def test_proxy_matrix_symmetric_zero_diagonal() -> None:
    grid = rn.VoxelGrid(dx=0.01, dy=0.01, dtheta=math.radians(10.0))
    # Two families that share a voxel, one that is isolated.
    a = rn.ReachableSet()
    a._log((1, 1, 1), 0.1, 0.1)
    b = rn.ReachableSet()
    b._log((1, 1, 1), 0.2, 0.1)
    c = rn.ReachableSet()
    c._log((5, 5, 5), 0.1, 0.1)
    fwd = [a, b, c]
    back = [rn.mirror_reachable_set(r, grid) for r in fwd]
    mat = rn.proxy_matrix(fwd, back, grid)
    assert np.allclose(np.diag(mat), 0.0)
    assert np.allclose(mat, mat.T, equal_nan=False)


# ---------------------------------------------------------------------------
# Centralities on a tiny hand-checkable graph.
# ---------------------------------------------------------------------------


def test_centralities_star_graph() -> None:
    # Star: node 0 connected to 1, 2, 3 with unit cost; leaves not interconnected.
    inf = math.inf
    w = np.array(
        [
            [0.0, 1.0, 1.0, 1.0],
            [1.0, 0.0, inf, inf],
            [1.0, inf, 0.0, inf],
            [1.0, inf, inf, 0.0],
        ]
    )
    c = rn.centralities(w)
    # Strength: hub has 3 unit edges -> 3; each leaf has 1 -> 1.
    assert c.strength[0] == pytest.approx(3.0)
    assert c.strength[1] == pytest.approx(1.0)
    # The hub is strictly the strongest node.
    assert c.strength[0] > c.strength[1]
    # Harmonic closeness: hub reaches all at distance 1 -> 3; a leaf reaches hub
    # at 1 and the other two leaves at 2 -> 1 + 0.5 + 0.5 = 2.
    assert c.harmonic_closeness[0] == pytest.approx(3.0)
    assert c.harmonic_closeness[1] == pytest.approx(2.0)
    # Betweenness: hub lies on all 3 leaf-leaf shortest paths -> 3; leaves -> 0.
    assert c.betweenness[0] == pytest.approx(3.0)
    assert c.betweenness[1] == pytest.approx(0.0)
    assert c.betweenness[2] == pytest.approx(0.0)


def test_centralities_path_graph_betweenness() -> None:
    # Path 0 - 1 - 2: middle node 1 is on the single 0-2 shortest path.
    inf = math.inf
    w = np.array(
        [
            [0.0, 1.0, inf],
            [1.0, 0.0, 1.0],
            [inf, 1.0, 0.0],
        ]
    )
    c = rn.centralities(w)
    assert c.betweenness[1] == pytest.approx(1.0)
    assert c.betweenness[0] == pytest.approx(0.0)
    assert c.betweenness[2] == pytest.approx(0.0)


def test_centralities_isolated_node_zero() -> None:
    inf = math.inf
    w = np.array(
        [
            [0.0, 1.0, inf],
            [1.0, 0.0, inf],
            [inf, inf, 0.0],
        ]
    )
    c = rn.centralities(w)
    assert c.strength[2] == pytest.approx(0.0)
    assert c.harmonic_closeness[2] == pytest.approx(0.0)
    assert c.betweenness[2] == pytest.approx(0.0)


def test_strength_prefers_cheaper_edges() -> None:
    # Reciprocal-cost: a node with cheap edges scores higher than one with
    # expensive edges, even at equal degree.
    inf = math.inf
    w = np.array(
        [
            [0.0, 0.1, inf, inf],
            [0.1, 0.0, 2.0, inf],
            [inf, 2.0, 0.0, 2.0],
            [inf, inf, 2.0, 0.0],
        ]
    )
    c = rn.centralities(w)
    # Node 1 has a 0.1 edge (1/0.1=10) + a 2.0 edge (0.5) = 10.5; node 2 has two
    # 2.0 edges = 1.0. Cheap-edge node dominates strength.
    assert c.strength[1] > c.strength[2]


def test_vu_ms_matches_earth_moon_velocity_unit() -> None:
    # VU = LU/TU = 384400 km / (4.34837740 d) ~ 1.0232 km/s (standard EM value).
    assert pytest.approx(1023.16, abs=1.0) == rn.VU_MS


def test_apply_budget_cap_drops_over_budget_edges() -> None:
    # An edge whose m/s cost exceeds dV_cap is dropped (set to inf); a cheap edge
    # is kept and converted to m/s.
    inf = math.inf
    # nondimensional: 0.01 -> ~10 m/s (kept); 1.0 -> ~1023 m/s (> 409.3, dropped).
    w_nd = np.array(
        [
            [0.0, 0.01, 1.0],
            [0.01, 0.0, inf],
            [1.0, inf, 0.0],
        ]
    )
    capped = rn.apply_budget_cap(w_nd, dv_cap_ms=409.3)
    assert capped[0, 1] == pytest.approx(0.01 * rn.VU_MS)
    assert math.isinf(capped[0, 2])  # over the cap -> dropped
    assert math.isinf(capped[1, 2])  # already missing
    assert capped[0, 0] == 0.0


def test_budget_cap_creates_relay_betweenness() -> None:
    # The cap can drop a direct edge and force relay routing -> nonzero betweenness
    # for the relay node. Without the cap (complete graph) betweenness is zero.
    # 0-1, 1-2 moderately cheap; 0-2 direct is cheaper than the 0-1-2 relay sum
    # but still over the cap. Uncapped: direct 0-2 wins (no relay). Capped: 0-2
    # dropped, forcing relay through 1.
    # direct 0-2 = 0.55*VU ~ 563 m/s; relay 0-1-2 = (0.4+0.4)*VU ~ 819 m/s.
    w_nd = np.array(
        [
            [0.0, 0.4, 0.55],
            [0.4, 0.0, 0.4],
            [0.55, 0.4, 0.0],
        ]
    )
    uncapped = rn.centralities(w_nd * rn.VU_MS)
    assert uncapped.betweenness[1] == pytest.approx(0.0)  # direct 0-2 cheaper: no relay
    capped = rn.apply_budget_cap(w_nd, dv_cap_ms=409.3)
    c = rn.centralities(capped)
    assert c.betweenness[1] > 0.0  # node 1 now relays 0<->2


def test_normalized_centralities_preserve_ranks() -> None:
    # The Nf normalization is a per-metric constant -> ranks unchanged vs raw.
    inf = math.inf
    w = np.array(
        [
            [0.0, 0.1, 2.0, inf],
            [0.1, 0.0, 0.5, 2.0],
            [2.0, 0.5, 0.0, 0.5],
            [inf, 2.0, 0.5, 0.0],
        ]
    )
    raw = rn.centralities(w)
    norm = rn.normalized_centralities(w, n_families=4)
    assert np.argsort(-raw.strength).tolist() == np.argsort(-norm.strength).tolist()
    assert (
        np.argsort(-raw.harmonic_closeness).tolist()
        == np.argsort(-norm.harmonic_closeness).tolist()
    )
