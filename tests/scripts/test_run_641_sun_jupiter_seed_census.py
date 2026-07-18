"""#641: unit tests for the new clustering/degenerate-equilibrium logic factored
out for the Sun-Jupiter generative-seed census script.

These are the "reusable code" pieces this task's own dispatch called out for
tests: :func:`is_degenerate_equilibrium` (the trivial-Lagrange-point-fixed-
point detector this task's own pilot run discovered was needed),
:func:`lagrange_point_label`, and :func:`cluster_genuine_orbits` (the light
Jacobi+period+geometry clustering pass).
"""

from __future__ import annotations

import numpy as np

import scripts.run_641_sun_jupiter_seed_census as run


def _seed(
    *,
    jacobi: float,
    period: float,
    state0: list[float],
    residual: float = 1e-12,
) -> run.GeneratedSeed:
    """Build a minimal converged+sane GeneratedSeed for testing cluster logic."""
    arr = np.asarray(state0, dtype=np.float64)
    return run.GeneratedSeed(
        state0_guess=arr,
        period_guess=period,
        converged=True,
        physically_sane=True,
        state0=arr,
        period=period,
        jacobi=jacobi,
        residual=residual,
        stability_index=1.0,
        stability_note="test",
    )


# ---------------------------------------------------------------------------
# is_degenerate_equilibrium
# ---------------------------------------------------------------------------


def test_zero_velocity_state_is_equilibrium() -> None:
    # Exactly the shape observed at L4 in this task's own pilot run: position
    # at (0.5-mu, +-sqrt(3)/2), velocity at machine-precision noise.
    state0 = [0.499046, 0.866025, 0.0, 1e-14, -2e-15, 3e-16]
    assert run.is_degenerate_equilibrium(state0) is True


def test_nonzero_velocity_state_is_not_equilibrium() -> None:
    # Observed genuine orbit from this task's own pilot run (cluster 1).
    state0 = [0.9455745, -0.1004817, 0.0, -0.1186207, 0.1284685, 0.0]
    assert run.is_degenerate_equilibrium(state0) is False


def test_threshold_is_a_hard_boundary() -> None:
    below = [0.5, 0.5, 0.0, 5e-7, 0.0, 0.0]
    above = [0.5, 0.5, 0.0, 5e-6, 0.0, 0.0]
    assert run.is_degenerate_equilibrium(below, vnorm_threshold=1e-6) is True
    assert run.is_degenerate_equilibrium(above, vnorm_threshold=1e-6) is False


# ---------------------------------------------------------------------------
# lagrange_point_label
# ---------------------------------------------------------------------------


def test_l4_l5_labels_from_y_sign() -> None:
    mu = 0.001
    assert run.lagrange_point_label([0.499, 0.866, 0, 0, 0, 0], mu) == "L4"
    assert run.lagrange_point_label([0.499, -0.866, 0, 0, 0, 0], mu) == "L5"


def test_collinear_labels_from_x_position() -> None:
    mu = 0.001
    # Secondary at x = 1 - mu = 0.999; beyond it (x > 1-mu) is L2.
    assert run.lagrange_point_label([1.07, 0.0, 0, 0, 0, 0], mu) == "L2"
    # Between primary (-mu) and secondary (1-mu) is L1.
    assert run.lagrange_point_label([0.9, 0.0, 0, 0, 0, 0], mu) == "L1"
    # Beyond the primary in the negative direction (x < -mu) is L3.
    assert run.lagrange_point_label([-1.0, 0.0, 0, 0, 0, 0], mu) == "L3"


# ---------------------------------------------------------------------------
# cluster_genuine_orbits
# ---------------------------------------------------------------------------


def test_identical_jacobi_period_geometry_merge_into_one_cluster() -> None:
    a = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])
    b = _seed(jacobi=3.0001, period=1.5001, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])
    clusters = run.cluster_genuine_orbits([a, b])
    assert len(clusters) == 1
    assert len(clusters[0].members) == 2


def test_different_jacobi_gives_distinct_clusters() -> None:
    a = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])
    b = _seed(jacobi=3.3, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])
    clusters = run.cluster_genuine_orbits([a, b])
    assert len(clusters) == 2


def test_different_z_sign_gives_distinct_clusters() -> None:
    # Same Jacobi/period, but a "northern" vs "southern" halo-like branch
    # (opposite z0 sign) is a qualitatively different geometry.
    north = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.05, 0.1, 0.1, 0.0])
    south = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, -0.05, 0.1, 0.1, 0.0])
    clusters = run.cluster_genuine_orbits([north, south])
    assert len(clusters) == 2


def test_clusters_sorted_largest_first() -> None:
    big = [_seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0]) for _ in range(3)]
    small = [_seed(jacobi=3.3, period=2.0, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])]
    clusters = run.cluster_genuine_orbits(big + small)
    assert len(clusters[0].members) == 3
    assert len(clusters[1].members) == 1


def test_representative_is_lowest_residual_member() -> None:
    lo = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0], residual=1e-14)
    hi = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0], residual=1e-9)
    clusters = run.cluster_genuine_orbits([hi, lo])
    assert clusters[0].representative is lo


# ---------------------------------------------------------------------------
# cluster_topology_label
# ---------------------------------------------------------------------------


def test_out_of_plane_seed_gets_halo_label() -> None:
    s = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.08, 0.1, 0.1, 0.0])
    assert run.cluster_topology_label(s) == frozenset({"halo"})


def test_planar_seed_gets_unrestricted_label() -> None:
    s = _seed(jacobi=3.0, period=1.5, state0=[1.0, 0.0, 0.0, 0.1, 0.1, 0.0])
    assert run.cluster_topology_label(s) == frozenset()
