"""M-ED Phase 1: post-hoc bend feasibility (plan Phase 1 Task 1.3)."""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.correct import _bend_deg, _max_bend_deg, _residual_vector


def test_max_bend_decreases_with_vinf() -> None:
    # Higher V_inf -> tighter max turn at the same body (Mars).
    assert _max_bend_deg(3.0, "M") > _max_bend_deg(8.0, "M")


def test_bend_zero_for_parallel_vectors() -> None:
    assert _bend_deg((1.0, 0.0, 0.0), (2.0, 0.0, 0.0)) == 0.0


def _mars_flyby_nodes(vmag: float, bend_deg: float) -> dict[str, np.ndarray]:
    """One-intermediate (E, M, E) node set: a Mars flyby with equal-magnitude
    in/out V_inf separated by ``bend_deg`` (so magnitude continuity is exact and
    only the feasibility hinge can be nonzero)."""
    theta = np.radians(bend_deg)
    straight = np.array([vmag, 0.0, 0.0])
    return {
        "b0_out": straight,
        "b1_in": straight,
        "b1_out": vmag * np.array([np.cos(theta), np.sin(theta), 0.0]),
        "b2_in": straight,
    }


def test_bend_hinge_is_zero_at_exactly_max_bend() -> None:
    """#140 review: the vector-mode feasibility hinge term is <= epsilon at the
    constructed boundary where the required Mars-flyby bend EXACTLY equals the
    V_inf-limited max single-flyby turn. The hinge is ``max(0, required - max)``,
    so it must vanish at required == max (no infeasibility yet)."""
    vmag = 5.0
    max_turn = _max_bend_deg(vmag, "M")
    nodes = _mars_flyby_nodes(vmag, max_turn)
    res = _residual_vector(nodes, n_encounters=3, mode="vector", sequence=("E", "M", "E"))
    # res = [magnitude-continuity, bend-hinge, closure]; the hinge is res[1].
    assert abs(res[0]) < 1e-12  # equal magnitudes -> exact continuity
    assert res[1] <= 1e-9  # hinge term at exactly max bend


def test_bend_hinge_is_positive_just_beyond_max_bend() -> None:
    """Teeth: a bend even slightly beyond the max turn yields a strictly positive
    hinge (the boundary test above is not vacuous)."""
    vmag = 5.0
    max_turn = _max_bend_deg(vmag, "M")
    nodes = _mars_flyby_nodes(vmag, max_turn + 1.0)
    res = _residual_vector(nodes, n_encounters=3, mode="vector", sequence=("E", "M", "E"))
    assert res[1] > 0.0
