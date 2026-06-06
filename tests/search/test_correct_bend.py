"""M-ED Phase 1: post-hoc bend feasibility (plan Phase 1 Task 1.3)."""

from __future__ import annotations

from cyclerfinder.search.correct import _bend_deg, _max_bend_deg


def test_max_bend_decreases_with_vinf() -> None:
    # Higher V_inf -> tighter max turn at the same body (Mars).
    assert _max_bend_deg(3.0, "M") > _max_bend_deg(8.0, "M")


def test_bend_zero_for_parallel_vectors() -> None:
    assert _bend_deg((1.0, 0.0, 0.0), (2.0, 0.0, 0.0)) == 0.0
