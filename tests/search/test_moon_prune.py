"""Phase 6 Phase 1: moon-pair leg admissibility prune (plan Phase 1 Tasks 1.0-1.2)."""

from __future__ import annotations

from cyclerfinder.search.moon_prune import moon_leg_admissible


def test_europa_ganymede_leg_returns_reasoned_verdict() -> None:
    ok, reason = moon_leg_admissible(
        "Europa",
        "Ganymede",
        vinf_kms=4.0,
        budget_kms=10.0,
        primary="Jupiter",
    )
    assert isinstance(ok, bool)
    assert reason  # non-empty: the prune must record WHY


def test_zero_budget_prunes_on_vilm_floor() -> None:
    ok, reason = moon_leg_admissible(
        "Europa",
        "Ganymede",
        vinf_kms=4.0,
        budget_kms=0.0,
        primary="Jupiter",
    )
    assert ok is False
    assert "vilm" in reason.lower() or "floor" in reason.lower()


def test_prune_topology_records_per_leg_reasons() -> None:
    from cyclerfinder.data.discover_novel import TopologySpec
    from cyclerfinder.search.moon_prune import prune_topology_legs

    spec = TopologySpec(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single",) * 3,
        period_k=1,
        period_sec=1.0e6,
        tof_seed_days=(3.5, 7.2),
        slack_leg=2,
    )
    survives, reasons = prune_topology_legs(
        spec,
        vinf_seed_kms=4.0,
        budget_kms=10.0,
        primary="Jupiter",
    )
    assert isinstance(survives, bool)
    assert len(reasons) == len(spec.sequence) - 1  # one per leg


def test_prune_keeps_the_known_closing_ieg_family() -> None:
    from cyclerfinder.data.discover_novel import TopologySpec
    from cyclerfinder.search.moon_prune import prune_topology_legs

    spec = TopologySpec(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single",) * 3,
        period_k=1,
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        tof_seed_days=(3.5, 7.2),
        slack_leg=2,
    )
    # Budget set generously above the VILM floor: the known-closing family must
    # NOT be pruned (Bellman: never discard a feasible candidate).
    survives, _ = prune_topology_legs(
        spec,
        vinf_seed_kms=4.0,
        budget_kms=50.0,
        primary="Jupiter",
    )
    assert survives is True
