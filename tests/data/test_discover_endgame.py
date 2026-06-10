"""Endgame-genome discovery path (plan 2026-06-09, Component 3)."""

from __future__ import annotations

from cyclerfinder.data.discover_novel import (
    discover_endgame_moon,
    saturnian_titan_tour_topologies,
)
from cyclerfinder.search.endgame_graph import InterMoonTransfer
from cyclerfinder.search.leveraging_leg import LeveragingLegResult


def test_discover_endgame_yields_powered_findings_or_clean_empty() -> None:
    findings = list(
        discover_endgame_moon(
            topologies=saturnian_titan_tour_topologies(),
            center="Saturn",
            target_vinf_floor_kms=6.0,
            dv_budget_kms=4.0,
        )
    )
    for f in findings:
        assert f.powered is True
        assert f.endgame_route is not None
        for leg in f.endgame_route.leveraging_legs:
            assert isinstance(leg, LeveragingLegResult)
            assert leg.gamma_floor_ok
        # steps are only legs or transfers.
        for s in f.endgame_route.steps:
            assert isinstance(s, (LeveragingLegResult, InterMoonTransfer))
