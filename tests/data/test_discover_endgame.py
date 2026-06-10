"""Endgame-genome discovery path (plan 2026-06-09, Component 3)."""

from __future__ import annotations

from cyclerfinder.data.discover_novel import (
    discover_endgame_moon,
    endgame_route_to_nbody_request,
    saturnian_titan_tour_topologies,
)
from cyclerfinder.search.endgame_graph import InterMoonTransfer, solve_endgame
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


def test_endgame_route_to_nbody_request_shape() -> None:
    route = solve_endgame(
        moon_system="Jupiter",
        entry_moon="Europa",
        target_moon="Europa",
        vinf_entry_kms=2.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=3.0,
        system_moons=("Europa",),
    )
    assert route is not None
    req = endgame_route_to_nbody_request(route, center="Jupiter", moon="Europa")
    assert req["center"] == "Jupiter"
    assert req["n_maneuvers"] == len(route.leveraging_legs)
