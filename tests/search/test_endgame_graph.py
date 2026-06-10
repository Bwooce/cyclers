# tests/search/test_endgame_graph.py
"""Multi-moon endgame route search (plan 2026-06-09, Component 2)."""

from __future__ import annotations

import pytest

import cyclerfinder.search.endgame_graph as eg
import cyclerfinder.search.leveraging_leg as ll  # noqa: F401 (used in Task 5+)
import cyclerfinder.search.vilm as vilm


def test_lower_bound_is_admissible_vs_vilm() -> None:
    bound = eg.route_lower_bound_kms("Ganymede", "Europa")
    full = vilm.vilm_dv_min("Ganymede", "Europa")
    assert bound <= full + 1e-9
    assert bound > 0.0


def test_intermoon_transfer_ganymede_europa() -> None:
    # Coplanar Hohmann between Ganymede and Europa orbits: depart/arrive V∞ match
    # vilm._hohmann_vinf; ballistic dv ~ 0; positive ToF.
    t = eg.evaluate_intermoon_transfer("Ganymede", "Europa")
    vinf_outer, vinf_inner = vilm._hohmann_vinf("Ganymede", "Europa")
    assert t.vinf_depart_kms == pytest.approx(vinf_outer, abs=1e-6)
    assert t.vinf_arrive_kms == pytest.approx(vinf_inner, abs=1e-6)
    assert t.dv_kms == pytest.approx(0.0, abs=1e-9)
    assert t.tof_days > 0.0
    # Reversed direction swaps depart/arrive.
    r = eg.evaluate_intermoon_transfer("Europa", "Ganymede")
    assert r.vinf_depart_kms == pytest.approx(vinf_inner, abs=1e-6)
    assert r.vinf_arrive_kms == pytest.approx(vinf_outer, abs=1e-6)


def test_solve_endgame_single_moon_descends_at_europa() -> None:
    # Pure endgame: lower V∞ at Europa from 2.0 -> floor 0.8 (entry==target moon).
    route = eg.solve_endgame(
        moon_system="Jupiter",
        entry_moon="Europa",
        target_moon="Europa",
        vinf_entry_kms=2.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=3.0,
        system_moons=("Europa",),
    )
    assert route is not None
    assert route.vinf_final_kms <= 0.8 + 1e-6
    assert route.total_dv_kms >= route.lower_bound_kms - 1e-9
    assert all(leg.gamma_floor_ok for leg in route.leveraging_legs)


def test_solve_endgame_two_moon_tour_ganymede_to_europa() -> None:
    # Multi-moon: enter high at Ganymede, transfer to Europa, capture-feasible.
    route = eg.solve_endgame(
        moon_system="Jupiter",
        entry_moon="Ganymede",
        target_moon="Europa",
        vinf_entry_kms=3.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert route is not None
    assert route.steps[-1] is not None
    # At least one intermoon transfer (Ganymede -> Europa) is in the route.
    assert any(isinstance(s, eg.InterMoonTransfer) for s in route.steps)
    assert route.vinf_final_kms <= 0.8 + 1e-6


def test_solve_endgame_no_route_within_budget_returns_none() -> None:
    route = eg.solve_endgame(
        moon_system="Jupiter",
        entry_moon="Europa",
        target_moon="Europa",
        vinf_entry_kms=2.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=1e-4,
        system_moons=("Europa",),
    )
    assert route is None


def test_dijkstra_matches_brute_force_on_two_moon_grid() -> None:
    bb = eg.solve_endgame(
        moon_system="Jupiter",
        entry_moon="Ganymede",
        target_moon="Europa",
        vinf_entry_kms=3.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    brute = eg._brute_force_optimum(
        moon_system="Jupiter",
        entry_moon="Ganymede",
        target_moon="Europa",
        vinf_entry_kms=3.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert bb is not None and brute is not None
    assert bb.total_dv_kms == pytest.approx(brute, abs=1e-6)


def test_a2_ganymede_europa_min_dv_golden() -> None:
    # A2 (Part-1 Table 1, no-GA): Ganymede->Europa ΔV_min = 1.71 km/s.
    assert vilm.vilm_dv_min("Ganymede", "Europa") == pytest.approx(1.71, abs=0.17)


def test_a6_europa_endgame_scalar_golden() -> None:
    # A6: Europa endgame continuous floor + 46-day phasing scalar.
    dv_ms, days = vilm.europa_endgame_dv()
    assert days == pytest.approx(46.0)
    assert dv_ms < 154.0  # continuous floor < published discrete 3-VILM design


def test_phasefull_route_realises_leverage_above_continuous_floor() -> None:
    # The endgame route is V∞-SHAPING ONLY (a cycler does not capture into a
    # parking orbit), so the right published lower bound is the *leverage* portion
    # of A2 = full ΔV_min MINUS the escape+capture insertion floor. A finite-leg
    # route cannot beat that continuous leverage minimum, and (leverage realised)
    # it stays within a small multiple of it — NOT the un-leveraged impulsive cost.
    route = eg.solve_endgame(
        moon_system="Jupiter",
        entry_moon="Ganymede",
        target_moon="Europa",
        vinf_entry_kms=3.0,
        target_vinf_floor_kms=0.8,
        dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert route is not None
    leverage_floor = vilm.vilm_dv_min("Ganymede", "Europa") - vilm.vilm_dv_floor(
        "Ganymede", "Europa"
    )
    assert route.total_dv_kms >= leverage_floor - 1e-9
    assert route.total_dv_kms < 5.0 * leverage_floor
    assert all(leg.gamma_floor_ok for leg in route.leveraging_legs)
