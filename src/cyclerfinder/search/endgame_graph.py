# src/cyclerfinder/search/endgame_graph.py
"""Multi-moon endgame route search (spec 2026-06-09, Component 2).

Chains phase-full VILM legs (:mod:`cyclerfinder.search.leveraging_leg`) and
ballistic intermoon transfers to walk a moon-tour cycler's encounter V∞ from a
high entry down to a bend-feasible target floor. Dijkstra over (moon, V∞) states
(non-negative edge costs -> optimal **over the discretised graph**: states are
bucketed at ``_VINF_BUCKET_KMS`` and the candidate V∞ targets are a finite set, so
this is the optimum of the bucketed graph, not the continuous problem). The
phase-free Γ-quadrature (:mod:`cyclerfinder.search.vilm`) supplies the admissible
lower bound. The route is a conservative V∞-shaping lower bound — phasing closure
is confirmed downstream by the n-body step, not asserted here. Pure.
"""

from __future__ import annotations

import heapq
import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search import leveraging_leg as ll
from cyclerfinder.search import vilm

DAY_S = 86400.0


def route_lower_bound_kms(entry_moon: str, target_moon: str) -> float:
    """Admissible ΔV lower bound (km/s): vilm escape+capture insertion floor.

    For a multi-moon transfer (entry != target) this is the escape + capture cost
    at the two endpoints.  For a same-moon endgame (entry == target) the floor is
    trivially 0: no maneuver is required in the degenerate case where V∞ already
    satisfies the floor.  Any actual realised route has ΔV ≥ 0, so 0 is admissible.
    """
    if entry_moon == target_moon:
        return 0.0
    return vilm.vilm_dv_floor(entry_moon, target_moon)


@dataclass(frozen=True)
class InterMoonTransfer:
    """A ballistic coplanar Hohmann transfer between two moons' circular orbits.

    The transfer is gravity-assist-linked (Campagnola endgame model): its ΔV is
    folded into the V∞ bounds, so the standalone ΔV is ~0; the cost is paid by the
    leveraging legs that shape V∞ to the Hohmann departure value.
    """

    moon_from: str
    moon_to: str
    vinf_depart_kms: float
    vinf_arrive_kms: float
    dv_kms: float
    tof_days: float
    gamma_floor_ok: bool = True  # ballistic — trivially satisfied


def evaluate_intermoon_transfer(moon_from: str, moon_to: str) -> InterMoonTransfer:
    """Coplanar Hohmann transfer between ``moon_from`` and ``moon_to`` (km/s, days)."""
    sat_a = SATELLITES[moon_from]
    sat_b = SATELLITES[moon_to]
    mu = PRIMARIES[sat_a.primary]
    outer = moon_from if sat_a.sma_km >= sat_b.sma_km else moon_to
    inner = moon_to if outer == moon_from else moon_from
    vinf_outer, vinf_inner = vilm._hohmann_vinf(outer, inner)
    a_t = 0.5 * (sat_a.sma_km + sat_b.sma_km)
    tof_s = math.pi * math.sqrt(float(a_t**3) / mu)
    if moon_from == outer:
        vdep, varr = vinf_outer, vinf_inner
    else:
        vdep, varr = vinf_inner, vinf_outer
    return InterMoonTransfer(
        moon_from=moon_from,
        moon_to=moon_to,
        vinf_depart_kms=vdep,
        vinf_arrive_kms=varr,
        dv_kms=0.0,
        tof_days=tof_s / DAY_S,
    )


# Candidate resonances tried per leg; the cheapest converged one is used.
_RESONANCES: tuple[tuple[int, int], ...] = ((1, 1), (2, 1), (3, 2), (2, 3), (1, 2))
_VINF_MATCH_TOL_KMS = 1e-3
_VINF_BUCKET_KMS = 1e-3


@dataclass(frozen=True)
class EndgameRoute:
    steps: tuple[object, ...]  # mix of LeveragingLegResult | InterMoonTransfer
    total_dv_kms: float
    vinf_entry_kms: float
    vinf_final_kms: float
    lower_bound_kms: float

    @property
    def leveraging_legs(self) -> tuple[ll.LeveragingLegResult, ...]:
        return tuple(s for s in self.steps if isinstance(s, ll.LeveragingLegResult))


def _system_moons(moon_system: str, override: Sequence[str] | None) -> tuple[str, ...]:
    if override is not None:
        return tuple(override)
    return tuple(m for m, s in SATELLITES.items() if s.primary == moon_system)


def _best_leg(moon: str, vinf_in: float, vinf_out: float) -> ll.LeveragingLegResult | None:
    """Cheapest converged, Γ-floor-respecting leg over the candidate resonances."""
    best: ll.LeveragingLegResult | None = None
    for n, m in _RESONANCES:
        leg = ll.evaluate_leveraging_leg(
            moon=moon,
            n_moon_revs=n,
            m_sc_revs=m,
            vinf_in_kms=vinf_in,
            vinf_out_target_kms=vinf_out,
            exterior=(vinf_out > vinf_in),
        )
        if (
            leg.converged
            and leg.gamma_floor_ok
            and (best is None or leg.dv_dsm_kms < best.dv_dsm_kms)
        ):
            best = leg
    return best


def _vinf_targets(moon: str, moons: Sequence[str], floor: float) -> set[float]:
    """V∞ values worth steering to at ``moon``: the floor + each transfer's depart V∞."""
    targets: set[float] = {floor}
    for other in moons:
        if other != moon:
            targets.add(evaluate_intermoon_transfer(moon, other).vinf_depart_kms)
    return targets


def solve_endgame(
    *,
    moon_system: str,
    entry_moon: str,
    vinf_entry_kms: float,
    target_vinf_floor_kms: float,
    dv_budget_kms: float,
    target_moon: str | None = None,
    system_moons: Sequence[str] | None = None,
) -> EndgameRoute | None:
    """Cheapest leg+transfer chain lowering V∞ to the floor at ``target_moon``.

    Dijkstra over (moon, V∞) states (edges: leveraging legs within a moon +
    ballistic intermoon transfers). Returns ``None`` if no route reaches the floor
    within ``dv_budget_kms`` (-> a method-versioned EMPTY-region report upstream).
    """
    target_moon = target_moon or entry_moon
    moons = _system_moons(moon_system, system_moons)
    lower_bound = route_lower_bound_kms(entry_moon, target_moon)
    counter = itertools.count()

    pq: list[tuple[float, int, str, float, tuple[object, ...]]] = [
        (0.0, next(counter), entry_moon, vinf_entry_kms, ())
    ]
    best: EndgameRoute | None = None
    seen: dict[tuple[str, int], float] = {}

    while pq:
        cost, _c, moon, vinf, steps = heapq.heappop(pq)
        if cost > dv_budget_kms + 1e-12:
            continue
        if moon == target_moon and vinf <= target_vinf_floor_kms + 1e-9:
            if best is None or cost < best.total_dv_kms:
                best = EndgameRoute(
                    steps=steps,
                    total_dv_kms=cost,
                    vinf_entry_kms=vinf_entry_kms,
                    vinf_final_kms=vinf,
                    lower_bound_kms=lower_bound,
                )
            continue
        key = (moon, round(vinf / _VINF_BUCKET_KMS))
        if key in seen and seen[key] <= cost:
            continue
        seen[key] = cost

        # Leveraging legs: steer V∞ to each useful target at this moon.
        for vt in _vinf_targets(moon, moons, target_vinf_floor_kms):
            if abs(vt - vinf) < _VINF_BUCKET_KMS:
                continue
            leg = _best_leg(moon, vinf, vt)
            if leg is None:
                continue
            ncost = cost + leg.dv_dsm_kms
            if ncost > dv_budget_kms + 1e-12:
                continue
            heapq.heappush(pq, (ncost, next(counter), moon, leg.vinf_out_kms, (*steps, leg)))

        # Intermoon transfers: fire when V∞ matches a transfer's departure value.
        for other in moons:
            if other == moon:
                continue
            t = evaluate_intermoon_transfer(moon, other)
            if abs(vinf - t.vinf_depart_kms) < _VINF_MATCH_TOL_KMS:
                ncost = cost + t.dv_kms
                if ncost > dv_budget_kms + 1e-12:
                    continue
                heapq.heappush(pq, (ncost, next(counter), other, t.vinf_arrive_kms, (*steps, t)))
    return best


def _brute_force_optimum(
    *,
    moon_system: str,
    entry_moon: str,
    vinf_entry_kms: float,
    target_vinf_floor_kms: float,
    dv_budget_kms: float,
    target_moon: str | None = None,
    system_moons: Sequence[str] | None = None,
) -> float | None:
    """Exhaustive min total ΔV over the same edge set (test oracle)."""
    target_moon = target_moon or entry_moon
    moons = _system_moons(moon_system, system_moons)
    best = [math.inf]
    seen: dict[tuple[str, int], float] = {}

    def rec(moon: str, vinf: float, cost: float) -> None:
        if cost > dv_budget_kms + 1e-12:
            return
        if moon == target_moon and vinf <= target_vinf_floor_kms + 1e-9:
            best[0] = min(best[0], cost)
            return
        key = (moon, round(vinf / _VINF_BUCKET_KMS))
        if key in seen and seen[key] <= cost:
            return
        seen[key] = cost
        for vt in _vinf_targets(moon, moons, target_vinf_floor_kms):
            if abs(vt - vinf) < _VINF_BUCKET_KMS:
                continue
            leg = _best_leg(moon, vinf, vt)
            if leg is not None:
                rec(moon, leg.vinf_out_kms, cost + leg.dv_dsm_kms)
        for other in moons:
            if other == moon:
                continue
            t = evaluate_intermoon_transfer(moon, other)
            if abs(vinf - t.vinf_depart_kms) < _VINF_MATCH_TOL_KMS:
                rec(other, t.vinf_arrive_kms, cost + t.dv_kms)

    rec(entry_moon, vinf_entry_kms, 0.0)
    return None if math.isinf(best[0]) else best[0]
