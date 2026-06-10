# src/cyclerfinder/search/endgame_graph.py
"""Multi-moon endgame route search (spec 2026-06-09, Component 2).

Chains phase-full VILM legs (:mod:`cyclerfinder.search.leveraging_leg`) and
ballistic intermoon transfers to walk a moon-tour cycler's encounter V∞ from a
high entry down to a bend-feasible target floor. Dijkstra over (moon, V∞) states
(non-negative edge costs -> optimal); the phase-free Γ-quadrature
(:mod:`cyclerfinder.search.vilm`) supplies the admissible lower bound. Pure.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search import vilm

DAY_S = 86400.0


def route_lower_bound_kms(entry_moon: str, target_moon: str) -> float:
    """Admissible ΔV lower bound (km/s): vilm escape+capture insertion floor."""
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
