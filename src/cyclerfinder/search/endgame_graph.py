# src/cyclerfinder/search/endgame_graph.py
"""Multi-moon endgame route search (spec 2026-06-09, Component 2).

Chains phase-full VILM legs (:mod:`cyclerfinder.search.leveraging_leg`) and
ballistic intermoon transfers to walk a moon-tour cycler's encounter V∞ from a
high entry down to a bend-feasible target floor. Dijkstra over (moon, V∞) states
(non-negative edge costs -> optimal); the phase-free Γ-quadrature
(:mod:`cyclerfinder.search.vilm`) supplies the admissible lower bound. Pure.
"""

from __future__ import annotations

from cyclerfinder.search import vilm


def route_lower_bound_kms(entry_moon: str, target_moon: str) -> float:
    """Admissible ΔV lower bound (km/s): vilm escape+capture insertion floor."""
    return vilm.vilm_dv_floor(entry_moon, target_moon)
