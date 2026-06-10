# tests/search/test_endgame_graph.py
"""Multi-moon endgame route search (plan 2026-06-09, Component 2)."""

from __future__ import annotations

from cyclerfinder.search import endgame_graph as eg
from cyclerfinder.search import vilm


def test_lower_bound_is_admissible_vs_vilm() -> None:
    bound = eg.route_lower_bound_kms("Ganymede", "Europa")
    full = vilm.vilm_dv_min("Ganymede", "Europa")
    assert bound <= full + 1e-9
    assert bound > 0.0
