# tests/search/test_endgame_graph.py
"""Multi-moon endgame route search (plan 2026-06-09, Component 2)."""

from __future__ import annotations

import pytest

from cyclerfinder.search import endgame_graph as eg
from cyclerfinder.search import leveraging_leg as ll  # noqa: F401 (used in Task 5+)
from cyclerfinder.search import vilm


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
