"""Tier-1 Phase 3: corrector is centre-agnostic — mu_central plumbed into Lambert
(plan Phase 3 Task 3.1). Heliocentric default stays byte-identical."""

from __future__ import annotations

import inspect

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.search.correct import _vinf_nodes, ballistic_correct


def test_vinf_nodes_accepts_mu_central_defaulting_to_sun() -> None:
    sig = inspect.signature(_vinf_nodes)
    assert "mu_central" in sig.parameters
    assert sig.parameters["mu_central"].default == MU_SUN_KM3_S2


def test_ballistic_correct_accepts_mu_central() -> None:
    sig = inspect.signature(ballistic_correct)
    assert "mu_central" in sig.parameters
    assert sig.parameters["mu_central"].default == MU_SUN_KM3_S2


def test_max_bend_resolves_a_moon_code() -> None:
    from cyclerfinder.search.correct import _max_bend_deg

    # Europa is in SATELLITES, not PLANETS; the bend lookup must resolve it.
    bend = _max_bend_deg(3.0, "Europa")
    assert bend > 0.0
    # Higher V_inf -> tighter turn at the same moon.
    assert _max_bend_deg(8.0, "Europa") < bend
