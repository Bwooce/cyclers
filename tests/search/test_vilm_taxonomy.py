"""Tier-1 Phase 5: VILM n:m_K± leg taxonomy (plan Phase 5 Task 5.0).
Physics-invariant classification; NOT gated by a paper cell."""

from __future__ import annotations

from cyclerfinder.search.vilm import classify_vilm_leg


def test_exterior_vilm_classification() -> None:
    leg = classify_vilm_leg(n=3, m=2, body="Europa", exterior=True)
    assert leg.resonance == (3, 2)
    assert leg.body == "Europa"
    assert leg.exterior is True
