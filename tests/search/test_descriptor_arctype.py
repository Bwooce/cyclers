"""M-ED Phase 2: arc_type -> leg topology (plan Phase 2; spec §16.7.7)."""

from __future__ import annotations

from cyclerfinder.search.descriptor import arc_to_leg_topology


def test_generic_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("generic", resonance=None) == (0, "single")


def test_half_rev_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("half-rev", resonance=None) == (0, "single")


def test_full_rev_arc_uses_resonance_revs() -> None:
    # "3:2" -> spacecraft does 3 revs (M:N, M = spacecraft revs, spec §16.7.7).
    assert arc_to_leg_topology("full-rev", resonance="3:2") == (3, "low")
