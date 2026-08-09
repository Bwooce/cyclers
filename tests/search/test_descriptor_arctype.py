"""M-ED Phase 2: arc_type -> leg topology (plan Phase 2; spec §16.7.7)."""

from __future__ import annotations

from cyclerfinder.search.descriptor import arc_to_leg_topology


def test_generic_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("generic", resonance=None) == (0, "single")


def test_half_rev_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("half-rev", resonance=None) == (0, "single")


def test_full_rev_arc_uses_resonance_revs() -> None:
    # "3:2" -> 3 Earth years, 2 spacecraft revs (McConaghy/Russell/Longuski
    # 2005, JSR 42(4): "The second parameter, N, is the number of spacecraft
    # revolutions"; Russell 2004 p.126: "the number following the colon ...
    # represent[s] the number of revolutions by the spacecraft"). #794 fixed
    # the previously-reversed reading.
    assert arc_to_leg_topology("full-rev", resonance="3:2") == (2, "low")


def test_full_rev_one_to_one() -> None:
    assert arc_to_leg_topology("full-rev", resonance="1:1") == (1, "low")
