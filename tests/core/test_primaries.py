"""Tier-1 Phase 1: PRIMARIES mu table is present + sourced (plan Phase 1)."""

from __future__ import annotations

from cyclerfinder.core.satellites import PRIMARIES


def test_jupiter_and_saturn_primaries_present() -> None:
    assert "Jupiter" in PRIMARIES
    assert "Saturn" in PRIMARIES
    # Jovicentric GM ~1.2669e8 km^3/s^2 (JPL SSD planetary GM table).
    assert 1.26e8 < PRIMARIES["Jupiter"] < 1.27e8
    # Saturnian GM ~3.7931e7 km^3/s^2.
    assert 3.79e7 < PRIMARIES["Saturn"] < 3.80e7
