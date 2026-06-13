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


def test_outer_and_mars_primaries_present_and_sourced() -> None:
    """Mars/Uranus/Neptune/Pluto system GMs (JPL DE440), added 2026-06-14."""
    for p in ("Mars", "Uranus", "Neptune", "Pluto"):
        assert p in PRIMARIES
    # Mars system GM ~4.2828e4 km^3/s^2 (JPL DE440).
    assert 4.28e4 < PRIMARIES["Mars"] < 4.29e4
    # Uranus system GM ~5.7946e6 km^3/s^2 (JPL DE440).
    assert 5.79e6 < PRIMARIES["Uranus"] < 5.80e6
    # Neptune system GM ~6.8365e6 km^3/s^2 (JPL DE440).
    assert 6.83e6 < PRIMARIES["Neptune"] < 6.84e6
    # Pluto system GM ~975.5 km^3/s^2 (JPL DE440 / Brozovic et al. 2015).
    assert 970.0 < PRIMARIES["Pluto"] < 980.0
