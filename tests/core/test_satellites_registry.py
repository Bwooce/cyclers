"""Tier-1 Phase 1: SATELLITES registry coverage + internal consistency (plan Phase 1)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES, mean_motion_deg_day_about

_GALILEAN = ("Io", "Europa", "Ganymede", "Callisto", "Amalthea")
_SATURNIAN = (
    "Mimas",
    "Enceladus",
    "Tethys",
    "Dione",
    "Rhea",
    "Titan",
    "Iapetus",
    "Hyperion",
)
_MARTIAN = ("Phobos", "Deimos")
_URANIAN = ("Miranda", "Ariel", "Umbriel", "Titania", "Oberon")
_NEPTUNIAN = ("Triton", "Proteus")
_PLUTONIAN = ("Charon", "Nix", "Hydra")

_ALL_MOONS = _GALILEAN + _SATURNIAN + _MARTIAN + _URANIAN + _NEPTUNIAN + _PLUTONIAN

# Each system's full expected moon set (what _registry_moons_for(primary) returns).
# Drives the "system enumerates ONLY its own moons" invariant below.
_SYSTEM_MOONS: dict[str, frozenset[str]] = {
    "Earth": frozenset({"Moon"}),
    "Mars": frozenset(_MARTIAN),
    "Jupiter": frozenset(_GALILEAN),
    "Saturn": frozenset(_SATURNIAN),
    "Uranus": frozenset(_URANIAN),
    "Neptune": frozenset(_NEPTUNIAN),
    "Pluto": frozenset(_PLUTONIAN),
}


@pytest.mark.parametrize("moon", _ALL_MOONS)
def test_moon_present_and_keyed_by_full_name(moon: str) -> None:
    assert moon in SATELLITES
    assert SATELLITES[moon].code == moon  # full-name scheme
    assert SATELLITES[moon].primary in PRIMARIES


@pytest.mark.parametrize("moon", _ALL_MOONS)
def test_mean_motion_is_derived_not_handcopied(moon: str) -> None:
    s = SATELLITES[moon]
    expected = mean_motion_deg_day_about(s.sma_km, mu_primary=PRIMARIES[s.primary])
    assert s.mean_motion_deg_day == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("primary,expected", sorted(_SYSTEM_MOONS.items()))
def test_system_enumerates_only_its_own_moons(primary: str, expected: frozenset[str]) -> None:
    """Each primary's registry slice is EXACTLY its own moons — no cross-leak.

    The registry is filtered by ``sat.primary``, the same predicate
    ``discovery_campaign._registry_moons_for`` uses (so a Uranus campaign sees
    only Uranian moons, never Jovian ones, etc.).
    """
    got = {name for name, sat in SATELLITES.items() if sat.primary == primary}
    assert got == set(expected)
    # And the slice is disjoint from every OTHER system's moon set.
    for other, others_moons in _SYSTEM_MOONS.items():
        if other != primary:
            assert got.isdisjoint(others_moons)


def test_every_satellite_belongs_to_a_known_system() -> None:
    """No moon orbits a primary missing from PRIMARIES or the system map."""
    for name, sat in SATELLITES.items():
        assert sat.primary in PRIMARIES, f"{name} primary {sat.primary} missing from PRIMARIES"
        assert sat.primary in _SYSTEM_MOONS, f"{name} primary {sat.primary} unmapped"
        assert name in _SYSTEM_MOONS[sat.primary]
