"""Tier-1 Phase 1: SATELLITES registry coverage + internal consistency (plan Phase 1)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES, mean_motion_deg_day_about

_GALILEAN = ("Io", "Europa", "Ganymede", "Callisto")
_SATURNIAN = ("Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan")


@pytest.mark.parametrize("moon", _GALILEAN + _SATURNIAN)
def test_moon_present_and_keyed_by_full_name(moon: str) -> None:
    assert moon in SATELLITES
    assert SATELLITES[moon].code == moon  # full-name scheme
    assert SATELLITES[moon].primary in PRIMARIES


@pytest.mark.parametrize("moon", _GALILEAN + _SATURNIAN)
def test_mean_motion_is_derived_not_handcopied(moon: str) -> None:
    s = SATELLITES[moon]
    expected = mean_motion_deg_day_about(s.sma_km, mu_primary=PRIMARIES[s.primary])
    assert s.mean_motion_deg_day == pytest.approx(expected, rel=1e-12)
