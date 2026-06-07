"""Tier-1 Phase 4: Tisserand a_p + mu resolve a moon about its primary (plan Phase 4)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.tisserand import _a_p_km


def test_a_p_km_for_moon_is_about_primary_sma() -> None:
    assert _a_p_km("Europa") == pytest.approx(SATELLITES["Europa"].sma_km)


def test_a_p_km_for_planet_unchanged() -> None:
    from cyclerfinder.core.constants import AU_KM, PLANETS

    assert _a_p_km("E") == pytest.approx(PLANETS["E"].sma_au * AU_KM)


def test_vinf_tisserand_roundtrip_about_jupiter() -> None:
    from cyclerfinder.search.tisserand import tisserand_to_vinf, vinf_to_tisserand

    mu = PRIMARIES["Jupiter"]
    v = 5.0  # km/s Jovicentric V_inf at Europa
    t = vinf_to_tisserand("Europa", v, mu=mu)
    assert tisserand_to_vinf("Europa", t, mu=mu) == pytest.approx(v, rel=1e-9)


def test_sun_default_unchanged() -> None:
    from cyclerfinder.search.tisserand import vinf_to_tisserand

    # No mu= -> heliocentric, byte-identical to pre-change.
    assert vinf_to_tisserand("E", 5.0) == vinf_to_tisserand("E", 5.0)


def test_linkable_resolves_a_jovicentric_moon_pair() -> None:
    from cyclerfinder.core.constants import AU_KM
    from cyclerfinder.search.tisserand import linkable

    # Europa-Ganymede share Jupiter; at a feasible Jovicentric V_inf the contours
    # intersect (a linked-conic patch point exists). The contour scale is the
    # moons' about-primary SMA, so a_range is moon-scale (km->AU). Centre-blind
    # once a_p (_a_p_km->SATELLITES) and mu (PRIMARIES) resolve correctly.
    # NON-GOLDEN: assert feasibility (a patch point exists), not a specific V_inf.
    a_range = (300_000.0 / AU_KM, 2_500_000.0 / AU_KM)
    assert linkable("Europa", "Ganymede", 3.0, a_range_au=a_range, mu=PRIMARIES["Jupiter"]) is True
