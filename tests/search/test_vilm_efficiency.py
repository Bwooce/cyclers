"""Tier-1 Phase 5 GOLDEN: VILM efficiency root vs published Part-1 Table 3
V_inf_bar E/I (mining note A1, lines 337-354).

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Table 3 V_inf_bar values (km/s),
exterior (E) column. The min_vinf_for_vilm output is the side under test."""

from __future__ import annotations

import pytest

from cyclerfinder.search.vilm import min_vinf_for_vilm

# (moon, V_inf_bar exterior km/s) — Part-1 Table 3, E column (note lines 339-348).
_VINF_BAR_E = [
    ("Io", 0.351),
    ("Europa", 0.277),
    ("Ganymede", 0.372),
    ("Callisto", 0.328),
    ("Titan", 0.283),
    ("Rhea", 0.085),
    ("Dione", 0.067),
    ("Tethys", 0.052),
    ("Enceladus", 0.029),
]

# (moon, V_inf_bar interior km/s) — Part-1 Table 3, I column (note lines 339-348).
_VINF_BAR_I = [
    ("Io", 0.368),
    ("Europa", 0.290),
    ("Ganymede", 0.404),
    ("Callisto", 0.361),
    ("Titan", 0.321),
    ("Rhea", 0.087),
    ("Dione", 0.068),
    ("Tethys", 0.052),
    ("Enceladus", 0.029),
]


@pytest.mark.parametrize("moon,vbar_e", _VINF_BAR_E)
def test_efficiency_root_matches_published_vinf_bar(moon: str, vbar_e: float) -> None:
    # Paper rounds to 3 dp; band absorbs the V_c assumption + mu digit choice.
    assert min_vinf_for_vilm(moon) == pytest.approx(vbar_e, abs=0.01)


@pytest.mark.parametrize("moon,vbar_i", _VINF_BAR_I)
def test_efficiency_root_interior_matches_published(moon: str, vbar_i: float) -> None:
    assert min_vinf_for_vilm(moon, exterior=False) == pytest.approx(vbar_i, abs=0.01)
