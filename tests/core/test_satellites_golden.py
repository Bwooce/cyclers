"""Tier-1 Phase 1 GOLDEN: SATELLITES reproduces the published Endgame Part-1
Table 3 a_M / V_M (mining note A1, lines 337-352).

GOLDEN DISCIPLINE: EXPECTED = the PUBLISHED Campagnola & Russell Part-1 Table 3
values, sourced INDEPENDENTLY of the JPL-SSD values that built the registry. The
circular velocity is V_M = sqrt(mu_primary / a_M). Tolerance bands documented
inline; the paper rounds a_M to km and V_M to 3 dp.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# (moon, a_M km, V_M km/s) — Campagnola & Russell, "Endgame Problem" Part 1,
# Table 3 (p.17), transcribed at
# docs/notes/2026-06-05-endgame-tisserand-mining.md:339-348.
_TABLE3 = [
    ("Io", 421800.0, 17.330),
    ("Europa", 671100.0, 13.739),
    ("Ganymede", 1070400.0, 10.879),
    ("Callisto", 1882700.0, 8.203),
    ("Enceladus", 238040.0, 12.624),
    ("Tethys", 294670.0, 11.346),
    ("Dione", 377420.0, 10.025),
    ("Rhea", 527070.0, 8.484),
    ("Titan", 1221870.0, 5.572),
]


@pytest.mark.parametrize("moon,a_km,v_kms", _TABLE3)
def test_registry_reproduces_published_table3(moon: str, a_km: float, v_kms: float) -> None:
    s = SATELLITES[moon]
    # a_M: registry (JPL SSD) vs paper (JPL SSD) — same upstream, expect <0.1%.
    assert s.sma_km == pytest.approx(a_km, rel=1e-3)
    v_circ = math.sqrt(PRIMARIES[s.primary] / s.sma_km)
    # V_M: circular speed at a_M; paper rounds to 3 dp -> 0.5% band absorbs the
    # mu_primary digit choice between JPL releases.
    assert v_circ == pytest.approx(v_kms, rel=5e-3)
