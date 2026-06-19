"""#391 — tests for the amplitude-vs-Hill-fraction pre-screen.

Frozen gates:

  1. ``branch_C32_b0`` classifies ``V4_DOOMED`` — cross-checked against the
     known #389 V4 verdict Hill fraction of ~0.766 (the orbit that actually
     failed V4 by 4-5 orders into a ~10⁹ km escape).
  2. A small reference orbit (the sourced C11a Braik-Ross parent, ~0.27 Hill)
     classifies ``PASS``.
  3. The sourced Earth-Sun Hill radius matches the #389 V4 verdict value
     (~1.50e6 km) to within the SMA-vs-1-au difference.
  4. The screen's max excursion is >= the IC-only distance (it sees apoapsis).
  5. Input validation raises on malformed state / period.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.genome.hill_screen import (
    CLASS_PASS,
    CLASS_V4_DOOMED,
    earth_sun_hill_radius_km,
    screen_orbit,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

# branch_C32_b0 IC (data/floquet_phase2_family_network.jsonl row branch_C32_b0).
_C32_B0_STATE0 = np.array(
    [
        -0.7033325748418664,
        -2.9123784605081626,
        1.730700297697496e-22,
        -2.3503749595840504,
        0.567571628434372,
        8.295924032385729e-24,
    ],
    dtype=np.float64,
)
_C32_B0_PERIOD_TU = 23.355184434547017

# C11a sourced Braik-Ross Table 2 parent (small-amplitude reference, ~0.27 Hill).
_C11A_PARENT_STATE0 = np.array(
    [-0.7892215761812638, 0.0, 0.0, 0.0, -0.18611872501457488, 0.0], dtype=np.float64
)
_C11A_PARENT_PERIOD_TU = 9.715104513292248

# #389 V4 verdict structural_diagnostic values (data/branch_c32_b0_v4_verdict.jsonl).
_V4_VERDICT_HILL_FRACTION = 0.766437003544452
_V4_VERDICT_HILL_RADIUS_KM = 1502669.4473520927


def test_branch_c32_b0_is_v4_doomed() -> None:
    """branch_C32_b0 (the orbit that actually failed V4) screens V4_DOOMED."""
    system = braik_ross_system()
    result = screen_orbit(system, _C32_B0_STATE0, _C32_B0_PERIOD_TU)
    assert result.classification == CLASS_V4_DOOMED
    # Cross-check the Hill fraction against the #389 V4 verdict (~0.766). The
    # screen uses the true max excursion (vs the verdict's IC-only norm), so it
    # is slightly larger but must agree to ~1%.
    assert result.hill_fraction >= _V4_VERDICT_HILL_FRACTION
    assert abs(result.hill_fraction - _V4_VERDICT_HILL_FRACTION) < 0.01
    # The verdict isolated the solar tide at ~30% of Earth gravity.
    assert abs(result.solar_tide_to_earth_gravity_ratio - 0.30) < 0.02


def test_small_reference_orbit_passes() -> None:
    """A small-amplitude sourced parent (C11a, ~0.27 Hill) screens PASS."""
    system = braik_ross_system()
    result = screen_orbit(system, _C11A_PARENT_STATE0, _C11A_PARENT_PERIOD_TU)
    assert result.classification == CLASS_PASS
    assert result.hill_fraction < 0.3


def test_hill_radius_matches_v4_verdict() -> None:
    """The sourced Earth-Sun Hill radius matches the #389 V4 verdict value."""
    r_hill = earth_sun_hill_radius_km()
    # Verdict used exact 1 au; we use the sourced Earth SMA (1.0000026 au) -> a
    # ~3 km (2e-6 relative) difference. Both are ~1.50e6 km.
    assert abs(r_hill - _V4_VERDICT_HILL_RADIUS_KM) / _V4_VERDICT_HILL_RADIUS_KM < 1e-5


def test_max_excursion_at_least_ic_distance() -> None:
    """The propagated max excursion is >= the IC distance from Earth."""
    system = braik_ross_system()
    result = screen_orbit(system, _C32_B0_STATE0, _C32_B0_PERIOD_TU)
    earth = np.array([-system.mu, 0.0, 0.0], dtype=np.float64)
    ic_dist_km = float(np.linalg.norm(_C32_B0_STATE0[:3] - earth)) * float(system.l_km)
    assert result.max_amplitude_km >= ic_dist_km


def test_screen_rejects_bad_inputs() -> None:
    """Malformed state / period raise ValueError, not a silent default."""
    system = braik_ross_system()
    with pytest.raises(ValueError, match="shape"):
        screen_orbit(system, np.zeros(5), 1.0)
    with pytest.raises(ValueError, match="period"):
        screen_orbit(system, _C32_B0_STATE0, -1.0)
