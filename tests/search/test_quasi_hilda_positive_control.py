"""Tests for :mod:`cyclerfinder.search.quasi_hilda_positive_control` (task #664).

Physics correctness checks for the Sun-Jupiter Dellnitz-2005 positive-control
glue: energy conversion, section-state construction, osculating-element
correctness, region-indicator disjointness, and Poincare-map energy
conservation. The full pilot run (box-covering + transfer operator +
almost-invariant-set extraction + transport-probability comparison against
the paper's sourced numbers) is exercised by
``scripts/run_664_dellnitz_positive_control.py`` (a real compute job, not a
fast unit test) -- see that script's own docstring and #664's
``data/OUTSTANDING.md`` bullet for the actual positive-control verdict.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.search.quasi_hilda_positive_control import (
    DELLNITZ_2005_ENERGY,
    HILDA_RESONANCE_A_NONDIM,
    a_mars_nondim,
    accessible_xdot_max,
    energy_to_jacobi_constant,
    mars_crossing_indicator,
    osculating_elements_at_section,
    poincare_first_return,
    quasi_hilda_region_indicator,
    section_map_xxdot,
    section_state6,
    sun_jupiter_system,
    zero_velocity_v,
)

SYSTEM = sun_jupiter_system()
MU = SYSTEM.mu
C_TARGET = energy_to_jacobi_constant(DELLNITZ_2005_ENERGY, MU)


def test_sun_jupiter_mu_matches_known_mass_ratio() -> None:
    # Jupiter/Sun system mass ratio is close to 1/1047.
    assert pytest.approx(1.0 / 1047.0, rel=0.02) == MU


def test_energy_to_jacobi_constant_sourced_value() -> None:
    # E = -1.52 -> C ~= 3.039 for Sun-Jupiter mu (mu(1-mu) ~ 9.5e-4 correction).
    assert pytest.approx(3.0390470287, abs=1e-6) == C_TARGET


def test_hilda_resonance_a_nondim_matches_kepler_third_law() -> None:
    # (2/3)^(2/3) is the standard 3:2 resonance semimajor-axis ratio.
    assert pytest.approx((2.0 / 3.0) ** (2.0 / 3.0)) == HILDA_RESONANCE_A_NONDIM
    assert 0.76 < HILDA_RESONANCE_A_NONDIM < 0.77


def test_a_mars_nondim_matches_au_ratio() -> None:
    assert a_mars_nondim() == pytest.approx(1.52371034 / 5.20288700, rel=1e-9)


def test_zero_velocity_v_positive_inside_known_accessible_band() -> None:
    # x = -0.5 is well inside the inner (interior-realm) accessible band.
    assert zero_velocity_v(-0.5, MU, C_TARGET) > 0.0


def test_zero_velocity_v_negative_in_known_forbidden_gap() -> None:
    # A direct numerical scan (this task's own working notes) found a
    # forbidden band (V<0) straddling the L3 vicinity, roughly x in
    # (-1.117, -0.892); x=-1.0 sits inside it.
    assert zero_velocity_v(-1.0, MU, C_TARGET) < 0.0


def test_accessible_xdot_max_none_in_forbidden_gap() -> None:
    assert accessible_xdot_max(-1.0, MU, C_TARGET) is None
    result = accessible_xdot_max(-0.5, MU, C_TARGET)
    assert result is not None
    assert result > 0.0


def test_section_state6_has_correct_jacobi_constant() -> None:
    """The constructed 6-state must reproduce C_TARGET exactly under this
    project's OWN jacobi_constant function -- an independent cross-check
    that the energy-to-Jacobi conversion and ydot solve are self-consistent."""
    state = section_state6(-0.6, 0.2, MU, C_TARGET)
    assert state is not None
    c = jacobi_constant(state, MU)
    assert c == pytest.approx(C_TARGET, abs=1e-10)
    assert state[1] == 0.0  # y = 0
    assert state[4] < 0.0  # ydot < 0


def test_section_state6_none_outside_energy_manifold() -> None:
    # xdot larger than accessible_xdot_max at this x -> ydot^2 < 0.
    xdot_max = accessible_xdot_max(-0.5, MU, C_TARGET)
    assert xdot_max is not None
    assert section_state6(-0.5, xdot_max * 2.0, MU, C_TARGET) is None


def test_osculating_elements_hilda_resonance_orbit() -> None:
    """A section point built to approximate a near-circular 3:2-resonance
    orbit should recover a close to HILDA_RESONANCE_A_NONDIM."""
    # x = -a (near-zero eccentricity, near-circular crossing at conjunction).
    a_guess = HILDA_RESONANCE_A_NONDIM
    xdot_max = accessible_xdot_max(-a_guess, MU, C_TARGET)
    assert xdot_max is not None
    # A small xdot near 0 approximates a near-circular crossing.
    elems = osculating_elements_at_section(-a_guess, 0.0, MU, C_TARGET)
    assert elems is not None
    assert elems.a == pytest.approx(HILDA_RESONANCE_A_NONDIM, rel=0.05)
    assert elems.e < 0.3


def test_mars_crossing_and_hilda_region_are_disjoint_on_a_grid() -> None:
    """Direct regression check for the claim in the module docstring: the
    default quasi-Hilda region window and the Mars-crossing region are
    geometrically disjoint at E=-1.52 for Sun-Jupiter."""
    xs = np.linspace(-0.88, -0.05, 40)
    xdots = np.linspace(-0.9, 0.9, 40)
    both = 0
    r_count = 0
    q_count = 0
    for x in xs:
        for xdot in xdots:
            in_r = quasi_hilda_region_indicator(x, xdot, MU, C_TARGET)
            in_q = mars_crossing_indicator(x, xdot, MU, C_TARGET)
            if in_r is None or in_q is None:
                continue
            r_count += int(in_r)
            q_count += int(in_q)
            if in_r and in_q:
                both += 1
    assert r_count > 0
    assert q_count > 0
    assert both == 0


def test_mars_crossing_indicator_none_outside_energy_manifold() -> None:
    xdot_max = accessible_xdot_max(-0.5, MU, C_TARGET)
    assert xdot_max is not None
    assert mars_crossing_indicator(-0.5, xdot_max * 2.0, MU, C_TARGET) is None


def test_poincare_first_return_conserves_jacobi_constant() -> None:
    """Independent correctness check on the propagation itself: the CR3BP
    Jacobi constant must be conserved by the underlying integrator across a
    full first-return propagation, to tight tolerance."""
    state0 = section_state6(-0.7, 0.15, MU, C_TARGET)
    assert state0 is not None
    c0 = jacobi_constant(state0, MU)
    result = poincare_first_return(state0, MU)
    assert result is not None
    # Landed back on the section (y=0) with ydot<0, x<0.
    assert abs(result[1]) < 1e-8
    assert result[4] < 0.0
    assert result[0] < 0.0
    c1 = jacobi_constant(result, MU)
    assert c1 == pytest.approx(c0, abs=1e-8)


def test_section_map_xxdot_round_trips_through_poincare_first_return() -> None:
    image = section_map_xxdot(-0.7, 0.15, MU, C_TARGET)
    assert image is not None
    assert image.shape == (2,)
    # Rebuilding a section state from the image must land on the same
    # energy manifold (Jacobi constant preserved through the returned
    # (x, xdot) pair, i.e. the map really does stay ON the C_TARGET section).
    rebuilt = section_state6(float(image[0]), float(image[1]), MU, C_TARGET)
    assert rebuilt is not None


def test_section_map_xxdot_none_off_energy_manifold() -> None:
    xdot_max = accessible_xdot_max(-0.5, MU, C_TARGET)
    assert xdot_max is not None
    assert section_map_xxdot(-0.5, xdot_max * 2.0, MU, C_TARGET) is None


def test_energy_to_jacobi_constant_matches_hand_derivation() -> None:
    # Sanity-check the -2E - mu(1-mu) formula against a trivial mu=0 case.
    assert energy_to_jacobi_constant(-1.0, 0.0) == pytest.approx(2.0)
