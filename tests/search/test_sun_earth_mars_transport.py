"""Tests for :mod:`cyclerfinder.search.sun_earth_mars_transport` (task #685).

Physics-correctness checks for the Sun-Earth PCR3BP glue that applies the #664
set-oriented transfer-operator (GAIO) pipeline to this project's own
Earth-Mars transport domain: energy self-consistency, section-state
construction, osculating-element correctness, exterior-branch Poincare-map
conservation, and R/Q region disjointness + banding. The full pilot run
(box-covering + transfer operator + almost-invariant-set extraction + R->Q
transport probability) is exercised by
``scripts/run_685_sun_earth_mars_transport.py`` (a real compute job, not a
fast unit test) -- see that script's docstring and #685's ``data/OUTSTANDING.md``
bullet for the actual search verdict.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.search.sun_earth_mars_transport import (
    MARS_REACHING_R_A_MAX,
    SUN_EARTH_MARS_TRANSFER_C,
    a_mars_nondim,
    accessible_xdot_max,
    earth_neighborhood_region_indicator,
    hohmann_transfer_elements,
    mars_reaching_indicator,
    osculating_elements_at_section,
    poincare_first_return_exterior,
    section_map_xxdot,
    section_state6,
    sun_earth_system,
    zero_velocity_v,
)

SYSTEM = sun_earth_system()
MU = SYSTEM.mu
C = SUN_EARTH_MARS_TRANSFER_C


def test_sun_earth_mu_matches_known_mass_ratio() -> None:
    # Earth/(Sun+Earth) mass ratio ~= 3.0e-6.
    assert pytest.approx(3.0e-6, rel=0.05) == MU


def test_a_mars_nondim_matches_au_ratio() -> None:
    assert a_mars_nondim() == pytest.approx(1.52371034 / 1.00000261, rel=1e-9)
    # Earth-SMA length unit -> Mars sits at ~1.5237 length units.
    assert 1.52 < a_mars_nondim() < 1.525


def test_hohmann_transfer_elements() -> None:
    a, e = hohmann_transfer_elements()
    # Perihelion at Earth (1.0), aphelion at Mars (a_Mars).
    assert a * (1.0 - e) == pytest.approx(1.0, rel=1e-9)
    assert a * (1.0 + e) == pytest.approx(a_mars_nondim(), rel=1e-9)


def test_transfer_jacobi_constant_sourced_value() -> None:
    # Tisserand C of the Earth-Mars Hohmann ellipse ~= 2.9902.
    assert pytest.approx(2.990225, abs=1e-4) == C


def test_transfer_energy_below_c_l2_open_neck() -> None:
    """The chosen energy must sit BELOW the Sun-Earth L2 Jacobi constant (open
    L2 neck -- the transport regime this task deliberately works in)."""
    from scipy.optimize import brentq

    def d_omega(x: float) -> float:
        r1 = abs(x + MU)
        r2 = abs(x - 1.0 + MU)
        return x - (1.0 - MU) * (x + MU) / r1**3 - MU * (x - 1.0 + MU) / r2**3

    x_l2 = brentq(d_omega, 1.0 - MU + 1e-4, 1.1)
    c_l2 = x_l2 * x_l2 + 2.0 * (1.0 - MU) / abs(x_l2 + MU) + 2.0 * MU / abs(x_l2 - 1.0 + MU)
    assert c_l2 > C
    assert c_l2 == pytest.approx(3.0009, abs=1e-3)


def test_zero_velocity_v_positive_inside_domain() -> None:
    # x = 1.3 is well inside the accessible exterior band at this energy.
    assert zero_velocity_v(1.3, MU, C) > 0.0


def test_accessible_xdot_max_positive_on_domain() -> None:
    # At this open-neck energy the x^2 term keeps V(x) > 0 everywhere on the
    # positive-x exterior domain, so accessibility is bounded purely by
    # |xdot| <= xdot_max (a positive float), not by a forbidden x band.
    for x in (1.03, 1.3, 1.6):
        result = accessible_xdot_max(x, MU, C)
        assert result is not None
        assert result > 0.0
    # Off-manifold (xdot > xdot_max) is exercised by the section_state6 /
    # indicator None-branch tests below.
    edge = accessible_xdot_max(1.3, MU, C)
    assert edge is not None
    assert section_state6(1.3, edge * 1.001, MU, C) is None


def test_section_state6_has_correct_jacobi_constant() -> None:
    """The constructed 6-state must reproduce C exactly under this project's OWN
    jacobi_constant -- an independent cross-check that the section-state ydot
    solve is self-consistent with the chosen energy shell."""
    state = section_state6(1.3, 0.2, MU, C)
    assert state is not None
    assert jacobi_constant(state, MU) == pytest.approx(C, abs=1e-12)
    assert state[0] > 0.0  # x > 0 (exterior branch)
    assert state[1] == 0.0  # y = 0
    assert state[4] < 0.0  # ydot < 0


def test_section_state6_none_outside_energy_manifold() -> None:
    xdot_max = accessible_xdot_max(1.3, MU, C)
    assert xdot_max is not None
    assert section_state6(1.3, xdot_max * 2.0, MU, C) is None


def test_osculating_elements_apsis_crossing_hand_check() -> None:
    """A section crossing with xdot = 0 is an apsis (radial velocity zero), so
    the heliocentric distance there equals an osculating apside; at x = 1.1 the
    crossing distance IS the aphelion of a small (a ~ 1) ellipse."""
    elems = osculating_elements_at_section(1.1, 0.0, MU, C)
    assert elems is not None
    assert elems.r_a == pytest.approx(1.1, abs=1e-3)
    assert 0.9 < elems.a < 1.1
    assert 0.0 <= elems.e < 0.3


def test_mars_reaching_region_is_banded_not_a_half_plane() -> None:
    """A near-parabolic orbit whose aphelion runs far past the Mars annulus
    must NOT count as Mars-reaching (the banding, not an r_a >= a_Mars
    half-plane)."""
    # Find a section point on-manifold with a huge aphelion.
    found_huge = False
    for x in np.linspace(1.05, 1.58, 40):
        xdot_max = accessible_xdot_max(x, MU, C)
        if xdot_max is None:
            continue
        elems = osculating_elements_at_section(x, 0.95 * xdot_max, MU, C)
        if elems is not None and elems.e < 1.0 and elems.r_a > MARS_REACHING_R_A_MAX:
            found_huge = True
            assert mars_reaching_indicator(x, 0.95 * xdot_max, MU, C) is False
    assert found_huge, "expected at least one on-manifold near-parabolic (r_a > cap) sample"


def test_earth_and_mars_regions_disjoint_and_nonempty_on_grid() -> None:
    """Regression check on the module's own disjointness claim across the pilot
    domain."""
    xs = np.linspace(1.03, 1.60, 45)
    xdots = np.linspace(-0.85, 0.85, 45)
    r_count = q_count = both = 0
    for x in xs:
        for xdot in xdots:
            in_r = earth_neighborhood_region_indicator(x, xdot, MU, C)
            in_q = mars_reaching_indicator(x, xdot, MU, C)
            if in_r is None or in_q is None:
                continue
            r_count += int(in_r)
            q_count += int(in_q)
            if in_r and in_q:
                both += 1
    assert r_count > 0
    assert q_count > 0
    assert both == 0


def test_mars_reaching_indicator_none_off_manifold() -> None:
    xdot_max = accessible_xdot_max(1.3, MU, C)
    assert xdot_max is not None
    assert mars_reaching_indicator(1.3, xdot_max * 2.0, MU, C) is None


def test_poincare_first_return_exterior_conserves_jacobi_and_branch() -> None:
    """Independent correctness check on the propagation: the CR3BP Jacobi
    constant is conserved across a full exterior first-return, and the return
    lands on the y=0/ydot<0/x>0 section branch."""
    state0 = section_state6(1.35, 0.15, MU, C)
    assert state0 is not None
    c0 = jacobi_constant(state0, MU)
    result = poincare_first_return_exterior(state0, MU)
    assert result is not None
    assert abs(result[1]) < 1e-7  # y = 0
    assert result[4] < 0.0  # ydot < 0
    assert result[0] > 0.0  # x > 0 (exterior branch)
    assert jacobi_constant(result, MU) == pytest.approx(c0, abs=1e-8)


def test_section_map_xxdot_round_trips_on_manifold() -> None:
    image = section_map_xxdot(1.35, 0.15, MU, C)
    assert image is not None
    assert image.shape == (2,)
    rebuilt = section_state6(float(image[0]), float(image[1]), MU, C)
    assert rebuilt is not None
    assert jacobi_constant(rebuilt, MU) == pytest.approx(C, abs=1e-10)


def test_section_map_xxdot_none_off_manifold() -> None:
    xdot_max = accessible_xdot_max(1.3, MU, C)
    assert xdot_max is not None
    assert section_map_xxdot(1.3, xdot_max * 2.0, MU, C) is None


def test_kickoff_prevents_t0_self_detection() -> None:
    """The #664 t=0 self-detection gotcha regression: a state exactly on the
    section must NOT be returned unchanged as its own first return -- the
    exterior map must genuinely advance it."""
    x0, xdot0 = 1.4, 0.1
    image = section_map_xxdot(x0, xdot0, MU, C)
    assert image is not None
    # The return must differ from the start (a real first return, not the
    # spuriously self-detected t=0 point).
    assert not (
        math.isclose(image[0], x0, abs_tol=1e-9) and math.isclose(image[1], xdot0, abs_tol=1e-9)
    )
