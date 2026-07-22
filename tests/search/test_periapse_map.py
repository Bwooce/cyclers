"""Physics-correctness tests for the #683 periapse Poincaré-map machinery
(``cyclerfinder.search.periapse_map``), the Davis-Howell 2011 lineage map.

These are the calibration checks that must hold before any lobe structure read
off the map is trusted: exact (not sampling-based) periapsis detection, exact
Jacobi-constant parametrisation, Jacobi conservation across the map, and
geometrically-sane fate classification. Cross-checks against Davis & Howell
2011's own published Saturn-Titan / Sun-Saturn search energies are included.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.periapse_map as pm

# Davis & Howell 2011 (Acta Astronautica 69, 1038-1049) published search energies.
J1_SUN_SATURN = 3.0173046596239
J2_SATURN_TITAN = 3.015311017945150


def test_l_points_bracket_secondary_and_match_paper_regime() -> None:
    """C_L1 > C_L2, both bracket the secondary, and the paper's published
    search energies J1 / J2 each sit just below C_L2 for their system (both
    gateways open) -- an independent cross-check of the paper's stated regime
    against this project's own jacobi_constant."""
    ss = pm.build_periapse_map_system("Sun", "Saturn")
    st = pm.build_periapse_map_system("Saturn", "Titan")
    for sysm, jval in [(ss, J1_SUN_SATURN), (st, J2_SATURN_TITAN)]:
        assert sysm.x_l1 < sysm.x_secondary < sysm.x_l2
        assert sysm.c_l1 > sysm.c_l2  # interior point is higher-energy
        # search energy just below C_L2 => both L1 and L2 necks open
        assert jval < sysm.c_l2 < sysm.c_l1
        assert 0.0 < sysm.c_l2 - jval < 5e-4  # "just below" (narrow open neck)


def test_construct_periapse_state_is_exact() -> None:
    """A constructed periapse state reproduces the target Jacobi constant to
    machine precision and sits exactly on the section (rdot == 0)."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    x = (1 - mu) + 0.3 * st.r_hill
    y = 0.2 * st.r_hill
    s = pm.construct_periapse_state(x, y, J2_SATURN_TITAN, mu)
    assert s is not None
    assert abs(cr3bp.jacobi_constant(s, mu) - J2_SATURN_TITAN) < 1e-13
    assert abs(pm.secondary_radial_rate(s, mu)) < 1e-13


def test_construct_outside_zero_velocity_curve_returns_none() -> None:
    """A point where the Jacobi constraint needs v^2 < 0 is off the energy
    manifold and yields no state."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    # A point in the forbidden annulus (between Titan's deep well and the
    # far-field x^2+y^2 rise) requires v^2 < 0 at this energy.
    s = pm.construct_periapse_state((1 - mu), 1.5 * st.r_hill, J2_SATURN_TITAN, mu)
    assert s is None


def test_prograde_periapse_has_positive_angular_momentum_about_secondary() -> None:
    """The prograde construction gives counter-clockwise motion about P2
    (positive z angular momentum of the relative state)."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    s = pm.construct_periapse_state(
        (1 - mu) + 0.25 * st.r_hill, 0.1 * st.r_hill, J2_SATURN_TITAN, mu
    )
    assert s is not None
    dx, dy = float(s[0]) - (1 - mu), float(s[1])
    lz = dx * float(s[4]) - dy * float(s[3])
    assert lz > 0.0


def test_periapsis_detection_exact_and_jacobi_conserved() -> None:
    """The next-periapsis return is a genuine minimum of distance to P2
    (rddot > 0), sits exactly on the section (rdot ~ 0, NOT sampling-based),
    and conserves the Jacobi constant across the map."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    s0 = pm.construct_periapse_state(
        (1 - mu) + 0.2 * st.r_hill, 0.15 * st.r_hill, J2_SATURN_TITAN, mu
    )
    assert s0 is not None
    assert pm.is_periapsis(s0, mu)
    nxt = pm.next_periapse(s0, st)
    assert nxt is not None
    assert abs(pm.secondary_radial_rate(nxt, mu)) < 1e-11
    assert pm.secondary_radial_accel(nxt, mu) > 0.0
    assert abs(cr3bp.jacobi_constant(nxt, mu) - J2_SATURN_TITAN) < 1e-9


def test_periapsis_vs_apoapsis_sign() -> None:
    """secondary_radial_accel distinguishes the two apses that both satisfy
    rdot == 0: positive (radial minimum) at a periapsis, negative (radial
    maximum) at an apoapsis. Both constructed section points sit exactly on
    rdot == 0, so the rddot sign is the only discriminant, as the map's
    periapsis-vs-apoapsis (is_periapsis) test relies on."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    # A deep near-Titan section point is a periapsis (radial minimum).
    peri = pm.construct_periapse_state((1 - mu) + 0.12 * st.r_hill, 0.0, J2_SATURN_TITAN, mu)
    assert peri is not None
    assert abs(pm.secondary_radial_rate(peri, mu)) < 1e-13
    assert pm.secondary_radial_accel(peri, mu) > 0.0
    assert pm.is_periapsis(peri, mu)
    # A section point out near the zero-velocity contour is an apoapsis
    # (radial maximum): same rdot == 0, opposite rddot sign.
    apo = pm.construct_periapse_state((1 - mu), 0.5 * st.r_hill, J2_SATURN_TITAN, mu)
    assert apo is not None
    assert abs(pm.secondary_radial_rate(apo, mu)) < 1e-13
    assert pm.secondary_radial_accel(apo, mu) < 0.0
    assert not pm.is_periapsis(apo, mu)


def test_map_coords_roundtrip() -> None:
    """periapse_map_coords returns the periapse position relative to P2 in
    Hill-radius units, inverting the construction offset."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    xp_in, yp_in = 0.37, -0.21
    s = pm.construct_periapse_state(
        (1 - mu) + xp_in * st.r_hill, yp_in * st.r_hill, J2_SATURN_TITAN, mu
    )
    assert s is not None
    xp, yp, r_p, omega = pm.periapse_map_coords(s, st)
    assert xp == pytest.approx(xp_in, abs=1e-12)
    assert yp == pytest.approx(yp_in, abs=1e-12)
    assert r_p == pytest.approx(math.hypot(xp_in, yp_in) * st.r_hill, rel=1e-12)
    assert omega == pytest.approx(math.atan2(yp_in, xp_in), abs=1e-12)


def test_grazing_periapsis_classified_as_impact() -> None:
    """A periapse constructed inside the secondary's physical radius impacts."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    s = pm.construct_periapse_state(
        (1 - mu) + 0.5 * st.secondary_radius_nd, 0.0, J2_SATURN_TITAN, mu
    )
    assert s is not None
    fate, _ = pm.classify_fate(s, st, max_revs=6)
    assert fate is pm.PeriapseFate.IMPACT


def test_escape_classes_are_directional() -> None:
    """Over a coarse grid, escape-L1 periapses sit on the L1 (interior, x_p<0)
    side and escape-L2 periapses on the L2 (exterior, x_p>0) side."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    l1_xps: list[float] = []
    l2_xps: list[float] = []
    for xp in np.linspace(-1.1, 1.1, 22):
        for yp in np.linspace(-0.8, 0.8, 16):
            s = pm.construct_periapse_state(
                (1 - mu) + xp * st.r_hill, yp * st.r_hill, J2_SATURN_TITAN, mu
            )
            if s is None or not pm.is_periapsis(s, mu):
                continue
            fate, _ = pm.classify_fate(s, st, max_revs=1)
            if fate is pm.PeriapseFate.ESCAPE_L1:
                l1_xps.append(xp)
            elif fate is pm.PeriapseFate.ESCAPE_L2:
                l2_xps.append(xp)
    assert l1_xps and l2_xps
    assert np.mean(l1_xps) < -0.3  # L1 lobe is on the interior side
    assert np.mean(l2_xps) > 0.3  # L2 lobe is on the exterior side


def test_collect_periapses_time_ordered_and_on_section() -> None:
    """collect_titan_periapses returns time-ordered periapses, each exactly on
    the section (rdot ~ 0) and a genuine minimum (rddot > 0)."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    s0 = pm.construct_periapse_state(
        (1 - mu) + 0.2 * st.r_hill, 0.1 * st.r_hill, J2_SATURN_TITAN, mu
    )
    assert s0 is not None
    events, _impacted = pm.collect_titan_periapses(s0, st, t_total=8.0 * math.pi)
    assert len(events) >= 2
    times = [t for t, _ in events]
    assert times == sorted(times)
    for _t, state in events:
        assert abs(pm.secondary_radial_rate(state, mu)) < 1e-8
        assert pm.secondary_radial_accel(state, mu) > 0.0
        assert abs(cr3bp.jacobi_constant(state, mu) - J2_SATURN_TITAN) < 1e-8


def test_impact_flag_set_for_collision_seed() -> None:
    """A near-collision periapse seed terminates collect_titan_periapses with
    the impact flag set."""
    st = pm.build_periapse_map_system("Saturn", "Titan")
    mu = st.mu
    s0 = pm.construct_periapse_state(
        (1 - mu) + 0.5 * st.secondary_radius_nd, 0.0, J2_SATURN_TITAN, mu
    )
    assert s0 is not None
    _events, impacted = pm.collect_titan_periapses(s0, st, t_total=20.0 * math.pi)
    assert impacted
