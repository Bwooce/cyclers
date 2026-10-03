"""Tests for the #885 Uranus tour-to-cycler bridge primitives.

Expected values are either published (AAS 25-668 Table 8) or physics identities
(two-body integration, the mirror-closure turn), never numbers this module produced.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cyclerfinder.core.lambert import lambert
from cyclerfinder.search import uranus_bridge_885 as ub


def _two_body(_t: float, y: np.ndarray) -> np.ndarray:
    r = y[:2]
    out: np.ndarray = np.r_[y[2:], -ub.MU_URANUS * r / np.linalg.norm(r) ** 3]
    return out


def test_ellison_two_moon_orbit_predicts_the_other_published_vinf() -> None:
    """AAS 25-668 Table 8: the Oberon 3.448 / Ariel 4.325 conic gives Umbriel and
    Titania V-infinity close to the tour's own Umbriel (4.392-4.396) and Titania
    (3.880-3.922) flybys, i.e. the tour is near one coplanar Tisserand point."""
    a, e = ub.orbit_from_vinf_pair("Oberon", 3.448, "Ariel", 4.325)
    assert ub.vinf_at("Umbriel", a, e) == pytest.approx(4.394, abs=0.03)
    assert ub.vinf_at("Titania", a, e) == pytest.approx(3.90, abs=0.04)
    # and it reproduces its own defining pair
    assert ub.vinf_at("Oberon", a, e) == pytest.approx(3.448, abs=1e-9)
    assert ub.vinf_at("Ariel", a, e) == pytest.approx(4.325, abs=1e-9)


def test_kepler_propagation_matches_numerical_integration() -> None:
    r0 = np.array([190929.0, 0.0])
    v0 = np.array([0.4, 6.6])
    dt = 9.3 * ub.DAY_S
    r1, v1 = ub.propagate(r0, v0, dt)
    sol = solve_ivp(_two_body, (0.0, dt), np.r_[r0, v0], method="DOP853", rtol=1e-13, atol=1e-6)
    assert np.linalg.norm(sol.y[:2, -1] - r1) < 1e-2
    assert np.linalg.norm(sol.y[2:, -1] - v1) < 1e-7


def test_landau_leg_hits_the_moon_under_independent_integration() -> None:
    ph = ub.MoonPhases(theta0={"Ariel": 0.0, "Umbriel": 1.0, "Titania": 2.0, "Oberon": 3.0})
    # depart Ariel on a tour-like orbit (pump about 80 deg, 4 km/s)
    _, vm = ph.state("Ariel", 0.0)
    that = vm / np.linalg.norm(vm)
    rhat = np.array([1.0, 0.0])
    vout = 4.0 * (math.cos(math.radians(80)) * that + math.sin(math.radians(80)) * rhat)
    legs = ub.landau_leg(
        dep_moon="Ariel",
        t_dep=0.0,
        vinf_out=vout,
        arr_moon="Oberon",
        t_arr=30 * ub.DAY_S,
        phases=ph,
    )
    assert legs, "expected at least one feasible leg"
    leg = legs[0]
    # 1) coast to the burn
    s1 = solve_ivp(
        _two_body,
        (leg.t_dep, leg.t_burn),
        np.r_[leg.r_dep, leg.v_dep],
        method="DOP853",
        rtol=1e-13,
        atol=1e-6,
    )
    assert np.linalg.norm(s1.y[:2, -1] - leg.r_burn) < 0.1
    assert np.linalg.norm(s1.y[2:, -1] - leg.v_burn_before) < 1e-6
    # 2) burn, coast to the arrival; the spacecraft must be at Oberon
    s2 = solve_ivp(
        _two_body,
        (leg.t_burn, leg.t_arr),
        np.r_[leg.r_burn, leg.v_burn_after],
        method="DOP853",
        rtol=1e-13,
        atol=1e-6,
    )
    r_moon, _ = ph.state("Oberon", leg.t_arr)
    assert np.linalg.norm(s2.y[:2, -1] - r_moon) < 1.0
    assert leg.min_radius_km >= ub.RP_FLOOR_KM


def test_mirror_closure_turn_identity() -> None:
    """Physics invariant: a (1,1) symmetric closure arrives at the middle moon on
    one crossing and leaves on the mirror crossing, so the demanded turn is
    2 * min(pump, 180 - pump) for the conic's pump angle there."""
    ph = ub.MoonPhases(theta0={"Ariel": 0.0, "Umbriel": 0.0, "Titania": 0.0, "Oberon": 0.0})
    cg = ub.CyclerGeometry("Ariel", "Oberon", 7.751820498940574 * ub.DAY_S, 0.0, 0)
    out0, in1 = cg.leg_vinf(ph, 0.0, "Ariel", "Oberon")
    out1, _ = cg.leg_vinf(ph, cg.leg_s, "Oberon", "Ariel")
    a, e = ub.orbit_from_vinf_pair(
        "Ariel", float(np.linalg.norm(out0)), "Oberon", float(np.linalg.norm(in1))
    )
    arr = ub.arrival_pump("Oberon", a, e)
    assert arr is not None
    pump = arr[1]
    expected = 2.0 * min(pump, math.pi - pump)
    assert abs(ub.signed_angle(in1, out1)) == pytest.approx(expected, abs=1e-6)


def test_same_conic_node_demands_no_turn() -> None:
    """Flyable control: Lambert legs on either side of a point of ONE conic recover
    that conic, so the demanded turn there is zero."""
    r0 = np.array([583511.0, 0.0])
    v0 = np.array([-0.6, 2.0])
    t1, t2 = 3.0 * ub.DAY_S, 5.5 * ub.DAY_S
    r1, _ = ub.propagate(r0, v0, t1)
    r2, _ = ub.propagate(r0, v0, t2)
    s_a = lambert(np.r_[r0, 0.0], np.r_[r1, 0.0], t1, mu=ub.MU_URANUS)[0]
    s_b = lambert(np.r_[r1, 0.0], np.r_[r2, 0.0], t2 - t1, mu=ub.MU_URANUS)[0]
    assert np.linalg.norm(np.asarray(s_a.v2[:2]) - np.asarray(s_b.v1[:2])) < 1e-8


def test_entry_correction_is_zero_inside_the_bend_and_the_chord_outside() -> None:
    vin = np.array([1.5, 0.0])
    inside = ub.rotate(vin, math.radians(3.0))
    assert ub.entry_correction_kms(vin, inside, math.radians(6.0)) == pytest.approx(0.0, abs=1e-12)
    outside = ub.rotate(vin, math.radians(40.0))
    chord = 2.0 * 1.5 * math.sin(math.radians(40.0 - 6.0) / 2.0)
    assert ub.entry_correction_kms(vin, outside, math.radians(6.0)) == pytest.approx(
        chord, rel=1e-12
    )


def test_ladder_finds_a_trivial_target_in_one_flyby() -> None:
    a, e = ub.orbit_from_vinf_pair("Oberon", 3.448, "Ariel", 4.325)
    st = ub.arrival_pump("Ariel", a, e)
    assert st is not None
    res = ub.tisserand_ladder(("Ariel", st[0], st[1]), {"Ariel": st}, max_depth=2)
    assert res.n_flybys == 1
