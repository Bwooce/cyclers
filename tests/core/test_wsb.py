"""Belbruno weak-stability-boundary (WSB) surface tests (#378 Phase 1).

The `core/wsb.py` module implements Belbruno (2004) Chapter 3:

  * Task 1.1 -- Kepler two-body energy E_2 (Def 3.10 / eq 3.6, P_2-centred
    inertial) and the periapsis predicate (sigma: r-dot_23 = 0). Golden: the
    SIGN of E_2 at the collinear libration points is negative (the L-neck is
    two-body bound to the Moon -- the basis of Lemma 3.30, transit-orbits-are-
    ballistic-captures).
  * Task 1.2 -- analytic W surface (Lemma 3.21 / eq 3.29) + the parabolic
    closed-form golden C = +/- sqrt(2) + O(mu) (Lemma 3.34 / eq 3.39) and the
    C_1 validity boundary (Def 3.22).
  * Task 1.3 -- numerical stability-class one-revolution labelling (§3.2.1).

Goldens are SOURCED (Belbruno 2004 equations / closed forms), never values
this code computed (feedback_golden_tests_sourced_only). Where Belbruno's
printed number is in Hill-rescaled units (the mu^(2/3) L-point expansion), the
test asserts the SIGN / categorical fact that is frame-independent, and the
unit difference is documented in the module.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.wsb as wsb
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x


def test_kepler_energy_moon_l2_negative() -> None:
    """E_2 at the L_2 (and L_1) state is < 0 (two-body bound to the Moon).

    Belbruno Lemma 3.30: E_2(L_2) < 0, so a transit orbit through the L-neck IS
    a ballistic-capture transfer. The printed value -1.20187 is in Hill-rescaled
    coordinates (the mu^(2/3) expansion); the frame-independent sourced fact is
    the negative sign, which we assert here in the raw synodic E_2.
    """
    system = bcr4bp.andreu_default()
    mu = system.mu
    for point in ("L1", "L2"):
        xl = lagrange_collinear_x(mu, point)
        state = np.array([xl, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
        e2 = wsb.kepler_energy_moon(state, system)
        assert e2 < 0.0, f"E_2({point}) should be < 0 (bound), got {e2}"


def test_kepler_energy_moon_unbound_positive() -> None:
    """A fast state far from the Moon has E_2 > 0 (two-body unbound)."""
    system = bcr4bp.andreu_default()
    # Near the Moon but moving fast in the rotating frame => large inertial speed.
    state = np.array([1.0 - system.mu + 0.05, 0.0, 0.0, 0.0, 3.0, 0.0], dtype=np.float64)
    assert wsb.kepler_energy_moon(state, system) > 0.0


def test_is_periapsis_detects_radial_velocity_zero() -> None:
    """sigma predicate: periapsis iff Moon-relative radial rate r-dot_23 ~ 0.

    Construct a Moon-relative state whose inertial velocity is purely tangential
    (perpendicular to the Moon-relative position) -- that is periapsis/apoapsis
    (r-dot_23 = 0). A purely radial inertial velocity is NOT periapsis.
    """
    system = bcr4bp.andreu_default()
    mu = system.mu
    moon_x = 1.0 - mu
    # Position offset +x from Moon by 0.05. Tangential inertial velocity is +y.
    # In the rotating frame, inertial Xdot = v_rot + omega x X. For Xdot to be
    # purely +y (tangential), set v_rot so that omega x X is absorbed: here
    # X = (0.05, 0, 0), omega x X = (0, 0.05, 0). Choose v_rot = (0, vy, 0).
    # Then Xdot = (0, vy + 0.05, 0): radial component (x) is zero => periapsis.
    state_peri = np.array([moon_x + 0.05, 0.0, 0.0, 0.0, 0.2, 0.0], dtype=np.float64)
    assert wsb.is_periapsis(state_peri, system, tol=1e-9)

    # Radial inertial velocity: Xdot has a nonzero x-component along X.
    state_radial = np.array([moon_x + 0.05, 0.0, 0.0, 0.3, 0.2, 0.0], dtype=np.float64)
    assert not wsb.is_periapsis(state_radial, system, tol=1e-9)


def test_wsb_analytic_parabolic_limit() -> None:
    """Parabolic golden: wsb_analytic_c at e_2 -> 1, mu -> 0 returns +/- sqrt(2).

    Belbruno Lemma 3.34 / eq 3.39: for the mu = 0 parabolic orbit with respect
    to P_1, the Jacobi value on W-tilde is C = +/- sqrt(2). The analytic-W
    formula (eq 3.29) must reproduce this closed form in the parabolic limit.
    """
    # Direct (+) and retrograde (-) branches; mu -> 0, e_2 -> 1, r_2 = 1
    # (parabolic reference radius), residual A captures the +/- sqrt(2) limit.
    c_plus = wsb.wsb_analytic_c(r2=1.0, theta2=math.pi, e2=1.0, mu=0.0, branch="direct")
    c_minus = wsb.wsb_analytic_c(r2=1.0, theta2=math.pi, e2=1.0, mu=0.0, branch="retrograde")
    assert c_plus == pytest.approx(math.sqrt(2.0), abs=1e-9)
    assert c_minus == pytest.approx(-math.sqrt(2.0), abs=1e-9)


def test_wsb_validity_c1_earth_moon() -> None:
    """Def 3.22 validity boundary: C_1(Earth-Moon) ~ 3.184 from the L_1 Jacobi.

    C_1 is the Jacobi constant at L_1 (sourced from the CR3BP L-point, NOT a
    hardcoded literal). wsb_validity_ok(C, C_1) is True for C < C_1, False for
    C >= C_1. Belbruno: W exists for C in [2.22, C_1=3.184] (Earth-Moon).
    """
    system = bcr4bp.andreu_default()
    c1 = wsb.earth_moon_c1(system.mu)
    # Sourced sanity: Belbruno prints C_1 ~ 3.184 for Earth-Moon.
    assert c1 == pytest.approx(3.184, abs=0.01)
    assert wsb.wsb_validity_ok(3.0, c1)
    assert wsb.wsb_validity_ok(2.5, c1)
    assert not wsb.wsb_validity_ok(c1 + 0.1, c1)


def test_stability_class_labels() -> None:
    """One-rev stability classification: a bound periapsis is stable; a clearly
    escaping state is escape (§3.2.1 categorical labels, not magic numbers)."""
    system = bcr4bp.andreu_default()
    mu = system.mu
    moon_x = 1.0 - mu
    # Tightly bound near-circular periapsis ~3000 km from the Moon: stays bound.
    r_close = 0.05  # nondim (~19,000 km) low lunar periapsis
    # circular speed about the Moon (two-body) at r_close, inertial:
    v_circ = math.sqrt(mu / r_close)
    # Build a periapsis state: position +x of Moon, tangential inertial velocity.
    # rotating-frame vy chosen so inertial tangential speed = v_circ:
    # Xdot_y = v_rot_y + omega*X_x = v_rot_y + r_close => v_rot_y = v_circ - r_close
    state_bound = np.array(
        [moon_x + r_close, 0.0, 0.0, 0.0, v_circ - r_close, 0.0], dtype=np.float64
    )
    label_bound = wsb.stability_class(state_bound, system)
    assert label_bound == "stable", f"expected stable, got {label_bound}"

    # Clearly escaping: high inertial speed away from the Moon => unbound (E_2>0).
    state_escape = np.array([moon_x + r_close, 0.0, 0.0, 0.0, 5.0, 0.0], dtype=np.float64)
    label_escape = wsb.stability_class(state_escape, system)
    assert label_escape in ("escape", "unstable"), f"expected escape, got {label_escape}"
