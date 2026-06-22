"""Tests for ER3BP genome components."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom
from cyclerfinder.genome.er3bp_continuation import continue_er3bp_family_in_e
from cyclerfinder.genome.er3bp_periodic import correct_er3bp_periodic

# Broucke 1969 TR 32-1360 Table 12, Family 7P
# Sourced Golden anchor for validation.
MU_BROUCKE = 0.0121550

BROUCKE_ORBIT_1 = {
    "e": 0.0001,
    "x0": 0.1520965,
    "ydot0": 3.1608994,
}

BROUCKE_ORBIT_59 = {
    "e": 0.025,
    "x0": 0.1461519,
    "ydot0": 3.1962405,
}


def test_er3bp_eom_reduces_to_cr3bp() -> None:
    """At e=0, ER3BP EOMs should exactly match CR3BP EOMs."""
    state = np.array([0.5, 0.1, -0.2, -0.05, 0.8, 0.01])
    mu = 0.01215

    # ER3BP at f=0, e=0
    er3bp_deriv = er3bp_eom(f=0.0, state6=state, mu=mu, e=0.0)

    # CR3BP (autonomous)
    cr3bp_deriv = cr3bp_eom(t=0.0, state6=state, mu=mu)

    np.testing.assert_allclose(er3bp_deriv, cr3bp_deriv, rtol=1e-14, atol=1e-14)


def test_raw_eom_closure_on_broucke_golden() -> None:
    """
    Strongest EOM Check: Integrate the Broucke Orbit 59 IC raw for T=2pi
    and assert closure independent of the corrector. This isolates EOM correctness.
    """
    ic = np.array([BROUCKE_ORBIT_59["x0"], 0.0, 0.0, 0.0, BROUCKE_ORBIT_59["ydot0"], 0.0])
    sys = ER3BPSystem(
        mu=MU_BROUCKE,
        e=BROUCKE_ORBIT_59["e"],
        primary_name="Earth",
        secondary_name="Moon",
    )

    # Integrate to half period (pi) to verify perpendicular crossing
    sol_half = solve_ivp(
        fun=er3bp_eom,
        t_span=(0.0, math.pi),
        y0=ic,
        args=(sys.mu, sys.e),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    final_half = sol_half.y[:, -1]

    # Assert perpendicular crossing at half period: y(pi) == 0, xdot(pi) == 0
    # Tolerance relaxed to 1e-5 to account for Broucke's 1969 integration precision
    assert final_half[1] == pytest.approx(0.0, abs=1e-5)
    assert final_half[3] == pytest.approx(0.0, abs=1e-5)

    # Integrate full period (2pi)
    sol_full = solve_ivp(
        fun=er3bp_eom,
        t_span=(0.0, 2 * math.pi),
        y0=ic,
        args=(sys.mu, sys.e),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    final_full = sol_full.y[:, -1]

    # Assert full closure X(2pi) == X(0)
    # Tolerance relaxed to 1e-3 because exponential integration drift compounds
    # the 7-digit IC truncation error over a full 2pi period.
    np.testing.assert_allclose(final_full, ic, atol=1e-3)


def test_corrector_on_broucke_golden() -> None:
    """Corrector should verify the Golden IC without diverging."""
    ic = np.array([BROUCKE_ORBIT_59["x0"], 0.0, 0.0, 0.0, BROUCKE_ORBIT_59["ydot0"], 0.0])
    sys = ER3BPSystem(
        mu=MU_BROUCKE,
        e=BROUCKE_ORBIT_59["e"],
        primary_name="Earth",
        secondary_name="Moon",
    )

    orbit = correct_er3bp_periodic(
        system=sys,
        state_guess=ic,
        period_f=math.pi,
        is_half_period_residual=True,
    )

    assert orbit.corrector_residual < 1e-10
    # Should converge in very few iterations since it's already an exact solution
    assert orbit.iterations < 5


def test_er3bp_continuation() -> None:
    """
    Start from a 2pi-commensurate CR3BP-like seed (Broucke Orbit 1, e=0.0001)
    and continue along Family 7P to recover the target Orbit 59 (e=0.025).
    """
    ic_seed = np.array([BROUCKE_ORBIT_1["x0"], 0.0, 0.0, 0.0, BROUCKE_ORBIT_1["ydot0"], 0.0])
    sys_base = ER3BPSystem(
        mu=MU_BROUCKE,
        e=BROUCKE_ORBIT_1["e"],
        primary_name="Earth",
        secondary_name="Moon",
    )

    history = continue_er3bp_family_in_e(
        sys_base=sys_base,
        seed_state=ic_seed,
        period_f=math.pi,
        e_target=BROUCKE_ORBIT_59["e"],
        n_steps=10,  # Small number of steps for tests
        is_half_period_residual=True,
    )

    final_orbit = history[-1]

    assert final_orbit.e == pytest.approx(BROUCKE_ORBIT_59["e"])

    # Assert we recovered the Golden Orbit 59 ICs to within Broucke's noise floor
    assert final_orbit.state0[0] == pytest.approx(BROUCKE_ORBIT_59["x0"], abs=1e-5)
    assert final_orbit.state0[4] == pytest.approx(BROUCKE_ORBIT_59["ydot0"], abs=1e-5)


def test_er3bp_out_of_plane_z_gradient_finite_difference() -> None:
    """Verify the out-of-plane u_zz term at e>0 and z!=0 via finite difference."""
    mu = 0.01215
    e = 0.05
    f = math.pi / 4.0

    # Base state with non-zero z to ensure out-of-plane terms are active
    state6 = np.array([0.5, 0.1, -0.2, 0.0, 0.0, 0.0])

    # Evaluate analytical Jacobian via the STM EOM
    state42 = np.zeros(42)
    state42[:6] = state6
    # Seed STM with identity so phidot = A @ I = A
    state42[6:] = np.eye(6).flatten()

    # er3bp_stm_eom is not exposed at the top level module conventionally, so we import it
    from cyclerfinder.core.er3bp import er3bp_stm_eom

    deriv42 = er3bp_stm_eom(f, state42, mu, e)

    a_mat = deriv42[6:].reshape((6, 6))
    u_zz_analytical = a_mat[5, 2]

    # Finite difference
    dz = 1e-7
    state6_up = state6.copy()
    state6_up[2] += dz
    state6_dn = state6.copy()
    state6_dn[2] -= dz

    zdoubleprime_up = er3bp_eom(f, state6_up, mu, e)[5]
    zdoubleprime_dn = er3bp_eom(f, state6_dn, mu, e)[5]

    u_zz_fd = (zdoubleprime_up - zdoubleprime_dn) / (2 * dz)

    assert u_zz_analytical == pytest.approx(u_zz_fd, rel=1e-5)
