"""Unit tests for the QBCP invariant torus corrector.

Per the design guidelines, all exclamation marks are avoided in comments and
docstrings.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.qbcp_torus import (
    correct_qbcp_torus,
    evaluate_qbcp_torus,
    se_lyapunov_to_qbcp_torus_seed,
    se_to_em_transform,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi


def test_se_to_em_transform_maps_sun_and_secondary_exactly() -> None:
    """#891/#892: the Sun-Earth primary lands on the QBCP's own Sun, (alpha_7, alpha_8),
    with its velocity, and the secondary (the Earth+Moon mass point) on the origin at
    rest. The Sun starts at angle pi in this frame and regresses. The reference Sun
    velocity is differentiated from the published series here, independently of the
    module's central difference."""
    q = qbcp.qbcp_default()
    mu_se = 1.0 / (q.mu_sun + 1.0)
    om = q.omega_sun_nondim
    c7 = np.array(qbcp._COEFFS_ALPHA7)
    c8 = np.array(qbcp._COEFFS_ALPHA8)
    k7 = np.arange(len(c7))
    k8 = np.arange(len(c8))
    al0 = qbcp.evaluate_alphas(0.0, q)
    assert math.atan2(al0[8], al0[7]) == pytest.approx(math.pi, abs=1e-12)
    for t in (0.0, 0.9, 3.7):
        al = qbcp.evaluate_alphas(t, q)
        th = q.theta_sun0 + om * t
        dxs = -float(np.sum(c7 * k7 * om * np.sin(k7 * th)))
        dys = float(np.sum(c8 * k8 * om * np.cos(k8 * th)))
        sun = se_to_em_transform(np.array([-mu_se, 0.0, 0.0, 0.0, 0.0, 0.0]), t, q, mu_se)
        assert np.max(np.abs(sun[:2] - al[7:9])) < 1e-10
        assert abs(sun[2]) < 1e-12
        assert abs(sun[3] - dxs) < 1e-7
        assert abs(sun[4] - dys) < 1e-7
        sec = se_to_em_transform(np.array([1.0 - mu_se, 0.0, 0.0, 0.0, 0.0, 0.0]), t, q, mu_se)
        assert np.max(np.abs(sec)) < 1e-12
    # The Sun regresses: its angle decreases from pi.
    al1 = qbcp.evaluate_alphas(0.1, q)
    assert math.atan2(al1[8], al1[7]) < math.pi - 0.05


def test_se_to_em_transform_flow_bound_at_zero_moon_mass() -> None:
    """Measured bound, not an identity (#891/#892, 2026-10-04): with the Moon's mass set
    to zero, a Sun-Earth state transformed to the QBCP frame and propagated there agrees
    with the same state propagated in the Sun-Earth CR3BP and then transformed, to
    7.5e-7 in position and 9.6e-7 in velocity after 6 time units (the coherent Sun
    motion is not the uniform circle of the CR3BP, and the time scale is the mean rate
    1 - omega_S). Bound set at 5e-6. The old time scale (1 + omega_S) misses by more
    than 1."""
    q = qbcp.qbcp_default()
    qs = qbcp.QBCPSystem(
        mu=0.0,
        mu_sun=q.mu_sun,
        a_sun_nondim=q.a_sun_nondim,
        omega_sun_nondim=q.omega_sun_nondim,
        theta_sun0=q.theta_sun0,
    )
    mu_se = 1.0 / (q.mu_sun + 1.0)
    s0 = np.array([1.0 - mu_se + 0.008, 0.004, 0.002, 0.001, -0.006, 0.0015])
    t0, dt = 1.3, 6.0
    u0 = se_to_em_transform(s0, t0, qs, mu_se)
    _, states = qbcp.propagate_qbcp_pv(u0, (t0, t0 + dt), qs, rtol=1e-12, atol=1e-12)

    def _se_then_transform(time_factor: float) -> np.ndarray:
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, time_factor * dt),
            s0,
            args=(mu_se,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        return se_to_em_transform(sol.y[:, -1], t0 + dt, qs, mu_se)

    miss = np.abs(states[-1] - _se_then_transform(1.0 - q.omega_sun_nondim))
    assert np.max(miss[:3]) < 5e-6
    assert np.max(miss[3:]) < 5e-6
    miss_old = np.abs(states[-1] - _se_then_transform(1.0 + q.omega_sun_nondim))
    assert np.max(miss_old[:3]) > 1.0


def test_correct_qbcp_torus_convergence() -> None:
    """Test that a Sun-Earth L2 Lyapunov orbit can be corrected as a QBCP torus."""
    qbcp_sys = qbcp.qbcp_default()
    mu_se = 1.0 / (qbcp_sys.mu_sun + 1.0)

    sys_se = cr3bp.CR3BPSystem(
        mu=mu_se,
        primary="Sun",
        secondary="Earth",
        l_km=qbcp_sys.a_sun_nondim * 384400.0,
        t_s=1.0,
    )

    x_earth = 1.0 - mu_se
    c_target = 3.0008
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=x_earth + 0.010,
        jacobi=c_target,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )

    assert orbit_se.converged

    # Generate seed
    n_samples = 5
    n_modes = 2
    x0, phase_pin_idx, amplitude_pin = se_lyapunov_to_qbcp_torus_seed(
        orbit_se, qbcp_sys, mu_se, n_samples=n_samples
    )

    # Correct at mu = 0.0
    sys_mu0 = qbcp.QBCPSystem(
        mu=0.0,
        mu_sun=qbcp_sys.mu_sun,
        a_sun_nondim=qbcp_sys.a_sun_nondim,
        omega_sun_nondim=qbcp_sys.omega_sun_nondim,
        theta_sun0=qbcp_sys.theta_sun0,
    )

    torus_mu0 = correct_qbcp_torus(
        sys_mu0, x0, n_modes, n_samples, phase_pin_idx, amplitude_pin, tol=1e-3
    )

    assert torus_mu0.converged
    assert torus_mu0.invariance_residual < 1e-3

    # Evaluate at theta_long = 0, theta_trans = 0
    state_00 = evaluate_qbcp_torus(torus_mu0, 0.0, 0.0)
    assert state_00.shape == (6,)
