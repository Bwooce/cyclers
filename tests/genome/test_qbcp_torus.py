"""Unit tests for the QBCP invariant torus corrector.

Per the design guidelines, all exclamation marks are avoided in comments and
docstrings.
"""

from __future__ import annotations

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.qbcp_torus import (
    correct_qbcp_torus,
    evaluate_qbcp_torus,
    se_lyapunov_to_qbcp_torus_seed,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi


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
