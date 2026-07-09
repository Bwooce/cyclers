"""Unit tests for the BCR4BP invariant torus corrector.

Per the design guidelines, all exclamation marks are avoided in comments and
docstrings.
"""

from __future__ import annotations

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_torus import (
    correct_bcr4bp_torus,
    evaluate_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi


def test_correct_bcr4bp_torus_convergence() -> None:
    """Test that a Sun-Earth L2 Lyapunov orbit can be corrected as a BCR4BP torus."""
    bcr_sys = bcr4bp.andreu_default()
    mu_SE = 1.0 / (bcr_sys.mu_sun + 1.0)

    sys_se = cr3bp.CR3BPSystem(
        mu=mu_SE, primary="Sun", secondary="Earth", l_km=bcr_sys.a_sun_nondim * 384400.0, t_s=1.0
    )

    x_earth = 1.0 - mu_SE
    C_target = 3.0008
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=x_earth + 0.010,
        jacobi=C_target,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )

    assert orbit_se.converged

    # Generate seed
    n_samples = 5
    n_modes = 2
    x0, phase_pin_idx, amplitude_pin = se_lyapunov_to_bcr4bp_torus_seed(
        orbit_se, bcr_sys, mu_SE, n_samples=n_samples
    )

    # Correct at mu = 0.0
    sys_mu0 = bcr4bp.BCR4BPSystem(
        mu=0.0,
        mu_sun=bcr_sys.mu_sun,
        a_sun_nondim=bcr_sys.a_sun_nondim,
        omega_sun_nondim=bcr_sys.omega_sun_nondim,
    )

    torus_mu0 = correct_bcr4bp_torus(
        sys_mu0, x0, n_modes, n_samples, phase_pin_idx, amplitude_pin, tol=1e-6
    )

    assert torus_mu0.converged
    assert torus_mu0.invariance_residual < 1e-6

    # Evaluate at theta_long = 0, theta_trans = 0
    state_00 = evaluate_bcr4bp_torus(torus_mu0, 0.0, 0.0)
    assert state_00.shape == (6,)
