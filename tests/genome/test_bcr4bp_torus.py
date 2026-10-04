"""Unit tests for the BCR4BP invariant torus corrector.

Per the design guidelines, all exclamation marks are avoided in comments and
docstrings.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_torus import (
    correct_bcr4bp_torus,
    evaluate_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
    se_to_em_transform,
)
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi


def _identity_miss(mu_sun: float, se_time_factor: float, dt: float) -> tuple[float, float]:
    """Max position and velocity miss of the mu = 0 identity (#891).

    With the Moon's mass zero, the BCR4BP is exactly the Sun-Earth CR3BP seen from a
    frame turning at rate 1. Transform a Sun-Earth state at t0, propagate it in the
    BCR4BP for dt, and compare with the same state propagated in the Sun-Earth CR3BP
    for ``se_time_factor * dt`` and transformed at t0 + dt. The state is off the x
    axis, out of plane, a few Earth-Moon distances from the Earth, with a nonzero Sun
    phase and start time.
    """
    base = bcr4bp.andreu_default()
    sys0 = bcr4bp.BCR4BPSystem(
        mu=0.0,
        mu_sun=mu_sun,
        a_sun_nondim=base.a_sun_nondim,
        omega_sun_nondim=base.omega_sun_nondim,
        theta_sun0=0.7,
    )
    mu_se = 1.0 / (mu_sun + 1.0)
    s0 = np.array([1.0 - mu_se + 0.008, 0.004, 0.002, 0.001, -0.006, 0.0015])
    t0 = 1.3
    u0 = se_to_em_transform(s0, t0, sys0, mu_se)
    arc = bcr4bp.propagate_bcr4bp(sys0, u0, dt, t0=t0, rtol=1e-13, atol=1e-13)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, se_time_factor * dt),
        s0,
        args=(mu_se,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    u1: NDArray[np.float64] = se_to_em_transform(sol.y[:, -1], t0 + dt, sys0, mu_se)
    miss = np.abs(arc.state_f - u1)
    return float(np.max(miss[:3])), float(np.max(miss[3:]))


def test_se_to_em_transform_mu0_identity() -> None:
    """#891: the Sun-Earth to Earth-Moon map commutes with the flows at mu = 0.

    Measured 2026-10-04: 3.0e-8 in position and 3.2e-8 in velocity after 6 time units
    with the module's Sun mass, which misses the Kepler relation
    ``n_S^2 a_S^3 = 1 + mu_sun`` (n_S = 1 - omega_S) by 2.3e-8 relative; with the Sun
    mass set to satisfy it, 3.7e-12 (the integrator floor). The shipped transform
    before #891 (prograde Sun, factor 1 + omega_S) missed by 1.9e2.
    """
    base = bcr4bp.andreu_default()
    n_s = 1.0 - base.omega_sun_nondim
    kepler_mu_sun = n_s**2 * base.a_sun_nondim**3 - 1.0
    assert abs(kepler_mu_sun / base.mu_sun - 1.0) < 3e-8

    dpos, dvel = _identity_miss(base.mu_sun, n_s, 6.0)
    assert dpos < 1e-7
    assert dvel < 1e-7
    dpos, dvel = _identity_miss(kepler_mu_sun, n_s, 6.0)
    assert dpos < 1e-10
    assert dvel < 1e-10
    # Control: the old time scale (Sun-Earth time advancing at 1 + omega_S) fails.
    dpos, _ = _identity_miss(kepler_mu_sun, 1.0 + base.omega_sun_nondim, 6.0)
    assert dpos > 1.0


def test_correct_bcr4bp_torus_convergence() -> None:
    """Test that a Sun-Earth L2 Lyapunov orbit can be corrected as a BCR4BP torus.

    At mu = 0 the orbit is exactly an invariant torus of the BCR4BP, with rotation
    number ``2 pi (1 - omega_S) T_s / P_SE`` (one Sun synodic period ``T_s`` advances
    the Sun-Earth phase by ``n_S T_s``), so the corrected rho is checked against that.
    #891 (2026-10-04): the test used 5 samples and 2 modes, which converged only in
    the model with the Sun moving the wrong way; in the corrected model the 2-mode
    truncation floor of this orbit is 7.2e-4 and 11 samples (5 modes) are needed to
    reach 1e-6 (measured: residual 7.8e-7, rho within 4e-7 of the exact value).
    """
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
    n_samples = 11
    n_modes = 5
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
    t_s = 2.0 * math.pi / bcr_sys.omega_sun_nondim
    rho_exact = 2.0 * math.pi * (1.0 - bcr_sys.omega_sun_nondim) * t_s / orbit_se.period
    assert torus_mu0.rho == pytest.approx(rho_exact, abs=1e-5)

    # Evaluate at theta_long = 0, theta_trans = 0
    state_00 = evaluate_bcr4bp_torus(torus_mu0, 0.0, 0.0)
    assert state_00.shape == (6,)
