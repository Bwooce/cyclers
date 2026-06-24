"""CR3BP analytic linear-seed generator (#435 Task 1).

Validates that the collinear-Lyapunov linear seed, after the fixed-Jacobi
symmetric corrector, is a genuine planar periodic orbit (it closes on
re-propagation over one full period). Golden = closure + the linear-regime
period band; the collinear linear-Lyapunov construction is textbook
(Szebehely / Koon-Lo-Marsden-Ross), so we pin behaviour, not an unsourced
exact initial condition.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import cr3bp_system, jacobi_constant, propagate
from cyclerfinder.search.cr3bp_seed_generator import dro_seed, lyapunov_seed, lyapunov_seed_3d


def test_lyapunov_seed_closes_on_earth_moon_l1() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-3)
    assert state0.shape == (6,)
    assert 2.0 < period < 3.5  # EM L1 linear-Lyapunov period band (TU)
    # Closure: re-propagate state0 one full period and compare to state0.
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-7)


def test_lyapunov_seed_jacobi_consistent() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-3)
    # Symmetric form: y=z=vx=vz=0 at the corrected start.
    assert abs(state0[1]) < 1e-12
    assert abs(state0[2]) < 1e-12
    assert abs(state0[3]) < 1e-12
    assert abs(state0[5]) < 1e-12
    # Jacobi conserved across the full period.
    c0 = jacobi_constant(state0, sysem.mu)
    arc = propagate(sysem, state0, period)
    assert abs(jacobi_constant(arc.state_f, sysem.mu) - c0) < 1e-9


def test_lyapunov_seed_arbitrary_mu_sun_earth_l1() -> None:
    # Arbitrary (small) mu path: Sun-Earth L1. At this very small mu (~3e-6) the
    # corrector converges on its |xdot(t_half)| residual (tol 1e-8) but the
    # full-period state-closure floors near ~4e-7 (the half-crossing residual
    # does not directly bound full-period closure); 1e-6 is still a tight,
    # meaningful closure and confirms the arbitrary-mu seed path.
    sysem = cr3bp_system("Sun", "Earth")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-4)
    assert state0.shape == (6,)
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-6)


def test_lyapunov_seed_bad_point_raises() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError):
        lyapunov_seed(sysem, point="L9", amplitude=1e-3)


def test_lyapunov_seed_3d_earth_moon_l1_closes_nonplanar() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed_3d(sysem, point="L1", amplitude_z=0.02)
    assert state0.shape == (6,)
    assert abs(state0[2]) > 0.01  # genuinely out-of-plane, not a planar collapse
    assert period > 0
    # Independent closure: propagate one full period, assert it returns near itself.
    arc = propagate(sysem, state0, period)
    assert np.linalg.norm(arc.state_f - state0) < 1e-6


def test_dro_seed_closes_on_earth_moon() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = dro_seed(sysem, amplitude=5e-2)
    assert state0.shape == (6,)
    assert period > 0.0
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-7)
