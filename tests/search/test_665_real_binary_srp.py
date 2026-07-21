"""Regression tests for task #665's cannonball-SRP-augmented real-binary
dynamics (`cyclerfinder.search.real_binary_srp`).

Locks in the claims the module's own docstring makes: exact beta=0 reduction
to the pre-existing gravity-only machinery, conservation of the augmented
Jacobi-like quantity, and the load-bearing phi0-mirror-symmetry finding (only
phi0 in {0, pi} gives a genuinely periodic, not just numerically-"converged",
orbit).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.real_binary_srp as srp
from cyclerfinder.search.pluto_charon_kk_sweep import (
    _PC_ANCHOR_C,
    _PC_ANCHOR_HC,
    _PC_ANCHOR_T,
    _PC_ANCHOR_X0,
    PC_MU,
    make_pluto_charon_system,
)


@pytest.fixture(scope="module")
def pc_anchor_orbit() -> tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit]:
    """A known-good, already-converged gravity-only periodic orbit (the
    PC(3,2) anchor row also used by #504's/#665's own positive controls)."""
    system = make_pluto_charon_system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        _PC_ANCHOR_X0,
        _PC_ANCHOR_C,
        _PC_ANCHOR_T,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    assert orbit.converged
    return system, orbit


def test_srp_eom_reduces_to_gravity_only_at_beta_zero() -> None:
    state = np.array([0.7, 0.1, 0.0, 0.05, -0.3, 0.0])
    base = cr3bp.cr3bp_eom(0.0, state, PC_MU)
    aug = srp.cr3bp_srp_eom(0.0, state, PC_MU, beta_nd=0.0, phi0=0.0)
    np.testing.assert_array_equal(base, aug)
    aug_pi = srp.cr3bp_srp_eom(0.0, state, PC_MU, beta_nd=0.0, phi0=math.pi)
    np.testing.assert_array_equal(base, aug_pi)


def test_jacobi_srp_reduces_to_plain_jacobi_at_beta_zero() -> None:
    state = np.array([0.7, 0.1, 0.0, 0.05, -0.3, 0.0])
    c0 = cr3bp.jacobi_constant(state, PC_MU)
    c_srp = srp.jacobi_constant_srp(state, PC_MU, beta_nd=0.0, phi0=0.7)
    assert c_srp == pytest.approx(c0, abs=1e-14)


def test_srp_stm_variational_block_matches_gravity_only() -> None:
    """The SRP force is state-independent, so the A-matrix (columns 6:) of
    the augmented variational EOM must be bit-identical to the gravity-only
    one -- only the raw state-derivative (columns 0:6) differs."""
    y42 = np.concatenate([np.array([0.7, 0.1, 0.0, 0.05, -0.3, 0.0]), np.eye(6).reshape(36)])
    base_full = cr3bp.cr3bp_stm_eom(0.0, y42, PC_MU)
    aug_full = srp.cr3bp_srp_stm_eom(0.0, y42, PC_MU, beta_nd=1e-4, phi0=0.0)
    np.testing.assert_array_equal(base_full[6:], aug_full[6:])


def test_jacobi_srp_conserved_along_augmented_trajectory() -> None:
    """The augmented Jacobi-like quantity must be conserved along a
    trajectory of the augmented EOM (the whole point of deriving it) at a
    genuinely nonzero beta_nd."""
    from scipy.integrate import solve_ivp

    beta_nd, phi0 = 5e-4, 0.0
    state0 = np.array([0.7, 0.0, 0.0, 0.0, -0.6, 0.0])
    c0 = srp.jacobi_constant_srp(state0, PC_MU, beta_nd, phi0)
    sol = solve_ivp(
        srp.cr3bp_srp_eom,
        (0.0, 3.0),
        state0,
        args=(PC_MU, beta_nd, phi0),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.success
    assert sol.sol is not None
    for t in np.linspace(0.0, 3.0, 11):
        c_t = srp.jacobi_constant_srp(sol.sol(t), PC_MU, beta_nd, phi0)
        assert c_t == pytest.approx(c0, abs=1e-9)


def test_corrector_rejects_off_axis_phi0(
    pc_anchor_orbit: tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit],
) -> None:
    system, orbit = pc_anchor_orbit
    with pytest.raises(ValueError, match="phi0"):
        srp.correct_symmetric_fixed_jacobi_srp(
            system, orbit.x0, orbit.jacobi, orbit.period, beta_nd=1e-4, phi0=0.3, ydot0_sign=-1.0
        )


def test_srp_corrector_matches_gravity_only_at_beta_zero(
    pc_anchor_orbit: tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit],
) -> None:
    system, orbit = pc_anchor_orbit
    orbit_srp = srp.correct_symmetric_fixed_jacobi_srp(
        system,
        orbit.x0,
        orbit.jacobi,
        orbit.period,
        beta_nd=0.0,
        phi0=0.0,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    assert orbit_srp.converged
    assert orbit_srp.x0 == pytest.approx(orbit.x0, abs=1e-8)
    assert orbit_srp.ydot0 == pytest.approx(orbit.ydot0, abs=1e-8)
    assert orbit_srp.period == pytest.approx(orbit.period, abs=1e-6)


@pytest.mark.parametrize("phi0", [0.0, math.pi])
def test_onaxis_phi0_gives_genuinely_periodic_orbit(
    pc_anchor_orbit: tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit],
    phi0: float,
) -> None:
    """The module's own load-bearing claim: at phi0 in {0, pi} (SRP along
    the primary-secondary axis) the half-period corrector's result is a
    GENUINELY periodic orbit, confirmed by an independent full-period
    Radau crosscheck closing tightly."""
    system, orbit = pc_anchor_orbit
    beta_nd = 2e-4
    orbit_srp = srp.correct_symmetric_fixed_jacobi_srp(
        system,
        orbit.x0,
        orbit.jacobi,
        orbit.period,
        beta_nd=beta_nd,
        phi0=phi0,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    assert orbit_srp.converged
    state0 = np.array([orbit_srp.x0, 0.0, 0.0, 0.0, orbit_srp.ydot0, 0.0])
    ok, dj = srp.crosscheck_periodic_srp(
        system, state0, orbit_srp.period, beta_nd, phi0, closure_tol=1e-6, jacobi_tol=1e-8
    )
    assert ok, f"on-axis phi0={phi0} should give a genuinely closing periodic orbit (dj={dj:.3e})"
    assert dj < 1e-9


def _full_period_closure_srp(
    system: cr3bp.CR3BPSystem,
    state0: np.ndarray,
    period: float,
    beta_nd: float,
    phi0: float,
) -> float:
    """Full-period state closure residual (position/velocity mismatch, NOT
    the Jacobi delta) under the SRP-augmented flow -- the quantity the
    module docstring's phi0-mirror-symmetry claim is actually about."""
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        srp.cr3bp_srp_eom,
        (0.0, period),
        state0,
        args=(system.mu, beta_nd, phi0),
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    assert sol.success
    return float(np.linalg.norm(sol.y[:, -1] - state0))


def test_offaxis_phi0_breaks_periodicity_by_orders_of_magnitude(
    pc_anchor_orbit: tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit],
) -> None:
    """Off-axis SRP (sin(phi0) != 0) breaks the half-period mirror-symmetry
    the corrector relies on: the half-period crossing residual can still go
    to zero, but the resulting orbit is NOT genuinely periodic -- the
    FULL-period state closure (not the Jacobi delta, which is conserved
    along any trajectory of the autonomous augmented flow regardless of
    periodicity) is orders of magnitude worse than the on-axis case. This
    test bypasses the corrector's own guard (which correctly refuses
    off-axis phi0) to directly demonstrate WHY that guard exists: build a
    "solution" the same way the on-axis case is, but with the half-period
    ydot0 solved from an off-axis Jacobi relation."""
    system, orbit = pc_anchor_orbit
    beta_nd, phi0_offaxis = 2e-4, 0.3
    ydot0 = srp.ydot0_from_jacobi_srp(
        orbit.x0, orbit.jacobi, system.mu, beta_nd, phi0_offaxis, sign=-1.0
    )
    state0_offaxis = np.array([orbit.x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    closure_offaxis = _full_period_closure_srp(
        system, state0_offaxis, orbit.period, beta_nd, phi0_offaxis
    )

    orbit_onaxis = srp.correct_symmetric_fixed_jacobi_srp(
        system,
        orbit.x0,
        orbit.jacobi,
        orbit.period,
        beta_nd=beta_nd,
        phi0=0.0,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    state0_onaxis = np.array([orbit_onaxis.x0, 0.0, 0.0, 0.0, orbit_onaxis.ydot0, 0.0])
    closure_onaxis = _full_period_closure_srp(
        system, state0_onaxis, orbit_onaxis.period, beta_nd, 0.0
    )
    assert closure_onaxis < 1e-6
    assert closure_offaxis > 1e3 * max(closure_onaxis, 1e-9), (
        f"off-axis phi0 should degrade full-period closure by orders of magnitude vs on-axis "
        f"(off={closure_offaxis:.3e}, on={closure_onaxis:.3e})"
    )
