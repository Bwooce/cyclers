"""CR3BP dynamics core (plan 2026-06-10, Phase 1)."""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp


def test_jacobi_constant_value() -> None:
    # At a sample state the Jacobi constant matches the closed-form convention
    # C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2.
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0])
    x, y, z, vx, vy, vz = s
    r1 = math.hypot(x + mu, y, z)
    r2 = math.hypot(x - 1 + mu, y, z)
    expect = (x * x + y * y) + 2 * (1 - mu) / r1 + 2 * mu / r2 - (vx * vx + vy * vy + vz * vz)
    assert np.isclose(cr3bp.jacobi_constant(s, mu), expect)


def test_eom_shape_and_coriolis_sign() -> None:
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0])
    d = cr3bp.cr3bp_eom(0.0, s, mu)
    assert d.shape == (6,)
    # d[0:3] == velocity; ax includes +2*vy Coriolis term.
    assert np.allclose(d[0:3], s[3:6])


def test_jacobi_conserved_over_propagation() -> None:
    # Jacobi is conserved to ~1e-10 over a propagation (the integrator self-check).
    mu = 0.012277471
    s0 = np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0])
    c0 = cr3bp.jacobi_constant(s0, mu)
    arc = cr3bp.propagate(
        cr3bp.CR3BPSystem(mu=mu, primary="test", secondary="test", l_km=1.0, t_s=1.0), s0, 5.0
    )
    c1 = cr3bp.jacobi_constant(arc.state_f, mu)
    assert abs(c1 - c0) < 1e-9


def test_stm_matches_finite_difference() -> None:
    # The propagated 6x6 STM matches a finite-difference of the flow.
    mu = 0.012277471
    s0 = np.array([0.9, 0.0, 0.0, 0.0, -0.5, 0.0])
    sysm = cr3bp.CR3BPSystem(mu=mu, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    arc = cr3bp.propagate(sysm, s0, 1.0, with_stm=True)
    eps = 1e-6
    col0 = (
        cr3bp.propagate(sysm, s0 + np.array([eps, 0, 0, 0, 0, 0]), 1.0).state_f
        - cr3bp.propagate(sysm, s0 - np.array([eps, 0, 0, 0, 0, 0]), 1.0).state_f
    ) / (2 * eps)
    assert arc.stm is not None
    assert np.allclose(arc.stm[:, 0], col0, atol=1e-4)
