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
