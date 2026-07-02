"""Tests for the generic deflated-Newton root-enumeration primitive (#524).

Two closed-form POSITIVE CONTROLS for the machinery (both are Farrell,
Birkisson & Funke 2015's own worked examples for this exact technique, so
the EXPECTED roots trace to the published method, not to this module):

1. A cubic polynomial with three known real roots {1, 2, 3} -- classic
   scalar (N=1) deflation demo: repeated Newton from the SAME seed recovers
   all three roots in turn as each is deflated away.
2. A 2D system (unit circle intersect the line y=x) with two known distinct
   roots -- exercises the N=2 Jacobian path (M*J + outer(F, grad_M)).

Then a REAL CR3BP application, tied to the sourced #504 Pluto-Charon (3,2)
positive-control anchor (Ross-RT 2026 Table I row 4, already used as
``pluto_charon_kk_sweep.sweep_32_positive_control``'s seed): deflated Newton
with zero known roots, applied to the SAME scalar fixed-Jacobi
perpendicular-crossing residual ``correct_symmetric_fixed_jacobi`` solves,
must independently reproduce that corrector's converged ``x0`` -- an
agreement check between two independent implementations of Newton's method
on the same physical residual, not a self-computed golden.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.deflated_newton import (
    NewtonStopReason,
    enumerate_roots,
    newton_deflated,
)
from cyclerfinder.search.pluto_charon_kk_sweep import (
    _PC_ANCHOR_C,
    _PC_ANCHOR_HC,
    _PC_ANCHOR_T,
    _PC_ANCHOR_X0,
    PC_MU,
    make_pluto_charon_system,
)


def _cubic_residual(z: np.ndarray) -> np.ndarray | None:
    x = z[0]
    return np.array([(x - 1.0) * (x - 2.0) * (x - 3.0)])


def test_enumerate_roots_recovers_all_three_cubic_roots_from_same_seed() -> None:
    """Farrell's own worked example: (x-1)(x-2)(x-3)=0 from seed x=0 three times."""
    roots = enumerate_roots(_cubic_residual, [np.array([0.0])] * 6, tol=1e-12, max_iter=100)
    values = sorted(float(r[0]) for r in roots)
    assert len(values) == 3
    for found, expected in zip(values, [1.0, 2.0, 3.0], strict=True):
        assert abs(found - expected) < 1e-8


def test_newton_deflated_repels_from_a_known_root() -> None:
    """Deflating root 1.0 and starting AT it must not silently reconverge to it."""
    result = newton_deflated(
        _cubic_residual, np.array([1.0 + 1e-3]), [np.array([1.0])], tol=1e-12, max_iter=100
    )
    if result.stop_reason == NewtonStopReason.CONVERGED:
        assert result.z is not None
        assert abs(float(result.z[0]) - 1.0) > 1e-6, "deflation failed to repel from the known root"


def _circle_line_residual(z: np.ndarray) -> np.ndarray | None:
    x, y = z
    return np.array([x * x + y * y - 1.0, y - x])


def test_enumerate_roots_finds_both_circle_line_intersections() -> None:
    """Unit circle x^2+y^2=1 meets y=x at (+-1/sqrt(2), +-1/sqrt(2)) -- closed form."""
    seeds = [np.array([0.5, 0.5])] * 4
    roots = enumerate_roots(_circle_line_residual, seeds, tol=1e-12, max_iter=100)
    assert len(roots) == 2
    expected = 1.0 / np.sqrt(2.0)
    xs = sorted(float(r[0]) for r in roots)
    assert abs(xs[0] - (-expected)) < 1e-8
    assert abs(xs[1] - expected) < 1e-8


def _pc_crossing_residual(x0: float) -> float | None:
    """The same scalar residual correct_symmetric_fixed_jacobi drives to zero:
    xdot at the half_crossings-th x-axis crossing, at fixed Jacobi C=_PC_ANCHOR_C.
    """
    system = make_pluto_charon_system()
    try:
        ydot0 = cp.ydot0_from_jacobi(x0, _PC_ANCHOR_C, PC_MU, sign=-1.0)
    except ValueError:
        return None
    state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    t_hi = 1.25 * _PC_ANCHOR_T
    times, states = cp._xaxis_crossings(
        system, state0, t_hi, with_stm=False, rtol=1e-12, atol=1e-12
    )
    if len(times) < _PC_ANCHOR_HC:
        return None
    yf = states[_PC_ANCHOR_HC - 1]
    return float(yf[3])


def _pc_residual_fn(z: np.ndarray) -> np.ndarray | None:
    r = _pc_crossing_residual(float(z[0]))
    if r is None:
        return None
    return np.array([r])


def test_deflated_newton_reproduces_pc32_anchor_independently() -> None:
    """Plain (undeflated) Newton on the PC fixed-Jacobi crossing residual, seeded
    at the sourced #504 Ross-RT anchor x0, must land on the SAME root the
    existing correct_symmetric_fixed_jacobi corrector finds from the same seed
    -- independent-implementation agreement on real CR3BP data, PC mu, sourced
    C and half_crossings, not a value this module invented.
    """
    reference = cp.correct_symmetric_fixed_jacobi(
        make_pluto_charon_system(),
        _PC_ANCHOR_X0,
        _PC_ANCHOR_C,
        _PC_ANCHOR_T,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    assert reference.converged

    result = newton_deflated(
        _pc_residual_fn,
        np.array([_PC_ANCHOR_X0]),
        (),
        tol=1e-9,
        max_iter=60,
        step_cap=0.2,
    )
    assert result.stop_reason == NewtonStopReason.CONVERGED
    assert result.z is not None
    assert abs(float(result.z[0]) - reference.x0) < 1e-6
