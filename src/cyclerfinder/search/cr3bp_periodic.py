"""CR3BP periodic-orbit differential corrector (spec 2026-06-10, Phase 2).

STM-based single-shooting: from an initial guess, drive the state's return-to-start
residual to zero with Newton steps using the monodromy/STM. Pure.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp


@dataclass(frozen=True)
class PeriodicOrbit:
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    converged: bool
    closure_residual: float


def correct_periodic(
    system: cr3bp.CR3BPSystem,
    state0_guess: NDArray[np.float64],
    period_guess: float,
    *,
    tol: float = 1e-10,
    max_iter: int = 30,
    max_period_factor: float = 100.0,
) -> PeriodicOrbit:
    """Single-shooting periodicity correction.

    Free variables: the 6 initial-state components + the period T. Constraint:
    X(T) - X(0) = 0 (full-state periodicity). Newton step uses the STM Phi(T) and
    the time derivative at T:  [Phi - I | f(X(T))] dz = -(X(T) - X(0)).
    A min-norm least-squares step (6 eqns, 7 unknowns) is taken each iteration.
    Converged iff |X(T) - X(0)| < tol.

    *max_period_factor*: abort if the period exceeds
    ``max_period_factor * |period_guess|`` (Newton divergence guard; prevents
    unbounded integration time on poorly conditioned seeds).
    """
    s = np.asarray(state0_guess, float).copy()
    period = float(period_guess)
    period_lo = 1e-6 * abs(float(period_guess))
    period_hi = max_period_factor * abs(float(period_guess))
    residual = float("inf")
    for _ in range(max_iter):
        arc = cr3bp.propagate(system, s, period, with_stm=True)
        diff = arc.state_f - s
        residual = float(np.linalg.norm(diff))
        if residual < tol:
            break
        assert arc.stm is not None
        f_end = cr3bp.cr3bp_eom(period, arc.state_f, system.mu)
        jac = np.column_stack([arc.stm - np.eye(6), f_end])  # 6 x 7
        dz, *_ = np.linalg.lstsq(jac, -diff, rcond=None)  # min-norm step
        s = s + dz[:6]
        period = period + float(dz[6])
        if period <= period_lo or period >= period_hi:
            break
    converged = residual < tol and period_lo < period < period_hi
    return PeriodicOrbit(
        state0=s,
        period=period,
        jacobi=cr3bp.jacobi_constant(s, system.mu),
        converged=converged,
        closure_residual=residual,
    )


def crosscheck_periodic(
    system: cr3bp.CR3BPSystem,
    orbit: PeriodicOrbit,
    *,
    method: str = "Radau",
    rtol: float = 1e-11,
    atol: float = 1e-11,
    closure_tol: float = 1e-8,
    jacobi_tol: float = 1e-8,
) -> tuple[bool, float]:
    """Independent integrator cross-check for a corrected periodic orbit.

    Re-propagates ``orbit.state0`` for ``orbit.period`` using a DIFFERENT
    ``solve_ivp`` method (default: "Radau", an implicit Runge-Kutta), which is
    independent from the DOP853 used by :func:`correct_periodic` / :func:`propagate`.

    Returns ``(ok, dj)`` where:
      - ``ok`` is ``True`` iff the orbit re-closes within *closure_tol* AND the
        Jacobi constant is preserved within *jacobi_tol*;
      - ``dj = |C(T) - C(0)|`` is the absolute Jacobi drift over one period.
    """
    s0 = np.asarray(orbit.state0, float)
    c0 = cr3bp.jacobi_constant(s0, system.mu)

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        s0,
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    sf = sol.y[:, -1]
    closure = float(np.linalg.norm(sf - s0))
    dj = abs(cr3bp.jacobi_constant(sf, system.mu) - c0)
    ok = closure < closure_tol and dj < jacobi_tol
    return ok, dj
