"""CR3BP periodic-orbit differential corrector (spec 2026-06-10, Phase 2).

STM-based single-shooting: from an initial guess, drive the state's return-to-start
residual to zero with Newton steps using the monodromy/STM. Pure.
"""

from __future__ import annotations

import math
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


@dataclass(frozen=True)
class SymmetricOrbit:
    """A perpendicular-x-axis-crossing symmetric periodic orbit (Ross Eq. 9-12).

    The IC has the planar form ``(x0, 0, 0, 0, ydot0, 0)`` with ``ydot0`` fixed
    algebraically from the Jacobi constant; the period is ``T = 2 * t_half``
    where ``t_half`` is the time of the perpendicular crossing
    (``y = 0, xdot = 0``).
    """

    x0: float
    ydot0: float
    jacobi: float
    t_half: float
    period: float
    converged: bool
    crossing_residual: float  # |xdot(t_half)| at the corrected half-period
    n_iter: int


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


def _ubar_x_at_axis(x0: float, mu: float) -> float:
    """-2*Ubar evaluated on the x-axis (y=z=0): x^2 + 2(1-mu)/r1 + 2mu/r2.

    With the IC velocity purely in y, the Jacobi constant reduces to
    C = (-2*Ubar) - ydot0^2, so ydot0 = sqrt((-2*Ubar) - C).
    """
    r1 = abs(x0 + mu)
    r2 = abs(x0 - 1.0 + mu)
    return x0 * x0 + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2


def ydot0_from_jacobi(x0: float, jacobi: float, mu: float, *, sign: float = 1.0) -> float:
    """Solve C(x0,0,0,0,ydot0,0) = ``jacobi`` for ydot0 (Ross Eq. 9).

    Raises
    ------
    ValueError
        If the radicand ``(-2*Ubar) - C`` is negative (the requested Jacobi
        constant is not attainable at ``x0`` for a perpendicular crossing).
    """
    rad = _ubar_x_at_axis(x0, mu) - jacobi
    if rad < 0.0:
        raise ValueError(
            f"ydot0_from_jacobi: negative radicand {rad:.3e} at x0={x0:.6f}, C={jacobi:.6f}"
        )
    return float(sign) * math.sqrt(rad)


def _ubar_grad_x_at_axis(x0: float, mu: float) -> float:
    """dUbar/dx at (x0, 0, 0) with Ubar = -1/2(x^2+y^2) - (1-mu)/r1 - mu/r2.

    dUbar/dx = -x + (1-mu)(x+mu)/r1^3 + mu(x-1+mu)/r2^3.
    """
    r1 = abs(x0 + mu)
    r2 = abs(x0 - 1.0 + mu)
    return -x0 + (1.0 - mu) * (x0 + mu) / r1**3 + mu * (x0 - 1.0 + mu) / r2**3


def _xaxis_crossings(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t_hi: float,
    *,
    with_stm: bool,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], list[NDArray[np.float64]]]:
    """Integrate ``state0`` to ``t_hi`` and return all y=0 crossings (t>0).

    Returns ``(times, states)`` where ``times`` is the array of crossing times
    (excluding the t=0 root) and ``states`` the corresponding full state (with
    STM appended when ``with_stm``).
    """

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    rhs = cr3bp.cr3bp_stm_eom if with_stm else cr3bp.cr3bp_eom
    if with_stm:
        y0 = np.concatenate([np.asarray(state0, float), np.eye(6).reshape(36)])
    else:
        y0 = np.asarray(state0, float)

    sol = solve_ivp(
        rhs,
        (0.0, t_hi),
        y0,
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_lo = 1e-6 * t_hi
    pairs = [(t, y) for t, y in zip(t_events, y_events, strict=True) if t > t_lo]
    if not pairs:
        return np.array([]), []
    times = np.array([t for t, _ in pairs])
    states = [y for _, y in pairs]
    return times, states


def correct_symmetric_fixed_jacobi(
    system: cr3bp.CR3BPSystem,
    x0_guess: float,
    jacobi: float,
    period_guess: float,
    *,
    ydot0_sign: float = 1.0,
    half_crossings: int | None = None,
    tol: float = 1e-8,
    max_iter: int = 30,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    x0_bounds: tuple[float, float] = (-2.0, 2.0),
) -> SymmetricOrbit:
    """Fixed-Jacobi symmetric single-shooting corrector (Ross & Roberts-Tsoukkas
    2025, AAS 25-621, Eqs. 9-12).

    Holds the Jacobi constant ``jacobi`` fixed and finds ``x0`` so that the
    perpendicular x-axis-crossing orbit ``(x0, 0, 0, 0, ydot0(x0,C), 0)``
    re-crosses the x-axis perpendicularly (``y = 0`` AND ``xdot = 0``) at the
    half period. The single unknown is ``x0``; ``ydot0`` is re-derived from the
    Jacobi constraint each iteration (Eq. 9). The Newton step uses the STM at the
    half-period crossing (Eq. 12, expressed via the standard symmetric-orbit
    differential-correction identity rather than the paper's compressed render):

        d(xdot1)/dx0 = [Phi[3,0] + Phi[3,4] dydot0/dx0]
                       - xddot1/ydot1 [Phi[1,0] + Phi[1,4] dydot0/dx0]

    where the second group projects out the variation in crossing time so that
    the y=0 condition is preserved, and the Jacobi constraint couples
    ``dydot0/dx0 = -(dUbar/dx)(x0) / ydot0`` into the initial-velocity column
    (column 4) of the STM.

    A (k1,k2) cycler crosses the x-axis many times per period; the half-period
    perpendicular crossing is identified by its *index* among the y=0 crossings.
    ``half_crossings`` fixes that index (1-based); if ``None`` it is determined
    once from ``period_guess`` (the crossing nearest ``T/2`` on the seed orbit)
    and then held fixed across iterations — this keeps the Newton target on a
    single continuous branch. Converged iff ``|xdot(t_half)| < tol``. Period
    ``T = 2 * t_half``.
    """
    x0 = float(x0_guess)
    t_half_guess = 0.5 * float(period_guess)
    t_hi = 1.25 * float(period_guess)
    n_target = half_crossings
    crossing_res = float("inf")
    n_iter = 0
    t_half = t_half_guess
    yf: NDArray[np.float64] = np.zeros(6)
    lo, hi = x0_bounds
    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as iteration count
        ydot0 = ydot0_from_jacobi(x0, jacobi, system.mu, sign=ydot0_sign)
        state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
        times, states = _xaxis_crossings(system, state0, t_hi, with_stm=True, rtol=rtol, atol=atol)
        if len(times) == 0:
            break
        if n_target is None:
            # Fix the crossing index on the first iteration from T/2 proximity.
            n_target = int(np.argmin(np.abs(times - t_half_guess))) + 1
        # If the target crossing disappeared (orbit changed topology) fall back to
        # the crossing nearest the running half-period estimate.
        idx = int(np.argmin(np.abs(times - t_half))) if n_target > len(times) else n_target - 1
        ystm = states[idx]
        yf = ystm[:6]
        stm = ystm[6:].reshape(6, 6)
        t_half = float(times[idx])
        xdot1 = float(yf[3])
        crossing_res = abs(xdot1)
        if crossing_res < tol:
            break
        ydot1 = float(yf[4])
        f1 = cr3bp.cr3bp_eom(t_half, yf, system.mu)
        xddot1 = float(f1[3])  # d(xdot)/dt at crossing
        dydot0_dx0 = -_ubar_grad_x_at_axis(x0, system.mu) / ydot0
        dy_dx0 = float(stm[1, 0]) + float(stm[1, 4]) * dydot0_dx0
        dxdot_dx0 = float(stm[3, 0]) + float(stm[3, 4]) * dydot0_dx0
        if abs(ydot1) < 1e-14:
            break
        dxdot_total = dxdot_dx0 - xddot1 * dy_dx0 / ydot1
        if abs(dxdot_total) < 1e-14:
            break
        dx0 = -xdot1 / dxdot_total
        # Damp very large steps to keep the iterate on a physical branch.
        max_step = 0.2
        if abs(dx0) > max_step:
            dx0 = math.copysign(max_step, dx0)
        x0 = x0 + dx0
        x0 = min(max(x0, lo), hi)
    ydot0_final = ydot0_from_jacobi(x0, jacobi, system.mu, sign=ydot0_sign)
    period = 2.0 * t_half
    converged = crossing_res < tol
    return SymmetricOrbit(
        x0=x0,
        ydot0=ydot0_final,
        jacobi=cr3bp.jacobi_constant(np.array([x0, 0.0, 0.0, 0.0, ydot0_final, 0.0]), system.mu),
        t_half=t_half,
        period=period,
        converged=converged,
        crossing_residual=crossing_res,
        n_iter=n_iter,
    )


def barden_stability(
    system: cr3bp.CR3BPSystem,
    orbit: SymmetricOrbit,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, complex]:
    """Barden half-period monodromy stability index (Ross Eqs. 13-15).

    For a perpendicular-crossing symmetric orbit the full-period monodromy
    factors as ``M = G @ inv(Phi(T/2)) @ G @ Phi(T/2)`` with
    ``G = diag(1, -1, -1, 1, ... )`` acting on the planar (x, y, xdot, ydot)
    block (Barden 1994). Only the half-period STM is integrated. The stability
    parameter is ``nu = 1/2 (lambda + 1/lambda)`` for the nontrivial eigenpair
    ``(lambda, 1/lambda)``; the orbit is linearly stable iff ``|nu| < 1``.

    Returns ``(nu, lambda)`` with ``nu`` real-cast (its imaginary part is
    numerically ~0 for a genuine eigenpair) and ``lambda`` the selected
    nontrivial eigenvalue.
    """
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    arc = cr3bp.propagate(system, state0, orbit.t_half, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    phi_half = arc.stm
    # Planar CR3BP: the (x,y,vx,vy) block decouples from z. Use the 4x4 planar
    # sub-STM (indices 0,1,3,4) and the planar G = diag(1,-1,-1,1).
    idx = [0, 1, 3, 4]
    phi4 = phi_half[np.ix_(idx, idx)]
    g4 = np.diag([1.0, -1.0, -1.0, 1.0])
    monodromy = g4 @ np.linalg.inv(phi4) @ g4 @ phi4
    eigs = np.linalg.eigvals(monodromy)
    # Discard the eigenvalues nearest 1 (the trivial pair); the remaining
    # reciprocal pair carries the stability information.
    order = np.argsort(np.abs(eigs - 1.0))
    nontrivial = eigs[order[2:]]
    lam = complex(nontrivial[np.argmax(np.abs(nontrivial))])
    nu = 0.5 * (lam + 1.0 / lam)
    return float(nu.real), lam


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
