"""Asymmetric (general) planar CR3BP periodic-orbit corrector at fixed Jacobi (#249).

Motivation (Braik-Ross C21)
---------------------------
The existing symmetric corrector
(:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`)
assumes a perpendicular x-axis crossing IC ``(x0, 0, 0, 0, ydot0, 0)`` -- i.e.
both ``y0 = 0`` AND ``xdot0 = 0`` -- and corrects ``x0`` so the orbit re-crosses
the x-axis perpendicularly (``y = 0`` AND ``xdot = 0``). An exhaustive symmetric
all-roots search caps the Earth-Moon (2,1) family at T ~= 70.6 d, so it cannot
represent Braik-Ross's C21 (2,1)-cycler (C_J = 3.1294, T = 84.533 d). The
leading hypothesis is that C21 is ASYMMETRIC: ``xdot0 != 0`` at its y=0
crossings, which the symmetric machinery structurally cannot capture.

This module DROPS the ``xdot0 = 0`` assumption. The IC still starts on the
``y = 0`` axis (``y0 = 0``) but allows ``xdot0 != 0``; ``ydot0`` is fixed
algebraically from the Jacobi constant. A genuine periodic orbit is found via
the ``y = 0`` return map: the full state must return to the seed at the
``half_crossings``-th y=0 crossing.

The math (planar, fixed C)
--------------------------
IC on the y=0 axis: ``state0 = (x0, 0, xdot0, ydot0)`` (planar (x, y, xdot,
ydot)). The repo's Jacobi convention (see
:func:`cyclerfinder.core.cr3bp.jacobi_constant`) is

    C = (x^2 + y^2) + 2(1-mu)/r1 + 2 mu/r2 - (xdot^2 + ydot^2)
      = Ubar_axis(x) - (xdot^2 + ydot^2)   at y = 0,

where ``Ubar_axis(x) = x^2 + 2(1-mu)/r1 + 2 mu/r2`` is
:func:`cyclerfinder.search.cr3bp_periodic._ubar_x_at_axis`. Holding C fixed,

    ydot0 = sign * sqrt( Ubar_axis(x0) - C - xdot0^2 ),

requiring the radicand ``>= 0`` (else the requested Jacobi constant is
infeasible at ``(x0, xdot0)``). ``sign`` follows the symmetric convention
(``ydot0_sign``, default -1, to match the family seeds).

Free variables: ``(x0, xdot0)`` (2 unknowns). Integrate to the
``half_crossings``-th y=0 crossing; call the state there
``(x_c, 0, xdot_c, ydot_c)``. The y=0 return-map fixed point requires

    R(x0, xdot0) = [ x_c - x0,  xdot_c - xdot0 ]  =  0   (2 equations).

``ydot`` closes automatically because the Jacobi constant is conserved by the
flow and held fixed at the IC, so once ``(x, xdot)`` match on the y=0 section
and y=0, the remaining velocity component ``ydot`` is fixed up to sign (and the
sign is fixed by the crossing direction); we additionally report
``ydot_residual = |ydot_c - ydot0|`` as a diagnostic. Newton on the 2x2 system
``R`` with a finite-difference Jacobian.

Period definition
-----------------
``period`` is the RETURN TIME: the time of the ``half_crossings``-th y=0
crossing at which the residual ``R`` is driven to zero (the y=0 return-map
fixed point). For a genuinely periodic orbit this is the full period (the full
state has returned to the IC). For the symmetric special case ``xdot0 = 0`` the
``half_crossings``-th crossing is the perpendicular re-crossing and this return
time equals the symmetric corrector's ``t_half`` -- so callers seeding a
symmetric member with ``half_crossings`` set to the symmetric corrector's value
recover ``period == t_half`` (NOT ``2 * t_half``); see the reproduce test. The
returned ``period`` is independently verified by a full-period re-propagation
that ``||X(period) - X(0)||`` is small (exposed as ``closure_residual``).

Symmetric orbits fall out as the ``xdot0 ~= 0`` special case (reproduce-before-
trust gate). Model: pure planar CR3BP (z = 0). Pure (math/numpy/scipy + the
CR3BP core).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.outcome_log import log_outcome


@dataclass(frozen=True)
class GeneralPeriodicOrbit:
    """An asymmetric (general) planar CR3BP periodic orbit at fixed Jacobi.

    The IC is ``(x0, 0, xdot0, ydot0)`` on the y=0 axis with ``ydot0`` fixed
    from the Jacobi constant. ``period`` is the y=0 return time at which the
    full state returns to the IC (see the module docstring). ``asymmetry`` is
    ``abs(xdot0)`` -- zero for a symmetric (perpendicular-crossing) orbit.
    """

    x0: float
    y0: float  # always 0.0 (IC on the y=0 axis)
    xdot0: float
    ydot0: float
    jacobi: float
    period: float
    converged: bool
    residual: float  # max|R| of the 2x2 return-map fixed-point system
    closure_residual: float  # ||X(period) - X(0)|| from an independent re-propagation
    ydot_residual: float  # |ydot_c - ydot0| diagnostic at the closing crossing
    n_iter: int
    asymmetry: float  # abs(xdot0)


def _ydot0_general(x0: float, xdot0: float, c_target: float, mu: float, sign: float) -> float:
    """ydot0 from the Jacobi constant at y=0 with a nonzero xdot0 (module math).

    ydot0 = sign * sqrt( Ubar_axis(x0) - C - xdot0^2 ).

    Raises
    ------
    ValueError
        If the radicand is negative (the requested Jacobi constant is not
        attainable at ``(x0, xdot0)``).
    """
    rad = cp._ubar_x_at_axis(x0, mu) - c_target - xdot0 * xdot0
    if rad < 0.0:
        raise ValueError(
            f"_ydot0_general: negative radicand {rad:.3e} at x0={x0:.6f}, "
            f"xdot0={xdot0:.6f}, C={c_target:.6f}"
        )
    return float(sign) * math.sqrt(rad)


def _return_crossing(
    system: cr3bp.CR3BPSystem,
    x0: float,
    xdot0: float,
    c_target: float,
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, NDArray[np.float64]] | None:
    """Integrate the y=0 IC and return ``(t, state6)`` at the ``half_crossings``-th y=0 crossing.

    Returns ``None`` if ``ydot0`` is infeasible or the crossing does not exist
    within ``t_hi``.
    """
    try:
        ydot0 = _ydot0_general(x0, xdot0, c_target, mu, sign)
    except ValueError:
        return None
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.direction = 0.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_hi),
        state0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
    )
    pairs = [
        (t, y) for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True) if t > 1e-6 * t_hi
    ]
    if len(pairs) < half_crossings:
        return None
    t_c, yf = pairs[half_crossings - 1]
    return float(t_c), np.asarray(yf, dtype=np.float64)


def _residual(
    system: cr3bp.CR3BPSystem,
    x0: float,
    xdot0: float,
    c_target: float,
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], float, NDArray[np.float64]] | None:
    """Return ``(R, t_c, state_c)`` for the 2x2 return-map fixed point.

    ``R = [x_c - x0, xdot_c - xdot0]``.
    """
    rc = _return_crossing(
        system, x0, xdot0, c_target, mu, sign, half_crossings, t_hi, rtol=rtol, atol=atol
    )
    if rc is None:
        return None
    t_c, state_c = rc
    res = np.array([state_c[0] - x0, state_c[3] - xdot0], dtype=np.float64)
    return res, t_c, state_c


def _return_crossing_stm(
    x0: float,
    xdot0: float,
    c_target: float,
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, NDArray[np.float64], NDArray[np.float64]] | None:
    """As :func:`_return_crossing` but also integrate the STM.

    Returns ``(t_c, state6_c, stm6x6_c)`` at the ``half_crossings``-th y=0
    crossing, or ``None`` if infeasible / no such crossing.
    """
    try:
        ydot0 = _ydot0_general(x0, xdot0, c_target, mu, sign)
    except ValueError:
        return None
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])
    y0 = np.concatenate([state0, np.eye(6).reshape(36)])

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _y_event.direction = 0.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom,
        (0.0, t_hi),
        y0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_y_event,
    )
    pairs = [
        (t, y) for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True) if t > 1e-6 * t_hi
    ]
    if len(pairs) < half_crossings:
        return None
    t_c, yf = pairs[half_crossings - 1]
    yf = np.asarray(yf, dtype=np.float64)
    return float(t_c), yf[:6], yf[6:].reshape(6, 6)


def _residual_jac_stm(
    x0: float,
    xdot0: float,
    c_target: float,
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]] | None:
    """Return ``(R, t_c, state_c, J)`` with an STM-based analytic 2x2 Jacobian.

    ``J = dR/d(x0, xdot0)`` where ``R = [x_c - x0, xdot_c - xdot0]`` at the
    ``half_crossings``-th y=0 crossing. The Jacobian accounts for:

      * the parameterised IC ``X0(p) = (x0, 0, xdot0, ydot0(p))`` with
        ``ydot0 = sign*sqrt(Ubar(x0) - C - xdot0^2)``, whose partials are
        ``d(ydot0)/d(x0) = (dUbar/dx)(x0) / (2*ydot0)`` and
        ``d(ydot0)/d(xdot0) = -xdot0 / ydot0`` (planar (x,y,xdot,ydot) cols
        0, 3, with ydot0 entering planar row 3);
      * the variation of the crossing TIME, projected out so that the y=0 event
        condition is held: for any component ``i``,
        ``dX_c^i/dp = Phi[i,:] dX0/dp - (Xdot_c^i / Xdot_c^y) * (Phi[y,:] dX0/dp)``
        with ``Xdot_c = f(X_c)`` and ``y`` the y-row (planar index 1).

    Pure analytic (STM) Jacobian -- avoids the finite-difference noise floor that
    caps the achievable residual over a long multi-crossing return arc.
    """
    rc = _return_crossing_stm(
        x0, xdot0, c_target, mu, sign, half_crossings, t_hi, rtol=rtol, atol=atol
    )
    if rc is None:
        return None
    t_c, state_c, stm = rc
    res = np.array([state_c[0] - x0, state_c[3] - xdot0], dtype=np.float64)

    try:
        ydot0 = _ydot0_general(x0, xdot0, c_target, mu, sign)
    except ValueError:
        return None
    if abs(ydot0) < 1e-14:
        return None
    # Planar indices in the full 6-vector: x=0, y=1, xdot=3, ydot=4.
    # ydot0 = sign*sqrt(Ubar_axis(x0) - C - xdot0^2), so
    #   d(ydot0)/dx0    = d(Ubar_axis)/dx / (2*ydot0),
    #   d(ydot0)/dxdot0 = -xdot0 / ydot0.
    # cp._ubar_grad_x_at_axis returns dU_pseudo/dx for the pseudo-potential
    # U_pseudo = -1/2(x^2+y^2) - (1-mu)/r1 - mu/r2; our Ubar_axis = -2*U_pseudo,
    # hence d(Ubar_axis)/dx = -2 * _ubar_grad_x_at_axis.
    dubar_axis_dx = -2.0 * cp._ubar_grad_x_at_axis(x0, mu)
    dydot0_dx0 = dubar_axis_dx / (2.0 * ydot0)
    dydot0_dxdot0 = -xdot0 / ydot0
    dx0_col = np.zeros(6)  # dX0/d(x0)
    dx0_col[0] = 1.0
    dx0_col[4] = dydot0_dx0
    dxdot0_col = np.zeros(6)  # dX0/d(xdot0)
    dxdot0_col[3] = 1.0
    dxdot0_col[4] = dydot0_dxdot0
    # State derivative at the crossing (for the crossing-time projection).
    fdot = cr3bp.cr3bp_eom(t_c, state_c, mu)
    ydotc = float(fdot[1])  # d(y)/dt at crossing
    if abs(ydotc) < 1e-14:
        return None
    jac = np.zeros((2, 2))
    for col, dx0_full in enumerate((dx0_col, dxdot0_col)):
        dxc = stm @ dx0_full  # naive (fixed-time) variation of the crossing state
        # Project out the crossing-time variation to hold y_c = 0.
        proj = float(dxc[1]) / ydotc
        dxc_constr = dxc - proj * fdot
        # R = [x_c - x0, xdot_c - xdot0]; subtract the direct IC dependence.
        jac[0, col] = float(dxc_constr[0]) - (1.0 if col == 0 else 0.0)
        jac[1, col] = float(dxc_constr[3]) - (1.0 if col == 1 else 0.0)
    return res, float(t_c), state_c, jac


def correct_general_periodic(
    system: cr3bp.CR3BPSystem,
    x0_guess: float,
    xdot0_guess: float,
    c_target: float,
    period_guess: float,
    *,
    half_crossings: int,
    ydot0_sign: float = -1.0,
    tol: float = 1e-11,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    t_hi_frac: float = 1.5,
    dx_cap: float = 0.1,
    max_backtrack: int = 30,
) -> GeneralPeriodicOrbit:
    """Asymmetric (general) fixed-Jacobi single-shooting corrector via the y=0 return map.

    Free variables ``(x0, xdot0)`` (``y0 = 0`` fixed; ``ydot0`` fixed
    algebraically from the Jacobi constant each iteration). Damped Newton on the
    2x2 return-map fixed-point residual ``R = [x_c - x0, xdot_c - xdot0]`` at the
    ``half_crossings``-th y=0 crossing, using an analytic STM Jacobian
    (:func:`_residual_jac_stm`) with a backtracking line search on ``||R||``.
    The line search is essential: a multi-crossing return arc can be very
    ill-conditioned (condition number ~1e7 for the C11a member), so the raw
    Newton step overshoots from a rough seed; backtracking restores the
    quadratic basin. See the module docstring for the math and the period
    definition.

    Parameters
    ----------
    half_crossings:
        1-based index of the y=0 crossing on which the FULL state returns to the
        IC (the return-map fixed-point crossing). For a symmetric seed this is
        TWICE the symmetric corrector's ``half_crossings`` (the perpendicular
        half-period crossing is at index ``H``; the full return is at ``2H``):
        the corrector then recovers ``xdot0 ~= 0`` and ``period == T`` (see the
        reproduce test).
    ydot0_sign:
        Branch sign of ``ydot0`` (default -1, matching the family seeds).
    period_guess:
        Used only to set the integration horizon ``t_hi = t_hi_frac *
        period_guess`` (large enough to contain ``half_crossings`` crossings).

    Returns
    -------
    GeneralPeriodicOrbit
        ``converged`` is True iff ``max|R| < tol`` at exit. ``closure_residual``
        is an INDEPENDENT full-``period`` re-propagation check; ``ydot_residual``
        a diagnostic on the auto-closing ``ydot`` component.
    """
    mu = float(system.mu)
    x0 = float(x0_guess)
    xdot0 = float(xdot0_guess)
    t_hi = t_hi_frac * abs(float(period_guess))

    def _rj(
        a: float, b: float
    ) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]] | None:
        return _residual_jac_stm(
            a, b, c_target, mu, ydot0_sign, half_crossings, t_hi, rtol=rtol, atol=atol
        )

    res_norm = float("inf")
    t_c = abs(float(period_guess))
    state_c: NDArray[np.float64] = np.zeros(6)
    n_iter = 0
    cur = _rj(x0, xdot0)
    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as iteration count
        if cur is None:
            break
        res, t_c, state_c, jac = cur
        res_norm = float(np.max(np.abs(res)))
        if res_norm < tol:
            break
        try:
            dz = np.linalg.solve(jac, -res)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jac, -res, rcond=None)
        dz = np.asarray(dz, dtype=np.float64)
        # Cap the raw step (per component) to keep the iterate on a physical branch.
        dz[0] = float(np.clip(dz[0], -dx_cap, dx_cap))
        dz[1] = float(np.clip(dz[1], -dx_cap, dx_cap))
        if float(np.max(np.abs(dz))) < 1e-15:
            break
        # Backtracking line search on ||R||: ill-conditioned return arcs make the
        # full Newton step overshoot far from the root; halve until it decreases.
        scale = 1.0
        improved = None
        for _ in range(max_backtrack):
            trial = _rj(x0 + scale * dz[0], xdot0 + scale * dz[1])
            if trial is not None and float(np.max(np.abs(trial[0]))) < res_norm:
                improved = (x0 + scale * dz[0], xdot0 + scale * dz[1], trial)
                break
            scale *= 0.5
        if improved is None:
            break  # no decrease found -- at the achievable floor
        x0, xdot0, cur = improved

    converged = res_norm < tol
    # Final ydot0 / state on the corrected IC.
    try:
        ydot0 = _ydot0_general(x0, xdot0, c_target, mu, ydot0_sign)
        feasible = True
    except ValueError:
        ydot0 = float("nan")
        feasible = False
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])
    jacobi = cr3bp.jacobi_constant(state0, mu) if feasible else float("nan")
    period = float(t_c)
    ydot_residual = abs(float(state_c[4]) - ydot0) if feasible else float("nan")

    # Independent full-period re-propagation closure check (Radau, different
    # integrator than the DOP853 corrector loop). Computed whenever the IC is
    # feasible -- it is a diagnostic on the returned orbit, independent of the
    # corrector's own ``converged`` flag.
    closure_residual = float("nan")
    if feasible:
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period),
            state0,
            args=(mu,),
            method="Radau",
            rtol=max(rtol, 1e-11),
            atol=max(atol, 1e-11),
        )
        if sol.success:
            closure_residual = float(np.linalg.norm(sol.y[:, -1] - state0))

    # Passive training-data capture (#210): NO-OP unless CYCLERFINDER_OUTCOME_LOG
    # is set. Side effect only; never a validation input.
    log_outcome(
        solver="cr3bp.correct_general_periodic",
        inputs={
            "x0_guess": float(x0_guess),
            "xdot0_guess": float(xdot0_guess),
            "c_target": float(c_target),
            "period_guess": float(period_guess),
            "half_crossings": int(half_crossings),
            "ydot0_sign": float(ydot0_sign),
            "mu": mu,
        },
        outcome={
            "converged": bool(converged),
            "residual": float(res_norm),
            "x0": float(x0),
            "xdot0": float(xdot0),
            "ydot0": float(ydot0),
            "jacobi": float(jacobi),
            "period": float(period),
            "closure_residual": float(closure_residual),
            "n_iter": int(n_iter),
        },
        meta={"primary": system.primary, "secondary": system.secondary},
    )
    return GeneralPeriodicOrbit(
        x0=float(x0),
        y0=0.0,
        xdot0=float(xdot0),
        ydot0=float(ydot0),
        jacobi=float(jacobi),
        period=period,
        converged=converged,
        residual=float(res_norm),
        closure_residual=closure_residual,
        ydot_residual=ydot_residual,
        n_iter=n_iter,
        asymmetry=abs(float(xdot0)),
    )
