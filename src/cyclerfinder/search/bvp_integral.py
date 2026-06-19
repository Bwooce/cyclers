"""BVP integral-constraint corrector core (#380 Steps 1-3).

Augmented-STM quadrature for periodic-orbit / two-point-BVP correction with
INTEGRAL constraints of the form ``q_i(T) = integral_0^T h_i(t, X(t)) dt``.

Approach (per ``docs/notes/2026-06-19-380-bvp-integral-corrector-blueprint.md``)
-------------------------------------------------------------------------------
Append ``n_q`` integral-accumulator rows and their ``6 * n_q`` variational rows
to the existing 42-dim ``[X(6), Phi(36)]`` ODE so that each integral value
``q_i(T)`` and its sensitivity ``dq_i/dX0`` come out of ONE ``solve_ivp`` call:

  * ``dq_i/dt   = h_i(X(t), t)``          with ``q_i(0) = 0``    (value row)
  * ``dpsi_i/dt = (dh_i/dX(X, t)) @ Phi``  with ``psi_i(0) = 0``  (1x6 sens. row)

so ``dq_i/dX0 = psi_i(T)``. Augmented dim = ``42 + 7 * n_q``. The Newton
residual stacks the existing POINT rows (state closure) with the integral rows
``(q_i(T) - target_i) * weight_i``; the Jacobian gains rows
``[dq_i/dX0 (free state cols) | h_i(X_f, T) * t_scale (T col)]``. Step via
``np.linalg.lstsq`` (min-norm / overdetermined), backtracking on the TOTAL
residual norm, per-component step caps mirroring
:mod:`cyclerfinder.search.cr3bp_general_periodic_3d`. An independent Radau
re-propagation cross-checks the POINT residual only.

Scope
-----
Steps 1-3 only: this is a self-contained, additive module. The wiring into the
existing correctors (Steps 4-6: #378 Phase 3.2, #388 multi-arc) is DEFERRED. No
existing file is modified by this module.

References
----------
  * Bond & Allman (2021), *Modern Astrodynamics*, Ch. 14 (augmented STM /
    sensitivity quadrature).
  * Shepperd (1985), quadrature appended to universal-variable propagation.
  * Pellegrini & Russell (2016), JGCD DOI 10.2514/1.G001920 (the variable-step
    STM bias caveat that also applies to the augmented rows here; the corrector
    is single-shooting-residual-dominated so the bias does not change the
    Newton-step direction, identical to the rationale in
    :func:`cyclerfinder.search.cr3bp_general_periodic_3d._propagate_with_stm`).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp

# State-vector index aliases (full 6D rotating-frame state); IDX_T marks the
# period as the 7th free variable. Mirrors
# :mod:`cyclerfinder.search.cr3bp_general_periodic_3d`.
IDX_X = 0
IDX_Y = 1
IDX_Z = 2
IDX_XDOT = 3
IDX_YDOT = 4
IDX_ZDOT = 5
IDX_T = 6


# ---------------------------------------------------------------------------
# Public dataclasses.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IntegralConstraint:
    """One scalar integral constraint ``q(T) = integral_0^T h(t, X) dt = target``.

    Attributes
    ----------
    label :
        Human-readable name (diagnostics / logging only).
    integrand :
        ``h(t, X, params) -> float``. ``X`` is the 6-vector state; ``params`` is
        the propagation system object (:class:`cr3bp.CR3BPSystem` or
        :class:`bcr4bp.BCR4BPSystem`). Must be finite along the trajectory.
    integrand_grad :
        ``dh/dX(t, X, params) -> NDArray`` of shape ``(6,)``, the gradient of
        ``h`` w.r.t. the instantaneous state. If ``None`` a central
        finite-difference fallback is used (slower, less accurate).
    target :
        The desired value of ``q(T)``.
    weight :
        Scalar applied to the integral residual row before stacking, to control
        ``lstsq`` conditioning relative to the point rows (default 1.0).
    """

    label: str
    integrand: Callable[[float, NDArray[np.float64], Any], float]
    integrand_grad: Callable[[float, NDArray[np.float64], Any], NDArray[np.float64]] | None
    target: float
    weight: float = 1.0


@dataclass(frozen=True)
class AugmentedArc:
    """Result of an augmented-STM propagation.

    Attributes
    ----------
    state_f :
        Final state ``X(T)`` (6,).
    stm :
        Monodromy block ``Phi(T, 0)`` (6, 6).
    q_values :
        Integral values ``q_i(T)`` (n_q,), in constraint order.
    dq_dx0 :
        Sensitivities ``dq_i/dX0`` (n_q, 6), in constraint order.
    """

    state_f: NDArray[np.float64]
    stm: NDArray[np.float64]
    q_values: NDArray[np.float64]
    dq_dx0: NDArray[np.float64]


@dataclass(frozen=True)
class IntegralCorrectorResult:
    """A converged (or attempted) integral-constrained correction.

    Attributes
    ----------
    state0 :
        Corrected initial state (6,).
    period :
        Corrected period (nondim TU).
    q_values :
        Integral values ``q_i(T)`` at the corrected orbit (n_q,).
    point_residual :
        L2 norm of the point (state-closure) residual rows.
    integral_residual :
        L2 norm of the weighted integral residual rows.
    total_residual :
        L2 norm of the full stacked residual (point + integral).
    converged :
        True iff ``total_residual < tol`` AND the independent closure check on
        the point residual is below ``independent_tol``.
    independent_closure_residual :
        L2 norm of the full-period state closure from an independent Radau
        re-propagation (point rows only).
    n_iter :
        Newton iterations consumed.
    """

    state0: NDArray[np.float64]
    period: float
    q_values: NDArray[np.float64]
    point_residual: float
    integral_residual: float
    total_residual: float
    converged: bool
    independent_closure_residual: float
    n_iter: int


# ---------------------------------------------------------------------------
# Augmented-STM propagators.
# ---------------------------------------------------------------------------


def _finite_diff_grad(
    integrand: Callable[[float, NDArray[np.float64], Any], float],
    t: float,
    state6: NDArray[np.float64],
    params: Any,
    *,
    eps: float = 1e-7,
) -> NDArray[np.float64]:
    """Central finite-difference of ``integrand`` w.r.t. the 6 state components."""
    grad = np.zeros(6, dtype=np.float64)
    for k in range(6):
        sp = state6.copy()
        sm = state6.copy()
        h = eps * max(1.0, abs(float(state6[k])))
        sp[k] += h
        sm[k] -= h
        grad[k] = (integrand(t, sp, params) - integrand(t, sm, params)) / (2.0 * h)
    return grad


def _augmented_eom(
    base_eom: Callable[[float, NDArray[np.float64], Any], NDArray[np.float64]],
    stm_eom: Callable[[float, NDArray[np.float64], Any], NDArray[np.float64]],
    constraints: Sequence[IntegralConstraint],
    params: Any,
) -> Callable[[float, NDArray[np.float64]], NDArray[np.float64]]:
    """Build the augmented RHS over ``[X(6), Phi(36), q(n_q), psi(6*n_q)]``.

    The first 42 rows are exactly the existing variational EOM ``stm_eom``; the
    next ``n_q`` rows accumulate ``q_i`` and the final ``6 * n_q`` rows
    accumulate ``psi_i = dq_i/dX0`` via ``(dh_i/dX) @ Phi``.
    """
    n_q = len(constraints)

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        d42 = stm_eom(t, y[:42], params)
        state6 = y[:6]
        phi = y[6:42].reshape(6, 6)
        out = np.empty(42 + 7 * n_q, dtype=np.float64)
        out[:42] = d42
        for i, c in enumerate(constraints):
            out[42 + i] = c.integrand(t, state6, params)
            if c.integrand_grad is not None:
                grad = np.asarray(c.integrand_grad(t, state6, params), dtype=np.float64)
            else:
                grad = _finite_diff_grad(c.integrand, t, state6, params)
            psi_row = grad @ phi  # (6,)
            base = 42 + n_q + 6 * i
            out[base : base + 6] = psi_row
        return out

    return rhs


def _propagate_augmented(
    base_eom: Callable[[float, NDArray[np.float64], Any], NDArray[np.float64]],
    stm_eom: Callable[[float, NDArray[np.float64], Any], NDArray[np.float64]],
    params: Any,
    state6: NDArray[np.float64],
    t: float,
    constraints: Sequence[IntegralConstraint],
    *,
    t0: float,
    rtol: float,
    atol: float,
) -> AugmentedArc:
    """Core augmented propagation shared by the CR3BP / BCR4BP entry points."""
    state_arr = np.asarray(state6, dtype=np.float64)
    if state_arr.shape != (6,):
        raise ValueError(f"propagate_augmented: state6 must have shape (6,); got {state_arr.shape}")
    n_q = len(constraints)

    y0 = np.zeros(42 + 7 * n_q, dtype=np.float64)
    y0[:6] = state_arr
    y0[6:42] = np.eye(6).reshape(36)
    # q(0) = 0 and psi(0) = 0 already (zeros init).

    rhs = _augmented_eom(base_eom, stm_eom, constraints, params)
    sol = solve_ivp(
        rhs,
        (t0, t0 + t),
        y0,
        rtol=rtol,
        atol=atol,
        method="DOP853",
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"augmented propagation failed at t={sol.t[-1]}: {sol.message}")
    yf = sol.y[:, -1]
    state_f = yf[:6]
    stm = yf[6:42].reshape(6, 6)
    q_values = yf[42 : 42 + n_q].copy()
    dq_dx0 = np.zeros((n_q, 6), dtype=np.float64)
    for i in range(n_q):
        base = 42 + n_q + 6 * i
        dq_dx0[i, :] = yf[base : base + 6]
    return AugmentedArc(state_f=state_f, stm=stm, q_values=q_values, dq_dx0=dq_dx0)


def propagate_augmented_cr3bp(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t: float,
    constraints: Sequence[IntegralConstraint],
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> AugmentedArc:
    """Augmented-STM CR3BP propagation with integral-constraint quadrature.

    When ``constraints`` is empty this falls back to the existing plain STM
    propagation (:func:`cr3bp.propagate` with ``with_stm=True``) — identical
    result and cost; the augmented ODE is not run.
    """
    if len(constraints) == 0:
        arc = cr3bp.propagate(system, state0, t, with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        return AugmentedArc(
            state_f=arc.state_f,
            stm=arc.stm,
            q_values=np.zeros(0, dtype=np.float64),
            dq_dx0=np.zeros((0, 6), dtype=np.float64),
        )
    return _propagate_augmented(
        cr3bp.cr3bp_eom,
        cr3bp.cr3bp_stm_eom,
        system.mu,
        state0,
        t,
        constraints,
        t0=0.0,
        rtol=rtol,
        atol=atol,
    )


def propagate_augmented_bcr4bp(
    system: bcr4bp.BCR4BPSystem,
    state0: NDArray[np.float64],
    t: float,
    constraints: Sequence[IntegralConstraint],
    *,
    t0: float = 0.0,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> AugmentedArc:
    """Augmented-STM BCR4BP propagation with integral-constraint quadrature.

    The BCR4BP EOM is time-dependent (the Sun phase advances), so ``t0`` (the
    absolute time origin) matters; supply it when stitching arcs. When
    ``constraints`` is empty this falls back to the existing plain STM
    propagation (:func:`bcr4bp.propagate_bcr4bp` with ``with_stm=True``).
    """
    if len(constraints) == 0:
        arc = bcr4bp.propagate_bcr4bp(system, state0, t, with_stm=True, t0=t0, rtol=rtol, atol=atol)
        assert arc.stm is not None
        return AugmentedArc(
            state_f=arc.state_f,
            stm=arc.stm,
            q_values=np.zeros(0, dtype=np.float64),
            dq_dx0=np.zeros((0, 6), dtype=np.float64),
        )
    return _propagate_augmented(
        bcr4bp.bcr4bp_eom,
        bcr4bp.bcr4bp_stm_eom,
        system,
        state0,
        t,
        constraints,
        t0=t0,
        rtol=rtol,
        atol=atol,
    )


# ---------------------------------------------------------------------------
# Newton corrector.
# ---------------------------------------------------------------------------


def correct_with_integral_constraints(
    propagate_fn: Callable[[NDArray[np.float64], float], AugmentedArc],
    state_guess: NDArray[np.float64] | Sequence[float],
    period_guess: float,
    *,
    free_vars: tuple[int, ...] = (IDX_X, IDX_Y, IDX_Z, IDX_XDOT, IDX_YDOT, IDX_ZDOT, IDX_T),
    point_residual_indices: tuple[int, ...] = (
        IDX_X,
        IDX_Y,
        IDX_Z,
        IDX_XDOT,
        IDX_YDOT,
        IDX_ZDOT,
    ),
    is_half_period_residual: bool = False,
    integral_constraints: Sequence[IntegralConstraint] = (),
    integrand_params: Any = None,
    independent_eom: Callable[[float, NDArray[np.float64]], NDArray[np.float64]] | None = None,
    tol: float = 1e-10,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    state_step_cap: float = 0.1,
    period_step_cap_frac: float = 0.2,
    independent_tol: float = 1e-6,
    max_backtrack: int = 20,
) -> IntegralCorrectorResult:
    """Newton corrector with point (state-closure) AND integral constraints.

    Parameters
    ----------
    propagate_fn :
        ``(state6, t) -> AugmentedArc``. A closure binding a system + the
        ``integral_constraints`` (use a ``lambda`` around
        :func:`propagate_augmented_cr3bp` / :func:`propagate_augmented_bcr4bp`).
        The arc's ``q_values`` / ``dq_dx0`` must be in the SAME order as
        ``integral_constraints``.
    state_guess, period_guess :
        Initial 6-vector IC and period (nondim TU).
    free_vars :
        Indices in ``{0..6}`` (6 = period) the corrector may update.
    point_residual_indices :
        Indices in ``{0..5}`` whose ``X(t_event) - X(0)`` must vanish.
    is_half_period_residual :
        If True the point residual is evaluated at ``T/2`` (the integral rows
        still use the FULL period — they are accumulated over the whole arc the
        propagator integrates).
    integral_constraints :
        The integral constraints. Their order MUST match ``propagate_fn``'s
        arc outputs.
    integrand_params :
        The system object (or bare mu) handed to ``IntegralConstraint.integrand``
        ONLY for the period-column Jacobian term ``h_i(X_f, T)``. Constant
        integrands (e.g. the time integral, h=1) ignore it, so ``None`` is fine;
        supply the propagation system when steering on a ``params``-dependent
        integrand with the period free.
    independent_eom :
        ``f(t, X) -> dX/dt`` for the independent Radau closure cross-check on
        the POINT residual. If ``None`` the cross-check is skipped (residual
        reported as NaN and ``converged`` falls back to the corrector residual).
    tol, max_iter, rtol, atol, state_step_cap, period_step_cap_frac,
    independent_tol, max_backtrack :
        Standard corrector knobs, mirroring
        :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`.

    Returns
    -------
    IntegralCorrectorResult
    """
    free_vars = tuple(sorted(set(int(v) for v in free_vars)))
    point_residual_indices = tuple(sorted(set(int(v) for v in point_residual_indices)))
    if not all(0 <= v <= IDX_T for v in free_vars):
        raise ValueError(f"free_vars must be in [0,6]; got {free_vars}")
    if not all(0 <= v <= IDX_ZDOT for v in point_residual_indices):
        raise ValueError(f"point_residual_indices must be in [0,5]; got {point_residual_indices}")
    if not free_vars:
        raise ValueError("free_vars must contain at least one index")

    state0 = np.asarray(state_guess, dtype=np.float64).copy()
    if state0.shape != (6,):
        raise ValueError(f"state_guess must have shape (6,); got {state0.shape}")
    period = float(period_guess)
    if period <= 0.0:
        raise ValueError(f"period_guess must be > 0; got {period}")

    n_point = len(point_residual_indices)
    n_q = len(integral_constraints)
    n_res = n_point + n_q
    n_free = len(free_vars)
    weights = np.array([c.weight for c in integral_constraints], dtype=np.float64)
    targets = np.array([c.target for c in integral_constraints], dtype=np.float64)

    def _event_time(period_in: float) -> float:
        return 0.5 * period_in if is_half_period_residual else period_in

    def _evaluate(
        state_in: NDArray[np.float64], period_in: float
    ) -> tuple[NDArray[np.float64], float, AugmentedArc, AugmentedArc | None] | None:
        """Propagate and assemble the stacked residual vector.

        Returns ``(residual_vec, residual_norm, arc_point, arc_full)`` where
        ``arc_point`` is the propagation to the point-residual event time and
        ``arc_full`` is the full-period propagation carrying the integral
        accumulators (None when the two coincide).
        """
        t_point = _event_time(period_in)
        try:
            arc_point = propagate_fn(state_in, t_point)
        except RuntimeError:
            return None
        # Integral rows are accumulated over the FULL period; if the point
        # event is the full period these coincide.
        if is_half_period_residual and n_q > 0:
            try:
                arc_full = propagate_fn(state_in, period_in)
            except RuntimeError:
                return None
        else:
            arc_full = None
        res = np.empty(n_res, dtype=np.float64)
        d = arc_point.state_f - state_in
        res[:n_point] = d[list(point_residual_indices)]
        if n_q > 0:
            q_arc = arc_full if arc_full is not None else arc_point
            res[n_point:] = (q_arc.q_values - targets) * weights
        if not np.all(np.isfinite(res)):
            return None
        return res, float(np.linalg.norm(res)), arc_point, arc_full

    def _residual_only(state_in: NDArray[np.float64], period_in: float) -> float:
        out = _evaluate(state_in, period_in)
        return float("inf") if out is None else out[1]

    def _build_trial(
        cur_state: NDArray[np.float64], cur_period: float, step: NDArray[np.float64], scale_: float
    ) -> tuple[NDArray[np.float64], float]:
        ts = cur_state.copy()
        tp = cur_period
        for col_, unknown_ in enumerate(free_vars):
            if unknown_ == IDX_T:
                tp = cur_period + scale_ * float(step[col_])
                if not math.isfinite(tp) or tp <= 0.0:
                    tp = max(cur_period * 0.5, 1e-6)
            else:
                ts[unknown_] = cur_state[unknown_] + scale_ * float(step[col_])
        return ts, tp

    residual_arr = np.full(n_res, np.inf, dtype=np.float64)
    res_norm = float("inf")
    n_iter = 0

    cur = _evaluate(state0, period)
    if cur is not None:
        residual_arr, res_norm, _ap, _af = cur

    for n_iter in range(1, max_iter + 1):  # noqa: B007
        if cur is None:
            break
        residual_arr, res_norm, arc_point, arc_full = cur
        if res_norm < tol:
            break

        t_point = _event_time(period)
        arc_q = arc_full if arc_full is not None else arc_point
        jac = np.zeros((n_res, n_free), dtype=np.float64)

        # Point rows: d(X(t_point) - X0)/d(unknown). The period column needs
        # dX/dt at t_point; use ``independent_eom`` when supplied, else a short
        # finite-difference propagation as a fallback.
        if independent_eom is not None:
            f_point = independent_eom(t_point, arc_point.state_f)
        else:
            f_point = _state_deriv_fd(propagate_fn, arc_point.state_f, t_point)
        t_scale = 0.5 if is_half_period_residual else 1.0
        for col, unknown in enumerate(free_vars):
            if unknown == IDX_T:
                for row, ridx in enumerate(point_residual_indices):
                    jac[row, col] = float(f_point[ridx]) * t_scale
            else:
                for row, ridx in enumerate(point_residual_indices):
                    val = float(arc_point.stm[ridx, unknown])
                    if ridx == unknown:
                        val -= 1.0
                    jac[row, col] = val

        # Integral rows: d(q_i)/d(unknown). For state cols: dq_dx0 (full
        # period). For the T col: h_i(X_f, T) * t_scale_q where the integral is
        # over the full period (t_scale_q = 1.0 regardless of the point event).
        if n_q > 0:
            t_full = period
            xf_full = arc_q.state_f
            for i, c in enumerate(integral_constraints):
                row = n_point + i
                for col, unknown in enumerate(free_vars):
                    if unknown == IDX_T:
                        h_end = float(c.integrand(t_full, xf_full, integrand_params))
                        jac[row, col] = h_end * weights[i]
                    else:
                        jac[row, col] = float(arc_q.dq_dx0[i, unknown]) * weights[i]

        period_cap = period_step_cap_frac * abs(period)
        dz = _newton_step(
            jac, residual_arr, free_vars, state_cap=state_step_cap, period_cap=period_cap
        )
        if not np.all(np.isfinite(dz)) or float(np.max(np.abs(dz))) < 1e-15:
            break

        # Backtracking line search on the TOTAL residual norm.
        scale = 1.0
        improved = None
        for _bt in range(max_backtrack + 1):
            trial_state, trial_period = _build_trial(state0, period, dz, scale)
            trial_norm = _residual_only(trial_state, trial_period)
            if math.isfinite(trial_norm) and trial_norm < res_norm:
                trial = _evaluate(trial_state, trial_period)
                if trial is not None:
                    improved = (trial_state, trial_period, trial)
                    break
            scale *= 0.5
        if improved is None:
            break
        state0, period, cur = improved

    # Recover q_values at the final orbit (full period).
    if n_q > 0:
        try:
            final_arc = propagate_fn(state0, period)
            q_values = final_arc.q_values.copy()
        except RuntimeError:
            q_values = np.full(n_q, np.nan, dtype=np.float64)
    else:
        q_values = np.zeros(0, dtype=np.float64)

    point_res = float(np.linalg.norm(residual_arr[:n_point])) if n_point else 0.0
    integral_res = float(np.linalg.norm(residual_arr[n_point:])) if n_q else 0.0

    # Independent point-residual closure (Radau, machine precision).
    independent_residual = float("nan")
    if independent_eom is not None:
        try:
            sol = solve_ivp(
                independent_eom,
                (0.0, period),
                state0,
                method="Radau",
                rtol=max(rtol, 1e-12),
                atol=max(atol, 1e-12),
            )
            if sol.success:
                independent_residual = float(np.linalg.norm(sol.y[:, -1] - state0))
        except (RuntimeError, ValueError):
            pass

    converged_corrector = res_norm < tol
    if independent_eom is not None:
        converged_independent = (
            math.isfinite(independent_residual) and independent_residual < independent_tol
        )
        converged = converged_corrector and converged_independent
    else:
        converged = converged_corrector

    return IntegralCorrectorResult(
        state0=state0,
        period=float(period),
        q_values=q_values,
        point_residual=point_res,
        integral_residual=integral_res,
        total_residual=float(res_norm),
        converged=bool(converged),
        independent_closure_residual=float(independent_residual),
        n_iter=int(n_iter),
    )


def _state_deriv_fd(
    propagate_fn: Callable[[NDArray[np.float64], float], AugmentedArc],
    state_f: NDArray[np.float64],
    t_point: float,
    *,
    dt: float = 1e-6,
) -> NDArray[np.float64]:
    """Finite-difference dX/dt at the event time via two short propagations.

    Only used when no ``independent_eom`` is supplied for the period-column term
    of the point rows. Two-point forward difference around ``state_f``.
    """
    arc_p = propagate_fn(state_f, dt)
    return (arc_p.state_f - state_f) / dt


def _newton_step(
    jac: NDArray[np.float64],
    residual: NDArray[np.float64],
    free_vars: tuple[int, ...],
    *,
    state_cap: float,
    period_cap: float,
) -> NDArray[np.float64]:
    """Damped least-squares Newton step with per-component caps.

    Mirrors :func:`cyclerfinder.search.cr3bp_general_periodic_3d._newton_step`:
    ``solve`` when square + nonsingular, else min-norm ``lstsq``; raw step then
    clipped componentwise (state vs period scale).
    """
    nres, nfree = jac.shape
    if nres == nfree:
        try:
            dz = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jac, -residual, rcond=None)
    else:
        dz, *_ = np.linalg.lstsq(jac, -residual, rcond=None)
    dz_out: NDArray[np.float64] = np.asarray(dz, dtype=np.float64)
    for col, unknown in enumerate(free_vars):
        cap = period_cap if unknown == IDX_T else state_cap
        dz_out[col] = float(np.clip(dz_out[col], -cap, cap))
    return dz_out


# ---------------------------------------------------------------------------
# Constraint factories.
# ---------------------------------------------------------------------------


def time_integral_constraint(target: float, *, weight: float = 1.0) -> IntegralConstraint:
    """Trivial time integral ``q = integral_0^T 1 dt = T`` (steers the period).

    The integrand ``h = 1`` has a non-trivial T-column Jacobian term
    (``h(X_f, T) = 1``) and zero state gradient, so it is the simplest
    constraint that actually drives the corrector — it pins the period to
    ``target``. Useful as a sanity / exercise constraint and as a building
    block for arc-duration-sum constraints.
    """

    def integrand(t: float, state6: NDArray[np.float64], params: Any) -> float:
        return 1.0

    def grad(t: float, state6: NDArray[np.float64], params: Any) -> NDArray[np.float64]:
        return np.zeros(6, dtype=np.float64)

    return IntegralConstraint(
        label="time_integral",
        integrand=integrand,
        integrand_grad=grad,
        target=float(target),
        weight=float(weight),
    )


def sun_commensurate_period_constraint(*_args: Any, **_kwargs: Any) -> IntegralConstraint:
    """Sun-commensurate-period / arc-duration-sum constraint (Step 6, DEFERRED).

    This constraint is part of the deferred multi-arc wiring (#388 routes
    through the heliocentric DSM corrector, not this rotating-frame module; see
    the blueprint's deferred open question 2). Building it correctly requires
    the multi-arc duration bookkeeping that Steps 4-6 introduce, so it is a
    stub here.
    """
    raise NotImplementedError(
        "sun_commensurate_period_constraint is Step 6 (deferred): it needs the "
        "multi-arc duration wiring (#388 heliocentric DSM corrector path). Use "
        "time_integral_constraint for a single-arc period constraint."
    )


# Keep ``field`` imported-usage explicit for ruff (re-exported convenience).
__all__ = [
    "IDX_T",
    "IDX_X",
    "IDX_XDOT",
    "IDX_Y",
    "IDX_YDOT",
    "IDX_Z",
    "IDX_ZDOT",
    "AugmentedArc",
    "IntegralConstraint",
    "IntegralCorrectorResult",
    "correct_with_integral_constraints",
    "propagate_augmented_bcr4bp",
    "propagate_augmented_cr3bp",
    "sun_commensurate_period_constraint",
    "time_integral_constraint",
]
