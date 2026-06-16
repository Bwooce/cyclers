"""BCR4BP periodic-orbit corrector (#292 Phase 1 Part B).

Strict-periodicity Newton corrector for the standard (incoherent) BCR4BP. The
BCR4BP is non-autonomous (the Sun phase advances with time), so strict
periodicity requires the orbit period to commensurate with the Sun's synodic
frequency: ``omega_sun * T = 2 * pi * n`` for some positive integer ``n``.
The corrector enforces this by treating ``T`` as a free variable (or fixing it)
and uses an analytic STM Jacobian on the closure residual ``X(T) - X(0)``.

Design parallels :mod:`cyclerfinder.search.cr3bp_general_periodic_3d` (Phase 1
of the 3D CR3BP capability build, #291). The key differences are:

  * non-autonomous EOM => the period MUST be near a Sun-commensurate value
    (the caller passes ``sun_commensurate_n``);
  * no Jacobi constant => closure is purely on the 6D state (no fixed-energy
    constraint to thread through ydot0);
  * independent closure check re-propagates with ``Radau`` (a different
    integrator) at the converged IC, per the orbit-closure discipline.

What this module DOES
---------------------
  * Symmetric / general single-shooting Newton corrector with configurable
    ``free_vars`` and ``residual_indices`` (the 3D Phase 1 pattern).
  * Sourced golden tests: the Andreu / Rosales-Jorba POL1 IC is used as a
    SEED (not as a closure golden value -- the published IC is for the QBCP,
    the implemented model is the BCR4BP, so the model gap is O(eps^2)).
    The genome tests assert CLOSURE (residual < tol AND independent
    re-propagation < independent_tol), not specific IC values.
  * Independent (Radau) re-propagation closure check on every returned
    orbit, per the orbit-closure discipline.

What this module does NOT
-------------------------
  * Family continuation (Phase 2; natural-parameter / pseudo-arclength along
    POL1 or mu_sun continuation from CR3BP).
  * Halo-BCR4BP extension (Phase 3).
  * KNOWN_CORPUS / CandidateSignature widening (Phase 4).
  * V-pipeline adaptation / catalogue admission (Phase 5).

Discipline
----------
  * Sourced goldens only (`feedback_golden_tests_sourced_only`).
  * Independent cross-check mandatory (`feedback_orbit_closure_discipline`).
  * No catalogue writeback. No novelty claims.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.search.outcome_log import log_outcome

# State-vector index aliases (full 6D rotating-frame state). Mirror the
# conventions in :mod:`cyclerfinder.search.cr3bp_general_periodic_3d`.
IDX_X = 0
IDX_Y = 1
IDX_Z = 2
IDX_XDOT = 3
IDX_YDOT = 4
IDX_ZDOT = 5
IDX_T = 6  # period in 7th slot when free


@dataclass(frozen=True)
class BCR4BPPeriodicOrbit:
    """A converged (or attempted) BCR4BP periodic orbit.

    Attributes
    ----------
    state_initial :
        Corrected initial state (6 components, rotating frame).
    period_nondim :
        Full nondim period T at the converged IC.
    sun_commensurate_n :
        Commensurability number: ``omega_sun * period = 2*pi * n`` (or as
        close as the corrector landed). Strict periodicity in the BCR4BP
        requires an integer ``n``.
    sun_phase_drift :
        ``|omega_sun * period - 2*pi * n|`` at the converged IC. Zero only
        if the period is exactly commensurate; small nonzero values mean
        the corrector picked a nearby NON-periodic invariant curve (the
        flow is quasi-periodic when ``T`` and ``2*pi/omega_sun`` are
        incommensurate).
    converged :
        True iff the corrector residual is below ``tol`` AND the independent
        Radau closure residual is below ``independent_tol``.
    corrector_residual :
        L2 norm of the masked closure residual at the corrector's return time
        under DOP853.
    independent_closure_residual :
        L2 norm of the FULL 6D state difference ``X(period) - X(0)`` from a
        re-propagation under Radau. Cross-integrator check per the orbit-
        closure discipline.
    n_iter :
        Newton iterations consumed.
    system :
        The BCR4BPSystem the orbit was corrected in.
    free_vars / residual_indices / is_half_period_residual :
        Echoed back so the caller knows what was actually optimised.
    notes :
        Free-form caller-supplied annotation.
    """

    state_initial: NDArray[np.float64]
    period_nondim: float
    sun_commensurate_n: int
    sun_phase_drift: float
    converged: bool
    corrector_residual: float
    independent_closure_residual: float
    n_iter: int
    system: bcr4bp.BCR4BPSystem
    free_vars: tuple[int, ...]
    residual_indices: tuple[int, ...]
    is_half_period_residual: bool
    notes: str = ""


def _propagate_with_stm(
    system: bcr4bp.BCR4BPSystem,
    state0: NDArray[np.float64],
    t: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Wrap :func:`bcr4bp.propagate_bcr4bp` with STM; returns ``(state_f, STM)``."""
    arc = bcr4bp.propagate_bcr4bp(system, state0, t, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    return arc.state_f, arc.stm


def _newton_step(
    jac: NDArray[np.float64],
    residual: NDArray[np.float64],
    free_vars: tuple[int, ...],
    *,
    state_cap: float,
    period_cap: float,
) -> NDArray[np.float64]:
    """Compute a damped Newton step on the (n_res x n_free) Jacobian.

    Uses ``solve`` when square + nonsingular, else ``lstsq`` (min-norm). Each
    component is then clipped: state components to ``[-state_cap, state_cap]``,
    the period (if in ``free_vars``) to ``[-period_cap, period_cap]``.
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


def correct_bcr4bp_periodic(
    system: bcr4bp.BCR4BPSystem,
    state_guess: NDArray[np.float64] | Sequence[float],
    period_guess: float,
    *,
    sun_commensurate_n: int,
    free_vars: tuple[int, ...] = (IDX_X, IDX_YDOT, IDX_T),
    residual_indices: tuple[int, ...] = (IDX_Y, IDX_XDOT, IDX_ZDOT),
    is_half_period_residual: bool = False,
    tol: float = 1e-10,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    state_step_cap: float = 0.1,
    period_step_cap_frac: float = 0.2,
    independent_tol: float = 1e-6,
    max_backtrack: int = 20,
    require_monotone_decrease: bool = True,
    notes: str = "",
) -> BCR4BPPeriodicOrbit:
    """Symmetric / general BCR4BP single-shooting periodic-orbit corrector.

    The Newton system is constructed on the closure residual at the return
    time ``T_event = T`` (full period) or ``T/2`` (perpendicular-crossing
    half period). The Jacobian is built from the BCR4BP variational STM (6x6)
    plus, if ``T`` is free, a column from the EOM evaluated at the endpoint
    (chain-rule on the integration horizon).

    Parameters
    ----------
    system :
        BCR4BPSystem (see :class:`cyclerfinder.core.bcr4bp.BCR4BPSystem`).
    state_guess :
        6-vector IC. Free-var mask selects which components the corrector
        updates; others are held FIXED.
    period_guess :
        Initial period guess. For strict BCR4BP periodicity this should be
        near ``sun_commensurate_period(omega_sun, sun_commensurate_n)``.
    sun_commensurate_n :
        The integer commensurability number tracked alongside the corrector;
        ECHOED into the result. Pass the closest integer ``n`` to the seed
        period: ``round(omega_sun * period_guess / (2 * pi))``. The corrector
        is also queried for the actual ``sun_phase_drift`` at the converged
        period -- small for a genuine commensurate periodic orbit; large
        means the orbit is on an invariant curve (quasi-periodic), not a
        strict periodic orbit.
    free_vars :
        Indices of free unknowns. Default ``(IDX_X, IDX_YDOT, IDX_T)`` is
        the NRHO-style symmetric pattern: x0 and ydot0 of a perpendicular-
        crossing IC plus the period. Use ``(IDX_X, IDX_Y, IDX_Z, IDX_XDOT,
        IDX_YDOT, IDX_ZDOT, IDX_T)`` for the full asymmetric mode.
    residual_indices :
        Indices in the 6D state whose endpoint values must match the IC.
        Default ``(IDX_Y, IDX_XDOT, IDX_ZDOT)`` is the perpendicular-crossing
        closure at ``T/2`` (paired with ``is_half_period_residual=True``).
        Use ``(IDX_X, ..., IDX_ZDOT)`` for the full asymmetric closure at
        ``T``.
    is_half_period_residual :
        True iff the residual is evaluated at ``T/2`` (symmetric tulip/NRHO
        case); False iff at ``T`` (full closure).
    tol :
        Newton convergence tolerance (L2 norm of the residual vector).
    max_iter :
        Newton iteration cap.
    rtol, atol :
        Integrator tolerances inside the Newton loop (DOP853).
    state_step_cap, period_step_cap_frac :
        Per-component caps on the Newton step.
    max_backtrack :
        Backtracking line-search iterations on ``||R||``.
    require_monotone_decrease :
        If True (default) the step is rejected unless ``||R||`` strictly
        decreases; ill-conditioned BCR4BP arcs benefit from this. Set False
        for very rough seeds where transient increases are expected (e.g.
        QBCP -> BCR4BP model-gap convergence from a QBCP IC).
    independent_tol :
        Tolerance for the Radau closure cross-check.

    Returns
    -------
    BCR4BPPeriodicOrbit
        See dataclass docstring. ``converged`` is True iff BOTH the corrector
        residual AND the independent (Radau) closure check pass.

    Notes
    -----
    A bogus seed (e.g. inside a primary, or far from any periodic basin)
    typically fails to converge: the Newton iteration either stalls or the
    backtracking line search exhausts its budget. The result then has
    ``converged = False`` and the caller can act on the structured failure
    rather than getting a silently-bad periodic-orbit claim.
    """
    free_vars = tuple(sorted(set(int(v) for v in free_vars)))
    residual_indices = tuple(sorted(set(int(v) for v in residual_indices)))
    if not all(0 <= v <= IDX_T for v in free_vars):
        raise ValueError(f"free_vars must be in [0,6]; got {free_vars}")
    if not all(0 <= v <= IDX_ZDOT for v in residual_indices):
        raise ValueError(f"residual_indices must be in [0,5]; got {residual_indices}")
    if not free_vars or not residual_indices:
        raise ValueError("free_vars and residual_indices must each be non-empty")
    if sun_commensurate_n <= 0:
        raise ValueError(f"sun_commensurate_n must be positive int; got {sun_commensurate_n}")

    state0 = np.asarray(state_guess, dtype=np.float64).copy()
    if state0.shape != (6,):
        raise ValueError(f"state_guess must have shape (6,); got {state0.shape}")
    period = float(period_guess)
    if period <= 0.0:
        raise ValueError(f"period_guess must be > 0; got {period}")

    n_res = len(residual_indices)
    n_free = len(free_vars)
    residual_arr = np.full(n_res, np.inf, dtype=np.float64)
    res_norm = float("inf")
    n_iter = 0

    def _evaluate(
        state_in: NDArray[np.float64], period_in: float
    ) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]] | None:
        """Propagate, build the residual + endpoint state + STM.

        Returns ``(residual_vec, residual_norm, state_f, stm)`` or None on
        integrator failure.
        """
        t_ev = 0.5 * period_in if is_half_period_residual else period_in
        try:
            sf, stm_local = _propagate_with_stm(system, state_in, t_ev, rtol=rtol, atol=atol)
        except RuntimeError:
            return None
        d = sf - state_in
        res = d[list(residual_indices)].astype(np.float64, copy=True)
        if not np.all(np.isfinite(res)):
            return None
        return res, float(np.linalg.norm(res)), sf, stm_local

    def _residual_only(state_in: NDArray[np.float64], period_in: float) -> float:
        """Cheap residual norm (no STM) for the line search."""
        t_ev = 0.5 * period_in if is_half_period_residual else period_in
        try:
            arc = bcr4bp.propagate_bcr4bp(
                system, state_in, t_ev, with_stm=False, rtol=rtol, atol=atol
            )
        except RuntimeError:
            return float("inf")
        d = arc.state_f - state_in
        return float(np.linalg.norm(d[list(residual_indices)]))

    def _build_trial(
        cur_state: NDArray[np.float64],
        cur_period: float,
        step: NDArray[np.float64],
        scale_: float,
    ) -> tuple[NDArray[np.float64], float]:
        """Apply a scaled Newton step to ``(cur_state, cur_period)``."""
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

    cur = _evaluate(state0, period)
    if cur is not None:
        residual_arr, res_norm, _state_f, _stm = cur

    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as iter count
        if cur is None:
            break
        residual_arr, res_norm, state_f, stm = cur
        if res_norm < tol:
            break

        # Build the Jacobian: J[i,j] = d(residual_i) / d(unknown_j).
        # State columns: d(diff_i)/d(state0_j) = STM[i,j] - delta_{i,j}.
        # Period column: d(diff_i)/dT = f_end_i * t_scale (half-period chain).
        t_event = 0.5 * period if is_half_period_residual else period
        t_scale = 0.5 if is_half_period_residual else 1.0
        jac = np.zeros((n_res, n_free), dtype=np.float64)
        f_end = bcr4bp.bcr4bp_eom(t_event, state_f, system)
        for col, unknown in enumerate(free_vars):
            if unknown == IDX_T:
                for row, ridx in enumerate(residual_indices):
                    jac[row, col] = float(f_end[ridx]) * t_scale
            else:
                for row, ridx in enumerate(residual_indices):
                    val = float(stm[ridx, unknown])
                    if ridx == unknown:
                        val -= 1.0
                    jac[row, col] = val

        period_cap = period_step_cap_frac * abs(period)
        dz = _newton_step(
            jac, residual_arr, free_vars, state_cap=state_step_cap, period_cap=period_cap
        )
        if not np.all(np.isfinite(dz)) or float(np.max(np.abs(dz))) < 1e-15:
            break

        if require_monotone_decrease:
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
        else:
            # Damped Newton without strict monotone decrease (rare for BCR4BP).
            new_state, new_period = _build_trial(state0, period, dz, 1.0)
            new_cur = _evaluate(new_state, new_period)
            if new_cur is None:
                break
            state0, period, cur = new_state, new_period, new_cur

    # Sun-phase drift at the converged period: |omega_sun * T - 2*pi*n|.
    sun_phase_drift = abs(system.omega_sun_nondim * period - 2.0 * math.pi * sun_commensurate_n)

    # Independent (Radau) full-period closure cross-check.
    independent_residual = float("nan")
    try:
        sol = solve_ivp(
            bcr4bp.bcr4bp_eom,
            (0.0, period),
            state0,
            args=(system,),
            method="Radau",
            rtol=max(rtol, 1e-12),
            atol=max(atol, 1e-12),
        )
        if sol.success:
            independent_residual = float(np.linalg.norm(sol.y[:, -1] - state0))
    except (RuntimeError, ValueError):
        pass

    converged_corrector = res_norm < tol
    converged_independent = (
        math.isfinite(independent_residual) and independent_residual < independent_tol
    )
    converged = converged_corrector and converged_independent

    log_outcome(
        solver="bcr4bp.correct_bcr4bp_periodic",
        inputs={
            "state_guess": np.asarray(state_guess, dtype=np.float64).tolist(),
            "period_guess": float(period_guess),
            "sun_commensurate_n": int(sun_commensurate_n),
            "free_vars": list(free_vars),
            "residual_indices": list(residual_indices),
            "is_half_period_residual": bool(is_half_period_residual),
            "mu": float(system.mu),
            "mu_sun": float(system.mu_sun),
        },
        outcome={
            "converged": bool(converged),
            "corrector_residual": float(res_norm),
            "independent_closure_residual": float(independent_residual),
            "sun_phase_drift": float(sun_phase_drift),
            "state_initial": state0.tolist(),
            "period_nondim": float(period),
            "n_iter": int(n_iter),
        },
        meta={"model": "BCR4BP", "notes": notes},
    )

    return BCR4BPPeriodicOrbit(
        state_initial=state0,
        period_nondim=float(period),
        sun_commensurate_n=int(sun_commensurate_n),
        sun_phase_drift=float(sun_phase_drift),
        converged=bool(converged),
        corrector_residual=float(res_norm),
        independent_closure_residual=float(independent_residual),
        n_iter=int(n_iter),
        system=system,
        free_vars=free_vars,
        residual_indices=residual_indices,
        is_half_period_residual=bool(is_half_period_residual),
        notes=notes,
    )
