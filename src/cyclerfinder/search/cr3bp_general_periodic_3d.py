"""Full-3D broken-plane CR3BP periodic-orbit corrector (#291 Phase 1).

Phase 1 of the Track-A 3D capability build (#291; rank-1 frontier from #286).
This module DROPS the planar restriction baked into
:mod:`cyclerfinder.search.cr3bp_general_periodic` -- whose IC is hardcoded as
``(x0, 0, 0, xdot0, ydot0, 0)`` (z=0, zdot=0) at line 437 of that module -- and
exposes a generic single-shooting corrector with CONFIGURABLE free variables and
CONFIGURABLE closure residuals on the FULL 6D state.

Background (from the #287 spike, commit ``9068aa0``)
----------------------------------------------------
The spike showed that the existing infrastructure was already 3D-capable below
the corrector layer:

  * :func:`cyclerfinder.core.cr3bp.cr3bp_eom` is 6D and has been since spec 2026-
    06-10 (z and zdot appear in the EOM lines 73-81).
  * :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom` computes ``uzz``, ``uxz``, and
    ``uyz`` (lines 114, 116, 117) -- the variational equations are fully
    coupled in z.
  * :func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho` is
    already 3D-aware for the perpendicular-x-z-plane-crossing IC family
    ``(x0, 0, z0, 0, ydot0, 0)``.

The planar restriction lived ONLY in the asymmetric general corrector. The
#287 spike successfully closed a non-trivial 3D Braik-Ross (1,1) family member
at z0 = -0.241 (nondim) with corrector residual ~1e-13 and independent closure
~1.4e-10 by feeding ``correct_symmetric_nrho`` a non-zero ``z0_guess``; an
80-member family was mapped via natural-parameter continuation in z0 (see
``data/spike_287.jsonl``, ``scripts/spike_287_3d_aldrin.py``).

Phase 1 deliverable
-------------------
A 3D-aware corrector that callers (e.g. the #284 asymmetric scan, future
broken-plane campaigns, and Phase 2's 3D family-tracer) can invoke. Two usage
patterns are supported by a single function:

  1. **Symmetric 3D (tulip / NRHO style)** -- free vars ``(z0, ydot0, T)``,
     residual ``(y, xdot, zdot)`` at ``T/2``. This recovers the
     ``correct_symmetric_nrho`` behaviour with a generic API.
  2. **Full 3D asymmetric (broken-plane)** -- free vars ``(x0, y0, z0, xdot0,
     ydot0, zdot0, T)``, residual ``(x-x0, y-y0, z-z0, xdot-xdot0,
     ydot-ydot0, zdot-zdot0)`` at ``T``. 6 equations in 7 unknowns;
     min-norm least-squares Newton step.

The math is single-shooting Newton on the STM, identical in structure to the
existing :func:`cyclerfinder.search.cr3bp_periodic.correct_periodic` (which
already handles 6D + T), but with the free-var / residual-index masks layered
on top so callers can ALSO express partial-symmetric closures cleanly.

Independent closure check
-------------------------
Every converged orbit is RE-PROPAGATED at the converged IC with ``Radau`` at
``rtol = atol = 1e-12``. The L2 norm of ``X(T) - X(0)`` (or of the requested
residual subset at ``T_event``) is returned as ``independent_closure_residual``
-- a different integrator than the DOP853 in the Newton loop, gating against
single-integrator artefacts. Per ``feedback_orbit_closure_discipline``: a
clean residual at the corrector layer is NEVER trusted alone.

What Phase 1 does NOT deliver
-----------------------------
This module is the corrector. Phase 2-N will add:

  * 3D family-tracer integration (likely
    :func:`cyclerfinder.search.nrho_continuation.continue_nrho_family` widened
    to broken-plane stepping in a chosen continuation parameter).
  * Bifurcation classification UI for genuinely 3D bifurcation brackets.
  * CandidateSignature / KNOWN_CORPUS widening for 3D orbit fingerprints.
  * V-pipeline updates for 3D acceptance criteria.
  * Sourced same-model goldens for new 3D family discoveries.

The honest limitation surfaced by the spike: at small ``z0_guess`` (|z0| < 0.01
in the Earth-Moon (1,1) family) the corrector collapses to the planar member,
returning ``z0 ~ 0``. This is correct -- the planar manifold is dynamically
invariant -- but a caller searching for "genuinely 3D" orbits must screen on
``|z0| > eps`` after correction (see the ``degenerate_planar`` flag in
:class:`Periodic3DOrbit`).

Discipline
----------
  * NO catalogue writeback (Phase 1 is corrector + tests + doc only).
  * NO novelty claims (the spike's 3D family is likely-rediscovery of
    Antoniadou-Voyatzis 2018, arXiv:1811.09442; KNOWN_CORPUS gained that
    anchor at commit ``568d8a4``).
  * Sourced-only golden tests: the test EXPECTED states trace to
    ``data/catalogue.yaml`` row ``braik-ross-c11a-cycler-2026``. The 3D
    extension's correctness is verified by TOPOLOGY (z!=0 confirmed, period
    range plausible) + CLOSURE (independent re-propagation ~< 1e-9), not by
    a number our own code computed (which would be circular per
    ``feedback_golden_tests_sourced_only``).

References
----------
  * #286 frontier-scoping doc:
    ``docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md``.
  * #287 scoping spike:
    ``data/spike_287.jsonl`` + ``scripts/spike_287_3d_aldrin.py``;
    verdict in ``docs/notes/2026-06-16-287-3d-spike-verdict.md``.
  * Existing planar asymmetric corrector (kept for backward compat):
    :mod:`cyclerfinder.search.cr3bp_general_periodic`.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.outcome_log import log_outcome

# State-vector index aliases (full 6D rotating-frame state).
IDX_X = 0
IDX_Y = 1
IDX_Z = 2
IDX_XDOT = 3
IDX_YDOT = 4
IDX_ZDOT = 5
IDX_T = 6  # the period is the 7th unknown when ``T`` is in ``free_vars``.

# Convenience free-var bundles for the two most common usage patterns.
FREE_VARS_FULL_ASYMMETRIC: tuple[int, ...] = (
    IDX_X,
    IDX_Y,
    IDX_Z,
    IDX_XDOT,
    IDX_YDOT,
    IDX_ZDOT,
    IDX_T,
)
"""All 6 IC components + period free. 6 closure residuals, 7 unknowns.

Use for full asymmetric 3D broken-plane orbits where no symmetry pins any
component. The min-norm least-squares Newton step lands on the closest periodic
orbit in IC-space.
"""

FREE_VARS_SYMMETRIC_TULIP: tuple[int, ...] = (IDX_Z, IDX_YDOT, IDX_T)
"""Symmetric x-z-plane-crossing free vars: (z0, ydot0, T).

Use for tulip / NRHO-style orbits with IC ``(x0, 0, z0, 0, ydot0, 0)`` where
``x0`` is FIXED (the family-curve parameter) and the half-period perpendicular
crossing residual is ``(y, xdot, zdot)`` at ``T/2``. Equivalent to
:func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`'s setup.
"""

# Default closure residuals (full 6D state at the return time).
RESIDUAL_FULL_STATE_AT_T: tuple[int, ...] = (
    IDX_X,
    IDX_Y,
    IDX_Z,
    IDX_XDOT,
    IDX_YDOT,
    IDX_ZDOT,
)
"""All 6 state components at the return time. Default residual.

The residual is ``X(T_event) - X(0)`` masked to these indices.
"""

RESIDUAL_PERPENDICULAR_HALF_PERIOD: tuple[int, ...] = (IDX_Y, IDX_XDOT, IDX_ZDOT)
"""Tulip-style perpendicular-crossing residual: ``(y, xdot, zdot)`` at ``T/2``.

Pairs with ``FREE_VARS_SYMMETRIC_TULIP``.
"""


@dataclass(frozen=True)
class Periodic3DOrbit:
    """A converged (or attempted) full-3D CR3BP periodic orbit.

    Attributes
    ----------
    state0 :
        Corrected initial state (6 components, rotating frame).
    T_TU :
        Full period in nondim TU.
    jacobi :
        Jacobi constant at the corrected IC.
    converged :
        True iff the corrector residual is below ``tol`` AND the independent
        closure check is below ``independent_tol`` (default 1e-6).
    corrector_residual :
        L2 norm of the masked closure residual at the corrector's return time
        (DOP853 inside the Newton loop).
    independent_closure_residual :
        L2 norm of the FULL 6D state difference ``X(T) - X(0)`` from a
        re-propagation with ``Radau`` at ``rtol = atol = 1e-12``. A different
        integrator gates against single-integrator artefacts.
    n_iter :
        Newton iterations consumed.
    degenerate_planar :
        True iff ``|z0|`` and ``|zdot0|`` at the corrected IC are both below
        ``1e-9``. The corrector landed on the PLANAR member -- mathematically
        valid (the planar manifold is invariant) but not a genuinely 3D orbit.
        Callers searching for 3D orbits should screen on this flag.
    free_vars :
        The free-variable mask used. Echoed back for downstream consumers.
    residual_indices :
        The residual mask used. Echoed back.
    is_half_period_residual :
        True iff the residual was evaluated at ``T/2`` (perpendicular-crossing
        symmetric closure); False iff at ``T`` (full closure). Affects how
        downstream callers interpret ``corrector_residual``.
    """

    state0: NDArray[np.float64]
    T_TU: float
    jacobi: float
    converged: bool
    corrector_residual: float
    independent_closure_residual: float
    n_iter: int
    degenerate_planar: bool
    free_vars: tuple[int, ...]
    residual_indices: tuple[int, ...]
    is_half_period_residual: bool


def _propagate_with_stm(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Wrap :func:`cr3bp.propagate` for clarity; returns ``(state_f, STM)``."""
    arc = cr3bp.propagate(system, state0, t, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None  # we asked for it
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

    Uses ``solve`` when square + nonsingular, else ``lstsq`` (min-norm). The
    raw step is then clipped componentwise: state components to
    ``[-state_cap, state_cap]``, the period (if in ``free_vars``) to
    ``[-period_cap, period_cap]`` (the period numerical scale is much larger
    than the state-component scale).
    """
    nres, nfree = jac.shape
    if nres == nfree:
        try:
            dz = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jac, -residual, rcond=None)
    else:
        # Under- or over-determined: min-norm least-squares.
        dz, *_ = np.linalg.lstsq(jac, -residual, rcond=None)
    dz_out: NDArray[np.float64] = np.asarray(dz, dtype=np.float64)
    # Per-component cap (state vs period).
    for col, unknown in enumerate(free_vars):
        cap = period_cap if unknown == IDX_T else state_cap
        dz_out[col] = float(np.clip(dz_out[col], -cap, cap))
    return dz_out


def _is_planar(state0: NDArray[np.float64], *, eps: float = 1e-9) -> bool:
    """True iff the IC has |z0| < eps AND |zdot0| < eps (planar collapse)."""
    return bool(abs(float(state0[IDX_Z])) < eps and abs(float(state0[IDX_ZDOT])) < eps)


def correct_general_periodic_3d(
    system: cr3bp.CR3BPSystem,
    state0_guess: NDArray[np.float64] | Sequence[float],
    period_guess: float,
    *,
    free_vars: tuple[int, ...] = FREE_VARS_FULL_ASYMMETRIC,
    residual_indices: tuple[int, ...] = RESIDUAL_FULL_STATE_AT_T,
    is_half_period_residual: bool = False,
    tol: float = 1e-10,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    step_cap: float = 0.1,
    period_step_cap_frac: float = 0.2,
    independent_tol: float = 1e-6,
    planar_eps: float = 1e-9,
    max_backtrack: int = 20,
    require_monotone_decrease: bool = False,
) -> Periodic3DOrbit:
    """Full-3D single-shooting CR3BP periodic-orbit corrector.

    Configurable free vars and residual mask allow this one function to handle
    both the perpendicular-crossing symmetric (NRHO/tulip) case and the full
    asymmetric broken-plane case.

    Parameters
    ----------
    system :
        CR3BP system (only ``system.mu`` is used by the EOM).
    state0_guess :
        6-vector IC ``(x, y, z, xdot, ydot, zdot)``. The free-var mask selects
        which of these the corrector is allowed to update; the others are
        held FIXED at the guess.
    period_guess :
        Initial guess for the FULL nondim period. If ``IDX_T`` is in
        ``free_vars`` the corrector also updates ``T``; otherwise the period
        is held fixed.
    free_vars :
        Tuple of indices in ``{0..6}`` listing the unknowns. Index 6 means the
        period ``T`` is free. Default is the full 7-parameter asymmetric set;
        :data:`FREE_VARS_SYMMETRIC_TULIP` is the symmetric-NRHO bundle.
    residual_indices :
        Tuple of indices in ``{0..5}`` listing which state components must
        match at the return time. Default is the full 6D closure. For the
        perpendicular-crossing symmetric case (Y=0 at half-period) the
        natural choice is :data:`RESIDUAL_PERPENDICULAR_HALF_PERIOD`.
    is_half_period_residual :
        If True the residual is evaluated at ``T/2`` (symmetric tulip/NRHO
        case); else at ``T`` (full asymmetric closure). Default False.
    tol :
        Newton convergence tolerance (L2 norm of the residual vector). The
        default 1e-10 reflects the DOP853 integrator floor at rtol=atol=1e-12
        for a full-period 6D closure; the half-period 3-component symmetric
        residual reaches 1e-12 cleanly.
    max_iter :
        Newton iteration cap.
    rtol, atol :
        Integrator tolerances inside the Newton loop (DOP853).
    step_cap :
        Per-component cap on the raw Newton step for STATE unknowns.
        Defensive against wild overshoots from rough seeds.
    period_step_cap_frac :
        Period step cap as a fraction of the current period (default 20%).
        The period is on a much larger numerical scale than the state
        components and benefits from a separate cap.
    max_backtrack :
        Backtracking line-search iterations on ``||R||``. Ill-conditioned 3D
        return arcs make the full Newton step overshoot from a rough seed;
        halve the step until ``||R||`` decreases. Mirrors the pattern in
        :mod:`cyclerfinder.search.cr3bp_general_periodic` (the planar
        asymmetric corrector hit the same ill-conditioning).
    require_monotone_decrease :
        If True, the backtracking line search REJECTS the step unless
        ``||R||`` strictly decreases. If False (default), the capped Newton
        step is always applied -- transient residual increases are tolerated
        (damped Newton can oscillate before settling into the quadratic
        basin). Set True for ill-conditioned long-arc closures (asymmetric
        full-6D-at-T mode); the default False matches
        :func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`,
        which converges cleanly without monotone-decrease enforcement on the
        symmetric half-period 3-residual.
    independent_tol :
        Tolerance for the independent full-period closure check
        (re-propagation under Radau at the same rtol/atol). If the IC closes
        below this threshold AND the corrector residual is below ``tol``, the
        orbit is flagged ``converged = True``. Default 1e-6 is roughly 384 m
        in EM units -- tight enough to reject silent-failure cases but loose
        enough to admit moderately unstable orbits whose chaos amplifies
        round-trip noise.
    planar_eps :
        Threshold for the ``degenerate_planar`` flag. Default 1e-9 nondim
        (~384 m in Earth-Moon distance units).

    Returns
    -------
    Periodic3DOrbit
        See dataclass docstring. ``converged`` is True iff the corrector
        landed AND the independent closure check held.

    Notes
    -----
    The full-asymmetric default (6 residuals, 7 unknowns) is intentionally
    under-determined; the min-norm least-squares step lands on the closest
    periodic orbit in IC-space (a 1-parameter local family). This matches the
    behaviour of :func:`cyclerfinder.search.cr3bp_periodic.correct_periodic`,
    which uses the same min-norm strategy. For a Jacobi-pinned correction the
    caller should compose this corrector with a Jacobi residual on top
    (left to Phase 2 / family-tracer code).

    Multi-shooting is NOT used here; single-shooting handles the spike's basin
    cleanly (residual 1e-13, independent closure 1.4e-10 at z0 = -0.241). Deep
    low-perilune members may eventually require the
    :mod:`cyclerfinder.core.cr3bp_regularized` propagator; that's a Phase 2
    integration, not a Phase 1 dependency.
    """
    free_vars = tuple(sorted(set(int(v) for v in free_vars)))
    residual_indices = tuple(sorted(set(int(v) for v in residual_indices)))
    if not all(0 <= v <= IDX_T for v in free_vars):
        raise ValueError(f"free_vars must be in [0,6]; got {free_vars}")
    if not all(0 <= v <= IDX_ZDOT for v in residual_indices):
        raise ValueError(f"residual_indices must be in [0,5]; got {residual_indices}")
    if not free_vars:
        raise ValueError("free_vars must contain at least one index")
    if not residual_indices:
        raise ValueError("residual_indices must contain at least one index")

    state0 = np.asarray(state0_guess, dtype=np.float64).copy()
    if state0.shape != (6,):
        raise ValueError(f"state0_guess must have shape (6,); got {state0.shape}")
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
        """Propagate, build the residual + Jacobian. Returns ``None`` on failure.

        Returns ``(residual_vec, residual_norm, state_f, stm)``.
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
            arc = cr3bp.propagate(system, state_in, t_ev, with_stm=False, rtol=rtol, atol=atol)
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

    for n_iter in range(1, max_iter + 1):  # noqa: B007
        if cur is None:
            break
        residual_arr, res_norm, state_f, stm = cur
        if res_norm < tol:
            break

        # Build the Jacobian J[i,j] = d(residual_i) / d(unknown_j).
        #
        # For state-component unknowns (j in state_free):
        #   d(diff_i)/d(state0_j) = stm[i,j] - delta_{i,j}
        #   (the stm is dX(t_event)/dX(0); the -delta accounts for the X(0) term
        #    in diff = X(t_event) - X(0).)
        # For the period unknown (j == IDX_T):
        #   d(diff_i)/dT = (df/dt at t_event) * (1.0 if not half-period else 0.5)
        #   where f = cr3bp_eom(t_event, X(t_event), mu); the chain rule on
        #   T -> t_event = scale * T inserts the 0.5 factor for the half-period
        #   residual case.
        t_event = 0.5 * period if is_half_period_residual else period
        jac = np.zeros((n_res, n_free), dtype=np.float64)
        f_end = cr3bp.cr3bp_eom(t_event, state_f, system.mu)
        t_scale = 0.5 if is_half_period_residual else 1.0
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
        dz = _newton_step(jac, residual_arr, free_vars, state_cap=step_cap, period_cap=period_cap)
        if not np.all(np.isfinite(dz)) or float(np.max(np.abs(dz))) < 1e-15:
            break

        # Apply the (capped) Newton step. When ``require_monotone_decrease`` is
        # True the step is rejected unless ``||R||`` strictly decreases (halve
        # until it does); else apply the full capped step blind (matches
        # :func:`correct_symmetric_nrho`'s damped-Newton pattern, which
        # converges cleanly through transient residual oscillation).
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
            # Blind step (matches `correct_symmetric_nrho`'s pattern).
            new_state, new_period = _build_trial(state0, period, dz, 1.0)
            new_cur = _evaluate(new_state, new_period)
            if new_cur is None:
                break
            state0, period, cur = new_state, new_period, new_cur

    # Final state derived (no further integration).
    jacobi = float(cr3bp.jacobi_constant(state0, system.mu))
    is_planar = _is_planar(state0, eps=planar_eps)

    # Independent full-period closure check (Radau, machine precision).
    independent_residual = float("nan")
    try:
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period),
            state0,
            args=(system.mu,),
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

    # Passive training-data capture (#210): NO-OP unless CYCLERFINDER_OUTCOME_LOG.
    log_outcome(
        solver="cr3bp.correct_general_periodic_3d",
        inputs={
            "state0_guess": np.asarray(state0_guess, dtype=np.float64).tolist(),
            "period_guess": float(period_guess),
            "free_vars": list(free_vars),
            "residual_indices": list(residual_indices),
            "is_half_period_residual": bool(is_half_period_residual),
            "mu": float(system.mu),
        },
        outcome={
            "converged": bool(converged),
            "corrector_residual": float(res_norm),
            "independent_closure_residual": float(independent_residual),
            "state0": state0.tolist(),
            "T_TU": float(period),
            "jacobi": float(jacobi),
            "n_iter": int(n_iter),
            "degenerate_planar": bool(is_planar),
        },
        meta={"primary": system.primary, "secondary": system.secondary},
    )

    return Periodic3DOrbit(
        state0=state0,
        T_TU=float(period),
        jacobi=float(jacobi),
        converged=bool(converged),
        corrector_residual=float(res_norm),
        independent_closure_residual=float(independent_residual),
        n_iter=int(n_iter),
        degenerate_planar=bool(is_planar),
        free_vars=free_vars,
        residual_indices=residual_indices,
        is_half_period_residual=bool(is_half_period_residual),
    )
