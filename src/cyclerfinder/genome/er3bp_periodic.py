"""ER3BP periodic-orbit corrector (Axis-1).

Strict-periodicity Newton corrector for the Elliptic Restricted 3-Body Problem.
The ER3BP is non-autonomous (the true anomaly `f` explicitly appears in the EOMs),
so strict periodicity requires the orbit period in `f` to be a multiple of 2*pi.
The corrector enforces this by fixing `f_f` (the integration duration in `f`)
and treating only the state components as free variables.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.er3bp as er3bp
from cyclerfinder.search.outcome_log import log_outcome

IDX_X = 0
IDX_Y = 1
IDX_Z = 2
IDX_XDOT = 3
IDX_YDOT = 4
IDX_ZDOT = 5


@dataclass(frozen=True)
class ER3BPPeriodicOrbit:
    """A converged ER3BP periodic orbit."""

    state0: NDArray[np.float64]
    """Converged initial state [x, y, z, x', y', z'] at f=0."""

    period_f: float
    """The converged full period in true anomaly `f`. Must be 2*pi*n."""

    mu: float
    """Mass parameter."""

    e: float
    """Eccentricity of the primaries."""

    corrector_residual: float
    """The L2 norm of the final residual vector during correction."""

    independent_residual: float
    """The L2 norm of the full-period state discontinuity X(T) - X(0)
    measured using an independent (Radau) re-propagation."""

    iterations: int
    """Newton iterations required to converge."""

    notes: str = ""
    """Diagnostic annotations."""


class ConvergenceError(Exception):
    """Raised when the Newton corrector fails to converge."""


def correct_er3bp_periodic(
    system: er3bp.ER3BPSystem,
    state_guess: NDArray[np.float64] | Sequence[float],
    period_f: float,
    *,
    free_vars: tuple[int, ...] = (IDX_X, IDX_YDOT),
    residual_indices: tuple[int, ...] = (IDX_Y, IDX_XDOT),
    is_half_period_residual: bool = True,
    tol: float = 1e-10,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    state_step_cap: float = 0.1,
    independent_tol: float = 1e-5,
    max_backtrack: int = 20,
    require_monotone_decrease: bool = True,
    notes: str = "",
) -> ER3BPPeriodicOrbit:
    """Symmetric / general ER3BP single-shooting periodic-orbit corrector.

    The period `period_f` is FIXED (must be 2*pi for full period or pi for half).
    No time-variable is included in `free_vars`.

    Args:
        system: The ER3BP system definition (mu, e).
        state_guess: 6D state guess at true anomaly f=0.
        period_f: The integration time in true anomaly. For a symmetric half-period
            correction, this should be `pi`. For a full-period correction, `2*pi`.
        free_vars: Indices of state components to vary.
        residual_indices: Indices of state components to zero out at `period_f`.
            For symmetric half-period: (IDX_Y, IDX_XDOT, IDX_ZDOT).
            For full-period: (IDX_X, IDX_Y, IDX_Z, IDX_XDOT, IDX_YDOT, IDX_ZDOT).
            (Note: Z and ZDOT can be omitted for planar orbits).
        is_half_period_residual: If True, tests X(T) = 0 for residual indices.
            If False, tests X(T) - X(0) = 0 for residual indices.
        tol: Convergence tolerance for the residual vector L2 norm.
        ...
    """
    if len(free_vars) != len(residual_indices):
        raise ValueError("Must have equal number of free vars and residual indices.")

    state = np.array(state_guess, dtype=np.float64)
    if state.shape != (6,):
        raise ValueError("state_guess must be 6D.")

    history = []

    for i in range(max_iter):
        # 1. Propagate state and STM
        _, state_hist, stm_final = er3bp.propagate_er3bp(
            state6=state,
            f_span=(0.0, period_f),
            sys=system,
            rtol=rtol,
            atol=atol,
            with_stm=True,
        )

        state_final = state_hist[:, -1]

        # 2. Build residual vector F(X)
        res_list = []
        for ri in residual_indices:
            if is_half_period_residual:
                res_list.append(state_final[ri])
            else:
                res_list.append(state_final[ri] - state[ri])

        f_mat = np.array(res_list)
        err = float(np.linalg.norm(f_mat))
        history.append(err)

        if err < tol:
            break

        # 3. Build Jacobian DF(X)
        # df_mat_ij = d(residual_i) / d(free_var_j)
        n_vars = len(free_vars)
        df_mat = np.zeros((n_vars, n_vars))
        for r_idx, ri in enumerate(residual_indices):
            for v_idx, vi in enumerate(free_vars):
                df_dv = stm_final[ri, vi]
                if not is_half_period_residual and ri == vi:
                    df_dv -= 1.0
                df_mat[r_idx, v_idx] = df_dv

        # 4. Newton step
        try:
            delta = np.linalg.solve(df_mat, -f_mat)
        except np.linalg.LinAlgError as e:
            raise ConvergenceError(f"Singular Jacobian at iter {i}") from e

        # Apply cap to state variables
        delta = np.clip(delta, -state_step_cap, state_step_cap)

        # Backtracking line search
        alpha = 1.0
        best_state = state.copy()
        step_accepted = False

        for _backtrack in range(max_backtrack):
            candidate_state = state.copy()
            for v_idx, vi in enumerate(free_vars):
                candidate_state[vi] += alpha * delta[v_idx]

            # Quick re-propagate candidate (without STM)
            _, cand_hist, _ = er3bp.propagate_er3bp(
                state6=candidate_state,
                f_span=(0.0, period_f),
                sys=system,
                rtol=rtol,
                atol=atol,
                with_stm=False,
            )
            cand_final = cand_hist[:, -1]

            cand_res = []
            for ri in residual_indices:
                if is_half_period_residual:
                    cand_res.append(cand_final[ri])
                else:
                    cand_res.append(cand_final[ri] - candidate_state[ri])
            cand_err = float(np.linalg.norm(cand_res))

            if not require_monotone_decrease or cand_err < err:
                best_state = candidate_state
                step_accepted = True
                break

            alpha *= 0.5

        if not step_accepted and require_monotone_decrease:
            raise ConvergenceError(
                f"Line search failed to decrease error from {err:.2e} at iter {i}"
            )

        state = best_state

    else:
        raise ConvergenceError(f"Failed to converge in {max_iter} iterations. Final err: {err:.2e}")

    # Independent closure check (full period)
    full_period = period_f * 2.0 if is_half_period_residual else period_f

    radau_sol = solve_ivp(
        fun=er3bp.er3bp_eom,
        t_span=(0.0, full_period),
        y0=state,
        args=(system.mu, system.e),
        method="Radau",
        rtol=1e-12,
        atol=1e-12,
    )
    final_radau = radau_sol.y[:, -1]
    independent_err = float(np.linalg.norm(final_radau - state))

    if independent_err > independent_tol:
        log_outcome(
            solver="er3bp_periodic",
            inputs={"mu": system.mu, "e": system.e, "state_guess": state.tolist()},
            outcome={
                "independent_err": independent_err,
                "independent_tol": independent_tol,
                "status": "WARN",
            },
            meta={"message": "Converged but independent closure > tol"},
        )

    return ER3BPPeriodicOrbit(
        state0=state,
        period_f=full_period,
        mu=system.mu,
        e=system.e,
        corrector_residual=err,
        independent_residual=independent_err,
        iterations=i,
        notes=notes,
    )
