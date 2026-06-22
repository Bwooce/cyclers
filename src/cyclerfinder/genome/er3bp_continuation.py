"""ER3BP Family Continuation.

Continues a periodic orbit family in the Elliptic Restricted 3-Body Problem
along the eccentricity parameter `e`.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_periodic import (
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    ER3BPPeriodicOrbit,
    correct_er3bp_periodic,
)


class ContinuationError(Exception):
    """Raised when family continuation fails to converge at a step."""


def continue_er3bp_family_in_e(
    sys_base: ER3BPSystem,
    seed_state: NDArray[np.float64],
    period_f: float,
    e_target: float,
    n_steps: int,
    *,
    is_half_period_residual: bool = True,
    tol: float = 1e-10,
) -> list[ER3BPPeriodicOrbit]:
    """Continue an ER3BP periodic orbit family in eccentricity.

    Starting from the initial eccentricity `sys_base.e` (usually 0.0) with
    the provided `seed_state`, steps towards `e_target` using `n_steps` steps.

    Args:
        sys_base: Starting ER3BP system definition.
        seed_state: Converged or near-converged 6D state at sys_base.e.
        period_f: Fixed integration time in true anomaly `f`.
        e_target: Target eccentricity to continue towards.
        n_steps: Number of steps (inclusive of the target).
        is_half_period_residual: Symmetry flag for the corrector.
        tol: Convergence tolerance for the corrector.

    Returns:
        List of converged ER3BPPeriodicOrbit objects along the continuation path.
    """
    history: list[ER3BPPeriodicOrbit] = []

    current_e = sys_base.e
    e_step = (e_target - current_e) / n_steps

    # First, converge the seed exactly at the starting eccentricity
    try:
        current_orbit = correct_er3bp_periodic(
            system=sys_base,
            state_guess=seed_state,
            period_f=period_f,
            is_half_period_residual=is_half_period_residual,
            free_vars=(IDX_X, IDX_YDOT),
            residual_indices=(IDX_Y, IDX_XDOT),
            tol=tol,
        )
    except Exception as e:
        raise ContinuationError(f"Failed to converge initial seed at e={current_e}: {e}") from e

    history.append(current_orbit)

    # Secant predictor states
    prev_state = current_orbit.state0
    prev_e = current_e

    for i in range(1, n_steps + 1):
        target_e = sys_base.e + i * e_step
        next_sys = ER3BPSystem(
            mu=sys_base.mu,
            e=target_e,
            primary_name=sys_base.primary_name,
            secondary_name=sys_base.secondary_name,
        )

        # Predictor (Secant if we have at least 2 points, otherwise zeroth-order)
        if len(history) >= 2:
            e_diff = current_e - prev_e
            if abs(e_diff) > 1e-14:
                dstate_de = (current_orbit.state0 - prev_state) / e_diff
                state_guess = current_orbit.state0 + dstate_de * (target_e - current_e)
            else:
                state_guess = current_orbit.state0.copy()
        else:
            state_guess = current_orbit.state0.copy()

        try:
            converged = correct_er3bp_periodic(
                system=next_sys,
                state_guess=state_guess,
                period_f=period_f,
                is_half_period_residual=is_half_period_residual,
                free_vars=(IDX_X, IDX_YDOT),
                residual_indices=(IDX_Y, IDX_XDOT),
                tol=tol,
            )
        except Exception as e:
            raise ContinuationError(
                f"Continuation failed at e={target_e:.5f} (step {i}/{n_steps}): {e}"
            ) from e

        prev_state = current_orbit.state0
        prev_e = current_e

        current_orbit = converged
        current_e = target_e
        history.append(current_orbit)

    return history
