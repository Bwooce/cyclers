"""ER3BP Family Continuation.

Continues a periodic orbit family in the Elliptic Restricted 3-Body Problem
along the eccentricity parameter `e`.

Two continuators are provided:

  * :func:`continue_er3bp_family_in_e` — secant predictor in `e`. Robust on
    smooth families but STALLS at folds (turning points where `de/ds` changes
    sign), because the natural parameter `e` is no longer monotone there and
    the predicted next-`e` system has no nearby family member.
  * :func:`continue_er3bp_family_in_e_arclength` — pseudo-arclength predictor
    that walks THROUGH folds in `e`. Mirrors the proven 3D CR3BP arclength
    walker :func:`cyclerfinder.search.cr3bp_3d_family_tracer.continue_general_3d_family`:
    augment the unknowns with `e`, find the null tangent of the residual
    Jacobian via SVD, step along it by `ds`, then correct back onto BOTH the
    periodicity residual AND the arclength constraint
    `<z - z_pred, tau> = 0`.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp
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


def continue_er3bp_family_in_e_partial(
    sys_base: ER3BPSystem,
    seed_state: NDArray[np.float64],
    period_f: float,
    e_target: float,
    n_steps: int,
    *,
    is_half_period_residual: bool = True,
    tol: float = 1e-10,
) -> tuple[list[ER3BPPeriodicOrbit], float | None]:
    """Continue an ER3BP periodic orbit family in eccentricity, non-raising.

    Mirrors :func:`continue_er3bp_family_in_e` exactly (same zeroth/secant
    predictor and corrector configuration) but, instead of raising
    :class:`ContinuationError` on a failed step, returns the orbits converged so
    far together with the eccentricity at which the step failed.

    Returns:
        ``(orbits, death_e)`` where ``orbits`` is the list of converged orbits
        (including the initial seed) and ``death_e`` is the target eccentricity
        of the first step that failed to converge, or ``None`` if the
        continuation reached ``e_target``.
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
    except Exception:
        return history, current_e

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
        except Exception:
            return history, target_e

        prev_state = current_orbit.state0
        prev_e = current_e

        current_orbit = converged
        current_e = target_e
        history.append(current_orbit)

    return history, None


# ---------------------------------------------------------------------------
# Pseudo-arclength continuation in e (fold-capable).
# ---------------------------------------------------------------------------
#
# Mirrors `continue_general_3d_family` (cr3bp_3d_family_tracer.py). The
# continuation parameter is `e` (the ER3BPSystem pulsating-frame eccentricity).
# Augmented unknowns z = (x0, ydot0, e); residual r = (y, xdot) at `period_f`
# from the symmetric IC (x0, 0, 0, 0, ydot0, 0). Two residual equations in
# three unknowns => a 1-D family curve whose tangent is the SVD null vector.
#
# Per-member `e` is carried natively: ER3BPPeriodicOrbit already stores `.e`
# (it is built from `ER3BPSystem(e=...)`), so the family is returned as a plain
# `list[ER3BPPeriodicOrbit]` with no wrapper needed.

_FREE_VARS = (IDX_X, IDX_YDOT)
_RESIDUAL_INDICES = (IDX_Y, IDX_XDOT)


def _arclength_residual_and_jac(
    z: NDArray[np.float64],
    sys_base: ER3BPSystem,
    period_f: float,
    *,
    is_half_period_residual: bool,
    de_fd: float = 1e-7,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """Return ``(r, J)`` for the symmetric ER3BP residual at ``z=(x0, ydot0, e)``.

    The residual ``r = (y, xdot)`` is measured at ``period_f`` from the IC
    ``(x0, 0, 0, 0, ydot0, 0)`` (``is_half_period_residual=True``) or the
    full-period discontinuity ``X(T)-X(0)`` on the residual indices otherwise.

    The 2x3 Jacobian has columns ``(x0, ydot0, e)``. Sensitivities to the state
    components come from the STM (``stm[ri, vi]``); sensitivity to ``e`` is a
    central finite difference of the residual (the ER3BP STM tracks
    ``dX/dX(0)`` only, not ``dX/de``). Returns ``None`` on integrator failure.
    """
    x0_, ydot0_, e_ = float(z[0]), float(z[1]), float(z[2])
    state = np.array([x0_, 0.0, 0.0, 0.0, ydot0_, 0.0], dtype=np.float64)
    sys_e = ER3BPSystem(
        mu=sys_base.mu,
        e=e_,
        primary_name=sys_base.primary_name,
        secondary_name=sys_base.secondary_name,
    )

    try:
        _, state_hist, stm = propagate_er3bp(
            state6=state,
            f_span=(0.0, period_f),
            sys=sys_e,
            rtol=rtol,
            atol=atol,
            with_stm=True,
        )
    except (RuntimeError, ValueError):
        return None
    state_final = state_hist[:, -1]

    def _residual(sf: NDArray[np.float64], st: NDArray[np.float64]) -> NDArray[np.float64]:
        vals = []
        for ri in _RESIDUAL_INDICES:
            if is_half_period_residual:
                vals.append(float(sf[ri]))
            else:
                vals.append(float(sf[ri] - st[ri]))
        return np.array(vals, dtype=np.float64)

    r = _residual(state_final, state)

    jac = np.zeros((len(_RESIDUAL_INDICES), 3), dtype=np.float64)
    for r_idx, ri in enumerate(_RESIDUAL_INDICES):
        for v_idx, vi in enumerate(_FREE_VARS):
            df_dv = float(stm[ri, vi])
            if not is_half_period_residual and ri == vi:
                df_dv -= 1.0
            jac[r_idx, v_idx] = df_dv

    # Sensitivity to e via central finite difference of the residual.
    r_perturbed = []
    for sign in (+1.0, -1.0):
        sys_pert = ER3BPSystem(
            mu=sys_base.mu,
            e=e_ + sign * de_fd,
            primary_name=sys_base.primary_name,
            secondary_name=sys_base.secondary_name,
        )
        try:
            _, hist_p, _ = propagate_er3bp(
                state6=state,
                f_span=(0.0, period_f),
                sys=sys_pert,
                rtol=rtol,
                atol=atol,
                with_stm=False,
            )
        except (RuntimeError, ValueError):
            return None
        r_perturbed.append(_residual(hist_p[:, -1], state))
    jac[:, 2] = (r_perturbed[0] - r_perturbed[1]) / (2.0 * de_fd)

    return r, jac


def _arclength_tangent(
    jac: NDArray[np.float64],
    prev: NDArray[np.float64] | None,
) -> NDArray[np.float64] | None:
    """Unit null tangent of ``jac`` (last right-singular vector of the SVD).

    Oriented for continuity with ``prev`` (positive dot product). Returns
    ``None`` if the SVD fails or the tangent degenerates.
    """
    try:
        _, _sv, vt = np.linalg.svd(jac, full_matrices=True)
    except np.linalg.LinAlgError:
        return None
    if vt.shape != (3, 3):
        return None
    tau = np.asarray(vt[-1], dtype=np.float64)
    norm = float(np.linalg.norm(tau))
    if not np.all(np.isfinite(tau)) or norm < 1e-12:
        return None
    tau = tau / norm
    if prev is not None and float(np.dot(tau, prev)) < 0.0:
        tau = -tau
    return tau


def _correct_arclength(
    z_pred: NDArray[np.float64],
    tau: NDArray[np.float64],
    sys_base: ER3BPSystem,
    period_f: float,
    *,
    is_half_period_residual: bool,
    tol: float,
    max_iter: int = 60,
    state_cap: float = 0.1,
    e_cap: float = 0.05,
) -> NDArray[np.float64] | None:
    """Newton onto ``{r(z) = 0, tau . (z - z_pred) = 0}``.

    3 unknowns ``z = (x0, ydot0, e)``; 3 equations (2 residual + 1 arclength).
    Returns the converged ``z`` or ``None`` on failure.
    """
    z = z_pred.copy()
    for _ in range(max_iter):
        rj = _arclength_residual_and_jac(
            z,
            sys_base,
            period_f,
            is_half_period_residual=is_half_period_residual,
        )
        if rj is None:
            return None
        r0, grad = rj
        arc = float(np.dot(tau, z - z_pred))
        if float(np.linalg.norm(r0)) < tol and abs(arc) < 1e-12:
            return z
        jmat = np.vstack([grad, tau.reshape(1, 3)])
        rhs = -np.concatenate([r0, np.array([arc])])
        try:
            dz = np.linalg.solve(jmat, rhs)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jmat, rhs, rcond=None)
        dz = np.asarray(dz, dtype=np.float64)
        dz[0] = float(np.clip(dz[0], -state_cap, state_cap))  # x0
        dz[1] = float(np.clip(dz[1], -state_cap, state_cap))  # ydot0
        dz[2] = float(np.clip(dz[2], -e_cap, e_cap))  # e
        z = z + dz
    return None


def continue_er3bp_family_in_e_arclength(
    sys_base: ER3BPSystem,
    seed_state: NDArray[np.float64],
    period_f: float,
    e_target: float,
    *,
    ds: float = 0.005,
    max_steps: int = 400,
    is_half_period_residual: bool = True,
    tol: float = 1e-10,
) -> list[ER3BPPeriodicOrbit]:
    """Pseudo-arclength continuation of an ER3BP family in eccentricity `e`.

    Walks the family curve through folds (turning points in `e` where
    `de/ds` changes sign) that the secant :func:`continue_er3bp_family_in_e`
    stalls at. Mirrors the proven 3D walker
    :func:`cyclerfinder.search.cr3bp_3d_family_tracer.continue_general_3d_family`:
    augmented unknowns ``z = (x0, ydot0, e)``, SVD null-tangent prediction by
    ``ds``, Newton correction onto the periodicity residual plus the arclength
    constraint ``<z - z_pred, tau> = 0``.

    Args:
        sys_base: Starting ER3BP system (its ``.e`` is the seed eccentricity).
        seed_state: Converged or near-converged 6D state at ``sys_base.e``.
        period_f: Fixed integration time in true anomaly `f` (forwarded to the
            corrector exactly as the secant continuator does).
        e_target: Target eccentricity. The walk stops once it reaches (or
            passes) this value, or after ``max_steps`` steps.
        ds: Pseudo-arclength step in the unit ``(x0, ydot0, e)`` tangent space.
        max_steps: Maximum continuation steps.
        is_half_period_residual: Symmetry flag for the corrector / residual.
        tol: Convergence tolerance for the residual L2 norm.

    Returns:
        List of converged :class:`ER3BPPeriodicOrbit` members along the
        family, each carrying its own ``.e``. The first element is the seed
        member at ``sys_base.e``; the last reaches (or brackets) ``e_target``.

    Raises:
        ContinuationError: If the initial seed fails to converge.
    """
    history: list[ER3BPPeriodicOrbit] = []

    try:
        seed_orbit = correct_er3bp_periodic(
            system=sys_base,
            state_guess=seed_state,
            period_f=period_f,
            is_half_period_residual=is_half_period_residual,
            free_vars=_FREE_VARS,
            residual_indices=_RESIDUAL_INDICES,
            tol=tol,
        )
    except Exception as exc:
        raise ContinuationError(
            f"Failed to converge initial seed at e={sys_base.e}: {exc}"
        ) from exc

    history.append(seed_orbit)

    e_seed = sys_base.e
    z_cur = np.array(
        [float(seed_orbit.state0[IDX_X]), float(seed_orbit.state0[IDX_YDOT]), e_seed],
        dtype=np.float64,
    )
    # Orient the first tangent so that e advances toward e_target.
    e_dir = 1.0 if e_target >= e_seed else -1.0
    tau_prev: NDArray[np.float64] | None = None

    for _ in range(max_steps):
        rj = _arclength_residual_and_jac(
            z_cur,
            sys_base,
            period_f,
            is_half_period_residual=is_half_period_residual,
        )
        if rj is None:
            break
        _, jac = rj
        tau = _arclength_tangent(jac, tau_prev)
        if tau is None:
            break
        # On the first step there is no previous tangent: orient by e-direction.
        if tau_prev is None and tau[2] * e_dir < 0.0:
            tau = -tau

        z_pred = z_cur + ds * tau
        z_next = _correct_arclength(
            z_pred,
            tau,
            sys_base,
            period_f,
            is_half_period_residual=is_half_period_residual,
            tol=tol,
        )
        if z_next is None:
            break

        e_next = float(z_next[2])

        # Land exactly on e_target if this arclength step overshot it (only
        # when `e` is locally monotone toward the target — never near a fold,
        # where `e` is no longer advancing toward the target). Re-correct at
        # `e = e_target` with the standard fixed-e corrector seeded from this
        # member, so the final member sits precisely at the requested e.
        overshot = (e_dir > 0.0 and e_next > e_target) or (e_dir < 0.0 and e_next < e_target)
        advancing = (e_dir > 0.0 and e_next > float(z_cur[2])) or (
            e_dir < 0.0 and e_next < float(z_cur[2])
        )
        if overshot and advancing:
            sys_target = ER3BPSystem(
                mu=sys_base.mu,
                e=e_target,
                primary_name=sys_base.primary_name,
                secondary_name=sys_base.secondary_name,
            )
            guess = np.array(
                [float(z_next[0]), 0.0, 0.0, 0.0, float(z_next[1]), 0.0],
                dtype=np.float64,
            )
            try:
                final_orbit = correct_er3bp_periodic(
                    system=sys_target,
                    state_guess=guess,
                    period_f=period_f,
                    is_half_period_residual=is_half_period_residual,
                    free_vars=_FREE_VARS,
                    residual_indices=_RESIDUAL_INDICES,
                    tol=tol,
                )
            except Exception:
                final_orbit = None
            if final_orbit is not None:
                history.append(final_orbit)
            break
        sys_next = ER3BPSystem(
            mu=sys_base.mu,
            e=e_next,
            primary_name=sys_base.primary_name,
            secondary_name=sys_base.secondary_name,
        )
        state_guess = np.array(
            [float(z_next[0]), 0.0, 0.0, 0.0, float(z_next[1]), 0.0],
            dtype=np.float64,
        )
        try:
            orbit = correct_er3bp_periodic(
                system=sys_next,
                state_guess=state_guess,
                period_f=period_f,
                is_half_period_residual=is_half_period_residual,
                free_vars=_FREE_VARS,
                residual_indices=_RESIDUAL_INDICES,
                tol=tol,
            )
        except Exception:
            break

        history.append(orbit)
        z_cur = np.array(
            [float(orbit.state0[IDX_X]), float(orbit.state0[IDX_YDOT]), e_next],
            dtype=np.float64,
        )
        tau_prev = tau

        # Stop once we reach / pass e_target (in the requested direction).
        if (e_dir > 0.0 and e_next >= e_target - 1e-12) or (
            e_dir < 0.0 and e_next <= e_target + 1e-12
        ):
            break

    return history
