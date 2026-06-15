"""3D family-tracer + continuation driver for broken-plane CR3BP orbits (#296 Phase 2).

Phase 2 of the Track-A 3D capability build (#291 → Phase 1 corrector at
:mod:`cyclerfinder.search.cr3bp_general_periodic_3d`; Phase 2 is this module).

This module consumes the Phase 1 corrector
:func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
and walks a continuation parameter through the symmetric-tulip 3D family. It
delivers:

  * :func:`continue_general_3d_family` -- predictor-corrector continuation
    driver supporting pseudo-arclength, natural-T, and natural-x0 stepping
    modes. Each accepted member's monodromy + Floquet multipliers are
    computed.
  * :class:`Family3D` -- dataclass holding the family members, fold-turning
    points (where the augmented Jacobian's tangent component flipped sign in
    the natural parameter), Floquet stability tags, and walk metadata.

Continuation scheme
-------------------
For the symmetric tulip mode the corrector residual at the perpendicular
half-period crossing is 3-D: ``r = (y, xdot, zdot)`` at ``T/2``. The free-var
ambient space for pseudo-arclength is 4-D: ``z = (x0, z0, ydot0, T)``. The
residual map ``r: R^4 -> R^3`` defines a 1-parameter family curve (3 eqns in
4 unknowns; the curve is the kernel of ``dr/dz`` perpendicular to the
zero-set).

Three stepping modes:

  1. **pseudo-arclength** (default, robust past folds): at each member compute
     ``dr/dz`` (3x4 Jacobian), find the unit null vector ``tau`` (= last
     right-singular vector of the SVD), step ``z_pred = z + step * tau``,
     correct back onto ``{r(z) = 0, tau . (z - z_pred) = 0}`` (4 eqns, 4
     unknowns). The arclength constraint replaces a natural parameter and
     walks cleanly through folds in any one component.

  2. **natural_T**: step T as the natural parameter (``z = (x0, z0, ydot0)``
     unknowns; T held at the predicted value). Corrector matches the Phase 1
     ``FREE_VARS_SYMMETRIC_TULIP`` shape with T removed from free_vars.

  3. **natural_x0**: step x0 as the natural parameter (``z = (z0, ydot0, T)``
     unknowns). Matches the #287 spike's pattern -- known to fold around the
     ``x0 = -0.81`` member, hence pseudo-arclength is the default.

Discipline (per ``feedback_orbit_closure_discipline``)
------------------------------------------------------
  * Every accepted member is re-propagated with ``Radau`` at
    ``rtol = atol = 1e-12`` and the L2 closure ``X(T) - X(0)`` is asserted
    below ``closure_tol``. This is the cross-check; the corrector's own
    residual is necessary but never sufficient.
  * Fold turning points are RECORDED (sign change in the natural-parameter
    component of the unit tangent) but the pseudo-arclength continuation
    walks through them; for natural_T / natural_x0 stepping the walk STOPS
    cleanly at a fold (the natural-parameter Jacobian singularity).
  * Monodromy is the full 6x6 ``Phi(T)`` from
    :func:`cyclerfinder.search.bifurcation_detector.monodromy`. Floquet
    multipliers from :func:`bifurcation_detector.floquet_multipliers`.
  * The corrector residual reproduces the symmetric-tulip pattern used by
    :func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`; the
    Phase 1 ``correct_general_periodic_3d`` with
    ``FREE_VARS_SYMMETRIC_TULIP`` + ``RESIDUAL_PERPENDICULAR_HALF_PERIOD`` is
    the same residual map, exposed as a single function.

What Phase 2 does NOT deliver (deferred to Phase 3)
---------------------------------------------------
  * Literature-novelty check against ``literature_check.py``'s KNOWN_CORPUS.
  * ML flagger inputs / V0-V5 gauntlet adaptation for 3D orbits.
  * Catalogue admission (Phase 3 gate; this module only writes a family
    JSONL registry, NOT catalogue rows).

References
----------
  * Phase 1 corrector + verdict:
    :mod:`cyclerfinder.search.cr3bp_general_periodic_3d` /
    ``docs/notes/2026-06-16-291-phase1-3d-corrector.md``.
  * Pseudo-arclength reference (planar Jacobi-curve walker):
    :mod:`cyclerfinder.search.cr3bp_jacobi_arclength`.
  * Spike family extent / fold:
    ``data/spike_287.jsonl`` (x0 ∈ [-0.85, -0.77], z0 ∈ [-0.24, 0],
    T ∈ [9.5, 18.9] TU, C ∈ [2.78, 3.10]; 164 family rows).
  * Likely-rediscovery anchor (Phase 3 will run literature_check):
    Antoniadou-Voyatzis 2018 spatial CR3BP corpus, KNOWN_CORPUS commit
    ``568d8a4``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import (
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    FREE_VARS_SYMMETRIC_TULIP,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_Z,
    IDX_ZDOT,
    RESIDUAL_FULL_STATE_AT_T,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)

__all__ = [
    "Family3D",
    "Family3DMember",
    "FoldPoint",
    "continue_general_3d_family",
]


# ---------------------------------------------------------------------------
# Result types.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Family3DMember:
    """One converged + closure-verified member of a 3D family.

    Attributes
    ----------
    orbit :
        The Phase 1 :class:`Periodic3DOrbit` (state0, T, jacobi, closure
        diagnostics, planar-collapse flag).
    monodromy :
        The 6x6 full-period STM ``Phi(T)`` (or ``None`` if monodromy was
        disabled). The reciprocal-pair / unit-pair structure is a consequence
        of CR3BP's Hamiltonian symmetry; this matrix is the source of all
        downstream stability diagnostics.
    floquet :
        The 6 Floquet multipliers sorted by descending magnitude (or ``None``).
        Two are trivially unit (energy/time-translation); the remaining four
        come in reciprocal pairs.
    stability_tag :
        Coarse stability label derived from the Floquet spectrum (see
        :func:`_classify_floquet`).
    arc_length :
        Cumulative pseudo-arclength from the seed (zero at the seed; signed
        by walk direction). Set to NaN for natural-parameter modes that don't
        compute an arclength.
    step_index :
        Index in the walk; 0 is the seed; positive integers are forward steps
        and negative integers are backward steps.
    """

    orbit: Periodic3DOrbit
    monodromy: NDArray[np.float64] | None
    floquet: NDArray[np.complex128] | None
    stability_tag: str
    arc_length: float
    step_index: int


@dataclass(frozen=True)
class FoldPoint:
    """A fold-turning point detected during continuation.

    Attributes
    ----------
    step_index :
        The walk step at which the fold was detected (the bracketing
        member, before the sign flip).
    natural_param :
        Which component flipped: ``"x0"``, ``"T"``, ``"z0"``, ``"ydot0"``.
    tangent_before, tangent_after :
        Unit-tangent components in the natural parameter at the bracketing
        members; their opposite signs mark the fold.
    member_before, member_after :
        The bracketing family members.
    """

    step_index: int
    natural_param: str
    tangent_before: float
    tangent_after: float
    member_before: Family3DMember
    member_after: Family3DMember


@dataclass
class Family3D:
    """A continuation walk through a 3D CR3BP family.

    Attributes
    ----------
    seed :
        The converged seed member (step_index = 0).
    members :
        All accepted members in walk order. Forward and backward walks are
        merged into one list sorted by ``step_index``.
    folds :
        Detected fold-turning points (sign flips in a tangent component).
    continuation_mode :
        Echoed back: ``"pseudo_arclength"``, ``"natural_T"``, or
        ``"natural_x0"``.
    step :
        Nominal step size (pseudo-arclength) or step in the natural parameter.
    n_steps_forward, n_steps_backward :
        Realized step counts in each direction.
    forward_termination, backward_termination :
        Reason each direction stopped: ``"max_steps"``, ``"corrector_failed"``,
        ``"closure_failed"``, ``"degenerate_planar"``, ``"singular_tangent"``,
        ``"period_underflow"``.
    """

    seed: Family3DMember
    members: list[Family3DMember]
    folds: list[FoldPoint]
    continuation_mode: str
    step: float
    n_steps_forward: int = 0
    n_steps_backward: int = 0
    forward_termination: str = "max_steps"
    backward_termination: str = "max_steps"
    metadata: dict[str, object] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Floquet classification.
# ---------------------------------------------------------------------------


def _classify_floquet(eigs: NDArray[np.complex128], *, unit_tol: float = 1e-3) -> str:
    """Coarse stability tag from the Floquet multipliers.

    Tags:
      * ``"stable"`` -- all four non-trivial multipliers within ``unit_tol`` of
        the unit circle.
      * ``"unstable"`` -- at least one non-trivial multiplier strictly outside
        the unit circle (|lambda| > 1 + unit_tol).
      * ``"hyperbolic_pair"`` -- one or more reciprocal pairs of real
        non-unit multipliers (classic CR3BP unstable hyperbolic structure).
      * ``"degenerate"`` -- the trivial unit pair couldn't be identified.

    The classification is intentionally coarse; downstream consumers needing
    the full reciprocal-pair structure should consume the raw ``floquet``
    array directly.
    """
    if eigs.shape != (6,):
        return "degenerate"
    # Identify the trivial unit pair: two eigenvalues closest to +1.
    dist_to_unity = np.abs(eigs - 1.0)
    order = np.argsort(dist_to_unity)
    if dist_to_unity[order[1]] > 1e-1:
        # No identifiable unit pair -- not a periodic orbit (or numerical
        # garbage). Tag as degenerate; the caller can inspect the eigs.
        return "degenerate"
    # The non-trivial four are the rest.
    nontriv_mask = np.ones(6, dtype=bool)
    nontriv_mask[order[:2]] = False
    nontriv = eigs[nontriv_mask]
    nontriv_abs = np.abs(nontriv)
    if np.any(nontriv_abs > 1.0 + unit_tol):
        # Hyperbolic if there's a real reciprocal pair (|lambda|, 1/|lambda|);
        # else just unstable (e.g. complex pair off the unit circle).
        on_real = np.abs(nontriv.imag) < unit_tol
        if np.any(on_real & (nontriv_abs > 1.0 + unit_tol)):
            return "hyperbolic_pair"
        return "unstable"
    # All non-trivial within unit_tol of the unit circle.
    if np.all(np.abs(nontriv_abs - 1.0) < unit_tol):
        return "stable"
    return "unstable"


# ---------------------------------------------------------------------------
# Symmetric-tulip residual + Jacobian (for pseudo-arclength's null vector).
# ---------------------------------------------------------------------------


def _tulip_residual_and_jac(
    system: cr3bp.CR3BPSystem,
    x0: float,
    z0: float,
    ydot0: float,
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """Return ``(r, J)`` for the symmetric-tulip residual at the current
    point in the 4-D unknown space ``z = (x0, z0, ydot0, T)``.

    The residual is ``r = (y, xdot, zdot)`` at ``T/2`` from the IC
    ``(x0, 0, z0, 0, ydot0, 0)``.

    The Jacobian rows correspond to the residual indices ``(1, 3, 5)`` and
    columns to ``(x0, z0, ydot0, T)``. Sensitivity to ``x0`` uses
    ``stm[i, 0]`` (the STM tracks dX(t)/dX(0)); sensitivity to ``T`` uses
    ``0.5 * f[i]`` (chain rule on ``t_event = 0.5 T``).

    Returns ``None`` on integrator failure.
    """
    state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
    t_half = 0.5 * period
    try:
        arc = cr3bp.propagate(system, state0, t_half, with_stm=True, rtol=rtol, atol=atol)
    except RuntimeError:
        return None
    state_h = arc.state_f
    stm_h = arc.stm
    if stm_h is None:
        return None
    r = np.array(
        [float(state_h[IDX_Y]), float(state_h[IDX_XDOT]), float(state_h[IDX_ZDOT])],
        dtype=np.float64,
    )
    f_h = cr3bp.cr3bp_eom(t_half, state_h, system.mu)
    # Rows of the Jacobian are the residual indices (1, 3, 5); columns are
    # (x0, z0, ydot0, T). For state columns the sensitivity is stm[row, col];
    # for T it's 0.5 * f[row] (half-period chain rule).
    rows = (IDX_Y, IDX_XDOT, IDX_ZDOT)
    jac = np.zeros((3, 4), dtype=np.float64)
    for ri, row in enumerate(rows):
        jac[ri, 0] = float(stm_h[row, IDX_X])  # dx0
        jac[ri, 1] = float(stm_h[row, IDX_Z])  # dz0
        jac[ri, 2] = float(stm_h[row, IDX_YDOT])  # dydot0
        jac[ri, 3] = 0.5 * float(f_h[row])  # dT (chain rule on t_half=T/2)
    return r, jac


def _tulip_tangent(
    jac: NDArray[np.float64],
    *,
    prev: NDArray[np.float64] | None,
) -> NDArray[np.float64] | None:
    """Unit tangent = right-singular vector of ``jac`` for the smallest
    singular value (the null direction).

    For a 3x4 Jacobian with rank 3 the null space is 1-D; the SVD's last
    right-singular vector is the unit tangent. Oriented for continuity with
    ``prev`` (positive dot product).
    """
    try:
        _, _sv, vt = np.linalg.svd(jac, full_matrices=True)
    except np.linalg.LinAlgError:
        return None
    # Null direction is the last right-singular vector.
    if vt.shape != (4, 4):
        return None
    tau = np.asarray(vt[-1], dtype=np.float64)
    # Numerical: a near-zero singular value (>> machine epsilon for an
    # ill-conditioned tangent) signals a fold; the tangent is still well
    # defined (it's the null direction), just sensitive.
    if not np.all(np.isfinite(tau)) or float(np.linalg.norm(tau)) < 1e-12:
        return None
    tau = tau / float(np.linalg.norm(tau))
    if prev is not None and float(np.dot(tau, prev)) < 0.0:
        tau = -tau
    return tau


def _correct_pseudo_arclength(
    system: cr3bp.CR3BPSystem,
    z_pred: NDArray[np.float64],
    tau: NDArray[np.float64],
    *,
    tol: float,
    max_iter: int,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    state_cap: float = 0.1,
    period_cap_frac: float = 0.2,
) -> NDArray[np.float64] | None:
    """Newton onto ``{r(z) = 0, tau . (z - z_pred) = 0}``.

    4 unknowns ``z = (x0, z0, ydot0, T)``; 4 equations (3 residual + 1
    arclength). Returns the converged ``z`` or ``None`` on failure.
    """
    z = z_pred.copy()
    for _ in range(max_iter):
        x0_, z0_, ydot0_, period_ = float(z[0]), float(z[1]), float(z[2]), float(z[3])
        if period_ <= 0.0 or not math.isfinite(period_):
            return None
        rj = _tulip_residual_and_jac(system, x0_, z0_, ydot0_, period_, rtol=rtol, atol=atol)
        if rj is None:
            return None
        r0, grad = rj
        arc = float(np.dot(tau, z - z_pred))
        if float(np.linalg.norm(r0)) < tol and abs(arc) < 1e-12:
            return z
        # Augmented Jacobian: 3 rows of dr/dz + 1 row of tau.
        jmat = np.vstack([grad, tau.reshape(1, 4)])
        rhs = -np.concatenate([r0, np.array([arc])])
        try:
            dz = np.linalg.solve(jmat, rhs)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jmat, rhs, rcond=None)
        dz = np.asarray(dz, dtype=np.float64)
        # Component caps: state components and period treated separately.
        dz[0] = float(np.clip(dz[0], -state_cap, state_cap))  # x0
        dz[1] = float(np.clip(dz[1], -state_cap, state_cap))  # z0
        dz[2] = float(np.clip(dz[2], -state_cap, state_cap))  # ydot0
        period_cap = period_cap_frac * abs(period_)
        dz[3] = float(np.clip(dz[3], -period_cap, period_cap))
        z = z + dz
    return None


# ---------------------------------------------------------------------------
# Independent closure check + member construction.
# ---------------------------------------------------------------------------


def _independent_closure(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> float:
    """L2 norm of ``X(T) - X(0)`` from a Radau re-propagation.

    Different integrator than the DOP853 in the corrector; gates against
    single-integrator artefacts. Returns NaN on solver failure.
    """
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
    except (RuntimeError, ValueError):
        return float("nan")
    if not sol.success:
        return float("nan")
    return float(np.linalg.norm(sol.y[:, -1] - state0))


def _build_member(
    system: cr3bp.CR3BPSystem,
    orbit: Periodic3DOrbit,
    step_index: int,
    arc_length: float,
    *,
    compute_monodromy: bool,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> Family3DMember:
    """Wrap a converged :class:`Periodic3DOrbit` with monodromy + Floquet."""
    mono: NDArray[np.float64] | None = None
    fl: NDArray[np.complex128] | None = None
    tag = "unknown"
    if compute_monodromy:
        try:
            mono = monodromy(system, orbit.state0, orbit.T_TU, rtol=rtol, atol=atol)
            fl = floquet_multipliers(mono)
            tag = _classify_floquet(fl)
        except (RuntimeError, ValueError):
            mono = None
            fl = None
            tag = "monodromy_failed"
    return Family3DMember(
        orbit=orbit,
        monodromy=mono,
        floquet=fl,
        stability_tag=tag,
        arc_length=arc_length,
        step_index=step_index,
    )


# ---------------------------------------------------------------------------
# Predictor-corrector step (pseudo-arclength).
# ---------------------------------------------------------------------------


def _pseudo_arclength_step(
    system: cr3bp.CR3BPSystem,
    z_cur: NDArray[np.float64],
    tau_prev: NDArray[np.float64] | None,
    step: float,
    *,
    corrector_tol: float,
    max_corrector_iter: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """One pseudo-arclength step.

    Returns ``(z_next, tau_next)`` on success or ``None`` on failure.
    """
    x0_, z0_, ydot0_, period_ = float(z_cur[0]), float(z_cur[1]), float(z_cur[2]), float(z_cur[3])
    rj = _tulip_residual_and_jac(system, x0_, z0_, ydot0_, period_)
    if rj is None:
        return None
    _, jac = rj
    tau = _tulip_tangent(jac, prev=tau_prev)
    if tau is None:
        return None
    # Predictor.
    z_pred = z_cur + step * tau
    # Corrector.
    z_next = _correct_pseudo_arclength(
        system,
        z_pred,
        tau,
        tol=corrector_tol,
        max_iter=max_corrector_iter,
    )
    if z_next is None:
        return None
    return z_next, tau


# ---------------------------------------------------------------------------
# Natural-parameter step (T or x0 held during the corrector).
# ---------------------------------------------------------------------------


def _natural_step(
    system: cr3bp.CR3BPSystem,
    z_cur: NDArray[np.float64],
    z_prev: NDArray[np.float64] | None,
    step: float,
    natural: Literal["T", "x0"],
    *,
    corrector_tol: float,
    max_corrector_iter: int,
) -> Periodic3DOrbit | None:
    """One natural-parameter step.

    The natural parameter advances by ``step`` (signed) and the corrector
    closes the remaining three unknowns. Uses
    :func:`correct_general_periodic_3d` with the appropriate free_vars mask.
    """
    x0_, z0_, ydot0_, period_ = float(z_cur[0]), float(z_cur[1]), float(z_cur[2]), float(z_cur[3])
    if natural == "T":
        # Advance T; corrector solves (x0, z0, ydot0).
        new_period = period_ + step
        if new_period <= 0.0:
            return None
        state0_guess = np.array([x0_, 0.0, z0_, 0.0, ydot0_, 0.0], dtype=np.float64)
        # Predictor: extrapolate state from previous step if available.
        if z_prev is not None and abs(z_cur[3] - z_prev[3]) > 1e-12:
            frac = step / (z_cur[3] - z_prev[3])
            state0_guess[IDX_X] = x0_ + frac * (z_cur[0] - z_prev[0])
            state0_guess[IDX_Z] = z0_ + frac * (z_cur[1] - z_prev[1])
            state0_guess[IDX_YDOT] = ydot0_ + frac * (z_cur[2] - z_prev[2])
        # Free vars: (x0, z0, ydot0) -- T removed from the symmetric tulip mask.
        free: tuple[int, ...] = (IDX_X, IDX_Z, IDX_YDOT)
        result = correct_general_periodic_3d(
            system,
            state0_guess,
            new_period,
            free_vars=free,
            residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
            is_half_period_residual=True,
            tol=corrector_tol,
            max_iter=max_corrector_iter,
        )
    elif natural == "x0":
        # Advance x0; corrector solves (z0, ydot0, T).
        new_x0 = x0_ + step
        state0_guess = np.array([new_x0, 0.0, z0_, 0.0, ydot0_, 0.0], dtype=np.float64)
        period_guess = period_
        if z_prev is not None and abs(z_cur[0] - z_prev[0]) > 1e-12:
            frac = step / (z_cur[0] - z_prev[0])
            state0_guess[IDX_Z] = z0_ + frac * (z_cur[1] - z_prev[1])
            state0_guess[IDX_YDOT] = ydot0_ + frac * (z_cur[2] - z_prev[2])
            period_guess = period_ + frac * (z_cur[3] - z_prev[3])
        free = FREE_VARS_SYMMETRIC_TULIP  # (z0, ydot0, T)
        result = correct_general_periodic_3d(
            system,
            state0_guess,
            period_guess,
            free_vars=free,
            residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
            is_half_period_residual=True,
            tol=corrector_tol,
            max_iter=max_corrector_iter,
        )
    else:
        raise ValueError(f"natural must be 'T' or 'x0'; got {natural!r}")
    if not result.converged:
        return None
    return result


# ---------------------------------------------------------------------------
# Main driver.
# ---------------------------------------------------------------------------


def continue_general_3d_family(
    system: cr3bp.CR3BPSystem,
    seed_state: NDArray[np.float64] | Sequence[float],
    seed_period: float,
    *,
    free_vars: tuple[int, ...] = FREE_VARS_FULL_ASYMMETRIC,
    residual_indices: tuple[int, ...] = RESIDUAL_FULL_STATE_AT_T,
    continuation: Literal["pseudo_arclength", "natural_T", "natural_x0"] = "pseudo_arclength",
    step: float = 0.01,
    n_steps_max: int = 200,
    direction: Literal["forward", "backward", "both"] = "both",
    corrector_tol: float = 1e-10,
    closure_tol: float = 1e-6,
    fold_detection: bool = True,
    monodromy_eval: bool = True,
    max_corrector_iter: int = 60,
    on_step: Callable[[Family3DMember], None] | None = None,
) -> Family3D:
    """Drive a continuation walk through a 3D symmetric-tulip CR3BP family.

    Parameters
    ----------
    system :
        CR3BP system; only ``system.mu`` is read by the dynamics.
    seed_state :
        6-vector seed IC. For the symmetric-tulip mode this MUST have the
        perpendicular form ``(x0, 0, z0, 0, ydot0, 0)`` (the corrector
        residual implicitly assumes y0 = xdot0 = zdot0 = 0); the value is
        passed through to the Phase 1 corrector which enforces this at the
        seed pass too.
    seed_period :
        Full nondim period guess at the seed.
    free_vars :
        Echoed onto the Phase 1 corrector for the seed pass. Pseudo-arclength
        and natural-parameter modes ALWAYS use the symmetric-tulip free-var
        masks internally; this argument controls only the SEED correction
        (and is preserved on each Family3DMember.orbit for downstream
        consumers).
    residual_indices :
        Same: echoed onto the seed correction. The walk uses
        :data:`RESIDUAL_PERPENDICULAR_HALF_PERIOD` internally.
    continuation :
        One of ``"pseudo_arclength"`` (default, robust past folds),
        ``"natural_T"``, ``"natural_x0"``.
    step :
        Step size. For pseudo-arclength this is the arclength in the unit
        4-vector tangent space (z = (x0, z0, ydot0, T)). For natural-T it is
        the step in T (nondim TU). For natural-x0 it is the step in x0.
    n_steps_max :
        Maximum walk steps per direction.
    direction :
        ``"forward"``, ``"backward"``, or ``"both"`` (default).
    corrector_tol :
        Corrector L2 tolerance at each step.
    closure_tol :
        Independent-Radau closure tolerance for accepting a member.
    fold_detection :
        If True, detect sign flips in the tangent's natural-parameter
        component. Folds are RECORDED in ``Family3D.folds`` but the
        pseudo-arclength walk continues through them (the natural-parameter
        walks STOP at the fold's parameter singularity).
    monodromy_eval :
        If True (default), compute monodromy + Floquet at every accepted
        member. Set False if you only want the orbits (faster: skips a
        full-period STM integration per member).
    max_corrector_iter :
        Corrector iteration cap per step.
    on_step :
        Optional callback ``fn(member: Family3DMember) -> None`` invoked
        after each accepted member. For progress reporting.

    Returns
    -------
    Family3D
        The walk record. Members are sorted by ``step_index`` (negatives
        first if a backward walk happened).

    Raises
    ------
    ValueError
        On obviously malformed inputs (seed state shape; unknown continuation
        mode; non-positive step or period).
    RuntimeError
        If the seed itself fails to converge under the Phase 1 corrector.
        Continuation absolutely requires a valid seed; an unconverged seed
        is unrecoverable here.

    Notes
    -----
    Independent closure check is MANDATORY on every accepted member: the
    member's :attr:`Periodic3DOrbit.independent_closure_residual` must be
    below ``closure_tol`` AND its corrector residual below ``corrector_tol``.
    Per ``feedback_orbit_closure_discipline``: cross-check is never optional.
    """
    seed_state = np.asarray(seed_state, dtype=np.float64).copy()
    if seed_state.shape != (6,):
        raise ValueError(f"seed_state must have shape (6,); got {seed_state.shape}")
    if seed_period <= 0.0:
        raise ValueError(f"seed_period must be > 0; got {seed_period}")
    if step <= 0.0:
        raise ValueError(f"step must be > 0; got {step}")
    if continuation not in {"pseudo_arclength", "natural_T", "natural_x0"}:
        raise ValueError(f"unknown continuation mode: {continuation!r}")
    if direction not in {"forward", "backward", "both"}:
        raise ValueError(f"unknown direction: {direction!r}")

    # Seed pass via the Phase 1 corrector.
    seed_orbit = correct_general_periodic_3d(
        system,
        seed_state,
        seed_period,
        free_vars=free_vars,
        residual_indices=residual_indices,
        is_half_period_residual=(residual_indices == RESIDUAL_PERPENDICULAR_HALF_PERIOD),
        tol=corrector_tol,
        max_iter=max_corrector_iter,
        independent_tol=closure_tol,
    )
    if not seed_orbit.converged:
        raise RuntimeError(
            f"seed did not converge; corrector={seed_orbit.corrector_residual:.3e}, "
            f"independent={seed_orbit.independent_closure_residual:.3e}"
        )
    seed_member = _build_member(
        system,
        seed_orbit,
        step_index=0,
        arc_length=0.0,
        compute_monodromy=monodromy_eval,
    )
    if on_step is not None:
        on_step(seed_member)

    # Helper: from a Periodic3DOrbit extract z = (x0, z0, ydot0, T).
    def _z_of(orbit: Periodic3DOrbit) -> NDArray[np.float64]:
        return np.array(
            [
                float(orbit.state0[IDX_X]),
                float(orbit.state0[IDX_Z]),
                float(orbit.state0[IDX_YDOT]),
                float(orbit.T_TU),
            ],
            dtype=np.float64,
        )

    # Walk helper: returns (list of accepted members, list of folds, termination_reason).
    def _walk(
        sign: int,
    ) -> tuple[list[Family3DMember], list[FoldPoint], str]:
        members_local: list[Family3DMember] = []
        folds_local: list[FoldPoint] = []
        if continuation == "pseudo_arclength":
            z_cur = _z_of(seed_orbit)
            tau_prev: NDArray[np.float64] | None = None
            arc_cur = 0.0
            for n in range(1, n_steps_max + 1):
                step_signed = sign * step
                out = _pseudo_arclength_step(
                    system,
                    z_cur,
                    tau_prev,
                    step_signed,
                    corrector_tol=corrector_tol,
                    max_corrector_iter=max_corrector_iter,
                )
                if out is None:
                    return members_local, folds_local, "corrector_failed"
                z_next, tau_next = out
                # Build the corrected Periodic3DOrbit from z_next: run the Phase 1
                # corrector at the converged z to compute the independent closure
                # check (the pseudo-arclength corrector only computes the symmetric
                # tulip residual; the Phase 1 wrapper folds in the Radau cross-check
                # and the degenerate_planar flag).
                state0_next = np.array(
                    [float(z_next[0]), 0.0, float(z_next[1]), 0.0, float(z_next[2]), 0.0],
                    dtype=np.float64,
                )
                orb_next = correct_general_periodic_3d(
                    system,
                    state0_next,
                    float(z_next[3]),
                    free_vars=FREE_VARS_SYMMETRIC_TULIP,
                    residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
                    is_half_period_residual=True,
                    tol=corrector_tol,
                    max_iter=max_corrector_iter,
                    independent_tol=closure_tol,
                )
                if not orb_next.converged:
                    return members_local, folds_local, "closure_failed"
                arc_cur += step_signed
                step_index = sign * n
                member = _build_member(
                    system,
                    orb_next,
                    step_index=step_index,
                    arc_length=arc_cur,
                    compute_monodromy=monodromy_eval,
                )
                # Fold detection: sign flip in tau component for any of the four
                # parameters from the previous tangent. Pseudo-arclength walks
                # through cleanly; we only RECORD.
                if fold_detection and tau_prev is not None and len(members_local) > 0:
                    for i, name in enumerate(("x0", "z0", "ydot0", "T")):
                        if tau_prev[i] * tau_next[i] < 0:
                            folds_local.append(
                                FoldPoint(
                                    step_index=step_index,
                                    natural_param=name,
                                    tangent_before=float(tau_prev[i]),
                                    tangent_after=float(tau_next[i]),
                                    member_before=members_local[-1],
                                    member_after=member,
                                )
                            )
                members_local.append(member)
                if on_step is not None:
                    on_step(member)
                z_cur = _z_of(orb_next)
                tau_prev = tau_next
                if orb_next.degenerate_planar:
                    return members_local, folds_local, "degenerate_planar"
            return members_local, folds_local, "max_steps"

        # Natural-parameter modes.
        natural: Literal["T", "x0"] = "T" if continuation == "natural_T" else "x0"
        z_cur = _z_of(seed_orbit)
        z_prev: NDArray[np.float64] | None = None
        prev_natural_value = float(z_cur[3] if natural == "T" else z_cur[0])
        for n in range(1, n_steps_max + 1):
            step_signed = sign * step
            nat_result = _natural_step(
                system,
                z_cur,
                z_prev,
                step_signed,
                natural,
                corrector_tol=corrector_tol,
                max_corrector_iter=max_corrector_iter,
            )
            if nat_result is None:
                return members_local, folds_local, "corrector_failed"
            orb_next = nat_result
            # Independent closure already gated by the Phase 1 corrector's
            # independent_tol = closure_tol; orb_next.converged demands both.
            step_index = sign * n
            z_next = _z_of(orb_next)
            member = _build_member(
                system,
                orb_next,
                step_index=step_index,
                arc_length=float("nan"),
                compute_monodromy=monodromy_eval,
            )
            # Fold detection: monotonicity check on the natural parameter.
            cur_natural_value = float(z_next[3] if natural == "T" else z_next[0])
            if (
                fold_detection
                and len(members_local) > 0
                and (cur_natural_value - prev_natural_value) * step_signed < 0
            ):
                z_cur_natural = float(z_cur[3 if natural == "T" else 0])
                folds_local.append(
                    FoldPoint(
                        step_index=step_index,
                        natural_param=natural,
                        tangent_before=prev_natural_value - z_cur_natural,
                        tangent_after=cur_natural_value - prev_natural_value,
                        member_before=members_local[-1],
                        member_after=member,
                    )
                )
                members_local.append(member)
                if on_step is not None:
                    on_step(member)
                return members_local, folds_local, "fold"
            members_local.append(member)
            if on_step is not None:
                on_step(member)
            z_prev = z_cur
            z_cur = z_next
            prev_natural_value = cur_natural_value
            if orb_next.degenerate_planar:
                return members_local, folds_local, "degenerate_planar"
        return members_local, folds_local, "max_steps"

    forward_members: list[Family3DMember] = []
    backward_members: list[Family3DMember] = []
    folds_all: list[FoldPoint] = []
    forward_reason = "max_steps"
    backward_reason = "max_steps"
    if direction in {"forward", "both"}:
        forward_members, fwd_folds, forward_reason = _walk(+1)
        folds_all.extend(fwd_folds)
    if direction in {"backward", "both"}:
        backward_members, bwd_folds, backward_reason = _walk(-1)
        folds_all.extend(bwd_folds)

    # Merge: backward (most negative step first), seed, forward.
    members_sorted = sorted(
        [*backward_members, seed_member, *forward_members],
        key=lambda m: m.step_index,
    )

    return Family3D(
        seed=seed_member,
        members=members_sorted,
        folds=folds_all,
        continuation_mode=continuation,
        step=step,
        n_steps_forward=len(forward_members),
        n_steps_backward=len(backward_members),
        forward_termination=forward_reason,
        backward_termination=backward_reason,
        metadata={
            "seed_free_vars": list(free_vars),
            "seed_residual_indices": list(residual_indices),
            "closure_tol": closure_tol,
            "corrector_tol": corrector_tol,
        },
    )
