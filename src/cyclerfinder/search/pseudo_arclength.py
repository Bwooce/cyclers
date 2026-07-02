"""Generic pseudo-arclength continuation primitive (#524).

This project has pseudo-arclength continuation in at least four SPECIALIZED
forms already (``search/cr3bp_jacobi_arclength.py`` for CR3BP cyclers in
``(x0, C)``, ``search/mu_continuation.py`` for ``(x0, C, mu)``,
``search/qp_tori_arclength.py``, ``search/narc_continuation.py``) -- each
reimplements the SAME predict-correct-tangent machinery for its own
hardcoded residual function. #524 was originally scoped as "not yet built"
before this was checked; it isn't a from-scratch build, it's extracting the
shared machinery into a reusable primitive so a FUTURE search (#530's
Sun-Earth co-orbital exploration, #531's manifold-branch search, or any
one-off ``scripts/run_*.py``) can call ``continue_curve()`` on its own
residual function instead of hand-rolling a predictor-corrector loop --
directly the failure mode that produced bugs in #523/#530/#531 this session
(a wrong propagation horizon, a t_eval sort-direction crash, non-converging
candidates costing more than converging ones) and the #496/#520 lesson that
FIXED-GRID sampling misses narrow solution curves a continuation walk would
follow directly.

Scope, deliberately: co-dimension-1 curves only (``M = N-1`` scalar
constraint equations in an ``N``-vector unknown ``z``), the well-posed
textbook case with a UNIQUE (up to sign) tangent direction -- the same case
``cr3bp_jacobi_arclength.py`` handles (2 unknowns, 1 equation). This does
NOT attempt to generalize ``mu_continuation.py``'s ambiguous-tangent
2-surface case (1 equation, 3 unknowns, tangent selected by continuity with
the previous step rather than uniquely determined); that is out of scope
here.

Pure: numpy/scipy only. No CR3BP-specific code -- the caller supplies
``residual_fn`` and (optionally) ``jacobian_fn``.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import NDArray

ResidualFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]
JacobianFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]


class ContinuationStopReason(StrEnum):
    """Why a continuation walk stopped."""

    MAX_STEPS = "max_steps"
    CORRECTOR_FAILED = "corrector_failed"
    RESIDUAL_UNAVAILABLE = "residual_unavailable"
    TARGET_REACHED = "target_reached"


@dataclass(frozen=True)
class ContinuationPoint:
    """One converged point on the curve."""

    z: NDArray[np.float64]
    tangent: NDArray[np.float64]
    residual_norm: float
    step_index: int
    arclength_s: float


@dataclass(frozen=True)
class ContinuationCurve:
    """The full walked curve plus how it ended."""

    points: list[ContinuationPoint]
    stop_reason: ContinuationStopReason
    notes: str = ""

    def z_values(self) -> NDArray[np.float64]:
        """Stack every point's ``z`` into an ``(n_points, N)`` array."""
        return np.array([p.z for p in self.points], dtype=np.float64)


def _numerical_jacobian(
    residual_fn: ResidualFn, z: NDArray[np.float64], r0: NDArray[np.float64], *, h: float = 1e-7
) -> NDArray[np.float64] | None:
    """Forward-difference Jacobian ``dr/dz`` (M x N) when none is supplied."""
    n = z.shape[0]
    m = r0.shape[0]
    jac = np.zeros((m, n), dtype=np.float64)
    for i in range(n):
        zp = z.copy()
        zp[i] += h
        rp = residual_fn(zp)
        if rp is None:
            return None
        jac[:, i] = (rp - r0) / h
    return jac


def compute_tangent(
    residual_fn: ResidualFn,
    z: NDArray[np.float64],
    *,
    jacobian_fn: JacobianFn | None = None,
    prev_tangent: NDArray[np.float64] | None = None,
) -> NDArray[np.float64] | None:
    """Unit tangent to the curve ``residual_fn(z) = 0`` at ``z``.

    ``residual_fn`` must return an ``(N-1,)`` vector for an ``(N,)`` z (the
    co-dimension-1 case this module targets) -- the null space of the
    ``(N-1, N)`` Jacobian is then exactly 1-D and the tangent is unique up to
    sign, oriented for continuity with ``prev_tangent`` (or, at the very
    first point with no ``prev_tangent``, toward increasing ``z[-1]`` by
    convention -- callers that care about a specific orientation should pass
    ``prev_tangent`` explicitly).

    Returns ``None`` if the residual/Jacobian is unavailable at ``z``.
    """
    r0 = residual_fn(z)
    if r0 is None:
        return None
    jac = jacobian_fn(z) if jacobian_fn is not None else _numerical_jacobian(residual_fn, z, r0)
    if jac is None:
        return None
    n = z.shape[0]
    if jac.shape != (n - 1, n):
        raise ValueError(
            f"compute_tangent: residual must have N-1 components for an N-vector z "
            f"(co-dimension 1); got jacobian shape {jac.shape} for N={n}."
        )
    _, _, vt = np.linalg.svd(jac)
    tan = np.asarray(vt[-1], dtype=np.float64)
    if prev_tangent is None:
        if tan[-1] < 0:
            tan = -tan
    elif float(tan @ prev_tangent) < 0:
        tan = -tan
    return tan


def _correct(
    residual_fn: ResidualFn,
    z_pred: NDArray[np.float64],
    tangent_pred: NDArray[np.float64],
    *,
    jacobian_fn: JacobianFn | None,
    tol: float,
    max_iter: int,
    step_caps: NDArray[np.float64] | None,
) -> tuple[NDArray[np.float64], float] | None:
    """Newton onto ``{residual(z)=0, tangent_pred . (z - z_pred) = 0}``.

    Square (N x N) system: N-1 residual equations + 1 arclength equation.
    Solved directly (not least-squares) since it is exactly determined in the
    co-dimension-1 case this module targets.
    """
    z = z_pred.copy()
    n = z.shape[0]
    for _ in range(max_iter):
        r0 = residual_fn(z)
        if r0 is None:
            return None
        arc = float(tangent_pred @ (z - z_pred))
        res_norm = float(np.hypot(np.linalg.norm(r0), arc))
        if res_norm < tol:
            return z, res_norm
        jac = jacobian_fn(z) if jacobian_fn is not None else _numerical_jacobian(residual_fn, z, r0)
        if jac is None:
            return None
        full_jac = np.vstack([jac, tangent_pred.reshape(1, n)])
        full_res = np.concatenate([r0, [arc]])
        try:
            dz = np.linalg.solve(full_jac, -full_res)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(full_jac, -full_res, rcond=None)
        if step_caps is not None:
            dz = np.clip(dz, -step_caps, step_caps)
        z = z + dz
    return None


def continue_curve(
    residual_fn: ResidualFn,
    z0: NDArray[np.float64],
    *,
    jacobian_fn: JacobianFn | None = None,
    step_size: float = 1e-2,
    max_steps: int = 200,
    tol: float = 1e-10,
    max_iter: int = 40,
    step_caps: Sequence[float] | None = None,
    target_index: int | None = None,
    target_value: float | None = None,
    initial_tangent: NDArray[np.float64] | None = None,
) -> ContinuationCurve:
    """Walk the co-dimension-1 curve ``residual_fn(z) = 0`` from ``z0``.

    Predictor-corrector pseudo-arclength continuation (Keller 1977): predict
    ``z + step_size * tangent``, correct back onto the curve with the
    arclength constraint. Turns folds (where a single coordinate's
    natural-parameter derivative diverges) because the constraint is on
    ARCLENGTH, not on any one coordinate -- the #496/#520 lesson that a
    FIXED-GRID sample in one coordinate can step clean over a narrow
    solution region does not apply here; the walk follows the curve however
    it bends.

    ``z0`` must already (approximately) satisfy ``residual_fn(z0) ~ 0`` --
    this function does not do the initial correction, only the walk. If
    ``target_index``/``target_value`` are given, the walk stops (with
    ``TARGET_REACHED``) once ``z[target_index]`` crosses ``target_value``
    between two consecutive corrected points (no attempt to land exactly on
    it -- pass the crossing bracket to a separate landing/projection solve if
    an exact target is needed, matching ``land_at_jacobi``'s pattern in
    ``cr3bp_jacobi_arclength.py``).
    """
    z = np.asarray(z0, dtype=np.float64)
    caps = np.asarray(step_caps, dtype=np.float64) if step_caps is not None else None

    tangent = initial_tangent
    if tangent is None:
        tangent = compute_tangent(residual_fn, z, jacobian_fn=jacobian_fn)
        if tangent is None:
            return ContinuationCurve(
                points=[], stop_reason=ContinuationStopReason.RESIDUAL_UNAVAILABLE
            )

    r0 = residual_fn(z)
    if r0 is None:
        return ContinuationCurve(points=[], stop_reason=ContinuationStopReason.RESIDUAL_UNAVAILABLE)
    points = [
        ContinuationPoint(
            z=z.copy(),
            tangent=tangent.copy(),
            residual_norm=float(np.linalg.norm(r0)),
            step_index=0,
            arclength_s=0.0,
        )
    ]

    s = 0.0
    for step in range(1, max_steps + 1):
        z_pred = z + step_size * tangent
        corrected = _correct(
            residual_fn,
            z_pred,
            tangent,
            jacobian_fn=jacobian_fn,
            tol=tol,
            max_iter=max_iter,
            step_caps=caps,
        )
        if corrected is None:
            return ContinuationCurve(
                points=points, stop_reason=ContinuationStopReason.CORRECTOR_FAILED
            )
        z_new, res_norm = corrected
        new_tangent = compute_tangent(
            residual_fn, z_new, jacobian_fn=jacobian_fn, prev_tangent=tangent
        )
        if new_tangent is None:
            return ContinuationCurve(
                points=points, stop_reason=ContinuationStopReason.RESIDUAL_UNAVAILABLE
            )
        s += step_size

        if target_index is not None and target_value is not None:
            prev_val = z[target_index]
            new_val = z_new[target_index]
            crossed = (prev_val - target_value) * (new_val - target_value) <= 0.0
        else:
            crossed = False

        z, tangent = z_new, new_tangent
        points.append(
            ContinuationPoint(
                z=z.copy(),
                tangent=tangent.copy(),
                residual_norm=res_norm,
                step_index=step,
                arclength_s=s,
            )
        )
        if crossed:
            return ContinuationCurve(
                points=points, stop_reason=ContinuationStopReason.TARGET_REACHED
            )

    return ContinuationCurve(points=points, stop_reason=ContinuationStopReason.MAX_STEPS)


__all__ = [
    "ContinuationCurve",
    "ContinuationPoint",
    "ContinuationStopReason",
    "JacobianFn",
    "ResidualFn",
    "compute_tangent",
    "continue_curve",
]
