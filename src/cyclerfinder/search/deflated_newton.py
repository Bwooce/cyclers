"""Generic deflated-Newton root-enumeration primitive (#524).

Farrell, Birkisson & Funke (2015, SIAM J. Sci. Comput. 37(4):A2026-A2045,
"Deflation techniques for finding distinct solutions of nonlinear partial
differential equations"): given an already-found root ``u_i`` of
``F(u) = 0``, the DEFLATED system

    G(u) = M(u; {u_i}) F(u) = 0,     M(u; {u_i}) = prod_i (1/||u-u_i||^p + shift)

has the SAME roots as ``F`` away from any ``u_i`` (``M`` is finite and
nonzero there) but ``M(u) -> infinity`` as ``u -> u_i``, which repels a
Newton iteration away from already-found roots without moving OTHER roots.
Running Newton on ``G`` from a range of seeds, deflating each new root as it
is found, enumerates DISTINCT roots of ``F`` one at a time instead of the
same root being re-found from every seed in its basin.

Confirmed absent from this codebase before this module (searched: the two
prior "deflat" string hits are ``core.lambert``'s deflation-ANGLE handling
and ``search.bifurcation_detector._deflated_determinant``, a Doedel-1991
bifurcation TEST FUNCTION for detecting folds along a KNOWN continuation
branch -- neither is basin-repulsion root enumeration). This module is a
distinct, genuinely new capability, unlike ``pseudo_arclength.py``
(#524's other half), which turned out to already exist in specialized form.

Motivating case: ``cr3bp_periodic.correct_symmetric_fixed_jacobi`` solves a
SCALAR residual in ``x0`` at fixed ``(C, half_crossings)``; #504's
Pluto-Charon (3,1) sweep independently found a stable-but-wrong-topology
root by hand-picked seed near the (3,2) family's parameter range -- exactly
the kind of co-existing root deflation is built to enumerate systematically.

Pure: numpy only, no CR3BP-specific code -- callers supply
``residual_fn``/``jacobian_fn``, same contract as ``pseudo_arclength.py``.
Supports both the scalar case (``N=1``, e.g. a fixed-Jacobi symmetric-orbit
corrector) and general square ``N``-dimensional systems.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import NDArray

ResidualFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]
JacobianFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]


class NewtonStopReason(StrEnum):
    """Why a single deflated-Newton solve stopped."""

    CONVERGED = "converged"
    MAX_ITER = "max_iter"
    RESIDUAL_UNAVAILABLE = "residual_unavailable"
    SINGULAR = "singular"


@dataclass(frozen=True)
class DeflatedNewtonResult:
    """Outcome of one ``newton_deflated`` solve."""

    z: NDArray[np.float64] | None
    stop_reason: NewtonStopReason
    residual_norm: float
    n_iter: int


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
        jac[:, i] = (np.atleast_1d(np.asarray(rp, dtype=np.float64)) - r0) / h
    return jac


def deflation_factor(
    z: NDArray[np.float64],
    known_roots: Sequence[NDArray[np.float64]],
    *,
    p: float = 2.0,
    shift: float = 1.0,
) -> tuple[float, NDArray[np.float64]]:
    """Farrell deflation factor ``M(z)`` and its gradient at ``z``.

    ``M(z) = prod_i (1/||z - u_i||^p + shift)``. Returns ``(inf, [inf,...])``
    if ``z`` coincides with a known root (caller must treat that as a
    singular/repelled state, never divide by the returned gradient).

    The gradient uses the log-derivative product rule
    (``grad(prod f_i) = (prod f_i) * sum(grad(f_i) / f_i)``), valid because
    every factor ``1/||z-u_i||^p + shift`` is strictly positive except
    exactly at ``u_i``.
    """
    if not known_roots:
        return 1.0, np.zeros_like(z)
    m_total = 1.0
    grad_log_sum = np.zeros_like(z)
    for u in known_roots:
        diff = z - np.asarray(u, dtype=np.float64)
        dist = float(np.linalg.norm(diff))
        if dist < 1e-12:
            return math.inf, np.full_like(z, math.inf)
        m_i = 1.0 / dist**p + shift
        grad_m_i = -p * dist ** (-p - 2) * diff
        m_total *= m_i
        grad_log_sum = grad_log_sum + grad_m_i / m_i
    return m_total, m_total * grad_log_sum


def newton_deflated(
    residual_fn: ResidualFn,
    z0: NDArray[np.float64],
    known_roots: Sequence[NDArray[np.float64]] = (),
    *,
    jacobian_fn: JacobianFn | None = None,
    tol: float = 1e-10,
    max_iter: int = 50,
    p: float = 2.0,
    shift: float = 1.0,
    step_cap: NDArray[np.float64] | float | None = None,
) -> DeflatedNewtonResult:
    """Newton-solve the deflated system ``M(z; known_roots) * residual_fn(z) = 0``.

    With ``known_roots=()`` this is plain Newton on ``residual_fn``. Each
    additional known root repels the iterate away from it (without shifting
    where any OTHER root of ``residual_fn`` sits), so calling this
    repeatedly with a growing ``known_roots`` list from the same or nearby
    seeds enumerates distinct roots (see :func:`enumerate_roots`).

    Square systems use a direct solve; non-square (over/under-determined)
    Jacobians fall back to least-squares, matching
    ``pseudo_arclength.py``'s corrector.
    """
    z = np.atleast_1d(np.asarray(z0, dtype=np.float64)).copy()
    r_last = np.full(1, np.inf)
    for it in range(1, max_iter + 1):
        r_raw = residual_fn(z)
        if r_raw is None:
            return DeflatedNewtonResult(
                z=None,
                stop_reason=NewtonStopReason.RESIDUAL_UNAVAILABLE,
                residual_norm=math.inf,
                n_iter=it,
            )
        r = np.atleast_1d(np.asarray(r_raw, dtype=np.float64))
        r_last = r
        raw_norm = float(np.linalg.norm(r))
        if raw_norm < tol:
            return DeflatedNewtonResult(
                z=z.copy(),
                stop_reason=NewtonStopReason.CONVERGED,
                residual_norm=raw_norm,
                n_iter=it,
            )

        m, grad_m = deflation_factor(z, known_roots, p=p, shift=shift)
        if not math.isfinite(m):
            return DeflatedNewtonResult(
                z=None, stop_reason=NewtonStopReason.SINGULAR, residual_norm=raw_norm, n_iter=it
            )

        jac = jacobian_fn(z) if jacobian_fn is not None else _numerical_jacobian(residual_fn, z, r)
        if jac is None:
            return DeflatedNewtonResult(
                z=None,
                stop_reason=NewtonStopReason.RESIDUAL_UNAVAILABLE,
                residual_norm=raw_norm,
                n_iter=it,
            )
        jac = np.atleast_2d(np.asarray(jac, dtype=np.float64))

        jac_g = m * jac + np.outer(r, grad_m)
        rhs = -m * r
        try:
            if jac_g.shape[0] == jac_g.shape[1]:
                dz = np.linalg.solve(jac_g, rhs)
            else:
                dz, *_ = np.linalg.lstsq(jac_g, rhs, rcond=None)
        except np.linalg.LinAlgError:
            return DeflatedNewtonResult(
                z=None, stop_reason=NewtonStopReason.SINGULAR, residual_norm=raw_norm, n_iter=it
            )

        if step_cap is not None:
            dz = np.clip(dz, -step_cap, step_cap)
        z = z + dz

    return DeflatedNewtonResult(
        z=None,
        stop_reason=NewtonStopReason.MAX_ITER,
        residual_norm=float(np.linalg.norm(r_last)),
        n_iter=max_iter,
    )


def enumerate_roots(
    residual_fn: ResidualFn,
    seeds: Sequence[NDArray[np.float64]],
    *,
    jacobian_fn: JacobianFn | None = None,
    tol: float = 1e-10,
    max_iter: int = 50,
    p: float = 2.0,
    shift: float = 1.0,
    step_cap: NDArray[np.float64] | float | None = None,
    dedup_tol: float = 1e-6,
) -> list[NDArray[np.float64]]:
    """Run deflated Newton from each seed, deflating every root found so far.

    Returns the list of DISTINCT converged roots (deduplicated at
    ``dedup_tol``), in the order first found. A seed whose iterate diverges,
    stalls, or reconverges to an already-known root (rejected by the
    dedup check) contributes nothing -- this is the expected behavior for
    seeds outside any new root's basin, not an error.
    """
    roots: list[NDArray[np.float64]] = []
    for seed in seeds:
        result = newton_deflated(
            residual_fn,
            seed,
            roots,
            jacobian_fn=jacobian_fn,
            tol=tol,
            max_iter=max_iter,
            p=p,
            shift=shift,
            step_cap=step_cap,
        )
        z_found = result.z
        if z_found is None or result.stop_reason != NewtonStopReason.CONVERGED:
            continue
        is_new = not any(float(np.linalg.norm(z_found - r)) < dedup_tol for r in roots)
        if is_new:
            roots.append(z_found)
    return roots


__all__ = [
    "DeflatedNewtonResult",
    "JacobianFn",
    "NewtonStopReason",
    "ResidualFn",
    "deflation_factor",
    "enumerate_roots",
    "newton_deflated",
]
