"""Free-Jacobi extension of the general periodic-orbit return-map residual (#523 rework).

``cr3bp_general_periodic.correct_general_periodic`` solves the 2x2 return-map
fixed-point residual ``R = [x_c - x0, xdot_c - xdot0]`` at a FIXED Jacobi
constant, with an analytic STM Jacobian ``dR/d(x0, xdot0)``
(``_residual_jac_stm``). This module extends that Jacobian with a third,
equally analytic column ``dR/dC`` -- free, since it reuses the SAME STM
already computed for the 2x2 case, no extra propagation needed -- turning
the residual into a co-dimension-1 curve in the 3-vector ``(x0, xdot0, C)``:
exactly ``search/pseudo_arclength.py``'s target case (``M=N-1``).

Motivation (#523): the original brute-force grid search re-certified every
coarse-enumerator candidate independently via 60-100s per-candidate chained
Newton passes (``scripts/run_523_earth_coorbital_search.py``), making the
full 1,890-point grid impractical (~32-52 CPU-hours estimated). Once a
SINGLE orbit is certified, continuation along this free-C curve reuses that
convergence: each new family member near an already-converged point needs
only a few cheap Newton iterations from a good tangent-based predictor,
instead of a cold multi-pass certification from a coarse grid guess.
Measured this session: ~1.4s/continuation-step (dominated by two STM
propagations over the ~80-nondim-TU integration horizon), vs. 60-100s/
candidate for the original brute-force certification -- a ~50-70x per-point
speedup, and a much larger effective speedup versus grid re-sampling since
continuation does not re-discover already-covered family members.

The extra derivative, from ``ydot0(x0, xdot0, C) = sign * sqrt(Ubar_axis(x0)
- C - xdot0^2)`` (the same relation ``_ydot0_general`` implements):

    d(ydot0)/dC = -1 / (2 * ydot0)

projected through the STM exactly as the existing ``(x0, xdot0)`` columns are
in ``_residual_jac_stm`` (crossing-time variation held out so ``y_c = 0``).

Pure: reuses ``cr3bp_general_periodic``'s private ``_return_crossing_stm`` /
``_ydot0_general`` (same package, not a public API duplication) plus
``core.cr3bp`` and ``search.cr3bp_periodic``'s pseudo-potential gradient.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.cr3bp_general_periodic import _return_crossing_stm, _ydot0_general

ResidualFreeCFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]
JacobianFreeCFn = Callable[[NDArray[np.float64]], NDArray[np.float64] | None]


def residual_and_jacobian_free_c(
    x0: float,
    xdot0: float,
    c_target: float,
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """Return ``(R, J)`` with ``R`` the 2-component return-map residual and
    ``J = dR/d(x0, xdot0, C)`` (2x3), or ``None`` if infeasible.
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

    dubar_axis_dx = -2.0 * cp._ubar_grad_x_at_axis(x0, mu)
    dydot0_dx0 = dubar_axis_dx / (2.0 * ydot0)
    dydot0_dxdot0 = -xdot0 / ydot0
    dydot0_dc = -1.0 / (2.0 * ydot0)

    dx0_col = np.zeros(6)
    dx0_col[0] = 1.0
    dx0_col[4] = dydot0_dx0
    dxdot0_col = np.zeros(6)
    dxdot0_col[3] = 1.0
    dxdot0_col[4] = dydot0_dxdot0
    dc_col = np.zeros(6)
    dc_col[4] = dydot0_dc

    fdot = cr3bp.cr3bp_eom(t_c, state_c, mu)
    ydotc = float(fdot[1])
    if abs(ydotc) < 1e-14:
        return None

    jac = np.zeros((2, 3))
    for col, dcol in enumerate((dx0_col, dxdot0_col, dc_col)):
        dxc = stm @ dcol
        proj = float(dxc[1]) / ydotc
        dxc_constr = dxc - proj * fdot
        jac[0, col] = float(dxc_constr[0]) - (1.0 if col == 0 else 0.0)
        jac[1, col] = float(dxc_constr[3]) - (1.0 if col == 1 else 0.0)
    return res, jac


def make_residual_and_jacobian_fns(
    mu: float,
    sign: float,
    half_crossings: int,
    t_hi_frac: float,
    initial_period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[ResidualFreeCFn, JacobianFreeCFn]:
    """Build ``(residual_fn, jacobian_fn)`` closures for
    ``search.pseudo_arclength.continue_curve``, given ``z = (x0, xdot0, C)``.

    Caches the last ``(z, result)`` pair: ``pseudo_arclength.py`` calls
    ``residual_fn(z)`` and ``jacobian_fn(z)`` SEPARATELY for the same ``z``
    within one Newton iteration (``_correct``'s loop) -- without caching,
    that would compute the SAME expensive STM propagation twice per
    iteration.

    ``t_hi = t_hi_frac * initial_period`` is FIXED for the life of this
    closure pair (derived once from the seed orbit's own converged period,
    not the original script's generic ``T_MAX`` guess) -- adequate as long
    as the family's period does not grow enough along the walked ``C``-range
    to exceed it (if it does, ``_return_crossing_stm`` returns ``None``
    cleanly rather than silently misidentifying a later crossing; restart
    continuation from a freshly re-certified point with an updated
    ``initial_period`` if that happens, mirroring
    ``run_523_earth_coorbital_search.py``'s own period-guess chaining).
    """
    t_hi = t_hi_frac * abs(initial_period)
    cache: dict[
        tuple[float, float, float], tuple[NDArray[np.float64], NDArray[np.float64]] | None
    ] = {}

    def _evaluate(z: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
        key = (float(z[0]), float(z[1]), float(z[2]))
        if key in cache:
            return cache[key]
        cache.clear()  # only ever need the single most-recent evaluation
        out = residual_and_jacobian_free_c(
            key[0], key[1], key[2], mu, sign, half_crossings, t_hi, rtol=rtol, atol=atol
        )
        cache[key] = out
        return out

    def residual_fn(z: NDArray[np.float64]) -> NDArray[np.float64] | None:
        out = _evaluate(z)
        if out is None:
            return None
        return out[0]

    def jacobian_fn(z: NDArray[np.float64]) -> NDArray[np.float64] | None:
        out = _evaluate(z)
        if out is None:
            return None
        return out[1]

    return residual_fn, jacobian_fn


__all__ = ["make_residual_and_jacobian_fns", "residual_and_jacobian_free_c"]
