"""Empirical transit-branch manifold grids for QP-tori (#548).

The linking-number heteroclinic screen (:mod:`genome.qp_torus_heteroclinic`)
builds its stable/unstable manifold-endpoint grids via
:func:`genome.qp_torus_manifold.torus_manifold_grid`, which selects WHICH of
the two manifold branches (``+`` / ``-`` along the local unstable/stable
eigenvector) to propagate using a ``vec[0] * sign`` heuristic -- i.e. it trusts
the SIGN of the x-component of the STM eigenvector to point toward the transit
direction. Task #547 flagged this as a plausible-but-unconfirmed weak point for
genuine 3D quasi-halo eigenvectors (the x-component can be small or flip across
the grid, silently mixing the two manifold sheets into one endpoint grid and
corrupting the closed torus-map curve the linking number needs).

This module removes that dependency. It adapts #547's EMPIRICAL transit
classification (``genome.transit_manifold.classify_unstable_branch``, built and
validated for the planar-Lyapunov case) to the 3D quasi-halo torus grid: at
EVERY torus point it propagates BOTH signed perturbations and keeps the one that
actually reaches the surface of section first -- the genuine transit branch,
decided by the dynamics rather than an eigenvector-component sign. The
non-transit branch (which falls back and never crosses, or crosses only after a
long excursion) is discarded. The result is a :class:`ManifoldGrid` in the same
shape the linking-number scan already consumes, so it is a drop-in replacement
for :func:`torus_manifold_grid` in that pipeline.

Pure CR3BP + qp_tori/qp_torus_manifold dependency only.
"""

from __future__ import annotations

import concurrent.futures
import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus
from cyclerfinder.genome.qp_torus_manifold import (
    ManifoldGrid,
    local_stability,
    torus_point_stm,
)


@dataclass(frozen=True)
class _Crossing:
    t_abs: float
    state: NDArray[np.float64]


def _first_crossing(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    *,
    direction: float,
    surface_x: float,
    t_max: float,
    rtol: float,
    atol: float,
) -> _Crossing | None:
    """First crossing of ``x = surface_x`` under forward (``direction>0``) or
    backward (``direction<0``) flow, or ``None`` if none within ``|t_max|``.
    Returns the crossing state and the ABSOLUTE crossing time (so the earliest
    crossing across two candidate branches can be compared directly).
    """

    def _x_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[0] - surface_x)

    _x_event.terminal = True  # type: ignore[attr-defined]
    _x_event.direction = 0.0  # type: ignore[attr-defined]

    t_span = (0.0, math.copysign(t_max, direction))
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        np.asarray(state0, float),
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_x_event,
        max_step=abs(t_max) / 20.0,
    )
    if sol.t_events is None or len(sol.t_events[0]) == 0:
        return None
    return _Crossing(
        t_abs=abs(float(sol.t_events[0][0])),
        state=np.asarray(sol.y_events[0][0], dtype=np.float64),
    )


def _transit_crossing(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    vec: NDArray[np.float64],
    *,
    eps: float,
    direction: float,
    surface_x: float,
    t_max: float,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], int] | None:
    """Empirically choose the transit branch at one torus point.

    Propagates BOTH ``state0 + eps*vec`` and ``state0 - eps*vec`` and returns
    the crossing state of whichever branch reaches ``x = surface_x`` FIRST
    (the direct transit trajectory), together with the chosen sign (``+1`` /
    ``-1``). Returns ``None`` if NEITHER branch crosses within ``|t_max|`` (a
    locally non-transit point). No eigenvector-component sign is trusted.
    """
    best_state: NDArray[np.float64] | None = None
    best_sign = 0
    best_t = math.inf
    for s in (+1, -1):
        seed = state0 + float(s) * eps * vec
        cr = _first_crossing(
            system,
            seed,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=rtol,
            atol=atol,
        )
        if cr is not None and cr.t_abs < best_t:
            best_t = cr.t_abs
            best_state = cr.state
            best_sign = s
    if best_state is None:
        return None
    return best_state, best_sign


def _transit_row(
    args: tuple[
        int, float, NDArray[np.float64], QPTorus, str, float, float, float, float, float, float
    ],
) -> tuple[int, NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_], NDArray[np.int8]]:
    (
        i,
        theta_long,
        thetas_lat,
        torus,
        branch,
        eps,
        surface_x,
        t_max,
        hyperbolicity_tol,
        rtol,
        atol,
    ) = args
    n_lat = len(thetas_lat)
    row_origins = np.zeros((n_lat, 2))
    row_endpoints = np.full((n_lat, 6), np.nan)
    row_hyperbolic = np.zeros(n_lat, dtype=bool)
    row_sign = np.zeros(n_lat, dtype=np.int8)
    direction = 1.0 if branch == "unstable" else -1.0
    for j, theta_trans in enumerate(thetas_lat):
        row_origins[j, 0] = theta_long
        row_origins[j, 1] = theta_trans
        state, stm = torus_point_stm(
            torus, float(theta_long), float(theta_trans), rtol=rtol, atol=atol
        )
        stab = local_stability(state, stm, hyperbolicity_tol=hyperbolicity_tol)
        vec = stab.vec_u if branch == "unstable" else stab.vec_s
        if vec is None:
            continue
        row_hyperbolic[j] = True
        res = _transit_crossing(
            torus.system,
            state,
            vec,
            eps=eps,
            direction=direction,
            surface_x=surface_x,
            t_max=t_max,
            rtol=rtol,
            atol=atol,
        )
        if res is not None:
            row_endpoints[j, :] = res[0]
            row_sign[j] = res[1]
    return i, row_origins, row_endpoints, row_hyperbolic, row_sign


def transit_torus_manifold_grid(
    torus: QPTorus,
    *,
    n_long: int,
    n_lat: int,
    branch: str,
    surface_x: float,
    t_max: float,
    eps: float = 1e-5,
    hyperbolicity_tol: float = 1e-4,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    max_workers: int = 12,
) -> tuple[ManifoldGrid, NDArray[np.int8]]:
    """Build a manifold-endpoint grid with EMPIRICAL per-point transit-branch
    selection (see module docstring).

    Same output contract as :func:`torus_manifold_grid` (a
    :class:`ManifoldGrid` whose ``endpoints`` are the surface-of-section
    crossing states, ``nan`` where no branch transited). Additionally returns
    the ``(n_long, n_lat)`` grid of the chosen branch signs (``+1`` / ``-1`` /
    ``0`` where non-transit) for diagnostics -- a large, coherent single-sign
    region indicates a clean transit sheet; a salt-and-pepper sign field warns
    that the two sheets are being mixed (the failure mode #547 flagged).
    """
    if branch not in ("stable", "unstable"):
        raise ValueError(f"branch must be 'stable' or 'unstable', got {branch!r}")
    thetas_long = 2.0 * math.pi * np.arange(n_long) / n_long
    thetas_lat = 2.0 * math.pi * np.arange(n_lat) / n_lat

    origins = np.zeros((n_long, n_lat, 2), dtype=np.float64)
    endpoints = np.full((n_long, n_lat, 6), np.nan, dtype=np.float64)
    hyperbolic = np.zeros((n_long, n_lat), dtype=np.bool_)
    signs = np.zeros((n_long, n_lat), dtype=np.int8)

    tasks = [
        (
            i,
            float(tl),
            thetas_lat,
            torus,
            branch,
            eps,
            surface_x,
            t_max,
            hyperbolicity_tol,
            rtol,
            atol,
        )
        for i, tl in enumerate(thetas_long)
    ]
    if n_long >= 4 and max_workers > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(max_workers, n_long)) as ex:
            results = list(ex.map(_transit_row, tasks))
    else:
        results = [_transit_row(t) for t in tasks]

    for i, ro, re, rh, rs in results:
        origins[i] = ro
        endpoints[i] = re
        hyperbolic[i] = rh
        signs[i] = rs

    return ManifoldGrid(origins=origins, endpoints=endpoints, hyperbolic=hyperbolic), signs


__all__ = ["transit_torus_manifold_grid"]
