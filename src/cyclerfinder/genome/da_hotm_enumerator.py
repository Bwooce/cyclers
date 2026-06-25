"""Global multi-rev fixed-point enumerator over the section domain (#450 Task 4).

Given a :class:`~cyclerfinder.genome.da_hotm_backend.SectionMap` backend, a target
revolution count ``n``, and a 2-D ``(x, xdot)`` domain box on the Poincare section,
sweep the section-map fixed-point residual ``||P^n(s) - s||`` over a grid, isolate
the sub-tolerance cells, cluster them into distinct basins, and emit one coarse
candidate IC per basin. This is the GLOBAL, domain-covering enumeration that
seed-local continuation structurally cannot do (design draft §2): it finds every
fixed point in the box regardless of whether a sourced seed or a continuation path
to it exists.

The enumerator is deliberately COARSE and cheap -- the section-residual gate
(design draft §3 step 1). The expensive 1e-12 certification stays in the existing
``correct_general_periodic`` corrector, which the driver (Task 6) hands each
emitted candidate.

Backend-agnostic: consumes only the ``SectionMap`` interface, so the sampling and
Taylor-map backends are swappable here.

Pure: math / numpy + the backend.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import (
    DASectionMap,
    SectionMap,
    SectionPoint,
)


@dataclass(frozen=True)
class DomainBox:
    """An axis-aligned ``(x, xdot)`` sub-rectangle of the Poincare section."""

    x_lo: float
    x_hi: float
    xdot_lo: float
    xdot_hi: float

    def __post_init__(self) -> None:
        if not (self.x_lo < self.x_hi and self.xdot_lo < self.xdot_hi):
            raise ValueError(f"DomainBox: degenerate box {self}")


@dataclass(frozen=True)
class FixedPointCandidate:
    """A coarse fixed-point candidate emitted by the enumerator.

    ``(x0, xdot0)`` is the section IC at fixed Jacobi ``c_target``; ``n`` the
    revolution count; ``residual`` the section-map residual ``||P^n(s) - s||`` at
    emission (the coarse gate value, NOT a certified closure).
    """

    x0: float
    xdot0: float
    c_target: float
    n: int
    residual: float
    moved: float = 0.0  # distance the Taylor finder moved this from its reference


def _is_local_min(grid: np.ndarray, i: int, j: int) -> bool:
    """True if cell (i, j) is <= its (up to 8) finite neighbours."""
    v = grid[i, j]
    if not np.isfinite(v):
        return False
    ni, nj = grid.shape
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            a, b = i + di, j + dj
            if 0 <= a < ni and 0 <= b < nj:
                w = grid[a, b]
                if np.isfinite(w) and w < v:
                    return False
    return True


def enumerate_fixed_points(
    backend: SectionMap,
    box: DomainBox,
    n: int,
    *,
    residual_tol: float,
    grid: tuple[int, int] = (41, 31),
    dedup_radius: float = 0.01,
) -> list[FixedPointCandidate]:
    """Enumerate ``P^n`` fixed points over ``box`` and emit coarse candidate ICs.

    Parameters
    ----------
    backend:
        Any :class:`SectionMap` (sampling or Taylor-map).
    box:
        The ``(x, xdot)`` sub-rectangle to sweep.
    n:
        Revolution count -- fixed points of the n-th iterate ``P^n``.
    residual_tol:
        Coarse section-residual gate: a cell is a candidate only if
        ``||P^n(s) - s|| <= residual_tol``.
    grid:
        ``(n_x, n_xdot)`` sample counts over the box.
    dedup_radius:
        Section-distance below which two sub-tol local minima are the SAME basin;
        the lower-residual representative is kept.

    Returns
    -------
    list[FixedPointCandidate]
        One coarse candidate per distinct basin, sorted by ascending residual.
        Empty if no cell is sub-tolerance (a legitimate negative for the box).
    """
    nx, nxd = grid
    xs = np.linspace(box.x_lo, box.x_hi, nx)
    xds = np.linspace(box.xdot_lo, box.xdot_hi, nxd)
    res = np.full((nx, nxd), np.inf, dtype=np.float64)
    for i, x in enumerate(xs):
        for j, xd in enumerate(xds):
            res[i, j] = backend.residual(SectionPoint(x=float(x), xdot=float(xd)), n)

    # Sub-tolerance local minima are candidate basin centres.
    raw: list[tuple[float, float, float]] = []  # (residual, x, xdot)
    for i in range(nx):
        for j in range(nxd):
            if res[i, j] <= residual_tol and _is_local_min(res, i, j):
                raw.append((float(res[i, j]), float(xs[i]), float(xds[j])))
    raw.sort(key=lambda t: t[0])  # lowest residual first

    # Greedy dedup: keep the lowest-residual rep per basin (radius dedup_radius).
    kept: list[tuple[float, float, float]] = []
    for r, x, xd in raw:
        if all(
            (x - kx) ** 2 + (xd - kxd) ** 2 > dedup_radius * dedup_radius for _, kx, kxd in kept
        ):
            kept.append((r, x, xd))

    return [
        FixedPointCandidate(x0=x, xdot0=xd, c_target=backend.c_target, n=int(n), residual=r)
        for r, x, xd in kept
    ]


def taylor_enumerate(
    backend: DASectionMap,
    box: DomainBox,
    n: int,
    *,
    order: int = 2,
    h: float = 3e-4,
    samples: int = 6,
    ref_grid: tuple[int, int] = (13, 25),
    dedup_radius: float = 2e-3,
    max_iter: int = 10,
    stop_when: Callable[[float, float, float], bool] | None = None,
) -> list[FixedPointCandidate]:
    """Enumerate ``P^n`` fixed points over ``box`` with the Taylor-map backend.

    Sweeps a COARSE grid of references over the box; from each, runs the iterated
    Taylor-map fixed-point finder (:meth:`DASectionMap.taylor_fixed_point`) and
    collects the converged section points; dedups them into distinct basins. This
    is the GLOBAL enumeration that reaches the strongly-unstable multi-rev families
    (P5g') a brute-force grid cannot (their section basin is narrower than any
    feasible grid spacing -- see the #450 decision note). Each emitted candidate is
    COARSE (the FD-Taylor floor ~3e-4 for P5g'); the corrector micro-multistart
    (``search.da_hotm_close.close_candidate``) certifies it to 1e-12.

    Returns one :class:`FixedPointCandidate` per distinct Taylor fixed-point basin,
    with ``residual`` set to the (fragile) section residual at the candidate, or
    ``+inf`` if the float compose walls off there (still a valid coarse candidate
    for the corrector). Empty if the Taylor finder never converges in the box.

    ``stop_when(x, xdot, moved)`` -- optional early-exit predicate evaluated on
    each converged candidate; the sweep returns immediately (with the candidates
    found so far) when it is True. The driver uses this to stop at the first
    closable candidate instead of sweeping the whole grid.
    """
    nx, nxd = ref_grid
    xs = np.linspace(box.x_lo, box.x_hi, nx)
    xds = np.linspace(box.xdot_lo, box.xdot_hi, nxd)
    found: list[tuple[float, float, float, float]] = []  # (res_or_inf, x, xdot, moved)
    stop = False
    for x in xs:
        if stop:
            break
        for xd in xds:
            ref = SectionPoint(x=float(x), xdot=float(xd))
            try:
                fp = backend.taylor_fixed_point(
                    ref, n=n, order=order, h=h, samples=samples, max_iter=max_iter
                )
            except (ValueError, RuntimeError):
                continue
            if not (box.x_lo <= fp.x <= box.x_hi and box.xdot_lo <= fp.xdot <= box.xdot_hi):
                continue
            moved = math.hypot(fp.x - ref.x, fp.xdot - ref.xdot)
            res = backend.residual(fp, n)
            found.append((res, fp.x, fp.xdot, moved))
            if stop_when is not None and stop_when(fp.x, fp.xdot, moved):
                stop = True
                break

    # Dedup into basins; keep the rep that MOVED most (converged, not a stalled
    # reference) within each basin, residual as a tiebreak.
    found.sort(key=lambda t: (-t[3], not math.isfinite(t[0]), t[0]))
    kept: list[tuple[float, float, float, float]] = []
    for r, x, xd, mv in found:
        if all(
            (x - kx) ** 2 + (xd - kxd) ** 2 > dedup_radius * dedup_radius for _, kx, kxd, _ in kept
        ):
            kept.append((r, x, xd, mv))
    return [
        FixedPointCandidate(
            x0=x, xdot0=xd, c_target=backend.c_target, n=int(n), residual=r, moved=mv
        )
        for r, x, xd, mv in kept
    ]


def recover_png_candidate(
    backend: DASectionMap,
    box: DomainBox,
    n: int,
    *,
    period_guess: float | None = None,
    ref_grid: tuple[int, int] = (13, 25),
    close_reach: float = 2e-4,
) -> FixedPointCandidate | None:
    """Emit the best coarse Png' candidate from ``box`` (lane-recovery, Task 5).

    Runs :func:`taylor_enumerate` over the band, then ranks the distinct coarse
    candidates by how far the Taylor finder MOVED them from their grid reference
    (a candidate that the iterated map pulled toward a true fixed point moved more
    than a stalled off-family reference) and returns the one nearest the band's
    converged cluster. Returns the COARSE candidate -- an OUTPUT of the global
    sweep, NOT a handed-in IC. The corrector micro-multistart
    (``search.da_hotm_close.close_candidate``) certifies it to the published P5g'.

    ``period_guess`` is unused here (the closure horizon is set by the caller); it
    is accepted so the driver can pass it through uniformly.

    Sweeps a FINE reference grid (mesh ~1e-3) over the band -- fine enough that a
    node lands within the Taylor convergence basin (~2e-3) of the family's fixed
    point, which a coarser grid misses (decision-note finding). Each node runs the
    iterated Taylor finder. There is NO smooth scalar residual near this
    strongly-unstable multi-rev fixed point (the crossing structure flips under any
    perturbation, so the section residual is +inf almost everywhere except exactly
    on it), so candidates are ranked by how far the Taylor finder MOVED them from
    their grid node: a stalled off-family reference barely moves, while the genuine
    fixed point is PULLED far (toward P5g'). The most-moved candidate is returned;
    the corrector micro-multistart certifies it to the published P5g'. The sweep
    early-exits once a candidate has moved past ``close_reach`` from MULTIPLE nodes
    (a robust convergence signal), bounding the wall time.
    """
    _ = period_guess
    nx, nxd = ref_grid
    xs = np.linspace(box.x_lo, box.x_hi, nx)
    xds = np.linspace(box.xdot_lo, box.xdot_hi, nxd)
    best: FixedPointCandidate | None = None
    best_moved = -1.0
    hits_at_cluster = 0
    for x in xs:
        for xd in xds:
            ref = SectionPoint(x=float(x), xdot=float(xd))
            try:
                fp = backend.taylor_fixed_point(ref, n=n, order=2, h=3e-4, samples=6, max_iter=12)
            except (ValueError, RuntimeError):
                continue
            if not (box.x_lo <= fp.x <= box.x_hi and box.xdot_lo <= fp.xdot <= box.xdot_hi):
                continue
            moved = math.hypot(fp.x - ref.x, fp.xdot - ref.xdot)
            if moved > best_moved:
                best_moved = moved
                best = FixedPointCandidate(
                    x0=fp.x,
                    xdot0=fp.xdot,
                    c_target=backend.c_target,
                    n=int(n),
                    residual=_robust_section_residual(backend, fp, n),
                    moved=moved,
                )
            # Cluster confirmation: count nodes whose Taylor fixed point lands near
            # the current best (the family attracts a basin of references). Early
            # exit once the cluster is confirmed AND the best has moved well past
            # the close reach -- the family is surfaced.
            if best is not None and math.hypot(fp.x - best.x0, fp.xdot - best.xdot0) < 1e-3:
                hits_at_cluster += 1
                if hits_at_cluster >= 3 and best_moved > 5.0 * close_reach:
                    return best
    return best


def _robust_section_residual(backend: DASectionMap, s: SectionPoint, n: int) -> float:
    """Return-map residual to the 2n-th y=0 crossing (finite where structure intact).

    The fragile float ``compose`` walls off to +inf for strongly-unstable orbits;
    this propagates the lifted IC once to the 2n-th y=0 crossing and measures
    ``||(x, xdot)_{2n} - (x, xdot)_0||``, which stays finite as long as 2n
    crossings exist -- a better ranking signal for picking the closable candidate.
    """
    try:
        state0 = backend.lift(s)
    except ValueError:
        return float("inf")
    mu = float(backend.system.mu)

    def _y_event(t: float, y: np.ndarray, _mu: float) -> float:
        return float(y[1])

    _y_event.direction = 0.0  # type: ignore[attr-defined]
    t_hi = 8.0 * n
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_hi),
        state0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        events=_y_event,
    )
    crossings = [
        (t, y) for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True) if t > 1e-6 * t_hi
    ]
    if len(crossings) < 2 * n:
        return float("inf")
    _, yf = crossings[2 * n - 1]
    return float(math.hypot(yf[0] - s.x, yf[3] - s.xdot))


__all__ = [
    "DomainBox",
    "FixedPointCandidate",
    "enumerate_fixed_points",
    "recover_png_candidate",
    "taylor_enumerate",
]
