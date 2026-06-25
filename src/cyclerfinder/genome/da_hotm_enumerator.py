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

from dataclasses import dataclass

import numpy as np

from cyclerfinder.genome.da_hotm_backend import SectionMap, SectionPoint


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


__all__ = ["DomainBox", "FixedPointCandidate", "enumerate_fixed_points"]
