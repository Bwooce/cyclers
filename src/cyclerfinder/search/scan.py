"""Parallel epoch x branch multi-start scan over the N-arc ballistic corrector.

The N-arc differential corrector (:mod:`cyclerfinder.search.correct`) is a
single-start root-find: from one ``(t0, ToF)`` seed it converges to the nearest
ballistic family. Reaching a *specific* family (e.g. the Jones VEM members, or
the ~6.4 km/s S1L1 Mars family the prototype's ``main()`` loop originally found)
requires sampling many starts across launch epoch and per-leg rev/branch
topology -- the density lever (plan / task #110, ``scripts/correct_s1l1_twoarc.py``
``main()``).

This module drives that outer grid in parallel. The documented constraint
(``optimize.py`` ~line 885): inner-DE ``workers=-1`` lost to per-evaluation
pickling of ``(Cell, Ephemeris)``. So we parallelise the **outer** grid instead:
each worker receives only primitives (body codes, epoch float, ToF-seed floats,
rev/branch tuples, the ephemeris *model string*) and constructs its own
:class:`~cyclerfinder.core.ephemeris.Ephemeris` inside the worker process. The
only objects crossing the process boundary are the small primitive payload and a
small frozen result dataclass -- never a live ``Ephemeris`` or ``Cell``.

Determinism: results are returned sorted by ``(max_residual_kms, t0_seed_sec)``
so a parallel run is bit-for-bit reproducible against a serial run of the same
grid.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import BallisticClosureResult, ballistic_correct


@dataclass(frozen=True)
class ScanPoint:
    """One multi-start grid point -- pure primitives, picklable to a worker.

    Carries everything a worker needs to construct its own ``Ephemeris`` and run
    :func:`~cyclerfinder.search.correct.ballistic_correct`. No live objects.
    """

    sequence: tuple[str, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    t0_seed_sec: float
    tof_seed_days: tuple[float, ...]
    period_sec: float
    vinf_cap: float
    slack_leg: int | None = None
    rp_factors: tuple[tuple[str, float], ...] | None = None
    tol_kms: float = 0.1
    residual_mode: str = "magnitude"


@dataclass(frozen=True)
class ScanResult:
    """A grid point paired with its ballistic-closure outcome."""

    point: ScanPoint
    result: BallisticClosureResult

    @property
    def closed(self) -> bool:
        return self.result.converged

    @property
    def max_residual_kms(self) -> float:
        return self.result.max_residual_kms

    @property
    def t0_seed_sec(self) -> float:
        return self.point.t0_seed_sec


def _evaluate_point(point: ScanPoint) -> ScanResult:
    """Worker entry point: build a fresh Ephemeris and run the corrector.

    Runs in a worker process. Only primitives crossed the boundary; the
    ``Ephemeris`` is constructed here, never pickled. The corrector already maps
    Lambert pathologies to a non-converged result, so this returns a result for
    every point rather than raising.
    """
    ephem = Ephemeris(model=_EPHEM_MODEL)
    rp_factors = dict(point.rp_factors) if point.rp_factors is not None else None
    result = ballistic_correct(
        sequence=point.sequence,
        per_leg_revs=point.per_leg_revs,
        per_leg_branch=point.per_leg_branch,
        t0_seed_sec=point.t0_seed_sec,
        tof_seed_days=point.tof_seed_days,
        period_sec=point.period_sec,
        ephem=ephem,
        vinf_cap=point.vinf_cap,
        rp_factors=rp_factors,
        slack_leg=point.slack_leg,
        tol_kms=point.tol_kms,
        residual_mode=point.residual_mode,
    )
    return ScanResult(point=point, result=result)


# The ephemeris model is process-global in each worker so it is initialised once
# per worker (via the pool initialiser) rather than pickled per point. The parent
# sets it through the initialiser argument; a module default keeps direct
# (serial) calls to ``_evaluate_point`` working in tests.
_EPHEM_MODEL: str = "astropy"


def _init_worker(model: str) -> None:
    """Pool initialiser: pin the ephemeris model for this worker process."""
    global _EPHEM_MODEL
    _EPHEM_MODEL = model


def _sorted(results: Iterable[ScanResult]) -> list[ScanResult]:
    """Deterministic order: ascending residual, then ascending t0 seed."""
    return sorted(results, key=lambda r: (r.max_residual_kms, r.t0_seed_sec))


def scan_serial(
    points: Sequence[ScanPoint],
    *,
    ephem_model: str = "astropy",
) -> list[ScanResult]:
    """Serial reference scan -- the determinism oracle for :func:`scan_parallel`.

    Evaluates every grid point in this process with one shared ``Ephemeris``.
    Returns results in the same deterministic order as the parallel path.
    """
    global _EPHEM_MODEL
    saved = _EPHEM_MODEL
    _EPHEM_MODEL = ephem_model
    try:
        results = [_evaluate_point(p) for p in points]
    finally:
        _EPHEM_MODEL = saved
    return _sorted(results)


def scan_parallel(
    points: Sequence[ScanPoint],
    *,
    ephem_model: str = "astropy",
    max_workers: int | None = None,
    closed_only: bool = False,
) -> list[ScanResult]:
    """Run the epoch x branch multi-start grid over a process pool.

    Each worker constructs its own ``Ephemeris(model=ephem_model)`` once (pool
    initialiser) and evaluates points with it; only the primitive
    :class:`ScanPoint` payload and the small :class:`ScanResult` cross the
    process boundary. Results are sorted by ``(max_residual_kms, t0_seed_sec)``
    so the output is reproducible and identical to :func:`scan_serial` on the
    same grid.

    Parameters
    ----------
    points:
        The multi-start grid (epoch x branch x ToF-seed).
    ephem_model:
        Ephemeris model string passed to every worker (``"astropy"`` for DE440).
    max_workers:
        Worker process count. ``None`` -> ``os.cpu_count()``.
    closed_only:
        When ``True``, drop non-converged points from the returned list (the
        full grid is still evaluated; this only filters the output).

    Returns
    -------
    list[ScanResult]
        Every grid point's outcome (or only the closed ones when
        ``closed_only``), deterministically ordered.
    """
    if not points:
        return []
    if max_workers is None:
        max_workers = os.cpu_count() or 1
    max_workers = max(1, min(max_workers, len(points)))

    if max_workers == 1:
        results = scan_serial(points, ephem_model=ephem_model)
    else:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_worker,
            initargs=(ephem_model,),
        ) as pool:
            results = _sorted(pool.map(_evaluate_point, points))

    if closed_only:
        results = [r for r in results if r.closed]
    return results


def build_epoch_branch_grid(
    *,
    sequence: tuple[str, ...],
    period_sec: float,
    vinf_cap: float,
    t0_seeds_sec: Sequence[float],
    branch_topologies: Sequence[tuple[tuple[int, ...], tuple[str, ...]]],
    tof_seed_days: Sequence[float],
    slack_leg: int | None = None,
    rp_factors: dict[str, float] | None = None,
    tol_kms: float = 0.1,
    residual_mode: str = "magnitude",
) -> list[ScanPoint]:
    """Materialise the epoch x branch multi-start grid as a list of ScanPoints.

    The Cartesian product of launch-epoch seeds and ``(per_leg_revs,
    per_leg_branch)`` topologies, each carrying the shared ToF seed and period.
    This mirrors the prototype's ``main()`` double loop (epoch offsets x EE
    branches, ``scripts/correct_s1l1_twoarc.py:162-181``) generalised to N arcs.

    ``tof_seed_days`` is the *free* (slack-eliminated) ToF seed the corrector
    expects, not the full per-leg seed; the caller drops the slack leg.
    """
    rp_tuple = tuple(sorted(rp_factors.items())) if rp_factors else None
    grid: list[ScanPoint] = []
    for revs, branch in branch_topologies:
        for t0 in t0_seeds_sec:
            grid.append(
                ScanPoint(
                    sequence=sequence,
                    per_leg_revs=revs,
                    per_leg_branch=branch,
                    t0_seed_sec=float(t0),
                    tof_seed_days=tuple(float(t) for t in tof_seed_days),
                    period_sec=float(period_sec),
                    vinf_cap=float(vinf_cap),
                    slack_leg=slack_leg,
                    rp_factors=rp_tuple,
                    tol_kms=float(tol_kms),
                    residual_mode=residual_mode,
                )
            )
    return grid
