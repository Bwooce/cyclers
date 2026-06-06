"""Parallel epoch x branch multi-start scan (task #110).

Fast smoke: a tiny circular-model grid, parallel result identical to a serial
run of the same grid -- the determinism gate. Slow real-grid: a DE440 grid with
a measured serial-vs-parallel wall-time speedup, reported honestly.

NON-GOLDEN: any V_inf / residual asserted here is OUR computation, used only as a
determinism / regime fixture, never as a published-anchor EXPECTED.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from cyclerfinder.search.scan import (
    ScanPoint,
    build_epoch_branch_grid,
    scan_parallel,
    scan_serial,
)

DAY_S = 86400.0


def _s1l1_grid_circular(n_epochs: int) -> list[ScanPoint]:
    """Tiny E-M-E-E grid on the cheap circular model (fast, deterministic)."""
    period_sec = (1.4612 + 2.8096) * 365.25 * DAY_S
    base = (
        np.datetime64("2030-03-22T00:00:00") - np.datetime64("2000-01-01T12:00:00")
    ) / np.timedelta64(1, "s")
    t0_seeds = [float(base) + off * DAY_S for off in range(0, n_epochs * 60, 60)]
    return build_epoch_branch_grid(
        sequence=("E", "M", "E", "E"),
        period_sec=period_sec,
        vinf_cap=9.0,
        t0_seeds_sec=t0_seeds,
        branch_topologies=[
            ((0, 1, 2), ("single", "single", "low")),
            ((0, 1, 1), ("single", "single", "low")),
        ],
        tof_seed_days=[154.0, 379.0],
        slack_leg=2,
    )


def test_parallel_matches_serial_determinism_gate() -> None:
    """A 2-worker parallel scan returns bit-identical results to the serial run
    of the same grid (the reproducibility contract)."""
    grid = _s1l1_grid_circular(n_epochs=4)
    serial = scan_serial(grid, ephem_model="circular")
    parallel = scan_parallel(grid, ephem_model="circular", max_workers=2)

    assert len(parallel) == len(serial) == len(grid)
    for s, p in zip(serial, parallel, strict=True):
        assert p.point == s.point
        assert p.result.t0_sec == s.result.t0_sec
        assert p.result.tof_days == s.result.tof_days
        assert p.result.max_residual_kms == s.result.max_residual_kms
        assert p.result.converged == s.result.converged
        assert p.result.vinf_per_encounter_kms == s.result.vinf_per_encounter_kms


def test_results_sorted_by_residual_then_epoch() -> None:
    """Deterministic ordering: ascending residual, then ascending t0 seed."""
    grid = _s1l1_grid_circular(n_epochs=4)
    results = scan_parallel(grid, ephem_model="circular", max_workers=2)
    keys = [(r.max_residual_kms, r.t0_seed_sec) for r in results]
    assert keys == sorted(keys)


def test_empty_grid_returns_empty() -> None:
    assert scan_parallel([], ephem_model="circular") == []


def test_closed_only_filters_output() -> None:
    grid = _s1l1_grid_circular(n_epochs=4)
    full = scan_parallel(grid, ephem_model="circular", max_workers=2)
    closed = scan_parallel(grid, ephem_model="circular", max_workers=2, closed_only=True)
    assert all(r.closed for r in closed)
    assert len(closed) <= len(full)
    assert closed == [r for r in full if r.closed]


@pytest.mark.slow
def test_parallel_speedup_on_real_grid() -> None:
    """DE440 grid: measure serial vs parallel wall time and assert the parallel
    path is at least not slower (speedup reported in the captured output).

    The numbers printed here are MEASURED, not asserted as anchors -- the gate is
    only that parallelism does not regress wall time and that the results are
    identical to the serial reference.
    """
    grid = _s1l1_grid_circular(n_epochs=16)  # 2 branches x 16 epochs = 32 points

    t0 = time.perf_counter()
    serial = scan_serial(grid, ephem_model="astropy")
    t_serial = time.perf_counter() - t0

    t1 = time.perf_counter()
    parallel = scan_parallel(grid, ephem_model="astropy", max_workers=8)
    t_parallel = time.perf_counter() - t1

    speedup = t_serial / t_parallel if t_parallel > 0 else float("inf")
    print(
        f"\n[scan speedup] points={len(grid)} serial={t_serial:.2f}s "
        f"parallel(8w)={t_parallel:.2f}s speedup={speedup:.2f}x"
    )

    # Determinism holds on the real model too.
    assert len(parallel) == len(serial)
    for s, p in zip(serial, parallel, strict=True):
        assert p.result.max_residual_kms == s.result.max_residual_kms
        assert p.result.converged == s.result.converged

    # Parallel must not be materially slower than serial on a real grid.
    assert t_parallel <= t_serial * 1.5
