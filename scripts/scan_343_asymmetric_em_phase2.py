"""Phase 2 asymmetric-corrector scan (Earth-Moon, #343).

Phase 2 of #284. Rewires #284's per-cell guard chain through:

* :func:`cyclerfinder.parallel.parallel_sweep.parallel_sweep` (#321 substrate)
  -- 5-8x measured speedup on embarrassingly-parallel single-cell sweeps.
* :func:`cyclerfinder.search.single_orbit_prioritizer.score_single_orbit`
  (#310 adapter) -- closes the "prioritizer needs pair-shaped input" gap that
  blocked #284 Phase 1 from ranking surviving candidates.

The grid expands #284's 320-cell coverage to 1,944 cells (6.1x more) under the
same wall-time budget, with prioritizer Tier-0 ΔV-to-surrogate as the rank
signal. (A first Phase 2 trial targeted ~5,800 cells but the corrector cost on
the higher (k1, k2) bands under live CPU contention ran past the harness
budget; the 1,944-cell grid keeps the (k1, k2) coverage expansion -- the axis
#284 most under-sampled -- and trims hc=8 + the +ydot0 sign, both of which
contributed little new converged territory in Phase 1.) Output:
``data/scan_343_asymmetric_em_phase2.jsonl``.

Discipline (per the task brief):
- NO catalogue writeback.
- NO novelty claims in commit messages.
- SILVER = pass corrector + topology + DOP853 closure + ML low-FP.
- Novelty-claimable = SILVER + literature-fresh (offline; necessary-not-sufficient).
- Per-cell processing lives in ``asymmetric_novel_scan_parallel.process_cell``
  (top-level closure, pickle-safe under loky).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
from cyclerfinder.search.asymmetric_novel_scan_parallel import (
    AsymCell,
    process_cell,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "data" / "scan_343_asymmetric_em_phase2.jsonl"

# Phase 2 expanded grid. Phase 1 was 5x4x4x2x1x2 = 320 cells. Phase 2 below is
# 9x6x6x3x1x2 = 1,944 cells -- 6.1x more coverage. Expansion targets axes #284
# flagged as under-covered: (k1, k2) bands missing 3:4 / 4:3 / 5:2 / 2:5, plus
# finer C/x0/xdot0 sampling.
DEFAULT_K_TARGETS: tuple[tuple[int, int], ...] = (
    (1, 1),
    (1, 2),
    (2, 1),
    (3, 2),
    (2, 3),
    (3, 4),
    (4, 3),
    (5, 2),
    (2, 5),
)
DEFAULT_C_GRID = (3.06, 3.10, 3.13, 3.15, 3.17, 3.19)  # 6 Jacobi levels
DEFAULT_X0_GRID = (-0.92, -0.85, -0.80, -0.75, -0.70, -0.65)  # 6 x0 seeds
DEFAULT_XDOT0_SEEDS = (-0.05, 0.0, 0.05)  # 3 xdot0 seeds (now includes 0)
# ydot0_sign kept at (-1.0,) per #284 Phase 1 convention. Both-signs was probed
# in an early Phase 2 trial; corrector basin showed no sign-asymmetry that the
# other axes' expansion didn't already cover. Reverting to one sign halves the
# wall time and preserves the comparison axis to Phase 1.
DEFAULT_YDOT0_SIGNS = (-1.0,)
DEFAULT_HALF_CROSSINGS = (4, 6)  # match #284 Phase 1; hc=8 slowest, mostly null

# Reference: serial baseline cost from #284 Phase 1 (135 converged / 320 cells
# in 4,810s ~ 15 s/converged). Phase 2 budget: 16-core, 8-12x speedup -> 5,832
# cells in 10-25 min wall.
PROGRESS_EVERY_SEC = 60.0


def _log(msg: str) -> None:
    ts = datetime.now(UTC).strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def _build_grid(
    k_targets: tuple[tuple[int, int], ...] = DEFAULT_K_TARGETS,
    c_grid: tuple[float, ...] = DEFAULT_C_GRID,
    x0_grid: tuple[float, ...] = DEFAULT_X0_GRID,
    xdot0_seeds: tuple[float, ...] = DEFAULT_XDOT0_SEEDS,
    ydot0_signs: tuple[float, ...] = DEFAULT_YDOT0_SIGNS,
    half_crossings: tuple[int, ...] = DEFAULT_HALF_CROSSINGS,
) -> list[AsymCell]:
    cells: list[AsymCell] = []
    for k1, k2 in k_targets:
        for c in c_grid:
            for x0 in x0_grid:
                for xd in xdot0_seeds:
                    for sign in ydot0_signs:
                        for hc in half_crossings:
                            cells.append(
                                AsymCell(
                                    k1=int(k1),
                                    k2=int(k2),
                                    jacobi_request=float(c),
                                    seed_x0=float(x0),
                                    seed_xdot0=float(xd),
                                    ydot0_sign=float(sign),
                                    half_crossings=int(hc),
                                )
                            )
    return cells


# Dedupe keys for cataloguing one row per "physical orbit" at grid resolution.
_DEDUP_X0_DP = 4
_DEDUP_C_DP = 4
_DEDUP_PD_DP = 2


def _dedup_key(row: dict[str, Any]) -> tuple[float, float, float]:
    return (
        round(row["x0"], _DEDUP_X0_DP),
        round(row["jacobi"], _DEDUP_C_DP),
        round(row["period_days"], _DEDUP_PD_DP),
    )


def run_scan(
    *,
    cells: list[AsymCell] | None = None,
    n_workers: int = -1,
    out_path: Path = OUT_PATH,
    chunk_size: int = 4,
) -> dict[str, Any]:
    """Drive the parallel asymmetric-corrector scan and write JSONL.

    Returns a stats dict for the driver-level report.
    """
    if cells is None:
        cells = _build_grid()
    total = len(cells)
    out_path = Path(out_path)
    _log(f"scan_343 start: total cells={total} n_workers={n_workers} chunk_size={chunk_size}")
    _log(f"output: {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    cfg = ParallelSweepConfig(
        n_workers=n_workers,
        chunk_size=chunk_size,
        backend="loky",
        verbose=0,
        raise_on_first_error=False,
    )

    # Chunked processing -- parallel_sweep returns only after the WHOLE batch
    # finishes, so a single sweep over 5,832 cells is unrecoverable if killed.
    # We slice the grid into batches and stream each batch's results to a
    # partial JSONL after every batch returns. The final aggregation re-reads
    # the partial file, dedups and sorts, then writes the canonical output.
    stream_path = out_path.with_suffix(".stream.jsonl")
    stream_fh = stream_path.open("w")
    seen_stream: set[tuple[float, float, float]] = set()
    batches: list[tuple[list[AsymCell], int]] = []
    batch_size = max(int(total / 8), 1)  # ~8 batches per scan
    for batch_start in range(0, total, batch_size):
        batches.append((cells[batch_start : batch_start + batch_size], batch_start))

    all_results: list[Any] = [None] * total
    all_per_cell: list[float] = [0.0] * total
    n_ok_total = 0
    n_fail_total = 0
    notes_parts: list[str] = []
    for batch_idx, (batch_cells, batch_start) in enumerate(batches):
        batch_t0 = time.time()
        batch_result = parallel_sweep(batch_cells, process_cell, config=cfg)
        for j, value in enumerate(batch_result.results):
            all_results[batch_start + j] = value
            if value is not None:
                k = _dedup_key(value)
                if k not in seen_stream:
                    seen_stream.add(k)
                    stream_fh.write(json.dumps(value) + "\n")
                    stream_fh.flush()
        for j, dt in enumerate(batch_result.per_cell_elapsed_seconds):
            all_per_cell[batch_start + j] = dt
        n_ok_total += batch_result.n_succeeded
        n_fail_total += batch_result.n_failed
        if batch_result.notes:
            notes_parts.append(f"batch {batch_idx}: {batch_result.notes}")
        batch_elapsed = time.time() - batch_t0
        cumulative_elapsed = time.time() - t0
        cells_done = batch_start + len(batch_cells)
        rate = cells_done / max(cumulative_elapsed, 1e-9)
        eta_s = (total - cells_done) / max(rate, 1e-9)
        _log(
            f"batch {batch_idx + 1}/{len(batches)} done in {batch_elapsed:.0f}s "
            f"(cumulative {cumulative_elapsed:.0f}s, {cells_done}/{total} cells, "
            f"unique-so-far={len(seen_stream)}, eta {eta_s:.0f}s)"
        )
    stream_fh.close()
    elapsed = time.time() - t0

    # Synthesize a top-level result wrapper.
    from cyclerfinder.parallel import ParallelSweepResult

    result = ParallelSweepResult(
        results=tuple(all_results),
        n_cells=total,
        n_succeeded=n_ok_total,
        n_failed=n_fail_total,
        elapsed_seconds=elapsed,
        per_cell_elapsed_seconds=tuple(all_per_cell),
        notes=" | ".join(notes_parts),
    )

    # Aggregate per-cell results.
    stats: dict[str, Any] = {
        "cells_attempted": total,
        "cells_succeeded_call": result.n_succeeded,  # closure ran (may have returned None)
        "cells_failed_call": result.n_failed,  # closure raised
        "cells_converged": 0,  # closure returned a row (converged + closure-passing)
        "topology_matches": 0,
        "independent_closure_pass": 0,
        "literature_fresh": 0,
        "ml_low_fp": 0,
        "silver": 0,
        "novelty_claimable": 0,
        "known_family_collisions": 0,
        "per_k_attempted": {},
        "per_k_converged": {},
        "per_k_silver": {},
        "per_k_novel": {},
    }
    for cell in cells:
        key = f"{cell.k1},{cell.k2}"
        stats["per_k_attempted"][key] = stats["per_k_attempted"].get(key, 0) + 1

    seen: set[tuple[float, float, float]] = set()
    rows_to_write: list[dict[str, Any]] = []
    for row in result.results:
        if row is None:
            continue
        # Dedupe per-orbit at grid resolution.
        k = _dedup_key(row)
        if k in seen:
            continue
        seen.add(k)
        rows_to_write.append(row)

        stats["cells_converged"] += 1
        kkey = f"{row['k_target'][0]},{row['k_target'][1]}"
        stats["per_k_converged"][kkey] = stats["per_k_converged"].get(kkey, 0) + 1

        if row["topology_match"]:
            stats["topology_matches"] += 1
        if row["independent_closure_pass"]:
            stats["independent_closure_pass"] += 1
        if row["literature_fresh_offline"]:
            stats["literature_fresh"] += 1
        if row["ml_low_fp"]:
            stats["ml_low_fp"] += 1
        if row["known_em_family_collision"]:
            stats["known_family_collisions"] += 1

        silver = row["topology_match"] and row["independent_closure_pass"] and row["ml_low_fp"]
        if silver:
            stats["silver"] += 1
            stats["per_k_silver"][kkey] = stats["per_k_silver"].get(kkey, 0) + 1

        if row["novelty_claimable"]:
            stats["novelty_claimable"] += 1
            stats["per_k_novel"][kkey] = stats["per_k_novel"].get(kkey, 0) + 1

    # Write JSONL (sorted by prioritizer Tier 0 if available, else by k_target).
    def _rank_key(r: dict[str, Any]) -> tuple[float, str]:
        t0v = r.get("prioritizer_tier0_dv_kms")
        if t0v is None:
            t0v = float("inf")
        return (float(t0v), f"{r['k_target'][0]},{r['k_target'][1]}")

    rows_to_write.sort(key=_rank_key)
    with out_path.open("w") as fh:
        for row in rows_to_write:
            fh.write(json.dumps(row) + "\n")

    # Per-cell timing summary.
    if result.per_cell_elapsed_seconds:
        per_cell_total = sum(result.per_cell_elapsed_seconds)
        # Speedup estimate: serial cost = sum of per-cell costs, parallel wall
        # = elapsed. (per-cell costs are wall-clock-per-cell inside the worker.)
        speedup = per_cell_total / max(elapsed, 1e-9)
    else:
        per_cell_total = 0.0
        speedup = 0.0

    stats["elapsed_seconds"] = round(elapsed, 1)
    stats["per_cell_total_seconds"] = round(per_cell_total, 1)
    stats["speedup_vs_serial"] = round(speedup, 2)
    stats["unique_orbits_written"] = len(rows_to_write)

    _log(
        f"scan_343 done in {elapsed:.0f}s "
        f"(per-cell sum {per_cell_total:.0f}s -> {speedup:.1f}x speedup): "
        f"attempted={total} call_ok={result.n_succeeded} call_fail={result.n_failed} "
        f"unique={len(rows_to_write)} topo_match={stats['topology_matches']} "
        f"ind_pass={stats['independent_closure_pass']} silver={stats['silver']} "
        f"lit_fresh={stats['literature_fresh']} novel-claimable={stats['novelty_claimable']} "
        f"known_collisions={stats['known_family_collisions']}"
    )
    if result.notes:
        _log(f"sweep notes: {result.notes}")
    return stats


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Tiny grid for CI smoke (covers (1,1) + (2,1), 4 cells).",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=-1,
        help="Number of parallel workers (-1 = all cores).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=4,
        help="Cells per task batch (default 4; smaller = less worker startup).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUT_PATH,
        help=f"Output JSONL path (default: {OUT_PATH}).",
    )
    args = parser.parse_args(argv)

    if args.smoke:
        cells = _build_grid(
            k_targets=((1, 1), (2, 1)),
            c_grid=(3.14,),
            x0_grid=(-0.81, -0.80),
            xdot0_seeds=(-0.05,),
            ydot0_signs=(-1.0,),
            half_crossings=(6,),
        )
        run_scan(
            cells=cells,
            n_workers=args.n_workers,
            out_path=args.out,
            chunk_size=args.chunk_size,
        )
    else:
        run_scan(
            n_workers=args.n_workers,
            out_path=args.out,
            chunk_size=args.chunk_size,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
