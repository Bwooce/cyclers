"""#474 cross-worker cache benchmark for the parallel moon-tour campaign.

INFORMATIONAL ONLY — NOT a CI gate (timing is machine-dependent and flaky; no
timing assertion lives in the test suite). It runs a #468-style nested moon-tour
sweep under the parallel substrate THREE ways and reports the wall-clock of each
plus the parallel speedup the #474 pre-warm buys:

  A. loky, no prewarm        — the STATUS QUO: every fresh-interpreter worker
                               rebuilds the #472 cache cold (speedup evaporates).
  B. multiprocessing, no prewarm — forked workers, but the parent cache is cold
                               at fork, so each worker still builds its own copy.
  C. multiprocessing, prewarm — the #474 fix: parent warms the #472 caches, the
                               forked workers inherit them copy-on-write.

The real point is the in-script PARITY assertion: all three passes (and a direct
serial reference) produce BIT-IDENTICAL results — the cache / backend / pre-warm
change NO value. The speedup is reported as a number, not gated.

Per the incremental-progress-reports requirement: one flushed JSONL line per
sweep unit is written to ``out/bench_474_runlog.jsonl`` —
``{mode, item_id, sub_step, elapsed_s, ts}`` — so this is never a black box.

Run: ``uv run python scripts/bench_474_cross_worker_cache.py``
"""

from __future__ import annotations

import functools
import json
import time
from datetime import UTC, datetime
from pathlib import Path

from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
from cyclerfinder.search.cache_warm import DEFAULT_LEGS, warm_moon_leg_caches
from cyclerfinder.search.moon_prune import moon_leg_admissible

REPO = Path(__file__).resolve().parents[1]
RUNLOG = REPO / "out" / "bench_474_runlog.jsonl"

# A #468-style moon-tour skeleton. Larger V∞ grid than the test so the per-cell
# cost is non-trivial and the cold-cache penalty is visible at the wall clock.
_LEGS = DEFAULT_LEGS
_VINF_GRID = [round(2.0 + 0.05 * k, 4) for k in range(40)]  # 40 discrete V∞
_BUDGET_GRID = [1.0, 2.0, 3.0, 4.0]

_CELLS = [
    (primary, a, b, vinf, budget)
    for (primary, a, b) in _LEGS
    for vinf in _VINF_GRID
    for budget in _BUDGET_GRID
]


def _admit_cell(cell: tuple[str, str, str, float, float]) -> tuple:
    primary, a, b, vinf, budget = cell
    ok, reason = moon_leg_admissible(a, b, vinf_kms=vinf, budget_kms=budget, primary=primary)
    return (primary, a, b, vinf, budget, ok, reason)


def _serial_reference() -> list[tuple]:
    return [_admit_cell(c) for c in _CELLS]


def _run(mode: str, cfg: ParallelSweepConfig, runlog) -> tuple[list[tuple], float]:
    t0 = time.perf_counter()
    result = parallel_sweep(_CELLS, _admit_cell, config=cfg)
    elapsed = time.perf_counter() - t0
    runlog.write(
        json.dumps(
            {
                "mode": mode,
                "item_id": "moon_tour_sweep",
                "sub_step": "complete",
                "n_cells": result.n_cells,
                "n_failed": result.n_failed,
                "elapsed_s": round(elapsed, 4),
                "ts": datetime.now(UTC).isoformat(),
            }
        )
        + "\n"
    )
    runlog.flush()
    if result.n_failed:
        raise RuntimeError(f"{mode}: {result.n_failed} cells failed: {result.notes}")
    return list(result.results), elapsed


def main() -> None:
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)

    prewarm = functools.partial(
        warm_moon_leg_caches,
        legs=_LEGS,
        vinf_grid=_VINF_GRID,
        budget_grid=_BUDGET_GRID,
    )

    with RUNLOG.open("w", encoding="utf-8") as runlog:
        reference = _serial_reference()

        res_loky, t_loky = _run(
            "loky_no_prewarm",
            ParallelSweepConfig(n_workers=4, backend="loky"),
            runlog,
        )
        res_mp_cold, t_mp_cold = _run(
            "multiprocessing_no_prewarm",
            ParallelSweepConfig(n_workers=4, backend="multiprocessing"),
            runlog,
        )
        res_mp_warm, t_mp_warm = _run(
            "multiprocessing_prewarm",
            ParallelSweepConfig(n_workers=4, backend="multiprocessing", prewarm=prewarm),
            runlog,
        )

    # PARITY: every pass is bit-identical to the direct serial reference.
    assert res_loky == reference, "PARITY FAILURE: loky != serial reference"
    assert res_mp_cold == reference, "PARITY FAILURE: mp(cold) != serial reference"
    assert res_mp_warm == reference, "PARITY FAILURE: mp(prewarm) != serial reference"

    speedup_vs_loky = t_loky / t_mp_warm if t_mp_warm > 0 else float("inf")
    speedup_vs_mp_cold = t_mp_cold / t_mp_warm if t_mp_warm > 0 else float("inf")

    print(
        json.dumps(
            {
                "task": "#474 cross-worker cache benchmark (informational)",
                "n_cells": len(_CELLS),
                "t_loky_no_prewarm_s": round(t_loky, 4),
                "t_multiprocessing_no_prewarm_s": round(t_mp_cold, 4),
                "t_multiprocessing_prewarm_s": round(t_mp_warm, 4),
                "speedup_prewarm_vs_loky": round(speedup_vs_loky, 2),
                "speedup_prewarm_vs_mp_cold": round(speedup_vs_mp_cold, 2),
                "parity_bit_identical": True,
                "runlog": str(RUNLOG.relative_to(REPO)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
