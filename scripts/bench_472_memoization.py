"""#472 cached-vs-uncached benchmark for the memoized per-leg cost functions.

INFORMATIONAL ONLY — this is NOT a CI gate (no timing assertion lives in the
test suite; timing is machine-dependent and flaky). It runs a small
#468-style nested moon-tour sweep twice (cache warm vs every cache cleared per
call) and reports the wall-clock speedup factor, plus a parity assertion that
the two passes produce BIT-IDENTICAL results (the real point: the cache changes
no value).

Per the incremental-progress-reports requirement (feedback_incremental_progress_
reports, refined 2026-06-26): the nested sweep writes one flushed JSONL line per
sweep unit to ``out/bench_472_runlog.jsonl`` — ``{item_id, sub_step, elapsed_s,
ts}`` — so this is never a black box that only prints at the end.

Run: ``uv run python scripts/bench_472_memoization.py``
"""

from __future__ import annotations

import json
import time
from datetime import UTC, datetime
from pathlib import Path

from cyclerfinder.core import flyby
from cyclerfinder.search import tisserand, vilm
from cyclerfinder.search.correct import _max_bend_deg_nominal
from cyclerfinder.search.moon_prune import moon_leg_admissible

REPO = Path(__file__).resolve().parents[1]
RUNLOG = REPO / "out" / "bench_472_runlog.jsonl"

# A small #468-style moon-tour skeleton sweep: adjacent Jovian/Saturnian moon
# pairs, swept over a V∞ grid and a budget grid (mirrors the repeated-evaluation
# pattern of the campaign's ToF/phase sweeps).
_LEGS = [
    ("Jupiter", "Io", "Europa"),
    ("Jupiter", "Europa", "Ganymede"),
    ("Jupiter", "Ganymede", "Callisto"),
    ("Saturn", "Titan", "Rhea"),
    ("Saturn", "Rhea", "Dione"),
    ("Saturn", "Dione", "Tethys"),
]
_VINF_GRID = [round(2.0 + 0.1 * k, 4) for k in range(20)]  # 20 discrete V∞
_BUDGET_GRID = [1.0, 2.0, 3.0, 4.0]
_REPEATS = 8  # the campaign re-hits the same (leg, V∞, budget) across phasings


_CACHED_FNS = (
    vilm._vc_adim,
    vilm._vbar_vinf_adim,
    vilm.min_vinf_for_vilm,
    vilm._v_m,
    vilm._leverage_dv_kms,
    vilm._vilm_dv_min_pair,
    vilm.vilm_dv_floor,
    vilm.europa_endgame_dv,
    tisserand._a_p_km,
    tisserand.vinf_to_tisserand,
    tisserand.tisserand_to_vinf,
    tisserand.linkable,
    flyby.max_bend,
    flyby.dv_from_turn_deficit,
    flyby.dv_powered_flyby_periapsis,
    _max_bend_deg_nominal,
)


def _clear_all() -> None:
    for fn in _CACHED_FNS:
        fn.cache_clear()


def _sweep(*, bust_cache: bool, runlog) -> list[tuple]:
    """Run the nested sweep once. When ``bust_cache`` clears every cache before
    each leg-eval (simulating the uncached cost), otherwise leave caches warm."""
    results: list[tuple] = []
    t0 = time.perf_counter()
    for primary, a, b in _LEGS:
        item_id = f"{primary}:{a}->{b}"
        for vinf in _VINF_GRID:
            for budget in _BUDGET_GRID:
                for _ in range(_REPEATS):
                    if bust_cache:
                        _clear_all()
                    ok, _reason = moon_leg_admissible(
                        a, b, vinf_kms=vinf, budget_kms=budget, primary=primary
                    )
                    results.append((primary, a, b, vinf, budget, ok))
            elapsed = time.perf_counter() - t0
            runlog.write(
                json.dumps(
                    {
                        "item_id": item_id,
                        "sub_step": f"vinf={vinf}",
                        "mode": "uncached" if bust_cache else "cached",
                        "elapsed_s": round(elapsed, 4),
                        "ts": datetime.now(UTC).isoformat(),
                    }
                )
                + "\n"
            )
            runlog.flush()
    return results


def main() -> None:
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    with RUNLOG.open("w", encoding="utf-8") as runlog:
        # Uncached pass first (clears cache before every call).
        _clear_all()
        t0 = time.perf_counter()
        res_uncached = _sweep(bust_cache=True, runlog=runlog)
        t_uncached = time.perf_counter() - t0

        # Cached pass (caches stay warm across the whole sweep).
        _clear_all()
        t0 = time.perf_counter()
        res_cached = _sweep(bust_cache=False, runlog=runlog)
        t_cached = time.perf_counter() - t0

    # PARITY: the two passes must be bit-identical (the cache changes nothing).
    assert res_cached == res_uncached, "PARITY FAILURE: cached != uncached results"

    speedup = t_uncached / t_cached if t_cached > 0 else float("inf")
    n_evals = len(res_cached)
    print(
        json.dumps(
            {
                "task": "#472 memoization benchmark (informational)",
                "n_leg_evals": n_evals,
                "t_uncached_s": round(t_uncached, 4),
                "t_cached_s": round(t_cached, 4),
                "speedup_factor": round(speedup, 2),
                "parity_bit_identical": True,
                "runlog": str(RUNLOG.relative_to(REPO)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
