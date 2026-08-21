"""#859: Resonant Atlas pilot, Stage A driver — family recovery + eigenvalue
survey, checkpointed via `campaign_runner` (`#788`/`#800`).

Builds the full Stage A cell grid
(:func:`cyclerfinder.search.resonant_atlas_stage_a.build_stage_a_cells`: 4
systems x every coprime ``p:q`` with ``p,q<=8``) and runs it through
:func:`cyclerfinder.search.campaign_runner.run_grid_campaign`, which
checkpoints every cell durably to ``results.jsonl`` the instant its batch
returns (kill-safe, resumable — re-invoking this script with the SAME
arguments picks up exactly where it left off).

**THIS SCRIPT IS THE HARNESS, NOT THE REAL RUN.** `#859`'s own registration
budgets the full Stage A sweep at ~50-100 CPU-hours / ~2-3 days wall at 4
workers — far beyond a single foreground dispatched-agent session. This
script was built and SMOKE-TESTED on small slices only (a handful of cells
per system, `--systems`/`--max-pq`/`--max-batches` all support narrowing for
exactly that purpose); the full run is left for the coordinating session to
launch and monitor as its own long-lived background process. See
``docs/notes/2026-08-21-859-resonant-atlas-pilot-harness.md`` for the
smoke-test's own measured results and a refined per-cell/per-system time
estimate (materially different from `#858`'s own inferred estimate — read
that note before dispatching the real run, not just this docstring).

Usage (small slice, matching this task's own smoke test):

    uv run python scripts/run_859_resonant_atlas_stage_a.py \\
        --systems uranus-oberon --max-pq 3 --n-c-steps 3 --n-workers 2

Full run (coordinating session only, NOT this task):

    uv run python scripts/run_859_resonant_atlas_stage_a.py --n-workers 4 \\
        --pause-seconds-per-batch 30 --thermal-backoff-seconds 120

Output: ``data/found/859_resonant_atlas_stage_a/results.jsonl`` (the
checkpoint itself, one JSON line per cell, per `campaign_runner`'s own
"results.jsonl is its own checkpoint" design) plus, on a full run,
``data/found/859_resonant_atlas_stage_a/census_report.json`` (the Stage A
deliverable: per-system in-band cell counts, written by ``--report``).

Foreground only — this project's own standing lesson
(`feedback_subagent_background_is_fatal`) is that a dispatched agent must
NEVER self-background a long-running process. Chunk long invocations via
``--max-batches``, or hand the un-narrowed command to a session that CAN
launch and monitor a genuine multi-day background process.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.search.campaign_runner import (
    CampaignRunnerConfig,
    CampaignRunnerRouting,
    run_grid_campaign,
)
from cyclerfinder.search.resonant_atlas_stage_a import (
    DEFAULT_D_JACOBI,
    DEFAULT_MAX_PQ,
    DEFAULT_N_C_STEPS,
    DEFAULT_X0_SIGN,
    STAGE_A_SYSTEMS,
    build_stage_a_cells,
    stage_a_worker,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = _REPO_ROOT / "data" / "found" / "859_resonant_atlas_stage_a"


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--systems",
        nargs="*",
        default=None,
        help=(
            "subset of system_key values to run (default: all 4 pilot systems). "
            f"Choices: {', '.join(s.system_key for s in STAGE_A_SYSTEMS)}"
        ),
    )
    ap.add_argument(
        "--max-pq",
        type=int,
        default=DEFAULT_MAX_PQ,
        help="coprime p:q ceiling (default: 8, per #859's own registration)",
    )
    ap.add_argument(
        "--n-c-steps",
        type=int,
        default=DEFAULT_N_C_STEPS,
        help="continuation steps per cell after the seed member (default: 9 -> 10 members)",
    )
    ap.add_argument(
        "--d-jacobi",
        type=float,
        default=DEFAULT_D_JACOBI,
        help="Jacobi-constant step size per continuation step (default: 5e-4)",
    )
    ap.add_argument("--x0-sign", type=int, default=DEFAULT_X0_SIGN, choices=(-1,))
    ap.add_argument("--n-workers", type=int, default=4)
    ap.add_argument("--checkpoint-batch-size", type=int, default=8)
    ap.add_argument(
        "--max-batches",
        type=int,
        default=None,
        help="cap batches run THIS invocation (default: exhaust pending cells)",
    )
    ap.add_argument("--pause-seconds-per-batch", type=float, default=0.0)
    ap.add_argument("--thermal-backoff-seconds", type=float, default=0.0)
    ap.add_argument("--timeout-seconds-per-cell", type=float, default=600.0)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=OUT_DIR,
        help="where results.jsonl (and --report's census_report.json) live",
    )
    ap.add_argument(
        "--report",
        action="store_true",
        help="after the run, write census_report.json (the Stage A in-band cell census)",
    )
    return ap.parse_args()


def _systems_for(keys: list[str] | None) -> tuple[Any, ...]:
    if keys is None:
        return STAGE_A_SYSTEMS
    by_key = {s.system_key: s for s in STAGE_A_SYSTEMS}
    unknown = [k for k in keys if k not in by_key]
    if unknown:
        raise SystemExit(f"unknown system_key(s): {unknown}; choices: {sorted(by_key)}")
    return tuple(by_key[k] for k in keys)


def _write_census_report(results_path: Path, out_path: Path, cells: list[dict[str, Any]]) -> None:
    """Stage A's own deliverable: per-system in-band (p:q, C) cell census.

    Reads ``results.jsonl`` directly (never trusts an in-memory running
    total) so the report is correct even if this invocation only ran a
    subset of batches — it reports on whatever is durably checkpointed so
    far, honestly labelled with how many of the grid's cells that covers.
    """
    by_index: dict[int, dict[str, Any]] = {}
    if results_path.exists():
        for line in results_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            by_index[int(row["index"])] = row

    per_system: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "n_cells_total": 0,
            "n_cells_evaluated": 0,
            "n_cells_in_band": 0,
            "n_cells_seed_no_converge": 0,
            "n_cells_seed_domain_invalid": 0,
            "n_cells_error": 0,
            "in_band_cells": [],
        }
    )
    for idx, cell in enumerate(cells):
        sys_key = cell["system_key"]
        stats = per_system[sys_key]
        stats["n_cells_total"] += 1
        row = by_index.get(idx)
        if row is None:
            continue
        stats["n_cells_evaluated"] += 1
        status = row["status"]
        payload = row.get("payload", {})
        if status == "error":
            stats["n_cells_error"] += 1
            continue
        reason = payload.get("reason")
        if reason == "seed_did_not_converge":
            stats["n_cells_seed_no_converge"] += 1
        elif reason == "seed_domain_invalid":
            stats["n_cells_seed_domain_invalid"] += 1
        n_in_band = int(payload.get("n_in_band", 0) or 0)
        if n_in_band > 0:
            stats["n_cells_in_band"] += 1
            stats["in_band_cells"].append(
                {
                    "p": cell["p"],
                    "q": cell["q"],
                    "n_in_band_members": n_in_band,
                    "max_abs_lambda": payload.get("max_abs_lambda"),
                    "min_abs_lambda": payload.get("min_abs_lambda"),
                }
            )

    report = {
        "task": "#859",
        "generated_utc": datetime.now(UTC).isoformat(),
        "in_band_low": 50.0,
        "in_band_high": 2500.0,
        "per_system": dict(per_system),
    }
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"census report written to {out_path}", flush=True)
    for sys_key, stats in per_system.items():
        print(
            f"  {sys_key}: {stats['n_cells_evaluated']}/{stats['n_cells_total']} evaluated, "
            f"{stats['n_cells_in_band']} in-band, "
            f"{stats['n_cells_seed_domain_invalid']} seed-domain-invalid, "
            f"{stats['n_cells_seed_no_converge']} seed-no-converge, "
            f"{stats['n_cells_error']} errors",
            flush=True,
        )


def main() -> None:
    args = _parse_args()
    systems = _systems_for(args.systems)
    cells = build_stage_a_cells(
        systems,
        max_pq=args.max_pq,
        n_c_steps=args.n_c_steps,
        d_jacobi=args.d_jacobi,
        x0_sign=args.x0_sign,
    )
    print(
        f"[{datetime.now(UTC).isoformat(timespec='seconds')}] Stage A grid: "
        f"{len(cells)} cells across {len(systems)} system(s) "
        f"({', '.join(s.system_key for s in systems)}), max_pq={args.max_pq}, "
        f"n_c_steps={args.n_c_steps}, d_jacobi={args.d_jacobi}",
        flush=True,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    routing = CampaignRunnerRouting(results_path=args.out_dir / "results.jsonl")
    config = CampaignRunnerConfig(
        n_workers=args.n_workers,
        checkpoint_batch_size=args.checkpoint_batch_size,
        timeout_seconds_per_cell=args.timeout_seconds_per_cell,
        max_batches=args.max_batches,
        pause_seconds_per_batch=args.pause_seconds_per_batch,
        thermal_backoff_seconds=args.thermal_backoff_seconds,
    )

    t0 = datetime.now(UTC)
    stats = run_grid_campaign(cells, stage_a_worker, routing=routing, config=config)
    dt_s = (datetime.now(UTC) - t0).total_seconds()

    print(
        f"[{datetime.now(UTC).isoformat(timespec='seconds')}] this invocation: "
        f"{stats.batches_run_this_invocation} batches in {dt_s:.1f}s "
        f"({dt_s / max(stats.batches_run_this_invocation, 1):.1f}s/batch of "
        f"{args.checkpoint_batch_size} cells); cumulative: "
        f"{stats.evaluated_total}/{stats.total_cells} evaluated, "
        f"hits={stats.hits} misses={stats.misses} errors={stats.errors}",
        flush=True,
    )
    if args.report:
        _write_census_report(routing.results_path, args.out_dir / "census_report.json", cells)


if __name__ == "__main__":
    main()
