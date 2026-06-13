"""Discovery-campaign daemon (#253, Track C) — the thin driver.

Generalizes ``scripts/search_campaign_daemon.py`` from a CR3BP corpus loop into a
DISCOVERY run: it drives the :mod:`cyclerfinder.search.discovery_campaign` engine
over a pluggable search target (the #254 repeated-moon multi-rev genome over
Jupiter, wired here), routing every outcome to the unchanged SILVER->gauntlet
artefacts and NEVER to the catalogue.

Like the corpus daemon this is a PLAIN, RESUMABLE, QUOTA-PROOF process:

* writes the progress checkpoint to the gitignored ``out/`` tree (never ``data/``
  or git, so it never collides with concurrent work),
* shards the deterministic candidate stream by ``--worker-id`` so N copies use N
  cores with no shared-file contention,
* resumes mid-stream on restart from the checkpoint.

SILVER survivors land in ``data/review_queue.jsonl`` (the gauntlet governs; a
human promotes — the daemon never does). A no-hit sweep appends one
method-versioned record to ``data/empty_regions.jsonl``. The catalogue is
read-only.

Launch one worker per core, e.g.::

    for w in 0 1 2 3; do
      uv run python scripts/discovery_campaign_daemon.py --worker-id "$w" --n-workers 4 &
    done

Each worker has its own checkpoint shard; ``data/review_queue.jsonl`` and
``data/empty_regions.jsonl`` are append-only and JSONL, so concurrent appends
interleave cleanly (one line per write).
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from cyclerfinder.search.discovery_campaign import (
    CampaignConfig,
    CampaignRouting,
    RepeatedMoonTarget,
    run_campaign,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--n-workers", type=int, default=1)
    ap.add_argument("--primary", default="Jupiter")
    ap.add_argument(
        "--seq-lengths",
        default="3,4",
        help="comma-separated repeated-moon sequence lengths to sweep",
    )
    ap.add_argument(
        "--max-rev", type=int, default=3, help="max Lambert n_rev per leg (grid 0..max)"
    )
    ap.add_argument("--phase-samples", type=int, default=24)
    ap.add_argument(
        "--gate-residual-kms",
        type=float,
        default=0.05,
        help="canonical residual a candidate must beat to route SILVER (km/s)",
    )
    ap.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="cap candidates this run (default: exhaust the enumeration)",
    )
    ap.add_argument(
        "--out-dir",
        default="out/discovery_campaign",
        help="gitignored checkpoint dir (out/ is gitignored)",
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    target = RepeatedMoonTarget(
        primary=args.primary,
        seq_lengths=tuple(int(x) for x in args.seq_lengths.split(",")),
        n_rev_grid=tuple(range(args.max_rev + 1)),
        n_phase_samples=args.phase_samples,
        git_sha=_git_sha(),
    )
    config = CampaignConfig(
        gate_residual_kms=args.gate_residual_kms,
        worker_id=args.worker_id,
        n_workers=args.n_workers,
        max_candidates=args.max_candidates,
    )
    # data/ registries (real SILVER + negative artefacts); checkpoint in out/.
    data_dir = Path(__file__).resolve().parent.parent / "data"
    routing = CampaignRouting(
        review_queue_path=data_dir / "review_queue.jsonl",
        empty_regions_path=data_dir / "empty_regions.jsonl",
        checkpoint_path=out_dir / f"checkpoint_{target.target_id}_w{args.worker_id}.txt",
    )

    print(
        f"[discovery] target={target.target_id} worker={args.worker_id}/{args.n_workers} "
        f"seq_lengths={target.seq_lengths} max_rev={args.max_rev} "
        f"gate={config.gate_residual_kms} km/s sha={target.method_capability().git_sha}",
        flush=True,
    )
    stats = run_campaign(target, config, routing)
    print(f"[discovery] DONE {stats.as_dict()}", flush=True)
    # Best-effort result line in out/ (gitignored).
    (out_dir / f"result_{target.target_id}_w{args.worker_id}.json").write_text(
        json.dumps(stats.as_dict(), indent=2)
    )


if __name__ == "__main__":
    main()
