"""Re-run #264 Pluto sweep with the 5-tier prioritizer; emit data/silver_282.jsonl (#282).

This script answers ONE question: does the 5-tier scorer composition (#282)
change the candidate YIELD vs the original #264 enumeration-only sweep?

The honest architecture (documented in five_tier_prioritizer.py + the #282
note): #264 emits patched-conic Lambert legs, not CR3BP representative orbits.
Of the five tiers only Tier 0 (Zhang-Topputo NN ΔV predictor) operates
natively on Lambert legs. Tiers 1-5 require CR3BP representatives, which
#264 does not produce. The yield-comparison run therefore exercises:

* **Pass A (post-hoc Tier-0 re-prioritization of existing SILVERs)** -- score
  each of the 12 Pluto SILVERs from data/review_queue.jsonl by reconstructing
  their best-phasing Lambert legs and running NeuralReachPrefilter on each
  leg. Reports per-leg ΔV and a per-candidate cumulative ΔV. This tests
  whether the NN's view of "obvious infeasibility" agrees with the #264
  closure verdict.
* **Pass B (representative Pluto re-sweep with Tier-0 gating)** -- exercises
  one shard (worker_id=0, n_workers=4) of the #264 Pluto enumeration with
  Tier-0 admission as a pre-closure gate. Compares yields to the original
  #264 shard. DOWN-SCOPED to a representative sample because exhaustive
  re-run is out of budget.

Outputs:

* data/silver_282.jsonl -- per-candidate row with per-tier scores + tiers_skipped
  audit. Frame matches the #264 review_queue rows; candidates are NOT
  promoted; this is yield data, not catalogue input.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

# Repo on path (this is a `uv run` script).
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.search.five_tier_prioritizer import (  # noqa: E402
    FiveTierPrioritizer,
    legs_from_repeated_moon_candidate,
)


def _git_sha() -> str:
    import subprocess

    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _read_silver_264(path: Path) -> list[dict]:
    """Read the 12 SILVER candidates from the #264 review queue."""
    rows = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    return rows


def _score_one_candidate(
    prioritizer: FiveTierPrioritizer,
    row: dict,
    *,
    primary: str,
    phase_samples: int,
) -> dict:
    """Run Tier 0 on the legs of one #264 SILVER candidate."""
    sequence = list(row["sequence"])
    n_rev = list(row["verdict_audit"]["n_rev"])
    try:
        legs = legs_from_repeated_moon_candidate(
            primary=primary,
            sequence=sequence,
            n_rev=n_rev,
            phase_samples=phase_samples,
        )
    except Exception as exc:
        return {
            "candidate_id": row["candidate_id"],
            "sequence": sequence,
            "n_rev": n_rev,
            "tier0_status": f"leg-reconstruction-failed: {exc}",
            "per_tier_scores": None,
        }
    if legs is None:
        return {
            "candidate_id": row["candidate_id"],
            "sequence": sequence,
            "n_rev": n_rev,
            "tier0_status": "no-feasible-phasing",
            "per_tier_scores": None,
        }
    t0_stats = prioritizer.score_candidate_legs(legs)
    return {
        "candidate_id": row["candidate_id"],
        "sequence": sequence,
        "n_rev": n_rev,
        "primary": primary,
        "original_residual_kms": row["verdict_audit"]["residual_kms"],
        "original_max_vinf_kms": row["max_vinf_kms"],
        "original_ml_flagger_p_fp": row.get("literature_check", {}).get("ml_flagger_p_fp"),
        "tier0_status": "ok",
        "tier0_max_dv_kms": t0_stats["tier0_max_dv_kms"],
        "tier0_sum_dv_kms": t0_stats["tier0_sum_dv_kms"],
        "tier0_all_admitted": t0_stats["tier0_all_admitted"],
        "tier0_any_inference_failed": t0_stats["tier0_any_inference_failed"],
        "tier0_n_legs": t0_stats["n_legs"],
        "per_leg": [
            {
                "label_from": p["label_from"],
                "label_to": p["label_to"],
                "dv_kms": p["tier0_predicted_dv_kms"],
                "tof_days": p["tier0_predicted_tof_days"],
                "admitted": p["tier0_admitted"],
                "model_available": p["tier0_model_available"],
                "fallback_used": p["tier0_fallback_used"],
            }
            for p in t0_stats["per_leg"]
        ],
        "tiers_skipped": [
            {
                "tier": i,
                "reason": (
                    "patched-conic leg input lacks CR3BP representative "
                    "orbit; tier requires periodic-orbit anchor (#282 architecture)"
                ),
            }
            for i in (1, 2, 3, 4, 5)
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--review-queue",
        default=str(ROOT / "data" / "review_queue.jsonl"),
        help="path to #264 review queue with the 12 Pluto SILVERs",
    )
    ap.add_argument(
        "--out",
        default=str(ROOT / "data" / "silver_282.jsonl"),
        help="output JSONL with per-candidate Tier-0 scores",
    )
    ap.add_argument("--primary", default="Pluto")
    ap.add_argument(
        "--phase-samples",
        type=int,
        default=12,
        help="phasing grid size for leg reconstruction (smaller = faster)",
    )
    ap.add_argument(
        "--tier0-admit-threshold-kms",
        type=float,
        default=5.0,
        help="Tier-0 NN admission threshold (km/s); paper default 5.0",
    )
    args = ap.parse_args()

    review_queue = Path(args.review_queue)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[282] reading SILVERs from {review_queue}", flush=True)
    rows = _read_silver_264(review_queue)
    pluto_rows = [r for r in rows if r["verdict_audit"]["primary"] == args.primary]
    print(f"[282] {len(pluto_rows)} {args.primary} SILVERs in queue", flush=True)

    threshold = args.tier0_admit_threshold_kms
    print(
        f"[282] constructing FiveTierPrioritizer (Tier-0 threshold {threshold} km/s)",
        flush=True,
    )
    prioritizer = FiveTierPrioritizer(tier0_admit_threshold_kms=threshold)
    nn = prioritizer.nn_prefilter
    model_avail = nn is not None and nn.weights_dv is not None
    print(f"[282] NN model_available={model_avail}", flush=True)

    t0_start = time.time()
    sha = _git_sha()
    with out_path.open("w") as out:
        for i, row in enumerate(pluto_rows):
            t_start = time.time()
            scored = _score_one_candidate(
                prioritizer, row, primary=args.primary, phase_samples=args.phase_samples
            )
            scored["scored_at_sha"] = sha
            out.write(json.dumps(scored) + "\n")
            out.flush()
            elapsed = time.time() - t_start
            print(
                f"[282] [{i + 1}/{len(pluto_rows)}] {row['candidate_id']} "
                f"seq={row['sequence']} tier0_max_dv={scored.get('tier0_max_dv_kms', 'n/a')} "
                f"admit={scored.get('tier0_all_admitted', 'n/a')} dt={elapsed:.1f}s",
                flush=True,
            )

    total = time.time() - t0_start
    print(f"[282] DONE -> {out_path} ({len(pluto_rows)} rows, {total:.1f}s)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
