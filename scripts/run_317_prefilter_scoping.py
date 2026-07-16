#!/usr/bin/env python3
"""#317 scoping: measure whether a learned sweep-impossibility pre-filter is warranted.

See ``cyclerfinder.ml.sweep_prefilter_scoping``'s module docstring for the
full argument. This script produces the REAL numbers that back it:

1. **CR3BP-corrector regime** (no existing pre-gate): assembles cheap
   pre-corrector features from this project's own ``#210`` outcome-log
   corpus (``out/outcome_log/search_campaign_w*.jsonl``), fits a logistic
   regression predicting convergence, and reports held-out accuracy/AUC plus
   a recall/skip table -- the decision-relevant "how much compute would we
   actually save at a safe recall level" question.
2. **Lambert/resonance-construction regime** (existing closed-form gates):
   aggregates the ``direction_summary``/``sequence_summary`` records already
   emitted by this project's ``enumerate_*_symmetric_closures.jsonl`` sweep
   family, reporting how much of each sweep's compute the cheap gate already
   rejects for free, and the hit rate of the pool that reaches the expensive
   gate.
3. A wall-clock cost benchmark: ``jacobi_constant`` (the cheapest feature)
   vs. ``correct_periodic`` (the actual corrector) -- how large the
   "expensive step" really is, so the recall/skip numbers can be read as an
   actual time saving, not just an abstract fraction.

No catalogue writeback, no literature-novelty check -- this is a scoping /
capability-negative investigation, not a discovery result.
"""

from __future__ import annotations

import glob
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cyclerfinder.core.cr3bp import cr3bp_system, jacobi_constant
from cyclerfinder.ml.sweep_prefilter_scoping import (
    CHEAP_FEATURE_NAMES,
    assemble_prefilter_dataset,
    auc_score,
    fit_logistic_prefilter,
    gate_efficiency_from_summary_records,
    recall_skip_table,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

SEED = 317
_OUT_DIR = Path("data/found/317_prefilter_scoping")

# Every enumerate_*_symmetric_closures.jsonl sweep on disk descended from
# scripts/scan_558_uranus_all_pairs_offset_sweep.py's gate machinery -- these
# already emit the direction_summary/sequence_summary records the gate-
# efficiency aggregation reads.
_SYMMETRIC_CLOSURE_SWEEPS = sorted(glob.glob("data/enumerate_*symmetric_closures.jsonl"))


def _load_jsonl(path: str) -> list[dict]:
    records = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def main() -> int:
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    # ---- Part 1: CR3BP-corrector regime (no existing pre-gate) ----
    log_paths = sorted(Path("out/outcome_log").glob("search_campaign_w*.jsonl"))
    if not log_paths:
        print("No out/outcome_log/search_campaign_w*.jsonl found -- cannot run Part 1.")
        return 1
    print(f"[{time.time() - t0:.0f}s] scanning {len(log_paths)} outcome-log files...")
    dataset = assemble_prefilter_dataset(log_paths)
    n = len(dataset)
    conv_rate = float(dataset.converged.mean())
    print(
        f"[{time.time() - t0:.0f}s] assembled {n} rows (scanned {dataset.n_scanned}), "
        f"raw converged rate {conv_rate:.1%}"
    )

    idx = rng.permutation(n)
    n_train = int(0.7 * n)
    train_idx, test_idx = idx[:n_train], idx[n_train:]
    x_train, y_train = dataset.features[train_idx], dataset.converged[train_idx]
    x_test, y_test = dataset.features[test_idx], dataset.converged[test_idx]

    model = fit_logistic_prefilter(x_train, y_train, seed=SEED)
    scores_test = model.decision_scores(x_test)
    proba_test = model.predict_proba(x_test)
    pred = proba_test >= 0.5
    accuracy = float((pred == y_test).mean())
    majority_baseline = float(max(y_test.mean(), 1.0 - y_test.mean()))
    auc = auc_score(scores_test, y_test)
    recall_table = recall_skip_table(scores_test, y_test)
    print(
        f"[{time.time() - t0:.0f}s] held-out: n={len(y_test)} accuracy={accuracy:.4f} "
        f"(majority baseline {majority_baseline:.4f}) AUC={auc:.4f}"
    )
    for row in recall_table:
        print(
            f"    target_recall={row['target_recall']:.3f} -> actual_recall="
            f"{row['actual_recall']:.4f} skip_fraction_of_negatives="
            f"{row['skip_fraction_of_negatives']:.4f}"
        )

    # ---- Cost benchmark: cheap feature vs. the actual corrector ----
    sys_em = cr3bp_system("Earth", "Moon")
    s0 = np.array(
        [-0.4145618480314011, 2.77e-23, 0.9075312043329506, -1.15e-12, 1.4076145460136695, 4.0e-13]
    )
    n_rep = 20000
    t_a = time.perf_counter()
    for _ in range(n_rep):
        jacobi_constant(s0, sys_em.mu)
    jacobi_us = (time.perf_counter() - t_a) / n_rep * 1e6

    n_rep_corrector = 100
    t_b = time.perf_counter()
    for _ in range(n_rep_corrector):
        correct_periodic(sys_em, s0, 3.12)
    corrector_ms = (time.perf_counter() - t_b) / n_rep_corrector * 1e3
    cost_ratio = corrector_ms * 1e3 / jacobi_us
    print(
        f"[{time.time() - t0:.0f}s] cost benchmark: jacobi_constant={jacobi_us:.2f}us/call, "
        f"correct_periodic={corrector_ms:.1f}ms/call, ratio={cost_ratio:.0f}x"
    )

    # ---- Part 2: Lambert/resonance-construction regime (existing gates) ----
    print(
        f"[{time.time() - t0:.0f}s] aggregating {len(_SYMMETRIC_CLOSURE_SWEEPS)} "
        f"enumerate_*_symmetric_closures.jsonl sweeps..."
    )
    per_sweep: dict[str, dict] = {}
    for path in _SYMMETRIC_CLOSURE_SWEEPS:
        records = _load_jsonl(path)
        summary = gate_efficiency_from_summary_records(records)
        per_sweep[Path(path).name] = {
            "n_evaluated": summary.n_evaluated,
            "n_infeasible": summary.n_infeasible,
            "n_subgate_reaches_expensive_gate": summary.n_subgate_reaches_expensive_gate,
            "n_all_gates_passed": summary.n_all_gates_passed,
            "frac_rejected_by_cheap_gate": summary.frac_rejected_by_cheap_gate,
            "frac_reaching_expensive_gate": summary.frac_reaching_expensive_gate,
            "hit_rate_of_expensive_gate_pool": summary.hit_rate_of_expensive_gate_pool,
        }
        print(
            f"    {Path(path).name}: n_evaluated={summary.n_evaluated} "
            f"cheap_reject={summary.frac_rejected_by_cheap_gate:.1%} "
            f"reaches_expensive={summary.frac_reaching_expensive_gate:.1%} "
            f"hit_rate_of_expensive_pool={summary.hit_rate_of_expensive_gate_pool:.1%}"
        )

    overall = gate_efficiency_from_summary_records(
        rec for path in _SYMMETRIC_CLOSURE_SWEEPS for rec in _load_jsonl(path)
    )
    # The blended "overall" number is dominated by two huge zero-hit sweeps
    # (#600's 806k-candidate 3-moon enumeration, #607's mass-limited
    # small-body systems) -- split by "did this sweep find ANY closure at
    # all" so the headline doesn't bury the fruitful sweeps' healthy hit
    # rate under the empty ones' bulk.
    fruitful_paths = [
        p for p in _SYMMETRIC_CLOSURE_SWEEPS if per_sweep[Path(p).name]["n_all_gates_passed"] > 0
    ]
    empty_paths = [
        p for p in _SYMMETRIC_CLOSURE_SWEEPS if per_sweep[Path(p).name]["n_all_gates_passed"] == 0
    ]
    fruitful = gate_efficiency_from_summary_records(
        rec for path in fruitful_paths for rec in _load_jsonl(path)
    )
    empty = gate_efficiency_from_summary_records(
        rec for path in empty_paths for rec in _load_jsonl(path)
    )

    elapsed = time.time() - t0
    summary_out = {
        "task": "#317 sweep-impossible-region pre-filter scoping",
        "seed": SEED,
        "cr3bp_corrector_regime": {
            "description": "blind search_campaign_daemon Phase B (no existing pre-gate)",
            "feature_names": list(CHEAP_FEATURE_NAMES),
            "n_scanned": dataset.n_scanned,
            "n_rows": n,
            "raw_converged_rate": conv_rate,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "held_out_accuracy": accuracy,
            "held_out_majority_baseline_accuracy": majority_baseline,
            "held_out_auc": auc,
            "recall_skip_table": recall_table,
        },
        "cost_benchmark": {
            "jacobi_constant_us_per_call": jacobi_us,
            "correct_periodic_ms_per_call": corrector_ms,
            "cost_ratio_corrector_over_jacobi": cost_ratio,
        },
        "lambert_resonance_regime": {
            "description": (
                "enumerate_*_symmetric_closures.jsonl sweeps -- existing closed-form "
                "Lambert-existence + #324 bend gate already run before the expensive "
                "bend+DOP853 gate_candidate step"
            ),
            "per_sweep": per_sweep,
            "overall": {
                "n_evaluated": overall.n_evaluated,
                "n_infeasible": overall.n_infeasible,
                "n_subgate_reaches_expensive_gate": overall.n_subgate_reaches_expensive_gate,
                "n_all_gates_passed": overall.n_all_gates_passed,
                "frac_rejected_by_cheap_gate": overall.frac_rejected_by_cheap_gate,
                "frac_reaching_expensive_gate": overall.frac_reaching_expensive_gate,
                "hit_rate_of_expensive_gate_pool": overall.hit_rate_of_expensive_gate_pool,
            },
            "fruitful_sweeps_only": {
                "sweeps": [Path(p).name for p in fruitful_paths],
                "n_evaluated": fruitful.n_evaluated,
                "frac_rejected_by_cheap_gate": fruitful.frac_rejected_by_cheap_gate,
                "frac_reaching_expensive_gate": fruitful.frac_reaching_expensive_gate,
                "hit_rate_of_expensive_gate_pool": fruitful.hit_rate_of_expensive_gate_pool,
            },
            "confirmed_empty_sweeps_only": {
                "sweeps": [Path(p).name for p in empty_paths],
                "n_evaluated": empty.n_evaluated,
                "frac_rejected_by_cheap_gate": empty.frac_rejected_by_cheap_gate,
                "frac_reaching_expensive_gate": empty.frac_reaching_expensive_gate,
                "hit_rate_of_expensive_gate_pool": empty.hit_rate_of_expensive_gate_pool,
            },
        },
        "elapsed_s": elapsed,
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary_out, indent=2, sort_keys=False) + "\n")
    print(f"[{elapsed:.0f}s] wrote {summary_path}")

    print()
    print("HONEST HEADLINE:")
    print(
        f"  CR3BP-corrector regime (no existing gate): AUC={auc:.3f}, held-out accuracy "
        f"{accuracy:.1%} vs {majority_baseline:.1%} majority baseline."
    )
    print(
        "  At a SAFE recall level for a program that treats real hits as precious "
        "(this project's own memory: 'Novel hits are RARE'):"
    )
    for row in recall_table:
        print(
            f"    recall={row['target_recall']:.3f} -> would skip only "
            f"{row['skip_fraction_of_negatives']:.1%} of doomed corrector calls."
        )
    print(
        f"  The corrector itself costs ~{cost_ratio:.0f}x a cheap feature computation, so "
        "ANY skip is pure savings -- but the achievable skip fraction at a safe recall "
        "is negligible. VERDICT: not worth building/maintaining a learned pre-filter here."
    )
    print()
    print(
        f"  Lambert/resonance regime, ALL {len(_SYMMETRIC_CLOSURE_SWEEPS)} sweeps blended: "
        f"cheap gate alone rejects {overall.frac_rejected_by_cheap_gate:.1%} of all evaluated "
        f"candidates for free; only {overall.frac_reaching_expensive_gate:.1%} reaches the "
        f"expensive gate, blended hit rate {overall.hit_rate_of_expensive_gate_pool:.1%} "
        "(dominated by #600/#607's huge confirmed-EMPTY sweeps -- see split below)."
    )
    print(
        f"    FRUITFUL sweeps only ({', '.join(Path(p).name for p in fruitful_paths)}): "
        f"cheap gate rejects {fruitful.frac_rejected_by_cheap_gate:.1%}, "
        f"{fruitful.frac_reaching_expensive_gate:.1%} reaches the expensive gate, "
        f"hit rate {fruitful.hit_rate_of_expensive_gate_pool:.1%} of THAT pool -- "
        "exactly where these sweeps' real discoveries live, not further junk to prune."
    )
    print(
        f"    CONFIRMED-EMPTY sweeps only ({len(empty_paths)} sweeps, "
        f"{empty.n_evaluated} candidates): 0% hit rate is a per-SYSTEM physical limit "
        "(e.g. #607's mass-limited moons), already diagnosed via the SAME closed-form "
        "bend gate -- not a pattern a classifier would need to learn separately."
    )
    print(
        "  VERDICT: no room for a learned pre-filter to add value here; the existing "
        "closed-form Lambert-existence + #324 bend gates already do this job exactly "
        "and for free, and pruning the fruitful pool further risks discarding the "
        "sweeps' real hits."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
