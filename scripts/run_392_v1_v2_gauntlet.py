"""#392 V1 and V2 Gauntlet for Floquet low-amplitude cycler candidates.

Reads the branches from the low amplitude sweep output, filters for those that
pass the Hill screen, and runs the V1 and V2 gates on the highest Jacobi
constant branch (lowest amplitude).

Usage:
    uv run python scripts/run_392_v1_v2_gauntlet.py
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from cyclerfinder.data.validation.v2_3d import V2_DRIFT_FLOOR_KMS, run_v2_3d
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    correct_general_periodic_3d,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

V1_CORRECTOR_RESIDUAL_FLOOR = 1.0e-10
V1_INDEPENDENT_CLOSURE_FLOOR = 1.0e-6
PHASE_LABEL = "392_gauntlet"
N_CYCLES_LIST = (3, 5, 10)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/scan_392_floquet_low_amp.jsonl"),
        help="Path to the phase 2 sweep output",
    )
    parser.add_argument(
        "--output-v1",
        type=Path,
        default=Path("data/branch_392_v1_verdict.jsonl"),
    )
    parser.add_argument(
        "--output-v2",
        type=Path,
        default=Path("data/branch_392_v2_verdict.jsonl"),
    )
    parser.add_argument(
        "--target",
        type=str,
        default=None,
        help="Filter by parent_label (e.g. C11a)",
    )
    args = parser.parse_args()

    branches = []
    with args.input.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("kind") == "branch_record":
                if args.target and row.get("parent_label") != args.target:
                    continue
                branches.append(row)

    if not branches:
        print("No branches found in input.")
        return

    # Filter for those that pass Hill screen (or MARGINAL)
    survivors = [b for b in branches if b.get("hill_classification") in ("PASS", "MARGINAL")]
    if not survivors:
        print("No branches passed the Hill screen.")
        return

    # Sort by Jacobi constant descending (highest C_J = lowest amplitude)
    survivors.sort(key=lambda b: float(b["branch_jacobi"]), reverse=True)
    best_candidate = survivors[0]
    candidate_id = best_candidate.get(
        "branch_id",
        f"branch_{best_candidate['parent_label']}_C_{best_candidate['branch_jacobi']:.4f}",
    )

    state0 = np.array(best_candidate["branch_state0"], dtype=np.float64)
    period_guess = float(best_candidate["branch_period_TU"])
    system = braik_ross_system()

    print(f"Top candidate: {candidate_id}")
    print(f"  Jacobi: {best_candidate['branch_jacobi']:.5f}")
    print(f"  Hill Fraction: {best_candidate.get('hill_fraction', 0.0):.3f}")

    # --- V1 ---
    print(f"\n[V1-verify] candidate_id={candidate_id}")
    t_start_v1 = time.time()

    result = correct_general_periodic_3d(
        system,
        state0,
        period_guess,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        residual_indices=RESIDUAL_FULL_STATE_AT_T,
        is_half_period_residual=False,
        tol=1e-12,
        max_iter=80,
        rtol=1e-12,
        atol=1e-12,
        independent_tol=V1_INDEPENDENT_CLOSURE_FLOOR,
        require_monotone_decrease=False,
    )
    elapsed_v1 = time.time() - t_start_v1

    passes_v1 = bool(
        result.converged
        and result.corrector_residual < V1_CORRECTOR_RESIDUAL_FLOOR
        and result.independent_closure_residual < V1_INDEPENDENT_CLOSURE_FLOOR
    )
    print(f"[V1-verify] passes_v1: {passes_v1}")

    iso_end_v1 = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output_v1.parent.mkdir(parents=True, exist_ok=True)
    with args.output_v1.open("w") as fh:
        verdict = {
            "kind": "v1_verdict_cr3bp_periodic",
            "candidate_id": candidate_id,
            "system": "Earth-Moon CR3BP",
            "mu": float(system.mu),
            "state0_corrected_nondim": [float(x) for x in result.state0],
            "period_corrected_TU": float(result.T_TU),
            "jacobi_corrected": float(result.jacobi),
            "corrector_residual": float(result.corrector_residual),
            "independent_closure_residual": float(result.independent_closure_residual),
            "converged": bool(result.converged),
            "passes_v1": passes_v1,
        }
        fh.write(json.dumps(verdict) + "\n")

    if not passes_v1:
        print("V1 failed. Halting before V2.")
        return

    # --- V2 ---
    print("\n[V2-verify] running V2 bounded cycle drift...")
    t_start_v2 = time.time()

    state0_v1 = result.state0
    period_v1 = result.T_TU

    verdicts_v2 = {}
    all_pass_v2 = True
    for n in N_CYCLES_LIST:
        print(f"[V2-verify] running n_cycles={n}...")
        v = run_v2_3d(
            candidate_id,
            state0_v1,
            period_v1,
            system,
            n_cycles=n,
            drift_floor_kms=V2_DRIFT_FLOOR_KMS,
            rtol=1e-12,
            atol=1e-12,
            notes=f"#392 Gauntlet V2 n_cycles={n} for {candidate_id}",
        )
        print(f"  passes_v2={v.passes_v2}, max_drift_km={v.max_drift_kms:.4e}")
        verdicts_v2[n] = {
            "kind": "v2_verdict_cr3bp_periodic",
            "candidate_id": candidate_id,
            "n_cycles_requested": int(v.n_cycles_requested),
            "n_cycles_propagated": int(v.n_cycles_propagated),
            "per_cycle_drift_kms": [float(x) for x in v.per_cycle_drift_kms],
            "max_drift_kms": float(v.max_drift_kms),
            "drift_floor_kms": float(v.drift_floor_kms),
            "converged_at_each_return": bool(v.converged_at_each_return),
            "passes_v2": bool(v.passes_v2),
        }
        all_pass_v2 = all_pass_v2 and bool(v.passes_v2)

    args.output_v2.parent.mkdir(parents=True, exist_ok=True)
    with args.output_v2.open("w") as fh:
        for n in N_CYCLES_LIST:
            fh.write(json.dumps(verdicts_v2[n]) + "\n")

    print(f"[V2-verify] all V2 gates pass: {all_pass_v2}")


if __name__ == "__main__":
    main()
