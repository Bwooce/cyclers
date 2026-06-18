"""#389 P389.3 — V2 bounded-cycle gate for branch_C32_b0 (CR3BP-adapted).

For a CR3BP periodic orbit, the V2 question is bounded-drift survival under
DOP853 propagation for n consecutive cycles WITHOUT recorrecting. The orbit's
own essentially-stable Floquet character (max_floquet_mag = 1.000000000000617,
sigma_d/day = 6.08e-15) predicts that the per-cycle drift should remain at
the level of the integrator's cumulative round-off — well below the spec §14
50,000 km same-model floor.

This is the bounded-cycle CR3BP V2 (``v2_3d.py``), NOT the moontour Lambert-
relegs V2 (``v2_moontour.py``). branch_C32_b0 is a closed CR3BP orbit, not a
patched-conic tour.

Spec §14 V2 floor: 50,000 km (same-model bar). Cycles: {3, 5, 10}.

Usage:
    uv run python scripts/branch_c32_b0_v2_verify.py [--output PATH]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.validation.v2_3d import (
    V2_DRIFT_FLOOR_KMS,
    V2_N_CYCLES_MIN,
    run_v2_3d,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
V1_VERDICT_PATH = Path("data/branch_c32_b0_v1_verdict.jsonl")
N_CYCLES_LIST = (3, 5, 10)
PHASE_LABEL = "389_p389_3"


def _load_corrected_state_and_period() -> tuple[np.ndarray, float]:
    """Use the V1-corrected state + period (the V1-converged closure)."""
    with V1_VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v1_verdict_cr3bp_periodic"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                state = np.array(row["state0_corrected_nondim"], dtype=np.float64)
                period = float(row["period_corrected_TU"])
                return state, period
    raise AssertionError(f"V1 verdict row for {CANDIDATE_ID!r} not found in {V1_VERDICT_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_c32_b0_v2_verdict.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    state0, period = _load_corrected_state_and_period()
    system = braik_ross_system()
    print(f"[V2-verify] candidate_id={CANDIDATE_ID}")
    print(f"[V2-verify] V1-corrected state0={state0}")
    print(f"[V2-verify] V1-corrected period_TU={period:.15f}")
    print(f"[V2-verify] system.l_km={system.l_km}")
    print(f"[V2-verify] V2 floor: {V2_DRIFT_FLOOR_KMS:.0f} km, n_cycles_min={V2_N_CYCLES_MIN}")

    verdicts: dict[int, dict[str, Any]] = {}
    all_pass = True
    for n in N_CYCLES_LIST:
        print(f"[V2-verify] running n_cycles={n}...")
        v = run_v2_3d(
            CANDIDATE_ID,
            state0,
            period,
            system,
            n_cycles=n,
            drift_floor_kms=V2_DRIFT_FLOOR_KMS,
            rtol=1e-12,
            atol=1e-12,
            notes=f"#389 P389.3 V2 n_cycles={n} for branch_C32_b0",
        )
        per_cycle_kms = [float(x) for x in v.per_cycle_drift_kms]
        per_cycle_vel = [float(x) for x in v.per_cycle_velocity_drift_kms]
        print(
            f"  passes_v2={v.passes_v2}, max_drift_km={v.max_drift_kms:.4e}, "
            f"converged_at_each_return={v.converged_at_each_return}"
        )
        print(f"  per_cycle_drift_km={per_cycle_kms}")
        verdicts[n] = {
            "kind": "v2_verdict_cr3bp_periodic",
            "candidate_id": CANDIDATE_ID,
            "n_cycles_requested": int(v.n_cycles_requested),
            "n_cycles_propagated": int(v.n_cycles_propagated),
            "per_cycle_drift_kms": per_cycle_kms,
            "per_cycle_velocity_drift_kms": per_cycle_vel,
            "max_drift_kms": float(v.max_drift_kms),
            "drift_floor_kms": float(v.drift_floor_kms),
            "n_cycles_min": int(v.n_cycles_min),
            "converged_at_each_return": bool(v.converged_at_each_return),
            "passes_v2": bool(v.passes_v2),
            "notes": v.notes,
        }
        all_pass = all_pass and bool(v.passes_v2)

    elapsed = time.time() - t_start
    iso_end = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(
            json.dumps(
                {
                    "kind": "header",
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "v2_drift_floor_kms": V2_DRIFT_FLOOR_KMS,
                    "n_cycles_min": V2_N_CYCLES_MIN,
                    "n_cycles_list": list(N_CYCLES_LIST),
                    "discipline": (
                        "Spec §14 V2 bounded-cycle CR3BP propagation: cumulative "
                        "position drift at each cycle boundary held under the "
                        f"{V2_DRIFT_FLOOR_KMS:.0f} km same-model floor across "
                        f"n_cycles in {tuple(N_CYCLES_LIST)}. NOT a moontour V2 "
                        "(branch_C32_b0 is a closed CR3BP periodic orbit, not "
                        "a patched-Lambert tour)."
                    ),
                }
            )
            + "\n"
        )
        for n in N_CYCLES_LIST:
            fh.write(json.dumps(verdicts[n]) + "\n")
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "candidate_id": CANDIDATE_ID,
                    "phase": PHASE_LABEL,
                    "iso_end": iso_end,
                    "all_pass": all_pass,
                }
            )
            + "\n"
        )
    print(f"[V2-verify] verdict written to {args.output}")
    print(f"[V2-verify] all V2 gates pass: {all_pass}")


if __name__ == "__main__":
    main()
