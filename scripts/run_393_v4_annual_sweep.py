"""#393 V4 annual launch-epoch sweep for the low-amplitude floquet cycler.

The #338 SILVER pattern applied to the new low-amplitude floquet candidate.
The V4 single-epoch run showed that the candidate FAILS V4 at the single
epoch 2000-01-15, but importantly, the EM-only control ALSO failed.
Because it survives V3 (circular Earth-Moon), its failure in V4 (eccentric
Earth-Moon) means it is highly sensitive to the true lunar eccentricity.

We scan over 100 launch epochs to see if there is any epoch where the lunar
anomaly aligns favorably and preserves the bounded drift.

Usage:
    uv run python scripts/run_393_v4_annual_sweep.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cyclerfinder.search.reachable_representatives import braik_ross_system
from scripts.run_393_v4_verify import (
    EPHEMERIS,
    V4_AGREEMENT_FLOOR_KMS,
    _load_v1,
    _load_v3_drifts,
    _run_ias15_n_cycles,
)

PHASE_LABEL = "393_v4_annual_sweep"
SWEEP_YEARS = range(2000, 2100)  # 100 epochs
DOY_LABEL = "06-21T00:00:00"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_393_v4_annual_sweep.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    candidate_id, state0, period = _load_v1()
    system = braik_ross_system()
    v3_drifts = _load_v3_drifts()

    passed_v3_cycles = [n for n, drifts in v3_drifts.items() if max(drifts) < 50000.0]
    passed_v3_cycles.sort()
    n_cycles = passed_v3_cycles[-1] if passed_v3_cycles else 1
    v3_n = v3_drifts.get(n_cycles, [])

    print(f"[V4-sweep] candidate_id={candidate_id}")
    print(f"[V4-sweep] {len(list(SWEEP_YEARS))} epochs, n_cycles={n_cycles}, ephemeris={EPHEMERIS}")

    rows: list[dict[str, Any]] = []
    n_pass = 0
    n_fail = 0
    for year in SWEEP_YEARS:
        epoch = f"{year}-{DOY_LABEL}"
        per_cycle_v4_km, per_cycle_conv, _ = _run_ias15_n_cycles(
            state0, period, n_cycles, system, epoch
        )
        agreement = [float(abs(v4 - v3)) for v4, v3 in zip(per_cycle_v4_km, v3_n, strict=False)]
        max_agreement = max(agreement) if agreement else float("inf")
        max_v4_drift = max(per_cycle_v4_km) if per_cycle_v4_km else float("inf")
        converged = all(per_cycle_conv)
        passes_v4 = bool(converged and max_agreement < V4_AGREEMENT_FLOOR_KMS)
        n_pass += int(passes_v4)
        n_fail += int(not passes_v4)
        rows.append(
            {
                "kind": "annual_sweep_row",
                "year": int(year),
                "launch_epoch_utc": epoch,
                "passes_v4": passes_v4,
                "n_cycles": n_cycles,
                "max_v4_drift_kms": float(max_v4_drift),
                "drift_agreement_kms": float(max_agreement),
                "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "converged_at_each_cycle": converged,
            }
        )

    drifts = np.array([r["max_v4_drift_kms"] for r in rows], dtype=np.float64)
    verdict_label = "STRUCTURAL_FAIL_ALL_EPOCHS" if n_pass == 0 else "MIXED"
    if n_pass == len(rows):
        verdict_label = "STRUCTURAL_PASS_ALL_EPOCHS"

    boundary = {
        "kind": "boundary_verdict",
        "candidate_id": candidate_id,
        "n_epochs": len(rows),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "verdict_label": verdict_label,
        "epoch_dependent": bool(0 < n_pass < len(rows)),
        "min_max_v4_drift_kms": float(drifts.min()),
        "max_max_v4_drift_kms": float(drifts.max()),
        "median_max_v4_drift_kms": float(np.median(drifts)),
    }

    elapsed = time.time() - t_start
    iso_end = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(
            json.dumps(
                {
                    "kind": "header",
                    "candidate_id": candidate_id,
                    "phase": PHASE_LABEL,
                    "iso_start": iso_start,
                    "iso_end": iso_end,
                    "elapsed_seconds": elapsed,
                    "ephemeris": EPHEMERIS,
                    "n_cycles_per_epoch": n_cycles,
                    "sweep_years": [SWEEP_YEARS.start, SWEEP_YEARS.stop - 1],
                    "doy_label": DOY_LABEL,
                    "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                }
            )
            + "\n"
        )
        for row in rows:
            fh.write(json.dumps(row) + "\n")
        fh.write(json.dumps(boundary) + "\n")
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "candidate_id": candidate_id,
                    "phase": PHASE_LABEL,
                    "iso_end": iso_end,
                    "verdict_label": verdict_label,
                    "n_pass": n_pass,
                    "n_fail": n_fail,
                }
            )
            + "\n"
        )
    print(f"[V4-sweep] {n_pass} PASS / {n_fail} FAIL across {len(rows)} epochs")
    print(f"[V4-sweep] verdict_label={verdict_label}")
    print(f"[V4-sweep] min max drift={drifts.min():.4e} km, max max drift={drifts.max():.4e} km")
    print(f"[V4-sweep] written to {args.output}")


if __name__ == "__main__":
    main()
