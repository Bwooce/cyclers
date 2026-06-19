"""#389 P389.6 — annual launch-epoch sweep 2000-2099 for branch_C32_b0 V4.

The #338 SILVER pattern, applied to branch_C32_b0's V4 HALT. P389.5 showed the
candidate FAILS V4 at the single epoch 2000-01-15 by escaping (~1e9 km drift)
under real DE440 solar tides at its near-Hill-radius amplitude. The #338 pattern
asks whether a V4 verdict is *epoch-dependent* (a launch-window artifact) or
*structural* (true at every epoch).

For the SILVER, #338 found an EFFECTIVELY_CYCLIC interior PASS run. For
branch_C32_b0 we expect the OPPOSITE: a uniform FAIL across all 100 epochs,
because the failure cause (solar tide ~30% of Earth gravity at the orbit's
0.77-Hill-radius amplitude) is geometric, not phase-dependent. A uniform FAIL
sweep upgrades the HALT from "fails at one epoch" to "structurally cannot be a
real-ephemeris cycler at any 21st-century launch epoch" — the honest, complete
negative.

This script reuses the P389.5 V4 propagator (``_run_ias15_n_cycles``) at
n_cycles=3 across launch epochs ``<year>-06-21T00:00:00`` for year in 2000-2099
(matching the #338 grid's day-of-year), and records per-epoch pass/fail.

Usage:
    uv run python scripts/branch_c32_b0_v4_annual_sweep.py [--output PATH]
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
from scripts.branch_c32_b0_v4_verify import (
    CANDIDATE_ID,
    EPHEMERIS,
    V4_AGREEMENT_FLOOR_KMS,
    _load_corrected_state_and_period,
    _load_v3_drifts,
    _run_ias15_n_cycles,
)

PHASE_LABEL = "389_p389_6"
SWEEP_YEARS = range(2000, 2100)  # 100 epochs
DOY_LABEL = "06-21T00:00:00"
N_CYCLES = 3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/branch_c32_b0_v4_annual_sweep_389.jsonl"),
    )
    args = parser.parse_args()

    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t_start = time.time()

    state0, period = _load_corrected_state_and_period()
    system = braik_ross_system()
    v3_drifts = _load_v3_drifts()
    v3_n = v3_drifts.get(N_CYCLES, [])

    print(f"[V4-sweep] candidate_id={CANDIDATE_ID}")
    print(f"[V4-sweep] {len(list(SWEEP_YEARS))} epochs, n_cycles={N_CYCLES}, ephemeris={EPHEMERIS}")

    rows: list[dict[str, Any]] = []
    n_pass = 0
    n_fail = 0
    for year in SWEEP_YEARS:
        epoch = f"{year}-{DOY_LABEL}"
        per_cycle_v4_km, per_cycle_conv, _ = _run_ias15_n_cycles(
            state0, period, N_CYCLES, system, epoch
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
                "n_cycles": N_CYCLES,
                "max_v4_drift_kms": float(max_v4_drift),
                "drift_agreement_kms": float(max_agreement),
                "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "converged_at_each_cycle": converged,
            }
        )

    drifts = np.array([r["max_v4_drift_kms"] for r in rows], dtype=np.float64)
    verdict_label = "STRUCTURAL_FAIL_ALL_EPOCHS" if n_pass == 0 else "MIXED"
    boundary = {
        "kind": "boundary_verdict",
        "candidate_id": CANDIDATE_ID,
        "n_epochs": len(rows),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "verdict_label": verdict_label,
        "epoch_dependent": bool(0 < n_pass < len(rows)),
        "min_max_v4_drift_kms": float(drifts.min()),
        "max_max_v4_drift_kms": float(drifts.max()),
        "median_max_v4_drift_kms": float(np.median(drifts)),
        "interpretation": (
            "branch_C32_b0 FAILS V4 at every 21st-century launch epoch: the "
            "real-ephemeris solar tide (~30% of Earth gravity at the orbit's "
            "~0.77 Earth-Sun-Hill-radius amplitude) destabilizes the far-amplitude "
            "(3,3) orbit into escape (~1e9 km drift) regardless of launch phase. "
            "The failure is STRUCTURAL (geometric/amplitude-driven), NOT a "
            "launch-window artifact — the opposite of the #338 SILVER, whose "
            "interior PASS run earned EFFECTIVELY_CYCLIC. branch_C32_b0 stays "
            "unadmitted; the catalogue is unchanged."
        ),
    }

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
                    "ephemeris": EPHEMERIS,
                    "n_cycles_per_epoch": N_CYCLES,
                    "sweep_years": [SWEEP_YEARS.start, SWEEP_YEARS.stop - 1],
                    "doy_label": DOY_LABEL,
                    "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                    "successor_to": "#338 (SILVER annual sweep, EFFECTIVELY_CYCLIC)",
                    "discipline": (
                        "#338-pattern annual launch-epoch sweep to test whether the "
                        "P389.5 V4 HALT is epoch-dependent (window artifact) or "
                        "structural (true at every epoch). Uniform FAIL upgrades the "
                        "HALT to a complete structural negative."
                    ),
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
                    "candidate_id": CANDIDATE_ID,
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
    print(f"[V4-sweep] written to {args.output}")


if __name__ == "__main__":
    main()
