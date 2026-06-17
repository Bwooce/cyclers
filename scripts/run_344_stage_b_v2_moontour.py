"""#344 Phase 2 Stage B - moontour V2 gauntlet on the Saturn Titan-Rhea-Titan SILVER.

The driver
----------
:func:`cyclerfinder.data.validation.v2_moontour.run_v2_moontour` (#306 Phase 2
Part A, the same driver used by #330 on the Umbriel-Oberon-Umbriel SILVER)
re-solves the SILVER's Lambert legs over ``n_cycles`` consecutive cycles,
advancing the moon ephemerides through their natural Keplerian motion each
cycle. The verdict gates: every cycle's Lambert converged AND per-cycle
V_inf-continuity residual <= 0.05 km/s AND inter-cycle rendezvous drift <=
50,000 km.

What this script writes
-----------------------
``data/silver_344_moontour_v2_verdicts.jsonl`` with one row per ``n_cycles``
sample (3, 5, 10) so the verdict + per-cycle drift trace are recorded for
audit + the discovery report. ``_meta``, ``per_cycle_trace``, and the final
``_meta verdict`` row complete the audit trail.

Discipline anchors
------------------
* READ-ONLY on ``data/silver_344_verified.jsonl``.
* READ-ONLY on ``src/cyclerfinder/data/validation/v2_moontour.py``.
* NO catalogue writeback. The verdict is recorded; what it means is
  documented in ``docs/notes/2026-06-17-344-phase2-stage-b-v2-moontour.md``.
* The verdict is whatever the math says - PASS / FAIL / quasi-cycler - no
  test-tuning.

Run as::

    uv run python scripts/run_344_stage_b_v2_moontour.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.validation.v2_moontour import (  # noqa: E402
    V2_MOONTOUR_CLOSURE_FLOOR_KMS,
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    V2_MOONTOUR_N_CYCLES_MIN,
    V2MoontourVerdict,
    run_v2_moontour,
)

SILVER_JSONL = ROOT / "data" / "silver_344_verified.jsonl"
OUT_JSONL = ROOT / "data" / "silver_344_moontour_v2_verdicts.jsonl"

# Stored SILVER row fields. Sourced from data/silver_344_verified.jsonl
# (#344 Stage A verification, commit 63809ec on main). Carried here as
# literal constants so the script is reproducible without re-running
# Stage A.
SILVER_ID = "repeated-moon-saturn-titan-rhea-titan-stage-b"
SILVER_SEQ: tuple[str, ...] = ("Titan", "Rhea", "Titan")
SILVER_VINF: tuple[float, ...] = (
    1.7375055995850324,
    1.6462740278228238,
    1.7273175030110421,
)
SILVER_TOF: tuple[float, ...] = (16.977266455394638, 16.977266455394638)
# Stage-A storage convention (from scripts/scan_344_saturn_titan_rhea_finer.py
# ``_sweep_one_cycle``): theta = {anchor: phase0, intermediate: phase0+rel_off}.
# For the (Titan, Rhea, Titan) cycle the ANCHOR is Titan and the INTERMEDIATE
# is Rhea, so Titan = 273.75 deg and Rhea = (273.75 + 288.75) mod 360 = 202.5 deg.
SILVER_STORED_PHASE0_DEG = 273.74999999999994
SILVER_STORED_REL_OFFSET_DEG = 288.75
# v2_moontour.run_v2_moontour convention (per its module docstring): the
# offset is applied to the SECOND distinct moon in registry-sorted order.
# distinct_moons = sorted({"Titan", "Rhea"}) = ("Rhea", "Titan"), so
# Rhea = phase0_v2 and Titan = phase0_v2 + rel_off_v2.
# To reproduce the Stage A SILVER geometry under this convention:
#   phase0_v2 = Rhea_pos = 202.5 deg
#   rel_off_v2 = (Titan_pos - Rhea_pos) mod 360 = (273.75 - 202.5) mod 360
#              = 71.25 deg
# The 0.025 km/s cycle-0 reproduction in #330 worked because Umbriel-Oberon
# had stored rel_off = 180 deg (palindromic under the swap); #344 has stored
# rel_off = 288.75 deg so the swap matters and the conversion above is
# load-bearing.
SILVER_PHASE0_DEG = 202.5
SILVER_REL_OFF_DEG = 71.25
SILVER_NREV: tuple[int, ...] = (1, 1)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _verdict_to_row(verdict: V2MoontourVerdict) -> dict[str, Any]:
    """Flatten a V2MoontourVerdict into a JSONL row."""
    return {
        "kind": "moontour_v2_verdict",
        "candidate_id": verdict.candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_requested": verdict.n_cycles_requested,
        "n_cycles_completed": verdict.n_cycles_completed,
        "drift_floor_kms": verdict.drift_floor_kms,
        "closure_floor_kms": verdict.closure_floor_kms,
        "n_cycles_min": verdict.n_cycles_min,
        "max_drift_kms": verdict.max_drift_kms,
        "max_drift_seconds": verdict.max_drift_seconds,
        "max_closure_residual_kms": verdict.max_closure_residual_kms,
        "passes_v2": verdict.passes_v2,
        "per_cycle": [
            {
                "cycle_index": c.cycle_index,
                "converged_legs": c.converged_legs,
                "n_legs": c.n_legs,
                "rendezvous_drift_kms": c.rendezvous_drift_kms,
                "rendezvous_drift_seconds": c.rendezvous_drift_seconds,
                "closure_residual_kms": c.closure_residual_kms,
                "notes": c.notes,
            }
            for c in verdict.per_cycle
        ],
        "notes": verdict.notes,
    }


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#344-stage-b] moontour V2 phase 2 -- sha={sha}", flush=True)
    print(f"[#344-stage-b] candidate = {SILVER_ID}", flush=True)
    print(f"[#344-stage-b] sequence = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#344-stage-b] V_inf = {SILVER_VINF}, ToF (days) = {SILVER_TOF}, "
        f"rel_off = {SILVER_REL_OFF_DEG} deg, phase0 = {SILVER_PHASE0_DEG} deg",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#344 Phase 2 Stage B - V2 moontour gauntlet on Titan-Rhea-Titan SILVER",
            "candidate_id": SILVER_ID,
            "source_jsonl": str(SILVER_JSONL.relative_to(ROOT)),
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg_stage_a_convention": SILVER_STORED_REL_OFFSET_DEG,
                "phase0_deg_stage_a_convention": SILVER_STORED_PHASE0_DEG,
                "rel_offset_deg_v2_convention": SILVER_REL_OFF_DEG,
                "phase0_deg_v2_convention": SILVER_PHASE0_DEG,
                "phase_convention_note": (
                    "Stage A uses _sweep_one_cycle's {anchor: phase0, "
                    "intermediate: phase0+rel_off}; v2_moontour uses "
                    "{registry-sorted moon 0: phase0, moon 1: phase0+rel_off}. "
                    "Conversion: Titan_pos=273.75, Rhea_pos=202.5 -> v2 "
                    "phase0=202.5 (Rhea, sorted-first), rel_off=71.25 (Titan-Rhea)."
                ),
                "n_rev": list(SILVER_NREV),
            },
            "driver_floors": {
                "drift_floor_kms": V2_MOONTOUR_DRIFT_FLOOR_KMS,
                "closure_floor_kms": V2_MOONTOUR_CLOSURE_FLOOR_KMS,
                "n_cycles_min": V2_MOONTOUR_N_CYCLES_MIN,
            },
            "git_sha": sha,
        }
    )

    n_cycles_grid = (3, 5, 10)
    headline: dict[int, V2MoontourVerdict] = {}
    for nc in n_cycles_grid:
        print(f"[#344-stage-b] run_v2_moontour(n_cycles={nc})...", flush=True)
        verdict = run_v2_moontour(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            None,  # primary auto-resolved from sequence (Saturn)
            n_cycles=nc,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0_DEG,
            notes=f"#344 phase-2 stage-b SILVER scan, n_cycles={nc}",
        )
        headline[nc] = verdict
        print(
            f"  -> passes_v2={verdict.passes_v2} | "
            f"completed={verdict.n_cycles_completed}/{verdict.n_cycles_requested} | "
            f"max_drift_kms={verdict.max_drift_kms:.3e} | "
            f"max_closure_kms={verdict.max_closure_residual_kms:.3e}",
            flush=True,
        )
        for c in verdict.per_cycle:
            print(
                f"     cycle {c.cycle_index}: drift={c.rendezvous_drift_kms:.3e} km, "
                f"residual={c.closure_residual_kms:.3e} km/s, "
                f"converged_legs={c.converged_legs}/{c.n_legs}",
                flush=True,
            )
        rows.append(_verdict_to_row(verdict))

    # Verdict label rules (mirror #330):
    # * PASS_STRICT_CYCLER iff strict V2 passes at n_cycles=3
    #   (per-cycle closure < 0.05 km/s AND drift < 50,000 km over 3 cycles).
    # * FAIL_QUASI_BOUNDED iff strict gate FAILS but every Lambert converges
    #   over 10 cycles AND closure residual stays below the v4.7 quasi
    #   envelope (< 0.5 km/s).
    # * FAIL otherwise (closure blows past 0.5 km/s, or Lambert fails to
    #   converge mid-tour).
    v3 = headline[3]
    v10 = headline[10]
    quasi_cycler_3cycle = (
        v3.n_cycles_completed >= 1
        and not v3.passes_v2
        and v10.n_cycles_completed == 10
        and v10.max_closure_residual_kms < 0.5
    )
    if v3.passes_v2:
        next_step = (
            "Stage C (V3 REBOUND IAS15) - strict V2 cycler verdict; the V3 "
            "real-ephemeris gate still applies before catalogue admission."
        )
        verdict_label = "PASS_STRICT_CYCLER"
    elif quasi_cycler_3cycle:
        next_step = (
            "Stage C (V3 REBOUND IAS15) recommended - Lambert legs all converge "
            "over 10 cycles with closure residual < 0.5 km/s; SILVER is a "
            "quasi-cycler under V2 but warrants a real-eph V3 to distinguish "
            "bounded drift from instability."
        )
        verdict_label = "PASS_QUASI_CYCLER"
    else:
        next_step = (
            "HALT - V2 fails on both drift AND closure beyond the v4.7 quasi "
            "envelope; admission as a quasi-cycler is not warranted at this "
            "phasing. Candidate retires to the negative-results registry."
        )
        verdict_label = "FAIL"

    rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "candidate_id": SILVER_ID,
            "n_cycles_3": {
                "passes_v2": v3.passes_v2,
                "max_drift_kms": v3.max_drift_kms,
                "max_closure_kms": v3.max_closure_residual_kms,
                "n_completed": v3.n_cycles_completed,
            },
            "n_cycles_5": {
                "passes_v2": headline[5].passes_v2,
                "max_drift_kms": headline[5].max_drift_kms,
                "max_closure_kms": headline[5].max_closure_residual_kms,
                "n_completed": headline[5].n_cycles_completed,
            },
            "n_cycles_10": {
                "passes_v2": v10.passes_v2,
                "max_drift_kms": v10.max_drift_kms,
                "max_closure_kms": v10.max_closure_residual_kms,
                "n_completed": v10.n_cycles_completed,
            },
            "verdict_label": verdict_label,
            "writeback_to_catalogue": False,
            "next_step": next_step,
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    print(f"[#344-stage-b] wrote {OUT_JSONL}", flush=True)
    print(
        f"[#344-stage-b] verdict: {verdict_label} (elapsed {time.time() - t0:.1f}s)",
        flush=True,
    )
    print(f"[#344-stage-b] next step: {next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
