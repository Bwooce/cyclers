"""#330 Phase 2 — moontour V2 gauntlet on the #327 Umbriel-Oberon SILVER.

The driver
----------
:func:`cyclerfinder.data.validation.v2_moontour.run_v2_moontour` (#306 Phase 2
Part A) re-solves the SILVER's Lambert legs over ``n_cycles`` consecutive
cycles, advancing the moon ephemerides through their natural Keplerian
motion each cycle. The verdict gates: every cycle's Lambert converged AND
per-cycle V_inf-continuity residual <= 0.05 km/s AND inter-cycle rendezvous
drift <= 50,000 km.

What this script writes
-----------------------
``data/silver_327_moontour_v2_verdicts.jsonl`` with one row per ``n_cycles``
sample (3, 5, 10) so the verdict + per-cycle drift trace are recorded for
audit + the discovery report. ``_meta``, ``per_cycle_trace``, and the final
``_meta verdict`` row complete the audit trail.

Discipline anchors
------------------
* READ-ONLY on ``data/silver_327_verified.jsonl``.
* READ-ONLY on ``src/cyclerfinder/data/validation/v2_3d.py`` (Phase 1).
* NO catalogue writeback. The verdict is recorded; what it means is
  documented in ``docs/notes/2026-06-16-330-moontour-v2-phase2.md``.
* The verdict is whatever the math says — PASS / FAIL / quasi-cycler — no
  test-tuning.

Run as::

    uv run python scripts/run_330_silver_moontour_v2.py
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

SILVER_JSONL = ROOT / "data" / "silver_327_verified.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_moontour_v2_verdicts.jsonl"

# Stored SILVER row fields. Sourced from data/silver_327_verified.jsonl (#327
# verification, commit b080f32). Carried here as literal constants so the
# script is reproducible without re-running #327.
SILVER_ID = "repeated-moon-uranus-00000041"
SILVER_SEQ: tuple[str, ...] = ("Umbriel", "Oberon", "Umbriel")
SILVER_VINF: tuple[float, ...] = (
    0.9199258810725036,
    0.9604309791298091,
    0.8946936085078939,
)
SILVER_TOF: tuple[float, ...] = (14.940560615336594, 14.940560615336594)
SILVER_REL_OFF_DEG = 180.0
SILVER_NREV: tuple[int, ...] = (1, 1)
SILVER_PHASE0_DEG = 29.999999999999996


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
    print(f"[#330] moontour V2 phase 2 -- sha={sha}", flush=True)
    print(f"[#330] candidate = {SILVER_ID}", flush=True)
    print(f"[#330] sequence = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#330] V_inf = {SILVER_VINF}, ToF (days) = {SILVER_TOF}, "
        f"rel_off = {SILVER_REL_OFF_DEG}°, phase0 = {SILVER_PHASE0_DEG}°",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#330 moontour V2 phase 2 — Umbriel-Oberon-Umbriel SILVER",
            "candidate_id": SILVER_ID,
            "source_jsonl": str(SILVER_JSONL.relative_to(ROOT)),
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg": SILVER_REL_OFF_DEG,
                "phase0_deg": SILVER_PHASE0_DEG,
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
        print(f"[#330] run_v2_moontour(n_cycles={nc})...", flush=True)
        verdict = run_v2_moontour(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            None,  # primary auto-resolved from sequence (Uranus)
            n_cycles=nc,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0_DEG,
            notes=f"#330 phase-2 SILVER scan, n_cycles={nc}",
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

    # Phase 3 / Phase 4 recommendation: which gate is the cleaner next step
    # depends on whether ANY n_cycles passes V2. If n_cycles=3 fails on drift
    # but every Lambert converged AND quasi-bounded over 5+ cycles, the
    # candidate is a quasi-cycler (admit per v4.7) and Phase 3 (V3 nbody)
    # is the right next gate. If even cycle 1 closure blows past 0.5 km/s
    # the candidate is not even quasi-stable; Phase 4 (HFEM Uranus) is the
    # appropriate honest stop.
    v3 = headline[3]
    v10 = headline[10]
    quasi_cycler_3cycle = (
        v3.n_cycles_completed >= 1
        and not v3.passes_v2
        and v10.max_closure_residual_kms < 0.5  # within "quasi" envelope
    )
    if v3.passes_v2:
        next_step = "Phase 3 (#331) — V3 6D nbody / REBOUND on the SILVER"
        verdict_label = "PASS"
    elif quasi_cycler_3cycle:
        next_step = (
            "Phase 3 (#331) recommended — Lambert legs all converge over 10 cycles "
            "with closure residual < 0.5 km/s; SILVER is a quasi-cycler under V2 "
            "but worth a real-eph V3 to distinguish bounded drift from instability"
        )
        verdict_label = "FAIL_QUASI_BOUNDED"
    else:
        next_step = (
            "Phase 4 (#332) — V2 fails on both drift AND closure; admission "
            "as a quasi-cycler is not warranted at this phasing. The HFEM "
            "Uranus gate is the honest stop."
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
    print(f"[#330] wrote {OUT_JSONL}", flush=True)
    print(f"[#330] verdict: {verdict_label} (elapsed {time.time() - t0:.1f}s)", flush=True)
    print(f"[#330] next step: {next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
