"""#331 Phase 3 — V3 6D n-body gauntlet on the #327 Umbriel-Oberon SILVER.

The driver
----------
:func:`cyclerfinder.data.validation.v3_3d.run_v3_3d` (#306 Phase 3 Part B)
re-propagates each cycle of the SILVER's Lambert tour under REBOUND IAS15
(Gauss-Radau) in the planet frame. The verdict gates: REBOUND IAS15 per-
cycle terminal positions agree with the V2 driver's (DOP853 + Lambert)
per-cycle terminal positions to within 100 km.

A V3 PASS means: the V2 quasi_cycler bounded-drift signature reported by
#330 (commit ff6a4ad) is a REAL property of the shared circular-coplanar
Keplerian model — NOT an artifact of the V2 driver's DOP853 + Lambert
internal numerics.

A V3 FAIL would mean: the V2 bounded-drift was integrator noise, and the
SILVER retires to the negative-results registry (#172).

What this script writes
-----------------------
``data/silver_327_v3_verdicts.jsonl`` with one row per ``n_cycles`` sample
(3, 5, 10) so the V3 / V2 per-cycle comparison is recorded for audit +
the Phase 3 discovery report. ``_meta``, per-cycle drift+agreement
traces, and the final ``_meta headline`` row complete the audit trail.

Discipline anchors
------------------
* READ-ONLY on ``data/silver_327_verified.jsonl``.
* READ-ONLY on ``data/silver_327_moontour_v2_verdicts.jsonl`` (#330).
* READ-ONLY on ``src/cyclerfinder/data/validation/v2_3d.py``,
  ``v2_moontour.py``, ``v1_3d.py`` (Phases 1+2).
* NO catalogue writeback. The verdict is recorded; what it means is
  documented in ``docs/notes/2026-06-16-331-v3-nbody-phase3.md``.
* The verdict is whatever the math says — PASS / FAIL — no test-tuning.

Run as::

    uv run python scripts/run_331_silver_v3_gauntlet.py
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
    V2MoontourVerdict,
    run_v2_moontour,
)
from cyclerfinder.data.validation.v3_3d import (  # noqa: E402
    V3_AGREEMENT_FLOOR_KMS,
    V3_N_CYCLES_MIN,
    V3Verdict3D,
    run_v3_3d,
)

SILVER_VERIFIED_JSONL = ROOT / "data" / "silver_327_verified.jsonl"
V2_VERDICTS_JSONL = ROOT / "data" / "silver_327_moontour_v2_verdicts.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_v3_verdicts.jsonl"

# Stored SILVER row fields. Sourced from data/silver_327_verified.jsonl
# (commit b080f32) and data/silver_327_moontour_v2_verdicts.jsonl
# (commit ff6a4ad). Carried as literal constants so the script is
# self-contained.
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


def _v3_verdict_to_row(verdict: V3Verdict3D) -> dict[str, Any]:
    """Flatten a V3Verdict3D into a JSONL row."""
    return {
        "kind": "moontour_v3_verdict",
        "candidate_id": verdict.candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "drift_agreement_kms": verdict.drift_agreement_kms,
        "v3_v2_agreement_floor_kms": verdict.v3_v2_agreement_floor_kms,
        "passes_v3": verdict.passes_v3,
        "per_cycle_drift_kms_v3": list(verdict.per_cycle_drift_kms_v3),
        "per_cycle_drift_kms_v2": list(verdict.per_cycle_drift_kms_v2),
        "per_cycle": [
            {
                "cycle_index": c.cycle_index,
                "converged_legs": c.converged_legs,
                "n_legs": c.n_legs,
                "rendezvous_drift_kms_v3": c.rendezvous_drift_kms_v3,
                "rendezvous_drift_kms_v2": c.rendezvous_drift_kms_v2,
                "agreement_kms": c.agreement_kms,
                "ias15_vs_analytic_kepler_kms": c.ias15_vs_analytic_kepler_kms,
                "notes": c.notes,
            }
            for c in verdict.per_cycle
        ],
        "notes": verdict.notes,
    }


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#331] V3 nbody phase 3 -- sha={sha}", flush=True)
    print(f"[#331] candidate = {SILVER_ID}", flush=True)
    print(f"[#331] sequence = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#331] V_inf = {SILVER_VINF}, ToF (days) = {SILVER_TOF}, "
        f"rel_off = {SILVER_REL_OFF_DEG}°, phase0 = {SILVER_PHASE0_DEG}°",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#331 V3 nbody phase 3 — Umbriel-Oberon-Umbriel SILVER",
            "candidate_id": SILVER_ID,
            "source_jsonl": str(SILVER_VERIFIED_JSONL.relative_to(ROOT)),
            "v2_verdicts_jsonl": str(V2_VERDICTS_JSONL.relative_to(ROOT)),
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg": SILVER_REL_OFF_DEG,
                "phase0_deg": SILVER_PHASE0_DEG,
                "n_rev": list(SILVER_NREV),
            },
            "driver_floors": {
                "agreement_floor_kms": V3_AGREEMENT_FLOOR_KMS,
                "n_cycles_min": V3_N_CYCLES_MIN,
            },
            "git_sha": sha,
        }
    )

    n_cycles_grid = (3, 5, 10)
    headline_v2: dict[int, V2MoontourVerdict] = {}
    headline_v3: dict[int, V3Verdict3D] = {}
    for nc in n_cycles_grid:
        print(f"[#331] run_v2_moontour(n_cycles={nc})...", flush=True)
        v2 = run_v2_moontour(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            None,
            n_cycles=nc,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0_DEG,
            notes=f"#331 phase-3 V2 input, n_cycles={nc}",
        )
        headline_v2[nc] = v2
        print(f"[#331] run_v3_3d(n_cycles={nc})...", flush=True)
        v3 = run_v3_3d(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            None,
            v2_verdict=v2,
            n_cycles=nc,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0_DEG,
            notes=f"#331 phase-3 V3 scan, n_cycles={nc}",
        )
        headline_v3[nc] = v3
        print(
            f"  -> passes_v3={v3.passes_v3} | "
            f"integrator={v3.integrator} | "
            f"completed={v3.n_cycles_propagated}/{nc} | "
            f"agreement_kms={v3.drift_agreement_kms:.3e}",
            flush=True,
        )
        for c in v3.per_cycle:
            print(
                f"     cycle {c.cycle_index}: V3_drift={c.rendezvous_drift_kms_v3:.3e} km, "
                f"V2_drift={c.rendezvous_drift_kms_v2:.3e} km, "
                f"|V3-V2|={c.agreement_kms:.3e} km, "
                f"IAS15-vs-Kepler={c.ias15_vs_analytic_kepler_kms:.3e} km",
                flush=True,
            )
        rows.append(_v3_verdict_to_row(v3))

    # V3 headline + phase-4 recommendation. The V3 PASS at all three cycle
    # counts means the V2 quasi_cycler bounded-drift signature survives an
    # independent integrator architecture: it is a REAL property of the
    # circular-coplanar Keplerian model. The next gate is V4 (#332 HFEM
    # Uranian-system real-ephemeris with SPICE) + #329 Heaton-Longuski
    # literature check.
    v3_3 = headline_v3[3]
    v3_5 = headline_v3[5]
    v3_10 = headline_v3[10]
    all_pass = v3_3.passes_v3 and v3_5.passes_v3 and v3_10.passes_v3
    any_fail = (not v3_3.passes_v3) or (not v3_5.passes_v3) or (not v3_10.passes_v3)

    if all_pass:
        verdict_label = "PASS"
        next_step = (
            "Phase 4 (#332) — V4 HFEM Uranian-system real-ephemeris via GMAT + SPICE. "
            "REBOUND IAS15 confirms the V2 bounded-drift signature is a real property "
            "of the shared model; the V4 question is whether it survives real-eph "
            "Uranian dynamics. Also run #329 Heaton-Longuski lit-check before "
            "claiming a catalogue novelty."
        )
    elif any_fail:
        verdict_label = "FAIL"
        next_step = (
            "Retire to negative-results registry (#172) — V2 bounded-drift was "
            "integrator artefact (REBOUND IAS15 disagrees with V2 DOP853+Lambert "
            "beyond the agreement floor). The SILVER is NOT a quasi-cycler under "
            "an integrator-independent reading of the shared model."
        )
    else:  # pragma: no cover — defensive
        verdict_label = "MIXED"
        next_step = "Investigate per-cycle disagreement before Phase 4."

    rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "candidate_id": SILVER_ID,
            "v2_max_drift_kms_at_n3": headline_v2[3].max_drift_kms,
            "v2_max_drift_kms_at_n5": headline_v2[5].max_drift_kms,
            "v2_max_drift_kms_at_n10": headline_v2[10].max_drift_kms,
            "v3_n_cycles_3": {
                "integrator": v3_3.integrator,
                "passes_v3": v3_3.passes_v3,
                "agreement_kms": v3_3.drift_agreement_kms,
                "n_completed": v3_3.n_cycles_propagated,
            },
            "v3_n_cycles_5": {
                "integrator": v3_5.integrator,
                "passes_v3": v3_5.passes_v3,
                "agreement_kms": v3_5.drift_agreement_kms,
                "n_completed": v3_5.n_cycles_propagated,
            },
            "v3_n_cycles_10": {
                "integrator": v3_10.integrator,
                "passes_v3": v3_10.passes_v3,
                "agreement_kms": v3_10.drift_agreement_kms,
                "n_completed": v3_10.n_cycles_propagated,
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
    print(f"[#331] wrote {OUT_JSONL}", flush=True)
    print(f"[#331] verdict: {verdict_label} (elapsed {time.time() - t0:.1f}s)", flush=True)
    print(f"[#331] next step: {next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
