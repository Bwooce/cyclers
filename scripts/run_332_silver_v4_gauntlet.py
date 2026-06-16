"""#332 Phase 4 — V4 HFEM-fallback Uranus gauntlet on the #327 SILVER.

The driver
----------
:func:`cyclerfinder.data.validation.v4_uranus.run_v4_uranus` (#306 Phase 4
Part B, commit 7d43631) re-propagates each cycle of the SILVER's Lambert tour
under scipy DOP853 + Uranus J2 + classical-moon third-body perturbations
(the documented scipy-fallback path — GMAT R2022a install lacks Uranian
satellite SPICE kernels; see ``src/cyclerfinder/data/validation/v4_uranus.py``
module docstring "GMAT-vs-fallback rationale" for the full rationale).

A V4 PASS means: the V3-confirmed bounded-drift signature of the SILVER
(#331 commit cae57ca) survives the dominant non-Keplerian perturbations of
the Uranian system at V4-class fidelity. A V4 FAIL means the V2/V3
Keplerian-coplanar bounded-drift was a model idealization that does not
survive J2 + other-moon n-body — retire to the negative-results registry
(#172) recording which order of perturbation breaks the signature.

Either outcome is the right session-closing result for #332.

Discipline anchors
------------------
* READ-ONLY on ``data/silver_327_verified.jsonl``, ``data/silver_327_v3_verdicts.jsonl``.
* READ-ONLY on ``src/cyclerfinder/data/validation/v[1-3]_*.py`` (Phases 1-3).
* NO catalogue writeback. The verdict is recorded; what it means is
  documented in ``docs/notes/2026-06-16-332-v4-uranus-phase4.md``.
* The verdict is whatever the math says — PASS / FAIL — no test-tuning.
* This is the scipy-fallback V4 (Part D — Phase 4.1 GMAT V4 with full
  Uranian satellite SPICE remains a downstream catalogue-admission gate
  if/when those kernels are installed).

Run as::

    uv run python scripts/run_332_silver_v4_gauntlet.py
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
    run_v2_moontour,
)
from cyclerfinder.data.validation.v3_3d import (  # noqa: E402
    V3Verdict3D,
    run_v3_3d,
)
from cyclerfinder.data.validation.v4_uranus import (  # noqa: E402
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4UranusVerdict,
    run_v4_uranus,
)

SILVER_VERIFIED_JSONL = ROOT / "data" / "silver_327_verified.jsonl"
V3_VERDICTS_JSONL = ROOT / "data" / "silver_327_v3_verdicts.jsonl"
OUT_JSONL = ROOT / "data" / "silver_327_v4_verdicts.jsonl"

# Stored SILVER row fields. Sourced from data/silver_327_verified.jsonl
# (commit b080f32) and data/silver_327_v3_verdicts.jsonl (commit cae57ca).
# Carried as literal constants so the script is self-contained.
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


def _v4_verdict_to_row(verdict: V4UranusVerdict) -> dict[str, Any]:
    """Flatten a V4UranusVerdict into a JSONL row."""
    return {
        "kind": "moontour_v4_verdict",
        "candidate_id": verdict.candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "drift_agreement_kms": verdict.drift_agreement_kms,
        "v4_v3_agreement_floor_kms": verdict.v4_v3_agreement_floor_kms,
        "bounded_drift_survives": verdict.bounded_drift_survives,
        "passes_v4": verdict.passes_v4,
        "per_cycle_drift_kms_v4": list(verdict.per_cycle_drift_kms_v4),
        "per_cycle_drift_kms_v3": list(verdict.per_cycle_drift_kms_v3),
        "per_cycle": [
            {
                "cycle_index": c.cycle_index,
                "converged_legs": c.converged_legs,
                "n_legs": c.n_legs,
                "rendezvous_drift_kms_v4": c.rendezvous_drift_kms_v4,
                "rendezvous_drift_kms_v3": c.rendezvous_drift_kms_v3,
                "agreement_kms": c.agreement_kms,
                "v4_terminal_offset_vs_moon_kms": c.v4_terminal_offset_vs_moon_kms,
                "notes": c.notes,
            }
            for c in verdict.per_cycle
        ],
        "notes": verdict.notes,
    }


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#332] V4 HFEM Uranus phase 4 -- sha={sha}", flush=True)
    print(f"[#332] candidate = {SILVER_ID}", flush=True)
    print(f"[#332] sequence = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#332] V_inf = {SILVER_VINF}, ToF (days) = {SILVER_TOF}, "
        f"rel_off = {SILVER_REL_OFF_DEG}deg, phase0 = {SILVER_PHASE0_DEG}deg",
        flush=True,
    )
    print(
        f"[#332] V4 model: scipy DOP853 + Uranus J2 ({URANUS_J2:.5e}, "
        f"R_eq={URANUS_R_EQ_KM:.0f} km) + classical-moon third-body "
        f"({', '.join(URANIAN_PERTURBER_MOONS)})",
        flush=True,
    )
    print(
        "[#332] GMAT-fallback rationale: Uranian satellite SPICE kernels "
        "(URA111/URA107) unavailable in GMAT R2022a install. Phase 4.1 "
        "(full GMAT V4) remains a downstream catalogue-admission gate.",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#332 V4 HFEM Uranus phase 4 - Umbriel-Oberon-Umbriel SILVER",
            "candidate_id": SILVER_ID,
            "source_jsonl": str(SILVER_VERIFIED_JSONL.relative_to(ROOT)),
            "v3_verdicts_jsonl": str(V3_VERDICTS_JSONL.relative_to(ROOT)),
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg": SILVER_REL_OFF_DEG,
                "phase0_deg": SILVER_PHASE0_DEG,
                "n_rev": list(SILVER_NREV),
            },
            "v4_model": {
                "integrator": "scipy DOP853",
                "uranus_j2": URANUS_J2,
                "uranus_r_eq_km": URANUS_R_EQ_KM,
                "perturber_moons": list(URANIAN_PERTURBER_MOONS),
                "j2_source": "Jacobson 2014, AJ 148:76, Table 4",
                "spice_status": "Uranian satellite SPICE kernels (URA111/URA107) NOT installed",
                "fallback_path": (
                    "documented scipy-fallback per task #332 spec; "
                    "Phase 4.1 full-GMAT V4 remains pending kernel install"
                ),
            },
            "driver_floors": {
                "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "n_cycles_min": V4_N_CYCLES_MIN,
            },
            "git_sha": sha,
        }
    )

    n_cycles_grid = (3, 5, 10)
    headline_v4: dict[int, V4UranusVerdict] = {}
    headline_v3: dict[int, V3Verdict3D] = {}
    for nc in n_cycles_grid:
        print(f"[#332] run_v2_moontour(n_cycles={nc})...", flush=True)
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
            notes=f"#332 phase-4 V2 input, n_cycles={nc}",
        )
        print(f"[#332] run_v3_3d(n_cycles={nc})...", flush=True)
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
            notes=f"#332 phase-4 V3 input, n_cycles={nc}",
        )
        headline_v3[nc] = v3
        print(f"[#332] run_v4_uranus(n_cycles={nc})...", flush=True)
        t_v4 = time.time()
        v4 = run_v4_uranus(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            None,
            v3_verdict=v3,
            n_cycles=nc,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0_DEG,
            notes=f"#332 phase-4 V4 scan, n_cycles={nc}",
        )
        headline_v4[nc] = v4
        print(
            f"  -> passes_v4={v4.passes_v4} | "
            f"bounded={v4.bounded_drift_survives} | "
            f"completed={v4.n_cycles_propagated}/{nc} | "
            f"agreement_kms={v4.drift_agreement_kms:.3e} | "
            f"elapsed={time.time() - t_v4:.1f}s",
            flush=True,
        )
        for c in v4.per_cycle:
            print(
                f"     cycle {c.cycle_index}: V4_drift={c.rendezvous_drift_kms_v4:.3e} km, "
                f"V3_drift={c.rendezvous_drift_kms_v3:.3e} km, "
                f"|V4-V3|={c.agreement_kms:.3e} km, "
                f"V4-vs-moon={c.v4_terminal_offset_vs_moon_kms:.3e} km",
                flush=True,
            )
        rows.append(_v4_verdict_to_row(v4))

    # V4 headline + catalogue-admission recommendation.
    v4_3 = headline_v4[3]
    v4_5 = headline_v4[5]
    v4_10 = headline_v4[10]
    all_pass = v4_3.passes_v4 and v4_5.passes_v4 and v4_10.passes_v4
    any_fail = (not v4_3.passes_v4) or (not v4_5.passes_v4) or (not v4_10.passes_v4)
    all_bounded = (
        v4_3.bounded_drift_survives and v4_5.bounded_drift_survives and v4_10.bounded_drift_survives
    )

    if all_pass:
        verdict_label = "PASS"
        next_step = (
            "9 of 10 catalogue-admission gates cleared (closure, DOP853 "
            "cross-check, physical-sanity, lit-novelty, ML flagger, V1 3D, "
            "V2 moontour, V3 IAS15, V4 scipy-fallback Uranus J2 + n-body). "
            "REMAINING: (1) #329 Heaton-Longuski 2003 JSR paywall human "
            "access; (2) Phase 4.1 -- full GMAT V4 once Uranian satellite "
            "SPICE kernels (URA111/URA107) are installed -- strictly the "
            "real-eph gate the task #332 originally specified. Catalogue "
            "admission as quasi_cycler row is gated on those two."
        )
    elif any_fail and not all_bounded:
        verdict_label = "FAIL"
        next_step = (
            "Bounded-drift signature did NOT survive V4-class perturbations "
            "(Uranus J2 + classical-moon third-body). Retire to negative-"
            "results registry (#172) recording which order of perturbation "
            "breaks the signature. The V3-confirmed bounded-drift was a "
            "Keplerian-coplanar model artifact under real Uranian dynamics."
        )
    elif any_fail and all_bounded:
        verdict_label = "MARGINAL"
        next_step = (
            "Bounded-drift survives qualitatively but the V4-vs-V3 agreement "
            "exceeds the 50,000 km floor at one or more n_cycles. The "
            "signature is robust to the dominant V4-class perturbations but "
            "the perturbations meaningfully shift the per-cycle drift. "
            "Phase 4.1 (full GMAT V4 with Uranian satellite SPICE) is the "
            "correct next gate; the marginal result is recorded but does "
            "not admit to catalogue."
        )
    else:  # pragma: no cover - defensive
        verdict_label = "MIXED"
        next_step = "Investigate per-cycle disagreement before catalogue admission."

    rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "candidate_id": SILVER_ID,
            "v3_n3_max_drift_kms": max(headline_v3[3].per_cycle_drift_kms_v3),
            "v3_n5_max_drift_kms": max(headline_v3[5].per_cycle_drift_kms_v3),
            "v3_n10_max_drift_kms": max(headline_v3[10].per_cycle_drift_kms_v3),
            "v4_n_cycles_3": {
                "passes_v4": v4_3.passes_v4,
                "bounded_drift_survives": v4_3.bounded_drift_survives,
                "agreement_kms": v4_3.drift_agreement_kms,
                "n_completed": v4_3.n_cycles_propagated,
                "v4_max_drift_kms": (
                    max(v4_3.per_cycle_drift_kms_v4) if v4_3.per_cycle_drift_kms_v4 else None
                ),
            },
            "v4_n_cycles_5": {
                "passes_v4": v4_5.passes_v4,
                "bounded_drift_survives": v4_5.bounded_drift_survives,
                "agreement_kms": v4_5.drift_agreement_kms,
                "n_completed": v4_5.n_cycles_propagated,
                "v4_max_drift_kms": (
                    max(v4_5.per_cycle_drift_kms_v4) if v4_5.per_cycle_drift_kms_v4 else None
                ),
            },
            "v4_n_cycles_10": {
                "passes_v4": v4_10.passes_v4,
                "bounded_drift_survives": v4_10.bounded_drift_survives,
                "agreement_kms": v4_10.drift_agreement_kms,
                "n_completed": v4_10.n_cycles_propagated,
                "v4_max_drift_kms": (
                    max(v4_10.per_cycle_drift_kms_v4) if v4_10.per_cycle_drift_kms_v4 else None
                ),
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
    print(f"[#332] wrote {OUT_JSONL}", flush=True)
    print(f"[#332] verdict: {verdict_label} (elapsed {time.time() - t0:.1f}s)", flush=True)
    print(f"[#332] next step: {next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
