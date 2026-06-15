"""Galileo VEEGA reproduction probe (#298 Phase 2 of #289).

Sweeps the Tisserand-Poincaré MGA enumerator over the published Galileo
launch window (1989-10-01..1989-11-15) with planet_set=(V,E,J), max_legs=4,
and writes each (V,E,E,J) candidate + its Phase-1 closure to
``data/scan_298_galileo_veega.jsonl``.

This is a reproduction PROBE, not a closure-grade reproduction. Galileo's
real flight walks across multiple V_inf shells (4 km/s at launch, ~9 km/s
by Earth-2, ~5.6 km/s at Jupiter); the Tisserand-Poincaré graph is a
single-V_inf-shell pre-screen. The probe confirms the *structural*
reproduction (sequence + launch window) and records the residual gap that
Phase 3 (DSM + per-leg TOF optimisation) must close.

See ``docs/notes/2026-06-16-298-289-phase2-tisserand-enumerator.md`` for
the full phase plan and the interpretation of the JSONL fields.

Sourcing: Diehl-Belbruno-Roberts 1986 AAS (Galileo trajectory) — KNOWN_CORPUS
anchor at corpus rev ``568d8a4``. Launch UTC: 1989-10-18T16:53:40
(public Galileo trajectory record); first body encounter (Venus): ~1990-02-10.

Run with:
    uv run python scripts/probe_298_galileo_veega.py
"""

from __future__ import annotations

import json
import time
from datetime import UTC, datetime
from pathlib import Path

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.tisserand_mga_window import find_mga_chains, validate_chain_candidate

OUTPUT_PATH = Path(__file__).resolve().parents[1] / "data" / "scan_298_galileo_veega.jsonl"

LAUNCH_WINDOW = ("1989-10-01T00:00:00", "1989-11-15T00:00:00")
PLANET_SET = ("V", "E", "J")
VINF_GRID_KMS = (8.0, 9.0, 10.0, 11.0)
TOF_BOX_DAYS = (60.0, 1200.0)
MAX_LEGS = 4
EPOCH_STEP_DAYS = 15.0

PUBLISHED_LAUNCH_UTC = "1989-10-18T16:53:40"  # Diehl-Belbruno-Roberts 1986


def main() -> None:
    """Enumerate, validate, and write JSONL."""
    eph = Ephemeris("astropy")
    print(f"[{datetime.now(UTC).isoformat()}] starting probe")
    t0 = time.time()

    n_total = 0
    n_veej = 0
    n_closed = 0
    rows: list[dict] = []

    for cand in find_mga_chains(
        launch_window=LAUNCH_WINDOW,
        planet_set=PLANET_SET,
        max_legs=MAX_LEGS,
        vinf_grid_kms=VINF_GRID_KMS,
        tof_box_days_per_leg=TOF_BOX_DAYS,
        epoch_step_days=EPOCH_STEP_DAYS,
    ):
        n_total += 1
        if cand.sequence != ("V", "E", "E", "J"):
            continue
        n_veej += 1
        # Try a loose closure tolerance — we already know the single-shell
        # seed will NOT meet the 0.5 km/s gate (Galileo is multi-shell).
        # Use 500 km/s ceilings so we always record the residual; the per-leg
        # V_inf disagreement and flyby continuity ΔV are the diagnostic, not
        # the converged/not-converged flag.
        closure = validate_chain_candidate(
            cand,
            eph,
            closure_tol_kms=500.0,
            flyby_continuity_tol_kms=500.0,
            independent_cross_check=False,
        )
        row = {
            "sequence": list(cand.sequence),
            "vinf_tuple_kms": list(cand.vinf_tuple_kms),
            "leg_tofs_days": list(cand.leg_tofs_days),
            "launch_epoch_utc": cand.launch_epoch_utc,
            "tisserand_parameter": cand.tisserand_parameter,
            "chain_score": cand.chain_score,
            "published_launch_utc": PUBLISHED_LAUNCH_UTC,
            "closure_residual_kms": (None if closure is None else closure.closure_residual_kms),
            "flyby_continuity_max_dv_kms": (
                None if closure is None else closure.flyby_continuity_max_dv_kms
            ),
            "per_encounter_vinf_kms": (
                None if closure is None else list(closure.per_encounter_vinf_kms)
            ),
            "diagnostic_recorded": closure is not None,
            "scan_id": "298",
            "scan_method_version": "tisserand_poincare_mga_window_v1",
            "scan_timestamp_utc": datetime.now(UTC).isoformat(),
            "notes": (
                "Phase-2 single-V_inf-shell pre-screen seed. Galileo VEEGA "
                "is structurally a multi-shell pump tour — residual far "
                "above the 0.5 km/s gate is EXPECTED. Phase 3 (DSM + "
                "per-leg TOF optimisation) must close the gap."
            ),
        }
        rows.append(row)
        if closure is not None:
            n_closed += 1

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")

    t1 = time.time()
    print(
        f"[{datetime.now(UTC).isoformat()}] done in {t1 - t0:.1f}s; "
        f"enumerated={n_total} VEEJ={n_veej} diagnostic_recorded={n_closed}"
    )
    print(f"Output: {OUTPUT_PATH}")

    # Honest residual summary.
    valid = [r for r in rows if r["closure_residual_kms"] is not None]
    if valid:
        residuals = [r["closure_residual_kms"] for r in valid]
        print(f"Closure residuals (km/s): min={min(residuals):.2f} max={max(residuals):.2f}")
        print("  --> single-shell seed gap; Phase 3 DSM + TOF-optimisation required")


if __name__ == "__main__":
    main()
