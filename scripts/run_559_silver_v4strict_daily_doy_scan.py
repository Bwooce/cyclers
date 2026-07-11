"""#559 — Daily launch-epoch DOY-sensitivity scan on the #327 SILVER under V4-strict.

Phase 2 follow-up to #338 (which found 2000-01-15 FAIL / 2000-06-21 PASS, same year).
That annual-resolution sweep confirmed sub-year DOY sensitivity but couldn't resolve
whether June 21 sits in a WIDE tolerant band or on a KNIFE-EDGE boundary.

This script samples EVERY calendar day across TWO specific year-long windows:

* 2000-01-01 through 2000-12-31 (366 days, leap year): brackets both the known
  FAIL point (2000-01-15) and the known PASS point (2000-06-21), allowing
  within-year boundary characterization.
* 2030-01-01 through 2030-12-31 (365 days): independent check on whether the
  PASS band repeats similarly in a different decade (decade-scale robustness).

Total ~731 epochs at ~1 sec each (~10-15 min wall-clock). Each runs V4-strict
at ``n_cycles=3`` (the minimum required by the driver).

Output: one row per daily epoch, showing launch_epoch, passes_v4_strict verdict,
and drift_agreement_kms_vs_v3 for direct visual/statistical boundary analysis.

NO catalogue writeback. READ-ONLY validation reporting.

Discipline anchors
------------------
* READ-ONLY on ``src/cyclerfinder/data/validation/v[1-4]_*.py`` modules.
* READ-ONLY on the existing #335/#338 outputs.
* NO catalogue writeback.
* NO `--no-verify`.

Run as::

    uv run python scripts/run_559_silver_v4strict_daily_doy_scan.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.data.validation.v2_moontour import run_v2_moontour  # noqa: E402
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d  # noqa: E402
from cyclerfinder.data.validation.v4_uranus import (  # noqa: E402
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4UranusVerdict,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (  # noqa: E402
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    run_v4_uranus_strict,
)

# Stored SILVER row fields — same constants as scripts/run_335_*.py
# (sourced from data/silver_327_verified.jsonl, READ-ONLY here).
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


# Daily sweep: generate every calendar day across TWO specific year-long windows
# to resolve whether the June 21 PASS from #338 sits in a wide band or on a knife-edge.
def _generate_daily_epochs() -> tuple[str, ...]:
    """Generate daily epochs for 2000-01-01 to 2000-12-31 and 2030-01-01 to 2030-12-31."""
    epochs = []
    for year in (2000, 2030):
        start = datetime(year, 1, 1)
        # 2000 is leap year (366 days), 2030 is not (365 days)
        end = datetime(2001, 1, 1) if year == 2000 else datetime(2031, 1, 1)
        current = start
        while current < end:
            epochs.append(current.strftime("%Y-%m-%dT00:00:00"))
            current += timedelta(days=1)
    return tuple(epochs)


EPOCHS = _generate_daily_epochs()

OUT_JSONL = ROOT / "data" / "silver_327_v4_strict_daily_sweep_559.jsonl"

_REGION_ID = "silver-327-umbriel-oberon-daily-doy-sensitivity-2026-07-11"
_METHOD = MethodCapability(
    genome=(
        "Daily-resolution launch-epoch sweep of catalogue row #312 "
        "(silver-327 Umbriel-Oberon-Umbriel) under the existing, frozen V4-strict "
        "validation pipeline -- read-only validation reporting, no genome/corrector change"
    ),
    corrector="existing V2->V3->V4-scipy->V4-strict chain (validation/v[1-4]_*.py), unmodified",
    capability_tags=frozenset(
        {"cr3bp", "real-ephemeris", "v4-strict", "doy-sensitivity", "validation-reporting"}
    ),
    git_sha="working-tree",
)

# V4-strict requires n_cycles >= V4_N_CYCLES_MIN (== 3). Stay at the floor
# for the sweep -- runtime is then ~1 s/epoch x 100 epochs ~= 100 s, well
# inside the 2-minute compute budget from the task brief.
N_CYCLES = V4_N_CYCLES_MIN

REPORT_EVERY = 50  # progress print cadence (731 epochs total, ~14x the annual sweep)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _build_v3_v4scipy(n_cycles: int) -> tuple[V3Verdict3D, V4UranusVerdict]:
    """V2 -> V3 -> V4-scipy chain.

    Epoch-blind (V2/V3/V4-scipy use circular-coplanar Kepler moons), so the
    chain is built ONCE and reused across all ~730 V4-strict daily epochs.
    """
    v2 = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#559 daily sweep V4-strict input chain, n_cycles={n_cycles}",
    )
    v3 = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v2_verdict=v2,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#559 daily sweep V4-strict input chain, n_cycles={n_cycles}",
    )
    v4 = run_v4_uranus(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v3_verdict=v3,
        n_cycles=n_cycles,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
        notes=f"#559 daily sweep V4-strict input chain, n_cycles={n_cycles}",
    )
    return v3, v4


def main() -> int:
    preflight_search(
        task_no=559,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(EPOCHS),
        override_reason=(
            "read-only validation-reporting sweep on an already-catalogued row (#312) "
            "using the existing frozen V4-strict pipeline unmodified -- not a discovery "
            "sweep; a timing pilot doesn't apply to a fixed, already-scoped 731-epoch "
            "daily scan bounded by the calendar (two named year-long windows)."
        ),
    )
    sha = _git_sha()
    t0 = time.time()
    print(f"[#559] daily DOY-sensitivity epoch sweep on #327 SILVER -- sha={sha}", flush=True)
    print(f"[#559] candidate = {SILVER_ID}", flush=True)
    print(f"[#559] sequence  = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#559] SPICE kernels = LSK + PCK + ura111.bsp ({DEFAULT_URA_PATH})",
        flush=True,
    )
    print(
        f"[#559] sweep = {len(EPOCHS)} daily epochs "
        f"(2000-01-01..2000-12-31 and 2030-01-01..2030-12-31), "
        f"n_cycles={N_CYCLES}",
        flush=True,
    )

    # Sanity-check kernels are present BEFORE running the long chain.
    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#559] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    # Build the V3 + V4-scipy chain ONCE (epoch-blind).
    print(f"[#559] V2->V3->V4-scipy chain at n_cycles={N_CYCLES}...", flush=True)
    t_chain = time.time()
    v3, v4_scipy = _build_v3_v4scipy(N_CYCLES)
    print(f"[#559]   chain ready (elapsed {time.time() - t_chain:.1f}s)", flush=True)

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": (
                "#559 daily DOY-sensitivity launch-epoch sweep on the #327 SILVER under V4-strict"
            ),
            "successor_to": "#338 annual sweep (found 2000-01-15 FAIL / 2000-06-21 PASS)",
            "candidate_id": SILVER_ID,
            "stored_silver": {
                "sequence": list(SILVER_SEQ),
                "vinf_per_encounter_kms": list(SILVER_VINF),
                "tof_days": list(SILVER_TOF),
                "rel_offset_deg": SILVER_REL_OFF_DEG,
                "phase0_deg": SILVER_PHASE0_DEG,
                "n_rev": list(SILVER_NREV),
            },
            "v4_strict_model": {
                "integrator": "scipy DOP853",
                "uranus_j2": URANUS_J2,
                "uranus_r_eq_km": URANUS_R_EQ_KM,
                "perturber_moons": list(URANIAN_PERTURBER_MOONS),
                "spice_kernel_source": (
                    "JPL/NAIF generic_kernels/spk/satellites/a_old_versions/ura111.bsp; "
                    "5 classical Uranian moons + Uranus + Earth + Sun, 1900-2099 ET"
                ),
            },
            "driver_floors": {
                "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "n_cycles_min": V4_N_CYCLES_MIN,
            },
            "sweep_params": {
                "windows": [
                    "2000-01-01 to 2000-12-31 (366 days, leap year)",
                    "2030-01-01 to 2030-12-31 (365 days)",
                ],
                "n_epochs": len(EPOCHS),
                "n_cycles_per_epoch": N_CYCLES,
                "sampling_note": (
                    "Daily resolution across two year-long windows to resolve whether "
                    "the June-21 PASS from #338 sits in a wide tolerant band or on a "
                    "knife-edge boundary. Windows bracket known FAIL (2000-01-15) and "
                    "PASS (2000-06-21 / 2030-06-21) reference points."
                ),
            },
            "git_sha": sha,
        }
    )

    n_pass = 0
    n_fail = 0
    n_completed_total = 0
    for i, epoch in enumerate(EPOCHS):
        t_run = time.time()
        try:
            v4s = run_v4_uranus_strict(
                SILVER_ID,
                SILVER_SEQ,
                SILVER_VINF,
                SILVER_TOF,
                SILVER_REL_OFF_DEG,
                epoch,
                None,
                v3_verdict=v3,
                v4_scipy_verdict=v4_scipy,
                n_cycles=N_CYCLES,
                n_revs=SILVER_NREV,
                notes=f"#559 daily sweep, epoch={epoch}, n_cycles={N_CYCLES}",
            )
            row = {
                "kind": "daily_sweep_row",
                "launch_epoch_utc": epoch,
                "passes_v4_strict": bool(v4s.passes_v4_strict),
                "bounded_drift_survives": bool(v4s.bounded_drift_survives),
                "n_cycles_propagated": int(v4s.n_cycles_propagated),
                "n_cycles_requested": N_CYCLES,
                "drift_agreement_kms_vs_v3": float(v4s.drift_agreement_kms_vs_v3),
                "drift_agreement_kms_vs_v4_scipy": float(v4s.drift_agreement_kms_vs_v4_scipy),
                "per_cycle_drift_kms_v4_strict": list(v4s.per_cycle_drift_kms_v4_strict),
                "per_cycle_drift_kms_v4_scipy": list(v4s.per_cycle_drift_kms_v4_scipy),
                "per_cycle_drift_kms_v3": list(v4s.per_cycle_drift_kms_v3),
                "eccentricity_used_e_umbriel": float(v4s.eccentricity_used_e_umbriel),
                "eccentricity_used_e_oberon": float(v4s.eccentricity_used_e_oberon),
                "inclination_used_deg_umbriel": float(v4s.inclination_used_deg_umbriel),
                "inclination_used_deg_oberon": float(v4s.inclination_used_deg_oberon),
                "wall_clock_s": float(time.time() - t_run),
            }
            if v4s.passes_v4_strict:
                n_pass += 1
            else:
                n_fail += 1
            n_completed_total += v4s.n_cycles_propagated
        except Exception as exc:
            row = {
                "kind": "daily_sweep_row",
                "launch_epoch_utc": epoch,
                "passes_v4_strict": False,
                "bounded_drift_survives": False,
                "n_cycles_propagated": 0,
                "n_cycles_requested": N_CYCLES,
                "drift_agreement_kms_vs_v3": float("inf"),
                "drift_agreement_kms_vs_v4_scipy": float("inf"),
                "error": f"{type(exc).__name__}: {exc}",
                "wall_clock_s": float(time.time() - t_run),
            }
            n_fail += 1

        rows.append(row)

        if (i + 1) % REPORT_EVERY == 0 or i == len(EPOCHS) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta_s = (len(EPOCHS) - (i + 1)) / rate if rate > 0 else 0.0
            print(
                f"[#559]   epoch {epoch} ({i + 1:3d}/{len(EPOCHS)}): "
                f"PASS={n_pass} FAIL={n_fail} | "
                f"elapsed={elapsed:.1f}s ETA={eta_s:.1f}s | "
                f"last drift_vs_v3={row.get('drift_agreement_kms_vs_v3', float('nan')):.3e} km",
                flush=True,
            )

    # Summary header (kept as a separate _meta row so the JSONL can be
    # consumed both row-by-row and headline-only).
    rows.append(
        {
            "_meta": True,
            "kind": "daily_sweep_summary",
            "candidate_id": SILVER_ID,
            "n_epochs_total": len(EPOCHS),
            "n_pass": n_pass,
            "n_fail": n_fail,
            "pass_fraction": n_pass / len(EPOCHS) if EPOCHS else 0.0,
            "n_cycles_completed_sum": n_completed_total,
            "n_cycles_completed_expected": len(EPOCHS) * N_CYCLES,
            "verdict_label": (
                "ALL_PASS" if n_fail == 0 else "ALL_FAIL" if n_pass == 0 else "MIXED"
            ),
            "next_step": (
                "analyze pass/fail pattern by day-of-year across both 2000 and 2030 windows"
            ),
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    print(f"\n[#559] wrote {OUT_JSONL}", flush=True)
    print(
        f"[#559] sweep complete: PASS={n_pass} FAIL={n_fail} (elapsed {time.time() - t0:.1f}s)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
