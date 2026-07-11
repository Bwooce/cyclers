"""#338 — Annual launch-epoch sweep 2000-2099 on the #327 SILVER under V4-strict.

Successor to #335 (commit ``011a289``), which ran V4-strict at 3 launch
epochs and produced a MIXED / EPOCH_DEPENDENT verdict:

* 2000-01-15 → FAIL (agreement-vs-V3 ≈ 90,620 km, exceeds 50,000 km floor)
* 2030-06-21 → PASS (agreement-vs-V3 ≈ 9,723 km)
* 2050-12-31 → PASS (agreement-vs-V3 ≈ 3,680 km)

3 epochs is too coarse to characterise the PASS/FAIL boundary.

This script samples ONE epoch per year — Y-06-21T00:00:00 (mid-range,
matches the 2030-06-21 anchor that PASSED) — over the full 2000-2099
window (100 epochs). Each epoch runs V4-strict at ``n_cycles=3``
(the minimum required by the driver, ~1 sec each → ~2 min wall-clock).

The output is a row per year, suitable for the boundary-analysis script
(``scripts/analyze_338_boundary.py``) to look for periodicity in
``passes_v4_strict`` and ``agreement_vs_v3`` vs year.

Honest scope
------------
Annual sampling cannot resolve sub-year periodicities (e.g. the
Umbriel-Oberon synodic period ~5.99 days, or Umbriel's orbital period
~4.14 days). What it CAN resolve:

  * decade-scale launch-epoch sensitivity
  * Uranus's heliocentric orbital phase (period ~84 yr)
  * aliasing into the annual grid is detectable as low-frequency beats
    against a known sub-year period

If the boundary analysis sees a CYCLIC PASS pattern with period close
to 84 yr (Uranus orbital) or fractions thereof, that's a real physical
signal. If the pattern looks irregular at annual resolution, future work
needs daily/weekly sampling near a transition (Part C of the doc).

NO catalogue writeback. The verdict is recorded; what to do with it is
the boundary analysis (#338 Part B) + doc (Part C).

Discipline anchors
------------------
* READ-ONLY on ``src/cyclerfinder/data/validation/v[1-4]_*.py`` modules.
* READ-ONLY on the existing #335 JSONL outputs.
* NO catalogue writeback.
* NO `--no-verify`.

Run as::

    uv run python scripts/run_338_silver_v4strict_annual_sweep.py
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

# Annual sweep: Y-06-21 (northern summer solstice; matches the #335 2030-06-21
# anchor that PASSED, so we're scanning around a known good operating point
# rather than the J2000-anchored FAIL).
YEARS = tuple(range(2000, 2100))
DOY_LABEL = "06-21T00:00:00"

OUT_JSONL = ROOT / "data" / "silver_327_v4_strict_annual_sweep_338.jsonl"

# V4-strict requires n_cycles >= V4_N_CYCLES_MIN (== 3). Stay at the floor
# for the sweep -- runtime is then ~1 s/epoch x 100 epochs ~= 100 s, well
# inside the 2-minute compute budget from the task brief.
N_CYCLES = V4_N_CYCLES_MIN

REPORT_EVERY = 10  # progress print cadence


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
    chain is built ONCE and reused across all 100 V4-strict epochs.
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
        notes=f"#338 annual sweep V4-strict input chain, n_cycles={n_cycles}",
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
        notes=f"#338 annual sweep V4-strict input chain, n_cycles={n_cycles}",
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
        notes=f"#338 annual sweep V4-strict input chain, n_cycles={n_cycles}",
    )
    return v3, v4


def main() -> int:
    sha = _git_sha()
    t0 = time.time()
    print(f"[#338] annual epoch sweep on #327 SILVER -- sha={sha}", flush=True)
    print(f"[#338] candidate = {SILVER_ID}", flush=True)
    print(f"[#338] sequence  = {SILVER_SEQ}, n_rev = {SILVER_NREV}", flush=True)
    print(
        f"[#338] SPICE kernels = LSK + PCK + ura111.bsp ({DEFAULT_URA_PATH})",
        flush=True,
    )
    print(
        f"[#338] sweep = {len(YEARS)} epochs (Y-{DOY_LABEL}, Y in {YEARS[0]}..{YEARS[-1]}), "
        f"n_cycles={N_CYCLES}",
        flush=True,
    )

    # Sanity-check kernels are present BEFORE running the long chain.
    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#338] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    # Build the V3 + V4-scipy chain ONCE (epoch-blind).
    print(f"[#338] V2->V3->V4-scipy chain at n_cycles={N_CYCLES}...", flush=True)
    t_chain = time.time()
    v3, v4_scipy = _build_v3_v4scipy(N_CYCLES)
    print(f"[#338]   chain ready (elapsed {time.time() - t_chain:.1f}s)", flush=True)

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#338 annual launch-epoch sweep 2000-2099 on the #327 SILVER under V4-strict",
            "successor_to": "#335 (MIXED / EPOCH_DEPENDENT verdict at 3 epochs)",
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
                "years": [YEARS[0], YEARS[-1]],
                "doy_label": DOY_LABEL,
                "n_epochs": len(YEARS),
                "n_cycles_per_epoch": N_CYCLES,
                "sampling_caveat": (
                    "Annual sampling cannot resolve sub-year periodicities; "
                    "if the boundary analysis suggests periods << 1 yr "
                    "(e.g. Umbriel-Oberon synodic ~5.99 days), the annual "
                    "grid is aliased and finer sampling is needed."
                ),
            },
            "git_sha": sha,
        }
    )

    n_pass = 0
    n_fail = 0
    n_completed_total = 0
    for i, year in enumerate(YEARS):
        epoch = f"{year:04d}-{DOY_LABEL}"
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
                notes=f"#338 annual sweep, epoch={epoch}, n_cycles={N_CYCLES}",
            )
            row = {
                "kind": "annual_sweep_row",
                "year": year,
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
                "eccentricity_used_e_umbriel": float(v4s.eccentricity_used_e_body1),
                "eccentricity_used_e_oberon": float(v4s.eccentricity_used_e_body2),
                "inclination_used_deg_umbriel": float(v4s.inclination_used_deg_body1),
                "inclination_used_deg_oberon": float(v4s.inclination_used_deg_body2),
                "wall_clock_s": float(time.time() - t_run),
            }
            if v4s.passes_v4_strict:
                n_pass += 1
            else:
                n_fail += 1
            n_completed_total += v4s.n_cycles_propagated
        except Exception as exc:
            row = {
                "kind": "annual_sweep_row",
                "year": year,
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

        if (i + 1) % REPORT_EVERY == 0 or i == len(YEARS) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta_s = (len(YEARS) - (i + 1)) / rate if rate > 0 else 0.0
            print(
                f"[#338]   year {year} ({i + 1:3d}/{len(YEARS)}): "
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
            "kind": "annual_sweep_summary",
            "candidate_id": SILVER_ID,
            "n_epochs_total": len(YEARS),
            "n_pass": n_pass,
            "n_fail": n_fail,
            "pass_fraction": n_pass / len(YEARS) if YEARS else 0.0,
            "n_cycles_completed_sum": n_completed_total,
            "n_cycles_completed_expected": len(YEARS) * N_CYCLES,
            "verdict_label": (
                "ALL_PASS" if n_fail == 0 else "ALL_FAIL" if n_pass == 0 else "MIXED"
            ),
            "next_step": "scripts/analyze_338_boundary.py for the boundary verdict",
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")

    print(f"\n[#338] wrote {OUT_JSONL}", flush=True)
    print(
        f"[#338] sweep complete: PASS={n_pass} FAIL={n_fail} (elapsed {time.time() - t0:.1f}s)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
