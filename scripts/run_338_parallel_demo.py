"""#321 Part C -- parallel demo of the #338 annual V4-strict sweep.

Side-by-side wall-clock comparison: the original ``scripts/run_338_silver_v4strict_annual_sweep.py``
runs 100 epochs serially through V4-strict at ``n_cycles=3`` (~1 s/epoch),
total ~100 s. This script reproduces the same 100 V4-strict runs through
the new ``cyclerfinder.parallel.parallel_sweep`` substrate and records both
wall times for the docs note.

The closure semantics MUST match the serial baseline byte-for-byte: the
``year -> dict`` mapping produced by each cell here is identical in shape
and field semantics to the per-year row that ``run_338_silver_v4strict_annual_sweep.py``
writes to ``data/silver_327_v4_strict_annual_sweep_338.jsonl``. The only
difference is the order in which rows are produced (parallel reorders by
worker completion; we restore input order before writing).

NO catalogue writeback. This is a demo + speedup measurement only.

The output goes to a DIFFERENT path
(``data/silver_327_v4_strict_annual_sweep_321_parallel_demo.jsonl``) so the
serial baseline's output is not overwritten. The script ALSO writes a small
``data/scan_321_parallel_demo_summary.jsonl`` describing the wall-clock
comparison.

Discipline anchors
------------------
* READ-ONLY on ``run_338_silver_v4strict_annual_sweep.py`` and on the v3/v4
  validation modules; we only import the published surfaces.
* NO catalogue writeback.
* NO `--no-verify`.

Run as::

    uv run python scripts/run_338_parallel_demo.py [--n-workers 4]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.validation.v2_moontour import run_v2_moontour  # noqa: E402
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d  # noqa: E402
from cyclerfinder.data.validation.v4_uranus import (  # noqa: E402
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
from cyclerfinder.parallel import (  # noqa: E402
    ParallelSweepConfig,
    parallel_sweep,
)

# Constants matched against scripts/run_338_silver_v4strict_annual_sweep.py.
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

YEARS = tuple(range(2000, 2100))
DOY_LABEL = "06-21T00:00:00"
N_CYCLES = V4_N_CYCLES_MIN  # == 3

OUT_JSONL = ROOT / "data" / "silver_327_v4_strict_annual_sweep_321_parallel_demo.jsonl"
SUMMARY_JSONL = ROOT / "data" / "scan_321_parallel_demo_summary.jsonl"


# ---------------------------------------------------------------------------
# Cell payload + closure (top-level for pickle safety)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Cell:
    """Per-epoch cell payload.

    The V3 + V4-scipy verdicts are epoch-blind in this sweep (V3/V4-scipy
    use circular-coplanar Kepler moons) so they are computed ONCE outside
    the closure and passed into every cell. They must be picklable; both
    dataclasses are frozen and contain only floats / tuples / lists, which
    pickle cleanly.
    """

    year: int
    epoch_utc: str
    v3: V3Verdict3D
    v4_scipy: V4UranusVerdict


def _run_one_epoch(cell: _Cell) -> dict[str, Any]:
    """Cell closure: run V4-strict at the given epoch and return a row dict.

    Matches the row schema produced by the serial baseline at
    ``scripts/run_338_silver_v4strict_annual_sweep.py`` (the post-`v4s` block).
    """
    t_run = time.time()
    try:
        v4s = run_v4_uranus_strict(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            cell.epoch_utc,
            None,
            v3_verdict=cell.v3,
            v4_scipy_verdict=cell.v4_scipy,
            n_cycles=N_CYCLES,
            n_revs=SILVER_NREV,
            notes=f"#321 parallel demo, epoch={cell.epoch_utc}, n_cycles={N_CYCLES}",
        )
        return {
            "kind": "annual_sweep_row",
            "year": cell.year,
            "launch_epoch_utc": cell.epoch_utc,
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
    except Exception as exc:
        return {
            "kind": "annual_sweep_row",
            "year": cell.year,
            "launch_epoch_utc": cell.epoch_utc,
            "passes_v4_strict": False,
            "bounded_drift_survives": False,
            "n_cycles_propagated": 0,
            "n_cycles_requested": N_CYCLES,
            "drift_agreement_kms_vs_v3": float("inf"),
            "drift_agreement_kms_vs_v4_scipy": float("inf"),
            "error": f"{type(exc).__name__}: {exc}",
            "wall_clock_s": float(time.time() - t_run),
        }


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _build_v3_v4scipy(n_cycles: int) -> tuple[V3Verdict3D, V4UranusVerdict]:
    """V2 -> V3 -> V4-scipy chain (epoch-blind). Built once."""
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
        notes=f"#321 parallel demo V4-strict input chain, n_cycles={n_cycles}",
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
        notes=f"#321 parallel demo V4-strict input chain, n_cycles={n_cycles}",
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
        notes=f"#321 parallel demo V4-strict input chain, n_cycles={n_cycles}",
    )
    return v3, v4


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--n-workers",
        type=int,
        default=-1,
        help="parallel workers (-1 = all cores; default: -1)",
    )
    ap.add_argument(
        "--skip-serial",
        action="store_true",
        help="skip the serial timing pass (use historical 100s estimate)",
    )
    args = ap.parse_args()

    sha = _git_sha()
    print(f"[#321] parallel-demo on #338 annual sweep -- sha={sha}", flush=True)
    print(f"[#321] n_workers={args.n_workers}, n_epochs={len(YEARS)}", flush=True)

    # Pre-flight: kernels.
    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#321] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    # Build the V3+V4-scipy chain ONCE.
    print(f"[#321] V2->V3->V4-scipy chain at n_cycles={N_CYCLES}...", flush=True)
    t_chain = time.time()
    v3, v4_scipy = _build_v3_v4scipy(N_CYCLES)
    print(f"[#321]   chain ready ({time.time() - t_chain:.1f}s)", flush=True)

    # Assemble cells.
    cells = [
        _Cell(year=year, epoch_utc=f"{year:04d}-{DOY_LABEL}", v3=v3, v4_scipy=v4_scipy)
        for year in YEARS
    ]

    # ---------------------------------------------------------------------
    # Serial pass (baseline timing)
    # ---------------------------------------------------------------------
    if args.skip_serial:
        serial_wall = float("nan")
        serial_rows: list[dict[str, Any]] = []
        print("[#321] serial pass SKIPPED (use --skip-serial=False to re-time)", flush=True)
    else:
        print(f"[#321] serial pass over {len(cells)} cells...", flush=True)
        t_serial = time.time()
        serial_rows = [_run_one_epoch(cell) for cell in cells]
        serial_wall = time.time() - t_serial
        n_pass_serial = sum(1 for r in serial_rows if r.get("passes_v4_strict"))
        print(
            f"[#321]   serial: {serial_wall:.1f}s, PASS={n_pass_serial}/{len(serial_rows)}",
            flush=True,
        )

    # ---------------------------------------------------------------------
    # Parallel pass
    # ---------------------------------------------------------------------
    cfg = ParallelSweepConfig(
        n_workers=args.n_workers,
        backend="loky",
        verbose=0,
        raise_on_first_error=False,
    )
    print(f"[#321] parallel pass over {len(cells)} cells ({cfg.backend})...", flush=True)
    t_par = time.time()
    par_result = parallel_sweep(cells, _run_one_epoch, config=cfg)
    par_wall = time.time() - t_par
    par_rows = [
        r if r is not None else {"kind": "annual_sweep_row", "error": "cell None"}
        for r in par_result.results
    ]
    n_pass_par = sum(1 for r in par_rows if r.get("passes_v4_strict"))
    print(
        f"[#321]   parallel: {par_wall:.1f}s, PASS={n_pass_par}/{len(par_rows)}, "
        f"n_failed={par_result.n_failed}",
        flush=True,
    )

    # ---------------------------------------------------------------------
    # Equivalence check (parallel vs serial)
    # ---------------------------------------------------------------------
    equivalence_notes = "skipped (no serial baseline)"
    if not args.skip_serial:
        # Compare passes_v4_strict + drift_agreement_kms_vs_v3 byte-for-byte
        # (these are the deterministic publish-relevant fields).
        mismatches = 0
        for s, p in zip(serial_rows, par_rows, strict=True):
            if s.get("passes_v4_strict") != p.get("passes_v4_strict"):
                mismatches += 1
                continue
            sd = s.get("drift_agreement_kms_vs_v3", float("nan"))
            pd = p.get("drift_agreement_kms_vs_v3", float("nan"))
            # Exact equality of floats is the contract: SAME inputs, SAME
            # numerical kernels, SAME single-thread per-cell scipy. Allow
            # NaN==NaN as equal (both inf rows).
            if not (sd == pd or (sd != sd and pd != pd)):
                mismatches += 1
        equivalence_notes = (
            f"{len(serial_rows) - mismatches}/{len(serial_rows)} rows match "
            f"({mismatches} mismatches)"
        )
        print(f"[#321] equivalence: {equivalence_notes}", flush=True)

    speedup = (serial_wall / par_wall) if (serial_wall > 0 and par_wall > 0) else float("nan")
    print(
        f"[#321] speedup: {speedup:.2f}x (serial {serial_wall:.1f}s / parallel {par_wall:.1f}s)",
        flush=True,
    )

    # ---------------------------------------------------------------------
    # Write rows (parallel ordering preserved via parallel_sweep input order)
    # ---------------------------------------------------------------------
    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    meta: dict[str, Any] = {
        "_meta": True,
        "task": "#321 Part C parallel-sweep demo on #338 annual V4-strict sweep",
        "candidate_id": SILVER_ID,
        "n_epochs": len(YEARS),
        "n_cycles_per_epoch": N_CYCLES,
        "serial_wall_seconds": serial_wall,
        "parallel_wall_seconds": par_wall,
        "speedup_ratio": speedup,
        "n_workers_requested": args.n_workers,
        "n_workers_used_note": "joblib -1 = all logical cores",
        "equivalence_check": equivalence_notes,
        "driver_floors": {
            "agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
            "n_cycles_min": V4_N_CYCLES_MIN,
        },
        "git_sha": sha,
    }
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps(meta) + "\n")
        for row in par_rows:
            fh.write(json.dumps(row) + "\n")
    print(f"[#321] wrote {OUT_JSONL}", flush=True)

    # Summary file (small, machine-parseable for the docs note).
    SUMMARY_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY_JSONL.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#321 parallel-sweep demo summary",
                    "git_sha": sha,
                    "n_epochs": len(YEARS),
                    "n_cycles_per_epoch": N_CYCLES,
                    "serial_wall_seconds": serial_wall,
                    "parallel_wall_seconds": par_wall,
                    "speedup_ratio": speedup,
                    "n_workers_requested": args.n_workers,
                    "parallel_succeeded": par_result.n_succeeded,
                    "parallel_failed": par_result.n_failed,
                    "equivalence_notes": equivalence_notes,
                }
            )
            + "\n"
        )
    print(f"[#321] wrote {SUMMARY_JSONL}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
