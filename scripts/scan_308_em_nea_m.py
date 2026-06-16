"""#308 Phase 1 — NEA-augmented cycler discovery sweep (Earth -> NEA -> Mars).

Per task #308 the fresh-ground asteroid-leveraging class is the sparsest of
the three #302 candidates (asteroid leveraging / low-thrust / non-E-M). This
script enumerates length-3 heliocentric chains ``E -> NEA -> M`` over the
:data:`cyclerfinder.search.asteroid_leveraging.LARGEST_NEAS` pool (10 anchors)
and writes one JSONL row per chain that clears the configured gates.

PHASE 1 EXPECTED VERDICT: most or all NEA candidates fail the physical-sanity
gate. NEAs are tiny; their gravity wells are essentially zero compared to a
planet. The patched-conic max bend at any NEA in the pool at typical cycler
V_inf (3-7 km/s) is well below the 5 deg useful floor. The script honestly
reports what survives — likely 0 candidates pass the gate at the default
V_inf grid; the gate's job is to reject these.

This run produces TWO JSONL summary buckets:

* **All Lambert-converging candidates** (no physical-sanity gate). This maps
  the geometric / closure region per NEA so a Phase 2 (low-thrust)
  re-evaluation can target the bend-required-but-NEA-can't-provide gap.
* **Gate-passing survivors** (with the physical-sanity gate ON). The real
  Phase 1 verdict: how many ballistic NEA flybys are physically meaningful.

NO catalogue writeback. NO novelty claims.

Run as::

    uv run python scripts/scan_308_em_nea_m.py

Outputs ``data/scan_308_em_nea_m.jsonl`` — one leading ``_meta`` row, then
per-candidate rows (gate-passing only, by default), then a trailing summary.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.search.asteroid_leveraging import (  # noqa: E402
    LARGEST_NEAS,
    NEACyclerCandidate,
    search_nea_augmented_cyclers,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------
# Sweep grid
# --------------------------------------------------------------------------
#
# V_inf grid covers the operational cycler regime (3-7 km/s) — the regime
# the search machinery exists to support. The gate-passing expectation is
# zero at these V_inf for every NEA in the pool, by construction of the gate
# (the Phase 1 honest verdict). The grid is symmetric on either side of the
# Aldrin cycler V_inf (~5.6 km/s at Earth, ~3.0 km/s at Mars) so the report
# spans both sides of the natural shell.
VINF_GRID_KMS: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0)

# TOF grid covers the broad Earth-Mars / Earth-NEA leg range. Aldrin's
# baseline outbound is 146 d; the NEAs span 0.92-2.66 AU semi-major axis so
# the relevant TOF range needs to cover both shorter (Apophis, perihelion
# near Earth) and longer (Ganymed, Amor group beyond Mars).
TOF_BOX_DAYS_PER_LEG: tuple[float, float] = (60.0, 500.0)
N_TOF_SAMPLES: int = 6  # 6x6 = 36 (TOF, TOF) pairs per (NEA, V_inf)

# Closure floors: match #302's defaults (the previous E-M structural sweep).
CLOSURE_FLOOR_KMS: float = 0.5
FLYBY_CONTINUITY_FLOOR_KMS: float = 0.5

# Physical-sanity floor: same 5 deg as #324 (the Umbriel-prompted gate).
MIN_USEFUL_BEND_DEG: float = 5.0


def _candidate_to_dict(cand: NEACyclerCandidate) -> dict[str, object]:
    """Render an :class:`NEACyclerCandidate` to a JSON-serialisable dict."""
    d = asdict(cand)
    # tuples -> lists for JSON.
    for key, value in list(d.items()):
        if isinstance(value, tuple):
            d[key] = list(value)
    return d


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_308_em_nea_m.jsonl"
    print(f"[308] starting -- sha={sha} -- out={out_path}", flush=True)

    eph = Ephemeris("circular")
    n_neas = len(LARGEST_NEAS)
    n_vinf = len(VINF_GRID_KMS)
    n_tof = N_TOF_SAMPLES
    cells_per_nea = n_vinf * n_tof * n_tof
    total_cells = n_neas * cells_per_nea
    print(
        f"[308] sweep grid: {n_neas} NEAs x {n_vinf} V_inf x {n_tof}x{n_tof} TOFs "
        f"= {total_cells} cells",
        flush=True,
    )

    meta_row = {
        "_meta": True,
        "kind": "config",
        "scan": "308_em_nea_m",
        "task": "#308 asteroid-leveraging Phase 1",
        "git_sha": sha,
        "primary_sequence": ["E", "M"],
        "nea_pool_names": [nea.name for nea in LARGEST_NEAS],
        "nea_pool_designations": [nea.designation for nea in LARGEST_NEAS],
        "vinf_grid_kms": list(VINF_GRID_KMS),
        "tof_box_days_per_leg": list(TOF_BOX_DAYS_PER_LEG),
        "n_tof_samples": N_TOF_SAMPLES,
        "closure_floor_kms": CLOSURE_FLOOR_KMS,
        "flyby_continuity_floor_kms": FLYBY_CONTINUITY_FLOOR_KMS,
        "min_useful_bend_deg": MIN_USEFUL_BEND_DEG,
        "ephemeris": "circular",
        "discipline": (
            "NO catalogue writeback. NO novelty claims. Frame: 'Phase 1 candidate; "
            "passes physical-sanity gate; awaits literature check'. Most NEA "
            "candidates expected to FAIL the physical-sanity gate at cycler V_inf "
            "— that IS the right answer for ballistic NEA flybys."
        ),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write(json.dumps(meta_row) + "\n")

    # ---- Pass 1: pre-gate (Lambert + closure only) — diagnostic counts ----
    print("[308] pass 1: pre-gate Lambert + closure (diagnostic counts only)", flush=True)
    pre_gate_candidates = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=LARGEST_NEAS,
            vinf_grid_kms=VINF_GRID_KMS,
            tof_box_days_per_leg=TOF_BOX_DAYS_PER_LEG,
            n_tof_samples=N_TOF_SAMPLES,
            closure_floor_kms=CLOSURE_FLOOR_KMS,
            flyby_continuity_floor_kms=FLYBY_CONTINUITY_FLOOR_KMS,
            min_useful_bend_deg=MIN_USEFUL_BEND_DEG,
            use_physical_sanity_gate=False,
            ephemeris=eph,
        )
    )
    print(f"[308] pre-gate Lambert/closure passers: {len(pre_gate_candidates)}", flush=True)

    # Per-NEA breakdown.
    pre_gate_by_nea: dict[str, int] = {}
    for cand in pre_gate_candidates:
        nea_name = cand.nea_in_sequence[0]
        pre_gate_by_nea[nea_name] = pre_gate_by_nea.get(nea_name, 0) + 1
    print(f"[308] pre-gate per-NEA: {pre_gate_by_nea}", flush=True)

    # ---- Pass 2: post-gate (physical-sanity ON) — the real Phase 1 verdict ----
    print("[308] pass 2: post-gate (physical-sanity ON)", flush=True)
    post_gate_candidates = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=LARGEST_NEAS,
            vinf_grid_kms=VINF_GRID_KMS,
            tof_box_days_per_leg=TOF_BOX_DAYS_PER_LEG,
            n_tof_samples=N_TOF_SAMPLES,
            closure_floor_kms=CLOSURE_FLOOR_KMS,
            flyby_continuity_floor_kms=FLYBY_CONTINUITY_FLOOR_KMS,
            min_useful_bend_deg=MIN_USEFUL_BEND_DEG,
            use_physical_sanity_gate=True,
            ephemeris=eph,
        )
    )
    print(f"[308] gate-passing survivors: {len(post_gate_candidates)}", flush=True)

    # Per-NEA post-gate breakdown.
    post_gate_by_nea: dict[str, int] = {}
    for cand in post_gate_candidates:
        nea_name = cand.nea_in_sequence[0]
        post_gate_by_nea[nea_name] = post_gate_by_nea.get(nea_name, 0) + 1
    if post_gate_by_nea:
        print(f"[308] gate-passing per-NEA: {post_gate_by_nea}", flush=True)
    else:
        print(
            "[308] no NEA in the pool passes the physical-sanity gate at any "
            "(V_inf, TOF) — the Phase 1 expected verdict.",
            flush=True,
        )

    # ---- Write rows: write ALL gate-passing survivors, plus a CAPPED sample
    # of pre-gate-only candidates (for downstream Phase 2 diagnostic). ----
    with out_path.open("a", encoding="utf-8") as out:
        for cand in post_gate_candidates:
            row = _candidate_to_dict(cand)
            row["_meta"] = {"git_sha": sha, "scan": "308_em_nea_m", "pass": "post_gate"}
            out.write(json.dumps(row) + "\n")

        # Cap the pre-gate diagnostic at 200 rows (sampled by even stride
        # across the iterator order) so the JSONL stays small.
        pre_only = [c for c in pre_gate_candidates if not c.physical_sanity_passed]
        if pre_only:
            stride = max(1, len(pre_only) // 200)
            sampled_pre_only = pre_only[::stride][:200]
            for cand in sampled_pre_only:
                row = _candidate_to_dict(cand)
                row["_meta"] = {
                    "git_sha": sha,
                    "scan": "308_em_nea_m",
                    "pass": "pre_gate_diagnostic_sample",
                }
                out.write(json.dumps(row) + "\n")
        else:
            sampled_pre_only = []

        summary = {
            "_meta": True,
            "kind": "summary",
            "scan": "308_em_nea_m",
            "git_sha": sha,
            "enumerated_cells": total_cells,
            "lambert_closure_passers": len(pre_gate_candidates),
            "physical_sanity_passers": len(post_gate_candidates),
            "physical_sanity_rejects": len(pre_only),
            "pre_gate_per_nea": pre_gate_by_nea,
            "post_gate_per_nea": post_gate_by_nea,
            "diagnostic_sample_written": len(sampled_pre_only),
            "phase_1_verdict": (
                f"Gate-passing survivors at default V_inf grid (3-7 km/s) are "
                f"expected to be ~0 — every NEA in the pool has a useful-V_inf "
                f"ceiling near its surface escape speed (~0.2-23 m/s), six "
                f"orders of magnitude below cycler V_inf. The {len(post_gate_candidates)} "
                f"non-zero count IS the Phase 1 answer."
            ),
            "phase_2_recommendation": (
                f"Phase 2 (low-thrust V_inf modification at the NEA encounter) "
                f"may open the bend-required-but-NEA-can't-provide gap. The "
                f"{len(pre_gate_candidates)} Lambert-converging pre-gate candidates "
                f"map the geometric region that Phase 2 should target."
            ),
            "notes": (
                "Offline run: no literature_check yet (Phase 1 only produces "
                "the JSONL). literature_check is downstream and mandatory "
                "before any 'novel' framing."
            ),
        }
        out.write(json.dumps(summary) + "\n")

    print(
        f"[308] DONE -- enumerated={total_cells}, "
        f"lambert_passers={len(pre_gate_candidates)}, "
        f"gate_passers={len(post_gate_candidates)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
