"""#311 Phase 1 Part A -- finer (k1, k2) sweep at Saturn Rhea-Dione.

The #285 Saturn campaign left Rhea-Dione-Rhea (1, 1) closing at residual
0.10688 km/s -- 2.1x above the 0.05 km/s gate (vs Uranus Oberon-Titania at
1.2x). Direct analog of the #312 Uranus Part A treatment: widen the
(k1, k2) grid + sweep phase resolution to confirm whether the 0.107 km/s
residual is the genome ceiling or a phase-grid artifact, AND surface any
sub-gate cells that the original (1, 2, 3) grid missed.

NO catalogue writeback. NO novelty claims. Sourced Saturn-system mu and
moon-SMA values come from JPL DE440 satellites registry already in
``src/cyclerfinder/core/satellites.py`` (Saturn system GM 3.7931207e7
km^3/s^2; Rhea mu 153.94 / a 527070 km; Dione mu 73.116 / a 377420 km --
all JPL SSD phys_par + sat441 mean elements, accessed 2026-06-14).

Outputs:
  * ``data/scan_311_saturn_rhea_dione_finer.jsonl`` (broader grid)
  * ``data/scan_311_saturn_rhea_dione_robustness.jsonl`` (phase sweep)

Run as::

    uv run python scripts/scan_311_saturn_rhea_dione_finer.py
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

from cyclerfinder.search.discovery_campaign import RepeatedMoonTarget  # noqa: E402
from cyclerfinder.search.saturn_uranus_campaign import (  # noqa: E402
    run_prioritized_scan,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# Part A.1: wider (k1, k2) grid over the Saturn near-miss moon set.
# CRITICAL: same moon set as #285 ((Titan, Rhea, Dione, Tethys)) so the
# phase-grid offsets are identical to that scan -- restricting to
# {Rhea, Dione} would change ``len(consts)`` and shift per-moon initial
# longitudes, breaking comparability with the 0.10688 km/s baseline.
# n_rev_grid covers (k1, k2) in [0, 5]^2 -- 36 cells per length-3 sequence,
# including the brief's targeted (2,3), (4,3), (3,4), (1,3), (3,1),
# (5,2), (2,5) cells.
SATURN_BAND_A_FINER = dict(
    primary="Saturn",
    moons=("Titan", "Rhea", "Dione", "Tethys"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),  # 36 (k1, k2) cells per sequence
)


def run_robustness_sweep(
    *,
    out_path: Path,
    git_sha: str,
    phase_samples_grid: tuple[int, ...] = (12, 24, 48, 96),
) -> list[dict[str, Any]]:
    """Re-close Rhea-Dione-Rhea (1, 1) at increasing phase resolution.

    The #285 Saturn near-miss was 0.10688 km/s at phase_samples=12; this
    sweep asks whether finer phase resolution moves the residual (=> the
    near-miss is a phase-grid artifact) or holds it (=> genome ceiling,
    Lambert single-shot coplanar can do no better).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": "#311 Phase 1 Part A.2 -- Rhea-Dione-Rhea (1,1) phase-resolution robustness",
            "primary": "Saturn",
            "sequence": ["Rhea", "Dione", "Rhea"],
            "n_rev": [1, 1],
            "phase_samples_grid": list(phase_samples_grid),
            "tof_resonance_grid": [0.5, 1.0, 1.5, 2.0],
            "reference_285_residual_kms": 0.10688173280775803,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        for ps in phase_samples_grid:
            # Same moon set as #285 (Titan, Rhea, Dione, Tethys) so the
            # per-moon longitude offsets match.
            target = RepeatedMoonTarget(
                primary="Saturn",
                moons=("Titan", "Rhea", "Dione", "Tethys"),
                seq_lengths=(3,),
                n_rev_grid=(1,),  # ONLY (1, 1) for the robustness check
                n_phase_samples=ps,
                tof_resonance_grid=(0.5, 1.0, 1.5, 2.0),
                git_sha=git_sha,
            )
            t0 = time.time()
            best_target_row: dict[str, Any] | None = None
            for cand in target.enumerate_candidates():
                if tuple(cand.sequence) != ("Rhea", "Dione", "Rhea"):
                    continue
                if tuple(cand.payload["n_rev"]) != (1, 1):
                    continue
                closure = target.close(cand)
                row = {
                    "phase_samples": ps,
                    "sequence": list(cand.sequence),
                    "n_rev": list(cand.payload["n_rev"]),
                    "converged": closure.converged,
                    "residual_kms": closure.residual_kms,
                    "vinf_per_encounter_kms": list(closure.vinf_per_encounter_kms),
                    "tof_days": list(closure.tof_days),
                    "elapsed_s": time.time() - t0,
                    "git_sha": git_sha,
                }
                best_target_row = row
                break
            if best_target_row is None:
                best_target_row = {
                    "phase_samples": ps,
                    "sequence": ["Rhea", "Dione", "Rhea"],
                    "n_rev": [1, 1],
                    "converged": False,
                    "residual_kms": None,
                    "vinf_per_encounter_kms": [],
                    "tof_days": [],
                    "elapsed_s": time.time() - t0,
                    "git_sha": git_sha,
                    "note": "candidate not enumerated (deterministic order may differ)",
                }
            fh.write(json.dumps(best_target_row) + "\n")
            fh.flush()
            rows.append(best_target_row)
            print(
                f"[311-A.2] phase_samples={ps} residual="
                f"{best_target_row['residual_kms']!r} km/s elapsed="
                f"{best_target_row['elapsed_s']:.1f}s",
                flush=True,
            )

        summary = {
            "_meta": True,
            "kind": "summary",
            "robustness_rows": rows,
            "verdict": _robustness_verdict(rows),
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()
    return rows


def _robustness_verdict(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Decide GENOME-CEILING vs PHASE-RESOLUTION-ARTIFACT from the rows."""
    residuals = [r["residual_kms"] for r in rows if r["residual_kms"] is not None]
    if not residuals:
        return {"label": "NO_DATA", "min_res": None, "max_res": None, "spread": None}
    rmin = float(min(residuals))
    rmax = float(max(residuals))
    spread = rmax - rmin
    # Same criteria as #312 Part A.2: if spread < 5e-3 km/s, GENOME CEILING.
    if spread < 5e-3:
        label = "GENOME_CEILING"
    elif rmin < 0.05:
        label = "PHASE_RESOLUTION_CLOSED_GAP"
    else:
        label = "PHASE_RESOLUTION_PARTIAL"
    return {
        "label": label,
        "min_res": rmin,
        "max_res": rmax,
        "spread": spread,
        "interpretation": {
            "GENOME_CEILING": "spread < 5e-3 km/s -- finer phase grid does "
            "not move the residual; the 0.107 km/s near-miss is the GENOME "
            "ceiling (Lambert single-shot coplanar). To close the gap, "
            "extend the genome (DSM legs, multi-arc, low-thrust, 3D).",
            "PHASE_RESOLUTION_CLOSED_GAP": "a finer phase grid found a "
            "sub-gate residual -- the original 0.107 was a phase-resolution "
            "artifact; the cell DOES close ballistically.",
            "PHASE_RESOLUTION_PARTIAL": "finer phase grid moved the residual "
            "but not below the gate -- partial phase-sensitivity, but the "
            "ballistic ceiling still stops short of the gate.",
        }.get(label, "(no interpretation)"),
    }


def main() -> int:
    sha = _git_sha()
    print(f"[311-A] Saturn near-miss extended sweep -- sha={sha}", flush=True)

    # Part A.1: wider (k1, k2) grid.
    out_finer = ROOT / "data" / "scan_311_saturn_rhea_dione_finer.jsonl"
    print(
        f"[311-A.1] Wider (k1, k2) grid n_rev=(0..5) at Saturn length-3 -- out={out_finer}",
        flush=True,
    )
    sum_a1 = run_prioritized_scan(
        out_path=out_finer,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **SATURN_BAND_A_FINER,
    )
    print(f"[311-A.1] Finer-grid summary: {sum_a1.as_dict()}", flush=True)

    with out_finer.open("a", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "summary",
                    "summary": sum_a1.as_dict(),
                    "git_sha": sha,
                }
            )
            + "\n"
        )

    # Part A.2: phase-resolution robustness on the (1, 1) cell.
    out_robust = ROOT / "data" / "scan_311_saturn_rhea_dione_robustness.jsonl"
    print(
        f"[311-A.2] Phase-resolution robustness at Rhea-Dione-Rhea (1,1) -- out={out_robust}",
        flush=True,
    )
    rows = run_robustness_sweep(
        out_path=out_robust,
        git_sha=sha,
        phase_samples_grid=(12, 24, 48, 96),
    )
    verdict = _robustness_verdict(rows)
    print(f"[311-A.2] Robustness verdict: {verdict}", flush=True)

    print("[311-A] DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
