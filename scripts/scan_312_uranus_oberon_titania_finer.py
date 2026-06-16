"""#312 Phase 1 Part A -- finer (k1, k2) sweep at Uranus Oberon-Titania.

The #285 Uranus campaign closed Oberon-Titania-Oberon (1,1) at
residual 0.0617 km/s -- only 25% above the 0.05 km/s gate. This is the
session's highest novelty-leverage near-miss because (a) it is close to the
gate, and (b) there is NO published Uranian-cycler prior in our KNOWN_CORPUS.

Part A asks two questions with the existing #254 repeated-moon multi-rev
genome (no genome changes; READ-ONLY on the search machinery):

1. **Broader (k1, k2) grid.** #285 swept ``n_rev_grid = (0, 1, 2, 3)`` --
   16 (k1, k2) cells per length-3 sequence. We widen to
   ``n_rev_grid = (0, 1, 2, 3, 4, 5)`` -- 36 cells per sequence -- adding
   (2,3), (4,3), (3,4), (1,3), (3,1), (5,2), (2,5) and intermediates that
   the brief lists. Does any of those close below the 0.05 km/s gate?

2. **Robustness of the 0.0617 km/s near-miss.** Is that residual the
   GENOME ceiling or a phase-grid resolution artifact? #285 ran at
   ``phase_samples = 12``. We re-run Oberon-Titania-Oberon (1,1) at
   ``phase_samples in {12, 24, 48, 96}`` and report the best residual at
   each resolution. The Saturn near-miss verdict (#285) was that 0.107 km/s
   held across {12, 24, 48}; we do the same robustness here.

NO catalogue writeback. NO novelty claims. Frame in commit: closure
quantified at best cell across the wider grid; if any cell closes <= gate,
the SILVER row is plumbed for #305/#306 gauntlet ingestion (not catalogue
admission here).

Sourced moon-GM and SMA values come from JPL DE440 satellites registry
``src/cyclerfinder/core/satellites.py`` (Uranus system GM 5.7945564e6
km^3/s^2; Titania mu 226.9 / a 436298 km; Oberon mu 205.3 / a 583511 km;
all from JPL SSD phys_par + sat441 mean elements).

Outputs:
  * ``data/scan_312_uranus_oberon_titania_finer.jsonl`` (broader grid)
  * ``data/scan_312_uranus_robustness.jsonl`` (phase-resolution sweep)

Run as::

    uv run python scripts/scan_312_uranus_oberon_titania_finer.py
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


# Part A.1: wider (k1, k2) grid over Uranus regulars (length-3 cycles).
# CRITICAL: use the EXACT same moon set as #285 ((Titania, Oberon, Umbriel))
# so the phase-grid offsets are identical to that scan -- restricting to
# {Oberon, Titania} changes ``len(consts)`` and shifts the per-moon initial
# longitudes, breaking comparability with the 0.0617 km/s baseline. The
# wider n_rev_grid covers all the (k1, k2) pairs the brief lists, including
# (2,3), (4,3), (3,4), (1,3), (3,1), (5,2), (2,5), plus the (1,1)-(3,3)
# block #285 had + (4,*) and (5,*) extensions.
URANUS_BAND_A_FINER = dict(
    primary="Uranus",
    moons=("Titania", "Oberon", "Umbriel"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),  # 36 (k1, k2) cells per sequence
)


def run_robustness_sweep(
    *,
    out_path: Path,
    git_sha: str,
    phase_samples_grid: tuple[int, ...] = (12, 24, 48, 96),
) -> list[dict[str, Any]]:
    """Re-close Oberon-Titania-Oberon (1,1) at increasing phase resolution.

    The brief: "rerun the Oberon-Titania-Oberon (1,1) cell at phase-grid
    {12, 24, 48, 96} samples to confirm the 0.062 is the genome ceiling or
    a phase-resolution artifact." Each row is one (phase_samples, residual)
    pair plus best vinf/ToF triplet.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": "#312 Phase 1 Part A.2 -- Oberon-Titania-Oberon (1,1) "
            "phase-resolution robustness",
            "primary": "Uranus",
            "sequence": ["Oberon", "Titania", "Oberon"],
            "n_rev": [1, 1],
            "phase_samples_grid": list(phase_samples_grid),
            "tof_resonance_grid": [0.5, 1.0, 1.5, 2.0],
            "reference_285_residual_kms": 0.0617,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        for ps in phase_samples_grid:
            # CRITICAL: same moon set as #285 (Titania, Oberon, Umbriel)
            # so the per-moon longitude offsets match. Restricting to
            # {Oberon, Titania} would change ``len(consts)`` and shift the
            # offsets, returning a non-comparable residual.
            target = RepeatedMoonTarget(
                primary="Uranus",
                moons=("Titania", "Oberon", "Umbriel"),
                seq_lengths=(3,),
                n_rev_grid=(1,),  # ONLY (1, 1) for the robustness check
                n_phase_samples=ps,
                tof_resonance_grid=(0.5, 1.0, 1.5, 2.0),
                git_sha=git_sha,
            )
            t0 = time.time()
            # Find the specific Oberon-Titania-Oberon (1,1) candidate
            # (skip the trivial Titania-Oberon-Titania (1,1) too -- both are
            # in the n_rev=(1,) enumeration).
            best_target_row: dict[str, Any] | None = None
            for cand in target.enumerate_candidates():
                if tuple(cand.sequence) != ("Oberon", "Titania", "Oberon"):
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
                    "sequence": ["Oberon", "Titania", "Oberon"],
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
                f"[312-A.2] phase_samples={ps} residual="
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
    # If all four phase-resolution residuals fall within +/- 5e-3 km/s (1/12
    # of the gate-overage of #285's 0.0617 km/s), call it a GENOME CEILING:
    # finer phase-grid does NOT close the gap, so the bottleneck is the
    # Lambert genome (single-ellipse coplanar single-shot per leg), not
    # how densely we sample initial moon longitudes.
    if spread < 5e-3:
        label = "GENOME_CEILING"
    elif rmin < 0.05:
        # The fine grid found a sub-gate residual -- phase resolution mattered!
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
            "not move the residual; the 0.062 km/s near-miss is the GENOME "
            "ceiling (Lambert single-shot coplanar). To close the gap, "
            "extend the genome (DSM legs, multi-arc, low-thrust, 3D).",
            "PHASE_RESOLUTION_CLOSED_GAP": "a finer phase grid found a "
            "sub-gate residual -- the original 0.062 was a phase-resolution "
            "artifact; the cell DOES close ballistically.",
            "PHASE_RESOLUTION_PARTIAL": "finer phase grid moved the residual "
            "but not below the gate -- partial phase-sensitivity, but the "
            "ballistic ceiling still stops short of the gate.",
        }.get(label, "(no interpretation)"),
    }


def main() -> int:
    sha = _git_sha()
    print(f"[312-A] Uranus near-miss extended sweep -- sha={sha}", flush=True)

    # Part A.1: wider (k1, k2) grid.
    out_finer = ROOT / "data" / "scan_312_uranus_oberon_titania_finer.jsonl"
    print(
        f"[312-A.1] Wider (k1, k2) grid n_rev=(0..5) at Oberon-Titania length-3 -- out={out_finer}",
        flush=True,
    )
    sum_a1 = run_prioritized_scan(
        out_path=out_finer,
        gate_residual_kms=0.05,
        phase_samples=12,
        git_sha=sha,
        **URANUS_BAND_A_FINER,
    )
    print(f"[312-A.1] Finer-grid summary: {sum_a1.as_dict()}", flush=True)

    # Append summary row to the finer scan.
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

    # Part A.2: phase-resolution robustness on the (1,1) cell.
    out_robust = ROOT / "data" / "scan_312_uranus_robustness.jsonl"
    print(
        f"[312-A.2] Phase-resolution robustness at Oberon-Titania-Oberon (1,1) -- out={out_robust}",
        flush=True,
    )
    rows = run_robustness_sweep(
        out_path=out_robust,
        git_sha=sha,
        phase_samples_grid=(12, 24, 48, 96),
    )
    verdict = _robustness_verdict(rows)
    print(f"[312-A.2] Robustness verdict: {verdict}", flush=True)

    print("[312-A] DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
