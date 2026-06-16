"""#341 Phase 1 Part A -- finer (k1, k2) sweep + phase robustness at Neptune
Proteus-Triton with the #320 Vector B basin-search convention.

#320 Vector B (commit 0e6f3f2) ran an initial epoch-aware sweep on Neptune
over moons ``(Triton, Proteus)`` with ``n_rev_grid = (0, 1, 2, 3)`` using a
24-sample global-phase x 24-sample relative-offset basin (the
``_sweep_one_cycle`` convention). Result:

  * 30 cells enumerated, 0 SILVER at the 0.05 km/s gate.
  * Proteus-Triton-Proteus (1, 1) best near-miss at 0.0584 km/s
    (phase0=240 deg, rel_offset=255 deg, tof_scale=2.0,
    vinfs [1.983, 1.427, 1.925]).
  * KNOWN_CORPUS anchor overlap = 0 (lit-fresh class).

Per the brief, the strongest signal of the whole #320 sweep was the
combination of low residual (16% above gate) + zero anchor overlap.

CRITICAL: this script re-uses the EXACT #320 Vector B ``_sweep_one_cycle``
convention (phase x rel_offset basin), NOT the ``RepeatedMoonTarget.close``
single-phase grid. ``RepeatedMoonTarget.close`` freezes the rel_offset at
``2 pi * j / len(consts)`` and only sweeps the global phase, missing the
rel_offset=255 deg basin where #320 found 0.058. A direct port keeps the
A.1 / A.2 results comparable to the #320 baseline.

Two questions, mirroring the #312 Phase 1 Part A workflow at Uranus that
produced the Umbriel-Oberon SILVER (#339):

1. **Broader (k1, k2) grid.** #320 swept ``(0, 1, 2, 3)``. Widen to
   ``(0, 1, 2, 3, 4, 5)`` per the brief, so the (3, 2), (2, 3), (4, 3),
   (3, 4), (1, 3), (3, 1), and (3, 4) etc. corners are covered. Does any
   of those close below the 0.05 km/s gate?

2. **Robustness of the 0.0584 km/s near-miss.** Is that residual the
   GENOME ceiling or a phase-grid resolution artifact? #320 sampled at
   (n_phase=24, n_offset=24). Re-run Proteus-Triton-Proteus (1, 1) at
   (n_phase=n_offset in {12, 24, 48, 96}) and report the best residual at
   each resolution.

NO catalogue writeback. NO novelty claims. Frame in commit: closure
quantified at best cell across the wider grid.

Sourced moon-GM and SMA values from JPL DE440 satellites registry
``src/cyclerfinder/core/satellites.py``:
  * Neptune system GM 6.836527100580e6 km^3/s^2 (JPL DE440 planetary)
  * Triton GM 1428.49546 km^3/s^2, a 354800 km, R 1352.6 km
  * Proteus GM 2.58342 km^3/s^2, a 117600 km, R 208.0 km
  (all from JPL SSD phys_par + sat441 mean elements, accessed 2026-06-14).

Outputs:
  * ``data/scan_341_neptune_proteus_triton_finer.jsonl`` (wider grid)
  * ``data/scan_341_neptune_robustness.jsonl`` (phase-resolution sweep)

Run as::

    uv run python scripts/scan_341_neptune_proteus_triton_finer.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.literature_check import (  # noqa: E402
    CandidateSignature,
    _candidate_anchors,
)
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    candidate_passes_physical_gate,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# #320 Vector B closure-residual machinery (direct port from
# scripts/scan_320_epoch_aware_moon_systems.py for comparability).
# ---------------------------------------------------------------------------


def _close_one(
    *,
    seq: tuple[str, ...],
    nrevs: tuple[int, ...],
    theta_per_moon: dict[str, float],
    tof_scale: float,
    consts: dict[str, tuple[float, float]],
    mu: float,
) -> tuple[bool, float, tuple[float, ...], tuple[float, ...]]:
    """Compute the closure residual of one (sequence, n_rev, phasing, tof_scale).

    Returns ``(feasible, residual_kms, vinf_per_encounter, tof_per_leg_days)``.
    Residual is the worst V_inf-magnitude continuity defect over all encounters
    INCLUDING the closed-cycle anchor wrap (matches the production genome).
    """
    n_legs = len(seq) - 1
    tofs: list[float] = []
    for k in range(n_legs):
        _, na = consts[seq[k]]
        _, nb = consts[seq[k + 1]]
        pa = 2.0 * math.pi / na
        pb = 2.0 * math.pi / nb
        tofs.append(tof_scale * math.sqrt(pa * pb))
    epochs = [0.0]
    for tof in tofs:
        epochs.append(epochs[-1] + tof)

    states = []
    for m, t in zip(seq, epochs, strict=True):
        sma, n = consts[m]
        states.append(_moon_state(theta_per_moon[m], n, t, sma, mu))

    vinf_in: list[float | None] = [None] * len(seq)
    vinf_out: list[float | None] = [None] * len(seq)
    for k in range(n_legs):
        r_a, v_a = states[k]
        r_b, v_b = states[k + 1]
        sols = lambert(r_a, r_b, tofs[k] * DAY_S, mu=mu, max_revs=max(0, nrevs[k]))
        wanted = [s for s in sols if s.n_revs == nrevs[k]]
        if not wanted:
            return (False, math.inf, (), ())
        best = min(wanted, key=lambda s, va=v_a: float(np.linalg.norm(s.v1 - va)))
        vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
        vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b))

    worst = 0.0
    per_enc: list[float] = []
    for k in range(len(seq)):
        vi = vinf_in[k]
        vo = vinf_out[k]
        if vi is not None and vo is not None:
            worst = max(worst, abs(vi - vo))
        rep = vi if vi is not None else vo
        per_enc.append(rep if rep is not None else 0.0)
    wo0 = vinf_out[0]
    wi_n = vinf_in[-1]
    if wo0 is not None and wi_n is not None:
        worst = max(worst, abs(wo0 - wi_n))
    return (True, worst, tuple(per_enc), tuple(tofs))


def _sweep_one_cycle(
    *,
    seq: tuple[str, ...],
    nrevs: tuple[int, ...],
    consts: dict[str, tuple[float, float]],
    mu: float,
    n_phase: int,
    n_offset: int,
    tof_scales: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0),
) -> dict[str, Any]:
    """Sweep global phase x relative offset on a closed length-3 cycle.

    Mirrors ``scan_320_epoch_aware_moon_systems._sweep_one_cycle`` exactly.
    """
    anchor = seq[0]
    distinct = sorted(set(seq))
    if len(distinct) != 2:
        raise ValueError(
            f"_sweep_one_cycle assumes length-3 closed cycle with 2 distinct moons; got {seq}"
        )
    intermediate = next(m for m in distinct if m != anchor)

    best = {
        "residual_kms": math.inf,
        "phase0_deg": 0.0,
        "rel_offset_deg": 0.0,
        "tof_scale": 1.0,
        "vinf_per_encounter_kms": [],
        "tof_days": [],
    }
    for i in range(n_offset):
        rel_off = 2.0 * math.pi * i / n_offset
        for j in range(n_phase):
            phase0 = 2.0 * math.pi * j / n_phase
            theta = {anchor: phase0, intermediate: phase0 + rel_off}
            for ts in tof_scales:
                ok, res, vinfs, tofs = _close_one(
                    seq=seq,
                    nrevs=nrevs,
                    theta_per_moon=theta,
                    tof_scale=ts,
                    consts=consts,
                    mu=mu,
                )
                if ok and res < best["residual_kms"]:
                    best = {
                        "residual_kms": res,
                        "phase0_deg": math.degrees(phase0),
                        "rel_offset_deg": math.degrees(rel_off),
                        "tof_scale": ts,
                        "vinf_per_encounter_kms": list(vinfs),
                        "tof_days": list(tofs),
                    }
    return best


def _closed_length3_cycles(moons: tuple[str, ...]) -> list[tuple[str, ...]]:
    cycles: list[tuple[str, ...]] = []
    for a in moons:
        for b in moons:
            if a == b:
                continue
            cycles.append((a, b, a))
    return cycles


# ---------------------------------------------------------------------------
# Part A.1: wider (k1, k2) grid using the #320 Vector B sweep convention.
# ---------------------------------------------------------------------------


def run_finer_grid(
    *,
    out_path: Path,
    git_sha: str,
    moons: tuple[str, str] = ("Triton", "Proteus"),
    nrev_grid: tuple[int, ...] = (0, 1, 2, 3, 4, 5),
    n_phase: int = 24,
    n_offset: int = 24,
    silver_gate_kms: float = 0.05,
    near_miss_band_kms: float = 1.0,
) -> dict[str, Any]:
    """Sweep all closed length-3 cycles on ``moons`` over the wider (k1, k2) grid."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    primary = "Neptune"
    mu = PRIMARIES[primary]
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    cycles = _closed_length3_cycles(moons)
    # (k1, k2) cells per sequence (exclude trivial (0, 0)).
    cells_per_seq = len(nrev_grid) ** 2 - 1
    n_cells = len(cycles) * cells_per_seq
    print(
        f"[341-A.1] moons={moons} cycles={len(cycles)} "
        f"n_rev_grid={nrev_grid} cells={n_cells} "
        f"n_phase={n_phase} n_offset={n_offset}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    near_miss: list[dict[str, Any]] = []
    silver: list[dict[str, Any]] = []
    t0 = time.time()

    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": "#341 Phase 1 Part A.1 -- Neptune Proteus-Triton wider (k1,k2) grid",
            "primary": primary,
            "moons": list(moons),
            "n_rev_grid": list(nrev_grid),
            "n_phase": n_phase,
            "n_offset": n_offset,
            "tof_scales": [0.5, 1.0, 1.5, 2.0],
            "silver_gate_kms": silver_gate_kms,
            "near_miss_band_kms": near_miss_band_kms,
            "reference_320_residual_kms": 0.0584,
            "reference_320_config": {
                "sequence": ["Proteus", "Triton", "Proteus"],
                "n_rev": [1, 1],
                "phase0_deg": 240.0,
                "rel_offset_deg": 255.0,
                "tof_scale": 2.0,
            },
            "system_mu_source": "JPL DE440 / src/cyclerfinder/core/satellites.py",
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        idx = 0
        for cycle in cycles:
            for nr1 in nrev_grid:
                for nr2 in nrev_grid:
                    if nr1 == 0 and nr2 == 0:
                        continue
                    nrevs = (nr1, nr2)
                    best = _sweep_one_cycle(
                        seq=cycle,
                        nrevs=nrevs,
                        consts=consts,
                        mu=mu,
                        n_phase=n_phase,
                        n_offset=n_offset,
                    )
                    if not math.isfinite(best["residual_kms"]):
                        idx += 1
                        continue
                    row: dict[str, Any] = {
                        "idx": idx,
                        "primary": primary,
                        "sequence": list(cycle),
                        "n_rev": list(nrevs),
                        **best,
                    }
                    res = float(best["residual_kms"])
                    vinfs_t: tuple[float, ...] = tuple(
                        float(v) for v in best["vinf_per_encounter_kms"]
                    )
                    if vinfs_t and all(v > 0.0 for v in vinfs_t):
                        passed, verdicts = candidate_passes_physical_gate(
                            tuple(cycle),
                            vinfs_t,
                            min_useful_bend_deg=5.0,
                        )
                        row["physical_gate_passed"] = bool(passed)
                        row["max_bend_deg_per_enc"] = [v.max_bend_deg for v in verdicts]
                    else:
                        row["physical_gate_passed"] = False
                        row["max_bend_deg_per_enc"] = []
                    sig = CandidateSignature(
                        primary=primary,
                        sequence=tuple(cycle),
                        vinf_per_encounter_kms=vinfs_t,
                        n_rev=tuple(int(x) for x in nrevs),
                    )
                    anchors = _candidate_anchors(sig)
                    row["n_lit_corpus_anchors_overlap"] = len(anchors)
                    row["lit_corpus_anchor_names"] = [a.name for a in anchors][:3]
                    is_silver = bool(res < silver_gate_kms and row["physical_gate_passed"])
                    row["is_silver"] = is_silver
                    is_lit_fresh_class = bool(len(anchors) == 0)
                    row["is_lit_fresh_class"] = is_lit_fresh_class
                    if res < near_miss_band_kms:
                        fh.write(json.dumps(row) + "\n")
                        fh.flush()
                        rows.append(row)
                        near_miss.append(row)
                        if is_silver:
                            silver.append(row)
                            print(
                                f"[341-A.1] SILVER seq={cycle} n_rev={nrevs} "
                                f"res={res:.4f} km/s vinfs={vinfs_t} "
                                f"physical={row['physical_gate_passed']} "
                                f"anchors={row['n_lit_corpus_anchors_overlap']}",
                                flush=True,
                            )
                    idx += 1

        near_miss.sort(key=lambda r: r["residual_kms"])
        silver.sort(key=lambda r: r["residual_kms"])
        elapsed = time.time() - t0

        summary = {
            "_meta": True,
            "kind": "summary",
            "primary": primary,
            "n_cells_evaluated": idx,
            "n_below_near_miss_band": len(near_miss),
            "n_silver": len(silver),
            # near_miss is sorted ascending by residual; rows is insertion order.
            "best_residual_kms": (near_miss[0]["residual_kms"] if near_miss else None),
            "best_record_brief": (
                {
                    "sequence": near_miss[0]["sequence"],
                    "n_rev": near_miss[0]["n_rev"],
                    "residual_kms": near_miss[0]["residual_kms"],
                    "vinf_per_encounter_kms": near_miss[0]["vinf_per_encounter_kms"],
                    "phase0_deg": near_miss[0]["phase0_deg"],
                    "rel_offset_deg": near_miss[0]["rel_offset_deg"],
                    "tof_scale": near_miss[0]["tof_scale"],
                    "is_silver": near_miss[0]["is_silver"],
                    "is_lit_fresh_class": near_miss[0]["is_lit_fresh_class"],
                }
                if near_miss
                else None
            ),
            "top5_near_miss_brief": [
                {
                    "sequence": r["sequence"],
                    "n_rev": r["n_rev"],
                    "residual_kms": r["residual_kms"],
                    "is_silver": r["is_silver"],
                    "is_lit_fresh_class": r["is_lit_fresh_class"],
                }
                for r in near_miss[:5]
            ],
            "elapsed_s": elapsed,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()
        print(
            f"[341-A.1] DONE -- {idx} cells, {len(near_miss)} near-miss, "
            f"{len(silver)} SILVER, best={summary['best_residual_kms']!r} "
            f"in {elapsed:.1f}s",
            flush=True,
        )
    return summary


# ---------------------------------------------------------------------------
# Part A.2: phase-resolution robustness on Proteus-Triton-Proteus (1, 1).
# ---------------------------------------------------------------------------


def run_robustness_sweep(
    *,
    out_path: Path,
    git_sha: str,
    moons: tuple[str, str] = ("Triton", "Proteus"),
    phase_samples_grid: tuple[int, ...] = (12, 24, 48, 96),
) -> list[dict[str, Any]]:
    """Re-close Proteus-Triton-Proteus (1, 1) at (n_phase=n_offset=ps) for each ps.

    Uses the same _sweep_one_cycle convention as Part A.1 / #320 Vector B,
    so the residuals are directly comparable.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    primary = "Neptune"
    mu = PRIMARIES[primary]
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    rows: list[dict[str, Any]] = []
    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": "#341 Phase 1 Part A.2 -- Proteus-Triton-Proteus (1,1) "
            "phase-resolution robustness via #320 Vector B basin convention",
            "primary": primary,
            "sequence": ["Proteus", "Triton", "Proteus"],
            "n_rev": [1, 1],
            "phase_samples_grid": list(phase_samples_grid),
            "tof_resonance_grid": [0.5, 1.0, 1.5, 2.0],
            "reference_320_residual_kms": 0.0584,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        for ps in phase_samples_grid:
            t0 = time.time()
            best = _sweep_one_cycle(
                seq=("Proteus", "Triton", "Proteus"),
                nrevs=(1, 1),
                consts=consts,
                mu=mu,
                n_phase=ps,
                n_offset=ps,
            )
            row = {
                "phase_samples": ps,
                "n_phase": ps,
                "n_offset": ps,
                "sequence": ["Proteus", "Triton", "Proteus"],
                "n_rev": [1, 1],
                "residual_kms": best["residual_kms"],
                "phase0_deg": best["phase0_deg"],
                "rel_offset_deg": best["rel_offset_deg"],
                "tof_scale": best["tof_scale"],
                "vinf_per_encounter_kms": list(best["vinf_per_encounter_kms"]),
                "tof_days": list(best["tof_days"]),
                "elapsed_s": time.time() - t0,
                "git_sha": git_sha,
            }
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            rows.append(row)
            print(
                f"[341-A.2] ps={ps} residual={row['residual_kms']!r} km/s "
                f"phase0={row['phase0_deg']:.1f} rel_off={row['rel_offset_deg']:.1f} "
                f"tof_scale={row['tof_scale']} elapsed={row['elapsed_s']:.1f}s",
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
    """GENOME-CEILING vs PHASE-RESOLUTION-ARTIFACT (same verdict logic as #312)."""
    residuals = [r["residual_kms"] for r in rows if r.get("residual_kms") is not None]
    if not residuals:
        return {"label": "NO_DATA", "min_res": None, "max_res": None, "spread": None}
    rmin = float(min(residuals))
    rmax = float(max(residuals))
    spread = rmax - rmin
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
            "GENOME_CEILING": "spread < 5e-3 km/s -- finer phase x rel_offset grid "
            "does not move the residual; the 0.058 km/s near-miss is the GENOME "
            "ceiling (Lambert single-shot coplanar). To close the gap, extend "
            "the genome (DSM legs, multi-arc, low-thrust, 3D).",
            "PHASE_RESOLUTION_CLOSED_GAP": "a finer phase x rel_offset grid found "
            "a sub-gate residual -- the original 0.058 was a phase-resolution "
            "artifact; the cell DOES close ballistically.",
            "PHASE_RESOLUTION_PARTIAL": "finer phase x rel_offset grid moved the "
            "residual but not below the gate -- partial phase-sensitivity, but "
            "the ballistic ceiling still stops short of the gate.",
        }.get(label, "(no interpretation)"),
    }


def main() -> int:
    sha = _git_sha()
    print(f"[341-A] Neptune Proteus-Triton extended sweep -- sha={sha}", flush=True)

    out_finer = ROOT / "data" / "scan_341_neptune_proteus_triton_finer.jsonl"
    print(
        f"[341-A.1] Wider (k1, k2) grid n_rev=(0..5) at Neptune (Triton, Proteus) "
        f"length-3 -- out={out_finer}",
        flush=True,
    )
    sum_a1 = run_finer_grid(out_path=out_finer, git_sha=sha)
    print(f"[341-A.1] Finer-grid summary: best={sum_a1['best_residual_kms']!r}", flush=True)

    out_robust = ROOT / "data" / "scan_341_neptune_robustness.jsonl"
    print(
        f"[341-A.2] Phase-resolution robustness at Proteus-Triton-Proteus (1,1) "
        f"-- out={out_robust}",
        flush=True,
    )
    rows = run_robustness_sweep(out_path=out_robust, git_sha=sha)
    verdict = _robustness_verdict(rows)
    print(f"[341-A.2] Robustness verdict: {verdict}", flush=True)

    print("[341-A] DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
