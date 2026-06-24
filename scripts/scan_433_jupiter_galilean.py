"""#433 Jupiter Galilean repeated-moon quasi_cycler sweep.

The Jupiter Galilean system was untouched by the #320 epoch-aware non-EM
moon-system sweep (which covered Saturn, Neptune, Pluto). This script mirrors
the #320 Saturn block structurally and extends it on two axes per the #344
lesson that finer phase resolution recovers aliased closures:

  * Phase grid 48x48 (vs #320's 24x24): the Galilean Laplace 1:2:4 resonance
    (Io : Europa : Ganymede mean-motion commensurability) phase-locks many
    cells; a coarse grid aliases the true basin floor.
  * Per-leg n_rev in (0, 1, 2, 3, 4, 5) (vs #320's 0..3): wider revolution
    window, excluding only the trivial all-zero direct-transfer corner.

For each closed length-3 repeated-moon cycle over the Galilean set
[Io, Europa, Ganymede, Callisto] (+ Amalthea, which carries the needed
GM/elements but whose ~0.16 km3/s2 GM makes it useless as a flyby body --
the #324 physical-sanity gate rejects cycles that try to use it), the sweep:

  * builds the epoch-aware (sma, mean-motion) constants from JPL DE440 /
    SAT441-sourced satellite elements,
  * runs the epoch-aware Lambert closure per phase x relative-offset x
    tof-scale cell,
  * records the basin-floor residual,
  * runs the #324 physical max-bend gate at each encounter,
  * runs the offline #328 KNOWN_CORPUS literature-novelty anchor overlap
    (the Galilean corpus IS populated: Hernandez/Jones/Jesick, Strange/
    Campagnola/Russell, Niehoff -- so an anchor overlap means the family is
    published; anchor-overlap = 0 is NECESSARY-not-sufficient for novelty).

The resonance/irrationality screening of #320 is kept EXACTLY: cells whose
closure collapses onto the resonance lattice (Lambert returns no requested-rev
solution, or the V_inf-magnitude continuity is a degenerate phase-locked
artefact) are excluded from the near-miss band. The Laplace resonance will
phase-lock many Galilean cells; the screening must filter them.

Output: ``data/scan_433_jupiter_galilean.jsonl`` (one leading ``_meta`` row,
one row per closed cycle below the near-miss band, one trailing ``summary``
row -- same schema as the #320 outputs).

NO catalogue writeback. SILVER candidates feed a #306-style follow-up task.

Discipline:
  * READ-ONLY on substrate modules.
  * Sourced golden discipline: all moon/system GM/element values from JPL
    DE440 / SAT441 (see ``src/cyclerfinder/core/satellites.py`` sourcing).
  * The JSONL is a clean negative or a flagged SILVER with provenance.

Run the full sweep as::

    uv run python scripts/scan_433_jupiter_galilean.py

Run the tiny end-to-end smoke as::

    uv run python scripts/scan_433_jupiter_galilean.py --smoke
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

# Module constants so the controller can drive the full sweep without editing
# the body. #433 deltas vs #320: finer phase grid + wider n_rev window.
PRIMARY = "Jupiter"
MOON_SET: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto", "Amalthea")
NREV_GRID: tuple[int, ...] = (0, 1, 2, 3, 4, 5)
N_PHASE = 48
N_OFFSET = 48
TOF_SCALES: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)
NEAR_MISS_BAND_KMS = 1.0
SILVER_GATE_KMS = 0.05
OUT_PATH = ROOT / "data" / "scan_433_jupiter_galilean.jsonl"


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


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
    INCLUDING the closed-cycle anchor wrap (matches the #320 / production genome).

    The ``feasible=False`` return is the resonance/irrationality screen: a
    requested-rev Lambert solution that does not exist on the resonance lattice
    is excluded (matches #320 ``_close_one``).
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
    # Closed-cycle anchor wrap (matches #285 fix C).
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
    n_phase: int = N_PHASE,
    n_offset: int = N_OFFSET,
    tof_scales: tuple[float, ...] = TOF_SCALES,
) -> dict[str, Any]:
    """Sweep global phase x relative offset on a closed length-3 cycle.

    The cycle has two distinct moons (anchor + 1 intermediate). The relative
    offset = phase(intermediate) - phase(anchor) is swept on ``n_offset``
    samples; the absolute anchor phase is swept on ``n_phase`` samples. This
    is the #312/#344 finer-phase basin probe -- raised to 48x48 for #433 so
    the Laplace-resonant Galilean cells are not aliased.

    Returns the best basin-floor record (residual + diagnostics).
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
    """Enumerate closed length-3 repeated-moon cycles over ``moons``.

    Returns ``(anchor, intermediate, anchor)`` for every ordered distinct pair.
    For the 5-moon Galilean set this is 5*4 = 20 closed cycles (same scaling
    rule as the #320 Saturn 5-moon enumeration).
    """
    cycles: list[tuple[str, ...]] = []
    for a in moons:
        for b in moons:
            if a == b:
                continue
            cycles.append((a, b, a))
    return cycles


def _per_system_sweep(
    *,
    primary: str,
    moons: tuple[str, ...],
    nrev_grid: tuple[int, ...],
    out_path: Path,
    sha: str,
    near_miss_band_kms: float = NEAR_MISS_BAND_KMS,
    silver_gate_kms: float = SILVER_GATE_KMS,
    n_phase: int = N_PHASE,
    n_offset: int = N_OFFSET,
    tof_scales: tuple[float, ...] = TOF_SCALES,
) -> dict[str, Any]:
    """Sweep one system's closed length-3 cycles over (n_rev1, n_rev2) x basin."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mu = PRIMARIES[primary]

    # Resolve epoch-aware constants; drop any moon that lacks the needed
    # elements rather than fabricating them (sourced golden discipline).
    consts: dict[str, tuple[float, float]] = {}
    usable_moons: list[str] = []
    for m in moons:
        sat = SATELLITES.get(m)
        if sat is None or sat.sma_km is None:
            print(
                f"[433 {primary}] DROP {m}: missing satellite element record",
                flush=True,
            )
            continue
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))
        usable_moons.append(m)
    moons = tuple(usable_moons)

    cycles = _closed_length3_cycles(moons)
    n_cells = len(cycles) * len(nrev_grid) * len(nrev_grid)
    print(
        f"[433 {primary}] moons={moons} cycles={len(cycles)} "
        f"n_rev_grid={nrev_grid} phase={n_phase}x{n_offset} cells={n_cells}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    near_miss: list[dict[str, Any]] = []
    silver: list[dict[str, Any]] = []
    # Per-cycle best tracking for the end-of-run report (keyed by tuple).
    per_cycle_best: dict[tuple[str, ...], dict[str, Any]] = {}
    n_resonance_screened = 0
    t0 = time.time()

    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": f"#433 Jupiter Galilean repeated-moon quasi_cycler sweep at {primary}",
            "primary": primary,
            "moons": list(moons),
            "n_rev_grid": list(nrev_grid),
            "n_phase": n_phase,
            "n_offset": n_offset,
            "tof_scales": list(tof_scales),
            "silver_gate_kms": silver_gate_kms,
            "near_miss_band_kms": near_miss_band_kms,
            "system_mu_source": "JPL DE440 SAT441 / src/cyclerfinder/core/satellites.py",
            "deltas_vs_320": "phase 48x48 (was 24x24); n_rev 0..5 (was 0..3)",
            "git_sha": sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        idx = 0
        for cycle in cycles:
            for nr1 in nrev_grid:
                for nr2 in nrev_grid:
                    if nr1 == 0 and nr2 == 0:
                        continue  # exclude trivial direct-transfer corner (#285 fix B)
                    nrevs = (nr1, nr2)
                    best = _sweep_one_cycle(
                        seq=cycle,
                        nrevs=nrevs,
                        consts=consts,
                        mu=mu,
                        n_phase=n_phase,
                        n_offset=n_offset,
                        tof_scales=tof_scales,
                    )
                    if not math.isfinite(best["residual_kms"]):
                        # No requested-rev Lambert solution anywhere on the
                        # phase x offset x tof lattice: resonance/irrationality
                        # screen rejected the whole cell (the Laplace lock).
                        n_resonance_screened += 1
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
                    # Physical-sanity gate (#324) at every encounter.
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
                    # Literature-novelty check (#328) -- offline anchor overlap.
                    sig = CandidateSignature(
                        primary=primary,
                        sequence=tuple(cycle),
                        vinf_per_encounter_kms=vinfs_t,
                        n_rev=tuple(int(x) for x in nrevs),
                    )
                    anchors = _candidate_anchors(sig)
                    row["n_lit_corpus_anchors_overlap"] = len(anchors)
                    row["lit_corpus_anchor_names"] = [a.name for a in anchors][:3]
                    # SILVER status: residual under gate + physical gate passed.
                    # Lit-fresh is COMPLEMENTARY signal (a populated anchor means
                    # the family is published; not-found via offline anchors is
                    # NECESSARY-not-sufficient for novelty).
                    is_silver = bool(res < silver_gate_kms and row["physical_gate_passed"])
                    row["is_silver"] = is_silver
                    is_lit_fresh_class = bool(len(anchors) == 0)
                    row["is_lit_fresh_class"] = is_lit_fresh_class

                    # Track per-cycle best (over all n_rev cells) for the report.
                    pcb = per_cycle_best.get(cycle)
                    if pcb is None or res < pcb["residual_kms"]:
                        per_cycle_best[cycle] = row

                    if res < near_miss_band_kms:
                        fh.write(json.dumps(row) + "\n")
                        fh.flush()
                        rows.append(row)
                        near_miss.append(row)
                        if is_silver:
                            silver.append(row)
                            print(
                                f"[433 {primary}] SILVER seq={cycle} n_rev={nrevs} "
                                f"res={res:.4f} km/s vinfs={vinfs_t} "
                                f"physical={row['physical_gate_passed']} "
                                f"anchors={row['n_lit_corpus_anchors_overlap']}",
                                flush=True,
                            )
                    idx += 1

        near_miss.sort(key=lambda r: r["residual_kms"])
        silver.sort(key=lambda r: r["residual_kms"])
        elapsed = time.time() - t0

        # Per-cycle best report (every cycle, not just near-miss ones).
        per_cycle_report = [
            {
                "sequence": list(cyc),
                "best_residual_kms": rec["residual_kms"],
                "n_rev": rec["n_rev"],
                "dv_band": _dv_band(float(rec["residual_kms"])),
                "physical_gate_passed": rec["physical_gate_passed"],
                "n_lit_corpus_anchors_overlap": rec["n_lit_corpus_anchors_overlap"],
                "is_lit_fresh_class": rec["is_lit_fresh_class"],
            }
            for cyc, rec in sorted(per_cycle_best.items(), key=lambda kv: kv[1]["residual_kms"])
        ]

        summary = {
            "_meta": True,
            "kind": "summary",
            "primary": primary,
            "n_cells_evaluated": idx,
            "n_resonance_screened": n_resonance_screened,
            "n_below_near_miss_band": len(near_miss),
            "n_silver": len(silver),
            "best_residual_kms": (rows[0]["residual_kms"] if rows else None),
            "best_record_brief": (
                {
                    "sequence": rows[0]["sequence"],
                    "n_rev": rows[0]["n_rev"],
                    "residual_kms": rows[0]["residual_kms"],
                    "vinf_per_encounter_kms": rows[0]["vinf_per_encounter_kms"],
                    "is_silver": rows[0]["is_silver"],
                    "is_lit_fresh_class": rows[0]["is_lit_fresh_class"],
                }
                if rows
                else None
            ),
            "per_cycle_best": per_cycle_report,
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
            "git_sha": sha,
            "discipline": (
                "SILVER here is residual < silver_gate_kms AND physical-gate "
                "passed; offline lit-corpus anchor-overlap = 0 is "
                "NECESSARY-not-sufficient for novelty (see literature_check "
                "module docstring). No catalogue writeback."
            ),
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()

        print(
            f"[433 {primary}] DONE -- {idx} cells "
            f"({n_resonance_screened} resonance-screened), {len(near_miss)} "
            f"near-miss, {len(silver)} SILVER, "
            f"best={summary['best_residual_kms']!r} in {elapsed:.1f}s",
            flush=True,
        )
        print(f"[433 {primary}] per-cycle best residual + dv_band + gate + lit:", flush=True)
        for rec in per_cycle_report:
            print(
                f"  seq={rec['sequence']} n_rev={rec['n_rev']} "
                f"res={rec['best_residual_kms']:.4f} km/s band={rec['dv_band']} "
                f"physical={'PASS' if rec['physical_gate_passed'] else 'FAIL'} "
                f"lit_anchors={rec['n_lit_corpus_anchors_overlap']} "
                f"lit_fresh={rec['is_lit_fresh_class']}",
                flush=True,
            )
        return summary


def _dv_band(residual_kms: float) -> str:
    """Coarse closure-residual band label for the end-of-run report."""
    if residual_kms < 0.05:
        return "silver_gate"
    if residual_kms < 0.1:
        return "sub_0.1"
    if residual_kms < 0.5:
        return "sub_0.5"
    if residual_kms < 1.0:
        return "near_miss"
    return "above_band"


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    smoke = "--smoke" in argv
    sha = _git_sha()

    if smoke:
        # Tiny end-to-end config: ONE Galilean 3-cycle, 8x8 phase, n_rev (0,1,2).
        print(f"[433] SMOKE -- Jupiter Galilean tiny end-to-end -- sha={sha}", flush=True)
        moons = ("Io", "Europa")
        out_path = ROOT / "data" / "scan_433_jupiter_galilean_smoke.jsonl"
        summary = _per_system_sweep(
            primary=PRIMARY,
            moons=moons,
            nrev_grid=(0, 1, 2),
            out_path=out_path,
            sha=sha,
            n_phase=8,
            n_offset=8,
        )
        print(
            f"[433] SMOKE best-closure residual = "
            f"{summary['best_residual_kms']!r} km/s "
            f"(resonance-screened {summary['n_resonance_screened']} cells)",
            flush=True,
        )
        return 0

    print(f"[433] Jupiter Galilean repeated-moon sweep -- sha={sha}", flush=True)
    summary = _per_system_sweep(
        primary=PRIMARY,
        moons=MOON_SET,
        nrev_grid=NREV_GRID,
        out_path=OUT_PATH,
        sha=sha,
    )
    print(
        f"[433] cross-system summary: {summary['primary']}: "
        f"cells={summary['n_cells_evaluated']}, "
        f"resonance_screened={summary['n_resonance_screened']}, "
        f"near_miss={summary['n_below_near_miss_band']}, "
        f"silver={summary['n_silver']}, best={summary['best_residual_kms']!r}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
