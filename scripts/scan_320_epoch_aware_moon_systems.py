"""#320 Vector B -- epoch-aware repeated-moon sweep at non-EM systems.

Phase 1 of #320 (first systematic quasi_cycler discovery sweep). For each of
three candidate non-EM moon systems:

  * Saturn Rhea-Dione  -- revisit the #285/#311 0.107 km/s near-miss with the
    finer (k1, k2) grid at higher phase resolution and the relative-offset
    sweep added (analogous to #312's Uranus extended sweep that surfaced
    Umbriel-Oberon (1,1) at 0.025 km/s -- the candidate that became #339).
  * Neptune Triton-Proteus  -- a speculative scout (no published cycler
    corpus for the Neptune system).
  * Pluto Charon-Hydra  -- similar topology to Umbriel-Oberon but at lower
    primary mu (Pluto system GM is ~6 orders below Saturn's).

For each system, the sweep enumerates closed repeated-moon length-3 cycles
with per-leg revolution counts in ``[0, 5]`` and finds the basin-floor
residual via a 24-sample relative-offset sweep (per the #312 lesson:
sweeping the relative offset between moon pairs reveals deeper basins than
the global-phase grid alone). Surviving candidates with residual below the
0.05 km/s gate are run through the #324 physical-sanity gate and the #328
literature-novelty check.

Output: ``data/scan_320_epoch_aware_<system>.jsonl`` per system (one leading
``_meta`` row, one row per closed cycle below the near-miss band, one trailing
``summary`` row).

NO catalogue writeback. SILVER candidates feed a #306-style follow-up task.

Discipline:
  * READ-ONLY on substrate modules.
  * Sourced golden discipline: all moon/system GM values from JPL DE440
    (see ``src/cyclerfinder/core/satellites.py`` sourcing comments).
  * Each per-system JSONL is a clean negative or a flagged SILVER with
    provenance; novelty claims defer to #328's KNOWN_CORPUS check.

Run as::

    uv run python scripts/scan_320_epoch_aware_moon_systems.py
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
    n_phase: int = 24,
    n_offset: int = 24,
    tof_scales: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0),
) -> dict[str, Any]:
    """Sweep global phase x relative offset on a closed length-3 cycle.

    The cycle has two distinct moons (anchor + 1 intermediate). The relative
    offset = phase(intermediate) - phase(anchor) is swept on ``n_offset``
    samples; the absolute anchor phase is swept on ``n_phase`` samples.

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
    near_miss_band_kms: float = 1.0,
    silver_gate_kms: float = 0.05,
    n_phase: int = 24,
    n_offset: int = 24,
) -> dict[str, Any]:
    """Sweep one system's closed length-3 cycles over (n_rev1, n_rev2) x basin."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mu = PRIMARIES[primary]
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    cycles = _closed_length3_cycles(moons)
    n_cells = len(cycles) * len(nrev_grid) * len(nrev_grid)
    print(
        f"[320-B {primary}] moons={moons} cycles={len(cycles)} "
        f"n_rev_grid={nrev_grid} cells={n_cells}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    near_miss: list[dict[str, Any]] = []
    silver: list[dict[str, Any]] = []
    t0 = time.time()

    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": f"#320 Vector B -- epoch-aware repeated-moon sweep at {primary}",
            "primary": primary,
            "moons": list(moons),
            "n_rev_grid": list(nrev_grid),
            "n_phase": n_phase,
            "n_offset": n_offset,
            "tof_scales": [0.5, 1.0, 1.5, 2.0],
            "silver_gate_kms": silver_gate_kms,
            "near_miss_band_kms": near_miss_band_kms,
            "system_mu_source": "JPL DE440 / src/cyclerfinder/core/satellites.py",
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
                    if res < near_miss_band_kms:
                        fh.write(json.dumps(row) + "\n")
                        fh.flush()
                        rows.append(row)
                        near_miss.append(row)
                        if is_silver:
                            silver.append(row)
                            print(
                                f"[320-B {primary}] SILVER seq={cycle} n_rev={nrevs} "
                                f"res={res:.4f} km/s vinfs={vinfs_t} "
                                f"physical={row['physical_gate_passed']} "
                                f"anchors={row['n_lit_corpus_anchors_overlap']}",
                                flush=True,
                            )
                    idx += 1

        # Sort near_miss + silver by residual.
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
            f"[320-B {primary}] DONE -- {idx} cells, {len(near_miss)} near-miss, "
            f"{len(silver)} SILVER, best={summary['best_residual_kms']!r} "
            f"in {elapsed:.1f}s",
            flush=True,
        )
        return summary


def main() -> int:
    sha = _git_sha()
    print(f"[320-B] non-EM moon-system sweeps -- sha={sha}", flush=True)

    targets = [
        # Saturn Rhea-Dione: the #285/#311 near-miss at 0.107 km/s ballistic
        # ceiling (genome confirmed per #311 robustness sweep). Include the
        # adjacent Tethys/Titan/Enceladus options so the (k1,k2) grid surfaces
        # any deeper basins among the regular moons.
        {
            "primary": "Saturn",
            "moons": ("Enceladus", "Tethys", "Dione", "Rhea", "Titan"),
            "nrev_grid": (0, 1, 2, 3),
            "tag": "saturn",
        },
        # Neptune: Triton (large, retrograde) + Proteus (small, prograde). The
        # GM ratio Triton/Proteus is ~550 so a Proteus flyby contributes negligible
        # bending; only Triton-only inner ring cycles are physically interesting,
        # but there's only ONE prograde sibling for triple-cycle topology, so the
        # length-3 closed cycle is just Triton-Proteus-Triton.
        {
            "primary": "Neptune",
            "moons": ("Triton", "Proteus"),
            "nrev_grid": (0, 1, 2, 3),
            "tag": "neptune",
        },
        # Pluto: Charon-Hydra-Charon and Charon-Nix-Charon. Nix/Hydra have GM
        # ~0.002 (negligible bending); the physical-sanity gate will reject
        # cycles that try to use them as flyby bodies. This is the speculative
        # leg of Vector B -- expected to be clean negatives on the physical
        # gate but the closure-residual surface is recorded for completeness.
        {
            "primary": "Pluto",
            "moons": ("Charon", "Hydra", "Nix"),
            "nrev_grid": (0, 1, 2, 3),
            "tag": "pluto",
        },
    ]

    summaries: list[dict[str, Any]] = []
    for t in targets:
        out_path = ROOT / "data" / f"scan_320_epoch_aware_{t['tag']}.jsonl"
        s = _per_system_sweep(
            primary=t["primary"],
            moons=t["moons"],
            nrev_grid=t["nrev_grid"],
            out_path=out_path,
            sha=sha,
        )
        summaries.append(s)

    print("[320-B] cross-system summary:", flush=True)
    for s in summaries:
        print(
            f"  {s['primary']}: cells={s['n_cells_evaluated']}, near_miss="
            f"{s['n_below_near_miss_band']}, silver={s['n_silver']}, "
            f"best={s['best_residual_kms']!r}",
            flush=True,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
