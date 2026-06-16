"""#344 Phase 1 Part B -- other Saturn-system pair / 3-body sweeps.

The brief enumerates five candidate Part B configurations:

  * Titan-Iapetus  (Iapetus highly inclined; NEW to the Vector B set --
    Iapetus is not in #320's moon set)
  * Titan-Dione   (in #320 at coarse (0..3) grid -- revisit at (0..5))
  * Titan-Tethys  (in #320 at coarse (0..3) grid -- revisit at (0..5))
  * Rhea-Dione    (#285's 0.107 km/s ceiling; revisit at finer methodology)
  * 3-body Titan-Rhea-Dione length-5 cycles (Takubo 2210.14996 prior, in
    KNOWN_CORPUS via Brinckerhoff/Lo/Marsden + Cassini anchors)

All use the same ``_sweep_one_cycle`` (basin-floor) convention as Part A
and #320 Vector B, so residuals are directly comparable to the Saturn
SILVER baseline (Titan-Rhea-Titan (1,1) at 0.0316 km/s).

NO catalogue writeback. NO novelty claims. The post-#334 KNOWN_CORPUS
anchor overlap is recorded per row; the verdict doc analyzes match depth
(body-set / topology / V_inf).

Sourced moon GMs from JPL DE440 / sat441 mean elements via
``src/cyclerfinder/core/satellites.py`` (accessed 2026-06-14):
  * Tethys GM 41.21, a 294670 km
  * Dione GM 73.116, a 377420 km
  * Rhea GM 153.94, a 527070 km
  * Titan GM 8978.14, a 1221870 km
  * Iapetus GM 120.51511, a 3561700 km

Iapetus is genuinely interesting because (a) it's outside the #320 set
(novel coverage), (b) its mean inclination relative to Saturn's
equatorial plane is ~15 deg (much higher than the inner regular
satellites), (c) its semi-major axis is ~3x Titan's, giving a long
period and large V_inf differentials. For the patched-conic Lambert
genome here that uses circular-coplanar moon orbits, the inclination
is averaged out; but the Iapetus inclination is a flag for Part C
follow-ups in any 3D extension.

Output: ``data/scan_344_saturn_other_pairs_<pair>.jsonl`` per pair, plus
a top-level index ``data/scan_344_saturn_other_pairs_index.jsonl``.

Run as::

    uv run python scripts/scan_344_saturn_other_pairs.py
"""

from __future__ import annotations

import itertools
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
# #320 Vector B closure-residual machinery (direct port -- same as Part A).
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
    """Compute the closure residual of one (sequence, n_rev, phasing, tof_scale)."""
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
    """Sweep global phase x relative offset on a closed length-N cycle.

    Generalised from the length-3 form to support length-5 3-body cycles
    (Titan-Rhea-Dione-Rhea-Titan etc.). The relative-offset sweep ranges
    over the SECOND moon (lowest-index non-anchor), and any additional
    distinct moons are phased proportionally so the per-moon initial-true-
    anomaly is consistent across the cycle. This is the same convention
    as ``_sweep_one_cycle`` in scan_320, which assumed exactly two
    distinct moons; here we keep the two-moon assumption for length-3 and
    extend to a coupled-offset form for the 3-body case.

    For 3-body cycles we set ``theta[moon_i] = phase0 + (i-1) * rel_off``,
    i.e. the second moon is at ``+rel_off`` and the third at ``+2*rel_off``.
    This couples the second-third phasing through a single offset
    parameter, which is the same dimensionality as the length-3 basin --
    it can MISS configurations where the second and third moons have
    independent offsets, but it keeps the cost the same order as Part A.
    """
    distinct = sorted(set(seq))
    anchor = seq[0]
    others = [m for m in distinct if m != anchor]

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
            theta = {anchor: phase0}
            for k_idx, m in enumerate(others, start=1):
                theta[m] = phase0 + k_idx * rel_off
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


def _closed_length5_cycles_3body(moons: tuple[str, ...]) -> list[tuple[str, ...]]:
    """Enumerate closed length-5 cycles on a 3-moon set.

    A closed length-5 cycle is ``(a, x, y, z, a)`` with: ``a != x``,
    ``x != y``, ``y != z``, ``z != a``. ``x``, ``y``, ``z`` are drawn from
    the moon set. The body_set used in the lit-check is the union of all
    distinct moons in the sequence.
    """
    seqs: list[tuple[str, ...]] = []
    for combo in itertools.product(moons, repeat=5):
        if combo[0] != combo[-1]:
            continue
        if any(combo[i] == combo[i + 1] for i in range(4)):
            continue
        if len(set(combo)) < 2:
            continue
        seqs.append(combo)
    return seqs


# ---------------------------------------------------------------------------
# Per-pair sweep (length-3, same as Part A but generic over moon pair).
# ---------------------------------------------------------------------------


def run_pair_sweep(
    *,
    out_path: Path,
    git_sha: str,
    moons: tuple[str, str],
    pair_label: str,
    nrev_grid: tuple[int, ...] = (0, 1, 2, 3, 4, 5),
    n_phase: int = 24,
    n_offset: int = 24,
    silver_gate_kms: float = 0.05,
    near_miss_band_kms: float = 1.0,
) -> dict[str, Any]:
    """Sweep closed length-3 cycles on ``moons`` at the wider (k1, k2) grid."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    primary = "Saturn"
    mu = PRIMARIES[primary]
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    cycles = _closed_length3_cycles(moons)
    cells_per_seq = len(nrev_grid) ** 2 - 1
    n_cells = len(cycles) * cells_per_seq
    print(
        f"[344-B] {pair_label} moons={moons} cycles={len(cycles)} "
        f"n_rev_grid={nrev_grid} cells={n_cells}",
        flush=True,
    )

    near_miss: list[dict[str, Any]] = []
    silver: list[dict[str, Any]] = []
    t0 = time.time()

    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": f"#344 Phase 1 Part B -- Saturn {pair_label} wider (k1,k2) grid",
            "primary": primary,
            "moons": list(moons),
            "n_rev_grid": list(nrev_grid),
            "n_phase": n_phase,
            "n_offset": n_offset,
            "tof_scales": [0.5, 1.0, 1.5, 2.0],
            "silver_gate_kms": silver_gate_kms,
            "near_miss_band_kms": near_miss_band_kms,
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
                    row["lit_corpus_anchor_names"] = [a.name for a in anchors][:5]
                    is_silver = bool(res < silver_gate_kms and row["physical_gate_passed"])
                    row["is_silver"] = is_silver
                    is_lit_fresh_class = bool(len(anchors) == 0)
                    row["is_lit_fresh_class"] = is_lit_fresh_class
                    if res < near_miss_band_kms:
                        fh.write(json.dumps(row) + "\n")
                        fh.flush()
                        near_miss.append(row)
                        if is_silver:
                            silver.append(row)
                            print(
                                f"[344-B] SILVER {pair_label} seq={cycle} n_rev={nrevs} "
                                f"res={res:.4f} km/s vinfs={vinfs_t} "
                                f"physical={row['physical_gate_passed']} "
                                f"anchors={row['n_lit_corpus_anchors_overlap']} "
                                f"lit_fresh={is_lit_fresh_class}",
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
            "pair_label": pair_label,
            "n_cells_evaluated": idx,
            "n_below_near_miss_band": len(near_miss),
            "n_silver": len(silver),
            "best_residual_kms": (near_miss[0]["residual_kms"] if near_miss else None),
            "best_record_brief": (
                {
                    "sequence": near_miss[0]["sequence"],
                    "n_rev": near_miss[0]["n_rev"],
                    "residual_kms": near_miss[0]["residual_kms"],
                    "vinf_per_encounter_kms": near_miss[0]["vinf_per_encounter_kms"],
                    "is_silver": near_miss[0]["is_silver"],
                    "is_lit_fresh_class": near_miss[0]["is_lit_fresh_class"],
                    "n_lit_corpus_anchors_overlap": near_miss[0]["n_lit_corpus_anchors_overlap"],
                    "lit_corpus_anchor_names": near_miss[0]["lit_corpus_anchor_names"],
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
                    "n_lit_corpus_anchors_overlap": r["n_lit_corpus_anchors_overlap"],
                }
                for r in near_miss[:5]
            ],
            "elapsed_s": elapsed,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()
        print(
            f"[344-B] {pair_label} DONE -- {idx} cells, {len(near_miss)} near-miss, "
            f"{len(silver)} SILVER, best={summary['best_residual_kms']!r} "
            f"in {elapsed:.1f}s",
            flush=True,
        )
    return summary


# ---------------------------------------------------------------------------
# Length-5 3-body sweep (Titan-Rhea-Dione).
# ---------------------------------------------------------------------------


def run_3body_sweep(
    *,
    out_path: Path,
    git_sha: str,
    moons: tuple[str, str, str],
    label: str,
    nrev_grid: tuple[int, ...] = (0, 1, 2, 3),
    n_phase: int = 12,
    n_offset: int = 12,
    silver_gate_kms: float = 0.05,
    near_miss_band_kms: float = 1.0,
) -> dict[str, Any]:
    """Sweep closed length-5 cycles on a 3-moon Saturn set.

    Length-5 cycles have 4 legs, so (k1..k4) over (0..3) is 4^4 = 256
    cells per sequence; capped at 3 to keep cost comparable to Part A's
    35 cells per length-3 sequence at (0..5). The basin uses a coarser
    n_phase = n_offset = 12 for the same reason -- ratio of 12^2 * 256 /
    24^2 * 35 = ~2.6x the per-sequence cost.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    primary = "Saturn"
    mu = PRIMARIES[primary]
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    cycles = _closed_length5_cycles_3body(moons)
    cells_per_seq = len(nrev_grid) ** 4 - 1  # exclude trivial (0,0,0,0)
    n_cells = len(cycles) * cells_per_seq
    print(
        f"[344-B-3b] {label} moons={moons} cycles={len(cycles)} "
        f"n_rev_grid={nrev_grid} cells={n_cells}",
        flush=True,
    )

    near_miss: list[dict[str, Any]] = []
    silver: list[dict[str, Any]] = []
    t0 = time.time()

    with out_path.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": f"#344 Phase 1 Part B -- Saturn 3-body {label} length-5 cycles",
            "primary": primary,
            "moons": list(moons),
            "seq_length": 5,
            "n_rev_grid": list(nrev_grid),
            "n_phase": n_phase,
            "n_offset": n_offset,
            "tof_scales": [0.5, 1.0, 1.5, 2.0],
            "silver_gate_kms": silver_gate_kms,
            "near_miss_band_kms": near_miss_band_kms,
            "phase_convention": (
                "coupled-offset: theta[moon_i] = phase0 + i * rel_off for "
                "i in 1..N-1 of distinct moons. Limited basin dimensionality "
                "for cost; may MISS configurations with independent offsets."
            ),
            "system_mu_source": "JPL DE440 / src/cyclerfinder/core/satellites.py",
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        idx = 0
        for cycle in cycles:
            for nrevs in itertools.product(nrev_grid, repeat=4):
                if all(nr == 0 for nr in nrevs):
                    idx += 1
                    continue
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
                vinfs_t: tuple[float, ...] = tuple(float(v) for v in best["vinf_per_encounter_kms"])
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
                row["lit_corpus_anchor_names"] = [a.name for a in anchors][:5]
                is_silver = bool(res < silver_gate_kms and row["physical_gate_passed"])
                row["is_silver"] = is_silver
                is_lit_fresh_class = bool(len(anchors) == 0)
                row["is_lit_fresh_class"] = is_lit_fresh_class
                if res < near_miss_band_kms:
                    fh.write(json.dumps(row) + "\n")
                    fh.flush()
                    near_miss.append(row)
                    if is_silver:
                        silver.append(row)
                        print(
                            f"[344-B-3b] SILVER {label} seq={cycle} n_rev={nrevs} "
                            f"res={res:.4f} km/s anchors={row['n_lit_corpus_anchors_overlap']}",
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
            "label": label,
            "n_cells_evaluated": idx,
            "n_below_near_miss_band": len(near_miss),
            "n_silver": len(silver),
            "best_residual_kms": (near_miss[0]["residual_kms"] if near_miss else None),
            "best_record_brief": (
                {
                    "sequence": near_miss[0]["sequence"],
                    "n_rev": near_miss[0]["n_rev"],
                    "residual_kms": near_miss[0]["residual_kms"],
                    "vinf_per_encounter_kms": near_miss[0]["vinf_per_encounter_kms"],
                    "is_silver": near_miss[0]["is_silver"],
                    "is_lit_fresh_class": near_miss[0]["is_lit_fresh_class"],
                    "n_lit_corpus_anchors_overlap": near_miss[0]["n_lit_corpus_anchors_overlap"],
                    "lit_corpus_anchor_names": near_miss[0]["lit_corpus_anchor_names"],
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
                    "n_lit_corpus_anchors_overlap": r["n_lit_corpus_anchors_overlap"],
                }
                for r in near_miss[:5]
            ],
            "elapsed_s": elapsed,
            "git_sha": git_sha,
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()
        print(
            f"[344-B-3b] {label} DONE -- {idx} cells, {len(near_miss)} near-miss, "
            f"{len(silver)} SILVER, best={summary['best_residual_kms']!r} "
            f"in {elapsed:.1f}s",
            flush=True,
        )
    return summary


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


PAIRS: list[tuple[str, str, tuple[str, str]]] = [
    ("titan_iapetus", "Titan-Iapetus", ("Titan", "Iapetus")),
    ("titan_dione", "Titan-Dione", ("Titan", "Dione")),
    ("titan_tethys", "Titan-Tethys", ("Titan", "Tethys")),
    ("rhea_dione", "Rhea-Dione", ("Rhea", "Dione")),
]


def main() -> int:
    sha = _git_sha()
    print(f"[344-B] Saturn other-pair extended sweep -- sha={sha}", flush=True)

    index_rows: list[dict[str, Any]] = []
    for slug, label, moons in PAIRS:
        out_path = ROOT / "data" / f"scan_344_saturn_{slug}.jsonl"
        print(f"[344-B] {label} -- out={out_path}", flush=True)
        summary = run_pair_sweep(
            out_path=out_path,
            git_sha=sha,
            moons=moons,
            pair_label=label,
        )
        index_rows.append(
            {
                "slug": slug,
                "label": label,
                "moons": list(moons),
                "best_residual_kms": summary["best_residual_kms"],
                "n_silver": summary["n_silver"],
                "best_record_brief": summary["best_record_brief"],
                "path": str(out_path.relative_to(ROOT)),
            }
        )

    # 3-body length-5 Titan-Rhea-Dione.
    out_3b = ROOT / "data" / "scan_344_saturn_titan_rhea_dione_3body.jsonl"
    summary_3b = run_3body_sweep(
        out_path=out_3b,
        git_sha=sha,
        moons=("Titan", "Rhea", "Dione"),
        label="Titan-Rhea-Dione (length-5, 3-body)",
    )
    index_rows.append(
        {
            "slug": "titan_rhea_dione_3body",
            "label": "Titan-Rhea-Dione (length-5, 3-body)",
            "moons": ["Titan", "Rhea", "Dione"],
            "best_residual_kms": summary_3b["best_residual_kms"],
            "n_silver": summary_3b["n_silver"],
            "best_record_brief": summary_3b["best_record_brief"],
            "path": str(out_3b.relative_to(ROOT)),
        }
    )

    # Top-level index.
    out_index = ROOT / "data" / "scan_344_saturn_other_pairs_index.jsonl"
    with out_index.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "index",
                    "task": "#344 Phase 1 Part B -- Saturn other-pair extended sweep",
                    "systems": index_rows,
                    "git_sha": sha,
                }
            )
            + "\n"
        )
    print(f"[344-B] DONE -- index: {out_index}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
