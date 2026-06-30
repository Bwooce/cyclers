"""#501 — Broadened real-ephemeris Galilean joint-search campaign.

Sweeps six Galilean moon-tour sequences beyond the CGCEC Phase-2b smoke (#318),
over a wider epoch window (2030-2040) and larger Sobol budget (512 cells per
sequence = 3072 total).  The CGCEC smoke found 0/256 in 2033-2035; this run
explores the wider sequence / epoch landscape.

**Sequences surveyed (all closed tours, seq[0]==seq[-1]):**

=======  ==========================================  ======  =====
Tag      Sequence                                    Legs    ToF/leg
=======  ==========================================  ======  =====
EGE      Europa-Ganymede-Europa                      2       5-40 d
GCG      Ganymede-Callisto-Ganymede                  2       10-50 d
EGCE     Europa-Ganymede-Callisto-Europa             3       5-40 d
IEI      Io-Europa-Io                                2       3-20 d
IEGI     Io-Europa-Ganymede-Io                       3       3-30 d
EGCGE    Europa-Ganymede-Callisto-Ganymede-Europa    4       10-50 d
=======  ==========================================  ======  =====

**Pipeline (per sequence):**

1. Positive control (Liang Member D, CGCEC) - abort entire run if prefilter fails.
2. Sobol sample N_SAMPLES_PER_SEQ cells (epoch + ToF seed + n_revs + branch).
3. Parallel evaluate_joint_cell prefilter (real-eph patched-conic surrogate).
4. N-body shoot top-K_SHOOT_PER_SEQ feasible survivors (jovian_shoot, jacobian=stm).
5. Literature check on any closers (post-hoc, necessary-not-sufficient).
6. Register compute-bounded empty regions per sequence (0 closes -> empty_regions.jsonl).
7. Write verdict note.

**Honest discipline:** EXPECT mostly empty.  "Not found" is a real result -- it maps
the explored region.  Any survivor that BOTH closes (n-body) AND clears lit-novelty
is FLAGGED for human adjudication -- no self-admission, no novelty claim.

Outputs::

    data/scan_501_broadened_<TAG>.jsonl        per-sequence prefilter results
    data/scan_501_shoot_results.jsonl          n-body shoot results (all sequences)
    data/empty_regions.jsonl                   appended (one entry per empty sequence)
    docs/notes/2026-07-01-501-broadened-joint-search-verdict.md
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from cyclerfinder.nbody import jovian  # noqa: E402
from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep  # noqa: E402
from cyclerfinder.search.joint_cell import (  # noqa: E402
    JointCell,
    evaluate_joint_cell,
    liang_member_d_cell,
)
from cyclerfinder.search.joint_sobol import make_sobol_cells  # noqa: E402

# ---------------------------------------------------------------------------
# Campaign configuration
# ---------------------------------------------------------------------------

EPOCH_WINDOW = ("2030-01-01", "2040-01-01")
N_REVS_RANGE = (1, 3)
POWERED_MIN_ALT_KM = 100.0
N_SAMPLES_PER_SEQ = 512
K_SHOOT_PER_SEQ = 5  # n-body shoot top-K survivors per sequence

POSITIVE_CONTROL_DEFECT_THRESHOLD_MS = 1.0

_KERNEL_FALLBACK = str(Path.home() / "dev" / "references" / "kernels" / "jup365.bsp")


@dataclass(frozen=True)
class SequenceConfig:
    tag: str
    sequence: tuple[str, ...]
    tof_seed_range: tuple[float, float]  # d/leg
    n_samples: int = N_SAMPLES_PER_SEQ
    n_revs_range: tuple[int, int] = N_REVS_RANGE


SEQUENCE_CONFIGS: list[SequenceConfig] = [
    SequenceConfig(
        tag="EGE",
        sequence=("Europa", "Ganymede", "Europa"),
        tof_seed_range=(5.0, 40.0),
    ),
    SequenceConfig(
        tag="GCG",
        sequence=("Ganymede", "Callisto", "Ganymede"),
        tof_seed_range=(10.0, 50.0),
    ),
    SequenceConfig(
        tag="EGCE",
        sequence=("Europa", "Ganymede", "Callisto", "Europa"),
        tof_seed_range=(5.0, 40.0),
    ),
    SequenceConfig(
        tag="IEI",
        sequence=("Io", "Europa", "Io"),
        tof_seed_range=(3.0, 20.0),
    ),
    SequenceConfig(
        tag="IEGI",
        sequence=("Io", "Europa", "Ganymede", "Io"),
        tof_seed_range=(3.0, 30.0),
    ),
    SequenceConfig(
        tag="EGCGE",
        sequence=("Europa", "Ganymede", "Callisto", "Ganymede", "Europa"),
        tof_seed_range=(10.0, 50.0),
    ),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

OUT_DIR = REPO_ROOT / "data"
SHOOT_JSONL = OUT_DIR / "scan_501_shoot_results.jsonl"
EMPTY_REGIONS_PATH = OUT_DIR / "empty_regions.jsonl"
VERDICT_NOTE = REPO_ROOT / "docs" / "notes" / "2026-07-01-501-broadened-joint-search-verdict.md"


def _resolve_kernel() -> str | None:
    k = jovian.jup365_kernel_path()
    if k is not None:
        return k
    if Path(_KERNEL_FALLBACK).exists():
        return _KERNEL_FALLBACK
    return None


def _log(msg: str) -> None:
    ts = datetime.now(UTC).strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Top-level picklable process_cell — module level for loky
# ---------------------------------------------------------------------------


def process_cell(cell: JointCell) -> dict[str, Any] | None:
    """Evaluate one Sobol joint-cell via the patched-conic surrogate (loky worker)."""
    kernel = _resolve_kernel()
    if kernel is None:
        return None
    try:
        ephem = jovian.JovianEphemeris(kernel)
        v = evaluate_joint_cell(cell, ephem)
        return {
            "epoch_iso": cell.epoch_iso,
            "sequence": list(cell.sequence),
            "n_revs": list(cell.n_revs),
            "branches": list(cell.branches),
            "tof_seed_days": list(cell.tof_seed_days),
            "powered_min_alt_km": cell.powered_min_alt_km,
            "closure_defect_ms": float(v.closure_defect_ms),
            "feasible": bool(v.feasible),
            "cycle_tof_days": float(v.cycle_tof_days),
            "vinf_kms": [float(x) for x in v.vinf_kms],
            "min_alt_km": float(v.min_alt_km),
            "altitudes_km": [float(x) for x in v.altitudes_km],
            "axis_c_z_km": float(v.axis_c_max_abs_z_km),
        }
    except Exception as exc:
        return {"error": str(exc), "feasible": False, "closure_defect_ms": float("inf")}


# ---------------------------------------------------------------------------
# N-body shoot helpers
# ---------------------------------------------------------------------------


def _build_jovian_seed(
    cell: JointCell,
    ephem: jovian.JovianEphemeris,
) -> Any | None:
    """Build a ShootingSeed from a JointCell's converged patched-conic solution."""
    from cyclerfinder.nbody.jovian import _solve_cycle_legs, optimize_cycle, tdb_sec_from_iso
    from cyclerfinder.nbody.shooter import ShootingSeed

    branch_plan = tuple(zip(cell.n_revs, cell.branches, strict=True))
    t0 = tdb_sec_from_iso(cell.epoch_iso)
    cyc, _ = optimize_cycle(
        t0,
        list(cell.tof_seed_days),
        ephem,
        cycle_index=1,
        vinf_in_prev=None,
        bound_days=cell.bound_days,
        min_alt_km=cell.powered_min_alt_km,
        sequence=cell.sequence,
        branch_plan=branch_plan,
    )
    if not cyc.converged or cyc.cycle_tof_days <= 0.0:
        return None

    legs = _solve_cycle_legs(list(cyc.epochs_sec), ephem, branch_plan, cell.sequence)
    if legs is None:
        return None
    sol_vinf_out, sol_vinf_in = legs
    n_legs = len(cell.sequence) - 1

    node_states = []
    seed_vinf_in = []
    seed_vinf_out = []
    for k in range(n_legs + 1):
        moon = cell.sequence[k]
        t_sec = float(cyc.epochs_sec[k])
        r_m, v_m = ephem.state(moon, t_sec)
        r_m = np.asarray(r_m, dtype=np.float64)
        v_m = np.asarray(v_m, dtype=np.float64)

        vi_out = np.asarray(sol_vinf_out[k] if k < n_legs else sol_vinf_out[0], dtype=np.float64)
        vi_in = np.asarray(
            sol_vinf_in[k - 1] if k > 0 else sol_vinf_in[n_legs - 1], dtype=np.float64
        )
        v_sc = v_m + (vi_in if k == n_legs else vi_out)
        node_states.append(np.concatenate([r_m, v_sc]))
        seed_vinf_in.append(vi_in)
        seed_vinf_out.append(vi_out)

    return ShootingSeed(
        node_states=node_states,
        epochs=[float(t) for t in cyc.epochs_sec],
        tofs=[float(t) for t in cyc.tofs_days],
        sequence=tuple(cell.sequence),
        slack_leg=0,
        period_days=float(cyc.cycle_tof_days),
        vinf_in=seed_vinf_in,
        vinf_out=seed_vinf_out,
    )


def _shoot_cell(
    cell: JointCell,
    kernel: str,
    ephem: jovian.JovianEphemeris,
) -> dict[str, Any]:
    """Build seed + run jovian_shoot(jacobian='stm'); return a verdict dict."""
    t_shoot = time.monotonic()
    seed = _build_jovian_seed(cell, ephem)
    if seed is None:
        return {"converged": False, "error": "seed_build_failed", "wall_s": 0.0}
    try:
        result = jovian.jovian_shoot(seed, kernel_path=kernel, jacobian="stm", max_nfev=40)
        return {
            "converged": bool(result.converged),
            "defect_norm": float(result.defect_norm),
            "seed_defect_norm": float(result.seed_defect_norm),
            "correction_dv_kms": float(result.correction_dv_kms),
            "bend_feasible": bool(result.bend_feasible),
            "vinf_kms": [float(x) for x in result.vinf_per_encounter_kms],
            "n_iterations": int(result.n_iterations),
            "wall_s": float(time.monotonic() - t_shoot),
        }
    except Exception as exc:
        return {
            "converged": False,
            "error": str(exc),
            "wall_s": float(time.monotonic() - t_shoot),
        }


# ---------------------------------------------------------------------------
# Per-sequence sweep
# ---------------------------------------------------------------------------


def _run_sequence(
    cfg: SequenceConfig,
    kernel: str,
    ephem: jovian.JovianEphemeris,
    git_sha: str,
    date_tag: str,
    shoot_fh: Any,
) -> dict[str, Any]:
    """Run the full prefilter + shoot pipeline for one sequence.  Returns a summary dict."""
    tag = cfg.tag
    seq = list(cfg.sequence)
    out_jsonl = OUT_DIR / f"scan_501_broadened_{tag}.jsonl"

    _log(f"[{tag}] building {cfg.n_samples} Sobol cells  sequence={'→'.join(seq)}")
    cells = make_sobol_cells(
        n_samples=cfg.n_samples,
        sequence=seq,
        epoch_window=EPOCH_WINDOW,
        n_revs_range=cfg.n_revs_range,
        tof_seed_range=cfg.tof_seed_range,
        powered_min_alt_km=POWERED_MIN_ALT_KM,
        seed=0,
    )
    _log(f"[{tag}] {len(cells)} cells built")

    _log(f"[{tag}] parallel prefilter sweep (loky, {cfg.n_samples} cells)...")
    t_sweep_start = time.time()
    par_cfg = ParallelSweepConfig(
        n_workers=-1,
        chunk_size=4,
        backend="loky",
        verbose=0,
        raise_on_first_error=False,
    )
    sweep = parallel_sweep(cells, process_cell, config=par_cfg)
    t_sweep_end = time.time()
    elapsed_sweep = t_sweep_end - t_sweep_start

    n_feasible = 0
    feasible_rows: list[tuple[float, int, dict[str, Any], JointCell]] = []
    with out_jsonl.open("w") as fh:
        for i, row in enumerate(sweep.results):
            if row is None:
                continue
            fh.write(json.dumps(row) + "\n")
            if row.get("feasible"):
                n_feasible += 1
                defect = float(row.get("closure_defect_ms", float("inf")))
                feasible_rows.append((defect, i, row, cells[i]))

    feasible_rows.sort(key=lambda x: x[0])
    _log(
        f"[{tag}] sweep done {elapsed_sweep:.1f}s: "
        f"total={cfg.n_samples} succeeded={sweep.n_succeeded} failed={sweep.n_failed} "
        f"feasible={n_feasible}"
    )

    # --- N-body shoot top-K
    n_close = 0
    shoot_rows: list[dict[str, Any]] = []
    top_k = feasible_rows[:K_SHOOT_PER_SEQ]
    t_shoot_start = time.time()
    if not top_k:
        _log(f"[{tag}] no feasible survivors — skip n-body shoot")
    else:
        _log(f"[{tag}] n-body shoot (jacobian=stm) on top {len(top_k)} feasible cells")
        for rank, (defect, idx, prefilter_row, cell) in enumerate(top_k):
            _log(
                f"[{tag}]   shooting rank {rank + 1}/{len(top_k)}: "
                f"epoch={cell.epoch_iso[:10]}  defect={defect:.3f} m/s"
            )
            shoot_result = _shoot_cell(cell, kernel, ephem)
            row: dict[str, Any] = {
                "tag": tag,
                "rank": rank + 1,
                "cell_idx": idx,
                "sequence": list(cell.sequence),
                "epoch_iso": cell.epoch_iso,
                "n_revs": list(cell.n_revs),
                "branches": list(cell.branches),
                "prefilter": prefilter_row,
                "shoot": shoot_result,
            }
            shoot_rows.append(row)
            shoot_fh.write(json.dumps(row) + "\n")
            shoot_fh.flush()
            if shoot_result.get("converged"):
                n_close += 1
                _log(
                    f"[{tag}]   CLOSED: defect_norm={shoot_result.get('defect_norm', '?'):.4e} "
                    f"correction_dv={shoot_result.get('correction_dv_kms', '?'):.4f} km/s"
                )
            else:
                _log(
                    f"[{tag}]   not closed: "
                    + (
                        f"error={shoot_result['error']!r}"
                        if "error" in shoot_result
                        else f"defect_norm={shoot_result.get('defect_norm', '?')}"
                    )
                )
    t_shoot_end = time.time()
    elapsed_shoot = t_shoot_end - t_shoot_start

    summary: dict[str, Any] = {
        "tag": tag,
        "sequence": list(cfg.sequence),
        "n_cells": cfg.n_samples,
        "n_succeeded_prefilter": int(sweep.n_succeeded),
        "n_failed_prefilter": int(sweep.n_failed),
        "n_feasible": n_feasible,
        "n_shot": len(top_k),
        "n_close": n_close,
        "elapsed_sweep_s": round(elapsed_sweep, 1),
        "elapsed_shoot_s": round(elapsed_shoot, 1),
        "feasible_rows": feasible_rows,
        "shoot_rows": shoot_rows,
        "top_k": top_k,
        "git_sha": git_sha,
        "date": date_tag,
    }

    # Register empty region if nothing closes
    if n_close == 0:
        entry = {
            "region_id": (
                f"jovian-{tag.lower()}-sobol-broadened-501-{datetime.now(UTC).strftime('%Y-%m-%d')}"
            ),
            "family": "planet-centric moon system (Jupiter)",
            "centre": "Jupiter",
            "topologies": [
                {
                    "sequence": list(cfg.sequence),
                    "per_leg_revs": list(cfg.n_revs_range),
                    "epoch_window": list(EPOCH_WINDOW),
                }
            ],
            "method_capability": {
                "genome": (
                    "Sobol joint-cell (epoch+tof+n_revs+branch) "
                    "real-eph patched-conic surrogate + jovian_shoot(jacobian=stm)"
                ),
                "corrector": "jovian_shoot (JovianRestrictedNBody + analytic STM Jacobian)",
                "capability_tags": [
                    "real-ephemeris",
                    "multi-rev",
                    "patched-conic-prefilter",
                    "n-body-shoot",
                    "analytic-stm",
                ],
                "git_sha": git_sha,
            },
            "search_extent": {
                "n_sobol_samples": cfg.n_samples,
                "sequence": list(cfg.sequence),
                "epoch_window": list(EPOCH_WINDOW),
                "n_revs_range": list(cfg.n_revs_range),
                "tof_seed_range_days": list(cfg.tof_seed_range),
                "top_k_shot": K_SHOOT_PER_SEQ,
                "ephem_model": "JUP365 real-ephemeris",
                "center": "Jupiter",
            },
            "prune_gates": [
                "evaluate_joint_cell prefilter (feasible + altitude + closure_defect)",
                "jovian_shoot convergence (continuity floor Jones AAS 17-577 §2.5)",
            ],
            "result": {
                "n_cells": cfg.n_samples,
                "n_succeeded_prefilter": int(sweep.n_succeeded),
                "n_failed_prefilter": int(sweep.n_failed),
                "n_feasible": n_feasible,
                "n_shot": len(top_k),
                "n_close": n_close,
                "elapsed_sweep_s": round(elapsed_sweep, 1),
                "elapsed_shoot_s": round(elapsed_shoot, 1),
            },
            "verdict": (
                f"EMPTY -- 0 cells closed under jovian_shoot(jacobian=stm) "
                f"for {tag} in the #501 broadened Sobol budget"
            ),
            "interpretation": (
                f"{cfg.n_samples}-cell Sobol scan of {tag} "
                f"({'→'.join(cfg.sequence)}) over epoch 2030-2040, "
                f"n_revs {cfg.n_revs_range}, ToF {cfg.tof_seed_range} d/leg, "
                f"{n_feasible} feasible prefilter survivors; top-{K_SHOOT_PER_SEQ} "
                f"n-body shot, {n_close} converged. Compute-bounded empty region."
            ),
            "source_anchors": "#501 broadened joint-search discovery campaign",
            "run": {
                "date": date_tag,
                "git_sha": git_sha,
                "elapsed_sweep_s": elapsed_sweep,
                "elapsed_shoot_s": elapsed_shoot,
            },
        }
        with EMPTY_REGIONS_PATH.open("a") as fh:
            fh.write(json.dumps(entry) + "\n")
        _log(f"[{tag}] EMPTY — appended to empty_regions.jsonl")
    else:
        _log(f"[{tag}] {n_close} CLOSE(S) — flagged for human adjudication")

    return summary


# ---------------------------------------------------------------------------
# Literature check
# ---------------------------------------------------------------------------


def _check_lit_for_closers(
    all_shoot_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Run check_literature on any converged n-body closers (post-hoc, necessary-not-sufficient)."""
    lit_results: list[dict[str, Any]] = []
    try:
        from cyclerfinder.search.literature_check import CandidateSignature, check_literature

        try:
            from claude_code_tools import web_search as _ws  # type: ignore[import-not-found]

            def _search_fn(q: str) -> list[dict[str, Any]]:
                return list(_ws(q))

        except ImportError:
            _search_fn = None  # type: ignore[assignment]

        for sr in all_shoot_rows:
            if not sr["shoot"].get("converged"):
                continue
            sig = CandidateSignature(
                primary="Jupiter",
                sequence=tuple(sr["sequence"]),
                n_rev=tuple(sr["n_revs"]),
                resonances=(),
            )
            if _search_fn is not None:
                lit = check_literature(sig, search=_search_fn)
                lit_results.append(
                    {
                        "tag": sr["tag"],
                        "cell_idx": sr["cell_idx"],
                        "epoch_iso": sr["epoch_iso"],
                        "sequence": sr["sequence"],
                        "status": lit.status,
                        "citation": lit.citation,
                    }
                )
                _log(f"  lit_check [{sr['tag']}]: status={lit.status} citation={lit.citation!r}")
            else:
                _log("  WebSearch not available — skipping literature check")
    except Exception as exc:
        _log(f"  literature check skipped: {exc}")
    return lit_results


# ---------------------------------------------------------------------------
# Verdict note
# ---------------------------------------------------------------------------


def _write_verdict_note(
    summaries: list[dict[str, Any]],
    lit_results: list[dict[str, Any]],
    pc_passed: bool,
    date_tag: str,
    git_sha: str,
    total_elapsed_s: float,
) -> None:
    total_cells = sum(s["n_cells"] for s in summaries)
    total_feasible = sum(s["n_feasible"] for s in summaries)
    total_shot = sum(s["n_shot"] for s in summaries)
    total_close = sum(s["n_close"] for s in summaries)
    verdict = "EMPTY (all sequences)" if total_close == 0 else f"{total_close} CLOSE(S)"

    lines = [
        "# #501 Broadened Galilean Joint-Search Verdict (2026-07-01)",
        "",
        f"**Date:** {date_tag[:10]}. **Campaign:** #501 broadened real-eph discovery. "
        f"**Verdict:** {verdict}.",
        "",
        "## Configuration",
        "",
        f"- Epoch window: {EPOCH_WINDOW[0]} to {EPOCH_WINDOW[1]}",
        f"- n_revs range: {N_REVS_RANGE}",
        f"- Powered flyby floor: {POWERED_MIN_ALT_KM} km",
        f"- Cells per sequence: {N_SAMPLES_PER_SEQ} (Sobol, seed=0, scrambled)",
        f"- Top-K shoot per sequence: {K_SHOOT_PER_SEQ}",
        f"- Total cells: {total_cells} ({len(summaries)} sequences x {N_SAMPLES_PER_SEQ})",
        f"- git: `{git_sha}`",
        "",
        "## Positive Control",
        "",
        f"- Liang Member D (CGCEC, 2033-09-25): **{'PASSED' if pc_passed else 'FAILED'}**",
        "",
        "## Per-Sequence Results",
        "",
        "| Tag | Sequence | Cells | Prefilter OK | Feasible | Shot | Closed | Sweep s | Shoot s |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for s in summaries:
        seq_str = "→".join(s["sequence"])
        lines.append(
            f"| {s['tag']} | `{seq_str}` | {s['n_cells']} | "
            f"{s['n_succeeded_prefilter']} | {s['n_feasible']} | "
            f"{s['n_shot']} | {s['n_close']} | "
            f"{s['elapsed_sweep_s']:.1f} | {s['elapsed_shoot_s']:.1f} |"
        )

    lines += [
        "",
        f"**Totals:** {total_cells} cells, {total_feasible} feasible, "
        f"{total_shot} shot, {total_close} closed.  "
        f"Total elapsed: {total_elapsed_s:.0f}s.",
        "",
    ]

    # Top feasible survivors per sequence
    lines += [
        "## Top Feasible Prefilter Survivors (per sequence, best by defect_ms)",
        "",
    ]
    for s in summaries:
        feasible_rows = s.get("feasible_rows", [])
        tag = s["tag"]
        if not feasible_rows:
            lines.append(f"**{tag}:** no feasible survivors.")
            lines.append("")
            continue
        lines += [
            f"**{tag}:**",
            "",
            "| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |",
            "|---|---|---|---|---|---|---|",
        ]
        for rank, (d_ms, _idx, row, cell) in enumerate(feasible_rows[:5], 1):
            lines.append(
                f"| {rank} | {cell.epoch_iso[:10]} | {list(cell.n_revs)} | "
                f"{list(cell.branches)} | {d_ms:.3f} | "
                f"{row.get('cycle_tof_days', '?'):.1f} | "
                f"{row.get('min_alt_km', '?'):.1f} |"
            )
        lines.append("")

    # N-body shoot results
    lines += [
        "## N-body Shoot Results (all sequences)",
        "",
        "| Tag | Rank | Epoch | Converged | defect_norm | correction_dv_kms | wall_s |",
        "|---|---|---|---|---|---|---|",
    ]
    for s in summaries:
        for sr in s.get("shoot_rows", []):
            sh = sr["shoot"]
            lines.append(
                f"| {sr['tag']} | {sr['rank']} | {sr['epoch_iso'][:10]} | "
                f"{sh.get('converged', False)} | "
                f"{sh.get('defect_norm', 'N/A')} | "
                f"{sh.get('correction_dv_kms', 'N/A')} | "
                f"{sh.get('wall_s', '?'):.1f} |"
            )
    lines.append("")

    if total_close == 0:
        lines += [
            "## Interpretation",
            "",
            f"0 / {total_shot} shot cells closed across all {len(summaries)} sequences.",
            "",
            "The broadened #501 campaign ("
            + ", ".join(s["tag"] for s in summaries)
            + f") over epoch {EPOCH_WINDOW[0]}-{EPOCH_WINDOW[1]}, "
            f"n_revs {N_REVS_RANGE}, is **compute-bounded empty** at the "
            f"{total_cells}-cell budget.",
            "",
            "Each sequence is registered as a compute-bounded empty region in "
            "`data/empty_regions.jsonl`.",
            "Future scale-up, surrogate-based importance sampling, or a qualitatively different "
            "global method is required to break the empty-region wall.",
            "",
        ]
    else:
        lines += [
            "## Flagged Closing Candidates",
            "",
            f"{total_close} candidate(s) closed under jovian_shoot(jacobian=stm).  "
            "**FLAGGED for human adjudication — no novelty claim, no catalogue writeback.**",
            "",
        ]
        for s in summaries:
            for sr in s.get("shoot_rows", []):
                if sr["shoot"].get("converged"):
                    lines += [
                        f"### {s['tag']} rank {sr['rank']}: {sr['epoch_iso'][:10]}",
                        f"- Sequence: `{'→'.join(sr['sequence'])}`",
                        f"- n_revs: {sr['n_revs']}, branches: {sr['branches']}",
                        f"- N-body correction ΔV: "
                        f"{sr['shoot'].get('correction_dv_kms', '?'):.4f} km/s",
                        f"- N-body defect_norm: {sr['shoot'].get('defect_norm', '?'):.4e}",
                        "",
                    ]
        if lit_results:
            lines += [
                "### Literature Check",
                "",
                "| Tag | cell_idx | Epoch | Status | Citation |",
                "|---|---|---|---|---|",
            ]
            for lr in lit_results:
                lines.append(
                    f"| {lr['tag']} | {lr['cell_idx']} | {lr['epoch_iso'][:10]} | "
                    f"{lr['status']} | {lr.get('citation', '—')} |"
                )
            lines.append("")
            lines.append(
                "**`not_found` is necessary-not-sufficient for novelty.  "
                "Human adjudication required before any claim.**"
            )
            lines.append("")

    lines += [
        "## Empty-Region Summary",
        "",
        f"Registered {len([s for s in summaries if s['n_close'] == 0])} compute-bounded "
        "empty regions in `data/empty_regions.jsonl`.",
        "",
        "| Tag | Sequence | Epoch Window | n_revs | cells | verdict |",
        "|---|---|---|---|---|---|",
    ]
    for s in summaries:
        verdict_str = "EMPTY" if s["n_close"] == 0 else f"{s['n_close']} CLOSE"
        lines.append(
            f"| {s['tag']} | `{'→'.join(s['sequence'])}` | "
            f"{EPOCH_WINDOW[0]}-{EPOCH_WINDOW[1]} | "
            f"{N_REVS_RANGE} | {s['n_cells']} | {verdict_str} |"
        )

    lines += [
        "",
        "## Method Notes",
        "",
        "- Positive control (Liang Member D, CGCEC 2033-09-25): PASSED prefilter.",
        "- Prefilter: `evaluate_joint_cell` (Nelder-Mead patched-conic, JUP365 real-eph).",
        "- N-body shoot: `jovian_shoot(jacobian='stm')` -- analytic block-bidiagonal STM Jacobian.",
        "- Convergence: Jones AAS 17-577 s2.5 floors (1e-3 km / 1e-6 km/s per component).",
        f"- Sequences swept: {', '.join(s['tag'] for s in summaries)}.",
        "- NO catalogue writeback. NO novelty claim.",
        "",
    ]

    VERDICT_NOTE.parent.mkdir(parents=True, exist_ok=True)
    VERDICT_NOTE.write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    t_run_start = time.time()
    date_tag = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    import subprocess

    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=REPO_ROOT
        ).strip()
    except Exception:
        git_sha = "unknown"

    _log(f"scan_501_broadened_joint_search start  git={git_sha}")
    _log(f"date -Iseconds: {datetime.now(UTC).isoformat()}")
    _log(
        f"sequences: {[c.tag for c in SEQUENCE_CONFIGS]}  "
        f"n_samples={N_SAMPLES_PER_SEQ}/seq  "
        f"total={len(SEQUENCE_CONFIGS) * N_SAMPLES_PER_SEQ}"
    )

    kernel = _resolve_kernel()
    if kernel is None:
        _log(
            "ERROR: JUP365 kernel not found "
            "(set CYCLERFINDER_JUP365 env var or place at "
            f"{_KERNEL_FALLBACK}). Aborting."
        )
        return 1
    _log(f"JUP365 kernel: {kernel}")

    ephem = jovian.JovianEphemeris(kernel)

    # --- Step 1: Positive control -------------------------------------------
    _log("Step 1: positive control (Liang Member D, CGCEC)")
    pc_cell = liang_member_d_cell()
    pc_verdict = evaluate_joint_cell(pc_cell, ephem)
    defect_ok = pc_verdict.closure_defect_ms < POSITIVE_CONTROL_DEFECT_THRESHOLD_MS
    if not pc_verdict.feasible or not defect_ok:
        _log(
            f"POSITIVE CONTROL FAILED: defect={pc_verdict.closure_defect_ms:.4f} m/s "
            f"feasible={pc_verdict.feasible} min_alt={pc_verdict.min_alt_km:.1f} km"
        )
        _log("ABORTING — do not proceed to broadened sweep with a failing positive control.")
        return 2
    _log(
        f"Positive control PASSED: defect={pc_verdict.closure_defect_ms:.4f} m/s "
        f"min_alt={pc_verdict.min_alt_km:.1f} km cycle_tof={pc_verdict.cycle_tof_days:.2f} d"
    )

    # --- Step 2: Multi-sequence sweep ----------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    shoot_fh = SHOOT_JSONL.open("w")

    summaries: list[dict[str, Any]] = []
    for cfg in SEQUENCE_CONFIGS:
        _log(f"--- Sequence {cfg.tag} {'→'.join(cfg.sequence)} ---")
        summary = _run_sequence(cfg, kernel, ephem, git_sha, date_tag, shoot_fh)
        summaries.append(summary)
        _log(
            f"[{cfg.tag}] done: "
            f"feasible={summary['n_feasible']} shot={summary['n_shot']} "
            f"closed={summary['n_close']}"
        )

    shoot_fh.close()

    # --- Step 3: Literature check for any closers ----------------------------
    all_shoot_rows = [sr for s in summaries for sr in s.get("shoot_rows", [])]
    n_close_total = sum(s["n_close"] for s in summaries)
    lit_results: list[dict[str, Any]] = []
    if n_close_total > 0:
        _log(f"Step 3: literature check for {n_close_total} closing cell(s)")
        lit_results = _check_lit_for_closers(all_shoot_rows)

    # --- Step 4: Final summary -----------------------------------------------
    total_elapsed = time.time() - t_run_start
    _log("=" * 70)
    _log(f"CAMPAIGN SUMMARY  git={git_sha}  date={date_tag}")
    _log(f"  total elapsed: {total_elapsed:.0f}s")
    _log(f"  {'Tag':<8} {'sequence':<40} {'cells':>6} {'feasible':>9} {'shot':>5} {'closed':>7}")
    for s in summaries:
        seq_str = "→".join(s["sequence"])
        _log(
            f"  {s['tag']:<8} {seq_str:<40} {s['n_cells']:>6} {s['n_feasible']:>9} "
            f"{s['n_shot']:>5} {s['n_close']:>7}"
        )
    total_cells = sum(s["n_cells"] for s in summaries)
    total_feasible = sum(s["n_feasible"] for s in summaries)
    total_shot = sum(s["n_shot"] for s in summaries)
    _log(
        f"  {'TOTAL':<8} {'':40} {total_cells:>6} {total_feasible:>9}"
        f" {total_shot:>5} {n_close_total:>7}"
    )
    if n_close_total == 0:
        _log("  verdict: EMPTY (all sequences) — registered in empty_regions.jsonl")
    else:
        _log(f"  verdict: {n_close_total} CLOSE(S) flagged for human adjudication")
        for s in summaries:
            for sr in s.get("shoot_rows", []):
                if sr["shoot"].get("converged"):
                    _log(
                        f"    CLOSE: [{s['tag']}] rank={sr['rank']} "
                        f"epoch={sr['epoch_iso'][:10]} "
                        f"n_revs={sr['n_revs']} branches={sr['branches']}"
                    )

    # --- Step 5: Verdict note ------------------------------------------------
    _write_verdict_note(
        summaries,
        lit_results,
        pc_passed=True,
        date_tag=date_tag,
        git_sha=git_sha,
        total_elapsed_s=total_elapsed,
    )
    _log(f"Verdict note: {VERDICT_NOTE}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
