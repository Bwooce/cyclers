"""#318 Phase 2b — Sobol joint-search smoke driver for Galilean CGCEC short cyclers.

256 Sobol cells (CGCEC, epoch window 2033-01-01 → 2035-01-01) evaluated via the
real-ephemeris patched-conic surrogate (:func:`evaluate_joint_cell`), then top-10
feasible survivors n-body shot via :func:`jovian_shoot` with the analytic STM
Jacobian.

**Positive control (Liang Member D)** is verified BEFORE the sweep. If it fails
the prefilter, the driver aborts — do not trust a gauntlet that can't find its
positive control ([[feedback_verify_gauntlet_with_positive_control]]).

Outputs::

    data/scan_318_sobol_smoke.jsonl       -- per-cell prefilter results
    data/empty_regions.jsonl              -- appended if 0 cells close under shoot
    docs/notes/2026-06-30-318-phase2b-smoke-verdict.md  -- verdict note

Constraints:
  - NO catalogue writeback
  - NO novelty claims
  - SMOKE scale only: 256 cells
"""

from __future__ import annotations

import json
import sys
import time
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

OUT_JSONL = REPO_ROOT / "data" / "scan_318_sobol_smoke.jsonl"
EMPTY_REGIONS_PATH = REPO_ROOT / "data" / "empty_regions.jsonl"
VERDICT_NOTE = REPO_ROOT / "docs" / "notes" / "2026-06-30-318-phase2b-smoke-verdict.md"

N_SAMPLES = 256
SEQUENCE = list(jovian.CGCEC)
EPOCH_WINDOW = ("2033-01-01", "2035-01-01")
N_REVS_RANGE = (1, 2)
TOF_SEED_RANGE = (15.0, 45.0)  # days per leg (CGCEC legs ~18-32 d; give extra slack)
POWERED_MIN_ALT_KM = 100.0  # Liang's flyby floor
TOP_K_SHOOT = 10  # n-body shoot the best K survivors

# Positive control: Liang Member D epoch (2033-09-25) is in [2033-01-01, 2035-01-01].
POSITIVE_CONTROL_DEFECT_THRESHOLD_MS = 1.0
POSITIVE_CONTROL_ALT_THRESHOLD_KM = 100.0

# A cell "closes" under n-body shoot if the corrector converges.
# No minimum correction_dv_kms threshold — converged IS the criterion.

_KERNEL_FALLBACK = str(Path.home() / "dev" / "references" / "kernels" / "jup365.bsp")


def _resolve_kernel() -> str | None:
    """Env var first, then well-known fallback — mirrors test_joint_sobol._resolve_kernel."""
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
# Top-level picklable process_cell — must be at module level for loky
# ---------------------------------------------------------------------------


def process_cell(cell: JointCell) -> dict[str, Any] | None:
    """Evaluate one Sobol joint-cell via the patched-conic surrogate.

    Creates a fresh :class:`JovianEphemeris` from the kernel path each call;
    loky workers inherit the env and furnish once per process (module-level cache).
    Returns ``None`` when no kernel is available or the cell raises.
    """
    kernel = _resolve_kernel()
    if kernel is None:
        return None
    try:
        ephem = jovian.JovianEphemeris(kernel)
        v = evaluate_joint_cell(cell, ephem)
        return {
            "epoch_iso": cell.epoch_iso,
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
# N-body shoot helpers (main-process only — not picklable, not called in workers)
# ---------------------------------------------------------------------------


def _build_jovian_seed(
    cell: JointCell,
    ephem: jovian.JovianEphemeris,
) -> Any | None:
    """Build a ShootingSeed from a JointCell's converged patched-conic solution.

    Re-runs :func:`optimize_cycle` to get the converged epochs, then re-solves
    the legs to obtain the per-leg V∞ vectors needed for the node states.
    Returns ``None`` when the conic doesn't converge.
    """
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

        # Outgoing V∞ at node k: sol_vinf_out[k] for k < n_legs; wrap uses [0].
        vi_out = np.asarray(sol_vinf_out[k] if k < n_legs else sol_vinf_out[0], dtype=np.float64)
        # Inbound V∞ at node k: node 0 uses wrap arrival = sol_vinf_in[n_legs-1].
        vi_in = np.asarray(
            sol_vinf_in[k - 1] if k > 0 else sol_vinf_in[n_legs - 1], dtype=np.float64
        )

        # Wrap node (k == n_legs) state uses inbound V∞ (arrived); others use outgoing.
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
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    date_tag = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    import subprocess

    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=REPO_ROOT
        ).strip()
    except Exception:
        git_sha = "unknown"

    _log(f"scan_318_joint_sobol_smoke start  git={git_sha}")
    _log(f"date -Iseconds: {datetime.now(UTC).isoformat()}")

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
    _log("Step 1: positive control (Liang Member D)")
    pc_cell = liang_member_d_cell()
    pc_verdict = evaluate_joint_cell(pc_cell, ephem)
    defect_ok = pc_verdict.closure_defect_ms < POSITIVE_CONTROL_DEFECT_THRESHOLD_MS
    if not pc_verdict.feasible or not defect_ok:
        _log(
            f"POSITIVE CONTROL FAILED: defect={pc_verdict.closure_defect_ms:.4f} m/s "
            f"feasible={pc_verdict.feasible} min_alt={pc_verdict.min_alt_km:.1f} km"
        )
        _log("ABORTING — do not proceed to smoke run with a failing positive control.")
        return 2
    _log(
        f"Positive control PASSED: defect={pc_verdict.closure_defect_ms:.4f} m/s "
        f"min_alt={pc_verdict.min_alt_km:.1f} km cycle_tof={pc_verdict.cycle_tof_days:.2f} d"
    )

    # --- Step 2: Build Sobol cells ------------------------------------------
    _log(f"Step 2: building {N_SAMPLES} Sobol cells  sequence={SEQUENCE}")
    _log(f"  epoch_window={EPOCH_WINDOW}  n_revs={N_REVS_RANGE}  tof={TOF_SEED_RANGE} d")
    cells = make_sobol_cells(
        n_samples=N_SAMPLES,
        sequence=SEQUENCE,
        epoch_window=EPOCH_WINDOW,
        n_revs_range=N_REVS_RANGE,
        tof_seed_range=TOF_SEED_RANGE,
        powered_min_alt_km=POWERED_MIN_ALT_KM,
        seed=0,
    )
    _log(f"  {len(cells)} cells built")

    # --- Step 3: Parallel patched-conic prefilter sweep ---------------------
    _log("Step 3: parallel patched-conic prefilter sweep (loky)")
    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    stream_fh = OUT_JSONL.open("w")

    t_sweep_start = time.time()
    cfg = ParallelSweepConfig(
        n_workers=-1,
        chunk_size=4,
        backend="loky",
        verbose=0,
        raise_on_first_error=False,
    )
    sweep = parallel_sweep(cells, process_cell, config=cfg)
    t_sweep_end = time.time()

    n_feasible = 0
    feasible_rows: list[tuple[float, int, dict[str, Any], JointCell]] = []
    for i, row in enumerate(sweep.results):
        if row is None:
            continue
        stream_fh.write(json.dumps(row) + "\n")
        if row.get("feasible"):
            n_feasible += 1
            defect = float(row.get("closure_defect_ms", float("inf")))
            feasible_rows.append((defect, i, row, cells[i]))
    stream_fh.close()
    elapsed_sweep = t_sweep_end - t_sweep_start

    feasible_rows.sort(key=lambda x: x[0])
    _log(
        f"Sweep done in {elapsed_sweep:.1f}s: "
        f"total={N_SAMPLES} succeeded={sweep.n_succeeded} failed={sweep.n_failed} "
        f"feasible={n_feasible}"
    )
    if sweep.notes:
        _log(f"  sweep notes: {sweep.notes}")

    # --- Step 4: N-body shoot top-K feasible survivors ----------------------
    shoot_rows: list[dict[str, Any]] = []
    n_close = 0
    t_shoot_start = time.time()
    top_k = feasible_rows[:TOP_K_SHOOT]
    if not top_k:
        _log("Step 4: no feasible survivors to shoot — skipping n-body stage")
    else:
        _log(
            f"Step 4: n-body shoot (jovian_shoot, jacobian=stm) on top {len(top_k)} feasible cells"
        )
        for rank, (defect, idx, prefilter_row, cell) in enumerate(top_k):
            _log(
                f"  shooting rank {rank + 1}/{len(top_k)}: "
                f"epoch={cell.epoch_iso}  defect={defect:.4f} m/s"
            )
            shoot_result = _shoot_cell(cell, kernel, ephem)
            row = {
                "rank": rank + 1,
                "cell_idx": idx,
                "prefilter": prefilter_row,
                "shoot": shoot_result,
            }
            shoot_rows.append(row)
            if shoot_result.get("converged"):
                n_close += 1
                _log(
                    f"    CLOSED: defect_norm={shoot_result.get('defect_norm', '?'):.4e} "
                    f"correction_dv={shoot_result.get('correction_dv_kms', '?'):.4f} km/s"
                )
            else:
                _log(
                    f"    not closed: defect_norm={shoot_result.get('defect_norm', '?')}"
                    + (f"  error={shoot_result['error']!r}" if "error" in shoot_result else "")
                )
    t_shoot_end = time.time()
    elapsed_shoot = t_shoot_end - t_shoot_start

    # --- Step 5: Literature check for any closers ---------------------------
    lit_results: list[dict[str, Any]] = []
    if n_close > 0:
        _log(f"Step 5: literature check for {n_close} closing cells")
        try:
            from cyclerfinder.search.literature_check import (
                CandidateSignature,
                check_literature,
            )

            try:
                from claude_code_tools import web_search as _ws  # type: ignore[import-not-found]

                def _search_fn(q: str) -> list[dict[str, Any]]:
                    return list(_ws(q))

            except ImportError:
                _search_fn = None  # type: ignore[assignment]

            for sr in shoot_rows:
                if not sr["shoot"].get("converged"):
                    continue
                cell = cells[sr["cell_idx"]]
                sig = CandidateSignature(
                    primary="Jupiter",
                    sequence=cell.sequence,
                    n_rev=cell.n_revs,
                    resonances=(),
                )
                if _search_fn is not None:
                    lit = check_literature(sig, search=_search_fn)
                    lit_results.append(
                        {"cell_idx": sr["cell_idx"], "status": lit.status, "citation": lit.citation}
                    )
                    _log(f"  lit_check status={lit.status} citation={lit.citation!r}")
                else:
                    _log("  WebSearch not available — skipping literature check")
        except Exception as exc:
            _log(f"  literature check skipped: {exc}")

    # --- Step 6: Empty-regions entry if nothing closes ----------------------
    summary = {
        "n_cells": N_SAMPLES,
        "n_succeeded_prefilter": int(sweep.n_succeeded),
        "n_failed_prefilter": int(sweep.n_failed),
        "n_feasible": n_feasible,
        "n_shot": len(top_k),
        "n_close": n_close,
        "elapsed_sweep_s": round(elapsed_sweep, 1),
        "elapsed_shoot_s": round(elapsed_shoot, 1),
        "git_sha": git_sha,
        "date": date_tag,
    }

    if n_close == 0:
        _log("Step 6: registering EMPTY region (0 closes under n-body shoot)")
        empty_entry = {
            "region_id": f"jovian-cgcec-sobol-smoke-318-{datetime.now(UTC).strftime('%Y-%m-%d')}",
            "family": "planet-centric moon system (Jupiter)",
            "centre": "Jupiter",
            "topologies": [
                {
                    "sequence": list(SEQUENCE),
                    "per_leg_revs": list(N_REVS_RANGE),
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
                "n_sobol_samples": N_SAMPLES,
                "sequence": list(SEQUENCE),
                "epoch_window": list(EPOCH_WINDOW),
                "n_revs_range": list(N_REVS_RANGE),
                "tof_seed_range_days": list(TOF_SEED_RANGE),
                "top_k_shot": TOP_K_SHOOT,
                "ephem_model": "JUP365 real-ephemeris",
                "center": "Jupiter",
            },
            "prune_gates": [
                "evaluate_joint_cell prefilter (feasible + altitude + closure_defect)",
                "jovian_shoot convergence (continuity floor Jones AAS 17-577 §2.5)",
            ],
            "result": summary,
            "verdict": (
                "EMPTY -- 0 cells closed under jovian_shoot(jacobian=stm) "
                "within the Sobol smoke budget"
            ),
            "interpretation": (
                "256-cell Sobol smoke scan in the CGCEC joint manifold "
                "(epoch 2033-2035, 1-2 revs/leg, 15-45 d/leg) produced "
                f"{n_feasible} feasible prefilter survivors; top-{TOP_K_SHOOT} "
                "n-body shot, 0 converged. This is a compute-bounded negative-region map "
                "(design doc: the realistic Phase-2 outcome is a compute-bounded empty-region map, "
                "not a discovery). Scale-up or surrogate-based sampling "
                "required to improve coverage."
            ),
            "run": {
                "date": date_tag,
                "git_sha": git_sha,
                "elapsed_sweep_s": elapsed_sweep,
                "elapsed_shoot_s": elapsed_shoot,
            },
        }
        with EMPTY_REGIONS_PATH.open("a") as fh:
            fh.write(json.dumps(empty_entry) + "\n")
        _log(f"  appended to {EMPTY_REGIONS_PATH}")

    # --- Step 7: Print summary -----------------------------------------------
    _log("=" * 60)
    _log(f"SUMMARY  git={git_sha}  date={date_tag}")
    _log(f"  cells:          {N_SAMPLES}")
    _log(f"  prefilter ok:   {sweep.n_succeeded}  (failed: {sweep.n_failed})")
    _log(f"  feasible:       {n_feasible}")
    _log(f"  shot (top-{TOP_K_SHOOT}): {len(top_k)}")
    _log(f"  closed (n-body): {n_close}")
    _log(f"  sweep wall:     {elapsed_sweep:.1f}s")
    _log(f"  shoot wall:     {elapsed_shoot:.1f}s")
    _log(f"  output:         {OUT_JSONL}")
    if n_close == 0:
        _log("  verdict:        EMPTY — registered in empty_regions.jsonl")
    else:
        _log(f"  verdict:        {n_close} CLOSING candidate(s) — see shoot_rows below")
        for sr in shoot_rows:
            if sr["shoot"].get("converged"):
                c = cells[sr["cell_idx"]]
                _log(
                    f"    CLOSE #{sr['rank']}: epoch={c.epoch_iso}  "
                    f"n_revs={c.n_revs}  branches={c.branches}"
                )

    # --- Step 8: Write verdict note -----------------------------------------
    _write_verdict_note(summary, feasible_rows, shoot_rows, top_k, n_close, lit_results)
    _log(f"Verdict note: {VERDICT_NOTE}")

    return 0


def _write_verdict_note(
    summary: dict[str, Any],
    feasible_rows: list[tuple[float, int, dict[str, Any], JointCell]],
    shoot_rows: list[dict[str, Any]],
    top_k: list[tuple[float, int, dict[str, Any], JointCell]],
    n_close: int,
    lit_results: list[dict[str, Any]],
) -> None:
    date_tag = summary["date"]
    git_sha = summary["git_sha"]
    verdict = "EMPTY" if n_close == 0 else f"{n_close} CLOSE(S)"

    lines = [
        f"# #318 Phase 2b — Sobol Smoke Verdict ({date_tag[:10]})",
        "",
        f"**Date:** {date_tag[:10]}. **Status:** Phase-2b smoke run complete. "
        f"**Verdict:** {verdict}.",
        "",
        "## Configuration",
        "",
        f"- Sequence: `{'→'.join(SEQUENCE)}`",
        f"- Epoch window: {EPOCH_WINDOW[0]} to {EPOCH_WINDOW[1]}",
        f"- n_revs range: {N_REVS_RANGE}",
        f"- ToF seed range: {TOF_SEED_RANGE} d/leg",
        f"- Sobol samples: {N_SAMPLES} (seed=0, scrambled)",
        f"- N-body shot top-K: {TOP_K_SHOOT}",
        f"- git: `{git_sha}`",
        "",
        "## Results",
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Cells submitted | {summary['n_cells']} |",
        f"| Prefilter succeeded | {summary['n_succeeded_prefilter']} |",
        f"| Prefilter failed | {summary['n_failed_prefilter']} |",
        f"| Feasible (surrogate) | {summary['n_feasible']} |",
        f"| N-body shot | {summary['n_shot']} |",
        f"| Closed (jovian_shoot) | {summary['n_close']} |",
        f"| Sweep wall time | {summary['elapsed_sweep_s']:.1f}s |",
        f"| Shoot wall time | {summary['elapsed_shoot_s']:.1f}s |",
        "",
    ]

    if feasible_rows:
        lines += [
            "## Top-10 Feasible Prefilter Survivors (by closure_defect_ms)",
            "",
            "| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |",
            "|---|---|---|---|---|---|---|",
        ]
        for rank, (d_ms, _idx, row, cell) in enumerate(feasible_rows[:10], 1):
            lines.append(
                f"| {rank} | {cell.epoch_iso[:10]} | {list(cell.n_revs)} | "
                f"{list(cell.branches)} | {d_ms:.3f} | "
                f"{row.get('cycle_tof_days', '?'):.1f} | "
                f"{row.get('min_alt_km', '?'):.1f} |"
            )
        lines.append("")

    if shoot_rows:
        lines += [
            f"## N-body Shoot Results (top-{TOP_K_SHOOT})",
            "",
            "| Rank | Epoch | Converged | defect_norm | correction_dv_kms | wall_s |",
            "|---|---|---|---|---|---|",
        ]
        for sr in shoot_rows:
            cell = None
            for _defect, idx, _row, c in top_k:
                if idx == sr["cell_idx"]:
                    cell = c
                    break
            epoch = cell.epoch_iso[:10] if cell else "?"
            sh = sr["shoot"]
            lines.append(
                f"| {sr['rank']} | {epoch} | {sh.get('converged', False)} | "
                f"{sh.get('defect_norm', 'N/A')} | "
                f"{sh.get('correction_dv_kms', 'N/A')} | "
                f"{sh.get('wall_s', '?'):.1f} |"
            )
        lines.append("")

    if n_close == 0:
        lines += [
            "## Interpretation",
            "",
            f"0 / {TOP_K_SHOOT} shot cells closed under `jovian_shoot(jacobian='stm')`.",
            f"The {N_SAMPLES}-cell Sobol smoke scan in the CGCEC joint manifold "
            f"(epoch {EPOCH_WINDOW[0]}-{EPOCH_WINDOW[1]}, "
            f"n_revs {N_REVS_RANGE}, ToF {TOF_SEED_RANGE} d/leg) is **compute-bounded empty**.",
            "",
            "This is consistent with the Phase-2 design doc's honest prior:",
            "> *The realistic Phase-2 outcome is a compute-bounded "
            "empty-region map, not a discovery*",
            "",
            "Registered in `data/empty_regions.jsonl`. Scale-up or surrogate-based sampling",
            "required to improve coverage beyond this smoke result.",
            "",
        ]
    else:
        lines += [
            "## Closing Candidates",
            "",
            f"{n_close} cell(s) closed under `jovian_shoot(jacobian='stm')`. "
            "See shoot results above.",
            "",
        ]
        if lit_results:
            lines += [
                "### Literature Check",
                "",
                "| cell_idx | status | citation |",
                "|---|---|---|",
            ]
            for lr in lit_results:
                lines.append(f"| {lr['cell_idx']} | {lr['status']} | {lr.get('citation', '—')} |")
            lines.append("")

    lines += [
        "## Method Notes",
        "",
        "- Positive control (Liang Member D, 2033-09-25): PASSED prefilter.",
        "- Prefilter: `evaluate_joint_cell` (Nelder-Mead patched-conic, JUP365 real-eph).",
        "- N-body shoot: `jovian_shoot(jacobian='stm')` — analytic block-bidiagonal STM Jacobian,",
        "  DOP853+STM co-integration per leg, avoids REBOUND variational-particle gotcha.",
        "- Convergence criterion: Jones AAS 17-577 §2.5 floors "
        "(1e-3 km / 1e-6 km/s per component).",
        "- NO catalogue writeback. NO novelty claim.",
        "",
    ]

    VERDICT_NOTE.parent.mkdir(parents=True, exist_ok=True)
    VERDICT_NOTE.write_text("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
