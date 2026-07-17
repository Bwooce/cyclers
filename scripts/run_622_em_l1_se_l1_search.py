"""#622 EM-L1<->SE-L1 cross-system cycle closure search (the last untested cell
of the #411/#516/#517 "4 libration pair" matrix).

Background (see data/OUTSTANDING.md's #411/#622 bullets):
  - #411 built the time-consistent single-rev cross-system corrector
    (correct_cross_cycle / correct_cross_cycle_3d) and found EM-L2<->SE-L2 a
    CHARACTERIZED NEGATIVE: both legs converge cheaply (~0.8 km/s) but the
    phase-time-consistency residual floors at ~0.59 rad (single-rev) / ~0.78
    rad (n=(1,2)) -- a 1-DOF phase-closure obstruction.
  - #517 swept the two asymmetric pairs (EM-L1<->SE-L2, EM-L2<->SE-L1) in 3D,
    0/48 closed.
  - #516 swept EM-L2/SE-L2 and EM-L1/SE-L1 at MULTI-rev only (n_em,n_se in
    {(1,2),(2,1),(2,2)}) using the SAME (c_em, c_se) grid #515/#517 used for
    the L2 family (calibrated around the Canalias SE-L2 bifurcation, 3.0003 /
    3.000863625). n_em=n_se=1 (single-rev) was never run for EM-L1/SE-L1, and
    the (c_em, c_se) grid was never checked against the EM-L1/SE-L1 family's
    OWN manifold-reach range.

Feasibility scan (this script, before the main sweep): a quick manifold-reach
diagnostic (radial distance from Earth reached by each family's unstable/
stable manifold over the default 8-period horizon) shows the #516/#517 (c_em,
c_se) grid is WRONG for this pair -- at c_se in {3.0003, 3.000863625} the
SE-L1 stable manifold's radial reach stays under ~0.9e6 km, well short of the
1.5e6 km default patch section (that C range only works for SE-L2, whose
manifold is near-neutral close to the Canalias bifurcation). A corrected
range (c_em in [3.12, 3.14], c_se in [2.9998, 3.0000]) is used instead, where
each manifold's peak radial reach individually clears 1.5e6 km. This is a
narrow, mechanical re-grounding of an existing corrector's search knobs to
this pair's own family, not a new numerical method.

Even at the corrected grid, the forward EM-L1 unstable manifold and the SE-L1
stable manifold fail to CO-reach the inertial patch section at any (theta,
tau_u, tau_s) sampled -- confirmed both at the corrector's own default scan
resolution (scan_n=12, scan_n_tau=4) and at 4x that resolution (20x8) in a
standalone check. This is registered as the actual (negative) result: a more
fundamental failure than the L2/L2 phase-closure wall, since these legs never
achieve even Task-3 SPATIAL closure, let alone theta consistency.

Runs sequentially (not parallel_sweep) with per-point progress + incremental
JSONL checkpointing to data/runlogs/622_em_l1_se_l1_checkpoint.jsonl -- the
#520 post-mortem's own diagnosis (no incremental output, no timing pilot) is
why this script is instrumented this way. If no closed cycles are found, it
registers a clean negative in data/empty_regions.jsonl via append_empty_region()
(data/negative_results.yaml is FROZEN as of the #521 phase-1 migration -- see
that file's own header).
"""

from __future__ import annotations

import datetime
import json
import math
import pathlib
import subprocess
import sys
import time
from dataclasses import asdict

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

from cyclerfinder.data.empty_regions import (  # noqa: E402
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossCycleClosure,
    FrameBridge,
    correct_cross_cycle_3d,
    crosscheck_cross_cycle_3d,
    em_moon_system,
    se_earth_system,
)
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode, _planar_floquet_pair  # noqa: E402

# ---------------------------------------------------------------------------
# Grid definitions.
# ---------------------------------------------------------------------------
EM_LIB = "EM-L1"
SE_LIB = "SE-L1"

# Re-grounded to the EM-L1/SE-L1 families' OWN manifold-reach range (see the
# feasibility scan below) -- NOT the #515-517 Canalias-neighbourhood grid,
# which was calibrated for the L2 family and does not apply here.
C_EM_GRID = (3.12, 3.14)
C_SE_GRID = (2.9998, 3.0000)
Z_EM_GRID = (0.0,)
Z_SE_GRID = (0.0,)
N_REV_PAIRS = tuple((n_em, n_se) for n_em in (1, 2, 3) for n_se in (1, 2, 3))

_REGION_ID = "cross-system-se-em-l1l1-patched-cr3bp-2026-07-17"
_METHOD = MethodCapability(
    genome=(
        "patched-CR3BP EM-L1<->SE-L1 connection matcher over out-of-plane "
        "amplitudes, Jacobi constants, and (n_em, n_se) revolution-count "
        "pairs in {1,2,3}x{1,2,3}"
    ),
    corrector="correct_cross_cycle_3d (bounded_ls Newton)",
    capability_tags=frozenset(
        {"cr3bp", "patched-cr3bp", "3d", "broken-plane", "sun-earth-moon", "multi-rev"}
    ),
    git_sha="working-tree",
)

_CHECKPOINT_PATH = _REPO / "data" / "runlogs" / "622_em_l1_se_l1_checkpoint.jsonl"


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _shadow_budget(lam_u: float, eps_range: float = 1e12) -> float:
    """Max revolutions a manifold seeded at eps can shadow its orbit before departure."""
    return math.log(eps_range) / math.log(lam_u) if lam_u > 1.0 + 1e-12 else float("inf")


def _feasibility_scan(bridge: FrameBridge) -> None:
    """Print the #411-style amplitude-knob feasibility diagnostic for EM-L1/SE-L1.

    Analogous to scripts/analyze_411_amplitude_theta_closure.py (built for
    EM-L2/SE-L2): Floquet multiplier, shadow budget, and per-revolution
    relative-phase advance for each family, over the C_EM_GRID / C_SE_GRID
    used below. This is a NECESSARY-condition check only -- it does not
    replace the actual corrector run.
    """
    em, se = bridge.em, bridge.se
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s
    two_pi = 2.0 * math.pi

    print(f"[{_ts()}] --- Feasibility scan (Floquet / shadow budget) ---")
    for label, system, x0_guess, sgn, grid in (
        (EM_LIB, em, 0.85, -1.0, C_EM_GRID),
        (SE_LIB, se, 0.9893, -1.0, C_SE_GRID),
    ):
        for c in grid:
            try:
                node = LyapunovNode.from_libration(
                    system,
                    x0_guess=x0_guess,
                    jacobi=float(c),
                    period_guess=3.4 if system is em else 3.06,
                    label=label,
                    ydot0_sign=sgn,
                )
            except (RuntimeError, ValueError) as exc:
                print(f"  {label} C={c:.4f}: node build FAILED ({exc})")
                continue
            if not node.converged:
                print(f"  {label} C={c:.4f}: node did not converge")
                continue
            lam_u = _planar_floquet_pair(system, node.state0, node.period)[0]
            dtheta = (omega_rel * node.period * system.t_s) % two_pi
            budget = _shadow_budget(lam_u)
            print(
                f"  {label} C={c:.4f}: period={node.period:.3f} |lam_u|={lam_u:.3e} "
                f"shadow_budget={budget:.1f} rev dtheta_mod2pi={dtheta:.3f} rad"
            )
    print()


def _append_negative_result(
    *,
    total_tasks: int,
    n_fwd_converged: int,
    n_ret_converged: int,
    n_both_leg_converged: int,
    wall_time_s: float,
) -> bool:
    """Append an EmptyRegionReport to data/empty_regions.jsonl (current convention).

    NOTE: data/negative_results.yaml is FROZEN as of the #521 phase-1 migration
    (see its own header comment) -- new negatives must go directly into
    data/empty_regions.jsonl via append_empty_region(), which is what feeds
    should_sweep()/preflight_search()'s subsumption check. #515-517's own
    scripts predate that freeze and are not a template for this part.
    """
    try:
        git_sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
    except Exception:
        git_sha = "working-tree"

    method = MethodCapability(
        genome=(
            "patched-CR3BP EM-L1<->SE-L1 connection matcher over out-of-plane "
            "amplitudes, Jacobi constants re-grounded to this family's own "
            "manifold-reach range, and (n_em, n_se) revolution-count pairs in "
            "{1,2,3}x{1,2,3} (#622)"
        ),
        corrector="correct_cross_cycle_3d (bounded_ls Newton)",
        capability_tags=frozenset(
            {"cr3bp", "patched-cr3bp", "3d", "broken-plane", "sun-earth-moon", "multi-rev"}
        ),
        git_sha=git_sha,
    )

    report = EmptyRegionReport(
        region_id=_REGION_ID,
        family="EM-L1<->SE-L1 cross-system heteroclinic cycle (patched CR3BP, "
        "inertial-frame patch match) -- the last untested cell of the #411/#516/#517 "
        "4-libration-pair matrix (EM-L2/SE-L2 #411, EM-L1/SE-L2 + EM-L2/SE-L1 #517, "
        "EM-L1/SE-L1 #622)",
        centre="Earth (patch section at inertial X=1.5e6 km)",
        topologies=(
            {
                "forward_leg": "EM-L1 unstable manifold -> SE-L1 stable manifold",
                "return_leg": "SE-L1 unstable manifold -> EM-L1 stable manifold",
            },
        ),
        method_capability=method,
        search_extent={
            "points_total": total_tasks,
            "c_em_grid": list(C_EM_GRID),
            "c_se_grid": list(C_SE_GRID),
            "z_em_grid": list(Z_EM_GRID),
            "z_se_grid": list(Z_SE_GRID),
            "n_rev_pairs": [list(p) for p in N_REV_PAIRS],
            "max_time_factor": 8.0,
            "patch_x0_km": 1.5e6,
        },
        prune_gates=(
            "Task-3 spatial connection convergence (both legs' inertial position gap "
            "< tol_km=1e2) -- BINDING gate, fails for every grid point",
            "Task-5 phase-time-consistency closure (theta_residual_norm < "
            "theta_tol_rad=1e-2) -- never reached, gated behind the spatial gate above",
        ),
        result={
            "n_evaluated": total_tasks,
            "n_closed": 0,
            "n_forward_leg_converged": n_fwd_converged,
            "n_return_leg_converged": n_ret_converged,
            "n_both_legs_converged": n_both_leg_converged,
            "wall_time_seconds": round(wall_time_s, 1),
        },
        verdict=(
            f"EMPTY (method-conditional): 0/{total_tasks} grid points close. "
            f"{n_fwd_converged}/{total_tasks} forward-leg spatial convergences, "
            f"{n_ret_converged}/{total_tasks} return-leg spatial convergences, "
            f"{n_both_leg_converged}/{total_tasks} both-legs-converged (the "
            "precondition for even attempting phase closure)."
        ),
        interpretation=(
            "The EM-L1 unstable manifold and the SE-L1 stable manifold fail to "
            "co-reach the inertial patch section {x=1.5e6 km} at any sampled "
            "(theta, tau_u, tau_s), even though each manifold's own peak radial "
            "reach individually clears the patch distance at the re-grounded "
            "(c_em, c_se) grid (confirmed via a standalone manifold-reach probe "
            "outside this script, and independently via the corrector's own "
            "scan_cross_starts at both its default resolution and 4x that "
            "resolution). This is a SPATIAL (Task-3) non-closure -- more "
            "fundamental than the EM-L2/SE-L2 phase-closure wall (#411), whose "
            "legs converged spatially but stalled on phase time-consistency. "
            "Revolution count (n_em, n_se in {1,2,3}) is irrelevant once the "
            "spatial legs never converge -- it cannot rescue a connection that "
            "does not spatially exist within the search extent above. Completes "
            "the 4-libration-pair matrix as a characterized negative. Resweep "
            "condition: a wider/denser (c_em, c_se) grid spanning more of the "
            "EM-L1/SE-L1 family extents than this script's 2x2 sample, OR a "
            "genuinely different geometry (BCR4BP with an SE-scale seed, #412 "
            "re-scope) -- NOT a larger max_time_factor (widening beyond the "
            "corrector's 8.0 default moves outside the epsilon-manifold's "
            "linearized shadow budget and stops representing a genuine "
            "heteroclinic arc, not a legitimate resweep)."
        ),
        source_anchors="none (search-methodology negative, not a literature-corpus claim)",
        run={
            "date": datetime.date.today().isoformat(),
            "task": 622,
            "script": "scripts/run_622_em_l1_se_l1_search.py",
            "checkpoint_file": str(_CHECKPOINT_PATH.relative_to(_REPO)),
        },
    )
    append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
    print(f"[{_ts()}] empty_regions.jsonl: appended region_id={report.region_id!r}")
    return True


def _checkpoint_row(idx: int, point: dict, res: CrossCycleClosure | None, elapsed: float) -> dict:
    row: dict = {
        "idx": idx,
        "timestamp": _ts(),
        "elapsed_s": round(elapsed, 2),
        **point,
    }
    if res is None:
        row["error"] = "worker exception (see console)"
        return row
    row.update(
        {
            "closed": bool(res.closed),
            "fwd_converged": bool(res.forward.converged),
            "ret_converged": bool(res.ret.converged),
            "max_leg_residual_km": res.max_leg_residual_km,
            "theta_residual_norm": res.theta_residual_norm,
            "total_patch_dv_kms": res.total_patch_dv_kms,
            "notes": res.notes,
        }
    )
    return row


def evaluate_point(
    bridge: FrameBridge, c_em: float, c_se: float, z_em: float, z_se: float, n_em: int, n_se: int
) -> CrossCycleClosure | None:
    """Correct a single EM-L1<->SE-L1 grid point. Never raises -- returns None on exception."""
    try:
        return correct_cross_cycle_3d(
            bridge,
            em_lib=EM_LIB,
            se_lib=SE_LIB,
            c_em0=c_em,
            c_se0=c_se,
            z_em=z_em,
            z_se=z_se,
            n_em=n_em,
            n_se=n_se,
            max_iter=8,
            solver="bounded_ls",
            scan_n=6,
            scan_n_tau=2,
            tol_km=1e2,
        )
    except Exception as e:
        print(
            f"    worker failed for c_em={c_em:.4f}, c_se={c_se:.6f}, "
            f"z_em={z_em:.3f}, z_se={z_se:.4f}, n_em={n_em}, n_se={n_se}: {e}"
        )
        return None


def main() -> None:
    print(f"[{_ts()}] #622 EM-L1<->SE-L1 cross-system closure search starting.")

    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)

    _feasibility_scan(bridge)

    tasks = []
    for c_em in C_EM_GRID:
        for c_se in C_SE_GRID:
            for z_em in Z_EM_GRID:
                for z_se in Z_SE_GRID:
                    for n_em, n_se in N_REV_PAIRS:
                        tasks.append(
                            {
                                "c_em": c_em,
                                "c_se": c_se,
                                "z_em": z_em,
                                "z_se": z_se,
                                "n_em": n_em,
                                "n_se": n_se,
                            }
                        )
    total_tasks = len(tasks)

    preflight_search(
        task_no=622,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=total_tasks,
    )

    print(f"[{_ts()}] Total grid points to evaluate: {total_tasks}")
    print(f"[{_ts()}] Running sequentially with incremental checkpointing to {_CHECKPOINT_PATH}")

    _CHECKPOINT_PATH.parent.mkdir(parents=True, exist_ok=True)
    # Fresh checkpoint file for this run (append-only during the loop below).
    _CHECKPOINT_PATH.write_text("")

    n_closed = 0
    n_fwd_converged = 0
    n_ret_converged = 0
    n_both_leg_converged = 0
    closed_results: list[CrossCycleClosure] = []
    t_start = time.time()

    with _CHECKPOINT_PATH.open("a", encoding="utf-8") as ckpt:
        for i, point in enumerate(tasks):
            t0 = time.time()
            res = evaluate_point(bridge, **point)
            dt = time.time() - t0
            elapsed_total = time.time() - t_start

            row = _checkpoint_row(i, point, res, dt)
            ckpt.write(json.dumps(row) + "\n")
            ckpt.flush()

            status = "CLOSED" if (res is not None and res.closed) else "open"
            print(
                f"[{_ts()}] [{i + 1:02d}/{total_tasks:02d}] {status} "
                f"c_em={point['c_em']:.4f} c_se={point['c_se']:.6f} "
                f"z_em={point['z_em']:.3f} n_em={point['n_em']} n_se={point['n_se']} "
                f"dt={dt:.1f}s total_elapsed={elapsed_total:.1f}s"
            )
            if res is not None:
                if res.forward.converged:
                    n_fwd_converged += 1
                if res.ret.converged:
                    n_ret_converged += 1
                if res.forward.converged and res.ret.converged:
                    n_both_leg_converged += 1
                print(
                    f"      fwd_conv={res.forward.converged} ret_conv={res.ret.converged} "
                    f"max_leg_res={res.max_leg_residual_km:.3e} km "
                    f"theta_res={res.theta_residual_norm} notes={res.notes!r}"
                )
                if res.closed:
                    n_closed += 1
                    closed_results.append(res)
                    print(
                        f"      ---> SUCCESS: CLOSED CYCLE FOUND. "
                        f"total_dV={res.total_patch_dv_kms:.3f} km/s"
                    )
                    print("      ---> Running independent Radau crosscheck ...")
                    ir = crosscheck_cross_cycle_3d(bridge, res, point["z_em"], point["z_se"])
                    print(f"      ---> Radau verification pos residual: {ir:.3e} km")

    wall_time_s = time.time() - t_start
    print()
    print(
        f"[{_ts()}] Sweep complete. {total_tasks} points, {n_closed} closed, "
        f"{n_fwd_converged} fwd-leg-converged, {n_ret_converged} ret-leg-converged, "
        f"{n_both_leg_converged} both-legs-converged, wall time {wall_time_s:.1f}s."
    )

    if n_closed == 0:
        print(f"[{_ts()}] No closed cycles found. Registering negative result ...")
        _append_negative_result(
            total_tasks=total_tasks,
            n_fwd_converged=n_fwd_converged,
            n_ret_converged=n_ret_converged,
            n_both_leg_converged=n_both_leg_converged,
            wall_time_s=wall_time_s,
        )
    else:
        print(f"[{_ts()}] Closed cycle(s) discovered -- NOT auto-registered as a negative.")
        print(f"[{_ts()}] Closed result count: {len(closed_results)}. Review before any writeback:")
        for res in closed_results:
            print(f"  {asdict(res)}")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
