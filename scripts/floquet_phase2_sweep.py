"""#347 Phase 2 P2.3 — Saddle-center sweep driver.

For each parent in
:data:`cyclerfinder.genome.floquet_phase2_parents.PHASE2_SWEEP_PARENTS`:

  1. Recover the parent's symmetric IC via the existing pipeline
     :func:`cyclerfinder.search.reachable_representatives.recover_all_cyclers_braik_ross`
     (the Braik-Ross 2026 Table 2 anchor).
  2. Run :func:`cyclerfinder.search.cr3bp_continuation.continue_family` in the
     parent's cj_window using the parent's dc + n_steps.
  3. Build a list of :class:`cyclerfinder.search.bifurcation_detector.FamilyMember`
     from the converged members; run
     :func:`cyclerfinder.search.bifurcation_detector.detect_saddle_center_bracket`
     with ``stm_mode="fixed_path"`` (#372 P372.3 finding: variable-step
     trivial-pair smear at 0.35% is comparable to k=1 saddle-center signal).
  4. For each saddle-center bracket found, hand the POST-bifurcation member to
     :func:`cyclerfinder.genome.asymmetric_branch.branch_at_saddle_center`
     (with the Phase 2 P2.1 Gram-Schmidt fix in place) at epsilon=5e-4.
  5. Emit one JSONL row per (parent, bracket, branched orbit) tuple to
     ``data/floquet_phase2_sweep_results.jsonl``.

Phase 2 is DISCOVERY only: each emitted row is a "discovery candidate", NOT a
catalogue admission. Catalogue writeback is Phase 4+ (gauntlet) work.

Per the discipline:
  * Sourced anchors only (Braik-Ross Table 2 + Ross-RT 2025 C_J_C21).
  * Independent topology cross-check via ``winding_topology``.
  * Independent closure cross-check is in the corrector's gate
    (Periodic3DOrbit.independent_closure_residual).
  * HONEST_NEGATIVES are first-class outcomes: a parent that does not yield
    a saddle-center bracket (e.g. C21's tiny family-extent) or a branch
    corrector that does not converge MUST emit a row recording the
    no-find / no-converge, not silently skip.

Usage:
    uv run python scripts/floquet_phase2_sweep.py [--output PATH]
                                                  [--epsilon EPS]
                                                  [--parents LABEL[,LABEL...]]
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from pathlib import Path

import numpy as np

from cyclerfinder.genome.asymmetric_branch import branch_at_saddle_center
from cyclerfinder.genome.floquet_phase2_parents import (
    PHASE2_SWEEP_PARENTS,
    SweepParent,
)
from cyclerfinder.search.bifurcation_detector import (
    FamilyMember,
    detect_saddle_center_bracket,
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_continuation import continue_family
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit
from cyclerfinder.search.reachable_representatives import (
    _CYCLER_SEEDS,
    TU_DAYS,
    braik_ross_system,
    recover_all_cyclers_braik_ross,
)


def _recover_parent_seed(parent: SweepParent) -> SymmetricOrbit:
    """Recover the parent's anchor symmetric orbit via the existing pipeline.

    For `_down`-suffixed labels, recover the base family (e.g. C32 for C32_down)
    — the only difference is the continuation direction the driver applies.
    """
    base_label = parent.label.removesuffix("_down")
    if base_label not in _CYCLER_SEEDS:
        raise ValueError(f"unknown parent label {parent.label} (base {base_label})")
    system = braik_ross_system()
    reps = recover_all_cyclers_braik_ross(system)
    rep = next(r for r in reps if r.label == base_label)
    # Anchor-stage gate: confirmed period within 1% of sourced.
    if not rep.confirmed:
        raise RuntimeError(
            f"parent {parent.label} anchor recovery NOT confirmed "
            f"(period={rep.period_days:.3f} d, sourced={parent.sourced_period_days:.3f} d)"
        )
    seed = SymmetricOrbit(
        x0=float(rep.state0[0]),
        ydot0=float(rep.state0[4]),
        jacobi=rep.jacobi,
        t_half=rep.period * 0.5,
        period=rep.period,
        converged=True,
        crossing_residual=1e-12,
        n_iter=0,
    )
    return seed


def _topology_check(system_mu: float, state0: np.ndarray, period: float) -> tuple[int, int, bool]:
    """Independent topology cross-check via ``winding_topology``.

    Returns (k1, k2, prograde). Wraps any RuntimeError / ValueError into a
    sentinel (-1, -1, False) so the sweep can record "topology check failed"
    without aborting the whole run.
    """
    try:
        topo = winding_topology(system_mu, state0, period)
        return int(topo.k1), int(topo.k2), bool(topo.prograde)
    except (RuntimeError, ValueError):
        return -1, -1, False


def _process_parent(
    parent: SweepParent,
    *,
    epsilon: float,
    output_fh,
    log_prefix: str,
) -> dict[str, int]:
    """Run the saddle-center sweep for one parent. Emits JSONL rows.

    Returns a dict of counters: members_walked, brackets_found, branches_converged,
    branches_topology_changed, errors.
    """
    counters = {
        "members_walked": 0,
        "brackets_found": 0,
        "branches_converged": 0,
        "branches_topology_changed": 0,
        "errors": 0,
    }
    t_parent_start = time.time()
    system = braik_ross_system()

    # Step 1: recover the anchor.
    print(f"{log_prefix} parent={parent.label}: recovering anchor...")
    try:
        seed = _recover_parent_seed(parent)
    except (RuntimeError, ValueError) as e:
        print(f"{log_prefix} parent={parent.label}: anchor recovery FAILED — {e}")
        output_fh.write(
            json.dumps(
                {
                    "kind": "parent_error",
                    "parent_label": parent.label,
                    "stage": "anchor_recovery",
                    "error": str(e),
                    "elapsed_seconds": time.time() - t_parent_start,
                }
            )
            + "\n"
        )
        counters["errors"] += 1
        return counters
    base_label = parent.label.removesuffix("_down")
    seed_meta = _CYCLER_SEEDS[base_label]

    # Step 2: natural-parameter continuation.
    direction = -1 if parent.label.endswith("_down") else +1
    c_min, c_max = parent.cj_window
    half_crossings = int(seed_meta["half_crossings"])
    ydot0_sign = float(seed_meta["ydot0_sign"])
    print(
        f"{log_prefix} parent={parent.label}: continuing direction={direction:+d}, dC={parent.dc}, "
        f"n_steps={parent.n_steps}, window={parent.cj_window}..."
    )
    t_continue_start = time.time()
    try:
        branch = continue_family(
            system,
            seed,
            direction=direction,
            d_jacobi=parent.dc,
            n_steps=parent.n_steps,
            min_jacobi=c_min,
            max_jacobi=c_max,
            half_crossings=half_crossings,
            ydot0_sign=ydot0_sign,
            seed_label=parent.label,
            radau_closure_tol=5e-2,  # Phase 1's loose closure (high-instability families)
            radau_jacobi_tol=1e-7,
            period_step_frac=0.05,
            period_floor_frac=0.5,
            period_ceiling_frac=1.5,
        )
    except (RuntimeError, ValueError) as e:
        print(f"{log_prefix} parent={parent.label}: continuation FAILED — {e}")
        output_fh.write(
            json.dumps(
                {
                    "kind": "parent_error",
                    "parent_label": parent.label,
                    "stage": "continuation",
                    "error": str(e),
                    "elapsed_seconds": time.time() - t_parent_start,
                }
            )
            + "\n"
        )
        counters["errors"] += 1
        return counters
    walk_dt = time.time() - t_continue_start
    counters["members_walked"] = len(branch.members)
    print(
        f"{log_prefix} parent={parent.label}: walk done — stop={branch.stop_reason}, "
        f"members={len(branch.members)}, walk_time={walk_dt:.1f}s"
    )

    # Step 3: build FamilyMembers + detect saddle-center brackets.
    fam_members = [
        FamilyMember(
            label=f"{parent.label}_i{i}",
            state0=np.array([m.x0, 0.0, 0.0, 0.0, m.ydot0, 0.0], dtype=np.float64),
            period=float(m.period),
            mu=float(system.mu),
            parameter=float(m.jacobi),
        )
        for i, m in enumerate(branch.members)
    ]
    if len(fam_members) < 2:
        print(
            f"{log_prefix} parent={parent.label}: "
            f"not enough members for detector ({len(fam_members)})"
        )
        output_fh.write(
            json.dumps(
                {
                    "kind": "parent_summary",
                    "parent_label": parent.label,
                    "stage": "detector_skipped_short_walk",
                    "members_walked": len(fam_members),
                    "stop_reason": str(branch.stop_reason),
                    "elapsed_seconds": time.time() - t_parent_start,
                }
            )
            + "\n"
        )
        return counters
    print(f"{log_prefix} parent={parent.label}: running detector (stm_mode=fixed_path)...")
    t_detect_start = time.time()
    try:
        brackets = detect_saddle_center_bracket(fam_members)
    except (RuntimeError, ValueError) as e:
        print(f"{log_prefix} parent={parent.label}: detector FAILED — {e}")
        output_fh.write(
            json.dumps(
                {
                    "kind": "parent_error",
                    "parent_label": parent.label,
                    "stage": "saddle_center_detector",
                    "error": str(e),
                    "elapsed_seconds": time.time() - t_parent_start,
                }
            )
            + "\n"
        )
        counters["errors"] += 1
        return counters
    detect_dt = time.time() - t_detect_start
    counters["brackets_found"] = len(brackets)
    print(
        f"{log_prefix} parent={parent.label}: detector found {len(brackets)} "
        f"saddle-center bracket(s) in {detect_dt:.1f}s"
    )

    # Step 4: emit a parent summary row regardless of bracket count.
    output_fh.write(
        json.dumps(
            {
                "kind": "parent_summary",
                "parent_label": parent.label,
                "k1": parent.k1,
                "k2": parent.k2,
                "jacobi_anchor": parent.jacobi_anchor,
                "cj_window": list(parent.cj_window),
                "direction": direction,
                "members_walked": len(fam_members),
                "stop_reason": str(branch.stop_reason),
                "n_brackets_found": len(brackets),
                "walk_time_seconds": walk_dt,
                "detect_time_seconds": detect_dt,
                "elapsed_seconds": time.time() - t_parent_start,
                "notes": parent.notes,
            }
        )
        + "\n"
    )

    if not brackets:
        print(f"{log_prefix} parent={parent.label}: NO BRACKETS — honest negative")
        return counters

    # Step 5: branch at each saddle-center bracket's post-side.
    for b_idx, bracket in enumerate(brackets):
        post_bif_member = bracket.members[1]
        parent_state0 = post_bif_member.state0
        parent_period = post_bif_member.period
        print(
            f"{log_prefix} parent={parent.label} bracket={b_idx}: "
            f"branching at C={bracket.extras.get('param_after', 'nan')}, eps={epsilon}..."
        )
        try:
            result = branch_at_saddle_center(system, parent_state0, parent_period, epsilon=epsilon)
        except (RuntimeError, ValueError) as e:
            print(
                f"{log_prefix} parent={parent.label} bracket={b_idx}: branch_at_saddle_center "
                f"raised — {e}"
            )
            print(traceback.format_exc())
            output_fh.write(
                json.dumps(
                    {
                        "kind": "branch_error",
                        "parent_label": parent.label,
                        "bracket_index": b_idx,
                        "stage": "branch_corrector",
                        "error": str(e),
                    }
                )
                + "\n"
            )
            counters["errors"] += 1
            continue

        if result is None:
            print(
                f"{log_prefix} parent={parent.label} bracket={b_idx}: "
                f"branch corrector DID NOT CONVERGE"
            )
            output_fh.write(
                json.dumps(
                    {
                        "kind": "branch_no_converge",
                        "parent_label": parent.label,
                        "parent_k1": parent.k1,
                        "parent_k2": parent.k2,
                        "bracket_index": b_idx,
                        "bracket_param_before": bracket.extras.get("param_before"),
                        "bracket_param_after": bracket.extras.get("param_after"),
                        "epsilon": epsilon,
                    }
                )
                + "\n"
            )
            continue

        counters["branches_converged"] += 1
        bo = result.branched_orbit

        # Independent topology cross-check on the branched orbit.
        b_k1, b_k2, _b_prograde = _topology_check(system.mu, bo.state0, bo.T_TU)
        if (b_k1, b_k2) != (result.branched_topology.k1, result.branched_topology.k2):
            # The corrector's internal topology computation should match the
            # re-run here; if it doesn't, log the divergence as a diagnostic.
            print(
                f"{log_prefix} parent={parent.label} bracket={b_idx}: topology "
                f"re-check disagrees with corrector: corrector "
                f"({result.branched_topology.k1}, {result.branched_topology.k2}) vs "
                f"re-run ({b_k1}, {b_k2})"
            )

        # Branched-orbit Floquet diagnostic.
        try:
            branched_mono = monodromy(system, bo.state0, bo.T_TU, stm_mode="fixed_path")
            branched_eigs = floquet_multipliers(branched_mono)
            max_floquet_mag = float(np.max(np.abs(branched_eigs)))
            sigma_tu = math.log(max_floquet_mag) / bo.T_TU if max_floquet_mag > 1.0 else 0.0
            sigma_d = sigma_tu / TU_DAYS
        except (RuntimeError, ValueError):
            max_floquet_mag = float("nan")
            sigma_tu = float("nan")
            sigma_d = float("nan")

        # Cycler-candidate flag (informational; Phase 4 gauntlet decides).
        cycler_candidate = (
            result.topology_changed and bo.corrector_residual < 1e-10 and (b_k1, b_k2) != (-1, -1)
        )
        if result.topology_changed:
            counters["branches_topology_changed"] += 1

        row = {
            "kind": "branch_record",
            "parent_label": parent.label,
            "parent_k1": parent.k1,
            "parent_k2": parent.k2,
            "parent_jacobi_anchor": parent.jacobi_anchor,
            "parent_state0": [float(x) for x in parent_state0],
            "parent_period_TU": float(parent_period),
            "parent_period_days": float(parent_period * TU_DAYS),
            "parent_topology": [
                int(result.parent_topology.k1),
                int(result.parent_topology.k2),
            ],
            "bracket_index": b_idx,
            "bracket_param_before": bracket.extras.get("param_before"),
            "bracket_param_after": bracket.extras.get("param_after"),
            "bracket_eig_before_real": float(bracket.eig_before.real),
            "bracket_eig_before_imag": float(bracket.eig_before.imag),
            "bracket_eig_after_real": float(bracket.eig_after.real),
            "bracket_eig_after_imag": float(bracket.eig_after.imag),
            "branch_state0": [float(x) for x in bo.state0],
            "branch_period_TU": float(bo.T_TU),
            "branch_period_days": float(bo.T_TU * TU_DAYS),
            "branch_jacobi": float(bo.jacobi),
            "branch_z0": float(bo.state0[2]),
            "branch_zdot0": float(bo.state0[5]),
            "branch_degenerate_planar": bool(bo.degenerate_planar),
            "branch_k1": int(result.branched_topology.k1),
            "branch_k2": int(result.branched_topology.k2),
            "branch_k1_recheck": b_k1,
            "branch_k2_recheck": b_k2,
            "topology_changed": bool(result.topology_changed),
            "corrector_residual": float(bo.corrector_residual),
            "independent_closure_residual": float(bo.independent_closure_residual),
            "max_floquet_mag_branched": max_floquet_mag,
            "sigma_TU_branched": sigma_tu,
            "sigma_d_per_day_branched": sigma_d,
            "epsilon_used": float(result.epsilon),
            "eigenvector_sign_used": int(result.sign),
            "eigenvalue_used_real": float(result.eigenvalue_used.real),
            "eigenvalue_used_imag": float(result.eigenvalue_used.imag),
            "eigenvector_used": [float(x) for x in result.eigenvector_used],
            "cycler_candidate_flag": bool(cycler_candidate),
        }
        output_fh.write(json.dumps(row) + "\n")
        print(
            f"{log_prefix} parent={parent.label} bracket={b_idx}: branched orbit "
            f"converged (residual={bo.corrector_residual:.3e}, T={bo.T_TU * TU_DAYS:.2f}d, "
            f"topology=({result.branched_topology.k1}, {result.branched_topology.k2}), "
            f"cycler_candidate={cycler_candidate})"
        )

    return counters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/floquet_phase2_sweep_results.jsonl"),
        help="output JSONL path",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=5e-4,
        help="perturbation amplitude for branch_at_saddle_center (default 5e-4, Phase 1 value)",
    )
    parser.add_argument(
        "--parents",
        type=str,
        default="",
        help="comma-separated list of parent labels to sweep (default: all 6)",
    )
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.parents:
        wanted = set(args.parents.split(","))
        parents = tuple(p for p in PHASE2_SWEEP_PARENTS if p.label in wanted)
        missing = wanted - {p.label for p in parents}
        if missing:
            raise SystemExit(f"unknown parent labels: {missing}")
    else:
        parents = PHASE2_SWEEP_PARENTS

    t_start = time.time()
    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    header = {
        "kind": "header",
        "phase": "347_phase2_p2_3",
        "iso_timestamp": iso_start,
        "epsilon": args.epsilon,
        "stm_mode": "fixed_path",
        "n_parents": len(parents),
        "parents": [p.label for p in parents],
        "discipline": (
            "Sourced anchors only (Braik-Ross 2026 Table 2 + Ross-RT 2025 C_J_C21). "
            "Phase 2 is DISCOVERY only — emitted rows are 'discovery candidates', "
            "NOT catalogue admissions. Catalogue writeback is Phase 4+ gauntlet work."
        ),
    }
    with args.output.open("w") as fh:
        fh.write(json.dumps(header) + "\n")
        fh.flush()
        total_counters = {
            "members_walked": 0,
            "brackets_found": 0,
            "branches_converged": 0,
            "branches_topology_changed": 0,
            "errors": 0,
        }
        for idx, parent in enumerate(parents):
            log_prefix = f"[{idx + 1}/{len(parents)}]"
            try:
                counters = _process_parent(
                    parent, epsilon=args.epsilon, output_fh=fh, log_prefix=log_prefix
                )
            except (RuntimeError, ValueError) as e:
                print(f"{log_prefix} parent={parent.label}: UNHANDLED — {e}")
                print(traceback.format_exc())
                fh.write(
                    json.dumps(
                        {
                            "kind": "parent_error",
                            "parent_label": parent.label,
                            "stage": "unhandled",
                            "error": str(e),
                        }
                    )
                    + "\n"
                )
                counters = {k: 0 for k in total_counters}
                counters["errors"] = 1
            for k, v in counters.items():
                total_counters[k] += v
            fh.flush()
        wall_dt = time.time() - t_start
        footer = {
            "kind": "footer",
            "phase": "347_phase2_p2_3",
            "iso_end": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "wall_time_seconds": wall_dt,
            "counters": total_counters,
        }
        fh.write(json.dumps(footer) + "\n")
        print(
            f"\nPhase 2 sweep done — wall_time={wall_dt:.1f}s, "
            f"parents={len(parents)}, "
            f"members_walked={total_counters['members_walked']}, "
            f"brackets_found={total_counters['brackets_found']}, "
            f"branches_converged={total_counters['branches_converged']}, "
            f"topology_changed={total_counters['branches_topology_changed']}, "
            f"errors={total_counters['errors']}"
        )
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
