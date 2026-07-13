"""#582 (stage 3a of #581): asymmetric/spatial-isolated 3D CR3BP niching-GA
search -- driver + positive-control harness.

Reuses ``search/niching_ga.py::run_deterministic_crowding`` UNMODIFIED with
the new ``search/isolated_3d_asymmetric_fitness.py`` objective, at mu=0.001,
one run per interior MMR (mirrors #581 stage 2's one-run-per-set pattern).

**This dispatch's deliverable is the positive control ONLY**
(``--positive-control``): reproduce an already-known #440 circular member
through the NEW fitness/bounds (not the old symmetric corrector), refine
through the EXISTING asymmetric 3D corrector
(``cr3bp_general_periodic_3d.correct_general_periodic_3d``), classify its
symmetry, populate its literature-matcher signature, and compare
``(x0, ydot0, T, C)`` against the known member to a stated tolerance. The
full 5-MMR asymmetric NOVELTY sweep (``--mode ga`` with the full paper-scale
GA budget) is a follow-up job, not run by this script's own invocation here.

Usage:
    uv run python scripts/run_582_asymmetric_3d_niching_search.py \\
        --mmr 3:2 --positive-control --workers 8

    # generic (checkpointed) GA run for one MMR, no downstream pipeline:
    uv run python scripts/run_582_asymmetric_3d_niching_search.py \\
        --mmr 3:2 --mode ga --max-gens 50 --workers 8

Checkpoints/results: data/found/582_niching_ga/
"""

from __future__ import annotations

import argparse
import functools
import json
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.preflight import MethodCapability, preflight_search
from cyclerfinder.search.er3bp_isolated_seeds import MMR_SEMI_MAJOR_AXES, resonant_po_seed
from cyclerfinder.search.isolated_3d_asymmetric_fitness import (
    IsolatedAsymmetricFitnessConfig,
    isolated_3d_asymmetric_fitness,
    mmr_bounds,
)
from cyclerfinder.search.isolated_3d_asymmetric_pipeline import (
    build_candidate_signature,
    classify_symmetry,
    literature_anchors_engaged,
    refine_ga_candidate,
)
from cyclerfinder.search.niching_ga import DeterministicCrowdingConfig, run_deterministic_crowding

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "582_niching_ga"

MU = 0.001
MMR_BY_LABEL: dict[str, tuple[int, int, float]] = {
    f"{p}:{q}": (p, q, a1) for p, q, a1 in MMR_SEMI_MAJOR_AXES
}

_METHOD = MethodCapability(
    genome=(
        "isolated_3d_asymmetric_fitness (x0,z0,xdot0,ydot0,zdot0,T genome, y0=0) "
        "evaluated by deterministic-crowding niching GA (search/niching_ga.py) "
        "at mu=0.001 over the 5 tabulated interior MMRs (er3bp_isolated_seeds.py)"
    ),
    corrector=(
        "cr3bp_general_periodic_3d.correct_general_periodic_3d "
        "(full-asymmetric free vars, independent Radau closure check)"
    ),
    capability_tags=frozenset(
        {"cr3bp", "3d", "asymmetric", "niching-ga", "isolated-mmr", "mu-0.001"}
    ),
    git_sha="working-tree",
)


def _system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=MU, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)


def _small_ga_config(seed: int) -> DeterministicCrowdingConfig:
    """Small, fast GA budget for the positive-control harness (not a paper-scale sweep)."""
    return DeterministicCrowdingConfig(population_size=40, generations=60, seed=seed)


def run_ga(
    p: int,
    q: int,
    a1: float,
    *,
    config: DeterministicCrowdingConfig,
    workers: int,
    max_gens: int | None,
    checkpoint: Path | None,
    runlog: Path | None,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Run (or resume) the niching GA for one MMR.

    Returns ``(phenotypes, fitness, x0_guess, ydot0_guess, t0)``.
    """
    bounds, x0_guess, ydot0_guess, t0 = mmr_bounds(a1, mu=MU)
    seed_state = np.array([x0_guess, 0.0, 0.0, 0.0, ydot0_guess, 0.0])
    jacobi_target = float(cr3bp.jacobi_constant(seed_state, MU))
    fit_config = IsolatedAsymmetricFitnessConfig(mu=MU, t0=t0, jacobi_target=jacobi_target)
    fit = functools.partial(isolated_3d_asymmetric_fitness, config=fit_config)

    t_start = time.monotonic()

    def progress(stats: dict[str, float]) -> None:
        rec = {
            "ts": datetime.now(UTC).isoformat(timespec="seconds"),
            "mmr": f"{p}:{q}",
            **stats,
            "elapsed_s": round(time.monotonic() - t_start, 1),
        }
        if runlog is not None:
            with runlog.open("a") as fh:
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
        gen = int(stats["generation"])
        if gen % 10 == 0 or gen == config.generations:
            print(
                f"[{p}:{q}] gen {gen}/{config.generations} "
                f"mean={stats['fitness_mean']:.6f} max={stats['fitness_max']:.6f} "
                f"elapsed={rec['elapsed_s']}s",
                flush=True,
            )

    result = run_deterministic_crowding(
        fit,
        bounds,
        config,
        workers=workers,
        checkpoint_path=checkpoint,
        max_generations_this_call=max_gens,
        progress_fn=progress,
    )
    print(
        f"[{p}:{q}] GA done through generation {result.generations_run}/{config.generations}; "
        f"best fitness={float(np.max(result.fitness)):.6f}",
        flush=True,
    )
    return result.phenotypes, result.fitness, x0_guess, ydot0_guess, t0


def positive_control(p: int, q: int, a1: float, *, workers: int) -> dict[str, object]:
    """End-to-end: GA -> refine -> classify symmetry -> literature signature -> compare.

    Reproduces the known #440 circular member for MMR ``p:q`` USING the new
    fitness function and bounds (not ``correct_symmetric_fixed_jacobi``), then
    checks the result against the known member's ``(x0, ydot0, T, C)``.
    """
    system = _system()
    label = f"{p}:{q}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / f"positive_control_{p}_{q}_checkpoint.npz"
    runlog = OUT_DIR / f"positive_control_{p}_{q}_runlog.jsonl"
    if ckpt.exists():
        ckpt.unlink()  # positive control is a fresh, self-contained run each time

    known = resonant_po_seed(system, p, q, a1)
    print(
        f"[{label}] known #440 member: x0={known.state0[0]:.6f} ydot0={known.state0[4]:.6f} "
        f"T={known.period:.6f} C={known.jacobi:.6f} (converged={known.converged})",
        flush=True,
    )

    config = _small_ga_config(seed=582000 + p * 100 + q)
    phen, fitness, _x0_guess, _ydot0_guess, _t0 = run_ga(
        p, q, a1, config=config, workers=workers, max_gens=None, checkpoint=ckpt, runlog=runlog
    )

    best_i = int(np.argmax(fitness))
    best_genome = phen[best_i]
    best_fit = float(fitness[best_i])
    print(
        f"[{label}] GA best genome: x0={best_genome[0]:.6f} z0={best_genome[1]:.6f} "
        f"xdot0={best_genome[2]:.6f} ydot0={best_genome[3]:.6f} zdot0={best_genome[4]:.6f} "
        f"T={best_genome[5]:.6f} fitness={best_fit:.6f}",
        flush=True,
    )

    orbit = refine_ga_candidate(system, best_genome)
    print(
        f"[{label}] corrector: converged={orbit.converged} "
        f"corrector_residual={orbit.corrector_residual:.3e} "
        f"independent_closure={orbit.independent_closure_residual:.3e} "
        f"state0={orbit.state0} T={orbit.T_TU:.6f} C={orbit.jacobi:.6f} "
        f"degenerate_planar={orbit.degenerate_planar}",
        flush=True,
    )

    symmetry = classify_symmetry(system, orbit.state0, orbit.T_TU)
    print(
        f"[{label}] symmetry: is_symmetric={symmetry.is_symmetric} "
        f"best_crossing_residual={symmetry.best_crossing_residual:.3e} "
        f"(n_crossings_checked={symmetry.n_crossings_checked})",
        flush=True,
    )

    sig = build_candidate_signature(system, orbit, p=p, q=q)
    anchors = literature_anchors_engaged(sig)
    print(f"[{label}] literature_check anchors engaged: {anchors}", flush=True)

    # Tolerances (stated explicitly; see the results note for the achieved margin).
    # The full-asymmetric corrector is 6 residuals / 7 unknowns (min-norm
    # least-squares, cr3bp_general_periodic_3d's own documented behaviour): it
    # lands on the CLOSEST periodic orbit to the GA seed, which need not be
    # bit-identical to the #440 seed's own symmetric-corrector member -- a
    # few-percent family-member shift is expected, not a bug.
    tol_x0_rel = 0.03
    tol_ydot0_rel = 0.05
    tol_t_rel = 0.05
    tol_jacobi_abs = 0.02

    x0_err = abs(orbit.state0[0] - known.state0[0]) / abs(known.state0[0])
    ydot0_err = abs(orbit.state0[4] - known.state0[4]) / abs(known.state0[4])
    t_err = abs(orbit.T_TU - known.period) / abs(known.period)
    jacobi_err = abs(orbit.jacobi - known.jacobi)

    matched = (
        orbit.converged
        and x0_err < tol_x0_rel
        and ydot0_err < tol_ydot0_rel
        and t_err < tol_t_rel
        and jacobi_err < tol_jacobi_abs
    )

    summary = {
        "mmr": label,
        "known_member": {
            "x0": float(known.state0[0]),
            "ydot0": float(known.state0[4]),
            "T": float(known.period),
            "C": float(known.jacobi),
            "converged": bool(known.converged),
        },
        "ga_best_genome": [float(v) for v in best_genome],
        "ga_best_fitness": best_fit,
        "refined_orbit": {
            "state0": [float(v) for v in orbit.state0],
            "T": float(orbit.T_TU),
            "C": float(orbit.jacobi),
            "converged": bool(orbit.converged),
            "corrector_residual": float(orbit.corrector_residual),
            "independent_closure_residual": float(orbit.independent_closure_residual),
            "degenerate_planar": bool(orbit.degenerate_planar),
        },
        "symmetry": {
            "is_symmetric": bool(symmetry.is_symmetric),
            "best_crossing_residual": float(symmetry.best_crossing_residual),
            "n_crossings_checked": symmetry.n_crossings_checked,
        },
        "literature_anchors_engaged": anchors,
        "tolerances": {
            "x0_rel": tol_x0_rel,
            "ydot0_rel": tol_ydot0_rel,
            "T_rel": tol_t_rel,
            "jacobi_abs": tol_jacobi_abs,
        },
        "errors": {
            "x0_rel": float(x0_err),
            "ydot0_rel": float(ydot0_err),
            "T_rel": float(t_err),
            "jacobi_abs": float(jacobi_err),
        },
        "matched": bool(matched),
    }

    out = OUT_DIR / f"positive_control_{p}_{q}_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    verdict = "PASS" if matched else "FAIL"
    print(
        f"[{label}] POSITIVE CONTROL {verdict}: "
        f"x0_err={x0_err:.4%} ydot0_err={ydot0_err:.4%} T_err={t_err:.4%} "
        f"jacobi_err={jacobi_err:.2e}",
        flush=True,
    )
    print(f"[{label}] summary written to {out}", flush=True)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mmr", required=True, choices=sorted(MMR_BY_LABEL), help="p:q interior MMR")
    ap.add_argument(
        "--mode",
        choices=("positive-control", "ga"),
        default="positive-control",
        help="'positive-control' runs the full GA->corrector->classify->signature "
        "pipeline and compares to the known #440 member; 'ga' just runs (or "
        "resumes) the checkpointed niching GA with no downstream pipeline",
    )
    ap.add_argument(
        "--positive-control", action="store_true", help="alias for --mode positive-control"
    )
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--max-gens", type=int, default=None, help="generation cap this call (--mode ga only)"
    )
    ap.add_argument(
        "--population", type=int, default=200, help="GA population size (--mode ga only)"
    )
    ap.add_argument(
        "--generations", type=int, default=400, help="GA generation budget (--mode ga only)"
    )
    args = ap.parse_args()

    p, q, a1 = MMR_BY_LABEL[args.mmr]
    mode = "positive-control" if args.positive_control else args.mode

    region_id = f"582-asymmetric-3d-isolated-mmr-{args.mmr.replace(':', '-')}-{mode}-2026-07-13"
    preflight_search(
        task_no=582,
        region_id=region_id,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=1,
        override_reason=(
            "positive-control reproduction of an already-known #440 circular "
            "MMR member (or its small-budget GA precursor) -- validates the "
            "new isolated_3d_asymmetric_fitness + pipeline against a known "
            "answer, not the full 5-MMR novelty sweep #582 explicitly defers "
            "to a separate, coordinator-owned dispatch"
        ),
    )

    if mode == "positive-control":
        positive_control(p, q, a1, workers=args.workers)
    else:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        ckpt = OUT_DIR / f"ga_{p}_{q}_checkpoint.npz"
        runlog = OUT_DIR / f"ga_{p}_{q}_runlog.jsonl"
        config = DeterministicCrowdingConfig(
            population_size=args.population, generations=args.generations, seed=582000 + p * 100 + q
        )
        phen, fitness, *_ = run_ga(
            p,
            q,
            a1,
            config=config,
            workers=args.workers,
            max_gens=args.max_gens,
            checkpoint=ckpt,
            runlog=runlog,
        )
        if fitness.shape[0] and np.max(fitness) >= 0:
            np.savez(OUT_DIR / f"ga_{p}_{q}_final.npz", phenotypes=phen, fitness=fitness)


if __name__ == "__main__":
    main()
