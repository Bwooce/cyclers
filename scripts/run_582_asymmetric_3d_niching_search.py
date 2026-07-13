"""#582 (stage 3a of #581): asymmetric/spatial-isolated 3D CR3BP niching-GA
search -- driver + positive-control harness + population analysis.

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

**``--mode analyze`` (added for the analysis follow-up to the above sweep)**
turns one finished ``ga_{p}_{q}_final.npz`` into an adjudicated per-MMR
summary, following the same pattern as
``run_583_widened_bounded_drift_search.py``'s ``analyze_partition``:

1. **Clustering.** The niching GA is designed to hold multiple co-existing
   basins in one final population -- taking only the single best-fitness
   member would silently collapse that structure. :func:`cluster_representatives`
   greedily selects cluster representatives: sort the population by fitness
   descending, then walk down the list adding a member as a NEW cluster
   representative only if it is more than ``distance_threshold`` away (in
   the SAME bounds-normalized Euclidean metric ``niching_ga.py``'s own
   parent-child pairing rule uses -- see that module's docstring) from every
   representative already picked. This is deterministic-crowding's own
   distance metric turned into a threshold-based greedy clustering, not a
   new metric invented for this step; it needs no new dependency (no
   scikit-learn -- checked ``pyproject.toml``, not a project dependency and
   this project prefers to avoid adding one for a step this simple).
2. **Fitness floor 0.9** (default, matching #583's own ``analyze_partition``
   choice for the exact same reason): ``isolated_3d_asymmetric_fitness`` is a
   bounded-reciprocal-penalty objective peaking at 1.0 for a near-periodic,
   near-target-Jacobi genome and falling off fast away from any basin (see
   that module's docstring); 0.9 keeps only genomes close enough to a real
   periodicity/Jacobi basin to be worth the corrector's cost, while still
   being well below the ``converged`` corrector's own much tighter closure
   tolerances -- the corrector step immediately downstream is the actual
   accept/reject gate, not this floor.
3. **Distance threshold 0.1** (default, in bounds-normalized ``[0,1]^6``
   units per free gene): coarse enough that near-duplicate members of the
   SAME niche (which a 200-member population typically produces many of
   near its best point) collapse to one representative, tight enough that
   two genuinely different resonant-family basins inside one MMR's box
   (which spans only +-15-50% around the analytic e=0 guess per
   ``mmr_bounds``) are not accidentally merged. Not paper-sourced -- there is
   no published clustering recipe for this new fitness landscape -- but
   documented and adjustable via ``--distance-threshold``.
4. **Pipeline per representative**: EXACTLY the positive control's own
   sequence -- ``isolated_3d_asymmetric_pipeline``'s ``refine_ga_candidate``
   -> check ``orbit.converged`` -> ``classify_symmetry`` ->
   ``build_candidate_signature`` -> ``literature_anchors_engaged`` -- no new
   corrector or classifier.
5. **Explicitly NOT done here**: the live literature-corpus search
   (``search/literature_check.py::check_literature``, which needs an
   injected ``SearchFn``) and any novelty adjudication -- per
   [[feedback_literature_novelty_check_baseline]] this analysis pass only
   gets a candidate to "converged, classified, matcher-engaged", which is
   necessary but not sufficient for any "novel" claim. The summary JSON says
   this explicitly (``literature_check_status`` field) rather than leaving
   it to be inferred.

**#585 follow-up (2026-07-14): resonance-scaled symmetry-breaking bounds.**
``mmr_bounds()``'s 3 flat absolute half-widths (``z0_abs=xdot0_abs=zdot0_abs=
0.05`` for every MMR) were replaced by a single resonance-scaled fraction
``s`` (Fable-reviewed GO design; see that function's docstring for the exact
formula). Every mode now requires ``--symmetry-breaking-s`` explicitly --
there is no default, so a caller can never silently inherit a changed box.
The Fable-reviewed 2-rung ladder is ``s=0.15`` and ``s=0.30``. Output
filenames are tagged with the ``s`` value used (e.g. ``s0p15``) so runs at
different rungs never collide or silently overwrite each other, and so
#582's own untagged ``ga_{p}_{q}_final.npz``/``{p}_{q}_analysis_summary.json``
files (produced under the old flat-0.05 box) are left untouched.

``--mode analyze`` also now runs a drift-detection check (#585 scope item
2): for each converged cluster representative, the converged period is
inverted back to an implied semi-major axis
(``isolated_3d_asymmetric_fitness.mmr_a1_from_t0``, the exact inverse of the
``mmr_t0`` formula ``mmr_bounds`` itself uses) and compared against every
tabulated interior MMR's own ``a1`` (``MMR_SEMI_MAJOR_AXES``) -- if the
NEAREST MMR by that implied ``a1`` is not this run's own target MMR, the
representative is flagged ``drifted_to_neighboring_mmr=true``. The summary
JSON reports a per-run ``n_drifted_to_neighboring_mmr`` count. This does not
reject or filter drifted candidates (a drifted-but-converged orbit may still
be a genuine member of ITS OWN resonance's family, just not this run's
target) -- it makes over-widening detectable rather than silently folded
into "this MMR's" cluster count.

Usage:
    uv run python scripts/run_582_asymmetric_3d_niching_search.py \\
        --mmr 3:2 --positive-control --symmetry-breaking-s 0.15 --workers 8

    # generic (checkpointed) GA run for one MMR, no downstream pipeline:
    uv run python scripts/run_582_asymmetric_3d_niching_search.py \\
        --mmr 3:2 --mode ga --symmetry-breaking-s 0.15 --max-gens 50 --workers 8

    # analyze one MMR's finished final population (SAME s used to run it):
    uv run python scripts/run_582_asymmetric_3d_niching_search.py \\
        --mmr 3:2 --mode analyze --symmetry-breaking-s 0.15

Checkpoints/results: data/found/582_niching_ga/
"""

from __future__ import annotations

import argparse
import functools
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.preflight import MethodCapability, preflight_search
from cyclerfinder.search.er3bp_isolated_seeds import MMR_SEMI_MAJOR_AXES, resonant_po_seed
from cyclerfinder.search.isolated_3d_asymmetric_fitness import (
    IsolatedAsymmetricFitnessConfig,
    isolated_3d_asymmetric_fitness,
    mmr_a1_from_t0,
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


def _s_tag(s: float) -> str:
    """Filesystem-safe tag for the #585 symmetry-breaking fraction, e.g. 0.15 -> 's0p15'."""
    return "s" + f"{s:g}".replace(".", "p").replace("-", "m")


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
    s: float,
    config: DeterministicCrowdingConfig,
    workers: int,
    max_gens: int | None,
    checkpoint: Path | None,
    runlog: Path | None,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Run (or resume) the niching GA for one MMR.

    ``s`` is #585's resonance-scaled symmetry-breaking fraction, threaded
    straight through to :func:`mmr_bounds` (see that function's docstring).

    Returns ``(phenotypes, fitness, x0_guess, ydot0_guess, t0)``.
    """
    bounds, x0_guess, ydot0_guess, t0 = mmr_bounds(a1, mu=MU, s=s)
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


def positive_control(p: int, q: int, a1: float, *, s: float, workers: int) -> dict[str, object]:
    """End-to-end: GA -> refine -> classify symmetry -> literature signature -> compare.

    Reproduces the known #440 circular member for MMR ``p:q`` USING the new
    fitness function and bounds (not ``correct_symmetric_fixed_jacobi``), then
    checks the result against the known member's ``(x0, ydot0, T, C)``.

    ``s`` is #585's resonance-scaled symmetry-breaking fraction (see
    ``mmr_bounds``'s docstring); the mandatory #585 gate is re-running this
    positive control at EACH new rung (``s=0.15``, ``s=0.30``) and confirming
    it still recovers the known member before any novelty sweep at that rung.
    """
    system = _system()
    label = f"{p}:{q}"
    tag = _s_tag(s)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / f"positive_control_{p}_{q}_{tag}_checkpoint.npz"
    runlog = OUT_DIR / f"positive_control_{p}_{q}_{tag}_runlog.jsonl"
    if ckpt.exists():
        ckpt.unlink()  # positive control is a fresh, self-contained run each time

    known = resonant_po_seed(system, p, q, a1)
    print(
        f"[{label}] s={s} known #440 member: x0={known.state0[0]:.6f} ydot0={known.state0[4]:.6f} "
        f"T={known.period:.6f} C={known.jacobi:.6f} (converged={known.converged})",
        flush=True,
    )

    config = _small_ga_config(seed=582000 + p * 100 + q)
    phen, fitness, _x0_guess, _ydot0_guess, _t0 = run_ga(
        p, q, a1, s=s, config=config, workers=workers, max_gens=None, checkpoint=ckpt, runlog=runlog
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
        "symmetry_breaking_s": s,
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

    out = OUT_DIR / f"positive_control_{p}_{q}_{tag}_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    verdict = "PASS" if matched else "FAIL"
    print(
        f"[{label}] s={s} POSITIVE CONTROL {verdict}: "
        f"x0_err={x0_err:.4%} ydot0_err={ydot0_err:.4%} T_err={t_err:.4%} "
        f"jacobi_err={jacobi_err:.2e}",
        flush=True,
    )
    print(f"[{label}] summary written to {out}", flush=True)
    return summary


# ---------------------------------------------------------------------------
# --mode analyze: cluster a finished GA population, run the mandatory
# refine -> classify -> signature -> anchor pipeline on each representative.
# ---------------------------------------------------------------------------

DEFAULT_ANALYZE_FITNESS_FLOOR = 0.9
DEFAULT_ANALYZE_DISTANCE_THRESHOLD = 0.1
DEFAULT_ANALYZE_MAX_CLUSTERS = 25


def cluster_representatives(
    phen: NDArray[np.float64],
    fitness: NDArray[np.float64],
    bounds: list[tuple[float, float]],
    *,
    fitness_floor: float,
    distance_threshold: float,
    max_clusters: int,
) -> list[int]:
    """Greedy deterministic-crowding-style niche clustering.

    Walks the population in DESCENDING fitness order (stopping once fitness
    drops below ``fitness_floor``) and accepts a member as a NEW cluster
    representative only if it is more than ``distance_threshold`` away, in
    the bounds-normalized ``[0,1]^6`` Euclidean metric ``niching_ga.py``'s
    own parent-child pairing rule uses (see that module's ``_dist``/
    ``_Encoding.decode_norm``), from EVERY representative already accepted.
    This keeps near-duplicate members of one niche from each producing their
    own "cluster" while still separating genuinely distinct basins. Returns
    population indices (into ``phen``/``fitness``), highest-fitness first.
    """
    lo = np.array([b[0] for b in bounds], dtype=np.float64)
    span = np.array([b[1] - b[0] for b in bounds], dtype=np.float64)
    span[span == 0.0] = 1.0  # defensive only -- this genome has no fixed genes
    norm = (phen - lo) / span

    order = np.argsort(-fitness)
    reps: list[int] = []
    rep_vecs: list[NDArray[np.float64]] = []
    for idx in order:
        if float(fitness[idx]) < fitness_floor:
            break  # order is descending, so nothing further qualifies either
        v = norm[idx]
        if all(float(np.linalg.norm(v - rv)) > distance_threshold for rv in rep_vecs):
            reps.append(int(idx))
            rep_vecs.append(v)
            if len(reps) >= max_clusters:
                break
    return reps


def nearest_mmr_by_implied_a1(period: float) -> tuple[str, float, float]:
    """#585 drift-detection: which tabulated MMR is closest to this period's implied a1?

    Inverts the converged period back to an implied semi-major axis via
    :func:`mmr_a1_from_t0` (the exact inverse of the ``mmr_t0`` formula
    ``mmr_bounds`` itself is centered on), then returns
    ``(nearest_label, a1_implied, dist_to_nearest)`` where ``nearest_label``
    is the ``p:q`` MMR whose own tabulated ``a1`` (``MMR_SEMI_MAJOR_AXES``)
    is closest to ``a1_implied``.
    """
    a1_implied = mmr_a1_from_t0(period)
    best_label = ""
    best_dist = math.inf
    for pp, qq, aa in MMR_SEMI_MAJOR_AXES:
        dist = abs(a1_implied - aa)
        if dist < best_dist:
            best_dist = dist
            best_label = f"{pp}:{qq}"
    return best_label, a1_implied, best_dist


def analyze_ga_population(
    p: int,
    q: int,
    a1: float,
    *,
    s: float,
    fitness_floor: float = DEFAULT_ANALYZE_FITNESS_FLOOR,
    distance_threshold: float = DEFAULT_ANALYZE_DISTANCE_THRESHOLD,
    max_clusters: int = DEFAULT_ANALYZE_MAX_CLUSTERS,
) -> dict[str, object]:
    """Cluster ``ga_{p}_{q}_{s-tag}_final.npz`` and run the mandatory pipeline on each niche.

    Per #582's own mandate (module docstring), a GA fitness peak is a basin
    INDICATOR, not a converged orbit: every cluster representative is routed
    through :func:`refine_ga_candidate` -> :func:`classify_symmetry` ->
    :func:`build_candidate_signature` -> :func:`literature_anchors_engaged`,
    the SAME sequence the positive control above already validated. Writes
    ``data/found/582_niching_ga/{p}_{q}_{s-tag}_analysis_summary.json`` and
    returns the same summary dict. Does NOT run any live literature search --
    see the module docstring and the summary's own ``literature_check_status``
    field.

    ``s`` MUST match the ``--symmetry-breaking-s`` the GA population being
    analyzed was actually run with (it re-derives ``mmr_bounds`` purely for
    the clustering normalization scale, so a mismatched ``s`` would silently
    use the wrong per-gene normalization). Every converged, symmetric-or-not
    representative also gets #585's drift-detection check (see
    :func:`nearest_mmr_by_implied_a1`): ``drifted_to_neighboring_mmr=true``
    means the representative's OWN converged period implies a semi-major
    axis nearer a DIFFERENT tabulated MMR than this run's target -- a sign
    the widened box let the GA wander into a neighboring resonance's basin
    rather than this one's (#440's own documented risk). Drifted candidates
    are still reported, not discarded -- they may be genuine members of
    THEIR OWN resonance, just not evidence about the one this run targeted.
    """
    label = f"{p}:{q}"
    tag = _s_tag(s)
    fpath = OUT_DIR / f"ga_{p}_{q}_{tag}_final.npz"
    if not fpath.exists():
        print(f"[{label}] final population missing ({fpath}); run --mode ga to completion first")
        return {"mmr": label, "error": "final population missing", "path": str(fpath)}

    data = np.load(fpath)
    phen, fitness = data["phenotypes"], data["fitness"]
    bounds, _x0_guess, _ydot0_guess, _t0 = mmr_bounds(a1, mu=MU, s=s)

    rep_idx = cluster_representatives(
        phen,
        fitness,
        bounds,
        fitness_floor=fitness_floor,
        distance_threshold=distance_threshold,
        max_clusters=max_clusters,
    )
    n_above_floor = int(np.sum(fitness >= fitness_floor))
    print(
        f"[{label}] population={phen.shape[0]} n_above_floor({fitness_floor})={n_above_floor} "
        f"clusters_selected={len(rep_idx)}",
        flush=True,
    )

    system = _system()
    clusters: list[dict[str, object]] = []
    n_converged = 0
    n_ready = 0
    n_drifted = 0
    for rank, idx in enumerate(rep_idx):
        genome = phen[idx]
        fit = float(fitness[idx])
        entry: dict[str, object] = {
            "cluster_rank": rank,
            "population_index": idx,
            "ga_genome": [float(v) for v in genome],
            "ga_fitness": fit,
        }
        try:
            orbit = refine_ga_candidate(system, genome)
        except Exception as exc:  # one bad seed must not kill the whole analysis pass
            entry["refine_error"] = str(exc)
            clusters.append(entry)
            print(f"[{label}] cluster {rank}: REFINE ERROR ({exc})", flush=True)
            continue

        entry["refined_orbit"] = {
            "state0": [float(v) for v in orbit.state0],
            "T": float(orbit.T_TU),
            "C": float(orbit.jacobi),
            "converged": bool(orbit.converged),
            "corrector_residual": float(orbit.corrector_residual),
            "independent_closure_residual": float(orbit.independent_closure_residual),
            "degenerate_planar": bool(orbit.degenerate_planar),
        }
        if not orbit.converged:
            entry["ready_for_literature_check"] = False
            print(
                f"[{label}] cluster {rank}: NOT CONVERGED fitness={fit:.6f} "
                f"corrector_residual={orbit.corrector_residual:.3e} "
                f"independent_closure={orbit.independent_closure_residual:.3e}",
                flush=True,
            )
            clusters.append(entry)
            continue

        symmetry = classify_symmetry(system, orbit.state0, orbit.T_TU)
        entry["symmetry"] = {
            "is_symmetric": bool(symmetry.is_symmetric),
            "best_crossing_residual": float(symmetry.best_crossing_residual),
            "n_crossings_checked": symmetry.n_crossings_checked,
        }

        # #585 scope item 2: drift-detection -- does this representative's own
        # converged period imply a semi-major axis nearer a NEIGHBORING MMR
        # than this run's own target? Widening the symmetry-breaking bounds
        # (via `s`) risks exactly #440's own documented neighboring-MMR /
        # exterior-1:2 family-selection trap; this makes it detectable rather
        # than silently folded into "this MMR's" cluster count.
        nearest_label, a1_implied, nearest_dist = nearest_mmr_by_implied_a1(orbit.T_TU)
        own_dist = abs(a1_implied - a1)
        drifted = nearest_label != label
        entry["drift_check"] = {
            "a1_implied": a1_implied,
            "own_mmr_a1_dist": own_dist,
            "nearest_mmr": nearest_label,
            "nearest_mmr_a1_dist": nearest_dist,
            "drifted_to_neighboring_mmr": drifted,
        }
        n_drifted += int(drifted)

        sig = build_candidate_signature(system, orbit, p=p, q=q)
        anchors = literature_anchors_engaged(sig)
        entry["candidate_signature"] = {
            "primary": sig.primary,
            "sequence": list(sig.sequence),
            "resonances": list(sig.resonances),
            "topology_label": sorted(sig.topology_label),
            "topology_3d": sig.topology_3d,
        }
        ready = bool(anchors)
        entry["literature_anchors_engaged"] = anchors
        entry["ready_for_literature_check"] = ready
        n_converged += 1
        n_ready += int(ready)
        classification = "symmetric" if symmetry.is_symmetric else "asymmetric"
        drift_note = f" DRIFTED->{nearest_label}" if drifted else ""
        print(
            f"[{label}] cluster {rank}: CONVERGED fitness={fit:.6f} {classification}{drift_note} "
            f"independent_closure={orbit.independent_closure_residual:.2e} "
            f"anchors={anchors}",
            flush=True,
        )
        clusters.append(entry)

    summary: dict[str, object] = {
        "mmr": label,
        "symmetry_breaking_s": s,
        "population_size": int(phen.shape[0]),
        "fitness_floor": fitness_floor,
        "distance_threshold": distance_threshold,
        "n_above_floor": n_above_floor,
        "n_clusters": len(clusters),
        "n_converged": n_converged,
        "n_ready_for_literature_check": n_ready,
        "n_drifted_to_neighboring_mmr": n_drifted,
        "clusters": clusters,
        "literature_check_status": (
            "NOT RUN. Per [[feedback_literature_novelty_check_baseline]], "
            "search/literature_check.py::check_literature() requires a live "
            "injected SearchFn and MUST be run against every cluster with "
            "ready_for_literature_check=true before any 'novel' claim -- this "
            "analysis pass only confirms the structural matcher engages "
            "(non-empty anchor pool via literature_anchors_engaged()), it does "
            "not itself search or adjudicate novelty. Not-found is necessary, "
            "not sufficient."
        ),
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{p}_{q}_{tag}_analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    print(
        f"[{label}] analysis: {len(clusters)} cluster(s), {n_converged} converged, "
        f"{n_ready} literature-matcher-ready, {n_drifted} drifted-to-neighboring-mmr; "
        f"summary -> {out}",
        flush=True,
    )
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mmr", required=True, choices=sorted(MMR_BY_LABEL), help="p:q interior MMR")
    ap.add_argument(
        "--mode",
        choices=("positive-control", "ga", "analyze"),
        default="positive-control",
        help="'positive-control' runs the full GA->corrector->classify->signature "
        "pipeline and compares to the known #440 member; 'ga' just runs (or "
        "resumes) the checkpointed niching GA with no downstream pipeline; "
        "'analyze' clusters a finished ga_{p}_{q}_final.npz and runs the same "
        "mandatory pipeline on each cluster representative",
    )
    ap.add_argument(
        "--positive-control", action="store_true", help="alias for --mode positive-control"
    )
    ap.add_argument(
        "--symmetry-breaking-s",
        type=float,
        required=True,
        dest="symmetry_breaking_s",
        help="#585's resonance-scaled symmetry-breaking fraction, threaded into "
        "mmr_bounds() as xdot0_abs=zdot0_abs=s*v_circ, z0_abs=max(0.05, s*a1). "
        "No default -- every invocation must pick a value explicitly (the "
        "Fable-reviewed ladder is 0.15 then 0.30). Required for all 3 modes, "
        "including --mode analyze (it must match the s the population being "
        "analyzed was actually run with).",
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
    ap.add_argument(
        "--fitness-floor",
        type=float,
        default=DEFAULT_ANALYZE_FITNESS_FLOOR,
        help="minimum GA fitness to be considered for clustering (--mode analyze only)",
    )
    ap.add_argument(
        "--distance-threshold",
        type=float,
        default=DEFAULT_ANALYZE_DISTANCE_THRESHOLD,
        help="bounds-normalized Euclidean distance a candidate must exceed from every "
        "existing representative to start a new cluster (--mode analyze only)",
    )
    ap.add_argument(
        "--max-clusters",
        type=int,
        default=DEFAULT_ANALYZE_MAX_CLUSTERS,
        help="cap on cluster representatives analyzed (--mode analyze only)",
    )
    args = ap.parse_args()

    p, q, a1 = MMR_BY_LABEL[args.mmr]
    mode = "positive-control" if args.positive_control else args.mode
    s = args.symmetry_breaking_s
    tag = _s_tag(s)

    region_id = (
        f"585-asymmetric-3d-isolated-mmr-{args.mmr.replace(':', '-')}-{mode}-{tag}-2026-07-14"
    )
    preflight_search(
        # task_no stays 582 (not 585): preflight's filename-consistency check
        # requires task_no to match this script's own run_582_*.py filename;
        # the #585 attribution for this specific invocation is carried in
        # region_id and override_reason instead.
        task_no=582,
        region_id=region_id,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=1,
        override_reason=(
            "#585 follow-up to #582's 0/104 asymmetric result: re-runs the SAME "
            "positive-control reproduction of an already-known #440 circular "
            "MMR member (or its small-budget GA precursor) at a resonance-"
            "scaled, wider symmetry-breaking bounds box (mmr_bounds(s=...)) -- "
            "validates the widened box didn't break basin containment before "
            "any novelty claim; 'analyze' mode is pure post-processing of an "
            "already-finished checkpoint (clustering + the same mandatory "
            "pipeline plus the new #585 drift-detection check), not a new "
            "search; the full 5-MMR novelty sweep itself (--mode ga at paper "
            "scale, at either s=0.15 or s=0.30) is a separate, coordinator-"
            "owned dispatch, not run by this invocation"
        ),
    )

    if mode == "positive-control":
        positive_control(p, q, a1, s=s, workers=args.workers)
    elif mode == "analyze":
        analyze_ga_population(
            p,
            q,
            a1,
            s=s,
            fitness_floor=args.fitness_floor,
            distance_threshold=args.distance_threshold,
            max_clusters=args.max_clusters,
        )
    else:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        ckpt = OUT_DIR / f"ga_{p}_{q}_{tag}_checkpoint.npz"
        runlog = OUT_DIR / f"ga_{p}_{q}_{tag}_runlog.jsonl"
        config = DeterministicCrowdingConfig(
            population_size=args.population, generations=args.generations, seed=582000 + p * 100 + q
        )
        phen, fitness, *_ = run_ga(
            p,
            q,
            a1,
            s=s,
            config=config,
            workers=args.workers,
            max_gens=args.max_gens,
            checkpoint=ckpt,
            runlog=runlog,
        )
        if fitness.shape[0] and np.max(fitness) >= 0:
            np.savez(OUT_DIR / f"ga_{p}_{q}_{tag}_final.npz", phenotypes=phen, fitness=fitness)


if __name__ == "__main__":
    main()
