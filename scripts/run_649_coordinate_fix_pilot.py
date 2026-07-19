#!/usr/bin/env python3
"""#649: coordinate-transform cross-mu rescue pilot for #628's generative seed model.

`#642` proved `#628`'s Earth-Moon-trained generative model's raw
``(state0, period)`` output collapses onto degenerate L4/L5 equilibria at
BOTH mu targets `#624` originally tested (mu=0.001, Sun-Earth mu~3.0e-6) --
a STRUCTURAL failure (the model has no mu-conditioning), not a sampling
artifact. `#645`'s shortlist item 4 proposed a CHEAP, non-ML test before
treating cross-mu transfer as permanently out of reach: reinterpret the raw
Earth-Moon-shaped output through a mu-INDEPENDENT coordinate (`#629`'s own
rho=(C-3)/(C_L1(mu)-3) scaled-energy quantity, plus a Hill-radius length
scale for position/velocity) and invert that coordinate at the target mu,
BEFORE handing the result to the corrector -- rather than feeding the raw
guess straight in as `#624` did.

This script reuses `#624`'s EXACT evaluation protocol: same N=100, same two
mu targets (mu=0.001 and Sun-Earth), the SAME already-trained #608 model
(re-derived bit-for-bit, not retrained), the SAME uniform-bounding-box
baseline construction, and the SAME ``is_physically_sane`` yardstick
(already `#642`-fixed to reject degenerate equilibria by default). The ONLY
change is that the "generated" arm's raw model samples are passed through
``cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu``
before refinement, instead of straight into ``correct_periodic`` like #624's
original script did. Baseline and transformed-generated draws share one rng
stream at each mu (mirroring #624's own ``_evaluate_mu``), so the comparison
is apples-to-apples within this single run -- #642's own already-recorded
corrected baseline numbers (see this script's own printed cross-reference at
the end) are reused for CONTEXT only, not as a substitute for this run's own
baseline measurement.

**This is a decisive test, not an optimization pass -- report the actual
measured numbers plainly, whichever way they go** (per this task's own
`data/OUTSTANDING.md` `#649` bullet, step 5). See
``src/cyclerfinder/ml/cross_mu_coordinate_transform.py``'s own module
docstring for the exact, explicitly-documented scaling choices.

No catalogue writeback, no literature-novelty check -- capability-transfer
evaluation only, same scope as `#624`.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.ml.cross_mu_coordinate_transform import transform_seed_to_target_mu
from cyclerfinder.ml.orbit_generative import (
    FEATURE_NAMES,
    assemble_corpus,
    fit_clustered_gaussian,
    is_physically_sane,
    uniform_bounding_box_sample,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUTCOME_LOG_DIR = _REPO_ROOT / "out" / "outcome_log"
_OUT_DIR = _REPO_ROOT / "data" / "found" / "649_coordinate_fix_pilot"

N_GENERATE = 100  # SAME N #624 used
SEED = 608  # SAME seed as #608/#624 -- bit-for-bit identical re-derived model
TRAINING_MU = 0.01215  # #628's TRAINING_MU

# #642's own already-recorded corrected baseline numbers, reused here for
# CONTEXT ONLY (this script measures its own fresh in-run baseline for the
# actual ratio -- see module docstring).
_642_CONTEXT = {
    "mu=0.001": {
        "archival_baseline_sane_rate": 0.0,  # 0/100
        "live_rerun_baseline_sane_rate": 2.0 / 60.0,  # 2/60
        "raw_generated_sane_rate": 0.0,  # 0/60 and 0/100, both trials
    },
    "Sun-Earth": {
        "archival_baseline_sane_rate": 1.0 / 100.0,  # 1/100
        "live_rerun_baseline_sane_rate": 2.0 / 60.0,  # 2/60
        "raw_generated_sane_rate": 0.0,  # 0/60 and 0/100, both trials
    },
}


def _corpus_paths() -> list[Path]:
    import glob

    patterns = ["search_campaign*.jsonl", "jpl_corpus.jsonl"]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(sorted(Path(p) for p in glob.glob(str(_OUTCOME_LOG_DIR / pattern))))
    if not paths:
        raise SystemExit(
            f"no outcome-log files found under {_OUTCOME_LOG_DIR} -- this script "
            "expects the project's own gitignored out/outcome_log/ corpus to be "
            "present locally (it is not committed)."
        )
    return paths


def _refine_batch(system: cr3bp.CR3BPSystem, rows: list[dict], *, tag: str) -> list[dict]:
    """Run each row's (state0, period) through the existing corrector (SAME as
    #608/#624/#642). ``rows`` are dicts already carrying whatever
    construction metadata (e.g. transform rho/scale, or None-construction
    failure) the caller wants preserved alongside the refinement outcome.
    """
    results = []
    for i, row in enumerate(rows):
        result = {"tag": tag, "index": i, **{k: v for k, v in row.items() if k != "state0_period"}}
        state0_period = row.get("state0_period")
        if state0_period is None:
            # construction failure upstream (e.g. transform returned None) --
            # counted as a non-converged attempt, matching how #624's own
            # loop counts a corrector exception: never silently dropped from
            # the denominator.
            result.update(converged=False, error="seed construction failed (see construction_note)")
            results.append(result)
            continue
        state0_guess, period_guess = state0_period
        try:
            if period_guess <= 1e-6:
                raise ValueError("non-positive period guess")
            orbit = correct_periodic(system, state0_guess, period_guess, tol=1e-10, max_iter=30)
        except Exception as exc:  # a bad guess must never crash the sweep
            result.update(converged=False, error=str(exc))
            results.append(result)
            continue
        result.update(
            converged=bool(orbit.converged),
            residual=orbit.closure_residual,
            period=orbit.period,
            jacobi=orbit.jacobi,
            state0=orbit.state0.tolist(),
        )
        results.append(result)
    return results


def _physically_sane_rate(results: list[dict]) -> tuple[float, int]:
    """Fraction of ALL rows (including construction/convergence failures in
    the denominator) that both converged AND landed on a physically-sane
    orbit -- SAME yardstick (`is_physically_sane`, default bounds, #642's
    equilibrium filter already applied by default) as #608/#624/#642.
    """
    n_sane = 0
    for r in results:
        if not r.get("converged"):
            continue
        if r.get("state0") is None or r.get("period") is None or r.get("jacobi") is None:
            continue
        if is_physically_sane(r["state0"], r["period"], r["jacobi"]):
            n_sane += 1
    return (n_sane / len(results) if results else 0.0), n_sane


def _evaluate_mu(
    *,
    label: str,
    system: cr3bp.CR3BPSystem,
    model,
    train: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    """Transformed-generated vs. fresh-in-run uniform-baseline, both at
    `system`'s mu, mirroring #624's own ``_evaluate_mu`` exactly except the
    generated arm is coordinate-transformed before refinement.
    """
    print(f"\n=== mu target: {label} (mu={system.mu:.10e}) ===")

    raw_generated = model.sample(N_GENERATE, rng=rng)
    n_construction_failed = 0
    transformed_rows: list[dict] = []
    for raw_row in raw_generated:
        state0_guess = raw_row[:6]
        period_guess = float(raw_row[6])
        transformed = transform_seed_to_target_mu(
            state0_guess, period_guess, TRAINING_MU, system.mu
        )
        if transformed is None:
            n_construction_failed += 1
            transformed_rows.append(
                {
                    "state0_period": None,
                    "construction_note": "transform_seed_to_target_mu returned None "
                    "(unrealizable velocity at rho-matched target Jacobi constant, "
                    "or zero-velocity guess)",
                }
            )
        else:
            transformed_rows.append(
                {
                    "state0_period": (transformed.state0, transformed.period),
                    "rho": transformed.rho,
                    "scale": transformed.scale,
                    "c_guess_earth_moon": transformed.c_guess,
                    "c_target": transformed.c_target,
                }
            )
    gen_results = _refine_batch(system, transformed_rows, tag=f"{label}_transformed_generated")
    gen_rate = sum(r["converged"] for r in gen_results) / N_GENERATE
    gen_sane_rate, gen_n_sane = _physically_sane_rate(gen_results)
    print(
        f"  [transformed-generated] n={N_GENERATE} "
        f"({n_construction_failed} construction failures): "
        f"raw convergence {gen_rate:.1%}; "
        f"converged+physically-sane {gen_sane_rate:.1%} ({gen_n_sane}/{N_GENERATE})"
    )

    baseline_raw = uniform_bounding_box_sample(train, N_GENERATE, rng=rng)
    baseline_rows = [{"state0_period": (row[:6], float(row[6]))} for row in baseline_raw]
    base_results = _refine_batch(system, baseline_rows, tag=f"{label}_baseline")
    base_rate = sum(r["converged"] for r in base_results) / N_GENERATE
    base_sane_rate, base_n_sane = _physically_sane_rate(base_results)
    print(
        f"  [baseline]  n={N_GENERATE}: raw convergence {base_rate:.1%}; "
        f"converged+physically-sane {base_sane_rate:.1%} ({base_n_sane}/{N_GENERATE})"
    )

    if base_sane_rate > 0:
        ratio: float | None = gen_sane_rate / base_sane_rate
        print(f"  ratio transformed-generated/baseline (physically-sane): {ratio:.2f}x")
    elif gen_sane_rate > 0:
        ratio = math.inf
        print("  baseline: 0 physically-sane; transformed-generated: >0 (infinite ratio)")
    else:
        ratio = None
        print(
            "  BOTH transformed-generated and baseline: 0 physically-sane -- no signal either way"
        )

    ctx = _642_CONTEXT[label]
    print(
        f"  [#642 context, NOT re-measured here] raw (untransformed) generated arm: "
        f"~{ctx['raw_generated_sane_rate']:.1%}; baseline archival "
        f"{ctx['archival_baseline_sane_rate']:.1%}, live re-run "
        f"{ctx['live_rerun_baseline_sane_rate']:.1%}"
    )

    return {
        "label": label,
        "mu": system.mu,
        "transformed_generated": {
            "n": N_GENERATE,
            "n_construction_failed": n_construction_failed,
            "raw_convergence_rate": gen_rate,
            "converged_and_physically_sane_rate": gen_sane_rate,
            "n_converged_and_physically_sane": gen_n_sane,
        },
        "baseline_uniform_fresh_in_run": {
            "n": N_GENERATE,
            "raw_convergence_rate": base_rate,
            "converged_and_physically_sane_rate": base_sane_rate,
            "n_converged_and_physically_sane": base_n_sane,
        },
        "ratio_transformed_generated_over_baseline": ratio,
        "context_642": ctx,
        "_raw_results": gen_results + base_results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--draw-seed",
        type=int,
        default=None,
        help=(
            "Optional independent rng seed for the GENERATION draw only (model "
            "training/split stays pinned to SEED=608 -- bit-for-bit identical model "
            "either way). Default: reuse SEED=608 for the draw too (this script's "
            "original single run). Use a different value to check the headline "
            "result isn't a lucky single draw -- writes to a suffixed summary/detail "
            "file so the original run's artifacts are never overwritten."
        ),
    )
    args = parser.parse_args()

    t_start = time.time()
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] #649 coordinate-fix pilot starting")

    paths = _corpus_paths()
    corpus = assemble_corpus(paths)  # Earth-Moon, SAME filters #608/#624 used
    print(
        f"  scanned {corpus.n_scanned} raw lines -> {corpus.n_converged_prefilter} "
        f"converged -> {len(corpus)} unique, physically-sane "
        f"({corpus.primary}-{corpus.secondary}, mu={corpus.mu:.10f})"
    )

    rng = np.random.default_rng(SEED)
    n = len(corpus)
    idx = rng.permutation(n)
    n_train = int(n * 0.8)
    train_idx = idx[:n_train]
    train = corpus.features[train_idx]
    print(f"  split: {len(train)} train (SAME seed/split as #608/#624 -> identical train set)")

    model = fit_clustered_gaussian(train, n_latent=5, n_clusters=8, seed=SEED)
    print(
        "  re-derived #608's PCA(5)+kmeans(8)-Gaussian model (SAME seed -> bit-for-bit identical)"
    )

    if args.draw_seed is not None:
        rng = np.random.default_rng(args.draw_seed)
        print(f"  using INDEPENDENT draw rng (--draw-seed={args.draw_seed}), model unchanged")

    mu001_system = cr3bp.CR3BPSystem(mu=0.001, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    sun_earth_system = cr3bp_system("Sun", "Earth")

    result_mu001 = _evaluate_mu(
        label="mu=0.001", system=mu001_system, model=model, train=train, rng=rng
    )
    result_sun_earth = _evaluate_mu(
        label="Sun-Earth", system=sun_earth_system, model=model, train=train, rng=rng
    )

    elapsed = time.time() - t_start
    summary = {
        "task": "#649",
        "corpus": {
            "primary": corpus.primary,
            "secondary": corpus.secondary,
            "mu": corpus.mu,
            "n_scanned_raw_lines": corpus.n_scanned,
            "n_converged_prefilter": corpus.n_converged_prefilter,
            "n_unique_physically_sane": len(corpus),
            "n_train": len(train),
            "feature_names": list(FEATURE_NAMES),
        },
        "model": "PCA(5) + kmeans(8) per-cluster Gaussian, re-derived from #608's own seed/split "
        "(bit-for-bit identical to #608's trained model, NOT retrained)",
        "transform": "cyclerfinder.ml.cross_mu_coordinate_transform.transform_seed_to_target_mu "
        "(rho=(C-3)/(C_L1(mu)-3) energy-matching + exact Hill-radius position/velocity-"
        "direction scaling + velocity-magnitude solve to hit target rho exactly; period "
        "left unchanged)",
        "mu_0_001": {k: v for k, v in result_mu001.items() if not k.startswith("_")},
        "sun_earth": {k: v for k, v in result_sun_earth.items() if not k.startswith("_")},
        "elapsed_s": elapsed,
        "seed": SEED,
        "draw_seed": args.draw_seed,
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_drawseed{args.draw_seed}" if args.draw_seed is not None else ""
    summary_path = _OUT_DIR / f"summary{suffix}.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    detail_path = _OUT_DIR / f"refine_results{suffix}.jsonl"
    with detail_path.open("w") as fh:
        for r in result_mu001["_raw_results"] + result_sun_earth["_raw_results"]:
            fh.write(json.dumps(r, default=str) + "\n")

    print(f"\n[{time.strftime('%Y-%m-%dT%H:%M:%S')}] done in {elapsed:.1f}s")
    print(f"  wrote {summary_path}")
    print(f"  wrote {detail_path}")

    print("\nHONEST HEADLINE (converged+physically-sane rate, the decisive metric):")
    for label, res in (("mu=0.001", result_mu001), ("Sun-Earth", result_sun_earth)):
        g = res["transformed_generated"]["converged_and_physically_sane_rate"]
        b = res["baseline_uniform_fresh_in_run"]["converged_and_physically_sane_rate"]
        ratio = res["ratio_transformed_generated_over_baseline"]
        ratio_str = (
            f"{ratio:.2f}x" if isinstance(ratio, float) and math.isfinite(ratio) else str(ratio)
        )
        print(
            f"  {label:10s}: transformed-generated {g:.1%} vs fresh baseline {b:.1%}  "
            f"(ratio: {ratio_str})"
        )
    print(
        "\n  #642's own recorded numbers, for context only (NOT re-measured this run):\n"
        "    Earth-Moon in-distribution (#608/#642 live re-run): generated 27% vs baseline 2% "
        "(~13.5x)\n"
        "    mu=0.001 RAW (untransformed) generated: 0% in two independent trials\n"
        "    Sun-Earth RAW (untransformed) generated: 0% in two independent trials"
    )


if __name__ == "__main__":
    main()
