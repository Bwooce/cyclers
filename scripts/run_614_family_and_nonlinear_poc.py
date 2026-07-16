#!/usr/bin/env python3
"""#614: family-tagging + nonlinear-encoder follow-up to `#608`'s generative

seed-model proof-of-concept.

`#608` built a linear-PCA + k-means-partitioned-Gaussian generative model
(``cyclerfinder.ml.orbit_generative.ClusteredGaussianLatentModel``) trained on
this project's own accumulated ``#210`` outcome-log corpus, and found it lifts
the existing corrector's physically-sane convergence rate ~12x over uniform
seeding (49% vs 4%). Its own honest verdict flagged two concrete gaps this
script tests, using the EXACT SAME corpus, split, and evaluation pipeline
`#608` already built (``assemble_corpus``, ``correct_periodic`` refinement,
``is_physically_sane``, ``nearest_neighbor_distances``) so the comparison is
apples-to-apples -- no new metric is invented here:

1. **Family-tagging** -- the corpus had no family label, only primary/
   secondary. `#614` adds ``heuristic_family_tag`` (Jacobi-band x period-band x
   z0-sign/deadband; see its docstring for the explicit "this is NOT a
   validated family classifier" caveat) and tests whether sampling from a
   per-tag-conditioned model (``fit_family_conditioned_gaussian``) beats the
   unsupervised k-means-partitioned model on the SAME linear PCA encoder.

2. **Nonlinear encoder** -- `#608` found generated candidates land unusually
   far (median NN distance ~40x a real orbit pair's spacing) from any single
   training example, evidence the LINEAR PCA encoder doesn't tightly track
   the true curved family manifolds. `#614` fits a small from-scratch
   numpy/scipy autoencoder (``fit_autoencoder`` / ``fit_autoencoder_
   clustered_gaussian`` -- no new heavy ML dependency, matching `#608`'s own
   established discipline) and reruns the SAME reconstruction / generate-
   then-refine / novelty checks.

Three models are compared against `#608`'s existing uniform baseline:
  (A) linear PCA + k-means(8)         -- `#608`'s baseline, REPRODUCED here
      bit-for-bit (same seed/split) as this run's regression floor.
  (B) linear PCA + family-tag groups  -- tests family-conditioning ALONE.
  (C) nonlinear autoencoder + k-means -- tests the nonlinear encoder ALONE.

No catalogue writeback, no literature-novelty check -- capability comparison,
not a discovery result (per `#614`'s explicit scope).
"""

from __future__ import annotations

import glob
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.ml.orbit_generative import (
    FEATURE_NAMES,
    assemble_corpus,
    fit_autoencoder_clustered_gaussian,
    fit_clustered_gaussian,
    fit_family_conditioned_gaussian,
    is_physically_sane,
    nearest_neighbor_distances,
    uniform_bounding_box_sample,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUTCOME_LOG_DIR = _REPO_ROOT / "out" / "outcome_log"
_OUT_DIR = _REPO_ROOT / "data" / "found" / "614_family_and_nonlinear_poc"

N_GENERATE = 100
N_HELDOUT_CHECK = 100
SEED = 608  # SAME seed as #608 so the split/baseline reproduce bit-for-bit.


def _corpus_paths() -> list[Path]:
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


def _refine_batch(system, rows: np.ndarray, *, tag: str) -> tuple[list[dict], np.ndarray]:
    """Run each row's (state0, period) through the existing corrector (SAME as #608)."""
    results = []
    converged_mask = np.zeros(len(rows), dtype=bool)
    for i, row in enumerate(rows):
        state0_guess = row[:6]
        period_guess = float(row[6])
        try:
            if period_guess <= 1e-6:
                raise ValueError("non-positive period guess")
            orbit = correct_periodic(system, state0_guess, period_guess, tol=1e-10, max_iter=30)
        except Exception as exc:  # a bad decoded guess must never crash the sweep
            results.append({"tag": tag, "index": i, "converged": False, "error": str(exc)})
            continue
        converged_mask[i] = orbit.converged
        results.append(
            {
                "tag": tag,
                "index": i,
                "converged": bool(orbit.converged),
                "residual": orbit.closure_residual,
                "period": orbit.period,
                "jacobi": orbit.jacobi,
                "state0": orbit.state0.tolist(),
            }
        )
    return results, converged_mask


def _physically_sane_rate(results: list[dict]) -> tuple[float, int]:
    """Fraction of ALL rows that both converged AND landed on a physically-sane
    orbit (see #608's script docstring for why the raw 'converged' flag alone
    is not the right metric)."""
    n_sane = 0
    for r in results:
        if not r.get("converged"):
            continue
        if r.get("state0") is None or r.get("period") is None or r.get("jacobi") is None:
            continue
        if is_physically_sane(r["state0"], r["period"], r["jacobi"]):
            n_sane += 1
    return (n_sane / len(results) if results else 0.0), n_sane


def _evaluate_model(
    *,
    name: str,
    model,
    system,
    corpus_features: np.ndarray,
    corpus_mean: np.ndarray,
    corpus_scale: np.ndarray,
    train: np.ndarray,
    heldout: np.ndarray,
    heldout_check_rows: np.ndarray,
    baseline_sane_rate: float,
    rng: np.random.Generator,
) -> dict:
    """Run the SAME #608 reconstruction / generate / novelty checks for one model."""
    print(f"\n  === model: {name} ===")

    # (a) reconstruction sanity check on held-out real orbits.
    latent_check = model.transform(heldout_check_rows)
    reconstructed = model.inverse_transform(latent_check)
    recon_results, recon_converged = _refine_batch(system, reconstructed, tag=f"{name}_recon")
    recon_rate = float(recon_converged.mean())
    recon_sane_rate, recon_n_sane = _physically_sane_rate(recon_results)
    recon_state_err = np.linalg.norm(reconstructed[:, :6] - heldout_check_rows[:, :6], axis=1)
    print(
        f"    [reconstruction] n={len(heldout_check_rows)}: raw convergence "
        f"{recon_rate:.1%}; converged+sane {recon_sane_rate:.1%} "
        f"({recon_n_sane}/{len(heldout_check_rows)}); "
        f"median ||state0_recon - state0_orig|| = {np.median(recon_state_err):.4e}"
    )

    # (b) generate NEW candidates.
    generated = model.sample(N_GENERATE, rng=rng)
    gen_results, gen_converged = _refine_batch(system, generated, tag=f"{name}_generated")
    gen_rate = float(gen_converged.mean())
    gen_sane_rate, gen_n_sane = _physically_sane_rate(gen_results)
    print(
        f"    [generated] n={N_GENERATE}: raw convergence {gen_rate:.1%}; "
        f"converged+sane {gen_sane_rate:.1%} ({gen_n_sane}/{N_GENERATE})"
    )
    if baseline_sane_rate > 0:
        ratio = gen_sane_rate / baseline_sane_rate
        print(f"    ratio vs uniform baseline (physically-sane): {ratio:.2f}x")

    # (c) novelty: NN distance of CONVERGED generated candidates to full corpus.
    def _standardize(rows: np.ndarray) -> np.ndarray:
        return (rows - corpus_mean) / corpus_scale

    reference_std = _standardize(corpus_features)
    gen_converged_rows = np.array(
        [
            r["state0"] + [r["period"], r["jacobi"]]
            for r, ok in zip(gen_results, gen_converged, strict=True)
            if ok
        ]
    )
    if len(gen_converged_rows) > 0:
        gen_nn_dist = nearest_neighbor_distances(_standardize(gen_converged_rows), reference_std)
        median_nn = float(np.median(gen_nn_dist))
        print(
            f"    [novelty] converged: {len(gen_converged_rows)}/{N_GENERATE}; "
            f"median NN-distance = {median_nn:.4f}"
        )
    else:
        median_nn = None
        print("    [novelty] no generated candidates converged -- skipped")

    return {
        "model": name,
        "reconstruction_check": {
            "n": len(heldout_check_rows),
            "raw_convergence_rate": recon_rate,
            "converged_and_physically_sane_rate": recon_sane_rate,
            "n_converged_and_physically_sane": recon_n_sane,
            "median_state0_error": float(np.median(recon_state_err)),
        },
        "generated": {
            "n": N_GENERATE,
            "raw_convergence_rate": gen_rate,
            "n_converged": int(gen_converged.sum()),
            "converged_and_physically_sane_rate": gen_sane_rate,
            "n_converged_and_physically_sane": gen_n_sane,
            "median_nn_distance_to_corpus": median_nn,
        },
    }


def main() -> None:
    t_start = time.time()
    print(
        f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] #614 family-tag + nonlinear-encoder POC starting"
    )

    paths = _corpus_paths()
    corpus = assemble_corpus(paths)
    n_unique_tags = len(set(corpus.family_tags.tolist()))
    print(
        f"  scanned {corpus.n_scanned} raw lines -> {corpus.n_converged_prefilter} "
        f"converged (any physical bounds) -> {len(corpus)} unique, physically-sane "
        f"({corpus.primary}-{corpus.secondary}, mu={corpus.mu:.10f}); "
        f"{n_unique_tags} distinct heuristic family tags"
    )

    rng = np.random.default_rng(SEED)
    n = len(corpus)
    idx = rng.permutation(n)
    n_train = int(n * 0.8)
    train_idx, heldout_idx = idx[:n_train], idx[n_train:]
    train = corpus.features[train_idx]
    heldout = corpus.features[heldout_idx]
    train_tags = corpus.family_tags[train_idx]
    print(f"  split: {len(train)} train / {len(heldout)} held-out (SAME seed/split as #608)")

    n_check = min(N_HELDOUT_CHECK, len(heldout))
    check_rows = heldout[rng.choice(len(heldout), size=n_check, replace=False)]

    system = cr3bp_system(corpus.primary, corpus.secondary)

    # --- Baseline: uniform-box sample, SAME as #608 ---
    baseline = uniform_bounding_box_sample(train, N_GENERATE, rng=rng)
    base_results, base_converged = _refine_batch(system, baseline, tag="baseline_uniform")
    base_rate = float(base_converged.mean())
    base_sane_rate, base_n_sane = _physically_sane_rate(base_results)
    print(
        f"\n  [baseline, uniform] n={N_GENERATE}: raw convergence {base_rate:.1%}; "
        f"converged+sane {base_sane_rate:.1%} ({base_n_sane}/{N_GENERATE})"
    )

    # --- Model A: linear PCA + k-means(8) -- #608's own baseline, reproduced ---
    model_a = fit_clustered_gaussian(train, n_latent=5, n_clusters=8, seed=SEED)

    # --- Model B: linear PCA + family-tag groups ---
    model_b = fit_family_conditioned_gaussian(train, train_tags, n_latent=5, seed=SEED)

    # --- Model C: nonlinear autoencoder + k-means(8) ---
    t_ae = time.time()
    model_c = fit_autoencoder_clustered_gaussian(
        train, n_latent=5, n_clusters=8, n_hidden=16, maxiter=300, seed=SEED
    )
    ae_mse = model_c.encoder.train_reconstruction_mse
    print(
        f"\n  fitted autoencoder in {time.time() - t_ae:.1f}s; "
        f"train recon MSE (standardized) = {ae_mse:.4f}"
    )

    # Linear-PCA reconstruction MSE for the SAME comparison (standardized units).
    lat_a_train = model_a.transform(train)
    recon_a_train = model_a.inverse_transform(lat_a_train)
    std_train = (train - model_a.mean) / model_a.scale
    std_recon_a = (recon_a_train - model_a.mean) / model_a.scale
    pca_train_mse = float(np.mean(np.sum((std_recon_a - std_train) ** 2, axis=1)))
    print(
        f"  encoder reconstruction MSE (standardized, TRAIN set): linear PCA = "
        f"{pca_train_mse:.4f}; autoencoder = {ae_mse:.4f}"
    )

    results = {}
    for name, model in (
        ("A_linear_kmeans", model_a),
        ("B_linear_family_tag", model_b),
        ("C_autoencoder_kmeans", model_c),
    ):
        results[name] = _evaluate_model(
            name=name,
            model=model,
            system=system,
            corpus_features=corpus.features,
            corpus_mean=model_a.mean,
            corpus_scale=model_a.scale,
            train=train,
            heldout=heldout,
            heldout_check_rows=check_rows,
            baseline_sane_rate=base_sane_rate,
            rng=rng,
        )

    # Reference: how far apart two REAL independent orbits typically sit (SAME
    # yardstick #608 used), standardized by model A's (linear PCA) mean/scale.
    def _standardize_a(rows: np.ndarray) -> np.ndarray:
        return (rows - model_a.mean) / model_a.scale

    heldout_vs_train_nn = nearest_neighbor_distances(_standardize_a(heldout), _standardize_a(train))
    median_real_nn = float(np.median(heldout_vs_train_nn))
    print(f"\n  reference: held-out-real-orbit-vs-train median NN distance = {median_real_nn:.4f}")

    elapsed = time.time() - t_start
    summary = {
        "task": "#614",
        "corpus": {
            "primary": corpus.primary,
            "secondary": corpus.secondary,
            "mu": corpus.mu,
            "n_scanned_raw_lines": corpus.n_scanned,
            "n_converged_prefilter": corpus.n_converged_prefilter,
            "n_unique_physically_sane": len(corpus),
            "n_train": len(train),
            "n_heldout": len(heldout),
            "n_distinct_heuristic_family_tags": n_unique_tags,
            "feature_names": list(FEATURE_NAMES),
        },
        "encoder_reconstruction_mse_standardized_train": {
            "linear_pca": pca_train_mse,
            "autoencoder": ae_mse,
        },
        "baseline_uniform": {
            "n": N_GENERATE,
            "raw_convergence_rate": base_rate,
            "n_converged": int(base_converged.sum()),
            "converged_and_physically_sane_rate": base_sane_rate,
            "n_converged_and_physically_sane": base_n_sane,
        },
        "models": results,
        "reference_real_orbit_nn_spacing_linear_pca_units": median_real_nn,
        "elapsed_s": elapsed,
        "seed": SEED,
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    print(f"\n[{time.strftime('%Y-%m-%dT%H:%M:%S')}] done in {elapsed:.1f}s")
    print(f"  wrote {summary_path}")

    print("\nHONEST HEADLINE (converged+physically-sane rate, the metric that matters):")
    print(f"  baseline (uniform, no model):      {base_sane_rate:.1%} ({base_n_sane}/{N_GENERATE})")
    for name in ("A_linear_kmeans", "B_linear_family_tag", "C_autoencoder_kmeans"):
        r = results[name]["generated"]
        ratio = (
            f"{r['converged_and_physically_sane_rate'] / base_sane_rate:.2f}x"
            if base_sane_rate > 0
            else "n/a"
        )
        print(
            f"  {name:22s}: {r['converged_and_physically_sane_rate']:.1%} "
            f"({r['n_converged_and_physically_sane']}/{N_GENERATE}), ratio vs baseline: {ratio}, "
            f"median NN dist: {r['median_nn_distance_to_corpus']}"
        )
    print(
        f"\n  encoder reconstruction MSE (standardized, lower=tighter manifold fit): "
        f"linear PCA={pca_train_mse:.4f} vs autoencoder={ae_mse:.4f}"
    )


if __name__ == "__main__":
    main()
