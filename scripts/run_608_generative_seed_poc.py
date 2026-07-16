#!/usr/bin/env python3
"""#608 proof-of-concept: statistical generative seed model for CR3BP periodic orbits.

Bounded feasibility test of the idea externally de-risked by Litteri, Gil,
Vasile, Rodriguez-Fernandez & Camacho, "Generation of periodic orbits in the
restricted three-body problem with a variational autoencoder," *Celestial
Mechanics and Dynamical Astronomy* 138:25 (June 2026) -- see
``cyclerfinder.ml.orbit_generative``'s module docstring for the architecture
comparison and the deliberate choice of a PCA + k-means-Gaussian model over a
new deep-learning dependency.

This script:

1. Assembles a REAL training corpus from this project's OWN accumulated #210
   outcome-log corpus (``out/outcome_log/*.jsonl``, ~540k raw lines from a
   dozen past Earth-Moon CR3BP corrector campaigns), applying physical-
   plausibility filters (see ``assemble_corpus``'s docstring).
2. Splits it 80/20 train/held-out.
3. Fits the clustered-Gaussian latent model on the train split.
4. SANITY CHECK: projects held-out real orbits through the model's
   encoder/decoder (PCA round-trip) and re-runs the existing
   ``cr3bp_periodic.correct_periodic`` corrector on the RECONSTRUCTED
   (state0, period) guess -- does the model's own compression still land in
   the right basin for orbits it was not shown?
5. GENERATE: samples brand-new latent points, decodes them to (state0,
   period) guesses, and refines them with the SAME existing corrector.
6. BASELINE: draws the same number of guesses uniformly from the training
   data's bounding box (no learned structure at all) and refines those too.
7. Reports honestly: convergence rate for reconstruction / generated /
   baseline, and how many CONVERGED generated candidates are "meaningfully
   new" (nearest-neighbor distance to the full assembled corpus, in
   standardized units) rather than trivial near-duplicates of a training
   point.

No catalogue writeback, no literature-novelty check -- this is a capability
proof-of-concept, not a discovery result (per #608's explicit scope).
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
    fit_clustered_gaussian,
    is_physically_sane,
    nearest_neighbor_distances,
    uniform_bounding_box_sample,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUTCOME_LOG_DIR = _REPO_ROOT / "out" / "outcome_log"
_OUT_DIR = _REPO_ROOT / "data" / "found" / "608_generative_seed_poc"

N_GENERATE = 100
N_HELDOUT_CHECK = 100
SEED = 608


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
    """Run each row's (state0, period) through the existing corrector.

    Returns (per-row result dicts, boolean converged mask).
    """
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
    """Fraction of ALL rows (not just converged ones) that both converged AND
    landed on a physically-sane orbit (see ``is_physically_sane``).

    The solver's raw ``converged`` flag alone is NOT a sufficient success
    metric here -- see this script's honest-headline discussion. A guess can
    satisfy the Newton periodicity tolerance while landing on a degenerate,
    non-physical solution (the exact contamination pattern found in the raw
    training corpus before ``assemble_corpus``'s filters were applied).
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


def main() -> None:
    t_start = time.time()
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] #608 generative-seed POC starting")

    paths = _corpus_paths()
    print(f"  scanning {len(paths)} outcome-log files under {_OUTCOME_LOG_DIR}")
    corpus = assemble_corpus(paths)
    print(
        f"  scanned {corpus.n_scanned} raw lines -> {corpus.n_converged_prefilter} "
        f"converged (any physical bounds) -> {len(corpus)} unique, physically-sane "
        f"({corpus.primary}-{corpus.secondary}, mu={corpus.mu:.10f})"
    )

    rng = np.random.default_rng(SEED)
    n = len(corpus)
    idx = rng.permutation(n)
    n_train = int(n * 0.8)
    train_idx, heldout_idx = idx[:n_train], idx[n_train:]
    train = corpus.features[train_idx]
    heldout = corpus.features[heldout_idx]
    print(f"  split: {len(train)} train / {len(heldout)} held-out")

    model = fit_clustered_gaussian(train, n_latent=5, n_clusters=8, seed=SEED)
    print(
        f"  fitted PCA(5)+kmeans(8)-Gaussian model; explained variance ratio "
        f"(top 5 PCs): {np.round(model.explained_variance_ratio, 4).tolist()}"
    )

    system = cr3bp_system(corpus.primary, corpus.secondary)

    # --- (a) Reconstruction sanity check on held-out real orbits ---
    n_check = min(N_HELDOUT_CHECK, len(heldout))
    check_rows = heldout[rng.choice(len(heldout), size=n_check, replace=False)]
    latent_check = model.transform(check_rows)
    reconstructed = model.inverse_transform(latent_check)
    # The decoded 8th column (jacobi) is a derived quantity, not independent;
    # recompute it from the decoded state so the corrector sees a physically
    # self-consistent guess rather than a PCA-smoothed-but-inconsistent one.
    recon_results, recon_converged = _refine_batch(system, reconstructed, tag="reconstruction")
    recon_rate = float(recon_converged.mean())
    recon_sane_rate, recon_n_sane = _physically_sane_rate(recon_results)
    # distance from the ORIGINAL held-out state to its reconstruction, in raw units.
    recon_state_err = np.linalg.norm((reconstructed[:, :6] - check_rows[:, :6]), axis=1)
    print(
        f"  [reconstruction check] {n_check} held-out orbits: raw convergence rate "
        f"{recon_rate:.1%}; converged+physically-sane rate {recon_sane_rate:.1%} "
        f"({recon_n_sane}/{n_check}); median ||state0_recon - state0_orig|| = "
        f"{np.median(recon_state_err):.4e}"
    )

    # --- (b) Generate NEW candidates from the fitted model ---
    generated = model.sample(N_GENERATE, rng=rng)
    gen_results, gen_converged = _refine_batch(system, generated, tag="generated")
    gen_rate = float(gen_converged.mean())
    gen_sane_rate, gen_n_sane = _physically_sane_rate(gen_results)
    print(
        f"  [generated] {N_GENERATE} latent samples -> refine: raw convergence rate "
        f"{gen_rate:.1%}; converged+physically-sane rate {gen_sane_rate:.1%} "
        f"({gen_n_sane}/{N_GENERATE})"
    )

    # --- (c) Naive baseline: uniform within the SAME training bounding box ---
    baseline = uniform_bounding_box_sample(train, N_GENERATE, rng=rng)
    base_results, base_converged = _refine_batch(system, baseline, tag="baseline_uniform")
    base_rate = float(base_converged.mean())
    base_sane_rate, base_n_sane = _physically_sane_rate(base_results)
    print(
        f"  [baseline]  {N_GENERATE} uniform-box samples -> refine: raw convergence rate "
        f"{base_rate:.1%}; converged+physically-sane rate {base_sane_rate:.1%} "
        f"({base_n_sane}/{N_GENERATE})"
    )

    # --- Novelty check on CONVERGED generated candidates ---
    # Standardize using the model's own fitted mean/scale so distances are
    # comparable across features of very different natural units.
    def _standardize(rows: np.ndarray) -> np.ndarray:
        return (rows - model.mean) / model.scale

    reference_std = _standardize(corpus.features)  # whole assembled corpus, not just train
    gen_converged_rows = np.array(
        [
            r["state0"] + [r["period"], r["jacobi"]]
            for r, ok in zip(gen_results, gen_converged, strict=True)
            if ok
        ]
    )
    if len(gen_converged_rows) > 0:
        gen_nn_dist = nearest_neighbor_distances(_standardize(gen_converged_rows), reference_std)
        # Reference scale: how far apart two REAL, independent orbits typically
        # sit (held-out orbit vs. its nearest TRAIN neighbor) -- an
        # apples-to-apples yardstick for "how novel" the generated candidates are.
        train_std = _standardize(train)
        heldout_std = _standardize(heldout)
        heldout_vs_train_nn = nearest_neighbor_distances(heldout_std, train_std)
        print(
            f"  [novelty] converged generated candidates: {len(gen_converged_rows)}/{N_GENERATE}; "
            f"median NN-distance-to-corpus = {np.median(gen_nn_dist):.4f} (standardized units); "
            f"reference: held-out-real-orbit-vs-train median NN distance = "
            f"{np.median(heldout_vs_train_nn):.4f} (i.e. how far apart two REAL, "
            f"independent orbits typically sit)"
        )
    else:
        gen_nn_dist = np.array([])
        heldout_vs_train_nn = np.array([])
        print("  [novelty] no generated candidates converged -- novelty check skipped")

    elapsed = time.time() - t_start
    summary = {
        "task": "#608",
        "corpus": {
            "primary": corpus.primary,
            "secondary": corpus.secondary,
            "mu": corpus.mu,
            "n_scanned_raw_lines": corpus.n_scanned,
            "n_converged_prefilter": corpus.n_converged_prefilter,
            "n_unique_physically_sane": len(corpus),
            "n_train": len(train),
            "n_heldout": len(heldout),
            "feature_names": list(FEATURE_NAMES),
        },
        "model": {
            "kind": "PCA(5) + kmeans(8) per-cluster Gaussian (linear-Gaussian VAE analog)",
            "explained_variance_ratio_top5": model.explained_variance_ratio.tolist(),
        },
        "reconstruction_check": {
            "n": n_check,
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
            "median_nn_distance_to_corpus": (
                float(np.median(gen_nn_dist)) if len(gen_nn_dist) else None
            ),
        },
        "baseline_uniform": {
            "n": N_GENERATE,
            "raw_convergence_rate": base_rate,
            "n_converged": int(base_converged.sum()),
            "converged_and_physically_sane_rate": base_sane_rate,
            "n_converged_and_physically_sane": base_n_sane,
        },
        "reference_real_orbit_nn_spacing": {
            "median_heldout_vs_train": float(np.median(heldout_vs_train_nn))
            if len(heldout_vs_train_nn)
            else None,
        },
        "elapsed_s": elapsed,
        "seed": SEED,
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    detail_path = _OUT_DIR / "refine_results.jsonl"
    with detail_path.open("w") as fh:
        for r in recon_results + gen_results + base_results:
            fh.write(json.dumps(r) + "\n")

    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] done in {elapsed:.1f}s")
    print(f"  wrote {summary_path}")
    print(f"  wrote {detail_path}")
    print()
    print("HONEST HEADLINE:")
    print("  Raw solver 'converged' rate (a Newton periodicity-tolerance flag only,")
    print("  NOT a physical-plausibility check -- see this script's docstring):")
    print(f"    reconstruction (sanity check): {recon_rate:.1%}")
    print(f"    generated (new samples):        {gen_rate:.1%}")
    print(f"    baseline (uniform, no model):   {base_rate:.1%}")
    print("  These raw rates are near-identical and are NOT the right comparison --")
    print("  a large fraction of 'converged' baseline solves land on degenerate,")
    print("  non-physical solutions (the same contamination the training corpus")
    print("  itself needed filtering for). The metric that matters is CONVERGED")
    print("  *AND* physically-sane (same bounds used to build the training corpus):")
    print(f"    reconstruction (sanity check): {recon_sane_rate:.1%} ({recon_n_sane}/{n_check})")
    print(f"    generated (new samples):        {gen_sane_rate:.1%} ({gen_n_sane}/{N_GENERATE})")
    print(f"    baseline (uniform, no model):   {base_sane_rate:.1%} ({base_n_sane}/{N_GENERATE})")
    if base_sane_rate > 0:
        ratio = gen_sane_rate / base_sane_rate
        print(f"    ratio generated/baseline (physically-sane): {ratio:.2f}x")
    elif gen_sane_rate > 0:
        print("    baseline: 0 physically-sane; generated: >0 physically-sane (infinite ratio)")
    else:
        print("    BOTH generated and baseline: 0 physically-sane -- no signal either way")


if __name__ == "__main__":
    main()
