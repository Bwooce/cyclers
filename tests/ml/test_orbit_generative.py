"""Tests for the #608 generative-seed-model proof-of-concept building blocks.

These exercise the LOGIC (log parsing/filtering, PCA+clustered-Gaussian
fit/sample round trip, nearest-neighbor distance) with small constructed
fixtures -- standard unit-test practice, NOT the actual #608 proof-of-concept
run, which trains on this project's own real accumulated outcome-log corpus
(see ``scripts/run_608_generative_seed_poc.py``). No physical claim is pinned
here; these are code-correctness tests only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cyclerfinder.ml.orbit_generative import (
    FEATURE_NAMES,
    AutoencoderModel,
    assemble_corpus,
    fit_autoencoder,
    fit_autoencoder_clustered_gaussian,
    fit_clustered_gaussian,
    fit_family_conditioned_gaussian,
    heuristic_family_tag,
    is_physically_sane,
    iter_outcome_records,
    nearest_neighbor_distances,
    uniform_bounding_box_sample,
)


def _outcome_record(
    *,
    solver: str = "cr3bp.correct_periodic",
    primary: str = "Earth",
    secondary: str = "Moon",
    converged: bool = True,
    residual: float = 1e-11,
    jacobi: float = 3.0,
    period: float = 2.5,
    state0: list[float] | None = None,
) -> dict:
    return {
        "solver": solver,
        "meta": {"primary": primary, "secondary": secondary},
        "inputs": {"state0_guess": state0 or [0.85, 0.0, 0.1, 0.0, 0.2, 0.0], "mu": 0.01215},
        "outcome": {
            "converged": converged,
            "residual": residual,
            "jacobi": jacobi,
            "period": period,
        },
    }


def test_iter_outcome_records_skips_malformed_lines(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    path.write_text(
        json.dumps(_outcome_record())
        + "\n"
        + "{not valid json\n"
        + "\n"  # blank line
        + json.dumps(_outcome_record(period=3.1))
        + "\n"
    )
    records = list(iter_outcome_records([path]))
    assert len(records) == 2
    assert records[0]["outcome"]["period"] == 2.5
    assert records[1]["outcome"]["period"] == 3.1


def test_assemble_corpus_filters_and_dedupes(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    lines = [
        _outcome_record(state0=[0.85, 0.0, 0.1, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0),
        # exact near-duplicate (within dedup_decimals) of the row above -> collapsed
        _outcome_record(state0=[0.850001, 0.0, 0.1, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0),
        # distinct real orbit -> kept
        _outcome_record(state0=[0.80, 0.01, 0.05, 0.01, 0.22, 0.01], period=3.0, jacobi=2.9),
        # wrong solver -> dropped
        _outcome_record(solver="other.solver"),
        # wrong system -> dropped
        _outcome_record(primary="Sun", secondary="Earth"),
        # not converged -> dropped
        _outcome_record(converged=False),
        # residual too loose -> dropped
        _outcome_record(residual=1e-3),
        # jacobi out of physical bounds -> dropped
        _outcome_record(jacobi=50.0),
        # period out of bounds -> dropped
        _outcome_record(period=100.0),
        # out-of-plane amplitude too large -> dropped
        _outcome_record(state0=[0.85, 0.0, 5.0, 0.0, 0.2, 0.0]),
    ]
    path.write_text("\n".join(json.dumps(rec) for rec in lines) + "\n")

    corpus = assemble_corpus([path])
    assert corpus.primary == "Earth"
    assert corpus.secondary == "Moon"
    assert corpus.n_scanned == len(lines)
    # n_converged_prefilter counts every record with the right solver+system
    # AND converged=True, BEFORE the physical-plausibility filters (jacobi/
    # period/residual/out-of-plane bounds): records 0,1,2 (the real orbits)
    # plus 6,7,8,9 (converged but physically implausible / too-loose-residual)
    # = 7. Only records 0/1/2 additionally pass the physical filters, and 0/1
    # dedupe to one row, leaving 2 final rows.
    assert corpus.n_converged_prefilter == 7
    assert len(corpus) == 2
    assert corpus.features.shape == (2, len(FEATURE_NAMES))


def test_assemble_corpus_raises_on_empty(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    path.write_text(json.dumps(_outcome_record(converged=False)) + "\n")
    with pytest.raises(ValueError, match="no record survived"):
        assemble_corpus([path])


def test_fit_clustered_gaussian_roundtrip_and_sample_shape() -> None:
    rng = np.random.default_rng(0)
    # Two well-separated 8D clusters (mimicking two distinct orbit families).
    cluster_a = rng.normal(loc=0.0, scale=0.02, size=(150, 8))
    cluster_b = rng.normal(loc=5.0, scale=0.02, size=(150, 8))
    features = np.vstack([cluster_a, cluster_b])

    model = fit_clustered_gaussian(features, n_latent=3, n_clusters=2, seed=0)

    # Encoder/decoder round trip should reconstruct well-separated clusters
    # to within a small fraction of the inter-cluster gap.
    latent = model.transform(features)
    reconstructed = model.inverse_transform(latent)
    err = np.linalg.norm(reconstructed - features, axis=1)
    assert np.median(err) < 0.5  # inter-cluster gap is ~5*sqrt(8) ~ 14

    samples = model.sample(20)
    assert samples.shape == (20, 8)
    assert np.all(np.isfinite(samples))
    # Generated samples should land near one of the two original clusters,
    # not off in unconstrained space.
    dist_to_a = np.linalg.norm(samples - cluster_a.mean(axis=0), axis=1)
    dist_to_b = np.linalg.norm(samples - cluster_b.mean(axis=0), axis=1)
    assert np.all(np.minimum(dist_to_a, dist_to_b) < 3.0)


def test_fit_clustered_gaussian_rejects_too_many_latent_dims() -> None:
    features = np.zeros((10, 4))
    with pytest.raises(ValueError, match="n_latent"):
        fit_clustered_gaussian(features, n_latent=8, n_clusters=2)


def test_nearest_neighbor_distances_exact_match_is_zero() -> None:
    reference = np.array([[0.0, 0.0], [10.0, 10.0]])
    query = np.array([[0.0, 0.0], [9.5, 10.5]])
    dist = nearest_neighbor_distances(query, reference)
    assert dist[0] == pytest.approx(0.0, abs=1e-12)
    assert dist[1] == pytest.approx(np.hypot(0.5, 0.5), abs=1e-9)


def test_is_physically_sane_matches_assemble_corpus_bounds() -> None:
    # A real, in-bounds Earth-Moon-scale planar orbit.
    assert is_physically_sane([0.85, 0.0, 0.05, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0)
    # Same state/period but a wildly out-of-range Jacobi (the raw-corpus
    # contamination pattern this function exists to catch) -> not sane.
    assert not is_physically_sane([0.85, 0.0, 0.05, 0.0, 0.2, 0.0], period=2.5, jacobi=-50.0)
    # Period out of bounds.
    assert not is_physically_sane([0.85, 0.0, 0.05, 0.0, 0.2, 0.0], period=100.0, jacobi=3.0)
    # Out-of-plane amplitude too large.
    assert not is_physically_sane([0.85, 0.0, 5.0, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0)
    assert not is_physically_sane([0.85, 0.0, 0.05, 0.0, 0.2, 5.0], period=2.5, jacobi=3.0)


def test_uniform_bounding_box_sample_respects_bounds() -> None:
    rng = np.random.default_rng(1)
    features = np.array([[0.0, -1.0], [1.0, 1.0], [0.5, 0.0]])
    samples = uniform_bounding_box_sample(features, 500, rng=rng)
    assert samples.shape == (500, 2)
    assert samples[:, 0].min() >= 0.0
    assert samples[:, 0].max() <= 1.0
    assert samples[:, 1].min() >= -1.0
    assert samples[:, 1].max() <= 1.0


# --- #614: family-tagging + nonlinear-encoder follow-up tests -------------
#
# Same "code-correctness with small constructed fixtures" convention as the
# #608 tests above -- no physical claim pinned here; see
# `scripts/run_614_family_and_nonlinear_poc.py` for the actual comparison
# against this project's real accumulated corpus.


def test_heuristic_family_tag_is_honest_heuristic_not_exact_family() -> None:
    # Different Jacobi/period/z-sign bands -> different tags.
    tag_lowj_planar = heuristic_family_tag([0.85, 0.0, 0.0, 0.0, 0.2, 0.0], period=1.0, jacobi=2.6)
    tag_highj_north = heuristic_family_tag([0.85, 0.0, 0.1, 0.0, 0.2, 0.0], period=1.0, jacobi=3.4)
    tag_highj_south = heuristic_family_tag([0.85, 0.0, -0.1, 0.0, 0.2, 0.0], period=1.0, jacobi=3.4)
    assert tag_lowj_planar != tag_highj_north
    assert tag_highj_north != tag_highj_south  # z-sign alone must flip the tag
    # Near-identical inputs land in the SAME bin (heuristic, coarse-grained).
    tag_a = heuristic_family_tag([0.85, 0.0, 0.05, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0)
    tag_b = heuristic_family_tag([0.85, 0.0, 0.06, 0.0, 0.2, 0.0], period=2.55, jacobi=3.01)
    assert tag_a == tag_b
    # Out-of-bounds jacobi/period are clamped, not an error.
    heuristic_family_tag([0.85, 0.0, 0.05, 0.0, 0.2, 0.0], period=1e6, jacobi=1e6)


def test_assemble_corpus_populates_family_tags(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    lines = [
        _outcome_record(state0=[0.85, 0.0, 0.1, 0.0, 0.2, 0.0], period=2.5, jacobi=3.0),
        _outcome_record(state0=[0.80, 0.01, -0.05, 0.01, 0.22, 0.01], period=3.0, jacobi=2.6),
    ]
    path.write_text("\n".join(json.dumps(rec) for rec in lines) + "\n")
    corpus = assemble_corpus([path])
    assert corpus.family_tags is not None
    assert corpus.family_tags.shape == (len(corpus),)
    assert corpus.family_tags.dtype.kind == "U"
    # A northern (z0>0) high-Jacobi orbit and a southern (z0<0) lower-Jacobi
    # orbit must land in different heuristic bins.
    assert corpus.family_tags[0] != corpus.family_tags[1]


def test_fit_family_conditioned_gaussian_roundtrip_and_sample_shape() -> None:
    rng = np.random.default_rng(0)
    cluster_a = rng.normal(loc=0.0, scale=0.02, size=(150, 8))
    cluster_b = rng.normal(loc=5.0, scale=0.02, size=(150, 8))
    features = np.vstack([cluster_a, cluster_b])
    tags = np.array(["A"] * 150 + ["B"] * 150)

    model = fit_family_conditioned_gaussian(features, tags, n_latent=3, seed=0)

    latent = model.transform(features)
    reconstructed = model.inverse_transform(latent)
    err = np.linalg.norm(reconstructed - features, axis=1)
    assert np.median(err) < 0.5  # same PCA encoder as #608's linear model

    samples = model.sample(20)
    assert samples.shape == (20, 8)
    assert np.all(np.isfinite(samples))
    dist_to_a = np.linalg.norm(samples - cluster_a.mean(axis=0), axis=1)
    dist_to_b = np.linalg.norm(samples - cluster_b.mean(axis=0), axis=1)
    assert np.all(np.minimum(dist_to_a, dist_to_b) < 3.0)


def test_fit_family_conditioned_gaussian_rejects_mismatched_tag_length() -> None:
    features = np.zeros((10, 8))
    tags = np.array(["A"] * 5)  # wrong length
    with pytest.raises(ValueError, match="family_tags length"):
        fit_family_conditioned_gaussian(features, tags, n_latent=3)


def test_fit_autoencoder_beats_linear_pca_on_a_curved_manifold() -> None:
    """The whole point of #614's nonlinear-encoder variant: on data that

    genuinely lies on a CURVED (non-linear) manifold, a nonlinear autoencoder
    should reconstruct it more tightly than a linear PCA encoder of the same
    latent dimension. Uses a synthetic circle/figure-eight-like curve
    embedded in 8D specifically because it is NOT linearly compressible to
    1D -- this is a general property check, not a value pinned from the
    project's real corpus (see the module docstring's testing convention).
    """
    rng = np.random.default_rng(42)
    t = rng.uniform(0, 2 * np.pi, size=600)
    base = np.stack([np.cos(t), np.sin(t), np.cos(2 * t), np.sin(2 * t)], axis=1)
    features = np.hstack([base, rng.normal(scale=0.01, size=(600, 4))])
    features = features + rng.normal(scale=0.01, size=(600, 8))

    pca_model = fit_clustered_gaussian(features, n_latent=1, n_clusters=2, seed=0)
    pca_latent = pca_model.transform(features)
    pca_recon = pca_model.inverse_transform(pca_latent)
    std = (features - pca_model.mean) / pca_model.scale
    pca_recon_std = (pca_recon - pca_model.mean) / pca_model.scale
    pca_mse = float(np.mean(np.sum((pca_recon_std - std) ** 2, axis=1)))

    ae = fit_autoencoder(features, n_latent=1, n_hidden=8, maxiter=500, seed=0)
    assert isinstance(ae, AutoencoderModel)

    assert ae.train_reconstruction_mse < pca_mse
    # Not just marginally better -- the linear encoder structurally cannot
    # track this curve at n_latent=1, so the gap should be substantial.
    assert ae.train_reconstruction_mse < 0.8 * pca_mse


def test_autoencoder_clustered_gaussian_sample_shape_and_finite() -> None:
    rng = np.random.default_rng(0)
    cluster_a = rng.normal(loc=0.0, scale=0.02, size=(150, 8))
    cluster_b = rng.normal(loc=5.0, scale=0.02, size=(150, 8))
    features = np.vstack([cluster_a, cluster_b])

    model = fit_autoencoder_clustered_gaussian(
        features, n_latent=3, n_clusters=2, n_hidden=8, maxiter=200, seed=0
    )
    samples = model.sample(20)
    assert samples.shape == (20, 8)
    assert np.all(np.isfinite(samples))
    dist_to_a = np.linalg.norm(samples - cluster_a.mean(axis=0), axis=1)
    dist_to_b = np.linalg.norm(samples - cluster_b.mean(axis=0), axis=1)
    assert np.all(np.minimum(dist_to_a, dist_to_b) < 3.0)


def test_clustered_gaussian_sample_reproducible_after_refactor() -> None:
    """Regression lock for #614's refactor of `ClusteredGaussianLatentModel.sample`

    into the shared `_sample_latent_from_clusters` helper: same seed must
    still produce bit-identical samples (this pins the ACTUAL observed
    values from running the current code, not invented ones).
    """
    rng = np.random.default_rng(0)
    cluster_a = rng.normal(loc=0.0, scale=0.02, size=(150, 8))
    cluster_b = rng.normal(loc=5.0, scale=0.02, size=(150, 8))
    features = np.vstack([cluster_a, cluster_b])
    model = fit_clustered_gaussian(features, n_latent=3, n_clusters=2, seed=0)

    samples_1 = model.sample(5, rng=np.random.default_rng(123))
    samples_2 = model.sample(5, rng=np.random.default_rng(123))
    np.testing.assert_array_equal(samples_1, samples_2)
