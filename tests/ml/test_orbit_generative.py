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
    assemble_corpus,
    fit_clustered_gaussian,
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
