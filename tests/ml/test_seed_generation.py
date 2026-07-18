"""Tests for the #628 productionized seed-generation API
(``cyclerfinder.ml.seed_generation``).

Follows the SAME "code-correctness with small constructed fixtures" testing
convention as ``tests/ml/test_orbit_generative.py`` -- no physical claim is
pinned via a synthetic fixture. Two tests below (marked
``skipif``) additionally reuse `#608`'s/`#624`'s own real-corpus evaluation
protocol as a REGRESSION CHECK against those tasks' headline lift numbers;
they are skipped automatically in any environment lacking this project's own
gitignored ``out/outcome_log/`` corpus (see
``cyclerfinder.ml.seed_generation.default_corpus_paths``'s docstring).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import CR3BPSystem
from cyclerfinder.ml.orbit_generative import (
    ClusteredGaussianLatentModel,
    fit_clustered_gaussian,
    is_physically_sane,
    uniform_bounding_box_sample,
)
from cyclerfinder.ml.seed_generation import (
    LIFT_ANCHORS,
    TRAINING_MU,
    GeneratedSeed,
    RecalibrationResult,
    SeedGenerationReport,
    calibrate_cluster_weights_for_mu,
    clear_model_cache,
    default_corpus_paths,
    delta_log10_mu,
    expected_lift_for_mu,
    generate_and_refine_seeds,
    get_default_model,
    resolve_system,
    stability_index,
)
from cyclerfinder.search.cr3bp_periodic import correct_periodic

try:
    default_corpus_paths()
    _HAS_REAL_CORPUS = True
except FileNotFoundError:
    _HAS_REAL_CORPUS = False

_skip_no_real_corpus = pytest.mark.skipif(
    not _HAS_REAL_CORPUS,
    reason="requires this project's own gitignored out/outcome_log/ corpus locally",
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
) -> dict[str, Any]:
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


def _synthetic_corpus_file(tmp_path: Path, *, n_per_cluster: int = 30, n_clusters: int = 8) -> Path:
    """Build a synthetic outcome-log file with enough distinct rows to fit
    ``fit_clustered_gaussian(n_latent=5, n_clusters=8)`` -- #608's exact
    hyperparameters -- without depending on the real project corpus.
    """
    rng = np.random.default_rng(42)
    lines: list[dict[str, Any]] = []
    for c in range(n_clusters):
        jacobi_center = 2.6 + 0.1 * c
        period_center = 1.0 + c
        for _ in range(n_per_cluster):
            state0 = [
                0.8 + rng.normal(scale=0.01),
                rng.normal(scale=0.01),
                0.05 * (1 if c % 2 == 0 else -1) + rng.normal(scale=0.005),
                rng.normal(scale=0.01),
                0.2 + rng.normal(scale=0.01),
                rng.normal(scale=0.01),
            ]
            lines.append(
                _outcome_record(
                    state0=state0,
                    jacobi=float(jacobi_center + rng.normal(scale=0.01)),
                    period=float(period_center + rng.normal(scale=0.05)),
                )
            )
    path = tmp_path / "synthetic_log.jsonl"
    path.write_text("\n".join(json.dumps(rec) for rec in lines) + "\n")
    return path


# --- resolve_system ---------------------------------------------------------


def test_resolve_system_bare_mu_builds_generic_system() -> None:
    system = resolve_system(mu=0.001)
    assert isinstance(system, CR3BPSystem)
    assert system.mu == pytest.approx(0.001)
    assert system.primary == "P1"
    assert system.secondary == "P2"


def test_resolve_system_named_pair_uses_cr3bp_system() -> None:
    system = resolve_system(primary="Earth", secondary="Moon")
    assert system.primary == "Earth"
    assert system.secondary == "Moon"
    assert system.mu == pytest.approx(TRAINING_MU, rel=1e-2)


def test_resolve_system_requires_mu_or_named_pair() -> None:
    with pytest.raises(ValueError, match="must supply"):
        resolve_system()


def test_resolve_system_rejects_partial_named_pair() -> None:
    with pytest.raises(ValueError, match="both primary and secondary"):
        resolve_system(primary="Earth")


# --- lift-vs-delta-log10-mu documentation -----------------------------------


def test_delta_log10_mu_zero_at_training_point() -> None:
    assert delta_log10_mu(TRAINING_MU, TRAINING_MU) == pytest.approx(0.0, abs=1e-12)


def test_expected_lift_for_mu_matches_anchors_at_anchor_points() -> None:
    for anchor in LIFT_ANCHORS:
        estimate = expected_lift_for_mu(anchor.mu)
        assert estimate.estimated_lift == pytest.approx(anchor.lift, rel=1e-6)
        assert estimate.beyond_validated_range is False


def test_expected_lift_for_mu_interpolates_between_anchors() -> None:
    # Halfway (in delta-space) between the mu=0.001 anchor and the Sun-Earth
    # anchor should land strictly between their lift values (30x and 3.5x) --
    # not equal to either endpoint, and not outside the [3.5, 30] envelope.
    mid_delta = 0.5 * (delta_log10_mu(0.001) + delta_log10_mu(LIFT_ANCHORS[-1].mu))
    mid_mu = 10 ** (np.log10(TRAINING_MU) - mid_delta)  # arbitrary side, magnitude only matters
    estimate = expected_lift_for_mu(mid_mu)
    assert 3.5 < estimate.estimated_lift < 30.0
    assert estimate.beyond_validated_range is False


def test_expected_lift_for_mu_flags_beyond_validated_range() -> None:
    # Far beyond the furthest tested anchor (Sun-Earth, delta~3.6).
    tiny_mu = 1e-12
    estimate = expected_lift_for_mu(tiny_mu)
    assert estimate.beyond_validated_range is True
    assert "unvalidated" in estimate.caveat.lower() or "UNVALIDATED" in estimate.method


def test_expected_lift_for_mu_never_returns_negative_lift() -> None:
    estimate = expected_lift_for_mu(1e-20)
    assert estimate.estimated_lift >= 0.0


# --- stability_index ---------------------------------------------------------


def test_stability_index_returns_finite_value_and_note() -> None:
    system = resolve_system(mu=0.01215)
    state0 = np.array([0.85, 0.0, 0.0, 0.0, 0.2, 0.0])
    idx, note = stability_index(system, state0, period=1e-3)
    assert idx is None or np.isfinite(idx)
    assert isinstance(note, str) and note


def test_stability_index_handles_propagation_failure_gracefully(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import cyclerfinder.ml.seed_generation as seed_generation_module

    def _boom(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(seed_generation_module.cr3bp, "propagate", _boom)
    system = resolve_system(mu=0.01215)
    idx, note = stability_index(system, np.zeros(6), period=1.0)
    assert idx is None
    assert "failed" in note


# --- get_default_model -------------------------------------------------------


def test_get_default_model_deterministic_given_same_inputs(tmp_path: Path) -> None:
    path = _synthetic_corpus_file(tmp_path)
    clear_model_cache()
    model_a, corpus_a, train_a = get_default_model(corpus_paths=[path], seed=608, cache=False)
    model_b, corpus_b, train_b = get_default_model(corpus_paths=[path], seed=608, cache=False)
    np.testing.assert_array_equal(train_a, train_b)
    np.testing.assert_array_equal(model_a.cluster_weights, model_b.cluster_weights)
    assert len(corpus_a) == len(corpus_b)


def test_get_default_model_caches_by_default(tmp_path: Path) -> None:
    path = _synthetic_corpus_file(tmp_path)
    clear_model_cache()
    model_a, _corpus_a, _train_a = get_default_model(corpus_paths=[path], seed=608)
    model_b, _corpus_b, _train_b = get_default_model(corpus_paths=[path], seed=608)
    assert model_a is model_b  # identical cached object, not just equal
    clear_model_cache()
    model_c, _corpus_c, _train_c = get_default_model(corpus_paths=[path], seed=608)
    assert model_c is not model_a  # cache cleared -> freshly re-derived object


def test_get_default_model_missing_corpus_raises_filenotfound(tmp_path: Path) -> None:
    empty_dir = tmp_path / "empty_outcome_log"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError, match="no outcome-log files"):
        default_corpus_paths(empty_dir)


# --- generate_and_refine_seeds: contract -------------------------------------


def _small_synthetic_model(tmp_path: Path) -> ClusteredGaussianLatentModel:
    path = _synthetic_corpus_file(tmp_path)
    clear_model_cache()
    model, _corpus, _train = get_default_model(corpus_paths=[path], seed=608, cache=False)
    return model


def test_generate_and_refine_seeds_shape_and_type_contract(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    report = generate_and_refine_seeds(
        5,
        mu=0.01215,
        model=model,
        rng=np.random.default_rng(0),
        compute_stability=False,
    )
    assert isinstance(report, SeedGenerationReport)
    assert report.n_requested == 5
    assert report.n_attempted == 5
    assert len(report.seeds) == 5
    assert 0.0 <= report.converged_and_physically_sane_rate <= 1.0
    assert report.n_physically_sane <= report.n_converged <= report.n_attempted
    for seed in report.seeds:
        assert isinstance(seed, GeneratedSeed)
        assert seed.state0_guess.shape == (6,)
        assert isinstance(seed.period_guess, float)
        assert isinstance(seed.converged, bool)
        assert isinstance(seed.physically_sane, bool)
        if seed.converged:
            assert seed.state0 is not None
            assert seed.period is not None
            assert seed.jacobi is not None
            assert seed.residual is not None
        else:
            assert seed.state0 is None


def test_generate_and_refine_seeds_is_mu_agnostic(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    for mu in (0.01215, 0.001, 3.0e-6):
        report = generate_and_refine_seeds(
            4, mu=mu, model=model, rng=np.random.default_rng(1), compute_stability=False
        )
        assert report.target_mu == pytest.approx(mu)
        assert "mu=" in report.target_label
        assert len(report.seeds) == 4


def test_generate_and_refine_seeds_named_pair_label(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    report = generate_and_refine_seeds(
        3,
        primary="Earth",
        secondary="Moon",
        model=model,
        rng=np.random.default_rng(2),
        compute_stability=False,
    )
    assert report.target_label == "Earth-Moon"


def test_generate_and_refine_seeds_populates_lift_estimate(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    report = generate_and_refine_seeds(
        3, mu=TRAINING_MU, model=model, rng=np.random.default_rng(3), compute_stability=False
    )
    assert report.lift_estimate.target_mu == pytest.approx(TRAINING_MU)
    assert report.lift_estimate.estimated_lift == pytest.approx(12.25, rel=1e-6)


def test_generate_and_refine_seeds_stability_computed_when_requested(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    report = generate_and_refine_seeds(
        20, mu=0.01215, model=model, rng=np.random.default_rng(4), compute_stability=True
    )
    # At least verify the contract: every converged seed either has a finite
    # stability index or an honest failure note, never a silent crash.
    for seed in report.seeds:
        if seed.converged:
            assert seed.stability_index is None or np.isfinite(seed.stability_index)
            assert isinstance(seed.stability_note, str) and seed.stability_note


# --- generate_and_refine_seeds: target_jacobi_bounds oversample/rank --------


@_skip_no_real_corpus
def test_generate_and_refine_seeds_target_jacobi_bounds_filters_and_ranks() -> None:
    model, _corpus, _train = get_default_model()
    bounds = (2.9, 3.1)
    report = generate_and_refine_seeds(
        5,
        primary="Earth",
        secondary="Moon",
        model=model,
        target_jacobi_bounds=bounds,
        max_oversample_factor=30.0,
        rng=np.random.default_rng(608),
        compute_stability=False,
    )
    assert report.n_matched_target_jacobi is not None
    assert len(report.seeds) <= 5
    for seed in report.seeds:
        assert seed.converged and seed.physically_sane
        # ``jacobi`` is only ``None`` for a non-converged seed (see
        # ``_refine_one``'s construction: ``jacobi`` is set in the SAME
        # ``if converged:`` block as ``state0``/``period``/``residual``) --
        # already excluded by the ``converged`` assert just above.
        assert seed.jacobi is not None
        assert bounds[0] <= seed.jacobi <= bounds[1]
    # Ranked by closeness to the bounds' midpoint -- monotonically
    # non-decreasing distance from center along the returned list.
    center = 0.5 * (bounds[0] + bounds[1])
    gaps = [abs(s.jacobi - center) for s in report.seeds if s.jacobi is not None]
    assert gaps == sorted(gaps)


def test_generate_and_refine_seeds_target_jacobi_bounds_honest_shortfall(tmp_path: Path) -> None:
    # An unreachable Jacobi window (far outside is_physically_sane's own
    # bounds) must never be silently padded -- the report must honestly
    # reflect 0 matches within the attempt budget, not error or fabricate.
    model = _small_synthetic_model(tmp_path)
    report = generate_and_refine_seeds(
        5,
        mu=0.01215,
        model=model,
        target_jacobi_bounds=(999.0, 1000.0),
        max_oversample_factor=3.0,
        rng=np.random.default_rng(5),
        compute_stability=False,
    )
    assert report.n_matched_target_jacobi == 0
    assert report.seeds == []


# --- calibrate_cluster_weights_for_mu ---------------------------------------


def test_calibrate_cluster_weights_for_mu_returns_valid_distribution(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    result = calibrate_cluster_weights_for_mu(
        model,
        mu=0.01215,
        n_probe_per_cluster=3,
        rng=np.random.default_rng(6),
    )
    assert isinstance(result, RecalibrationResult)
    k = len(model.cluster_weights)
    assert len(result.calibrated_weights) == k
    assert len(result.per_cluster_probe_sane_rate) == k
    assert np.sum(result.calibrated_weights) == pytest.approx(1.0, rel=1e-9)
    assert all(0.0 <= w <= 1.0 for w in result.calibrated_weights)
    # Means/covariances/encoder must be untouched -- only weights change.
    for a, b in zip(model.cluster_means, result.calibrated_model.cluster_means, strict=True):
        np.testing.assert_array_equal(a, b)
    for a, b in zip(model.cluster_covs, result.calibrated_model.cluster_covs, strict=True):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(model.components, result.calibrated_model.components)


def test_calibrate_cluster_weights_for_mu_sample_shape_still_valid(tmp_path: Path) -> None:
    model = _small_synthetic_model(tmp_path)
    result = calibrate_cluster_weights_for_mu(
        model, mu=0.01215, n_probe_per_cluster=3, rng=np.random.default_rng(7)
    )
    samples = result.calibrated_model.sample(10, rng=np.random.default_rng(8))
    assert samples.shape == (10, 8)
    assert np.all(np.isfinite(samples))


# --- regression checks against #608's/#624's own real-corpus headline lift --


@_skip_no_real_corpus
def test_generate_and_refine_seeds_reproduces_608_lift_floor_in_distribution() -> None:
    """Fast-path regression check (task step 8), floors updated by `#642`.

    `#642` found #608's/#624's original headline lift numbers were measured
    BEFORE ``is_physically_sane`` rejected degenerate Lagrange-point
    equilibria, and re-derived the true numbers from #608's own saved raw
    converged states: the Earth-Moon in-distribution GENERATED arm's real
    (non-equilibrium) physically-sane rate is ~22-27% (not #608's original
    49%, which counted 27/49 equilibria as real orbits), and the BASELINE
    arm's true rate is much lower and noisier than #608's original 4%
    (0-2 real hits per 100 draws in two independent #642 re-derivations, not
    4) -- both floors below widened accordingly to avoid flakiness against
    this now-understood noise, while still catching a silent regression in
    the ported generate-then-refine logic.
    """
    model, _corpus, train = get_default_model()
    n = 100
    report = generate_and_refine_seeds(
        n,
        primary="Earth",
        secondary="Moon",
        model=model,
        rng=np.random.default_rng(608),
        compute_stability=False,
    )
    gen_rate = report.converged_and_physically_sane_rate

    system = resolve_system(primary="Earth", secondary="Moon")
    baseline = uniform_bounding_box_sample(train, n, rng=np.random.default_rng(609))
    n_sane_base = 0
    for row in baseline:
        try:
            if row[6] <= 1e-6:
                continue
            orbit = correct_periodic(system, row[:6], float(row[6]), tol=1e-10, max_iter=30)
            if orbit.converged and is_physically_sane(orbit.state0, orbit.period, orbit.jacobi):
                n_sane_base += 1
        except Exception:
            pass
    base_rate = n_sane_base / n

    assert gen_rate > 0.15  # #642-corrected real rate is ~22-27%; a wide floor avoids flakiness
    if base_rate > 0:
        assert gen_rate / base_rate > 3.0  # #642-corrected real ratio is ~11-27x
    else:
        assert gen_rate > 0.0


@_skip_no_real_corpus
def test_generate_and_refine_seeds_cross_mu_lift_no_longer_assumed_positive() -> None:
    """`#642` FALSIFIED #624's headline cross-mu claim -- this is no longer a
    "must stay positive" regression test.

    #624's original claim (60%/30x at mu=0.001, 7%/3.5x at Sun-Earth,
    "lift TRANSFERS off-distribution") was measured before
    ``is_physically_sane`` rejected degenerate Lagrange-point equilibria.
    `#642` re-derived the true numbers from #624's own saved raw converged
    states and found the GENERATED arm's real (non-equilibrium) hits were
    ZERO at BOTH tested mu (every one of the original 60%/7% "physically
    sane" generated hits was actually an L4/L5 fixed point), while the
    uniform BASELINE retained a handful of genuine real orbits at both --
    i.e. the corrected lift is undefined-or-negative, not merely weaker.
    `#642`'s own live re-run of THIS module's API (different rng draws,
    N=60) independently reproduced the same collapse: 0/60 real generated
    hits at both mu, vs 2/60 real baseline hits at both. This is therefore
    NOT currently a "positive lift holds" regression test -- asserting that
    would re-encode the falsified claim. It stays as a CONTRACT/smoke test
    (the API runs cleanly end-to-end at both cross-mu targets and returns
    well-formed rates) until a proper re-pilot (flagged as follow-on work by
    `#642`) either confirms this collapse is durable or finds it was itself
    an N-too-small artifact. See `#642`'s ``data/OUTSTANDING.md`` bullet
    (AWAITING OPUS ADJUDICATION as of this test's writing) before re-adding
    any directional assertion here.
    """
    model, _corpus, train = get_default_model()
    n = 60
    for mu in (0.001, 3.0034805950690393e-06):
        report = generate_and_refine_seeds(
            n, mu=mu, model=model, rng=np.random.default_rng(624), compute_stability=False
        )
        gen_rate = report.converged_and_physically_sane_rate
        assert 0.0 <= gen_rate <= 1.0

        system = resolve_system(mu=mu)
        baseline = uniform_bounding_box_sample(train, n, rng=np.random.default_rng(625))
        n_sane_base = 0
        for row in baseline:
            try:
                if row[6] <= 1e-6:
                    continue
                orbit = correct_periodic(system, row[:6], float(row[6]), tol=1e-10, max_iter=30)
                if orbit.converged and is_physically_sane(orbit.state0, orbit.period, orbit.jacobi):
                    n_sane_base += 1
            except Exception:
                pass
        base_rate = n_sane_base / n
        assert 0.0 <= base_rate <= 1.0


# --- sanity: fit_clustered_gaussian import still works as expected ----------
# (guards against a stale re-export if orbit_generative's API ever changes)


def test_fit_clustered_gaussian_still_importable_for_get_default_model() -> None:
    features = np.random.default_rng(0).normal(size=(200, 8))
    model = fit_clustered_gaussian(features, n_latent=5, n_clusters=8, seed=608)
    assert model.cluster_weights.shape == (8,)
