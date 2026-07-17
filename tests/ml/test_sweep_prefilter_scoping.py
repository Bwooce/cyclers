"""Tests for the #317 scoping investigation's building blocks.

Same convention as ``tests/ml/test_orbit_generative.py``: these exercise the
LOGIC (feature extraction, dataset assembly, classifier fit, recall/skip
table, gate-efficiency aggregation) with small constructed fixtures -- code
correctness only. The actual #317 REAL-corpus numbers (AUC ~0.67, ~66% vs.
~63% held-out accuracy, negligible safe-recall compute savings, and the
Lambert/bend-gate efficiency stats) are pinned in
``scripts/run_317_prefilter_scoping.py``'s own report, not here -- see that
script and ``data/found/317_prefilter_scoping/summary.json``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import _r1_r2, jacobi_constant
from cyclerfinder.ml.sweep_prefilter_scoping import (
    CHEAP_FEATURE_NAMES,
    assemble_prefilter_dataset,
    auc_score,
    extract_cheap_features,
    fit_logistic_prefilter,
    gate_efficiency_from_summary_records,
    recall_skip_table,
)


def _corrector_record(
    *,
    solver: str = "cr3bp.correct_periodic",
    mu: float = 0.01215,
    state0: list[float] | None = None,
    period_guess: float = 3.0,
    converged: bool = True,
) -> dict[str, Any]:
    return {
        "solver": solver,
        "inputs": {
            "mu": mu,
            "state0_guess": state0 or [0.85, 0.0, 0.1, 0.0, 0.2, 0.0],
            "period_guess": period_guess,
        },
        "outcome": {
            "converged": converged,
            "jacobi": 3.0,
            "period": period_guess,
            "residual": 1e-10,
        },
    }


def test_cheap_feature_names_length_matches_extraction() -> None:
    rec = _corrector_record()
    result = extract_cheap_features(rec)
    assert result is not None
    feats, conv = result
    assert feats.shape == (len(CHEAP_FEATURE_NAMES),)
    assert conv is True


def test_extract_cheap_features_matches_hand_computed_values() -> None:
    mu = 0.01215
    state0 = [0.85, 0.01, 0.02, 0.03, 0.2, 0.04]
    rec = _corrector_record(mu=mu, state0=state0, period_guess=2.5, converged=False)
    result = extract_cheap_features(rec)
    assert result is not None
    feats, conv = result
    assert conv is False
    expected_jacobi = jacobi_constant(np.array(state0), mu)
    expected_r1, expected_r2 = _r1_r2(state0[0], state0[1], state0[2], mu)
    names = dict(zip(CHEAP_FEATURE_NAMES, feats, strict=True))
    assert names["jacobi_guess"] == pytest.approx(expected_jacobi)
    assert names["period_guess"] == pytest.approx(2.5)
    assert names["r1_guess"] == pytest.approx(expected_r1)
    assert names["r2_guess"] == pytest.approx(expected_r2)
    assert names["speed_guess"] == pytest.approx(np.linalg.norm(state0[3:]))
    assert names["posn_guess"] == pytest.approx(np.linalg.norm(state0[:3]))
    assert names["abs_z0_guess"] == pytest.approx(abs(state0[2]))
    assert names["abs_vz0_guess"] == pytest.approx(abs(state0[5]))


def test_extract_cheap_features_rejects_wrong_solver_and_malformed() -> None:
    assert extract_cheap_features(_corrector_record(solver="other.solver")) is None
    assert extract_cheap_features({"solver": "cr3bp.correct_periodic"}) is None
    bad_shape = _corrector_record(state0=[0.85, 0.0, 0.0])
    assert extract_cheap_features(bad_shape) is None


def test_assemble_prefilter_dataset_scans_and_skips_malformed(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    lines = [
        _corrector_record(converged=True),
        _corrector_record(converged=False, period_guess=4.0),
        _corrector_record(solver="other.solver"),  # dropped
        {"not": "even close"},  # dropped (missing solver key entirely)
    ]
    path.write_text("\n".join(json.dumps(x) for x in lines) + "\n" + "{not valid json\n" + "\n")
    dataset = assemble_prefilter_dataset([path])
    # n_scanned counts records `iter_outcome_records` successfully parses as
    # JSON (it already silently skips the invalid-JSON line and the blank
    # line before yielding) -- so this is exactly `len(lines)`, not more.
    assert dataset.n_scanned == len(lines)
    assert len(dataset) == 2
    assert dataset.features.shape == (2, len(CHEAP_FEATURE_NAMES))
    assert list(dataset.converged) == [True, False]


def test_assemble_prefilter_dataset_raises_on_empty(tmp_path: Path) -> None:
    path = tmp_path / "log.jsonl"
    path.write_text(json.dumps(_corrector_record(solver="other.solver")) + "\n")
    with pytest.raises(ValueError, match="no record survived"):
        assemble_prefilter_dataset([path])


def test_fit_logistic_prefilter_recovers_a_clearly_separable_signal() -> None:
    """Code-correctness only: a synthetic dataset where feature 0 alone

    perfectly determines the label should be fit near-perfectly by logistic
    regression (sanity check that the optimizer/gradient are wired up
    correctly), regardless of what the REAL corpus's signal strength turns
    out to be.
    """
    rng = np.random.default_rng(0)
    n = 400
    x = rng.normal(size=(n, len(CHEAP_FEATURE_NAMES)))
    y = x[:, 0] > 0.0
    model = fit_logistic_prefilter(x, y, seed=0)
    proba = model.predict_proba(x)
    acc = ((proba >= 0.5) == y).mean()
    assert acc > 0.95
    assert auc_score(model.decision_scores(x), y) > 0.95


def test_fit_logistic_prefilter_rejects_single_class() -> None:
    x = np.zeros((10, 3))
    y = np.zeros(10, dtype=bool)
    with pytest.raises(ValueError, match="both classes"):
        fit_logistic_prefilter(x, y)


def test_auc_score_random_scores_near_half_perfect_scores_is_one() -> None:
    rng = np.random.default_rng(1)
    labels = np.array([True] * 500 + [False] * 500)
    random_scores = rng.normal(size=1000)
    assert auc_score(random_scores, labels) == pytest.approx(0.5, abs=0.05)

    perfect_scores = np.concatenate([np.full(500, 1.0), np.full(500, 0.0)])
    assert auc_score(perfect_scores, labels) == pytest.approx(1.0)


def test_recall_skip_table_monotonic_in_target_recall() -> None:
    rng = np.random.default_rng(2)
    n = 2000
    labels = rng.random(n) < 0.4
    # Score correlates with label but imperfectly (mimics the real corpus's
    # weak-but-real signal), so there is a genuine recall/skip trade-off.
    scores = np.where(labels, rng.normal(1.0, 1.0, n), rng.normal(0.0, 1.0, n))
    table = recall_skip_table(scores, labels, target_recalls=(0.999, 0.99, 0.9, 0.5))
    skips = [row["skip_fraction_of_negatives"] for row in table]
    # Lower target recall (willing to lose more real hits) must never SKIP
    # fewer negatives than a higher target recall -- monotonic trade-off.
    assert skips == sorted(skips)
    for row in table:
        assert row["actual_recall"] >= row["target_recall"] - 1e-6


def test_gate_efficiency_from_summary_records_aggregates_correctly() -> None:
    records: list[dict[str, Any]] = [
        {"_meta": True, "unrelated": "record"},  # skipped: missing required keys
        {
            "kind": "direction_summary",
            "n_evaluated": 192,
            "n_infeasible": 122,
            "n_subgate_residual_only": 70,
            "n_all_gates_passed": 0,
        },
        {
            "kind": "direction_summary",
            "n_evaluated": 100,
            "n_infeasible": 60,
            "n_subgate_residual_only": 40,
            "n_all_gates_passed": 10,
        },
    ]
    summary = gate_efficiency_from_summary_records(records)
    assert summary.n_evaluated == 292
    assert summary.n_infeasible == 182
    assert summary.n_subgate_reaches_expensive_gate == 110
    assert summary.n_all_gates_passed == 10
    assert summary.frac_rejected_by_cheap_gate == pytest.approx(182 / 292)
    assert summary.frac_reaching_expensive_gate == pytest.approx(110 / 292)
    assert summary.hit_rate_of_expensive_gate_pool == pytest.approx(10 / 110)


def test_gate_efficiency_handles_empty_input() -> None:
    summary = gate_efficiency_from_summary_records([])
    assert summary.n_evaluated == 0
    assert summary.frac_rejected_by_cheap_gate == 0.0
    assert summary.hit_rate_of_expensive_gate_pool == 0.0
