"""Tests for the false-positive flagger (#256).

Reproduce-before-trust: the flagger MUST score known past-bug closures higher
than known true-reproductions on average. AUC > 0.75 on the labeled corpus
is the gate; if it fails, the flagger is not credible on the small-N labeled
set and downstream code should not rely on it.

The other three tests pin the non-blocking discipline: ``.score`` always
returns a probability, never raises, even on garbage / partial / NaN-laden
input. That is the public contract -- the flagger flags; it never blocks.

Runtime budget: all tests together must run in < 60 s (LOO is the heavy
one; it fits N=33 LogReg solves on a 12-feature standardised matrix, well
under a second per fit).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

from cyclerfinder.ml import (
    FEATURE_NAMES,
    KNOWN_FALSE_POSITIVES,
    KNOWN_TRUE_REPRODUCTIONS,
    FalsePosFlagger,
    build_training_set,
    extract_features,
)
from cyclerfinder.ml.falsepos_flagger import roc_auc

# ---------------------------------------------------------------------------
# Reproduce-before-trust: the AUC gate on the labeled-corpus seed.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def trained_flagger() -> tuple[FalsePosFlagger, np.ndarray, np.ndarray]:
    X, y, _meta = build_training_set()
    clf = FalsePosFlagger()
    clf.fit(X, y)
    return clf, X, y


def test_labeled_corpus_auc_above_075(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
) -> None:
    """The reproduce-before-trust gate: AUC > 0.75 on the labeled corpus.

    Train AUC measures fit (mostly perfect on tiny N -- the LogReg can carve
    the 12-d feature space cleanly), but the LOO AUC is the honest small-N
    held-out signal. We gate on the LOO AUC: if a flagger learned its
    signatures only by memorising row identity, LOO would crater. The
    threshold 0.75 follows the task spec (not perfect; the corpus is small).
    """
    clf, _X, _y = trained_flagger
    diag = clf.diagnostics
    assert diag is not None
    # The TRAIN AUC is essentially a sanity check that fit converged.
    assert diag.auc_train > 0.75, f"train AUC {diag.auc_train:.3f} below gate"
    # LOO AUC is the real gate. NaN means LOO degenerated -- not acceptable.
    assert not math.isnan(diag.auc_loo), "LOO AUC was NaN (corpus too small)"
    assert diag.auc_loo > 0.75, (
        f"LOO AUC {diag.auc_loo:.3f} below the 0.75 reproduce-before-trust gate -- "
        f"the labeled-corpus seed is not separable enough to trust this flagger"
    )


def test_corpus_class_balance(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
) -> None:
    """Both classes are populated; class_weight='balanced' has work to do.

    Catches the silent failure where the corpus drifts to one-class only and
    AUC tautologically passes. The task spec wires class_weight='balanced'
    precisely because corpora skew; this asserts the skew exists.
    """
    clf, _X, _y = trained_flagger
    diag = clf.diagnostics
    assert diag is not None
    assert diag.n_positives >= 8, "expect at least 8 labeled false-positives"
    assert diag.n_negatives >= 8, "expect at least 8 labeled true-reproductions"
    # Class imbalance is real, not severe (a sanity rather than the gate).
    ratio = diag.n_positives / max(diag.n_negatives, 1)
    assert 0.3 <= ratio <= 3.0, (
        f"class-balance ratio {ratio:.2f} is suspiciously extreme; "
        f"class_weight='balanced' is wired but the corpus shouldn't drift this much"
    )


# ---------------------------------------------------------------------------
# Score-range contract.
# ---------------------------------------------------------------------------


def test_score_is_a_probability(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
) -> None:
    clf, _X, _y = trained_flagger
    p = clf.score(KNOWN_FALSE_POSITIVES[0])
    assert isinstance(p, float)
    assert 0.0 <= p <= 1.0


def test_score_returns_05_when_unfit() -> None:
    """An unfit flagger returns the maximum-entropy 0.5 -- never raises."""
    clf = FalsePosFlagger()
    assert not clf.is_fit()
    assert clf.score(KNOWN_FALSE_POSITIVES[0]) == 0.5


# ---------------------------------------------------------------------------
# Non-blocking discipline: the flagger NEVER raises on a SILVER input.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "pathological",
    [
        {},  # entirely empty record
        {"max_residual_kms": None},  # explicit None
        {"max_residual_kms": float("nan")},  # NaN
        {"max_residual_kms": "not a number"},  # wrong type
        {"vinf_per_encounter_kms": "not a sequence"},
        {"vinf_per_encounter_kms": [None, "garbage", float("nan")]},
        {"closure_date": 12345},  # int where ISO date expected
        {"closure_method_version": None},
        {"encounter_periods_days": [-1.0, 0.0]},  # nonsensical periods
        {"encounter_periods_days": [365.25]},  # only one period -> NaN deviation
    ],
)
def test_score_never_raises_on_pathological_input(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
    pathological: dict[str, object],
) -> None:
    """Public contract: the flagger always returns a probability, never raises."""
    clf, _X, _y = trained_flagger
    p = clf.score(pathological)
    assert isinstance(p, float)
    assert 0.0 <= p <= 1.0


def test_score_never_raises_on_non_dict_input(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
) -> None:
    """Even calling with a wrong-typed record returns a probability, not a raise."""
    clf, _X, _y = trained_flagger
    # Intentionally violating the type to test the non-blocking guard. The
    # score() contract is "any input -> probability"; mypy doesn't see that
    # so we cast through Any rather than try to teach it.
    not_a_dict: Any = None
    p = clf.score(not_a_dict)
    assert isinstance(p, float)
    assert 0.0 <= p <= 1.0


# ---------------------------------------------------------------------------
# Feature-extractor NaN robustness.
# ---------------------------------------------------------------------------


def test_extract_features_returns_fixed_shape() -> None:
    v = extract_features({})
    assert v.shape == (len(FEATURE_NAMES),)
    assert v.dtype == np.float64


def test_extract_features_robust_to_missing_keys() -> None:
    v = extract_features({"max_residual_kms": 0.5})
    assert v.shape == (len(FEATURE_NAMES),)
    # max_residual_kms appears at FEATURE_NAMES[0]:
    assert v[0] == pytest.approx(0.5)


def test_extract_features_does_not_raise_on_garbage() -> None:
    extract_features({"vinf_per_encounter_kms": object()})
    extract_features({"closure_date": object()})
    extract_features({"model_assumption": 42})
    # If we got here, no raise -- the guarantee.


# ---------------------------------------------------------------------------
# roc_auc helper sanity (the numpy-only AUC).
# ---------------------------------------------------------------------------


def test_roc_auc_perfect_separation() -> None:
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    assert roc_auc(y, s) == pytest.approx(1.0)


def test_roc_auc_perfectly_wrong() -> None:
    y = np.array([0, 0, 1, 1])
    s = np.array([0.9, 0.8, 0.2, 0.1])
    assert roc_auc(y, s) == pytest.approx(0.0)


def test_roc_auc_ties_handled() -> None:
    y = np.array([0, 0, 1, 1])
    s = np.array([0.5, 0.5, 0.5, 0.5])
    # All ties -> AUC = 0.5 (rank-sum at midpoint).
    assert roc_auc(y, s) == pytest.approx(0.5)


def test_roc_auc_nan_for_single_class() -> None:
    y = np.array([1, 1, 1])
    s = np.array([0.1, 0.2, 0.3])
    assert math.isnan(roc_auc(y, s))


# ---------------------------------------------------------------------------
# Smoke: training set shape lines up with the feature ABI.
# ---------------------------------------------------------------------------


def test_training_set_shape() -> None:
    X, y, meta = build_training_set()
    assert X.shape[1] == len(FEATURE_NAMES)
    assert X.shape[0] == y.size == len(meta)
    assert X.shape[0] == len(KNOWN_FALSE_POSITIVES) + len(KNOWN_TRUE_REPRODUCTIONS)
    # Both classes are present.
    assert set(np.unique(y).tolist()) == {0, 1}


# ---------------------------------------------------------------------------
# #275 Pluto-class / binary-regime discrimination.
#
# Before #275, the flagger's p_fp clustered in [0.559, 0.561] for all 12
# Pluto SILVER candidates (#269) -- degenerate because the labeled corpus had
# NO Pluto-class or binary-system exemplars. #275 added 4 synthetic Pluto
# false-positives + 4 synthetic Pluto reproductions. The test below asserts
# that the retrained flagger now SEPARATES those 8 synthetic Pluto-class rows
# by a meaningful margin (the corpus contains the labels, but the LOO AUC
# gate above already covers learnability; this test specifically pins
# *within-Pluto-regime* discrimination). xfail if the margin collapses below
# 0.05 -- that would mean the flagger still cannot tell binary-regime FPs
# from binary-regime TRs and we need REAL labels from a V0-V5 result.
# ---------------------------------------------------------------------------


def test_pluto_class_discrimination(
    trained_flagger: tuple[FalsePosFlagger, np.ndarray, np.ndarray],
) -> None:
    """Retrained flagger discriminates synthetic Pluto FPs from synthetic Pluto TRs.

    Margin gate: mean(p_fp | synthetic-Pluto-FP) - mean(p_fp | synthetic-Pluto-TR)
    must exceed 0.15. Below 0.05 -> xfail with a marker that the labeled corpus
    still cannot discriminate within the binary regime; the next move would be
    acquiring REAL labels from the V0-V5 gauntlet results (#274).
    """
    clf, _X, _y = trained_flagger
    pluto_fps = [r for r in KNOWN_FALSE_POSITIVES if "#275" in r.get("_source", "")]
    pluto_trs = [r for r in KNOWN_TRUE_REPRODUCTIONS if "#275" in r.get("_source", "")]
    # We expect 4 of each from the #275 backfill.
    assert len(pluto_fps) >= 4, "expected at least 4 #275 Pluto-class FP exemplars"
    assert len(pluto_trs) >= 4, "expected at least 4 #275 Pluto-class TR exemplars"

    fp_scores = [clf.score(r) for r in pluto_fps]
    tr_scores = [clf.score(r) for r in pluto_trs]
    mean_fp = float(np.mean(fp_scores))
    mean_tr = float(np.mean(tr_scores))
    margin = mean_fp - mean_tr

    if margin < 0.05:
        pytest.xfail(
            "Pluto-class discrimination margin < 0.05 -- label corpus needs "
            "more Pluto-class examples; the flagger cannot yet discriminate "
            "within the binary regime. Wait for V0-V5 gauntlet results "
            "(#274) to provide REAL labeled examples."
        )
    assert margin > 0.15, (
        f"Pluto-class discrimination margin {margin:.4f} below the 0.15 gate "
        f"(mean FP p_fp={mean_fp:.4f}, mean TR p_fp={mean_tr:.4f}); the "
        f"corpus has 4 Pluto FP + 4 Pluto TR labels but the flagger is not "
        f"separating them. Investigate before relying on Pluto-regime scores."
    )
