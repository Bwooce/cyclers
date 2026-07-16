"""#317 scoping: is a learned pre-filter for "sweep-impossible regions" warranted?

``#317`` started as a one-line backlog idea ("PINN-based pre-filter for
sweep-impossible regions") with no scoping detail. This module is the
scoping + bounded-feasibility investigation, built on this project's OWN
history rather than synthetic data, following the ``#608``/``#614`` precedent
(``cyclerfinder.ml.orbit_generative``): no new heavy ML dependency (no
torch/jax/sklearn -- checked ``pyproject.toml``), a from-scratch numpy/scipy
model where a model is warranted at all.

**Two concrete "sweep-impossible region" regimes exist in this codebase, and
they turned out to need OPPOSITE verdicts:**

1. **Lambert/resonance-construction sweeps** (``scripts/enumerate_563_*``,
   ``enumerate_575_*``, ``enumerate_576_*``, ``enumerate_599_*``,
   ``enumerate_600_*``, ``enumerate_607_*``, ``enumerate_609_*`` and their
   ``scan_558_*`` ancestor): these ALREADY have a fast, closed-form,
   EXACT physical pre-gate two layers deep -- Lambert-solution existence
   (``residual_at_point`` returning ``None``) and the ``#324`` patched-conic
   bend gate (``cyclerfinder.search.physical_sanity.candidate_passes_physical_gate``,
   an ``arcsin`` formula). :func:`gate_efficiency_from_summary_records` reads
   the ``direction_summary``/``sequence_summary`` records these scripts
   already emit and shows the cheap Lambert-existence gate alone rejects
   ~70-87% of every evaluated candidate essentially for free, and the
   candidates that DO reach the remaining ~10-30% expensive step
   (bend + DOP853 cross-check) are exactly where these sweeps found their
   real closures (a ~25% hit rate in the systems with any positive
   closures at all -- see ``scripts/run_317_prefilter_scoping.py``'s report).
   A learned pre-filter here would have to out-perform gates that are
   already EXACT and effectively free to compute -- there is no honest room
   for it to add value, and pruning that ~10-30% pool further risks
   discarding the real hits it exists to find (exactly the "corrector could
   succeed but usually doesn't" trap ``#317``'s own scoping brief warns
   against).

2. **Blind CR3BP periodic-orbit correction** (``scripts/search_campaign_daemon.py``'s
   Phase B, the ``#210`` outcome-log corpus's dominant contributor): this
   regime has NO existing pre-gate at all -- a JPL seed is perturbed and
   fed straight to ``cyclerfinder.search.cr3bp_periodic.correct_periodic``
   (a Newton/STM shooter, ~250 ms/call, roughly 1.7e5x the cost of the
   cheapest feature below) with no check first, and ~62-63% of calls do not
   converge. This IS the regime ``#317`` was written for: an expensive
   solver invoked blindly with a real non-convergence rate. This module's
   own measurement (below) is nonetheless a NEGATIVE: cheap features
   computable before the corrector runs (the guess's own Jacobi constant,
   period, primary/secondary proximities, out-of-plane amplitude) carry only
   WEAK signal about eventual convergence (AUC ~0.67, ~66% held-out accuracy
   vs. a ~63% majority baseline), and at any recall level safe enough for a
   discovery program that treats real hits as precious (this project's own
   `MEMORY.md`: "Novel hits are RARE"), the achievable compute savings is
   negligible (well under 2% at 99% recall on the converged class; near-zero
   at 99.9%). This is consistent with most non-convergence here being Newton
   basin-of-attraction sensitivity to a random perturbation, not genuine
   physical infeasibility -- exactly the "search sparsity, not infeasibility"
   case the scoping brief flags as NOT the right thing to learn to reject.

**Why a classical classifier, not a literal PINN.** A PINN specifically means
a network trained with a loss term penalizing violation of a differential
equation via autodiff through both the network and the physics residual.
Nothing here needs that: the candidate labels (converged / did-not-converge,
gate-pass / gate-fail) are already observed outcomes of running the real
solver or the real closed-form gate, not points where we need to enforce an
ODE residual through a differentiable surrogate. A plain classifier over
physically-meaningful, cheaply-computed features is the right tool for
"predict this observed binary outcome from these inputs" -- which is all
either regime above actually needs. See ``fit_logistic_prefilter`` below: a
hand-rolled scipy L-BFGS-B logistic regression, the simplest classifier that
could show a real effect; if THIS shows no real effect, a deeper model is not
the missing ingredient (the features themselves carry too little signal).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from cyclerfinder.core.cr3bp import _r1_r2, jacobi_constant
from cyclerfinder.ml.orbit_generative import iter_outcome_records

__all__ = [
    "CHEAP_FEATURE_NAMES",
    "GateEfficiencySummary",
    "LogisticPrefilterModel",
    "PrefilterDataset",
    "assemble_prefilter_dataset",
    "auc_score",
    "extract_cheap_features",
    "fit_logistic_prefilter",
    "gate_efficiency_from_summary_records",
    "recall_skip_table",
]

# Cheap, closed-form, pre-corrector features for one `cr3bp.correct_periodic`
# guess: the guess's own Jacobi constant, its period guess, its distance to
# each primary in the rotating frame, its speed/position magnitude, and its
# out-of-plane amplitude/rate. Every one of these is O(1) arithmetic on the
# guess itself -- no propagation, no STM, nothing the corrector does.
CHEAP_FEATURE_NAMES: tuple[str, ...] = (
    "jacobi_guess",
    "period_guess",
    "r1_guess",
    "r2_guess",
    "speed_guess",
    "posn_guess",
    "abs_z0_guess",
    "abs_vz0_guess",
)


def extract_cheap_features(
    record: dict[str, Any],
    *,
    solver: str = "cr3bp.correct_periodic",
) -> tuple[NDArray[np.float64], bool] | None:
    """Compute :data:`CHEAP_FEATURE_NAMES` + the converged label from ONE outcome-log record.

    Returns ``None`` if ``record`` is not a well-formed record for ``solver``
    (wrong solver, missing fields, non-finite feature) -- callers scanning a
    real, occasionally-torn log must skip these, never crash on them.
    """
    if record.get("solver") != solver:
        return None
    try:
        inputs = record["inputs"]
        outcome = record["outcome"]
        mu = float(inputs["mu"])
        state0 = np.asarray(inputs["state0_guess"], dtype=float)
        period_guess = float(inputs["period_guess"])
        converged = bool(outcome["converged"])
    except (KeyError, TypeError, ValueError):
        return None
    if state0.shape != (6,):
        return None
    try:
        jacobi_guess = jacobi_constant(state0, mu)
        r1, r2 = _r1_r2(state0[0], state0[1], state0[2], mu)
    except Exception:
        return None
    speed_guess = float(np.linalg.norm(state0[3:]))
    posn_guess = float(np.linalg.norm(state0[:3]))
    feats = np.array(
        [
            jacobi_guess,
            period_guess,
            r1,
            r2,
            speed_guess,
            posn_guess,
            abs(float(state0[2])),
            abs(float(state0[5])),
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(feats)):
        return None
    return feats, converged


@dataclass(frozen=True)
class PrefilterDataset:
    """Assembled cheap-feature dataset for the CR3BP-corrector pre-filter question."""

    features: NDArray[np.float64]  # (n, len(CHEAP_FEATURE_NAMES))
    converged: NDArray[np.bool_]  # (n,)
    n_scanned: int

    def __len__(self) -> int:
        return int(self.features.shape[0])


def assemble_prefilter_dataset(
    paths: Iterable[Path],
    *,
    solver: str = "cr3bp.correct_periodic",
) -> PrefilterDataset:
    """Scan #210 outcome-log JSONL ``paths`` into a :class:`PrefilterDataset`.

    Every well-formed ``solver`` record contributes one row, regardless of
    ``converged`` -- unlike ``orbit_generative.assemble_corpus`` (which keeps
    only physically-sane CONVERGED rows to build a generative seed corpus),
    this dataset's whole point is the converged/not-converged LABEL itself.
    """
    feats: list[NDArray[np.float64]] = []
    labels: list[bool] = []
    n_scanned = 0
    for record in iter_outcome_records(paths):
        n_scanned += 1
        result = extract_cheap_features(record, solver=solver)
        if result is None:
            continue
        f, conv = result
        feats.append(f)
        labels.append(conv)
    if not feats:
        raise ValueError(f"no record survived feature extraction for solver={solver!r}")
    return PrefilterDataset(
        features=np.vstack(feats),
        converged=np.array(labels, dtype=bool),
        n_scanned=n_scanned,
    )


@dataclass(frozen=True)
class LogisticPrefilterModel:
    """Standardize -> logistic regression: the simplest classifier that could show an effect."""

    mean: NDArray[np.float64]
    scale: NDArray[np.float64]
    weights: NDArray[np.float64]  # (1 + n_features,): [bias, *coeffs]
    feature_names: tuple[str, ...]

    def decision_scores(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        """Raw logit scores (higher -> more likely to converge)."""
        standardized = (np.asarray(features, dtype=float) - self.mean) / self.scale
        return standardized @ self.weights[1:] + self.weights[0]

    def predict_proba(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        z = np.clip(self.decision_scores(features), -30.0, 30.0)
        return 1.0 / (1.0 + np.exp(-z))


def _sigmoid(z: NDArray[np.float64]) -> NDArray[np.float64]:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30.0, 30.0)))


def _logistic_nll_grad(
    w: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    l2: float,
) -> tuple[float, NDArray[np.float64]]:
    z = x @ w[1:] + w[0]
    p = _sigmoid(z)
    eps = 1e-9
    nll = -np.mean(y * np.log(p + eps) + (1.0 - y) * np.log(1.0 - p + eps))
    nll += l2 * float(np.sum(w[1:] ** 2))
    grad_w = x.T @ (p - y) / len(y) + 2.0 * l2 * w[1:]
    grad_b = np.mean(p - y)
    return float(nll), np.concatenate([[grad_b], grad_w])


def fit_logistic_prefilter(
    features: NDArray[np.float64],
    converged: NDArray[np.bool_],
    *,
    l2: float = 1e-3,
    seed: int = 0,
) -> LogisticPrefilterModel:
    """Fit a standardized logistic regression predicting ``converged`` from ``features``.

    scipy L-BFGS-B, analytic gradient -- same "from-scratch, no new ML
    dependency" discipline as ``orbit_generative.fit_autoencoder``. ``seed``
    only seeds the (deterministic, zero) initial weights' tie-break; included
    for interface symmetry with the rest of this project's ``fit_*`` helpers.
    """
    x = np.asarray(features, dtype=float)
    y = np.asarray(converged, dtype=float)
    if x.ndim != 2 or x.shape[0] != y.shape[0]:
        raise ValueError(f"features shape {x.shape} incompatible with converged shape {y.shape}")
    if x.shape[0] < 2 or len(np.unique(y)) < 2:
        raise ValueError("need at least 2 rows and both classes present to fit a classifier")
    mean = x.mean(axis=0)
    scale = x.std(axis=0)
    scale[scale == 0.0] = 1.0
    x_std = (x - mean) / scale

    rng = np.random.default_rng(seed)
    w0 = rng.normal(scale=1e-3, size=x.shape[1] + 1)
    result = minimize(
        _logistic_nll_grad,
        w0,
        args=(x_std, y, l2),
        jac=True,
        method="L-BFGS-B",
    )
    return LogisticPrefilterModel(
        mean=mean,
        scale=scale,
        weights=result.x,
        feature_names=CHEAP_FEATURE_NAMES,
    )


def auc_score(scores: NDArray[np.float64], labels: NDArray[np.bool_]) -> float:
    """Mann-Whitney-U-based AUC: P(score(positive) > score(negative))."""
    from scipy.stats import rankdata

    labels = np.asarray(labels, dtype=bool)
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUC requires both classes present")
    ranks = rankdata(scores)
    return float((ranks[labels].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def recall_skip_table(
    scores: NDArray[np.float64],
    labels: NDArray[np.bool_],
    *,
    target_recalls: tuple[float, ...] = (0.999, 0.99, 0.95, 0.90),
) -> list[dict[str, float]]:
    """For each target recall on the converged class, report the achievable skip fraction.

    A pre-filter that skips the corrector below some score threshold trades
    lost real convergences (recall < 1) for skipped-negative compute savings.
    This is the decision-relevant table for #317: "if we're only willing to
    silently lose X% of real hits, how much compute do we actually save?" --
    NOT raw accuracy/AUC, which don't answer that question directly.
    """
    labels = np.asarray(labels, dtype=bool)
    p_pos = scores[labels]
    p_neg = scores[~labels]
    if len(p_pos) == 0 or len(p_neg) == 0:
        raise ValueError("recall_skip_table requires both classes present")
    rows: list[dict[str, float]] = []
    for target in target_recalls:
        # `method="lower"` picks an actual data point rather than
        # interpolating between two -- guarantees `actual_recall >= target`
        # (interpolation can otherwise land the threshold a hair above the
        # exact target quantile and under-shoot it).
        threshold = float(np.quantile(p_pos, 1.0 - target, method="lower"))
        actual_recall = float((p_pos >= threshold).mean())
        skip_frac = float((p_neg < threshold).mean())
        rows.append(
            {
                "target_recall": float(target),
                "threshold": threshold,
                "actual_recall": actual_recall,
                "skip_fraction_of_negatives": skip_frac,
            }
        )
    return rows


@dataclass(frozen=True)
class GateEfficiencySummary:
    """Aggregate cheap-gate-vs-expensive-gate cost split over a sweep's own summary records."""

    n_evaluated: int
    n_infeasible: int
    n_subgate_reaches_expensive_gate: int
    n_all_gates_passed: int

    @property
    def frac_rejected_by_cheap_gate(self) -> float:
        return self.n_infeasible / self.n_evaluated if self.n_evaluated else 0.0

    @property
    def frac_reaching_expensive_gate(self) -> float:
        return self.n_subgate_reaches_expensive_gate / self.n_evaluated if self.n_evaluated else 0.0

    @property
    def hit_rate_of_expensive_gate_pool(self) -> float:
        """Of the candidates that reach the expensive gate, what fraction actually close?

        A LOW hit rate here means the expensive-gate pool is mostly further
        junk (room for a pre-filter, in principle); a HIGH hit rate means
        that pool is exactly where real discoveries live and should NOT be
        pruned further by an approximate model.
        """
        if self.n_subgate_reaches_expensive_gate == 0:
            return 0.0
        return self.n_all_gates_passed / self.n_subgate_reaches_expensive_gate


def gate_efficiency_from_summary_records(
    records: Iterable[dict[str, Any]],
) -> GateEfficiencySummary:
    """Aggregate ``direction_summary``/``sequence_summary`` records from an enumerate/scan sweep.

    These sweeps (``scripts/enumerate_563_*`` and siblings, descended from
    ``scripts/scan_558_uranus_all_pairs_offset_sweep.py``) already emit
    per-direction/per-sequence records with ``n_evaluated``, ``n_infeasible``
    (rejected by the cheap closed-form Lambert-existence/residual gate),
    ``n_subgate_residual_only`` (passed that cheap gate, reaches the
    expensive bend+DOP853 ``gate_candidate`` step), and
    ``n_all_gates_passed``. Records missing any of these keys (``_meta``,
    other ``kind`` values) are skipped.
    """
    n_evaluated = n_infeasible = n_subgate = n_passed = 0
    required = ("n_evaluated", "n_infeasible", "n_subgate_residual_only", "n_all_gates_passed")
    for rec in records:
        if not all(k in rec for k in required):
            continue
        n_evaluated += int(rec["n_evaluated"])
        n_infeasible += int(rec["n_infeasible"])
        n_subgate += int(rec["n_subgate_residual_only"])
        n_passed += int(rec["n_all_gates_passed"])
    return GateEfficiencySummary(
        n_evaluated=n_evaluated,
        n_infeasible=n_infeasible,
        n_subgate_reaches_expensive_gate=n_subgate,
        n_all_gates_passed=n_passed,
    )
