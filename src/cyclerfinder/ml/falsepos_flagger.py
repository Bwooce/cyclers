"""The false-positive flagger -- a NON-BLOCKING guard rail (#256).

Public API: :class:`FalsePosFlagger`. Given a SILVER closure record, it
returns a probability in ``[0, 1]`` that the record resembles one of our own
past bugs (see ``falsepos_features`` for the catalogued signatures). It
NEVER raises on a SILVER input, and it NEVER auto-rejects -- a downstream
caller routes high-probability scores to human re-check, never to auto-
promotion. ``.fit`` reports cross-val-like ROC AUC and a feature-importance
diagnostic; ``.score`` is the one-call probability.

Implementation choice (#256, Option 2): an L2-regularised LOGISTIC REGRESSION
fit by ``scipy.optimize.minimize`` (BFGS on the cross-entropy + L2 loss), with
``class_weight='balanced'`` (inverse-frequency weights) and a numpy-only
``roc_auc_score`` equivalent (rank-sum statistic). No new dependencies -- the
project pattern is lean deps + targeted modules, and this is the smallest
defensible ML item per the task spec. RandomForest was the alternative the
spec accepts; LogReg is sufficient at this corpus size and stays linear-
inspectable via ``coef * std(X)`` feature importance.

Discipline:

* **Reproduce before trust.** The trained flagger must score known false-
  positives higher than known true-reproductions on average -- AUC > 0.75 on
  the labeled corpus is the gate the test enforces. If the gate fails, the
  flagger MUST NOT be relied on.
* **Imputation is deterministic.** NaN features are filled with the training
  median (then 0.0 if a feature column is all-NaN); standardisation uses the
  training mean / std (std floored at 1e-8 to avoid divide-by-zero).
* **Non-blocking discipline.** ``.score`` returns ``0.5`` for any failure
  path (extractor raise, fit-time bug, pathological input). It never raises.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import minimize

from cyclerfinder.ml.falsepos_features import FEATURE_NAMES, extract_features


def roc_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Rank-sum AUC equivalent (numpy-only, ties-aware).

    Mathematically AUC = ``(U / (n_pos * n_neg))`` where ``U`` is the Mann-
    Whitney U statistic. With ``R_pos = sum of ranks of positives`` (ties
    averaged), ``U = R_pos - n_pos * (n_pos + 1) / 2``. Returns ``nan`` if
    either class is empty.
    """
    y_true = np.asarray(y_true, dtype=np.int64).ravel()
    y_score = np.asarray(y_score, dtype=np.float64).ravel()
    if y_true.shape != y_score.shape or y_true.size == 0:
        return float("nan")
    n_pos = int((y_true == 1).sum())
    n_neg = int((y_true == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Average ranks handle ties (the "method='average'" branch).
    order = np.argsort(y_score, kind="mergesort")
    sorted_scores = y_score[order]
    ranks = np.empty_like(order, dtype=np.float64)
    n = y_score.size
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_scores[j] == sorted_scores[i]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-based ranks, average of [i+1..j]
        ranks[order[i:j]] = avg_rank
        i = j
    r_pos = ranks[y_true == 1].sum()
    u = r_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


@dataclass(frozen=True)
class FitDiagnostics:
    """Diagnostic block reported by :meth:`FalsePosFlagger.fit`.

    Captures the AUC on the training set, the leave-one-out AUC (the
    small-N stand-in for held-out generalisation), feature importance, and
    the class-balance counts. ``auc_train`` measures fit, ``auc_loo`` is the
    one the reproduce-before-trust gate reads.
    """

    n_samples: int
    n_positives: int
    n_negatives: int
    auc_train: float
    auc_loo: float
    feature_importance: tuple[tuple[str, float], ...]


class FalsePosFlagger:
    """Non-blocking probabilistic flag for SILVER closures resembling past bugs.

    Use::

        clf = FalsePosFlagger()
        X, y, _meta = build_training_set()
        diag = clf.fit(X, y)
        p = clf.score(silver_record_dict)  # in [0, 1]; never raises

    The score is ``P(false_positive | features)``. Probabilities above ~0.5
    suggest routing the candidate to human re-check; the threshold belongs
    to the caller (the flagger only reports the probability).
    """

    # ------------------------------------------------------------------
    # Construction / state.
    # ------------------------------------------------------------------
    def __init__(
        self,
        *,
        l2: float = 1.0,
        class_weight_balanced: bool = True,
        random_state: int = 0,
    ) -> None:
        self._l2: float = float(l2)
        self._class_weight_balanced: bool = bool(class_weight_balanced)
        self._random_state: int = int(random_state)
        # Set by .fit (None == not yet fit).
        self._beta: np.ndarray | None = None  # shape (n_features + 1,) -- last is bias.
        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._x_median: np.ndarray | None = None
        self._x_std_unstd: np.ndarray | None = None  # std on imputed (un-standardised) X.
        self._fit_diagnostics: FitDiagnostics | None = None

    # ------------------------------------------------------------------
    # Core fit/predict math.
    # ------------------------------------------------------------------
    @staticmethod
    def _sigmoid(z: np.ndarray) -> np.ndarray:
        # Numerically-stable sigmoid (clip in exp space).
        out = np.empty_like(z, dtype=np.float64)
        pos = z >= 0
        neg = ~pos
        out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
        ez = np.exp(z[neg])
        out[neg] = ez / (1.0 + ez)
        return out

    @staticmethod
    def _impute(
        X: np.ndarray, *, median: np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Replace NaN with the column median (training-time computed)."""
        X = np.asarray(X, dtype=np.float64).copy()
        if median is None:
            with np.errstate(all="ignore"):
                med = np.nanmedian(X, axis=0)
            med = np.where(np.isnan(med), 0.0, med)
        else:
            med = median
        mask = np.isnan(X)
        if mask.any():
            for j in range(X.shape[1]):
                col_mask = mask[:, j]
                if col_mask.any():
                    X[col_mask, j] = med[j]
        return X, med

    @staticmethod
    def _standardise(
        X: np.ndarray,
        *,
        mean: np.ndarray | None = None,
        std: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mu = X.mean(axis=0) if mean is None else mean
        if std is None:
            sigma = X.std(axis=0)
            sigma = np.where(sigma < 1e-8, 1.0, sigma)
        else:
            sigma = std
        Xs = (X - mu) / sigma
        return Xs, mu, sigma

    def _loss_and_grad(
        self,
        beta: np.ndarray,
        Xb: np.ndarray,  # (N, p+1), bias col appended
        y: np.ndarray,  # (N,) in {0, 1}
        w: np.ndarray,  # (N,) per-sample weights
    ) -> tuple[float, np.ndarray]:
        """Weighted cross-entropy + L2 (no penalty on the bias)."""
        z = Xb @ beta
        p = self._sigmoid(z)
        # log-loss, with numerical guard:
        eps = 1e-15
        ll = -(w * (y * np.log(p + eps) + (1.0 - y) * np.log(1.0 - p + eps))).sum()
        # L2 (skip last entry = bias):
        l2 = 0.5 * self._l2 * float(np.dot(beta[:-1], beta[:-1]))
        loss = ll / max(w.sum(), 1e-12) + l2
        # grad
        grad_ll = (Xb.T @ (w * (p - y))) / max(w.sum(), 1e-12)
        reg = np.concatenate([self._l2 * beta[:-1], np.array([0.0])])
        grad = grad_ll + reg
        return float(loss), grad

    def _fit_beta(self, X_std: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
        n, p = X_std.shape
        Xb = np.concatenate([X_std, np.ones((n, 1))], axis=1)
        beta0 = np.zeros(p + 1, dtype=np.float64)

        def _objective(b: np.ndarray) -> tuple[float, np.ndarray]:
            return self._loss_and_grad(b, Xb, y, w)

        # scipy-stubs wants the (fun + jac=True) overload signed (ndarray) ->
        # (float, ndarray); the explicit local def hands the stubs that shape.
        res = minimize(
            _objective,
            beta0,
            jac=True,
            method="BFGS",
            options={"maxiter": 500, "gtol": 1e-8},
        )
        # BFGS may return non-success on tiny corpora; the iterate is still useful
        # as long as it's finite. If it isn't, fall back to zeros (uniform 0.5).
        beta = res.x if res.x is not None and np.all(np.isfinite(res.x)) else beta0
        return np.asarray(beta, dtype=np.float64)

    def _balanced_weights(self, y: np.ndarray) -> np.ndarray:
        n = y.size
        if not self._class_weight_balanced or n == 0:
            return np.ones(n, dtype=np.float64)
        n_pos = max(int((y == 1).sum()), 1)
        n_neg = max(int((y == 0).sum()), 1)
        w_pos = n / (2.0 * n_pos)
        w_neg = n / (2.0 * n_neg)
        return np.where(y == 1, w_pos, w_neg).astype(np.float64)

    # ------------------------------------------------------------------
    # Public fit / score.
    # ------------------------------------------------------------------
    def fit(self, X: np.ndarray, y: np.ndarray) -> FitDiagnostics:
        """Fit on the labeled corpus and return diagnostics.

        Reports train AUC + leave-one-out AUC + feature importance
        ``|coef| * std(X)`` so the caller can verify the
        reproduce-before-trust gate (AUC > 0.75 on the labeled set).
        """
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.int64).ravel()
        if X.ndim != 2 or X.shape[0] != y.size or X.shape[1] != len(FEATURE_NAMES):
            raise ValueError(
                f"X must be (N, {len(FEATURE_NAMES)}) and y (N,); "
                f"got X.shape={X.shape}, y.shape={y.shape}"
            )

        # Impute (median) then standardise; remember training stats.
        X_imp, med = self._impute(X)
        X_std, mu, sigma = self._standardise(X_imp)
        w = self._balanced_weights(y)
        beta = self._fit_beta(X_std, y, w)

        # Cache training state.
        self._beta = beta
        self._x_mean = mu
        self._x_std = sigma
        self._x_median = med
        self._x_std_unstd = X_imp.std(axis=0)

        # AUC on training set.
        p_train = self._predict_proba_from_std(X_std)
        auc_train = roc_auc(y, p_train)

        # Leave-one-out AUC -- the honest small-N gate.
        auc_loo = self._leave_one_out_auc(X, y)

        # Feature importance ~ |coef| * std(X) on the unstandardised scale.
        # (We have standardised beta -- multiply by inverse-std of standardisation
        # to get the per-unit coefficient, then by std on the original scale to
        # get the contribution. The two cancel because we standardised with
        # exactly that std, so the standardised coef IS the contribution.)
        importance = np.abs(beta[:-1])
        order = np.argsort(-importance)
        feat_imp = tuple((FEATURE_NAMES[i], float(importance[i])) for i in order)

        diag = FitDiagnostics(
            n_samples=int(y.size),
            n_positives=int((y == 1).sum()),
            n_negatives=int((y == 0).sum()),
            auc_train=float(auc_train),
            auc_loo=float(auc_loo),
            feature_importance=feat_imp,
        )
        self._fit_diagnostics = diag
        return diag

    def _predict_proba_from_std(self, X_std: np.ndarray) -> np.ndarray:
        assert self._beta is not None, "fit before score"
        Xb = np.concatenate([X_std, np.ones((X_std.shape[0], 1))], axis=1)
        return self._sigmoid(Xb @ self._beta)

    def _leave_one_out_auc(self, X: np.ndarray, y: np.ndarray) -> float:
        """Refit-each-row LOO probabilities, then AUC on them."""
        n = X.shape[0]
        if n < 4:
            return float("nan")
        loo_probs = np.zeros(n, dtype=np.float64)
        for i in range(n):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            Xi, yi = X[mask], y[mask]
            # Skip if the LOO fold lost a class entirely.
            if (yi == 1).sum() == 0 or (yi == 0).sum() == 0:
                loo_probs[i] = 0.5
                continue
            X_imp, med = self._impute(Xi)
            X_std, mu, sigma = self._standardise(X_imp)
            w = self._balanced_weights(yi)
            beta = self._fit_beta(X_std, yi, w)
            x_held = X[i : i + 1].copy()
            # Impute the held-out row using the LOO-fold median.
            x_imp_held, _ = self._impute(x_held, median=med)
            x_std_held = (x_imp_held - mu) / sigma
            xb = np.concatenate([x_std_held, np.ones((1, 1))], axis=1)
            loo_probs[i] = float(self._sigmoid(xb @ beta)[0])
        return roc_auc(y, loo_probs)

    def score(self, record: dict[str, Any]) -> float:
        """Probability this SILVER record resembles a past false-positive.

        NEVER raises. Returns:

        * a calibrated probability in ``[0, 1]`` when the flagger is fit and
          the record extracts cleanly;
        * ``0.5`` for an unfit flagger, an extractor failure, or any
          numerical degeneracy (the non-blocking contract).
        """
        if self._beta is None:
            return 0.5
        try:
            feats = extract_features(record if isinstance(record, dict) else {})
            x = feats.reshape(1, -1)
            x_imp, _ = self._impute(x, median=self._x_median)
            assert self._x_mean is not None and self._x_std is not None
            x_std = (x_imp - self._x_mean) / self._x_std
            xb = np.concatenate([x_std, np.ones((1, 1))], axis=1)
            p = float(self._sigmoid(xb @ self._beta)[0])
        except Exception:
            return 0.5
        if not np.isfinite(p):
            return 0.5
        # Clamp belt-and-braces.
        return float(min(max(p, 0.0), 1.0))

    # ------------------------------------------------------------------
    # Introspection.
    # ------------------------------------------------------------------
    @property
    def diagnostics(self) -> FitDiagnostics | None:
        return self._fit_diagnostics

    def is_fit(self) -> bool:
        return self._beta is not None


def _main() -> None:
    """CLI-ish entry: train on the labeled corpus and print diagnostics."""
    from cyclerfinder.ml.falsepos_labels import build_training_set

    X, y, _meta = build_training_set()
    clf = FalsePosFlagger()
    diag = clf.fit(X, y)
    print("=" * 72)
    print("FalsePosFlagger trained on labeled-corpus seed (#256)")
    print(
        f"  N={diag.n_samples}  positives(FP)={diag.n_positives}  negatives(TR)={diag.n_negatives}"
    )
    print(f"  AUC (train)  = {diag.auc_train:.4f}")
    print(f"  AUC (LOO)    = {diag.auc_loo:.4f}   [reproduce-before-trust gate: > 0.75]")
    print("  Feature importance (|coef| on standardised features):")
    for name, imp in diag.feature_importance:
        print(f"    {imp:7.4f}  {name}")
    print("=" * 72)


if __name__ == "__main__":  # pragma: no cover
    _main()


__all__ = [
    "FalsePosFlagger",
    "FitDiagnostics",
    "roc_auc",
]
