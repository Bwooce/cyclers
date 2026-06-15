"""Hand-curated labeled corpus seeds for the false-positive flagger (#256).

This file is the **labeled-corpus seed**, not a learned model on real run
data. Each entry is a dict shaped like a SILVER record (the keys
``falsepos_features.extract_features`` reads); fields marked ``_mocked=True``
have values inferred from the relevant bug-fix note rather than recovered
from a real run logfile. The hand-curation traces to documented sources --
each entry's ``_source`` field names the commit / docs note / catalogue row
the signature was extracted from.

DISCIPLINE: this is a SMALL N corpus (a few dozen rows). The trained flagger
is a guard rail, not a classifier we hand the keys to. ``FalsePosFlagger``
enforces ``class_weight='balanced'`` and reports cross-val-like AUC so the
caller can see when the corpus is too thin to trust.

Adding a new labeled example: append a dict with at least ``_label``
(``'false_positive'`` or ``'true_reproduction'``) and ``_source``, plus any of
the ``SILVER_RECORD_KEYS`` that are recoverable. Missing keys are imputed by
the extractor at NaN-safe defaults.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from cyclerfinder.ml.falsepos_features import extract_features

# ---------------------------------------------------------------------------
# NEGATIVE class -- closures we initially called wins, that turned out to be
# bugs. Each row's signature is traced to its bug-fix commit / docs note.
# ---------------------------------------------------------------------------
KNOWN_FALSE_POSITIVES: list[dict[str, Any]] = [
    # --- #181 ToF-artifact rows (Stage-B Lambert branch-ToF bug) ----------------
    # Pre-fix, the russell-ch4 G2 closure looked good but inflated Mars V_inf
    # ~1.6-2.1x. Signature: large max_residual_kms, vinf_max above sourced
    # floor, pre-fix solver tag.
    {
        "_label": "false_positive",
        "_source": "#181 / docs/notes/2026-06-10-dsm-tof-artifact-correction.md",
        "_mocked": True,
        "max_residual_kms": 2.4,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [6.4, 9.1],  # Mars inflated
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 1041.6],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-3b09614",
        "closure_date": "2026-06-08",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "false_positive",
        "_source": "#181 / Stage-B branch-ToF on russell-ch4-G_f2",
        "_mocked": True,
        "max_residual_kms": 1.9,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.8, 8.4],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 1090.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-3b09614",
        "closure_date": "2026-06-07",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "false_positive",
        "_source": "#181 / 6.44Gg3-like pre-fix Stage-B",
        "_mocked": True,
        "max_residual_kms": 3.1,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [7.2, 11.0],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 1100.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-3b09614",
        "closure_date": "2026-06-05",
        "model_assumption": "real_eph_DE440",
    },
    # --- #197 shared cross-check non-independence -------------------------------
    {
        "_label": "false_positive",
        "_source": "#197 / ba55b2e verify: independent Lambert endpoints",
        "_mocked": True,
        "max_residual_kms": 0.15,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.3, 7.1],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 780.0],
        "cross_check_shared_with_primary": True,  # the bug
        "closure_method_version": "pre-ba55b2e",
        "closure_date": "2026-06-08",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "false_positive",
        "_source": "#197 / shared upstream module",
        "_mocked": True,
        "max_residual_kms": 0.08,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.2, 6.9],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 800.0],
        "cross_check_shared_with_primary": True,
        "closure_method_version": "pre-ba55b2e",
        "closure_date": "2026-06-07",
        "model_assumption": "real_eph_DE440",
    },
    # --- #198 63s UTC/TDB epoch offset ------------------------------------------
    {
        "_label": "false_positive",
        "_source": "#198 / 439d279 J2000 TDB JD fix",
        "_mocked": True,
        "max_residual_kms": 0.32,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.4, 7.2],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 770.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-439d279",
        "closure_date": "2026-06-08",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "false_positive",
        "_source": "#198 / 63s offset before TDB fix",
        "_mocked": True,
        "max_residual_kms": 0.41,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.5, 7.4],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 750.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-439d279",
        "closure_date": "2026-06-07",
        "model_assumption": "real_eph_DE440",
    },
    # --- #212a cr3bp_system mu double-count (-1.2 %) ----------------------------
    {
        "_label": "false_positive",
        "_source": "#212a / 3cec84c mu double-count fix",
        "_mocked": True,
        "max_residual_kms": 0.02,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.6],  # CR3BP scaled
        "vinf_floors_kms": [0.8],
        "encounter_periods_days": [13.4, 27.2],  # off by ~1.2 %
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-3cec84c",
        "closure_date": "2026-06-06",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "false_positive",
        "_source": "#212a / RossEM Lyapunov pre-mu-fix",
        "_mocked": True,
        "max_residual_kms": 0.01,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.5],
        "vinf_floors_kms": [0.8],
        "encounter_periods_days": [13.2, 26.7],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "pre-3cec84c",
        "closure_date": "2026-06-05",
        "model_assumption": "CR3BP_EarthMoon",
    },
    # --- Bend-infeasible closure (Braik-Ross scorer note) -----------------------
    {
        "_label": "false_positive",
        "_source": "docs/notes/2026-06-13-braik-ross-reachable-set-scorer-results.md "
        "(evaluate_closure rejection rule)",
        "_mocked": True,
        "max_residual_kms": 0.04,
        "bend_feasible": False,  # the signature
        "topology_match": True,
        "vinf_per_encounter_kms": [5.5, 7.1],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 800.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "23b980e",
        "closure_date": "2026-06-14",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "false_positive",
        "_source": "bend-infeasible / impossible flyby bend at second body",
        "_mocked": True,
        "max_residual_kms": 0.06,
        "bend_feasible": False,
        "topology_match": True,
        "vinf_per_encounter_kms": [6.1, 8.4],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 820.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "23b980e",
        "closure_date": "2026-06-14",
        "model_assumption": "real_eph_DE440",
    },
    # --- #249 period-impostor (matched period, wrong winding) -------------------
    {
        "_label": "false_positive",
        "_source": "#249 / docs/notes/2026-06-14-249-unstable-member-recovery-plan.md",
        "_mocked": True,
        "max_residual_kms": 1e-6,  # numerically clean
        "bend_feasible": True,
        "topology_match": False,  # the signature
        "vinf_per_encounter_kms": [0.9],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [11.2],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "23b980e",
        "closure_date": "2026-06-14",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "false_positive",
        "_source": "#249 / period-only match before all-roots classifier",
        "_mocked": True,
        "max_residual_kms": 5e-7,
        "bend_feasible": True,
        "topology_match": False,
        "vinf_per_encounter_kms": [1.0],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [12.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "f608c6b",
        "closure_date": "2026-06-13",
        "model_assumption": "CR3BP_EarthMoon",
    },
    # --- "0.163 km/s" scratch closure -------------------------------------------
    {
        "_label": "false_positive",
        "_source": "MEMORY / non-reproducible scratch run",
        "_mocked": True,
        "max_residual_kms": 0.163,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [4.8, 6.7],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 820.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "scratch",
        "closure_date": None,
        "model_assumption": "FBS_analytic",
    },
    {
        "_label": "false_positive",
        "_source": "MEMORY / unreproducible single-run E-E-M-M trial",
        "_mocked": True,
        "max_residual_kms": 0.21,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.0, 7.0],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 800.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "wip",
        "closure_date": None,
        "model_assumption": "FBS_analytic",
    },
]


# ---------------------------------------------------------------------------
# POSITIVE class -- known-good reproductions (catalogue V1/V2/V3 rows, the
# 4 Braik-Ross #249 members, Liang CGE members #222, Ross 5-row EM #216).
# Numeric values trace to the catalogue / the named sourced anchors; we keep
# the corpus shape minimal -- residuals near zero, V_inf within floor, current
# fix-SHA. ``_mocked=True`` flags the SILVER-record-shape representation
# (we are not re-running each catalogue closure here; the FEATURES are the
# representation).
# ---------------------------------------------------------------------------
KNOWN_TRUE_REPRODUCTIONS: list[dict[str, Any]] = [
    # --- S1L1 russell-ch4 V3 rows (post-#181 fix) ------------------------------
    {
        "_label": "true_reproduction",
        "_source": "catalogue / russell-ch4-4.991gG2 V3 (cec9b90)",
        "_mocked": True,
        "max_residual_kms": 0.08,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [4.991, 5.55],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 779.96],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-10",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / russell-ch4-8.049gGf2 V3 (#188 powered)",
        "_mocked": True,
        "max_residual_kms": 0.06,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [8.049, 7.2],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 779.96],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-10",
        "model_assumption": "real_eph_DE440",
    },
    # --- Aldrin cycler outbound (V2) -------------------------------------------
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Aldrin cycler outbound V2",
        "_mocked": True,
        "max_residual_kms": 0.04,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [5.6, 8.6],
        "vinf_floors_kms": [4.0, 5.0],
        "encounter_periods_days": [365.25, 780.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-10",
        "model_assumption": "real_eph_DE440",
    },
    # --- Ross 5-row stable EM cyclers V2 (#216) --------------------------------
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Ross EM stable cycler family-1 V2 (#216)",
        "_mocked": True,
        "max_residual_kms": 1e-8,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.9],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.2, 26.4],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "4be2375",
        "closure_date": "2026-06-11",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Ross EM stable cycler family-2 V2 (#216)",
        "_mocked": True,
        "max_residual_kms": 2e-8,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.0],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.5, 27.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "4be2375",
        "closure_date": "2026-06-11",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Ross EM stable cycler family-3 V2 (#216)",
        "_mocked": True,
        "max_residual_kms": 3e-8,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.1],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [12.9, 25.8],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "4be2375",
        "closure_date": "2026-06-11",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Ross EM stable cycler family-4 V2 (#216)",
        "_mocked": True,
        "max_residual_kms": 1.5e-8,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.95],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.0, 26.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "4be2375",
        "closure_date": "2026-06-11",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Ross EM stable cycler family-5 V2 (#216)",
        "_mocked": True,
        "max_residual_kms": 2.5e-8,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.05],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.6, 27.2],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "4be2375",
        "closure_date": "2026-06-11",
        "model_assumption": "CR3BP_EarthMoon",
    },
    # --- Braik-Ross #249 unstable members C11a/C11b/C21/C32 --------------------
    {
        "_label": "true_reproduction",
        "_source": "#249 / Braik-Ross C11a V1 (a19eb24)",
        "_mocked": True,
        "max_residual_kms": 1e-7,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.92],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.3],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "a19eb24",
        "closure_date": "2026-06-14",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "#249 / Braik-Ross C11b V1 (a19eb24)",
        "_mocked": True,
        "max_residual_kms": 1e-7,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.94],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [13.5],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "a19eb24",
        "closure_date": "2026-06-14",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "#249 / Braik-Ross C21 V1 (a19eb24, asymmetric)",
        "_mocked": True,
        "max_residual_kms": 1e-7,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.97],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [84.533],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "a19eb24",
        "closure_date": "2026-06-14",
        "model_assumption": "CR3BP_EarthMoon",
    },
    {
        "_label": "true_reproduction",
        "_source": "#249 / Braik-Ross C32 V1 (a19eb24)",
        "_mocked": True,
        "max_residual_kms": 1e-7,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [1.02],
        "vinf_floors_kms": [0.5],
        "encounter_periods_days": [70.6],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "a19eb24",
        "closure_date": "2026-06-14",
        "model_assumption": "CR3BP_EarthMoon",
    },
    # --- Liang CGE V1 members (#222) -------------------------------------------
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Liang CGE member 1 V1 (#222)",
        "_mocked": True,
        "max_residual_kms": 0.02,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.88, 0.91],
        "vinf_floors_kms": [0.5, 0.5],
        "encounter_periods_days": [13.5, 27.0],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-12",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Liang CGE member 2 V1 (#222)",
        "_mocked": True,
        "max_residual_kms": 0.03,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.90, 0.93],
        "vinf_floors_kms": [0.5, 0.5],
        "encounter_periods_days": [13.4, 26.8],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-12",
        "model_assumption": "real_eph_DE440",
    },
    {
        "_label": "true_reproduction",
        "_source": "catalogue / Liang CGE member 3 V1 (#222)",
        "_mocked": True,
        "max_residual_kms": 0.025,
        "bend_feasible": True,
        "topology_match": True,
        "vinf_per_encounter_kms": [0.89, 0.92],
        "vinf_floors_kms": [0.5, 0.5],
        "encounter_periods_days": [13.6, 27.2],
        "cross_check_shared_with_primary": False,
        "closure_method_version": "cec9b90",
        "closure_date": "2026-06-12",
        "model_assumption": "real_eph_DE440",
    },
]


def build_training_set() -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Return ``(X, y, meta)`` for the labeled corpus.

    ``X`` is an ``(N, len(FEATURE_NAMES))`` feature matrix (may contain
    NaN -- the flagger imputes), ``y`` is ``(N,)`` int 1 = false-positive
    / 0 = true-reproduction (the flagger predicts ``P(false_positive)``),
    ``meta`` is the original list of dicts in row order for provenance.
    """
    rows: list[dict[str, Any]] = []
    labels: list[int] = []
    for r in KNOWN_FALSE_POSITIVES:
        rows.append(r)
        labels.append(1)
    for r in KNOWN_TRUE_REPRODUCTIONS:
        rows.append(r)
        labels.append(0)
    X = np.stack([extract_features(r) for r in rows], axis=0)
    y = np.array(labels, dtype=np.int64)
    return X, y, rows


__all__ = [
    "KNOWN_FALSE_POSITIVES",
    "KNOWN_TRUE_REPRODUCTIONS",
    "build_training_set",
]
