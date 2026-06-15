"""Feature extraction for the false-positive flagger (#256).

The flagger reads a SILVER-record-like dict and emits a fixed-shape numeric
feature vector that captures the **signatures of our own past bugs**. The
feature set is intentionally small and hand-crafted from documented past
failures, not a generic catalogue dump:

* **#181** (ToF artifact, `cec9b90` / `3b09614`): Stage-B Lambert used the
  G-arc *branch* ToF instead of the row's tabulated *signature transit*,
  inflating emerged Mars V_inf ~1.6-2.1x. Signature: large residual at one
  body, V_inf well above its sourced floor, candidate produced before the fix.
* **#197** (`ba55b2e`, shared-cross-check non-independence): the verifier and
  the primary solver shared an upstream module, so "independent" agreement was
  bogus. Signature: ``cross_check_shared_with_primary_flag``.
* **#198** (`439d279`, 63 s UTC/TDB epoch offset): the J2000 TDB epoch was
  built in UTC seconds, drifting all sourced anchors by ~63 s. Signature:
  candidate produced before the fix-SHA.
* **#212a** (`3cec84c`, ``cr3bp_system`` Earth-Moon mu double-count): the
  CR3BP mu was scaled twice (~-1.2 %), producing wrong period/Jacobi pairs
  that LOOKED like good closures. Signature: candidate produced before the
  fix in a CR3BP closure record.
* **Bend-infeasible closure** (Braik-Ross scorer note): a candidate with a
  numerically-clean residual but a geometrically impossible flyby bend. Signature:
  ``bend_infeasible_flag`` true.
* **#249 period-impostor** (2026-06-14 plan note): a (k1,k2) that matched the
  period yet had the wrong winding -- caught by the topology classifier.
  Signature: ``topology_match_flag`` false despite small ``max_residual_kms``.
* **"0.163 km/s" scratch closure** (MEMORY): non-reproducible single-run
  result. Signature: ``closure_method_version`` is empty / a scratch tag.

All extractors are NaN-safe (missing key -> sentinel feature; never raises) so
the flagger upstream contract ("non-blocking; always returns a probability")
holds for partial / malformed records. The feature ORDER is the public ABI of
the trained model; appending to ``FEATURE_NAMES`` is safe, reordering is not.
"""

from __future__ import annotations

import math
from itertools import pairwise
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Known bug-fix commit SHAs. A candidate's ``closure_method_version`` (a git
# SHA, solver version tag, or ISO-date) is "before the fix" if it does not
# match a known good SHA and -- if it carries a date -- predates the fix date.
# These dates are the COMMIT dates of the fixes in the recent log; they are the
# epoch-artifact signal, not a precise calibration.
# ---------------------------------------------------------------------------
_FIX_SHA_TO_DATE: dict[str, str] = {
    "3b09614": "2026-06-10",  # #181 Stage-B ToF artifact fix
    "cec9b90": "2026-06-10",  # #181 ToF writeback
    "ba55b2e": "2026-06-10",  # #197 independent cross-check
    "439d279": "2026-06-10",  # #198 63s UTC/TDB epoch offset
    "3cec84c": "2026-06-10",  # #212a cr3bp_system mu double-count
    "4be2375": "2026-06-11",  # #212b Ross-Roberts-Tsoukkas adoption
    "23b980e": "2026-06-14",  # #249 asymmetric CR3BP corrector
    "a19eb24": "2026-06-14",  # #249 C21 recovery
}

# Scratch / non-reproducible solver tags. A closure produced under one of
# these is by construction unreviewable -- this is the "0.163 km/s scratch"
# signature from MEMORY.
_SCRATCH_VERSION_TAGS: frozenset[str] = frozenset({"", "scratch", "wip", "unknown", "none"})

# Known SILVER-record top-level keys the extractor will *read*. Recorded here
# so downstream callers can verify completeness; missing keys never raise.
SILVER_RECORD_KEYS: frozenset[str] = frozenset(
    {
        "max_residual_kms",
        "bend_feasible",
        "topology_match",
        "vinf_per_encounter_kms",
        "vinf_floors_kms",
        "period_days",
        "encounter_periods_days",
        "cross_check_shared_with_primary",
        "closure_method_version",
        "closure_date",
        "model_assumption",
    }
)

# Public feature ORDER -- this IS the model ABI. Append-only.
FEATURE_NAMES: tuple[str, ...] = (
    "max_residual_kms",
    "bend_infeasible_flag",
    "topology_match_flag",
    "vinf_max",
    "vinf_min",
    "vinf_floor_gap_kms",
    "period_resonance_deviation",
    "cross_check_shared_with_primary_flag",
    "epoch_artifact_flag",
    "scratch_solver_flag",
    "cr3bp_pre_mu_fix_flag",
    "residual_to_vinf_ratio",
)


def _to_float(value: Any, default: float = math.nan) -> float:
    """Coerce ``value`` to ``float`` without raising. Missing -> ``default``."""
    if value is None:
        return default
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def _to_float_seq(values: Any) -> list[float]:
    """Coerce a sequence to a list of finite floats (drops bad entries)."""
    if values is None:
        return []
    try:
        iterable = list(values)
    except TypeError:
        return []
    out: list[float] = []
    for v in iterable:
        f = _to_float(v)
        if not math.isnan(f):
            out.append(f)
    return out


def _bool_flag(value: Any, *, default: float = 0.0) -> float:
    """A 1.0 / 0.0 flag from a truthy / falsy field; missing -> ``default``."""
    if value is None:
        return default
    return 1.0 if bool(value) else 0.0


def _resonance_deviation(periods_days: list[float]) -> float:
    """Closeness of consecutive period ratios to a small-integer p:q.

    For each consecutive ratio ``r = P_i+1 / P_i`` we find the small-integer
    p:q (p,q in 1..6) that minimises ``|r - p/q|``, then return the MAX such
    distance across the sequence. A real resonant cycler should be near-zero;
    a "looks-like-a-cycler-but-isn't" closure typically drifts. NaN if fewer
    than two periods are provided.
    """
    if len(periods_days) < 2:
        return math.nan
    candidates = [p / q for p in range(1, 7) for q in range(1, 7)]
    worst = 0.0
    for a, b in pairwise(periods_days):
        if a <= 0 or b <= 0:
            return math.nan
        r = b / a
        if r <= 0:
            return math.nan
        # Normalise to >= 1 so 2:1 and 1:2 hit the same lattice point.
        r_norm = r if r >= 1.0 else 1.0 / r
        dev = min(abs(r_norm - c) for c in candidates if c >= 1.0)
        if dev > worst:
            worst = dev
    return worst


def _epoch_artifact_flag(version: str | None, date_iso: str | None) -> float:
    """Was the candidate produced BEFORE the relevant fix-SHAs?

    1.0 means the candidate predates at least one bug-fix SHA from
    ``_FIX_SHA_TO_DATE`` (matched by SHA prefix or by ISO date), and is
    therefore eligible to carry an inherited bug. 0.0 if it matches a known
    fix-SHA exactly. 0.5 if neither tag is present (genuinely unknown -- the
    flagger learns this is mildly suspicious).
    """
    if version is None and date_iso is None:
        return 0.5
    version = (version or "").strip().lower()
    # Matches a known fix-SHA prefix -> presumed POST-fix.
    for fix_sha in _FIX_SHA_TO_DATE:
        if version.startswith(fix_sha):
            return 0.0
    # Has a date that's >= the latest fix? Treat as POST-fix.
    if date_iso is not None:
        latest_fix = max(_FIX_SHA_TO_DATE.values())
        if str(date_iso) >= latest_fix:
            return 0.0
    # Has a date strictly before any fix -> PRE-fix.
    if date_iso is not None:
        earliest_fix = min(_FIX_SHA_TO_DATE.values())
        if str(date_iso) < earliest_fix:
            return 1.0
    # Tag present but unrecognised, no informative date -> mildly suspicious.
    return 0.5


def _scratch_solver_flag(version: str | None) -> float:
    """1.0 if the solver tag matches a known non-reproducible scratch marker."""
    if version is None:
        return 1.0
    tag = str(version).strip().lower()
    return 1.0 if tag in _SCRATCH_VERSION_TAGS else 0.0


def _cr3bp_pre_mu_fix_flag(model: str | None, version: str | None, date_iso: str | None) -> float:
    """The #212a signature: a CR3BP closure produced before the mu fix.

    1.0 iff the model_assumption mentions CR3BP / cr3bp AND the epoch-artifact
    flag fires (version not matching a post-fix SHA and date < fix date).
    """
    if model is None:
        return 0.0
    is_cr3bp = "cr3bp" in str(model).lower()
    if not is_cr3bp:
        return 0.0
    return _epoch_artifact_flag(version, date_iso)


def extract_features(record: dict[str, Any]) -> np.ndarray:
    """Map a SILVER record dict to the fixed-length feature vector.

    NEVER raises. Missing / malformed fields fall back to a deterministic
    default (median-like values for continuous fields, ``nan`` is replaced
    downstream by ``FalsePosFlagger`` with the training median).

    Parameters
    ----------
    record:
        A dict that *resembles* a ``ReviewQueueEntry`` payload but is treated
        as fully optional. See ``SILVER_RECORD_KEYS`` for the keys read.

    Returns
    -------
    numpy.ndarray
        Shape ``(len(FEATURE_NAMES),)``, dtype ``float64``. Entries may be
        ``nan`` -- the flagger imputes at score time.
    """
    rec = record if isinstance(record, dict) else {}

    max_residual = _to_float(rec.get("max_residual_kms"))

    bend_feasible_raw = rec.get("bend_feasible")
    if bend_feasible_raw is None:
        bend_infeasible_flag = 0.5  # unknown
    else:
        bend_infeasible_flag = 0.0 if bool(bend_feasible_raw) else 1.0

    topology_raw = rec.get("topology_match")
    topology_match_flag = 0.5 if topology_raw is None else (1.0 if bool(topology_raw) else 0.0)

    vinfs = _to_float_seq(rec.get("vinf_per_encounter_kms"))
    vinf_max = max(vinfs) if vinfs else math.nan
    vinf_min = min(vinfs) if vinfs else math.nan

    floors = _to_float_seq(rec.get("vinf_floors_kms"))
    if vinfs and floors and len(floors) == len(vinfs):
        gaps = [v - f for v, f in zip(vinfs, floors, strict=False)]
        vinf_floor_gap = max(gaps)
    else:
        vinf_floor_gap = math.nan

    encounter_periods = _to_float_seq(rec.get("encounter_periods_days"))
    if not encounter_periods:
        # Fall back to a single period if the per-encounter list isn't given.
        period_days = _to_float(rec.get("period_days"))
        encounter_periods = [period_days] if not math.isnan(period_days) else []
    period_resonance_dev = _resonance_deviation(encounter_periods)

    cross_check_shared_flag = _bool_flag(rec.get("cross_check_shared_with_primary"), default=0.5)

    version = rec.get("closure_method_version")
    date_iso = rec.get("closure_date")
    version_str = str(version) if version is not None else None
    date_str = str(date_iso) if date_iso is not None else None

    epoch_artifact = _epoch_artifact_flag(version_str, date_str)
    scratch_flag = _scratch_solver_flag(version_str)
    pre_mu_fix = _cr3bp_pre_mu_fix_flag(rec.get("model_assumption"), version_str, date_str)

    if not math.isnan(max_residual) and vinfs:
        denom = max(max(vinfs), 1e-3)
        residual_to_vinf = max_residual / denom
    else:
        residual_to_vinf = math.nan

    out = np.array(
        [
            max_residual,
            bend_infeasible_flag,
            topology_match_flag,
            vinf_max,
            vinf_min,
            vinf_floor_gap,
            period_resonance_dev,
            cross_check_shared_flag,
            epoch_artifact,
            scratch_flag,
            pre_mu_fix,
            residual_to_vinf,
        ],
        dtype=np.float64,
    )
    assert out.shape == (len(FEATURE_NAMES),), "feature vector ABI broken"
    return out


__all__ = [
    "FEATURE_NAMES",
    "SILVER_RECORD_KEYS",
    "extract_features",
]
