"""Forge phase 3 — verdict-tier census ratchet.

The Axis-C provenance-tier census lives in
``tests/data/test_validation_tier_census.py`` (it freezes how the live catalogue
rows distribute across ``classify_validation`` tiers). This is its phase-3
sibling: it freezes the *gauntlet verdict* surface — the four
:class:`~cyclerfinder.verify.gauntlet.VerdictTier` values, their ordinal
confidence mapping, and the canonical axis-input → tier decision matrix the
combiner emits.

Two surfaces are frozen here:

1. The combiner's *decision surface* — a fixed table of representative axis
   inputs and the tier each must produce. A change to the combination rules
   (e.g. accidentally letting a single source promote to GOLD, or letting a
   failing axis pass) breaks this ratchet, forcing a reviewed diff.
2. The *live gauntlet sweep* (#125 Part 3) — the verdict-tier distribution
   ``scripts/sweep_gauntlet_ledger.py`` writes to ``data/gauntlet_ledger.jsonl``
   over all 237 catalogue rows (Axis C for all; recorded Axis A only for the
   Aldrin pair). Frozen as a coverage ratchet over our own verdict bookkeeping.

Every verdict (matrix and live ledger) is additionally round-tripped through
:func:`~cyclerfinder.verify.gauntlet.validate_verdict` so the census cannot
contain a verdict the validator would reject.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from cyclerfinder.data.provenance import Corroboration, Tier
from cyclerfinder.verify.agreement import (
    AgreementReport,
    ConstructionOptimiserPathResult,
    KeplerRepropPathResult,
    LamberthubPathResult,
)
from cyclerfinder.verify.fidelity import PersistenceClass, PersistenceReport
from cyclerfinder.verify.gauntlet import (
    ValidationVerdict,
    VerdictTier,
    run_gauntlet,
    validate_verdict,
)

# ---------------------------------------------------------------------------
# Frozen confidence mapping (ordinal label per tier; NOT a probability).
# ---------------------------------------------------------------------------
EXPECTED_CONFIDENCE: dict[str, str] = {
    "gold": "high",
    "silver": "medium",
    "bronze": "low",
    "rejected": "none",
}


def _agreement(*, agreed: bool, n_available: int, n_passed: int) -> AgreementReport:
    empty_a = LamberthubPathResult(available=False, per_leg=(), max_diff_mps=0.0, passed=False)
    empty_b = ConstructionOptimiserPathResult(
        available=False,
        resonant_vinf_kms={},
        cycler_vinf_kms={},
        construction_max_diff_kms=float("inf"),
        optimiser_available=False,
        optimiser_vinf_kms={},
        optimiser_max_diff_kms=None,
        max_diff_kms=float("inf"),
        passed=False,
    )
    empty_c = KeplerRepropPathResult(
        available=False, per_leg_residual_km=(), max_residual_km=float("inf"), passed=False
    )
    return AgreementReport(
        lamberthub=empty_a,
        construction_optimiser=empty_b,
        kepler_reprop=empty_c,
        n_paths_available=n_available,
        n_paths_passed=n_passed,
        agreed=agreed,
    )


def _persist(cls: PersistenceClass) -> PersistenceReport:
    return PersistenceReport(
        classification=cls,
        quantity="q",
        low_value=0.0,
        high_value=0.0,
        delta=0.0,
        within_tol=cls is PersistenceClass.PERSISTS,
    )


# Canonical decision matrix: (label, kwargs) -> expected tier. Frozen ratchet.
def _matrix() -> list[tuple[str, dict, VerdictTier]]:  # type: ignore[type-arg]
    return [
        (
            "machine-confirmed + independent source",
            dict(
                agreement=_agreement(agreed=True, n_available=2, n_passed=2),
                persistence_reports=(_persist(PersistenceClass.PERSISTS),),
                provenance_tier=Tier.CROSS_VALIDATED,
                corroboration=Corroboration.STRONGLY_SOURCED,
            ),
            VerdictTier.GOLD,
        ),
        (
            "machine-confirmed, single-sourced (novel holding)",
            dict(
                agreement=_agreement(agreed=True, n_available=3, n_passed=3),
                provenance_tier=Tier.CONSISTENCY_CHECKED,
                corroboration=Corroboration.SINGLE_SOURCED,
            ),
            VerdictTier.SILVER,
        ),
        (
            "thin coverage (one path, no veto)",
            dict(
                agreement=_agreement(agreed=False, n_available=1, n_passed=1),
                corroboration=Corroboration.SINGLE_SOURCED,
            ),
            VerdictTier.BRONZE,
        ),
        (
            "nothing ran",
            dict(),
            VerdictTier.BRONZE,
        ),
        (
            "falsified",
            dict(
                agreement=_agreement(agreed=True, n_available=3, n_passed=3),
                provenance_tier=Tier.CROSS_VALIDATED,
                corroboration=Corroboration.STRONGLY_SOURCED,
                falsified=True,
            ),
            VerdictTier.REJECTED,
        ),
        (
            "available-but-failing Axis A (veto)",
            dict(
                agreement=_agreement(agreed=False, n_available=2, n_passed=1),
                provenance_tier=Tier.CROSS_VALIDATED,
                corroboration=Corroboration.STRONGLY_SOURCED,
            ),
            VerdictTier.REJECTED,
        ),
        (
            "undocumented fidelity shift (S1L1 class)",
            dict(
                agreement=_agreement(agreed=True, n_available=2, n_passed=2),
                persistence_reports=(_persist(PersistenceClass.SHIFTS_UNDOCUMENTED),),
                provenance_tier=Tier.CROSS_VALIDATED,
                corroboration=Corroboration.STRONGLY_SOURCED,
            ),
            VerdictTier.REJECTED,
        ),
        (
            "cross-fidelity dispute (S1L1 class)",
            dict(
                agreement=_agreement(agreed=True, n_available=2, n_passed=2),
                corroboration=Corroboration.DISPUTED,
                corroboration_cross_fidelity=True,
            ),
            VerdictTier.REJECTED,
        ),
    ]


# Frozen verdict-tier census over the canonical decision matrix.
EXPECTED_VERDICT_CENSUS: dict[str, int] = {
    "gold": 1,
    "silver": 1,
    "bronze": 2,
    "rejected": 4,
}


def test_confidence_mapping_frozen() -> None:
    """Every VerdictTier maps to its frozen ordinal confidence label."""
    # Drive the mapping through the public combiner for each canonical case.
    seen: dict[str, str] = {}
    for _label, kwargs, _expected in _matrix():
        v = run_gauntlet("probe", **kwargs)
        seen[v.tier.value] = v.confidence
    assert seen == {k: EXPECTED_CONFIDENCE[k] for k in seen}, seen
    assert set(EXPECTED_CONFIDENCE) == {t.value for t in VerdictTier}


def test_decision_matrix_matches_expected_tiers() -> None:
    """Each canonical axis-input combination yields its frozen tier (teeth)."""
    for label, kwargs, expected in _matrix():
        v = run_gauntlet("probe", **kwargs)
        assert v.tier is expected, f"{label}: got {v.tier.value}, expected {expected.value}"
        # And every emitted verdict survives the validator.
        validate_verdict(v)


def test_verdict_tier_census_frozen() -> None:
    """The verdict-tier distribution over the canonical matrix is frozen.

    Update ``EXPECTED_VERDICT_CENSUS`` in the same commit as any change to the
    combiner's decision rules or the matrix above.
    """
    counts = Counter(run_gauntlet("probe", **kwargs).tier.value for _l, kwargs, _e in _matrix())
    assert dict(counts) == EXPECTED_VERDICT_CENSUS, dict(counts)


# ---------------------------------------------------------------------------
# #125 Part 3 — LIVE gauntlet sweep over the 237 catalogue rows.
#
# scripts/sweep_gauntlet_ledger.py runs run_gauntlet over every catalogue row to
# the extent its data supports (Axis C provenance for all rows; recorded Axis A
# agreement only for the Aldrin pair) and writes the verdicts to
# data/gauntlet_ledger.jsonl. This freezes that live distribution.
#
# Honest verdicts (gauntlet rules, verify/gauntlet.py):
#   * inbound  — Axis A available + agreed (2/2) AND cross_validated -> GOLD.
#   * outbound — Axis A available-but-FAILING (3 paths, 2 pass; path b vetoes,
#     agreed=False) -> REJECTED. The multi-axis gauntlet is stricter than the
#     single-axis §14 V1 gate the row still legitimately holds (Part 1).
#   * every other row — Axis A unavailable, no failing axis -> BRONZE.
# ---------------------------------------------------------------------------
_LEDGER_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "gauntlet_ledger.jsonl"

# Frozen live-ledger verdict census (2026-06-06). Update in the same commit as
# any provenance-tag change, gauntlet-rule change, or Axis-A evidence change.
EXPECTED_LIVE_VERDICT_CENSUS: dict[str, int] = {
    "gold": 1,  # aldrin-classic-em-k1-inbound
    "rejected": 1,  # aldrin-classic-em-k1-outbound (path-b veto)
    "bronze": 235,
}


def _load_ledger_tiers() -> dict[str, str]:
    tiers: dict[str, str] = {}
    for raw in _LEDGER_PATH.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        rec = json.loads(line)
        tiers[rec["cell_id"]] = rec["verdict_tier"]
    return tiers


def test_live_gauntlet_ledger_census_frozen() -> None:
    """The verdict-tier distribution in the committed ledger matches the ratchet."""
    tiers = _load_ledger_tiers()
    counts = Counter(tiers.values())
    assert dict(counts) == EXPECTED_LIVE_VERDICT_CENSUS, dict(counts)
    assert sum(counts.values()) == 237


def test_live_gauntlet_aldrin_pair_verdicts() -> None:
    """The Aldrin pair carries the honest, asymmetric gauntlet verdicts."""
    tiers = _load_ledger_tiers()
    assert tiers["aldrin-classic-em-k1-inbound"] == "gold"
    assert tiers["aldrin-classic-em-k1-outbound"] == "rejected"


def test_live_gauntlet_ledger_verdicts_validate() -> None:
    """Every recorded verdict round-trips through validate_verdict (no over-claim)."""
    for raw in _LEDGER_PATH.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        rec = json.loads(line)
        audit = rec["verdict_audit"]
        tier = VerdictTier(rec["verdict_tier"])
        conf = {"gold": "high", "silver": "medium", "bronze": "low", "rejected": "none"}[tier.value]
        verdict = ValidationVerdict(
            tier=tier,
            confidence=conf,
            axis_results=audit["axis_results"],
            provenance=audit["provenance"],
        )
        validate_verdict(verdict)  # raises if the recorded axes don't entail the tier
