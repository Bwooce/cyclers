"""Phase 3 — the gauntlet combiner: Axes A-D → :class:`ValidationVerdict`.

Unit tests for :func:`cyclerfinder.verify.gauntlet.run_gauntlet`, the pure
combiner that folds the four axis reports into a tiered verdict. These tests
are fast (no physics): they synthesise the axis reports and assert the tier
decision, the supersession-aware provenance, and the "teeth" (a failing
available axis demotes; a fabricated/invalid verdict is rejected by the
validator).
"""

from __future__ import annotations

import pytest

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
# Synthetic axis-report builders (no physics; shape-only)
# ---------------------------------------------------------------------------


def _agreement(*, agreed: bool, n_available: int, n_passed: int | None = None) -> AgreementReport:
    """Build a shape-only AgreementReport with the combiner-relevant fields."""
    if n_passed is None:
        n_passed = n_available if agreed else max(n_available - 1, 0)
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


def _persistence(cls: PersistenceClass) -> PersistenceReport:
    return PersistenceReport(
        classification=cls,
        quantity="outbound_tof_days",
        low_value=146.0,
        high_value=165.0,
        delta=19.0,
        within_tol=False,
    )


# ---------------------------------------------------------------------------
# Tier-decision matrix
# ---------------------------------------------------------------------------


def test_gold_requires_independent_source() -> None:
    """A and B clean, D ok, AND an independent source ⇒ GOLD."""
    verdict = run_gauntlet(
        "cand-gold",
        agreement=_agreement(agreed=True, n_available=2),
        persistence_reports=(_persistence(PersistenceClass.PERSISTS),),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
    )
    assert verdict.tier is VerdictTier.GOLD
    assert verdict.confidence == "high"


def test_silver_when_machine_confirmed_but_unsourced() -> None:
    """Machine-confirmed (A) but NO independent source ⇒ SILVER (capped)."""
    verdict = run_gauntlet(
        "cand-silver",
        agreement=_agreement(agreed=True, n_available=2),
        persistence_reports=(_persistence(PersistenceClass.SHIFTS_DOCUMENTED),),
        provenance_tier=Tier.CONSISTENCY_CHECKED,
        corroboration=Corroboration.SINGLE_SOURCED,
    )
    assert verdict.tier is VerdictTier.SILVER
    assert verdict.confidence == "medium"


def test_silver_not_promoted_to_gold_by_single_source() -> None:
    """SINGLE_SOURCED never counts as the independent source GOLD needs."""
    verdict = run_gauntlet(
        "cand-novel",
        agreement=_agreement(agreed=True, n_available=3),
        provenance_tier=Tier.CONSISTENCY_CHECKED,
        corroboration=Corroboration.SINGLE_SOURCED,
    )
    assert verdict.tier is VerdictTier.SILVER


def test_bronze_when_axis_a_unavailable_not_failing() -> None:
    """Fewer than two paths ran (unavailable, not failing) ⇒ BRONZE."""
    verdict = run_gauntlet(
        "cand-bronze",
        agreement=_agreement(agreed=False, n_available=1, n_passed=1),
        provenance_tier=Tier.UNVALIDATED,
        corroboration=Corroboration.SINGLE_SOURCED,
    )
    assert verdict.tier is VerdictTier.BRONZE
    assert verdict.confidence == "low"


def test_bronze_when_no_axes_run() -> None:
    """Nothing ran at all (no agreement report) ⇒ BRONZE, not REJECTED."""
    verdict = run_gauntlet("cand-empty")
    assert verdict.tier is VerdictTier.BRONZE


def test_rejected_on_falsification() -> None:
    """Falsification dominates: REJECTED even with everything else clean."""
    verdict = run_gauntlet(
        "cand-bogus",
        agreement=_agreement(agreed=True, n_available=3),
        persistence_reports=(_persistence(PersistenceClass.PERSISTS),),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
        falsified=True,
    )
    assert verdict.tier is VerdictTier.REJECTED
    assert verdict.confidence == "none"


def test_failing_available_axis_a_demotes_to_rejected() -> None:
    """An available-but-failing Axis A path vetoes ⇒ REJECTED (teeth)."""
    verdict = run_gauntlet(
        "cand-axisfail",
        # Two paths ran, did NOT agree → available-but-failing.
        agreement=_agreement(agreed=False, n_available=2, n_passed=1),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
    )
    assert verdict.tier is VerdictTier.REJECTED


def test_undocumented_shift_demotes_to_rejected() -> None:
    """A SHIFTS_UNDOCUMENTED persistence report (S1L1 class) ⇒ REJECTED."""
    verdict = run_gauntlet(
        "cand-shift",
        agreement=_agreement(agreed=True, n_available=2),
        persistence_reports=(_persistence(PersistenceClass.SHIFTS_UNDOCUMENTED),),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
    )
    assert verdict.tier is VerdictTier.REJECTED


def test_disputed_cross_fidelity_demotes_to_rejected() -> None:
    """A cross-fidelity DISPUTED corroboration (the S1L1 bug class) ⇒ REJECTED."""
    verdict = run_gauntlet(
        "cand-disputed",
        agreement=_agreement(agreed=True, n_available=2),
        provenance_tier=Tier.UNVALIDATED,
        corroboration=Corroboration.DISPUTED,
        corroboration_cross_fidelity=True,
    )
    assert verdict.tier is VerdictTier.REJECTED


def test_disputed_single_fidelity_does_not_reject() -> None:
    """A pure single-fidelity numeric dispute is surfaced, not auto-killed."""
    verdict = run_gauntlet(
        "cand-disputed-1f",
        agreement=_agreement(agreed=True, n_available=2),
        provenance_tier=Tier.UNVALIDATED,
        corroboration=Corroboration.DISPUTED,
        corroboration_cross_fidelity=False,
    )
    # Not rejected; machine-confirmed but no independent source ⇒ SILVER.
    assert verdict.tier is VerdictTier.SILVER


# ---------------------------------------------------------------------------
# Supersession-aware matching (Forge R1 delta 3)
# ---------------------------------------------------------------------------


def test_supersession_chain_in_provenance() -> None:
    """A verdict touching a superseded row carries the chain, not a clean known."""
    verdict = run_gauntlet(
        "cand-super",
        agreement=_agreement(agreed=True, n_available=2),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
        known_id="vem-emeeve-3syn",
        superseded_by=("jones-2017-vem-emevve-outbound",),
    )
    assert verdict.provenance["superseded_by"] == ["jones-2017-vem-emevve-outbound"]
    assert verdict.provenance["match_status"] == "superseded"
    # The tier is NOT downgraded on supersession alone (the superseding rows
    # may themselves be strong); the invalidated premise is just made visible.
    assert verdict.tier is VerdictTier.GOLD


def test_clean_known_match_when_not_superseded() -> None:
    verdict = run_gauntlet(
        "cand-known",
        agreement=_agreement(agreed=True, n_available=2),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
        known_id="aldrin-classic-em-k1-outbound",
    )
    assert verdict.provenance["match_status"] == "known"
    assert verdict.provenance["superseded_by"] == []


def test_unmatched_when_no_known_id() -> None:
    verdict = run_gauntlet("cand-unmatched", agreement=_agreement(agreed=True, n_available=2))
    assert verdict.provenance["match_status"] == "unmatched"


# ---------------------------------------------------------------------------
# Verdict validator — teeth: a fabricated/inconsistent verdict is rejected
# ---------------------------------------------------------------------------


def test_validate_accepts_a_real_verdict() -> None:
    verdict = run_gauntlet(
        "cand-ok",
        agreement=_agreement(agreed=True, n_available=2),
        provenance_tier=Tier.CROSS_VALIDATED,
        corroboration=Corroboration.STRONGLY_SOURCED,
    )
    validate_verdict(verdict)  # must not raise


def test_validate_rejects_gold_without_independent_source() -> None:
    """A hand-built GOLD verdict whose axis_results lack an independent source
    is internally inconsistent and must be refused (the fabrication guard)."""
    fake = ValidationVerdict(
        tier=VerdictTier.GOLD,
        confidence="high",
        axis_results={
            "A": {"available": True, "agreed": True},
            "B": {"ran": False, "undocumented_shift": False},
            "C": {"has_independent_source": False, "disputed_cross_fidelity": False},
            "D": {"falsified": False},
        },
        provenance={"candidate_id": "x", "superseded_by": []},
    )
    with pytest.raises(ValueError, match="GOLD requires an independent source"):
        validate_verdict(fake)


def test_validate_rejects_confidence_tier_mismatch() -> None:
    fake = ValidationVerdict(
        tier=VerdictTier.BRONZE,
        confidence="high",  # wrong: BRONZE is "low"
        axis_results={
            "A": {"available": False, "agreed": None},
            "B": {"ran": False, "undocumented_shift": False},
            "C": {"has_independent_source": False, "disputed_cross_fidelity": False},
            "D": {"falsified": False},
        },
        provenance={"candidate_id": "x", "superseded_by": []},
    )
    with pytest.raises(ValueError, match="confidence"):
        validate_verdict(fake)


def test_validate_rejects_rejected_tier_with_clean_axes() -> None:
    """A REJECTED verdict must have a refuting cause in its axis_results."""
    fake = ValidationVerdict(
        tier=VerdictTier.REJECTED,
        confidence="none",
        axis_results={
            "A": {"available": True, "agreed": True},
            "B": {"ran": True, "undocumented_shift": False},
            "C": {"has_independent_source": True, "disputed_cross_fidelity": False},
            "D": {"falsified": False},
        },
        provenance={"candidate_id": "x", "superseded_by": []},
    )
    with pytest.raises(ValueError, match="REJECTED requires"):
        validate_verdict(fake)
