"""Task 4: validation-tier classifier tests (synthetic inputs).

``classify_validation`` is a pure, source-independent function, so it is
unit-testable over synthetic source/fidelity combinations without touching
the catalogue. The live-row census ratchet (plan Task 4 bullet 3) is
deferred with Task 3 (it needs the YAML source/fidelity back-fill).
"""

from __future__ import annotations

from cyclerfinder.data.provenance import Tier, classify_validation


def test_cross_validated_two_independent_same_fidelity() -> None:
    """Two different independent sources at the same fidelity -> CROSS_VALIDATED."""
    assert (
        classify_validation("rogers-2012-t1", "russell-2004-t34", same_fidelity=True)
        is Tier.CROSS_VALIDATED
    )


def test_consistency_checked_same_source() -> None:
    """Same source on both sides (same fidelity) -> CONSISTENCY_CHECKED."""
    assert (
        classify_validation("russell-2004-t34", "russell-2004-t34", same_fidelity=True)
        is Tier.CONSISTENCY_CHECKED
    )


def test_cross_fidelity_caps_at_unvalidated() -> None:
    """Two independent sources but DIFFERENT fidelity -> UNVALIDATED (S1L1 bug class)."""
    assert (
        classify_validation("rogers-2012-t1", "russell-2004-t34", same_fidelity=False)
        is Tier.UNVALIDATED
    )


def test_missing_orbit_source_unvalidated() -> None:
    assert classify_validation(None, "russell-2004-t34", same_fidelity=True) is Tier.UNVALIDATED


def test_missing_vinf_source_unvalidated() -> None:
    assert classify_validation("rogers-2012-t1", None, same_fidelity=True) is Tier.UNVALIDATED


def test_both_missing_unvalidated() -> None:
    assert classify_validation(None, None, same_fidelity=True) is Tier.UNVALIDATED


def test_pseudo_source_cannot_cross_validate() -> None:
    """A 'derived' value paired with a real source is NOT cross-validated."""
    assert (
        classify_validation("derived", "russell-2004-t34", same_fidelity=True) is Tier.UNVALIDATED
    )


def test_computed_pseudo_source_unvalidated() -> None:
    """'computed' (our own optimiser) can never validate -> UNVALIDATED."""
    assert classify_validation("computed", "computed", same_fidelity=True) is Tier.UNVALIDATED
