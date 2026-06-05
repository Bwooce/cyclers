"""Task 5: multi-source corroboration scoring.

Exercises :func:`cyclerfinder.data.provenance.score_corroboration`, which
consumes the Task-3 provenance vocabulary (source keys + fidelity tiers) to
classify a single quantity as ``strongly-sourced`` / ``single-sourced`` /
``disputed`` and to *record* (never hide) the agreement spread.

The S1L1 V∞ case is the plan's motivating example: spec.md §9 anchors the
Earth V∞ at 5.65 km/s while the literature sources reachable for that cycler
give 4.99 km/s (Russell 2004 Table 4.9, the 4.991gG2/S1L1 parent cycler) and
4.7 km/s (McConaghy 2006 abstract). All three are values DOCUMENTED in the
catalogue prose / spec — this test consumes those documented numbers, it does
not invent any. The spread (~0.95 km/s) exceeds the V∞ corroboration tolerance,
so the quantity is correctly classified DISPUTED, documenting the known
fidelity-mismatch rather than averaging it away.
"""

from __future__ import annotations

from cyclerfinder.data.provenance import (
    Corroboration,
    SourcedValue,
    score_corroboration,
)


def test_strongly_sourced_when_two_sources_agree() -> None:
    """Two distinct independent sources within tolerance -> strongly-sourced."""
    score = score_corroboration(
        [
            SourcedValue(4.99, "russell-2004-t49_413", "circular-coplanar"),
            SourcedValue(5.10, "mcconaghy-2002", "circular-coplanar"),
        ]
    )
    assert score.classification is Corroboration.STRONGLY_SOURCED
    assert score.independent_source_count == 2
    assert score.spread < 0.5
    assert score.cross_fidelity is False


def test_single_sourced_when_one_source() -> None:
    score = score_corroboration([SourcedValue(5.65, "spec-9", "circular-coplanar")])
    assert score.classification is Corroboration.SINGLE_SOURCED
    assert score.independent_source_count == 1
    assert score.spread == 0.0


def test_single_sourced_when_same_source_repeated() -> None:
    """Two values from the SAME source are not corroboration."""
    score = score_corroboration(
        [
            SourcedValue(4.99, "russell-2004-t49_413", "circular-coplanar"),
            SourcedValue(5.10, "russell-2004-t49_413", "circular-coplanar"),
        ]
    )
    assert score.classification is Corroboration.SINGLE_SOURCED
    assert score.independent_source_count == 1


def test_pseudo_sources_do_not_corroborate() -> None:
    """derived/computed never count as independent corroboration."""
    score = score_corroboration(
        [
            SourcedValue(4.99, "russell-2004-t49_413", "circular-coplanar"),
            SourcedValue(4.99, "derived", "circular-coplanar"),
            SourcedValue(4.99, "computed", "circular-coplanar"),
        ]
    )
    assert score.classification is Corroboration.SINGLE_SOURCED
    assert score.independent_source_count == 1


def test_s1l1_earth_vinf_is_disputed() -> None:
    """The known S1L1 Earth-V∞ disagreement is classified DISPUTED, with the
    spread surfaced (plan Task 5: document, don't hide, the fidelity mismatch)."""
    score = score_corroboration(
        [
            SourcedValue(5.65, "spec-9", "circular-coplanar"),
            SourcedValue(4.99, "russell-2004-t49_413", "circular-coplanar"),
            SourcedValue(4.7, "mcconaghy-2006", "circular-coplanar"),
        ]
    )
    assert score.classification is Corroboration.DISPUTED
    assert score.independent_source_count == 3
    assert score.spread > 0.5  # 5.65 - 4.7 = 0.95


def test_cross_fidelity_dispute_flagged() -> None:
    """A dispute spanning fidelities sets cross_fidelity (the S1L1 bug class:
    a coplanar value compared against an analytic-ephemeris one)."""
    score = score_corroboration(
        [
            SourcedValue(5.65, "spec-9", "circular-coplanar"),
            SourcedValue(4.99, "russell-2004-t49_413", "analytic-ephemeris"),
        ]
    )
    assert score.classification is Corroboration.DISPUTED
    assert score.cross_fidelity is True


def test_agreement_at_one_fidelity_not_cross_fidelity() -> None:
    score = score_corroboration(
        [
            SourcedValue(4.99, "russell-2004-t49_413", "circular-coplanar"),
            SourcedValue(5.10, "mcconaghy-2002", "circular-coplanar"),
        ]
    )
    assert score.cross_fidelity is False
