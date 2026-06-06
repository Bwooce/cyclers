"""Phase 5 — the human-review queue artifact (Forge plan §35).

The review queue holds SILVER candidates (machine-confirmed, unsourced) with a
full audit trail, pending human approval. GOLDEN DISCIPLINE: the queue is
explicitly NON-catalogue — the loader/validator must refuse to treat a queued
row as a catalogue entry (no row creation without human approval).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cyclerfinder.data.review_queue import (
    ReviewQueueEntry,
    append_review_entry,
    is_catalogue_source,
    load_review_queue,
    validate_review_entry,
)
from cyclerfinder.verify.gauntlet import VerdictTier


def _entry(tier: str = "silver") -> ReviewQueueEntry:
    return ReviewQueueEntry(
        candidate_id="novel|E-M-E-E|k2|r001|bssl|t974531381",
        signature_hash="sha1:" + "a" * 40,
        verdict_tier=tier,
        match_outcome="novel",
        known_id=None,
        superseded_by=(),
        vinf_per_encounter_kms=(9.75, 13.01, 9.76, 9.75),
        tof_days=(165.7, 564.7, 829.5),
        bend_feasible=True,
        max_vinf_kms=13.01,
        sequence=("E", "M", "E", "E"),
        period_k=2,
        model_assumption="analytic-ephemeris",
        verdict_audit={"A": {"agreed": True}},
        panel={"n_verifiers": 3, "n_refuted": 0, "majority_refute": False},
        t_added="2026-06-06T00:00:00+00:00",
    )


def test_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "review_queue.jsonl"
    e = _entry()
    append_review_entry(path, e)
    loaded = list(load_review_queue(path))
    assert len(loaded) == 1
    assert loaded[0] == e


def test_queue_is_not_a_catalogue_source() -> None:
    """The queue is NON-catalogue: no row creation without human approval."""
    assert is_catalogue_source() is False


def test_validate_rejects_non_silver_with_clean_panel() -> None:
    """Only SILVER (or GOLD-pending) candidates belong in the queue; a REJECTED
    candidate must never be queued for human promotion."""
    with pytest.raises(ValueError):
        validate_review_entry(_entry(tier="rejected"))


def test_validate_rejects_panel_majority_refute() -> None:
    """A candidate the adversarial panel majority-refuted must not be queued."""
    e = _entry()
    bad = ReviewQueueEntry(
        **{
            **e.__dict__,
            "panel": {"n_verifiers": 3, "n_refuted": 2, "majority_refute": True},
        }
    )
    with pytest.raises(ValueError):
        validate_review_entry(bad)


def test_validate_accepts_clean_silver() -> None:
    validate_review_entry(_entry(tier="silver"))  # no raise


def test_jsonl_is_one_object_per_line(tmp_path: Path) -> None:
    path = tmp_path / "q.jsonl"
    append_review_entry(path, _entry())
    append_review_entry(path, _entry())
    lines = path.read_text().strip().splitlines()
    assert len(lines) == 2
    for ln in lines:
        json.loads(ln)  # each line is a standalone JSON object


def test_silver_tier_enum_value_matches() -> None:
    """Guard against tier-string drift between the queue and the gauntlet enum."""
    assert _entry().verdict_tier == VerdictTier.SILVER.value
