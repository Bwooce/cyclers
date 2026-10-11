"""#1056: the closing-repeat flag on CandidateSignature.

Pre-registration: docs/notes/2026-10-11-1056-closing-repeat-flag.md. The probe is gc6-0's shape: a
one-cycle list (Ganymede, Callisto, Ganymede) whose cycle has two consecutive Ganymede encounters.
"""

from __future__ import annotations

from typing import Any

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    _has_consecutive_same_body,
    check_literature,
    offline_corpus_search,
)

G, C = "Ganymede", "Callisto"


def _sig(seq: tuple[str, ...], closing_repeat: bool | None) -> CandidateSignature:
    kw: dict[str, Any] = {
        "primary": "Jupiter",
        "sequence": seq,
        "period_k": None,
        "vinf_per_encounter_kms": (),
        "topology_label": frozenset({"repeated-moon"}),
        "working_bodies": "two",
        "closing_repeat": closing_repeat,
    }
    return CandidateSignature(**kw)


def test_flag_semantics() -> None:
    assert not _has_consecutive_same_body((G, C, G))  # legacy: closing repeat inferred
    assert _has_consecutive_same_body((G, C, G), closing_repeat=False)  # one cycle: G-G wrap
    assert not _has_consecutive_same_body((G, C, G, C, G), closing_repeat=True)


def test_one_cycle_gcg_is_not_alternating() -> None:
    r = check_literature(_sig((G, C, G), False), search=offline_corpus_search)
    assert r.status != "published", r
    assert "Campagnola" not in str(r.citation), r


def test_explicit_closing_repeat_still_matches_gcgc() -> None:
    r = check_literature(_sig((G, C, G, C, G), True), search=offline_corpus_search)
    assert r.status == "published" and "Campagnola" in str(r.citation), r
