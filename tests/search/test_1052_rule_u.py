"""#1052 rule U, U2 and the Aldrin single-body fix (lead ruling on gate semantics, 2026-10-11).

Pre-registration: docs/notes/2026-10-11-1052-rule-u-untopologied-signatures.md (with addenda A1-A2).
"""

from __future__ import annotations

from typing import Any

import pytest

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)


def _check(**kw: Any) -> Any:
    base: dict[str, Any] = {
        "period_k": None,
        "vinf_per_encounter_kms": (),
        "topology_label": frozenset(),
    }
    base.update(kw)
    return check_literature(CandidateSignature(**base), search=offline_corpus_search)


@pytest.mark.parametrize(
    ("kw", "author"),
    [
        ({"primary": "Sun", "sequence": ("E", "M")}, "Byrnes"),  # Aldrin
        (
            {"primary": "Earth", "sequence": ("Moon",), "resonances": ("7:3",),
             "working_bodies": "one"},
            "Casoliva",
        ),
        (
            {"primary": "Jupiter", "sequence": ("Ganymede", "Callisto", "Ganymede"),
             "working_bodies": "one"},
            "Russell",
        ),
        (
            {"primary": "Saturn", "sequence": ("Titan", "Enceladus", "Titan"),
             "working_bodies": "one"},
            "Russell",
        ),
    ],
)  # fmt: skip
def test_positives_via_repeated_moon_anchors_stay_published(
    kw: dict[str, Any], author: str
) -> None:
    r = _check(**kw)
    assert r.status == "published" and str(r.citation).startswith(author), r


@pytest.mark.parametrize(
    ("primary", "seq"),
    [
        ("Saturn", ("Titan",)),
        ("Neptune", ("Triton",)),
        ("Mars", ("Phobos",)),
        ("Sun", ("Jupiter",)),
        ("Uranus", ("Umbriel", "Titania")),
        ("Pluto", ("Charon", "Nix")),
        ("Sun", ("E",)),
        ("Sun", ("M",)),
        ("Sun", ("E", "E")),
    ],
)
def test_untopologied_non_cycler_hits_capped(primary: str, seq: tuple[str, ...]) -> None:
    assert _check(primary=primary, sequence=seq).status == "inconclusive"


def test_aldrin_single_body_topologied_goes_to_a_human_not_fresh() -> None:
    r = _check(primary="Sun", sequence=("M",), topology_label=frozenset({"repeated-moon"}))
    assert r.status == "inconclusive", r


def test_two_body_ev_not_escalated() -> None:
    r = _check(
        primary="Sun",
        sequence=("V", "V", "V", "E"),
        working_bodies="two",
        return_types=frozenset({"FR", "GEN"}),
        topology_label=frozenset({"repeated-moon"}),
    )
    assert r.status == "not-found", r


def test_u2_unlabelled_anchor_cannot_publish_a_labelled_signature() -> None:
    r = _check(primary="Pluto", sequence=("Charon", "Nix"), topology_label=frozenset({"pump-tour"}))
    assert r.status != "published", r


def test_jupiter_negative_unchanged() -> None:
    r = _check(
        primary="Jupiter",
        sequence=("Io", "Callisto"),
        working_bodies="one",
        topology_label=frozenset({"repeated-moon"}),
    )
    assert r.status == "known-architecture-new-system", r
