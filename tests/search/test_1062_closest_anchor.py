"""#1062: the literature gate cites the closest in-footprint anchor, not the first in corpus order.

Pre-registration: docs/notes/2026-10-11-1062-closest-anchor-citation.md sec. 2. The expected
citations follow the sourced resonances_scope tags (Casoliva 2010 Table 3; Braik-Ross 2026 Table 2).
"""

from __future__ import annotations

from typing import Any

import pytest

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)


def _sig(resonances: tuple[str, ...], period_k: int | None = 3) -> CandidateSignature:
    kw: dict[str, Any] = {
        "primary": "Earth",
        "sequence": ("Moon",),
        "period_k": period_k,
        "vinf_per_encounter_kms": (),
        "resonances": resonances,
        "topology_label": frozenset({"resonant"}),
        "working_bodies": "one",
    }
    return CandidateSignature(**kw)


@pytest.mark.parametrize(
    ("reso", "author"),
    [
        (("7:3",), "Casoliva"),  # the #1047 Casoliva 7-3b control
        (("3:1",), "Braik"),
        (("2:1",), "Braik"),  # both declare 2:1; corpus order
        (("4:3",), "Kumar"),  # declared by neither; the first undeclared anchor
        ((), "Braik"),  # no resonances: unchanged from v3-A1
    ],
)
def test_1062_cites_closest_anchor(reso: tuple[str, ...], author: str) -> None:
    r = check_literature(_sig(reso), search=offline_corpus_search)
    assert r.status == "published", r
    assert str(r.citation).startswith(author), (reso, r.citation)
