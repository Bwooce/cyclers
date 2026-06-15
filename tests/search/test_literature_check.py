"""Tests for the literature-novelty check (#261).

The check is only trustworthy if it reliably flags KNOWN-PUBLISHED cyclers as
``published``. The SELF-VALIDATION block here feeds it the signatures of (a) a
Liang CGE member, (b) an Aldrin Earth-Mars cycler, and (c) a Russell/McConaghy
SnLm row, against a deterministic fake corpus mirroring the *real* published
search hits (titles/authors/dois as they appear in the literature), and asserts
each returns ``status="published"`` with a plausible citation. A deliberately-
fabricated nonsense signature must return ``not-found``.

The fake corpus is deterministic so CI is reproducible; the *live* WebSearch
self-validation (real queries + verdicts) is recorded in the build note
``docs/notes/2026-06-14-literature-novelty-check.md`` per task discipline.
"""

from __future__ import annotations

from collections.abc import Sequence

import pytest

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    LiteratureCheckResult,
    SearchResult,
    build_queries,
    check_literature,
    is_novelty_claimable,
)

# ---------------------------------------------------------------------------
# A deterministic fake corpus mirroring the REAL published search hits.
# Each entry's title/snippet carries the publication's actual authors + venue
# so the structural matcher scores it exactly as it would a live WebSearch hit.
# ---------------------------------------------------------------------------

_REAL_CORPUS: list[SearchResult] = [
    SearchResult(
        title="The Aldrin Earth-Mars cycler trajectory",
        url="https://doi.org/10.2514/6.2002-4420",
        snippet="Byrnes, McConaghy and Longuski analyse the Aldrin cycler, a "
        "ballistic Earth-Mars cycler trajectory using gravity assists.",
    ),
    SearchResult(
        title="Ballistic Earth-Mars cycler trajectories (systematic catalog)",
        url="https://doi.org/10.2514/1.10078",
        snippet="Russell and Ocampo present a systematic catalog of ballistic "
        "Earth-Mars cycler trajectories including two-synodic cyclers; "
        "McConaghy, Landau, Longuski extend the SnLm family.",
    ),
    SearchResult(
        title="Callisto-Ganymede-Europa Triple Cyclers",
        url="https://doi.org/10.2514/1.G008387",
        snippet="Liang, Yang, Bai and Qin present Callisto-Ganymede-Europa (CGE) "
        "triple cycler trajectories in the Jovian moon system, a ballistic "
        "moon cycler.",
    ),
    SearchResult(
        title="One Class of Io-Europa-Ganymede Triple Cyclers",
        url="https://example.org/aas-2017-ieg",
        snippet="Hernandez, Jones and Jesick describe Io-Europa-Ganymede triple "
        "cycler trajectories exploiting the Laplace resonance (Jovian moon cycler).",
    ),
    # Some noise hits that must NOT be mistaken for a structural match.
    SearchResult(
        title="A review of bicycle gear ratios",
        url="https://example.org/bikes",
        snippet="Cyclists and gear ratios for road bikes.",
    ),
    SearchResult(
        title="Lunar gateway station keeping",
        url="https://example.org/gateway",
        snippet="NRHO station keeping for the lunar Gateway, no cycler content.",
    ),
]


def fake_search(query: str) -> Sequence[SearchResult]:
    """Return corpus hits whose text overlaps the query tokens (like a search).

    A crude term-overlap ranker over the fixed corpus: any corpus row sharing
    >=2 lowercased alphabetic tokens with the query is "returned". Deterministic,
    no network. Mirrors how a real WebSearch would surface the published hit for
    a structurally-specific query.
    """
    q_terms = {t for t in _tokenise(query) if len(t) > 2}
    out: list[tuple[int, SearchResult]] = []
    for r in _REAL_CORPUS:
        text_terms = set(_tokenise(r.title + " " + r.snippet))
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def _tokenise(s: str) -> list[str]:
    return [t for t in "".join(c.lower() if c.isalnum() else " " for c in s).split()]


# ---------------------------------------------------------------------------
# Known-published self-validation: each MUST come back "published".
# ---------------------------------------------------------------------------

ALDRIN_SIG = CandidateSignature(
    primary="Sun",
    sequence=("E", "M"),
    period_k=1,
    period_years=2.135,
    vinf_per_encounter_kms=(6.5, 9.7),
)

RUSSELL_SIG = CandidateSignature(
    primary="Sun",
    sequence=("E", "E", "M", "M"),
    period_k=2,
    period_years=4.27,
    vinf_per_encounter_kms=(4.1, 2.0),
)

LIANG_SIG = CandidateSignature(
    primary="Jupiter",
    sequence=("Callisto", "Ganymede", "Callisto", "Europa", "Callisto"),
    vinf_per_encounter_kms=(5.673, 6.992, 4.668),
    n_rev=(1, 1, 1, 1),
)


@pytest.mark.parametrize(
    ("sig", "expect_token"),
    [
        (ALDRIN_SIG, "aldrin"),
        (RUSSELL_SIG, "russell"),
        (LIANG_SIG, "liang"),
    ],
)
def test_known_published_flagged(sig: CandidateSignature, expect_token: str) -> None:
    result = check_literature(sig, search=fake_search)
    assert result.status == "published", f"known-published cycler not flagged: {result}"
    assert result.citation, "published verdict must carry a citation"
    # The citation OR the matched url must reference the right publication line.
    blob = (result.citation + " " + (result.matched_url or "")).lower()
    assert expect_token in blob or result.doi, (
        f"citation does not point at the expected publication: {result.citation!r}"
    )
    assert result.confidence >= 0.70
    assert not is_novelty_claimable(result.to_review_block())


def test_aldrin_doi_extracted() -> None:
    result = check_literature(ALDRIN_SIG, search=fake_search)
    assert result.doi is not None and result.doi.startswith("10.2514")


def test_liang_is_moon_tour_published() -> None:
    result = check_literature(LIANG_SIG, search=fake_search)
    assert result.status == "published"
    assert "1.g008387" in (result.doi or "").lower() or "liang" in (result.citation or "").lower()


# ---------------------------------------------------------------------------
# Fabricated nonsense signature: MUST be "not-found".
# ---------------------------------------------------------------------------


def test_fabricated_signature_not_found() -> None:
    # A nonsense primary + invented bodies that no published cycler uses.
    bogus = CandidateSignature(
        primary="Neptune",
        sequence=("Triton", "Nereid", "Triton", "Proteus"),
        vinf_per_encounter_kms=(3.3, 7.1, 2.2),
        n_rev=(2, 1, 0, 1),
    )
    result = check_literature(bogus, search=fake_search)
    assert result.status == "not-found", f"fabricated sig should not match: {result}"
    assert result.citation is None
    assert result.doi is None
    # A clean not-found IS novelty-claimable -- it passed the rediscovery filter
    # (necessary-not-sufficient; the human + gauntlet still decide).
    assert is_novelty_claimable(result.to_review_block())
    # The query trail must show real searches were attempted.
    assert len(result.query_trail) >= 3


def test_fabricated_with_empty_search_is_inconclusive() -> None:
    """A search that returns nothing at all => inconclusive, NOT not-found.

    We must never emit a clean not-found when no search actually ran.
    """

    def empty_search(_query: str) -> Sequence[SearchResult]:
        return []

    bogus = CandidateSignature(primary="Sun", sequence=("E", "M"), period_k=1)
    result = check_literature(bogus, search=empty_search)
    assert result.status == "inconclusive"
    assert not is_novelty_claimable(result.to_review_block())


# ---------------------------------------------------------------------------
# The novelty gate.
# ---------------------------------------------------------------------------


def test_gate_rejects_unpopulated() -> None:
    assert is_novelty_claimable(None) is False
    assert is_novelty_claimable({}) is False
    assert is_novelty_claimable({"checked": False}) is False


def test_gate_rejects_published() -> None:
    block = LiteratureCheckResult(
        status="published", citation="X", doi=None, confidence=0.9
    ).to_review_block()
    assert is_novelty_claimable(block) is False


def test_gate_rejects_inconclusive() -> None:
    block = LiteratureCheckResult(
        status="inconclusive", citation=None, doi=None, confidence=0.5
    ).to_review_block()
    assert is_novelty_claimable(block) is False


def test_gate_allows_clean_not_found() -> None:
    block = LiteratureCheckResult(
        status="not-found", citation=None, doi=None, confidence=0.1
    ).to_review_block()
    assert is_novelty_claimable(block) is True
    # And the legacy is_promotion_eligible mapping agrees (result == no-match).
    assert block["result"] == "no-match"


# ---------------------------------------------------------------------------
# Query construction sanity.
# ---------------------------------------------------------------------------


def test_build_queries_specificity_order() -> None:
    qs = build_queries(ALDRIN_SIG)
    assert qs, "must build at least one query"
    # Most-specific (tour + period) before the generic fallback.
    assert any("Earth" in q and "Mars" in q for q in qs)
    assert qs == list(dict.fromkeys(qs)), "queries must be de-duplicated"
    # Named-corpus anchor surfaced for an E-M heliocentric candidate.
    assert any("Aldrin" in q for q in qs)


def test_build_queries_moon_tour_branch() -> None:
    qs = build_queries(LIANG_SIG)
    assert any("moon tour" in q for q in qs)
    assert any("Callisto" in q for q in qs)


# ---------------------------------------------------------------------------
# #272 corpus expansion: Pluto-Charon + recent capability papers.
#
# Each new anchor in KNOWN_CORPUS should be findable via a synthetic search
# result whose title/snippet carries the published authors + keywords. We
# feed three independent signatures (one Pluto-system, one Earth-Moon cycler
# family, one tulip-orbit) and verify check_literature flags each "published"
# with a citation that names the right line.
# ---------------------------------------------------------------------------


_NEW_CORPUS_HITS: list[SearchResult] = [
    SearchResult(
        title="Persephone: A Pluto-system Orbiter and Kuiper Belt Explorer",
        url="https://doi.org/10.3847/PSJ/abf837",
        snippet="Howard, Stern et al. design Persephone, a Pluto orbiter with "
        "Pluto-Charon CR3BP periodic orbits and Nix Hydra encounters; "
        "binary rotating frame cycler science orbits.",
    ),
    SearchResult(
        title="Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics",
        url="https://vsgc.odu.edu/wp-content/uploads/2026/04/Roberts-Tsoukkas_Michael_Cycler-Journal-Paper.pdf",
        snippet="Roberts-Tsoukkas and Ross present stable prograde Earth-Moon "
        "cycler families across mass parameters including a universal stable "
        "subfamily; multi-orbiter cycler trajectories in the three-body problem.",
    ),
    SearchResult(
        title="Novel Tulip-Shaped Three-body Orbits for Cislunar SDA Missions",
        url="https://doi.org/10.1007/s40295-025-00510-w",
        snippet="Koblick and Kelly construct tulip-shaped three-body cycler "
        "orbits in the Earth-Moon CR3BP for cislunar SDA missions; "
        "petal-count periodic orbit families.",
    ),
]


def _new_corpus_search(query: str) -> Sequence[SearchResult]:
    """Deterministic ranker over the #272 expansion corpus (and the old one)."""
    q_terms = {t for t in _tokenise(query) if len(t) > 2}
    out: list[tuple[int, SearchResult]] = []
    for r in (*_REAL_CORPUS, *_NEW_CORPUS_HITS):
        text_terms = set(_tokenise(r.title + " " + r.snippet))
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


PLUTO_PERSEPHONE_SIG = CandidateSignature(
    primary="Pluto",
    sequence=("Charon", "Nix", "Hydra", "Charon"),
    period_k=2,
    vinf_per_encounter_kms=(0.18, 0.07, 0.05),
)

ROBERTS_TSOUKKAS_SIG = CandidateSignature(
    primary="Earth",
    sequence=("Moon",),
    vinf_per_encounter_kms=(0.6,),
)

KOBLICK_TULIP_SIG = CandidateSignature(
    primary="Earth",
    sequence=("Moon",),
    vinf_per_encounter_kms=(0.4,),
    resonances=("tulip",),
)


# Tokens that any of the #272 new anchors could legitimately surface in a
# citation/url/doi. Per signature class, the matcher's structural overlap +
# author/keyword override should return at least ONE of these -- proving the
# expansion is reachable through check_literature. The matcher picks the
# first overlap in KNOWN_CORPUS insertion order; with multiple Earth-Moon
# anchors that all overlap a generic E-M cycler signature, ANY of them is a
# valid "new corpus" hit -- the point is the matcher resolves to a #272 line,
# not the pre-existing Aldrin/Russell/Liang/Strange/Jones/Hernandez entries.
_NEW_CORPUS_TOKENS: dict[str, tuple[str, ...]] = {
    "pluto-system": (
        "persephone",
        "showalter",
        "brozovic",
        "stern",
        "pluto",
        "10.3847/psj/abf837",
        "10.1038/nature14469",
    ),
    # Any of the Earth-Moon CR3BP #272 anchors (Braik-Ross orbital networks;
    # Roberts-Tsoukkas multi-orbiter; Kumar resonant transport; Koblick tulip;
    # Zhang tulip; Hiraiwa lobe-dynamics; Chinese J. tulip; Sanaga fidelity).
    "earth-moon-new": (
        "roberts-tsoukkas",
        "braik",
        "kumar",
        "rosengren",
        "koblick",
        "kelly",
        "tulip",
        "hiraiwa",
        "lobe dynamics",
        "orbital networks",
        "multi-orbiter",
        "cislunar resonant transport",
        "2605.31543",
        "2509.12675",
        "2602.17444",
        "10.1007/s40295-025-00510-w",
        "10.1007/s11071-026-12465-0",
        # Roberts-Tsoukkas journal manuscript URL host:
        "vsgc.odu.edu",
    ),
}


@pytest.mark.parametrize(
    ("sig", "token_class"),
    [
        (PLUTO_PERSEPHONE_SIG, "pluto-system"),
        (ROBERTS_TSOUKKAS_SIG, "earth-moon-new"),
        (KOBLICK_TULIP_SIG, "earth-moon-new"),
    ],
)
def test_new_corpus_entries_flagged_published(sig: CandidateSignature, token_class: str) -> None:
    """At least 3 of the #272 additions must be findable as ``published``.

    Mirrors the existing Aldrin/Russell/Liang self-validation: a synthetic
    search result with the publication's real authors + venue must score above
    MATCH_THRESHOLD and surface a citation/url referencing the right line.

    Multiple corpus anchors structurally overlap an Earth-Moon candidate
    (Braik-Ross, Roberts-Tsoukkas, Kumar, Koblick, Hiraiwa, ...); the matcher
    picks the first author/keyword overlap in KNOWN_CORPUS insertion order.
    The test passes if the verdict's citation OR matched_url surfaces ANY
    token from the expected publication class -- the point is that the new
    anchors are reachable through ``check_literature``, not that any one of
    them is uniquely picked from a multi-anchor structural overlap.
    """
    result = check_literature(sig, search=_new_corpus_search)
    assert result.status == "published", (
        f"new-corpus entry {token_class!r} not flagged as published: {result}"
    )
    assert result.citation, "published verdict must carry a citation"
    blob = (result.citation + " " + (result.matched_url or "") + " " + (result.doi or "")).lower()
    expected = _NEW_CORPUS_TOKENS[token_class]
    assert any(tok in blob for tok in expected), (
        f"verdict does not surface any expected token from {token_class!r}: "
        f"citation={result.citation!r} matched_url={result.matched_url!r} "
        f"doi={result.doi!r}"
    )
    assert result.confidence >= 0.70
    assert not is_novelty_claimable(result.to_review_block())


def test_new_corpus_anchors_registered() -> None:
    """Direct registration check for the #272 KNOWN_CORPUS additions.

    Independent of the live-matcher path, the corpus expansion must register
    the named publication anchors so downstream code (the discovery daemon's
    candidate-anchor walker, the offline citation override) can see them. We
    verify the canonical primaries + authors / keywords are present.
    """
    from cyclerfinder.search.literature_check import KNOWN_CORPUS

    pluto_anchors = [a for a in KNOWN_CORPUS if a.primary == "Pluto"]
    assert len(pluto_anchors) >= 3, (
        f"expected >=3 Pluto-system anchors after #272 expansion; got {len(pluto_anchors)}"
    )
    pluto_authors = {a for anchor in pluto_anchors for a in anchor.authors}
    assert "Howard" in pluto_authors  # Persephone
    assert "Showalter" in pluto_authors  # Styx-Nix-Hydra resonance

    earth_moon_anchors = [a for a in KNOWN_CORPUS if a.primary == "Earth" and "Moon" in a.body_set]
    em_authors = {a for anchor in earth_moon_anchors for a in anchor.authors}
    # Recent capability papers (Track A).
    assert "Braik" in em_authors  # Braik-Ross orbital networks
    assert "Roberts-Tsoukkas" in em_authors  # Roberts-Tsoukkas multi-orbiter
    assert "Kumar" in em_authors  # Kumar-Rawat-Rosengren-Ross
    assert "Koblick" in em_authors  # Koblick tulip
    assert "Hiraiwa" in em_authors  # Hiraiwa lobe dynamics


def test_review_entry_signature_roundtrip() -> None:
    from cyclerfinder.search.literature_check import signature_from_review_entry

    class _FakeEntry:
        sequence = ("Callisto", "Ganymede", "Europa")
        vinf_per_encounter_kms = (5.6, 6.9, 4.6)
        period_k = 4

        def __init__(self) -> None:
            self.verdict_audit = {"primary": "Jupiter", "n_rev": [1, 1, 1, 1]}

    sig = signature_from_review_entry(_FakeEntry())
    assert sig.primary == "Jupiter"
    assert sig.sequence == ("Callisto", "Ganymede", "Europa")
    assert sig.n_rev == (1, 1, 1, 1)
    assert sig.is_moon_tour
