"""Literature-novelty check for discovery-daemon SILVER candidates (#261).

The #253 discovery daemon's "novel" verdict only deduplicates a candidate
against our 277-row sourced catalogue plus the method-versioned negative
registry -- a *subset* of the published cycler literature. A daemon SILVER could
therefore be a cycler published *elsewhere* (a rediscovery mislabelled novel).
This module is the missing filter: given a candidate's defining structural
SIGNATURE, it searches the published record (WebSearch / WebFetch over arXiv,
NASA ADS / NTRS, AIAA, and the named cycler corpus) and returns a structured
verdict that populates the review-queue ``literature_check`` field.

DISCIPLINE (read before trusting any output):

* **"not-found" is NECESSARY-NOT-SUFFICIENT for novelty.** Absence of a search
  hit is *not* evidence of absence in the literature -- web search is partial,
  paywalled venues and conference proceedings are under-indexed, and structural
  fingerprints do not always survive into searchable text. This check is a
  FILTER against obvious rediscoveries, not a novelty proof. The human reviewer
  and the V0-V5 gauntlet still govern; a "not-found" only *clears* a candidate
  to continue, it never *certifies* it novel.
* **Cite by publication, never by repo path.** Citations returned here point at
  DOIs / arXiv ids / public URLs only.
* **No catalogue writeback.** This populates the review-queue
  ``literature_check`` block only (a human promotes).

The check is only trustworthy if it reliably flags KNOWN-PUBLISHED cyclers as
"published". The self-validation in ``tests/search/test_literature_check.py``
feeds it the signatures of a Liang CGE member, an Aldrin Earth-Mars cycler, and
a Russell/McConaghy SnLm row and asserts each returns ``status="published"``
with a plausible citation; a fabricated nonsense signature must return
"not-found". If those reproductions fail, the check is too weak to trust and
``not-found`` must not be believed.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

# A WebSearch callable: query -> list of {"title": str, "url": str} result dicts.
# Injectable so the engine drives the live tool while tests pin a deterministic
# fake corpus. The default raises -- a caller MUST wire the real WebSearch (the
# daemon does) so the module never silently "passes" with no search performed.
SearchFn = Callable[[str], Sequence["SearchResult"]]
FetchFn = Callable[[str, str], str]

Status = Literal["published", "not-found", "inconclusive"]


@dataclass(frozen=True)
class SearchResult:
    """One web-search hit (title + URL), the minimal shape the matcher reads."""

    title: str
    url: str
    snippet: str = ""


@dataclass(frozen=True)
class CandidateSignature:
    """The defining structural fingerprint of a candidate cycler.

    This is the literature-search key: the *physics identity* of the trajectory,
    not its print-precision numbers. Two trajectories with the same primary, the
    same body/moon tour (rotation-canonical), the same per-leg resonance and
    revolution structure, and the same V_inf regime are the same cycler family
    regardless of which paper printed which trailing digit.
    """

    primary: str
    """Central body (``Sun`` for heliocentric, ``Jupiter`` for a moon tour)."""

    sequence: tuple[str, ...]
    """Body / moon encounter tour (e.g. ``("E", "M")`` or the Galilean moons)."""

    period_k: int | None = None
    """Cycle period in synodic / repeat counts (``k`` in the SnLm naming)."""

    period_years: float | None = None
    """Cycle repeat period in years (heliocentric) -- ``None`` for moon tours."""

    vinf_per_encounter_kms: tuple[float, ...] = ()
    """V_inf magnitude at each encounter (km/s); the regime fingerprint."""

    resonances: tuple[str, ...] = ()
    """Per-leg p:q resonance labels if known (e.g. ``("3:2", "4:3")``)."""

    n_rev: tuple[int, ...] = ()
    """Per-leg full-revolution counts (the multi-rev structure)."""

    @property
    def is_moon_tour(self) -> bool:
        """A non-solar primary => a planetary-satellite (moon-tour) cycler."""
        return self.primary not in ("Sun", "", None)


@dataclass(frozen=True)
class LiteratureCheckResult:
    """Structured verdict for one candidate's literature search.

    ``status``:
      * ``published``  -- a published cycler matches the structural fingerprint
        (a rediscovery; NOT novelty-claimable).
      * ``not-found``  -- no published match surfaced (NECESSARY-not-sufficient
        for novelty; the candidate may continue to the human + gauntlet).
      * ``inconclusive`` -- the search itself could not be trusted (no results
        at all / search error); treat as NOT novelty-claimable until rerun.
    """

    status: Status
    citation: str | None
    doi: str | None
    confidence: float
    query_trail: list[str] = field(default_factory=list)
    matched_url: str | None = None
    notes: str = ""

    def to_review_block(self, *, reviewer: str = "literature_check.py (#261)") -> dict[str, Any]:
        """Render as the review-queue ``literature_check`` block (machine-written).

        Shape matches :class:`~cyclerfinder.data.review_queue.ReviewQueueEntry`'s
        documented contract: ``{checked, reviewer, date, sources_searched,
        result}`` plus the structured detail. ``result`` is mapped so the
        existing :func:`~cyclerfinder.data.review_queue.is_promotion_eligible`
        gate (``result == "no-match"``) only fires on a clean not-found.
        """
        from datetime import UTC, datetime

        result = {
            "published": "match",
            "not-found": "no-match",
            "inconclusive": "inconclusive",
        }[self.status]
        return {
            "checked": True,
            "reviewer": reviewer,
            "date": datetime.now(UTC).isoformat(),
            "sources_searched": list(self.query_trail),
            "result": result,
            "status": self.status,
            "citation": self.citation,
            "doi": self.doi,
            "confidence": self.confidence,
            "matched_url": self.matched_url,
            "notes": self.notes,
            "discipline": (
                "not-found is NECESSARY-NOT-SUFFICIENT for novelty: absence of a "
                "search hit is not evidence of absence; the human + V0-V5 "
                "gauntlet still govern. This is a rediscovery filter, not a "
                "novelty proof."
            ),
        }


# ---------------------------------------------------------------------------
# Named cycler corpus: structural fingerprints -> search anchors
# ---------------------------------------------------------------------------

_BODY_LONG = {
    "E": "Earth",
    "M": "Mars",
    "V": "Venus",
    "S": "Saturn",
    "J": "Jupiter",
    "Me": "Mercury",
}


def _expand_body(token: str) -> str:
    return _BODY_LONG.get(token, token)


@dataclass(frozen=True)
class CorpusAnchor:
    """A known published cycler line: how to recognise + cite it.

    ``primary`` + ``body_set`` (the set of tour members) + an optional
    ``period_k`` band define the structural footprint; ``citation`` / ``doi``
    are what we return when a search confirms a hit. ``authors`` and ``keywords``
    drive the structural-fingerprint queries.
    """

    name: str
    primary: str
    body_set: frozenset[str]
    authors: tuple[str, ...]
    keywords: tuple[str, ...]
    citation: str
    doi: str | None
    domains: tuple[str, ...] = ()


# Hand-curated from the catalogue's published rows + the task's named corpus.
# These are PUBLICATION facts (author/venue/doi), not values our code computed.
KNOWN_CORPUS: tuple[CorpusAnchor, ...] = (
    CorpusAnchor(
        name="Aldrin Earth-Mars cycler",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        authors=("Aldrin", "Byrnes", "McConaghy", "Longuski"),
        keywords=("Aldrin cycler", "Earth-Mars cycler", "cycler trajectory"),
        citation="Byrnes, McConaghy & Longuski, AIAA 2002-4420 (Aldrin cycler)",
        doi="10.2514/6.2002-4420",
    ),
    CorpusAnchor(
        name="Russell-Ocampo / McConaghy Earth-Mars SnLm cyclers",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        authors=("Russell", "Ocampo", "McConaghy", "Landau", "Longuski"),
        keywords=(
            "ballistic Earth-Mars cycler",
            "S1L1 cycler",
            "two-synodic cycler",
            "systematic cycler catalog",
        ),
        citation="Russell & Ocampo, J. Spacecraft & Rockets 41(1) 2004; "
        "McConaghy et al., J. Spacecraft & Rockets 43(2) 2006",
        doi="10.2514/1.10078",
    ),
    CorpusAnchor(
        name="Liang et al. Callisto-Ganymede-Europa triple cyclers",
        primary="Jupiter",
        body_set=frozenset({"Callisto", "Ganymede", "Europa"}),
        authors=("Liang", "Yang", "Bai", "Qin"),
        keywords=(
            "Callisto-Ganymede-Europa triple cycler",
            "CGE triple cycler",
            "Jovian moon cycler",
        ),
        citation="Liang, Yang, Li, Bai & Qin, JGCD Engineering Note 2024 "
        "(Callisto-Ganymede-Europa Triple Cyclers)",
        doi="10.2514/1.G008387",
    ),
    CorpusAnchor(
        name="Hernandez/Jones/Jesick Io-Europa-Ganymede triple cyclers",
        primary="Jupiter",
        body_set=frozenset({"Io", "Europa", "Ganymede"}),
        authors=("Hernandez", "Jones", "Jesick"),
        keywords=("Io-Europa-Ganymede triple cycler", "Jovian triple cycler"),
        citation="Hernandez, Jones & Jesick, 'One Class of Io-Europa-Ganymede "
        "Triple Cyclers' (AAS/AIAA 2017)",
        doi=None,
    ),
    CorpusAnchor(
        name="Strange/Campagnola/Russell moon-tour & V-infinity-leveraging",
        primary="Jupiter",
        body_set=frozenset({"Io", "Europa", "Ganymede", "Callisto"}),
        authors=("Strange", "Campagnola", "Russell", "Landau"),
        keywords=(
            "V-infinity leveraging",
            "Tisserand graph moon tour",
            "ballistic cycler Jovian system",
        ),
        citation="Strange, Russell & Buffington, 'Mapping the V-infinity globe' "
        "(AAS 07-277); Campagnola & Russell, JGCD 2010 (V-inf leveraging)",
        doi="10.2514/1.45645",
    ),
    CorpusAnchor(
        name="Jones et al. VEM triple cyclers (Venus-Earth-Mars)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M"}),
        authors=("Jones", "Hernandez", "Jesick"),
        keywords=("VEM triple cycler", "Venus-Earth-Mars cycler"),
        citation="Jones, Hernandez & Jesick, AAS 17-577 (VEM triple cyclers)",
        doi=None,
    ),
)


def _candidate_anchors(sig: CandidateSignature) -> list[CorpusAnchor]:
    """Corpus anchors whose structural footprint overlaps the signature.

    Match on PRIMARY (same dynamical system) and a non-trivial body-set overlap
    -- the structural fingerprint, not a keyword. A heliocentric Earth-Mars
    candidate cannot collide with a Jovian moon anchor and vice versa.
    """
    seq_set = frozenset(sig.sequence)
    anchors: list[CorpusAnchor] = []
    for anchor in KNOWN_CORPUS:
        if anchor.primary != sig.primary:
            continue
        overlap = seq_set & anchor.body_set
        # Require the candidate's tour to be a subset-ish of the anchor's body
        # set (every encountered body is one the anchor's family visits), which
        # is the structural-fingerprint test, not a single shared body.
        if overlap and seq_set <= anchor.body_set:
            anchors.append(anchor)
    return anchors


# ---------------------------------------------------------------------------
# Query construction (structural fingerprint -> search strings)
# ---------------------------------------------------------------------------


def build_queries(sig: CandidateSignature) -> list[str]:
    """Build the ordered query trail for a signature (most-specific first).

    Strategy: combine the structural fingerprint (primary + body/moon tour +
    resonance/n_rev + V_inf regime) with the named-corpus author/keyword anchors,
    then fall back to generic cycler queries over arXiv / ADS / NTRS. The order
    is most-specific -> most-generic so an early high-confidence hit short-
    circuits the search but a weak signature still casts a wide net.
    """
    bodies = [_expand_body(b) for b in sig.sequence]
    tour = "-".join(bodies)
    uniq_tour = "-".join(dict.fromkeys(bodies))  # dedup, order-preserving
    queries: list[str] = []

    # 1. Most specific: exact tour + period + cycler.
    if sig.period_k:
        queries.append(f"{uniq_tour} cycler {sig.period_k} synodic trajectory")
    queries.append(f"{tour} cycler trajectory")
    queries.append(f"{uniq_tour} cycler")

    # 2. Resonance / multi-rev structure (the SnLm / multi-rev fingerprint).
    if sig.resonances:
        reso = " ".join(sig.resonances)
        queries.append(f"{uniq_tour} cycler resonance {reso}")
    if any(nr > 0 for nr in sig.n_rev):
        queries.append(f"{uniq_tour} multi-revolution cycler trajectory")

    # 3. V_inf regime (leveraging / ballistic).
    if sig.vinf_per_encounter_kms:
        vmax = max(sig.vinf_per_encounter_kms)
        regime = "high V-infinity" if vmax > 6.0 else "ballistic low V-infinity"
        queries.append(f"{uniq_tour} {regime} cycler")

    # 4. Named-corpus author / keyword anchors for overlapping families.
    for anchor in _candidate_anchors(sig):
        for kw in anchor.keywords[:2]:
            queries.append(kw)
        queries.append(f"{anchor.authors[0]} {uniq_tour} cycler")

    # 5. Generic literature-index fallback.
    if sig.is_moon_tour:
        queries.append(f"{sig.primary} moon tour cycler ballistic patched-conic")
    else:
        queries.append(f"{sig.primary}-centered cycler trajectory periodic")

    # De-dup, preserve order.
    seen: set[str] = set()
    out: list[str] = []
    for q in queries:
        if q not in seen:
            seen.add(q)
            out.append(q)
    return out


# ---------------------------------------------------------------------------
# Result matching (structural fingerprint vs search-hit text)
# ---------------------------------------------------------------------------

_DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:a-z0-9]+", re.IGNORECASE)
_ARXIV_RE = re.compile(r"arxiv\.org/abs/(\d{4}\.\d{4,5})", re.IGNORECASE)


def _extract_doi(url: str, text: str) -> str | None:
    for blob in (url, text):
        m = _DOI_RE.search(blob)
        if m:
            return m.group(0).rstrip(".")
    return None


def _result_matches_fingerprint(sig: CandidateSignature, result: SearchResult) -> float:
    """Confidence (0-1) that a search hit is the candidate's cycler family.

    Scores on STRUCTURAL overlap with the hit's title+snippet: the word
    "cycler" (mandatory floor), then each tour body named, then a corpus
    author / keyword, then the primary. A hit that says "cycler" and names the
    full tour and a corpus author is a near-certain rediscovery.
    """
    hay = (result.title + " " + result.snippet).lower()
    if "cycler" not in hay and "cyclic" not in hay:
        return 0.0

    score = 0.30  # has "cycler"
    bodies = [_expand_body(b).lower() for b in sig.sequence]
    named = sum(1 for b in set(bodies) if b in hay)
    if bodies:
        score += 0.35 * (named / len(set(bodies)))

    if sig.primary.lower() in hay or (
        sig.is_moon_tour and ("moon" in hay or "satellite" in hay or "jovian" in hay)
    ):
        score += 0.10

    for anchor in _candidate_anchors(sig):
        if any(a.lower() in hay for a in anchor.authors):
            score += 0.20
            break
        if any(kw.lower() in hay for kw in anchor.keywords):
            score += 0.15
            break

    return min(score, 1.0)


# ---------------------------------------------------------------------------
# The check
# ---------------------------------------------------------------------------

MATCH_THRESHOLD: float = 0.70
"""Confidence at/above which a hit is treated as a published rediscovery."""

INCONCLUSIVE_FLOOR: float = 0.45
"""Best confidence below MATCH_THRESHOLD but >= this => 'inconclusive', not a
clean not-found: the search found cycler-adjacent material it could not rule out
as the same family, so a human must look (we do not certify novelty on it)."""

MAX_QUERIES: int = 8
"""Cap the live-search fan-out per candidate (cost + politeness)."""


def check_literature(
    sig: CandidateSignature,
    *,
    search: SearchFn,
    max_queries: int = MAX_QUERIES,
) -> LiteratureCheckResult:
    """Search the published record for a cycler matching ``sig``.

    Drives ``search`` (the injected WebSearch) over the query trail, scores each
    hit's structural overlap with the signature, and returns the best verdict:

    * any hit at/above :data:`MATCH_THRESHOLD` -> ``published`` (+ citation/doi),
    * a best hit in ``[INCONCLUSIVE_FLOOR, MATCH_THRESHOLD)`` -> ``inconclusive``,
    * no cycler-adjacent hit at all, with searches that DID return results
      -> ``not-found`` (necessary-not-sufficient for novelty),
    * a search that returned nothing/erroring everywhere -> ``inconclusive``
      (we cannot trust a not-found we never actually searched for).
    """
    queries = build_queries(sig)[:max_queries]
    trail: list[str] = []
    best_conf = 0.0
    best_hit: SearchResult | None = None
    any_results = False
    anchors = _candidate_anchors(sig)

    for q in queries:
        trail.append(q)
        try:
            results = list(search(q))
        except Exception as exc:  # a flaky search is inconclusive, not novel
            trail[-1] = f"{q}  [ERROR: {exc!r}]"
            continue
        if results:
            any_results = True
        for r in results:
            conf = _result_matches_fingerprint(sig, r)
            if conf > best_conf:
                best_conf = conf
                best_hit = r
        if best_conf >= MATCH_THRESHOLD:
            break  # short-circuit on a confident hit

    if best_hit is not None and best_conf >= MATCH_THRESHOLD:
        doi = _extract_doi(best_hit.url, best_hit.snippet)
        citation = best_hit.title
        # Prefer a curated corpus citation/doi when an overlapping anchor exists
        # and the hit corroborates it (named author/keyword present in the hit).
        hay = (best_hit.title + " " + best_hit.snippet).lower()
        for anchor in anchors:
            if any(a.lower() in hay for a in anchor.authors) or any(
                kw.lower() in hay for kw in anchor.keywords
            ):
                citation = anchor.citation
                doi = doi or anchor.doi
                break
        return LiteratureCheckResult(
            status="published",
            citation=citation,
            doi=doi,
            confidence=round(best_conf, 3),
            query_trail=trail,
            matched_url=best_hit.url,
            notes="Structural fingerprint matched a published cycler -- treat as "
            "a rediscovery; NOT novelty-claimable.",
        )

    if not any_results:
        return LiteratureCheckResult(
            status="inconclusive",
            citation=None,
            doi=None,
            confidence=0.0,
            query_trail=trail,
            notes="No search results returned at all -- the search could not be "
            "performed/trusted; rerun before believing a not-found.",
        )

    if best_conf >= INCONCLUSIVE_FLOOR:
        return LiteratureCheckResult(
            status="inconclusive",
            citation=best_hit.title if best_hit else None,
            doi=None,
            confidence=round(best_conf, 3),
            query_trail=trail,
            matched_url=best_hit.url if best_hit else None,
            notes="Cycler-adjacent literature surfaced but could not be confirmed "
            "as the same family; a human must adjudicate (not certified novel).",
        )

    return LiteratureCheckResult(
        status="not-found",
        citation=None,
        doi=None,
        confidence=round(best_conf, 3),
        query_trail=trail,
        notes="No published cycler matched the structural fingerprint. "
        "NECESSARY-NOT-SUFFICIENT for novelty: absence of a hit is not evidence "
        "of absence; the human + V0-V5 gauntlet still govern.",
    )


# ---------------------------------------------------------------------------
# The novelty GATE
# ---------------------------------------------------------------------------


def is_novelty_claimable(literature_check: dict[str, Any] | None) -> bool:
    """Is a candidate novelty-claimable given its ``literature_check`` block?

    The gate the daemon / review flow MUST consult before claiming novel:
    novelty is claimable iff the literature block is POPULATED **and** its status
    is not ``published``. An unpopulated (``None`` / not ``checked``) block or a
    ``published`` status => NOT claimable. ``inconclusive`` is also not
    claimable (we never searched cleanly enough to clear it).

    This is necessary-not-sufficient: a ``True`` here only means "no published
    rediscovery was found"; the human + V0-V5 gauntlet still decide novelty.
    """
    if not literature_check:
        return False
    if not literature_check.get("checked"):
        return False
    status = literature_check.get("status") or literature_check.get("result")
    # Accept either the structured status or the legacy result mapping.
    return status not in ("published", "match", "inconclusive", None)


def signature_from_review_entry(entry: Any) -> CandidateSignature:
    """Build a :class:`CandidateSignature` from a review-queue entry.

    Pulls primary from ``verdict_audit['primary']`` (the daemon stamps it),
    the encounter sequence + V_inf + period from the entry's own fields, and the
    per-leg ``n_rev`` from the audit block. Keeps the literature check wired to
    exactly the fields the daemon already records.
    """
    audit = getattr(entry, "verdict_audit", {}) or {}
    n_rev = audit.get("n_rev") or []
    return CandidateSignature(
        primary=audit.get("primary", "Sun"),
        sequence=tuple(entry.sequence),
        period_k=getattr(entry, "period_k", None),
        vinf_per_encounter_kms=tuple(entry.vinf_per_encounter_kms),
        n_rev=tuple(int(x) for x in n_rev),
    )


__all__ = [
    "INCONCLUSIVE_FLOOR",
    "KNOWN_CORPUS",
    "MATCH_THRESHOLD",
    "CandidateSignature",
    "CorpusAnchor",
    "LiteratureCheckResult",
    "SearchFn",
    "SearchResult",
    "build_queries",
    "check_literature",
    "is_novelty_claimable",
    "signature_from_review_entry",
]
