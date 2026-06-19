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

    period_band_tu: tuple[float, float] | None = None
    """Optional CR3BP nondim-period band (T_min, T_max) the candidate occupies.

    When set AND an overlapping anchor declares its own ``period_band_tu``, the
    candidate is treated as out-of-anchor-scope iff the two bands are disjoint.
    ``None`` (the historical default) preserves the old behaviour: no period
    filter, the anchor's structural-fingerprint footprint dominates. Used by
    #301 to escape the Antoniadou-Voyatzis 2018 anchor's low-integer scope when
    the candidate sits at a period-multiplied (k>1) sub-family.
    """

    topology_label: frozenset[str] = frozenset()
    """Optional topology classifier set (e.g. ``{"repeated-moon"}``).

    When BOTH the candidate AND an overlapping anchor declare non-empty
    ``topology_label`` sets, the matcher requires a non-empty intersection
    -- a candidate's topology must be one the anchor's published scope
    actually covers. Empty (the default) means "topology unrestricted":
    preserves the historical body-set-only matching for all anchors not
    yet annotated. Used by #349 to discriminate the Cassini-Huygens
    Titan-pump tour from a (k1, k2) repeated-moon cycler candidate
    despite sharing the {Titan, Rhea} body subset.

    Standard labels: ``"repeated-moon"`` (Aldrin / k1,k2 repeating-encounter
    cyclers), ``"pump-tour"`` (V-infinity-leveraging Titan-pump style),
    ``"mga-tour"`` (non-repeating multi-flyby tours, e.g. Galileo VEEGA),
    ``"tulip"`` (Sundman/petal Np-petal periodic orbits), ``"halo"``,
    ``"nrho"``, ``"resonant"``, ``"binary-coorbital"``.
    """

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
    period_band_tu: tuple[float, float] | None = None
    """Optional CR3BP nondim-period band the anchor's published scope covers.

    When set, a candidate whose ``period_band_tu_min`` lies outside this band is
    treated as out-of-scope for the anchor (the anchor is NOT used to flag it
    a rediscovery). ``None`` (default) means the anchor has no period-scope
    restriction. Used by #301 to filter the Antoniadou-Voyatzis 2018 anchor's
    low-integer-resonance scope away from the period-multiplied Neimark-Sacker
    sub-families at T_TU ~ 20-44.
    """

    topology_label: frozenset[str] = frozenset()
    """Optional topology classifier set this anchor's published scope covers.

    Sibling of :attr:`CandidateSignature.topology_label`; the matcher requires
    a non-empty intersection only when BOTH sides are non-empty. Empty (the
    default) means the anchor doesn't restrict by topology -- body-set-only
    matching, the historical behaviour. Added by #349 to discriminate the
    Cassini-Huygens Titan-pump tour anchor from a candidate (k1, k2)
    repeated-moon cycler despite shared {Titan, Rhea} body subset. See the
    sibling field on ``CandidateSignature`` for the standard labels.
    """


# Hand-curated from the catalogue's published rows + the task's named corpus.
# These are PUBLICATION facts (author/venue/doi), not values our code computed.
KNOWN_CORPUS: tuple[CorpusAnchor, ...] = (
    CorpusAnchor(
        name="Aldrin Earth-Mars cycler",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: Aldrin's canonical paradigm IS the (k1, k2) repeated-encounter
        # cycler -- the defining instance of the class.
        topology_label=frozenset({"repeated-moon"}),
        # #364: combined 1993 JSR (Byrnes-Longuski-Aldrin) + 2002 AIAA
        # (McConaghy-Longuski-Byrnes) author rosters.
        authors=("Aldrin", "Byrnes", "Longuski", "McConaghy"),
        keywords=("Aldrin cycler", "Earth-Mars cycler", "cycler trajectory"),
        # #364 errata fix: the canonical archival Aldrin cycler paper is the
        # 1993 JSR (Byrnes-Longuski-Aldrin, DOI 10.2514/3.25519); AIAA 2002-4420
        # is a follow-up. Per docs/notes/2026-06-17-digest-byrnes-longuski-
        # aldrin-1993.md the JSR paper publishes the circular-coplanar (a, e)
        # derivation and the 15-year DE405 numerical results (Tables 1-2).
        citation=(
            "Byrnes, D. V., Longuski, J. M. & Aldrin, B., 'Cycler Orbit Between "
            "Earth and Mars,' J. Spacecraft & Rockets 30(3):334-336 (1993), "
            "DOI 10.2514/3.25519 (canonical archival Aldrin cycler paper); "
            "McConaghy, T. T., Longuski, J. M. & Byrnes, D. V., 'Analysis of a "
            "Broad Class of Earth-Mars Cycler Trajectories,' AIAA-2002-4420 "
            "(Aldrin-class follow-up)"
        ),
        doi="10.2514/3.25519",
    ),
    CorpusAnchor(
        name="Russell-Ocampo / McConaghy Earth-Mars SnLm cyclers",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: SnLm is the systematic catalogue of repeated-encounter
        # Earth-Mars cyclers (S = synodic returns, L = sidereal returns).
        topology_label=frozenset({"repeated-moon"}),
        authors=("Russell", "Ocampo", "McConaghy", "Landau", "Longuski", "Byrnes"),
        keywords=(
            "ballistic Earth-Mars cycler",
            "S1L1 cycler",
            "two-synodic cycler",
            "systematic cycler catalog",
        ),
        # #364 errata fix: Russell-Ocampo is JGCD 27(3):321-335 (2004) DOI
        # 10.2514/1.1909, NOT JSR 41(1) DOI 10.2514/1.10078. The peer-reviewed
        # version stamp on the AAS-03-145 preprint p.1 establishes the JGCD
        # citation unambiguously. Per docs/notes/2026-06-17-digest-russell-
        # ocampo-2003.md.
        citation=(
            "Russell, R. P. & Ocampo, C. A., 'Systematic Method for "
            "Constructing Earth-Mars Cyclers Using Free-Return Trajectories,' "
            "J. Guidance, Control, and Dynamics 27(3):321-335 (2004), "
            "DOI 10.2514/1.1909 (preprint AAS-03-145); "
            "McConaghy et al., J. Spacecraft & Rockets 43(2) 2006; "
            "Byrnes, D. V., McConaghy, T. T. & Longuski, J. M., 'Analysis of "
            "Various Two Synodic Period Earth-Mars Cycler Trajectories,' "
            "AIAA/AAS Astrodynamics Specialist Conf., Monterey CA, Aug 2002 "
            "(#384; S1L1-B / Case-3 two-synodic precedence + real-eph V_inf "
            "envelope E 4.15-7.44 / M 2.97-7.83 km/s that brackets the spec-9 "
            "S1L1 anchor 5.65/3.05)"
        ),
        doi="10.2514/1.1909",
    ),
    CorpusAnchor(
        name="Liang et al. Callisto-Ganymede-Europa triple cyclers",
        primary="Jupiter",
        body_set=frozenset({"Callisto", "Ganymede", "Europa"}),
        # #350: 'triple cycler' = repeated CGE encounter sequence.
        topology_label=frozenset({"repeated-moon"}),
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
        # #350: 'triple cycler' = repeated IEG encounter sequence.
        topology_label=frozenset({"repeated-moon"}),
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
        # #350: V-infinity-leveraging is the same-body Tisserand-pump
        # methodology (AAS 07-277 abstract: 'transfers between the same
        # gravity-assist body'). Same family as the Cassini pump-tour.
        topology_label=frozenset({"pump-tour", "mga-tour"}),
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
        # #385: Campagnola-Buffington-Petropoulos 2014 EHM Europa orbiter +
        # lander Jovian tour (Acta Astro 100:68-81). Distinct paper / DOI from
        # the Strange/Campagnola/Russell anchor above (Campagnola & Russell
        # JGCD 2010). Tours 11-O3 / 12-L1 / 12-L4 with full per-encounter
        # V_inf tables (digest docs/notes/2026-06-17-digest-campagnola-2014.md).
        # topology mapped to the standard vocabulary: Tisserand-Poincare graph
        # design IS a pump-tour / mga-tour, not a separate topology label.
        name="Campagnola-Buffington-Petropoulos EHM Europa orbiter+lander Jovian tour (2014)",
        primary="Jupiter",
        body_set=frozenset({"Europa", "Ganymede", "Callisto"}),
        topology_label=frozenset({"pump-tour", "mga-tour"}),
        authors=("Campagnola", "Buffington", "Petropoulos"),
        keywords=(
            "Europa orbiter lander Jovian tour",
            "Tisserand-Poincare graph moon tour",
            "EHM Europa Habitability Mission tour",
        ),
        citation="Campagnola, Buffington & Petropoulos, 'Jovian tour design "
        "for orbiter and lander missions to Europa,' Acta Astronautica "
        "100:68-81 (2014); tours 11-O3 / 12-L1 / 12-L4.",
        doi="10.1016/j.actaastro.2014.02.005",
    ),
    CorpusAnchor(
        # #385: Niehoff 1970 'Touring the Galilean Satellites' (AIAA 70-1070)
        # -- the earliest sourced multi-flyby Galilean-moon tour paper in the
        # corpus, predating Strange-Russell-Buffington 2007 by 37 yr and the
        # Niehoff VISIT cyclers (already catalogued) by 15 yr. Foundational
        # resonance-locked Jovian tour paradigm (Mode 3, 14-day orbit, 73
        # encounters / 170 days). Digest
        # docs/notes/2026-06-17-digest-niehoff-1970.md. The digest's free-text
        # 'resonance-locked / laplace-syzygy' labels map to the standard
        # 'resonant' topology vocabulary; the tour structure is pump-tour /
        # mga-tour.
        name="Niehoff Galilean multi-flyby tour paradigm (1970, foundational)",
        primary="Jupiter",
        body_set=frozenset({"Io", "Europa", "Ganymede", "Callisto"}),
        topology_label=frozenset({"pump-tour", "mga-tour", "resonant"}),
        authors=("Niehoff",),
        keywords=(
            "Galilean satellite tour",
            "resonance-locked Jupiter orbit multi-flyby tour",
            "syzygy phase-locked satellite encounter sequence",
        ),
        citation="Niehoff, 'Touring the Galilean Satellites,' AAS/AIAA "
        "Astrodynamics Conference, Santa Barbara, Aug 19-21, 1970; AIAA "
        "Paper 70-1070. Foundational Galilean multi-flyby tour paradigm.",
        doi=None,
    ),
    CorpusAnchor(
        name="Jones et al. VEM triple cyclers (Venus-Earth-Mars)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M"}),
        # #350: 'VEM triple cycler' = repeated Venus-Earth-Mars encounter
        # sequence (the Jones-Hernandez-Jesick AAS 17-577 family).
        topology_label=frozenset({"repeated-moon"}),
        authors=("Jones", "Hernandez", "Jesick"),
        keywords=("VEM triple cycler", "Venus-Earth-Mars cycler"),
        citation="Jones, Hernandez & Jesick, AAS 17-577 (VEM triple cyclers)",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # Pluto-Charon dynamical literature (#272 expansion).
    #
    # The 12-row Pluto SILVER pass (#269 / 2026-06-15-pluto-silver-review.md)
    # identified three structural classes the offline corpus was missing. The
    # entries below pre-register the published anchors so the offline matcher
    # surfaces them deterministically without needing live web search; web
    # search remains authoritative for novelty.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Howard et al. Persephone Pluto-system orbiter (CR3BP periodic orbits)",
        primary="Pluto",
        body_set=frozenset({"Charon", "Nix", "Hydra", "Styx", "Kerberos"}),
        authors=("Howard", "Stern", "McKinnon"),
        keywords=(
            "Persephone Pluto orbiter",
            "Pluto-Charon CR3BP periodic orbit",
            "Pluto-Charon binary rotating frame orbit",
        ),
        citation="Howard, Stern et al., 'Persephone: A Pluto-system Orbiter "
        "and Kuiper Belt Explorer,' Planetary Science Journal 2(2):56 (2021); "
        "arXiv:2102.08282",
        doi="10.3847/PSJ/abf837",
    ),
    CorpusAnchor(
        name="Stern/SwRI Pluto Game-Changer (Charon gravity-assist tour)",
        primary="Pluto",
        body_set=frozenset({"Charon", "Nix", "Hydra", "Styx", "Kerberos"}),
        authors=("Stern", "Tapley", "Zangari"),
        keywords=(
            "Charon gravity assist tour",
            "Pluto orbiter Game-Changer",
            "Pluto system flyby tour",
        ),
        citation="Stern, Tapley, Zangari et al., 'Game-Changer Pluto Orbiter' "
        "concept, DPS 2018 (SwRI)",
        doi=None,
    ),
    CorpusAnchor(
        name="Showalter-Hamilton Styx-Nix-Hydra three-body resonance",
        primary="Pluto",
        body_set=frozenset({"Styx", "Nix", "Hydra"}),
        # #350: natural Laplace-like three-body resonance (Nature 2015) -- not
        # an engineered trajectory, observational dynamics paper.
        topology_label=frozenset({"resonant"}),
        authors=("Showalter", "Hamilton"),
        keywords=(
            "Styx Nix Hydra three-body resonance",
            "Pluto small moons Laplace-like resonance",
        ),
        citation="Showalter & Hamilton, 'Resonant interactions and chaotic "
        "rotation of Pluto's small moons,' Nature 522:45-49 (2015)",
        doi="10.1038/nature14469",
    ),
    CorpusAnchor(
        name="Brozovic et al. Pluto satellite orbit determination",
        primary="Pluto",
        body_set=frozenset({"Charon", "Styx", "Nix", "Kerberos", "Hydra"}),
        authors=("Brozovic", "Showalter", "Jacobson"),
        keywords=(
            "Pluto small satellite orbits",
            "Pluto satellite system orbit determination",
        ),
        citation="Brozovic et al., 'The orbits and masses of satellites of "
        "Pluto,' Icarus 246:317-329 (2015); Brozovic et al., AJ 163:241 (2022)",
        doi="10.1016/j.icarus.2014.03.015",
    ),
    CorpusAnchor(
        name="Pluto-Charon CR3BP tadpole/horseshoe at L3/L4/L5 (2025)",
        primary="Pluto",
        body_set=frozenset({"Charon"}),
        # #350: Tadpole + horseshoe = coorbital binary topology (the L3/L4/L5
        # equilateral-Lagrange family for the Pluto-Charon binary).
        topology_label=frozenset({"binary-coorbital"}),
        authors=("Pluto-Charon", "Charon"),  # arXiv preprint; authors not pinned in our notes
        keywords=(
            "Pluto-Charon tadpole horseshoe orbit",
            "Pluto-Charon L4 L5 L3 periodic orbit",
            "Pluto-Charon CR3BP coorbital",
        ),
        citation="'Tadpole and horseshoe orbits in the Pluto-Charon CR3BP at "
        "L3/L4/L5,' arXiv:2510.13479 (2025)",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # Recent capability papers (Track A, surfaced by the #265 capability sweep
    # 2026-06-13-discovery-capability-paper-sweep.md / forward-citation-sweep-2).
    # These are anchored here so the offline corpus + structural matcher can
    # cite the published genome / family-network / tulip / lobe-dynamics
    # papers when a candidate's fingerprint overlaps them.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Braik-Ross orbital networks in the three-body problem (2026)",
        primary="Earth",  # Earth-Moon CR3BP family network
        body_set=frozenset({"Moon"}),
        authors=("Braik", "Ross"),
        keywords=(
            "orbital network three-body problem",
            "reachable-set family accessibility",
            "Earth-Moon CR3BP family network",
        ),
        citation="Braik & Ross, 'Orbital Networks in the Three-Body Problem,' "
        "arXiv:2605.31543 (2026)",
        doi=None,
    ),
    CorpusAnchor(
        name="Roberts-Tsoukkas & Ross stable prograde Earth-Moon multi-orbiter cyclers",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        # #350: Earth-Moon multi-orbiter cyclers = repeated lunar-encounter
        # family (the 5 stable EM cyclers reproduced as V1 in catalogue).
        topology_label=frozenset({"repeated-moon"}),
        authors=("Roberts-Tsoukkas", "Ross"),
        keywords=(
            "stable prograde Earth-Moon cycler",
            "multi-orbiter cycler three-body dynamics",
            "binary-star mass parameter cycler family",
        ),
        citation="Roberts-Tsoukkas & Ross, 'Stable Prograde Earth-Moon "
        "Multi-Orbiter Cyclers via Three-Body Dynamics,' journal extension of "
        "AAS 25-621 (2026, VSGC manuscript)",
        doi=None,
    ),
    CorpusAnchor(
        name="Kumar-Rawat-Rosengren-Ross cislunar resonant transport (2026)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        authors=("Kumar", "Rawat", "Rosengren", "Ross"),
        keywords=(
            "cislunar resonant transport",
            "heteroclinic pathway resonant manifold",
            "Earth-Moon resonant family network",
        ),
        citation="Kumar, Rawat, Rosengren & Ross, 'Cislunar Resonant Transport "
        "and Heteroclinic Pathways,' Advances in Space Research (2026); "
        "arXiv:2509.12675",
        doi=None,
    ),
    CorpusAnchor(
        name="Koblick novel tulip-shaped three-body orbits (cislunar SDA)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        # #350: tulip-shaped Np-petal three-body orbit family (the #266
        # genome's source paper).
        topology_label=frozenset({"tulip"}),
        authors=("Koblick", "Kelly"),
        keywords=(
            "tulip-shaped three-body orbit",
            "tulip cislunar SDA orbit",
            "petal-count periodic orbit Earth-Moon",
        ),
        citation="Koblick, 'Novel Tulip-Shaped Three-body Orbits for Cislunar "
        "SDA Missions,' AMOSTECH (2023); Koblick & Kelly, J. Astronaut. Sci. "
        "(2025) DOI 10.1007/s40295-025-00510-w",
        doi="10.1007/s40295-025-00510-w",
    ),
    CorpusAnchor(
        name="Zhang-Jiang-Yuan tulip time-regularized bifurcation framework (2026)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        # #350: tulip orbit period-multiplying bifurcation methodology.
        topology_label=frozenset({"tulip"}),
        authors=("Zhang", "Jiang", "Yuan"),
        keywords=(
            "tulip orbit time-regularized bifurcation",
            "period-multiplying bifurcation cislunar tulip",
        ),
        citation="Zhang, Jiang & Yuan, 'Time-regularized bifurcation framework "
        "for tulip-shaped orbits,' Nonlinear Dynamics (2026) DOI "
        "10.1007/s11071-026-12465-0",
        doi="10.1007/s11071-026-12465-0",
    ),
    CorpusAnchor(
        name="Cislunar tulip robust construction (Chinese J. Aeronautics 2026)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        # #350: tulip orbit family construction methodology.
        topology_label=frozenset({"tulip"}),
        authors=("Chinese Journal of Aeronautics",),
        keywords=(
            "cislunar tulip robust construction",
            "tulip orbit cislunar mission",
        ),
        citation="'Robust construction of cislunar tulip-shaped orbits,' "
        "Chinese Journal of Aeronautics (2026), Elsevier "
        "S1000936126001755",
        doi=None,
    ),
    CorpusAnchor(
        name="Hiraiwa et al. lobe-dynamics low-energy cislunar transfers (2026)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        authors=("Hiraiwa", "Bando", "Sato", "Hokamoto"),
        keywords=(
            "lobe dynamics low-energy cislunar transfer",
            "weighted directed graph lobe sequence",
            "resonant orbit lobe-dynamics transfer",
        ),
        citation="Hiraiwa, Bando, Sato & Hokamoto, 'Design of low-energy "
        "transfers in cislunar space using sequences of lobe dynamics,' "
        "Acta Astronautica 248 (2026); arXiv:2602.17444",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # Ancillary anchors -- adjacent literature surfaced by the discovery
    # campaign that the matcher may bump into structurally.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Davis-Phillips-McCarthy Saturnian Ocean Worlds orbiters (2018)",
        primary="Saturn",
        # Body-set tightened by #346 PDF deep-read: paper's documented orbit
        # families are Saturn-Titan TULIP orbits (Fig. 8-11, ~6:1 resonance,
        # periapses over Titan's poles) and Saturn-Enceladus L1/L2 NRHO halos
        # + a heteroclinic transfer Enceladus-halo -> Titan polar orbit. Rhea
        # and Dione are NOT orbit-family targets in the paper, so a
        # repeated-moon (k1,k2) Titan-Rhea or Titan-Dione sequence is NOT
        # captured by this anchor and remains lit-fresh.
        body_set=frozenset({"Titan", "Enceladus"}),
        # #349: topology_label set per the same #346 deep-read -- the paper's
        # actual orbit families are tulip-shaped resonant orbits at Titan and
        # NRHO halo families at Enceladus. NOT a repeated-moon cycler, NOT a
        # pump tour. Belt-and-suspenders alongside the tightened body_set.
        topology_label=frozenset({"tulip", "halo", "nrho"}),
        authors=("Davis", "Phillips", "McCarthy"),
        keywords=(
            "Saturnian Ocean Worlds Poincare map",
            "Saturn-Titan tulip-shaped orbit polar coverage",
            "Saturn-Enceladus NRHO halo heteroclinic",
        ),
        citation="Davis, Phillips & McCarthy, 'Trajectory design for "
        "Saturnian Ocean Worlds orbiters using multidimensional Poincare "
        "maps,' Acta Astronautica 143:16-28 (2018)",
        doi="10.1016/j.actaastro.2017.11.004",
    ),
    CorpusAnchor(
        name="Sanaga-Park-Howell fidelity-transition framework (2026)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        authors=("Sanaga", "Park", "Howell"),
        keywords=(
            "fidelity transition framework CR3BP ephemeris",
            "CR3BP to high-fidelity continuation",
        ),
        citation="Sanaga, Park & Howell, J. Astronaut. Sci. (2026) DOI 10.1007/s40295-026-00571-5",
        doi="10.1007/s40295-026-00571-5",
    ),
    CorpusAnchor(
        # #357: Singh-Anderson-Taheri-Junkins 2021 (Acta Astro 183:255-272)
        # end-to-end GTO -> lunar low-thrust via Earth-Moon L1 halo manifolds.
        # 3 sourced 15-digit L1 halo ICs at C = 3.128 / 3.143 / 3.158 (Table 2
        # p.257). Digest docs/notes/2026-06-17-digest-singh-2021-L1-halo.md.
        name="Singh-Anderson-Taheri-Junkins Earth-Moon L1 halo low-thrust manifolds (2021)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        topology_label=frozenset({"halo"}),
        authors=("Singh", "Anderson", "Taheri", "Junkins"),
        keywords=(
            "Earth-Moon L1 halo manifold low-thrust transfer",
            "GTO to lunar polar orbit low-thrust",
            "bi-circular problem halo manifold",
        ),
        citation="Singh, Anderson, Taheri & Junkins, 'Exploiting manifolds of "
        "L1 halo orbits for end-to-end Earth-Moon low-thrust trajectory "
        "design,' Acta Astronautica 183:255-272 (2021); Table 2 p.257 "
        "(3 L1 halo ICs at C = 3.128 / 3.143 / 3.158).",
        doi="10.1016/j.actaastro.2021.03.017",
    ),
    CorpusAnchor(
        # #357: Singh-Anderson-Taheri-Junkins 2021 (JOTA 191(2-3)) low-thrust
        # transfers to 3 Southern L2 NRHOs via invariant manifolds. 3 sourced
        # NRHO ICs (9:2 Gateway / 24:5 / 4:1) with stability indices (Table 1
        # p.5). Digest docs/notes/2026-06-17-digest-singh-2021-NRHO.md. The
        # 9:2 Gateway NRHO stability indices (-1.3753 / 0.6626) are also #347
        # Floquet-multiplier ground truth.
        name="Singh-Anderson-Taheri-Junkins Earth-Moon Southern L2 NRHO manifolds (2021)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        topology_label=frozenset({"nrho", "halo"}),
        authors=("Singh", "Anderson", "Taheri", "Junkins"),
        keywords=(
            "Southern L2 NRHO low-thrust transfer invariant manifold",
            "9:2 Gateway NRHO",
            "near-rectilinear halo orbit stability index",
        ),
        citation="Singh, Anderson, Taheri & Junkins, 'Low-Thrust Transfers to "
        "Southern L2 Near-Rectilinear Halo Orbits Facilitated by Invariant "
        "Manifolds,' J. Optimization Theory & Applications 191(2-3) (2021); "
        "Table 1 p.5 (9:2 Gateway / 24:5 / 4:1 NRHO ICs + stability indices); "
        "preliminary AAS 20-565.",
        doi="10.1007/s10957-021-01898-9",
    ),
    # -----------------------------------------------------------------------
    # MGA / pump-tour / cycler-precursor literature (#294 scope expansion).
    #
    # Catalogue scope expanded 2026-06-15 from cyclers-only to a four-class
    # taxonomy: cycler / quasi_cycler / precursor_mga / mga_tour. The cluster
    # below anchors the published epoch-locked mission-design literature so
    # candidates in the new classes can be literature-checked before any
    # novelty claim. Galileo VEEGA, Cassini VVEJGA, Petropoulos pump tours,
    # Heaton-Longuski resonance hopping, etc.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Petropoulos-Longuski Jovian pump-tour combinatorics (2000)",
        primary="Jupiter",
        body_set=frozenset({"Io", "Europa", "Ganymede", "Callisto"}),
        # #350: Tisserand-graph pump-tour design (foundational paper).
        topology_label=frozenset({"pump-tour", "mga-tour"}),
        authors=("Petropoulos", "Longuski"),
        keywords=(
            "Jovian moon tour pump-up",
            "Tisserand graph pump tour",
            "satellite-to-satellite gravity-assist tour",
            "ballistic moon tour combinatorics",
        ),
        citation="Petropoulos & Longuski, 'A shape-based algorithm for "
        "automated design of low-thrust, gravity-assist trajectories,' "
        "J. Spacecraft & Rockets / AAS context (2000-2004); pump-tour "
        "enumeration foundational paper",
        doi=None,
    ),
    CorpusAnchor(
        name="Strange-Russell Tisserand pump-tour graph (2007 AAS 07-277)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M", "Jupiter"}),
        # #350: Heliocentric Tisserand-graph MGA-tour graph search (same
        # AAS-07-277 paper that grounds the Cassini-Huygens anchor's
        # same-body-pump methodology, here applied to heliocentric pump
        # tours).
        topology_label=frozenset({"pump-tour", "mga-tour"}),
        authors=("Strange", "Russell", "Buffington"),
        keywords=(
            "Tisserand graph MGA tour",
            "V-infinity globe mapping",
            "Tisserand pump tour design",
            "patched-conic flyby tour graph",
        ),
        citation="Strange, Russell & Buffington, 'Mapping the V-infinity "
        "Globe,' AAS 07-277 (2007); complements the line-246 entry by "
        "covering the heliocentric MGA-tour graph search",
        doi=None,
    ),
    CorpusAnchor(
        name="Heaton-Strange-Longuski resonance-hopping pump tours (2002)",
        primary="Jupiter",
        body_set=frozenset({"Io", "Europa", "Ganymede", "Callisto"}),
        # #350: resonance-hopping is same-body pump methodology applied to
        # Jovian moon tour design (Galileo / JIMO context).
        topology_label=frozenset({"pump-tour", "mga-tour"}),
        authors=("Heaton", "Strange", "Longuski"),
        keywords=(
            "resonance hopping moon tour",
            "Jovian resonance pump",
            "pump-down pump-up cycler tour",
        ),
        citation="Heaton, Strange & Longuski, 'Automated Design of the "
        "Europa Orbiter Tour' (2002 AAS / 2002-4727 AIAA Astrodynamics "
        "Specialist context); resonance-hop pump tour methodology",
        doi=None,
    ),
    CorpusAnchor(
        name="Vasile-Conway MGA-DSM optimisation (2006)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M", "Jupiter"}),
        # #350: MGA-DSM = multi-gravity-assist with deep-space maneuvers,
        # the canonical mga_tour optimisation paradigm.
        topology_label=frozenset({"mga-tour"}),
        authors=("Vasile", "Conway", "De Pascale"),
        keywords=(
            "MGA-DSM multi-gravity-assist deep-space maneuver",
            "global optimisation interplanetary trajectory",
            "GTOC MGA-DSM",
        ),
        citation="Vasile, Conway et al., MGA-DSM global-optimisation "
        "framework (Acta Astronautica / Conway book ch.) 2006-2009",
        doi=None,
    ),
    CorpusAnchor(
        name="Hughes-Edelman-Longuski VEM cycler extensions (2014)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M"}),
        # #350: extends Jones-Hernandez-Jesick AAS 17-577 VEM cycler family
        # (repeated Venus-Earth-Mars encounter sequence).
        topology_label=frozenset({"repeated-moon"}),
        authors=("Hughes", "Edelman", "Longuski"),
        keywords=(
            "Venus-Earth-Mars cycler extension",
            "VEM tour extension",
            "outbound-inbound Venus-Earth-Mars",
        ),
        citation="Hughes, Edelman & Longuski, AAS 14-822 / 'Venus-Earth-Mars "
        "Cyclers' extension paper (2014)",
        doi=None,
    ),
    CorpusAnchor(
        name="Genova-Aldrin purple Earth-Mars cycler precursors (2015)",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: Aldrin/purple cyclers + their precursor insertion trajectories
        # are repeated-encounter Earth-Mars families.
        topology_label=frozenset({"repeated-moon"}),
        authors=("Genova", "Aldrin"),
        keywords=(
            "purple Earth-Mars cycler",
            "Aldrin cycler precursor insertion",
            "Mars cycler insertion trajectory",
        ),
        citation="Genova & Aldrin, 'Mars Human Exploration: Aldrin and Purple "
        "Cyclers' (2015 AIAA / AAS context); cycler precursor insertion "
        "trajectories",
        doi=None,
    ),
    CorpusAnchor(
        name="McConaghy Earth-Mars cycler dissertation (2004)",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: McConaghy's PhD is the dissertation source for the SnLm
        # repeated-encounter Earth-Mars cycler catalogue.
        topology_label=frozenset({"repeated-moon"}),
        authors=("McConaghy",),
        keywords=(
            "Earth-Mars ballistic cycler dissertation",
            "Aldrin S1L1 cycler member",
            "low-thrust cycler insertion",
        ),
        citation="McConaghy, T.T., 'Design and Optimization of Interplanetary "
        "Trajectories' PhD dissertation, Purdue University (2004); the "
        "single-author Earth-Mars cycler dissertation underpinning the "
        "SnLm catalogue rows in data/catalogue.yaml",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # #364 — McConaghy 2004 JSR + McConaghy 2005 JSR + Rogers 2015 Acta Astro
    # added as Mars-cycler corpus anchors (errata fix wave from Agent A
    # digest pass). The 2004 JSR introduces the nPr family-tag nomenclature
    # and publishes the S1L1 DE405 itinerary (Table 6, V1-grade ground truth).
    # The 2005 JSR introduces the formal per-leg g/f/h descriptor
    # nomenclature (Table 2 Rosetta stone). The 2015 Acta Astro is the
    # establishment-cycler / precursor_mga class anchor — Tables 1-9 give
    # ~21 V1-grade insertion trajectories for eight cycler families.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name=("McConaghy-Longuski-Byrnes Earth-Mars cycler trajectory class analysis (2004)"),
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: nPr family + S1L1 DE405 itinerary are repeated-encounter
        # Earth-Mars cyclers.
        topology_label=frozenset({"repeated-moon"}),
        authors=("McConaghy", "Longuski", "Byrnes"),
        keywords=(
            "nPr nomenclature Earth-Mars cycler",
            "S1L1 ballistic cycler DE405 itinerary",
            "Earth-Mars cycler trajectory class",
            "n=7 ballistic cycler family",
        ),
        citation=(
            "McConaghy, T. T., Longuski, J. M. & Byrnes, D. V., 'Analysis "
            "of a Class of Earth-Mars Cycler Trajectories,' J. Spacecraft "
            "& Rockets 41(4):622-628 (2004), DOI 10.2514/1.11939. "
            "Introduces the nPr family-tag nomenclature (Aldrin = 1L1, "
            "Case 1 = 2L3, VISIT 1 = 7(R_p)12, VISIT 2 = 7(R_p)10). Table 4 "
            "lists 21 most-promising cyclers (1 <= n <= 6); Table 5 lists "
            "14 all-ballistic n=7 cyclers; Table 6 gives the 22-encounter "
            "outbound ballistic S1L1 cycler DE405 itinerary (launch "
            "9 June 2008) as V1-grade ground truth"
        ),
        doi="10.2514/1.11939",
    ),
    CorpusAnchor(
        name=("McConaghy-Russell-Longuski Earth-Mars cycler standard nomenclature (2005)"),
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: per-leg formal labels for repeated-encounter Earth-Mars cyclers.
        topology_label=frozenset({"repeated-moon"}),
        authors=("McConaghy", "Russell", "Longuski"),
        keywords=(
            "standard nomenclature Earth-Mars cycler",
            "g f h leg descriptor cycler label",
            "Lambert solution subtype U L cycler",
            "Ballistic S1L1 formal label",
        ),
        citation=(
            "McConaghy, T. T., Russell, R. P. & Longuski, J. M., 'Toward a "
            "Standard Nomenclature for Earth-Mars Cycler Trajectories,' "
            "J. Spacecraft & Rockets 42(4):694-698 (2005), DOI "
            "10.2514/1.8123. Introduces the formal per-leg [(body-seq)] n "
            "d_1...d_K label with d in {g(t_f, theta, eps), f(M:N, phi, "
            "lambda), h(t_f, N, eps, i')}; EBNF grammar in Table 3; Table 2 "
            "tabulates Aldrin = 1g(2-1/7, 1-1/7 rev, L), Ballistic S1L1 = "
            "2g(2.8277, 657.97 deg, U) g(1.4508, 522.29 deg, L) (confirms "
            "S1L1's two-arc structure), VISIT-1/2, Byrnes' Case 3, Russell "
            "Cycler-2.5.1.+0 and Cycler-4.3.1.-5 as Rosetta-stone formal "
            "labels"
        ),
        doi="10.2514/1.8123",
    ),
    CorpusAnchor(
        name=("Rogers-Hughes-Longuski-Aldrin Earth-Mars cycler establishment trajectories (2015)"),
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #364: establishment trajectories are precursor_mga (mga-tour) that
        # insert into the repeated-encounter cycler families (repeated-moon).
        # Both topology classes apply -- the cycler matcher should bump into
        # this anchor for either class. Per docs/notes/2026-06-17-digest-
        # rogers-2015.md section 4.
        topology_label=frozenset({"mga-tour", "repeated-moon"}),
        # Capped at 3 authors per memory rule; Aldrin (4th) named only in
        # the citation text.
        authors=("Rogers", "Hughes", "Longuski"),
        keywords=(
            "establishment cycler trajectory Earth-Mars",
            "K:L(M) V-infinity leveraging cycler insertion",
            "precursor cycler MGA low-thrust establishment",
            "cycler LEO insertion Sims-Longuski leveraging",
        ),
        citation=(
            "Rogers, B. A., Hughes, K. M., Longuski, J. M. & Aldrin, B., "
            "'Establishing cycler trajectories between Earth and Mars,' "
            "Acta Astronautica 112:114-125 (2015), DOI "
            "10.1016/j.actaastro.2015.03.002 (conf precursors AIAA "
            "Minneapolis 2012 + AAS Hilton Head 2013). Eight cycler families "
            "(Aldrin, VISIT-1, VISIT-2, Case 1, Case 2, Case 3, S1L1, U0L1) "
            "with per-cycler establishment Delta-v via Sims-Longuski K:L(M)+- "
            "V-infinity leveraging (Tables 3-5) and JPL MALTO low-thrust "
            "spiral (Tables 6-9); Table 1 confirms the Aldrin a=1.60 AU / "
            "e=0.393 / aphelion=2.23 AU / perihelion=0.97 AU values; ~21 "
            "V1-grade precursor_mga establishment trajectories"
        ),
        doi="10.1016/j.actaastro.2015.03.002",
    ),
    CorpusAnchor(
        name="Ceriotti GTOC-style MGA chain optimisation (2010)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M", "Jupiter"}),
        # #350: GTOC-style global optimisation of MGA chains (mga_tour class).
        topology_label=frozenset({"mga-tour"}),
        authors=("Ceriotti",),
        keywords=(
            "GTOC MGA tour optimisation",
            "global trajectory optimisation",
            "multi-gravity-assist tour design",
        ),
        citation="Ceriotti, M., 'Global Optimisation of Multiple Gravity "
        "Assist Trajectories,' University of Glasgow PhD (2010)",
        doi=None,
    ),
    CorpusAnchor(
        name="Vasile-Campagnola MGA-DSM tour optimisation (2009)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M", "Jupiter"}),
        # #350: MGA-DSM tour optimisation via Tisserand-Poincare graphs.
        topology_label=frozenset({"mga-tour"}),
        authors=("Vasile", "Campagnola"),
        keywords=(
            "MGA-DSM tour optimisation",
            "Jovian moon tour Tisserand-Poincaré",
            "multi-objective MGA-DSM",
        ),
        citation="Vasile & Campagnola, 'Design of Low-Energy Multi-Gravity "
        "Assist Trajectories Using Tisserand-Poincaré Graphs' (2009)",
        doi=None,
    ),
    CorpusAnchor(
        name="Diehl-Kaplan-Penzo / D'Amario-Byrnes Galileo design (1983, pre-Challenger)",
        primary="Sun",
        body_set=frozenset({"E", "Jupiter"}),
        # #356/#384 (2026-06-19): the Belbruno misattribution is now CORRECTED.
        # The acquired+digested 1983 design papers establish the true lineage
        # (docs/notes/2026-06-19-digest-damario-byrnes-1983-galileo-interplanetary.md
        # + 2026-06-17-digest-diehl-1983.md): the pre-Challenger Galileo concept
        # was a DIRECT Earth->Jupiter trajectory (NOT VEEGA), designed in two
        # companion AIAA-83 papers -- Diehl, Kaplan & Penzo (AIAA-83-0101,
        # Jovian satellite tour) and D'Amario & Byrnes (AIAA-83-0099,
        # interplanetary leg). "Belbruno" (weak-stability-boundary / Hiten lunar
        # capture, not Galileo design) and "Roberts" were NOT authors of these
        # papers -- both removed. The flown post-Challenger VEEGA is the separate
        # D'Amario-Bright-Wolf 1992 anchor below (body_set {V,E,Jupiter}).
        topology_label=frozenset({"mga-tour"}),
        authors=("Diehl", "Kaplan", "Penzo", "D'Amario", "Byrnes"),
        keywords=(
            "Galileo direct Earth-Jupiter trajectory",
            "Galileo satellite tour design",
            "pre-Challenger Galileo 1986-launch concept",
        ),
        citation="Diehl, Kaplan & Penzo, 'Satellite Tour Design for the Galileo "
        "Mission' (AIAA-83-0101) + D'Amario & Byrnes, 'Interplanetary Trajectory "
        "Design for the Galileo Mission' (AIAA-83-0099), AIAA 21st Aerospace "
        "Sciences Meeting, Reno, Jan 1983 -- the pre-Challenger DIRECT Earth-"
        "Jupiter 1986-launch concept (#384). The flown VEEGA is the "
        "D'Amario-Bright-Wolf 1992 anchor below.",
        doi=None,
    ),
    CorpusAnchor(
        name="D'Amario-Bright-Wolf Galileo VEEGA flown trajectory (1992)",
        primary="Sun",
        body_set=frozenset({"V", "E", "Jupiter"}),
        # #356 (2026-06-17): the post-launch / flown-trajectory canonical
        # Galileo VEEGA reference. Distinct paper from the Diehl 1986 anchor
        # above. Catalogue row damario-1992-galileo-veega (mga_tour, V0).
        topology_label=frozenset({"mga-tour"}),
        authors=("D'Amario", "Bright", "Wolf"),
        keywords=(
            "VEEGA Venus-Earth-Earth gravity assist",
            "Galileo Jupiter mission flown trajectory",
            "1989 VEEGA opportunity",
            "Earth-Earth 2-year resonance",
            "Gaspra Ida asteroid flyby",
        ),
        citation="D'Amario, Bright & Wolf, 'Galileo trajectory design,' "
        "Space Science Reviews 60(1-4), 23-78 (May 1992); DOI 10.1007/bf00216849; "
        "the canonical post-launch reference for the flown 1989 VEEGA "
        "(launch 18 Oct 1989, V 10 Feb 1990, E1 8 Dec 1990, E2 8 Dec 1992, "
        "JOI 8 Dec 1995). Catalogue row damario-1992-galileo-veega (mga_tour, V0).",
        doi="10.1007/bf00216849",
    ),
    CorpusAnchor(
        name="Tito-MacCallum 2018 Mars free-return mission design (2013)",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        # #350: Single-window Earth-Mars-Earth ballistic free-return
        # (catalogue row tito-2018-mars-free-return, mga_tour class).
        topology_label=frozenset({"mga-tour"}),
        authors=("Tito", "Anderson", "Carrico", "Hopkins", "Loucks", "Voels"),
        keywords=(
            "2018 Mars free-return",
            "Inspiration Mars manned flyby",
            "Earth-Mars-Earth ballistic free-return tour",
        ),
        citation="Tito, Anderson, Carrico, Hopkins, Loucks & Voels, "
        "'Feasibility Analysis for a Manned Mars Free-Return Mission in "
        "2018,' IEEE Aerospace Conference (2013); Tables III/IV. Catalogue "
        "row tito-2018-mars-free-return (mga_tour, V0)",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # #287 follow-up — spatial CR3BP corpus (rediscovered by the 3D-Aldrin
    # scoping spike at z0 = -0.241 nondim, ~93,000 km out-of-plane).
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Antoniadou-Voyatzis spatial resonant periodic orbits in CR3BP (2018)",
        primary="Earth",  # spike work; catalogue applies more generally
        body_set=frozenset({"Moon"}),
        # #350: 'Spatial Resonant Periodic Orbits' = resonant family.
        topology_label=frozenset({"resonant"}),
        authors=("Antoniadou", "Voyatzis"),
        keywords=(
            "spatial resonant periodic orbit",
            "3D CR3BP family",
            "out-of-plane Lyapunov-vertical family",
            "spatial three-body resonant periodic orbit",
        ),
        citation="Antoniadou & Voyatzis, 'Spatial Resonant Periodic Orbits "
        "in the Restricted Three-Body Problem,' (2018); arXiv:1811.09442. "
        "Anchor for #287's 3D Braik-Ross (1,1) family extension (likely "
        "rediscovery target).",
        doi=None,
        # The paper's published catalogue covers low-integer p:q resonant
        # orbits (typically 1:1, 2:1, 3:2) in the spatial CR3BP. Their families
        # sit at T_TU under ~15. The #299 Neimark-Sacker sub-families at T_TU
        # 20-44 are period-multiplied (k=3-6) derivatives that are outside the
        # paper's scope -- treat them as not-anchored on AV-2018 alone.
        period_band_tu=(0.0, 15.0),
    ),
    # -----------------------------------------------------------------------
    # #328 — Uranian dynamics + mission-design literature.
    #
    # Added after the #327 verified SILVER (Umbriel-Oberon (1,1) repeated-
    # encounter cycler at V_inf ~ 0.9 km/s) deep-dive lit pass
    # (docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md).
    #
    # The candidate cleared as CLEAN LITERATURE-FRESH at its specific topology
    # (no DIRECT MATCH found), but the surrounding Uranian published record
    # is broad enough that the offline corpus needed deterministic anchors so
    # any FUTURE Uranian candidate's structural fingerprint surfaces the
    # appropriate paper. None of these six anchors covers the (1,1) two-moon
    # Umbriel-Oberon cycler topology; they cover (a) one-shot satellite tours,
    # (b) single-moon MMR resonant orbits, (c) moon-pair halo-to-halo
    # one-shot manifold transfers, and (d) Decadal mission concepts.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Heaton-Longuski Galileo-style Uranian satellite tour (2003)",
        primary="Uranus",
        body_set=frozenset({"Miranda", "Ariel", "Umbriel", "Titania", "Oberon"}),
        # #350: Catalogue row heaton-longuski-2003-uranian-tour-u00-01
        # (mga_tour class). One-shot 40+-flyby tour, NOT periodic cycler.
        topology_label=frozenset({"mga-tour"}),
        authors=("Heaton", "Longuski"),
        keywords=(
            "Galileo-style Uranian satellite tour",
            "Uranian moon gravity-assist tour",
            "STOUR patched-conic Uranus tour",
            "Uranus orbiter moon flyby tour",
        ),
        citation="Heaton & Longuski, 'The Feasibility of a Galileo-Style Tour "
        "of the Uranian Satellites,' J. Spacecraft & Rockets 40(4):591-596 "
        "(2003); AIAA 2001-3859 / NTRS 20020021945. Foundational one-shot "
        "Uranus moon-tour anchor: 811-day three-phase tour with 40+ flybys "
        "across Miranda/Ariel/Umbriel/Titania/Oberon. NOT a periodic cycler.",
        doi="10.2514/2.3981",
    ),
    CorpusAnchor(
        name="Sims et al. polar Uranus orbiter & satellite tour (2014)",
        primary="Uranus",
        body_set=frozenset({"Miranda", "Ariel", "Umbriel", "Titania", "Oberon"}),
        # #350: Citation explicit: 'one-shot insertion tour, NOT cycler.'
        topology_label=frozenset({"mga-tour"}),
        authors=("Sims", "Finlayson", "Rinderle", "Vavrina", "Kawalkowski"),
        keywords=(
            "polar Uranus orbiter satellite tour",
            "conceptual mission design Uranus orbiter",
            "two-flyby-per-moon Uranus tour",
            "Uranus inclination reduction sequence",
        ),
        citation="Sims, Finlayson, Rinderle, Vavrina & Kawalkowski et al., "
        "'Conceptual mission design of a polar Uranus orbiter and satellite "
        "tour' (2014). Baseline 424-day, 619 m/s tour with two targeted "
        "flybys of each major moon. One-shot insertion tour, NOT cycler.",
        doi=None,
    ),
    CorpusAnchor(
        name="Kumar Uranus-Oberon PCRTBP MMR study (2025)",
        primary="Uranus",
        body_set=frozenset({"Oberon", "Titania"}),
        # #350: MMR = mean-motion resonance. Citation explicit: 'Single-moon
        # MMR topology, NOT moon-pair cycler.'
        topology_label=frozenset({"resonant"}),
        authors=("Kumar",),
        keywords=(
            "Uranus-Oberon PCRTBP mean motion resonance",
            "Uranus-Oberon unstable resonant periodic orbit",
            "Uranus-Titania-Oberon CCR4BP secondary resonance",
            "Uranian system heteroclinic resonance transition",
        ),
        citation="Kumar, 'Multi-shooting parameterization methods for "
        "invariant manifolds and heteroclinics of 2-DOF Hamiltonian Poincare "
        "maps, with applications to celestial resonant dynamics,' "
        "arXiv:2509.03655 (2025). Section 6.2 studies Uranus-Oberon PCRTBP "
        "3:4/4:5/5:6 exterior and 4:3/5:4/6:5 interior MMR unstable periodic "
        "orbits plus heteroclinic connections; extends to Uranus-Titania-"
        "Oberon CCR4BP secondary resonances. Single-moon MMR topology, NOT "
        "moon-pair cycler.",
        doi=None,
    ),
    CorpusAnchor(
        name="Canales-Howell-Fantino moon-to-moon analytical transfer (Titania-Oberon, 2021)",
        primary="Uranus",
        body_set=frozenset({"Titania", "Oberon"}),
        # #350: Citation explicit: 'L2 halo at Uranus-Titania -> L1 halo at
        # Uranus-Oberon... One-shot transfer between halo orbits, NOT
        # repeated-encounter cycler.'
        topology_label=frozenset({"halo"}),
        authors=("Canales", "Howell", "Fantino"),
        keywords=(
            "Titania-Oberon halo-to-halo transfer",
            "Uranian moon-to-moon analytical transfer MMAT",
            "Uranus-Titania L2 to Uranus-Oberon L1 manifold transfer",
            "2BP-CR3BP patched moon-to-moon transfer Uranus",
        ),
        citation="Canales, Howell & Fantino, 'Transfer design between "
        "neighborhoods of planetary moons in the circular restricted "
        "three-body problem: the moon-to-moon analytical transfer method,' "
        "Celest. Mech. Dyn. Astron. 133:36 (2021); arXiv:2110.03683. "
        "Uranian case study: L2 halo at Uranus-Titania -> L1 halo at "
        "Uranus-Oberon via unstable/stable manifolds. One-shot transfer "
        "between halo orbits, NOT repeated-encounter cycler. Companion "
        "FTLE-map extensions in Canales-Howell-Fantino JGCD 2023 "
        "(arXiv:2308.10029).",
        doi=None,
    ),
    CorpusAnchor(
        name="Jarmak QUEST Uranus orbiter New Frontiers concept (2020)",
        primary="Uranus",
        body_set=frozenset({"Miranda", "Ariel", "Umbriel", "Titania", "Oberon"}),
        # #350: Citation explicit: 'Polar orbiter, no satellite tour; not
        # cycler-relevant.' Mission concept = mga-tour class even when polar
        # orbiter has no extended satellite tour.
        topology_label=frozenset({"mga-tour"}),
        authors=("Jarmak", "Brinckerhoff"),
        keywords=(
            "QUEST Uranus orbiter New Frontiers",
            "Uranus polar orbit mission concept",
            "Jupiter-gravity-assist Uranus orbiter",
        ),
        citation="Jarmak, Brinckerhoff et al., 'QUEST: A New Frontiers "
        "Uranus orbiter mission concept study,' Acta Astronautica 170:6-26 "
        "(2020); ADS 2020AcAau.170....6J. Polar orbiter, no satellite "
        "tour; not cycler-relevant but anchors the New-Frontiers Uranus "
        "mission concept literature.",
        doi="10.1016/j.actaastro.2020.01.030",
    ),
    CorpusAnchor(
        name="UOP Decadal Flagship Uranus Orbiter & Probe + Aerocapture variants (2022-2025)",
        primary="Uranus",
        body_set=frozenset({"Miranda", "Ariel", "Umbriel", "Titania", "Oberon"}),
        # #350: Citation explicit: 'One-shot insertion tour, NOT periodic
        # cycler.' Decadal-class mga_tour mission concept.
        topology_label=frozenset({"mga-tour"}),
        authors=("Cohen", "Simon", "Saikia", "Hofstadter"),
        keywords=(
            "Uranus Orbiter and Probe UOP Flagship",
            "UOP Decadal trajectory satellite tour",
            "Uranus aerocapture flagship mission",
            "Uranus equatorial moon tour 4.5 year",
        ),
        citation="UOP Decadal Mission Concept Study (2023); Saikia et al. "
        "'A Flagship-class Uranus Orbiter and Probe mission concept using "
        "aerocapture,' Acta Astronautica (2022) DOI 10.1016/"
        "j.actaastro.2022.10.026; 'Uranus Orbiter and Probe: Mission "
        "Challenges and Concept Updates Since the OWL Decadal Survey,' "
        "Planet. Sci. J. 6:ae680c (2025) DOI 10.3847/PSJ/ae680c; Simon et "
        "al. Flagship Science-Driven Tour Design Community Input Poll, "
        "arXiv:2505.05514 (2025). Repeated-Titania-flyby inclination "
        "reduction + 4.5-year equatorial moon tour. One-shot insertion "
        "tour, NOT periodic cycler.",
        doi="10.3847/PSJ/ae680c",
    ),
    # -----------------------------------------------------------------------
    # #334 — BCR4BP system-swap literature anchors (Phase 4 Part D).
    #
    # The system-swap sweep characterised L1 Lyapunov mu_sun-continuation in
    # the Sun-Saturn-Titan, Sun-Mars-Phobos, Sun-Neptune-Triton, and
    # Sun-Pluto-Charon BCR4BP triples (data/scan_334_bcr4bp_system_swap.jsonl,
    # docs/notes/2026-06-17-334-bcr4bp-system-swap.md). The anchors below pin
    # the published CR3BP / mission-design literature for those systems so a
    # FUTURE candidate landing in any of these regimes has its structural
    # fingerprint surfaced deterministically before any novelty claim. None of
    # these papers studies the specific mu_sun-continuation that #334 does
    # (which is OUR computation, not catalogue-promotable per discipline);
    # they cover the surrounding moon-tour / CR3BP literature that the
    # corpus matcher should bump into. Pluto-Charon already has Persephone /
    # Game-Changer / Showalter-Hamilton / Brozovic / arXiv:2510.13479 anchors
    # above; not duplicated here.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Brinckerhoff-Lo-Marsden Saturn-Titan CR3BP libration-point orbits",
        primary="Saturn",
        body_set=frozenset({"Titan"}),
        # #350: Halo / libration-point orbit family at Saturn-Titan L1/L2/L3.
        topology_label=frozenset({"halo"}),
        authors=("Brinckerhoff", "Lo", "Marsden", "Howell"),
        keywords=(
            "Saturn-Titan CR3BP halo orbit",
            "Saturn-Titan libration-point family",
            "Titan moon CR3BP periodic orbit",
        ),
        citation="Brinckerhoff & Howell, 'Three-dimensional periodic orbits in "
        "the Saturn-Titan CR3BP' (AAS 09-129 family); Davis & Howell, halo / "
        "Lyapunov families at Saturn-Titan, J. Astronaut. Sci. (2011-2018 series). "
        "Anchors the published Saturn-Titan CR3BP periodic-orbit corpus that any "
        "#334 BCR4BP Saturn-Titan candidate must be compared against.",
        doi=None,
    ),
    CorpusAnchor(
        name="Cassini-Huygens Saturn-Titan satellite tour design",
        primary="Saturn",
        # #360 (2026-06-17): body_set extended to include Tethys + Mimas per
        # Wolf-Smith 1995 Table 2 verbatim (Tethys 21 targeted; Tethys 46N/49N
        # nontargeted; Mimas 31N nontargeted). Phoebe is NOT included: it is
        # mentioned in §3 for a pre-SOI approach flyby but does NOT appear in
        # the 1995 sample tour Table 2.
        body_set=frozenset({"Titan", "Enceladus", "Rhea", "Dione", "Iapetus", "Tethys", "Mimas"}),
        # #349: topology_label is decisive per Strange-Russell-Buffington
        # AAS-07-277 ("graphical method for the design of transfers between
        # the same gravity-assist body... used with great success in the
        # Cassini extended mission design") + Yam-Davis-Longuski-Howell-
        # Buffington JSR 2009 ("successive Titan flybys... Tisserand graphs").
        # Cassini's tour is fundamentally a SAME-BODY Titan-pump + Tisserand-
        # graph multi-body tour, NOT a repeated-moon (k1, k2) cycler. Rhea /
        # Dione / Enceladus / Iapetus were single-visit science targets
        # reached during the Titan-pump phase, NOT repeating tour members
        # in the cycler sense. A candidate signature with
        # topology_label={"repeated-moon"} is therefore correctly EXCLUDED
        # by this anchor (disjoint sets), unblocking #344 Phase 2 Stages
        # B-E for the Titan-Rhea-Titan (1, 1) candidate.
        # #360 (2026-06-17): topology_label retained as {pump-tour, mga-tour};
        # Wolf-Smith 1995 confirms both (sample tour has 33 Titan pump flybys
        # + 5 single-visit other-moon flybys, structurally a pump-tour wrapped
        # in an mga-tour).
        topology_label=frozenset({"pump-tour", "mga-tour"}),
        # #360 (2026-06-17): authors tuple extended with Wolf + Smith (the
        # 1995 pre-launch base-tour-design authors). Strange/Russell/Buffington/
        # Yam/Davis/Longuski retained as post-launch refinement authors.
        authors=(
            "Wolf",
            "Smith",
            "Strange",
            "Russell",
            "Buffington",
            "Yam",
            "Davis",
            "Longuski",
        ),
        keywords=(
            "Cassini Saturn tour design",
            "Titan same-body V-infinity leveraging pump tour",
            "Tisserand graph Cassini multi-body tour",
            "Cassini Equinox Solstice mission Titan tour",
        ),
        # #360 (2026-06-17): Wolf-Smith 1995 is now the LEADING pre-launch
        # citation (the canonical pre-launch base-tour reference: 63 orbits,
        # 38 targeted flybys, 33 of Titan + 5 of (Enceladus, Tethys, Dione,
        # Rhea, Iapetus); Wolf-Smith 1995 sample tour also includes a Mimas
        # 31N nontargeted flyby per Table 2). Strange/Yam/Valerino retained
        # as post-launch refinement references. See deep-read digest
        # docs/notes/2026-06-17-digest-wolf-smith-1995-cassini.md.
        citation="Wolf & Smith, 'Design of the Cassini Tour Trajectory in the "
        "Saturnian System,' Control Engineering Practice 3(11):1611-1619 "
        "(1995) DOI 10.1016/0967-0661(95)00172-7 -- pre-launch sample tour: "
        "63 orbits, 38 targeted flybys (33 of Titan + 5 of Enceladus/Tethys/"
        "Dione/Rhea/Iapetus); tour ID not given by paper (project-internal "
        "designation only); Strange, Russell & Buffington, 'Mapping the "
        "V-infinity globe' (AAS 07-277, JPL/Caltech, 2007) -- same-body "
        "Titan-pump method used in Cassini extended mission; Yam, Davis, "
        "Longuski, Howell & Buffington, 'Saturn Impact Trajectories for "
        "Cassini End-of-Mission,' JSR DOI 10.2514/1.38760 (2009) -- "
        "successive Titan flybys + Tisserand graphs for Saturn impact; "
        "Valerino, 'Updating the Reference Trajectory for the Cassini "
        "Solstice Mission,' SpaceOps 2014 DOI 10.2514/6.2014-1880 -- "
        "trajectory-update process for the Titan-flyby tour pattern.",
        # #360 (2026-06-17): DOI updated to lead with Wolf-Smith 1995
        # (the pre-launch canonical reference).
        doi="10.1016/0967-0661(95)00172-7",
    ),
    CorpusAnchor(
        name="Wallace Mars-Phobos CR3BP rendezvous trajectory (NASA TM)",
        primary="Mars",
        body_set=frozenset({"Phobos", "Deimos"}),
        authors=("Wallace", "Sims", "Bell"),
        keywords=(
            "Mars-Phobos CR3BP rendezvous",
            "Phobos sample return trajectory",
            "Mars-Phobos libration orbit",
        ),
        citation="Wallace et al., 'Mission concepts for Mars-Phobos exploration' "
        "(JPL / AAS context, 2000-2018); Genova et al., Mars Express extended "
        "Phobos-Deimos flyby campaign. The Sun-Mars-Phobos BCR4BP family in #334 "
        "is the dynamical idealisation of the Phobos-rendezvous design space. "
        "(GAP: no peer-reviewed source we have located studies the Sun-Mars-"
        "Phobos BCR4BP specifically; the published Phobos work is patched-conic "
        "and ephemeris-shooting, not bicircular-restricted-4-body.)",
        doi=None,
    ),
    CorpusAnchor(
        name="Voyager 2 Triton encounter + Trident / Triton-Hopper concept tour",
        primary="Neptune",
        body_set=frozenset({"Triton", "Proteus"}),
        authors=("Stone", "Miner", "Prockter", "Pappalardo"),
        keywords=(
            "Voyager Neptune Triton encounter trajectory",
            "Trident Neptune Triton mission concept",
            "Triton flyby Discovery-class concept",
            "Neptune-Triton CR3BP retrograde orbit",
        ),
        citation="Stone & Miner, 'The Voyager 2 Encounter with Neptune' Science "
        "246:1417 (1989); Prockter et al., 'Trident: a Discovery-class mission to "
        "Triton' (PSJ 2021). The Sun-Neptune-Triton BCR4BP family in #334 is the "
        "planar prograde idealisation of a system whose real Triton orbit is "
        "RETROGRADE + inclined -- carry this caveat into any matching claim. "
        "(GAP: no peer-reviewed Sun-Neptune-Triton bicircular-restricted-4-body "
        "study identified.)",
        doi=None,
    ),
    # -----------------------------------------------------------------------
    # #314 / #403 — Heteroclinic-cycle / Oterma literature anchors.
    # -----------------------------------------------------------------------
    CorpusAnchor(
        name="Koon-Lo-Marsden-Ross dynamical systems and space mission design (2000-2011)",
        primary="Sun",
        body_set=frozenset({"Jupiter"}),
        authors=("Koon", "Lo", "Marsden", "Ross"),
        keywords=(
            "Oterma heteroclinic connection",
            "resonance transition celestial mechanics",
            "Sun-Jupiter-Oterma PCR3BP",
            "dynamical systems space mission design",
        ),
        citation="Koon, W.S., Lo, M.W., Marsden, J.E. & Ross, S.D., 'Dynamical "
        "Systems, the Three-Body Problem and Space Mission Design' (2006/2011); "
        "Chaos 10(2):427-469 (2000) DOI 10.1063/1.166509. Foundational "
        "connection-only conceptual framework for Oterma.",
        doi="10.1063/1.166509",
    ),
    CorpusAnchor(
        name="Wilczak-Zgliczyński Sun-Jupiter-Oterma heteroclinic cycle proof",
        primary="Sun",
        body_set=frozenset({"Jupiter"}),
        authors=("Wilczak", "Zgliczyński", "Zgliczynski"),
        keywords=(
            "Oterma heteroclinic cycle",
            "computer-assisted proof",
            "L1 L2 Lyapunov heteroclinic",
            "Sun-Jupiter-Oterma PCR3BP",
        ),
        citation="Wilczak, D. & Zgliczyński, P., 'Heteroclinic Connections "
        "between Periodic Orbits in Planar Restricted Circular Three-Body "
        "Problem', Comm. Math. Phys. (Part I: arXiv:math/0201278; Part II: "
        "DOI 10.1007/s00220-005-1374-x). The sourced golden dataset for "
        "closed L1<->L2 Lyapunov heteroclinic cycles.",
        doi="10.1007/s00220-005-1374-x",
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
        if not (overlap and seq_set <= anchor.body_set):
            continue
        # Optional CR3BP period-band filter (#301): if BOTH sides declare a
        # ``period_band_tu``, drop the anchor when the bands are disjoint -- a
        # candidate at a period-multiplied sub-family is structurally out-of-
        # scope for an anchor that catalogues only the low-integer base family.
        if sig.period_band_tu is not None and anchor.period_band_tu is not None:
            c_min, c_max = sig.period_band_tu
            a_min, a_max = anchor.period_band_tu
            if c_max < a_min or c_min > a_max:
                continue
        # Optional topology-label filter (#349): if BOTH sides declare a
        # non-empty ``topology_label`` set, drop the anchor when the sets are
        # disjoint -- a (k1, k2) repeated-moon candidate is not the same family
        # as a Titan-pump tour anchor even when they share a body subset.
        # Empty on either side falls through to the historical body-set-only
        # match, preserving prior behaviour for un-annotated anchors.
        if (
            sig.topology_label
            and anchor.topology_label
            and not (sig.topology_label & anchor.topology_label)
        ):
            continue
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
