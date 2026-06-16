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
        name="Davis-Phillips-McCarthy Saturn tulip-shaped orbits",
        primary="Saturn",
        body_set=frozenset({"Titan", "Enceladus", "Rhea", "Dione"}),
        authors=("Davis", "Phillips", "McCarthy"),
        keywords=(
            "Saturn tulip-shaped orbit",
            "Titan period-multiplying orbit",
            "Saturnian periodic orbit Np petal",
        ),
        citation="Davis, Phillips & McCarthy, 'Tulip-shaped orbits in the "
        "Saturn system,' Acta Astronautica 143:16-28 (2018)",
        doi="10.1016/j.actaastro.2017.11.011",
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
    CorpusAnchor(
        name="Ceriotti GTOC-style MGA chain optimisation (2010)",
        primary="Sun",
        body_set=frozenset({"V", "E", "M", "Jupiter"}),
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
        name="Diehl-Belbruno-Roberts Galileo VEEGA design (1986)",
        primary="Sun",
        body_set=frozenset({"V", "E", "Jupiter"}),
        authors=("Diehl", "Belbruno", "Roberts", "D'Amario"),
        keywords=(
            "VEEGA Venus-Earth-Earth gravity assist",
            "Galileo Jupiter mission trajectory",
            "Earth-Earth gravity assist tour",
        ),
        citation="Diehl, Belbruno & Roberts et al., 'Galileo VEEGA Mission "
        "Design' (1986-1990 JPL / AAS); the canonical mga_tour archetype "
        "(October 1989 launch window once-per-~13yr alignment)",
        doi=None,
    ),
    CorpusAnchor(
        name="Tito-MacCallum 2018 Mars free-return mission design (2013)",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
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
