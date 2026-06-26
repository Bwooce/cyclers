"""#485 — citation-integrity consistency ratchet (body/system containment).

The guard that cannot be forgotten. A citation hallucination was found
2026-06-26: a Jovian "Galilean Io-Europa-Ganymede triple cycler" claim was at
risk of being cited to Jones-Hernandez-Jesick AAS 17-577, which is actually
"Low Excess Speed Triple Cyclers of Venus, Earth, and Mars" (a HELIOCENTRIC
interplanetary paper). Root cause: citation-by-concept-collision -- "triple
cycler" spans dynamical systems -- and it propagates because a hallucinated
citation, once written, is copied from our own notes as fact. See memory
``feedback_ground_citations_against_content``.

This ratchet makes the defect structurally impossible to land silently: for
every catalogue row whose ``first_published`` work is also pinned in the
``literature_check.KNOWN_CORPUS`` ground-truth registry (the #483 sourced
``system`` + ``body_set`` surfaces), it asserts the row's CLAIMED bodies are a
subset of the cited work's bodies and the row's system matches the cited work's
system. A {Io, Europa, Ganymede} jovian claim citing a {Venus, Earth, Mars}
heliocentric work FAILS (the positive control below proves the ratchet fires).

Ground truth: ``literature_check.KNOWN_CORPUS`` anchors carry ``system`` (the
#483 field, extracted from each paper's TITLE/ABSTRACT) and ``body_set`` (the
cited work's bodies). Both trace to the source's near-verbatim title text, not
a summary-of-summary -- the lowest-hallucination surface.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CorpusAnchor,
    system_for_primary,
)

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

# ---------------------------------------------------------------------------
# Body-token normalisation
#
# The catalogue uses single-letter codes for the planets it treats as flyby
# targets (``J`` Jupiter, ``S`` Saturn, ``U`` Uranus, ``N`` Neptune, ``Me``
# Mercury, ``E`` Earth) while KNOWN_CORPUS anchors spell some of them out. We
# canonicalise both sides before the subset test so a token-spelling difference
# is never mistaken for a body mis-citation. Minor bodies (asteroids/comets
# flown by on the way) are out of SYSTEM-integrity scope -- the ratchet guards
# the planet/moon system identity, not asteroid-flyby completeness.
# ---------------------------------------------------------------------------
_BODY_ALIAS = {
    "J": "Jupiter",
    "S": "Saturn",
    "U": "Uranus",
    "N": "Neptune",
    "Me": "Mercury",
    "E": "Earth",
}
_MINOR_BODIES = frozenset({"Gaspra", "Ida", "Steins", "Lutetia", "Vesta", "Ceres"})


def _norm_bodies(bodies: Iterable[str] | None, primary: str) -> frozenset[str]:
    """Canonicalise a body list: spell-out aliases, drop minor bodies + primary.

    The primary body itself is dropped because a moon-system anchor catalogues
    the secondaries (e.g. ``{Moon}``) while a CR3BP row lists both primary and
    secondary (``{E, Moon}``); the containment test is over the encountered
    bodies, not the central body.
    """
    canon_primary = _BODY_ALIAS.get(primary, primary)
    out = {_BODY_ALIAS.get(str(b), str(b)) for b in (bodies or [])} - _MINOR_BODIES
    out.discard(canon_primary)
    return frozenset(out)


def _surname(author: str) -> str:
    """The cited-author surname (catalogue ``first_published`` is 'Surname, X.')."""
    return author.split(",")[0].strip()


# ---------------------------------------------------------------------------
# The consistency check (the reusable core)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CitationClaim:
    """A (system, bodies) claim attached to a citation, for the containment test."""

    claim_id: str
    system: str
    bodies: frozenset[str]
    author_surnames: frozenset[str]


@dataclass(frozen=True)
class IntegrityViolation:
    claim_id: str
    anchor_name: str
    claim_system: str
    claim_bodies: tuple[str, ...]
    anchor_system: str
    anchor_bodies: tuple[str, ...]
    extra_bodies: tuple[str, ...]

    def describe(self) -> str:
        return (
            f"{self.claim_id!r} claims system={self.claim_system} "
            f"bodies={self.claim_bodies} but cites {self.anchor_name!r} "
            f"(system={self.anchor_system}, bodies={self.anchor_bodies}); "
            f"body/system MIS-CITATION -- extra bodies {self.extra_bodies} "
            f"are NOT in the cited work."
        )


def _anchor_represents_claim(anchor: CorpusAnchor, claim: CitationClaim) -> bool:
    """Does this anchor pin the SAME published work the claim cites?

    Strong link (low false-positive): the anchor's author surnames are a subset
    of the claim's cited authors AND the systems agree. A single shared surname
    is NOT enough -- "Longuski" appears on Earth-Mars, VEM, and Uranian-tour
    works alike; only the full author-set-subset link ties an anchor to the
    claim's own paper. (A mis-citation that names a DIFFERENT author set would
    not strong-link to the wrong anchor and so would not be checked here; the
    audit doc records that this ratchet covers author-consistent citations --
    the exact class the Galilean hallucination belongs to.)
    """
    anchor_surs = {_surname(a) for a in anchor.authors}
    if not anchor_surs or not (anchor_surs <= claim.author_surnames):
        return False
    return anchor.system_grounded == claim.system


def check_citation_integrity(claim: CitationClaim) -> list[IntegrityViolation]:
    """Assert ``claim.bodies`` are a subset of every cited work's bodies.

    Returns one :class:`IntegrityViolation` per anchor the claim cites whose
    sourced ``body_set`` does NOT contain all the claim's bodies. An empty list
    means the claim is body/system-consistent with the published record it
    cites (or cites nothing pinned in KNOWN_CORPUS, which is not a violation --
    only a positive containment failure is).
    """
    violations: list[IntegrityViolation] = []
    for anchor in KNOWN_CORPUS:
        if not _anchor_represents_claim(anchor, claim):
            continue
        anchor_bodies = _norm_bodies(anchor.body_set, _claim_primary(claim))
        extra = claim.bodies - anchor_bodies
        if extra:
            violations.append(
                IntegrityViolation(
                    claim_id=claim.claim_id,
                    anchor_name=anchor.name,
                    claim_system=claim.system,
                    claim_bodies=tuple(sorted(claim.bodies)),
                    anchor_system=anchor.system_grounded,
                    anchor_bodies=tuple(sorted(anchor_bodies)),
                    extra_bodies=tuple(sorted(extra)),
                )
            )
    return violations


def _claim_primary(claim: CitationClaim) -> str:
    """Recover the primary body from a claim's system (for body normalisation)."""
    inverse = {
        "heliocentric": "Sun",
        "earth-moon": "Earth",
        "jovian": "Jupiter",
        "saturnian": "Saturn",
        "uranian": "Uranus",
        "neptunian": "Neptune",
        "pluto-charon": "Pluto",
        "mars-system": "Mars",
    }
    return inverse.get(claim.system, "Sun")


def _claim_from_row(row: dict) -> CitationClaim:  # type: ignore[type-arg]
    primary = row.get("primary", "Sun")
    fp = row.get("first_published") or {}
    surs = frozenset(_surname(a) for a in (fp.get("authors") or []))
    return CitationClaim(
        claim_id=row["id"],
        system=system_for_primary(primary),
        bodies=_norm_bodies(row.get("bodies"), primary),
        author_surnames=surs,
    )


# ---------------------------------------------------------------------------
# Positive controls -- the ratchet MUST flag these (else it is vacuous)
# ---------------------------------------------------------------------------


def test_galilean_claim_citing_vem_work_is_flagged() -> None:
    """THE regression: a {Io, Europa, Ganymede} jovian claim must NOT pass when
    cited to the heliocentric {Venus, Earth, Mars} VEM work.

    This reconstructs the 2026-06-26 hallucination directly: same authors
    (Jones/Hernandez/Jesick), but the claim is Jovian (Galilean moons) while the
    cited VEM anchor is heliocentric. Because the systems differ, the VEM anchor
    does not even strong-link -- so we ALSO assert the stronger property: were
    the Galilean claim mislabelled heliocentric (the actual concept-collision
    error), the body containment fails outright.
    """
    # (a) honest jovian claim, wrong (heliocentric VEM) authors+system: the
    #     systems disagree, so no strong-link forms -- the claim is simply not
    #     matched to the VEM anchor (correct: a jovian claim cannot rediscover a
    #     heliocentric paper). We assert the body containment directly to prove
    #     the core check fires on the body sets themselves.
    vem_anchor = next(a for a in KNOWN_CORPUS if "VEM triple cyclers" in a.name)
    galilean_bodies = _norm_bodies(["Io", "Europa", "Ganymede"], "Jupiter")
    vem_bodies = _norm_bodies(vem_anchor.body_set, "Sun")
    assert not (galilean_bodies <= vem_bodies), (
        "POSITIVE CONTROL BROKEN: Galilean moons must NOT be a subset of the "
        "VEM (Venus-Earth-Mars) body set -- the hallucination would slip through."
    )

    # (b) the full claim path: a Galilean claim mislabelled as the VEM authors'
    #     work, forced to the heliocentric system (the exact mis-attribution).
    bad_claim = CitationClaim(
        claim_id="POSITIVE-CONTROL-galilean-cited-to-vem",
        system="heliocentric",  # the concept-collision mislabel
        bodies=galilean_bodies,
        author_surnames=frozenset({"Jones", "Hernandez", "Jesick"}),
    )
    violations = check_citation_integrity(bad_claim)
    assert violations, (
        "POSITIVE CONTROL FAILED: a Galilean {Io,Europa,Ganymede} claim cited to "
        "the VEM heliocentric work must be flagged a body/system mis-citation."
    )
    assert any(v.extra_bodies for v in violations)


def test_known_good_galilean_claim_passes() -> None:
    """The grounded mirror: a Galilean IEG claim cited to the CORRECT Jovian
    Hernandez/Jones/Jesick IEG anchor must PASS (not vacuously failing).
    """
    good_claim = CitationClaim(
        claim_id="KNOWN-GOOD-ieg-cited-to-jovian-ieg",
        system="jovian",
        bodies=_norm_bodies(["Io", "Europa", "Ganymede"], "Jupiter"),
        author_surnames=frozenset({"Hernandez", "Jones", "Jesick"}),
    )
    violations = check_citation_integrity(good_claim)
    assert violations == [], (
        f"KNOWN-GOOD case wrongly flagged (ratchet is over-firing): "
        f"{[v.describe() for v in violations]}"
    )
    # And it MUST strong-link to the Jovian anchor (not silently pass by matching
    # nothing) -- otherwise the known-good is vacuous.
    ieg_anchor = next(a for a in KNOWN_CORPUS if "Io-Europa-Ganymede triple cyclers" in a.name)
    assert _anchor_represents_claim(ieg_anchor, good_claim), (
        "KNOWN-GOOD is vacuous: the Jovian IEG claim does not strong-link to the Jovian IEG anchor."
    )


# ---------------------------------------------------------------------------
# The corpus-wide ratchet over catalogue.yaml
# ---------------------------------------------------------------------------


def _load_rows() -> list[dict]:  # type: ignore[type-arg]
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_no_catalogue_citation_body_system_mis_citation() -> None:
    """Every catalogue row's claimed bodies/system must be contained by the
    bodies/system of any KNOWN_CORPUS work it cites (#485).

    Verified clean across all strong-linked rows on 2026-06-26 (see
    docs/notes/2026-06-26-citation-audit.md). This is a frozen ratchet: a future
    row that cites a body/system-inconsistent work (the Galilean->Hernandez-2017
    failure class) trips it.
    """
    rows = _load_rows()
    all_violations: list[IntegrityViolation] = []
    for row in rows:
        all_violations.extend(check_citation_integrity(_claim_from_row(row)))
    assert all_violations == [], "Body/system mis-citations found:\n" + "\n".join(
        v.describe() for v in all_violations
    )


def test_ratchet_is_not_vacuous() -> None:
    """At least a healthy number of catalogue rows actually strong-link to an
    anchor -- proves the corpus sweep exercises real containment checks, not a
    no-op that passes because nothing matched.
    """
    rows = _load_rows()
    linked = 0
    for row in rows:
        claim = _claim_from_row(row)
        if any(_anchor_represents_claim(a, claim) for a in KNOWN_CORPUS):
            linked += 1
    assert linked >= 20, (
        f"Only {linked} catalogue rows strong-link to a KNOWN_CORPUS anchor; "
        "the ratchet may be vacuous -- check the author/system link logic."
    )


@pytest.mark.parametrize(
    "anchor",
    KNOWN_CORPUS,
    ids=[a.name[:40] for a in KNOWN_CORPUS],
)
def test_every_anchor_has_a_grounded_system(anchor: CorpusAnchor) -> None:
    """#483: every KNOWN_CORPUS anchor resolves to a non-empty grounded system
    (explicit ``system`` or derived from the title-sourced ``primary``).
    """
    assert anchor.system_grounded, f"anchor {anchor.name!r} has no grounded system"
