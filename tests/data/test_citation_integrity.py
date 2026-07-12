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
    anchor_for_key,
    build_citation_registry,
    can_anchor_decision,
    citation_key_exists,
    derive_citation_key,
    resolve_citation_key,
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
    is_family_seed: bool = False
    """True for a catalogue "family seed" placeholder row (id/name convention
    ``"... (family seed)"`` / ``"... (family seed entry)"``): an intentionally
    PAPER-LEVEL summary row (``legs: []``, all-null ``vinf_kms_at_encounters``)
    documenting a whole work pending per-member ingestion, not a specific
    cycler's claim. Per #578: Russell & Strange 2009's Saturnian family seed
    (``russell-strange-2009-saturnian-multimoon-family``) honestly claims all
    6 Saturnian moons the paper's multi-moon Titan-flyby cyclers pass through
    as science targets (digest ``docs/notes/2026-06-30-digest-russell-
    strange-2009-planetary-moon-cyclers.md``), but #578 registered only ONE
    Saturnian anchor (Titan-Enceladus -- R-S's only Table-1-enumerated
    repeated-moon PAIR; #578's own scope explicitly excludes adding a
    per-pair anchor for every science-target moon). A family-seed row's
    author+system correctness is still checked (the actual hallucination-
    prevention core -- the Galilean-vs-VEM regression this ratchet exists
    for is a SYSTEM mismatch, not a body-count mismatch); only the strict
    body-superset check is relaxed, since by the row's own documented intent
    it is a broader-than-any-single-anchor placeholder, not a point claim.
    """


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
    """Assert ``claim.bodies`` are covered by the strong-linked cited work(s).

    Returns a non-empty list iff the claim names a body no strong-linked
    anchor attests. Coverage is checked against the UNION of every anchor
    that strong-links to the claim (same author-surname-subset + system),
    not each anchor in isolation: #578 established the pattern of registering
    SEVERAL narrow per-pair ``CorpusAnchor``s for a single paper (Russell &
    Strange 2009's Jovian per-pair anchors -- deliberately NOT one body_set
    union, so a candidate cycler search never falsely collides an
    unenumerated pair like Io-Callisto against the whole family; see
    ``search/literature_check.py``'s ``_candidate_anchors`` docstring). A
    catalogue "family seed" row legitimately describing the WHOLE paper
    (e.g. ``russell-strange-2009-jovian-multimoon-family``, claiming all 4
    Galilean moons) is correctly covered by the union of its several
    per-pair anchors even though no SINGLE one of them spans the full body
    set -- that is not a mis-citation, it is the intended multi-anchor
    pattern. Only a body genuinely absent from every strong-linked anchor's
    body_set is a real violation. An empty list also covers "cites nothing
    pinned in KNOWN_CORPUS", which is not a violation -- only a positive
    containment failure is.
    """
    primary = _claim_primary(claim)
    matched = [anchor for anchor in KNOWN_CORPUS if _anchor_represents_claim(anchor, claim)]
    if not matched:
        return []
    if claim.is_family_seed:
        # A family-seed row's author+system correctness is exactly what
        # ``matched`` (non-empty) already proves; its body list is an
        # intentionally paper-level placeholder, not a specific claim to
        # hold to strict per-anchor(-union) containment. See
        # CitationClaim.is_family_seed.
        return []
    anchor_bodies_by_name = {
        anchor.name: _norm_bodies(anchor.body_set, primary) for anchor in matched
    }
    union_bodies: frozenset[str] = frozenset().union(*anchor_bodies_by_name.values())
    uncovered = claim.bodies - union_bodies
    if not uncovered:
        return []
    violations: list[IntegrityViolation] = []
    for anchor in matched:
        violations.append(
            IntegrityViolation(
                claim_id=claim.claim_id,
                anchor_name=anchor.name,
                claim_system=claim.system,
                claim_bodies=tuple(sorted(claim.bodies)),
                anchor_system=anchor.system_grounded,
                anchor_bodies=tuple(sorted(anchor_bodies_by_name[anchor.name])),
                extra_bodies=tuple(sorted(uncovered)),
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
        is_family_seed="family seed" in str(row.get("name", "")).lower(),
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


# ---------------------------------------------------------------------------
# #484 — structured citation-key registry: close the WRONG-AUTHOR gap.
#
# #485 (above) catches AUTHOR-CONSISTENT mis-citations (right author, wrong
# system/bodies -- the Galilean-vs-VEM class). A mis-citation that ALSO names
# the WRONG authors does not strong-link to any anchor, so #485 never sees it.
# #484 closes that: a citation resolves to a grounded work THROUGH a stable
# registry KEY, and a key that does not resolve (a wrong-author / fabricated
# reference) is FLAGGED.
# ---------------------------------------------------------------------------


def test_citation_registry_builds_with_unique_keys() -> None:
    """Every KNOWN_CORPUS anchor is registry-addressable by a UNIQUE key.

    ``build_citation_registry`` raises on a duplicate key, so a green build is
    itself the uniqueness assertion; we also assert the registry covers every
    anchor (no anchor is unreachable by key).
    """
    registry = build_citation_registry()
    assert len(registry) == len(KNOWN_CORPUS), (
        "registry key count does not match anchor count -- a key collision "
        "silently dropped an anchor."
    )
    for anchor in KNOWN_CORPUS:
        assert anchor.key_resolved in registry
        assert resolve_citation_key(anchor.key_resolved) is anchor


def test_ground_truthed_anchors_carry_explicit_stable_keys() -> None:
    """The two works we ground-truthed this session resolve by their explicit,
    stable, human-meaningful keys (not just a derived slug).
    """
    ieg = resolve_citation_key("hernandez-2017-ieg-608")
    assert ieg.body_set == frozenset({"Io", "Europa", "Ganymede"})
    assert ieg.system_grounded == "jovian"

    vem = resolve_citation_key("jones-2017-vem-577")
    assert vem.body_set == frozenset({"V", "E", "M"})
    assert vem.system_grounded == "heliocentric"

    # The two same-author sibling papers must NOT collide on one key (the exact
    # concept-collision the registry disambiguates).
    assert ieg.key_resolved != vem.key_resolved


def test_wrong_author_fabricated_key_is_flagged() -> None:
    """POSITIVE CONTROL (#484): a fabricated / wrong-author citation key does
    NOT resolve and is flagged -- the wrong-author gap #485 cannot reach.

    A mis-citation that names authors absent from the grounded record produces
    a key that is not in the registry; resolving it must raise, and the boolean
    existence check must report False. (This is the failure #485 misses: it only
    checks author-consistent citations, so a wrong-author reference slips past
    it -- here it is caught at the key.)
    """
    fabricated_key = "smith-2099-fabricated-galilean-vem-cycler"
    assert not citation_key_exists(fabricated_key)
    with pytest.raises(KeyError):
        resolve_citation_key(fabricated_key)

    # And a plausible-looking but WRONG-AUTHOR variant of a real work (e.g. the
    # IEG paper mis-attributed to "Longuski") derives a key that does not match
    # the grounded IEG key -- so it cannot silently bind to the real anchor.
    wrong_author_key = derive_citation_key(
        ("Longuski", "Byrnes"), "One Class of Io-Europa-Ganymede Triple Cyclers"
    )
    assert not citation_key_exists(wrong_author_key), (
        "a wrong-author variant of the IEG paper resolved to a registry key -- "
        "the wrong-author gap is NOT closed."
    )
    assert wrong_author_key != "hernandez-2017-ieg-608"


def test_every_catalogue_strong_link_resolves_to_a_registry_key() -> None:
    """Every catalogue row that strong-links to a KNOWN_CORPUS anchor resolves
    that anchor to a registry key -- the keyed surface the ratchet now guards.

    This wires the #484 registry into the corpus sweep: the same author-set +
    system strong-link #485 uses to bind a row to its cited work also yields a
    resolvable registry key for that work. A future anchor added without a
    resolvable key (or a duplicate that collides) trips this.
    """
    rows = _load_rows()
    registry = build_citation_registry()
    linked_keys: set[str] = set()
    for row in rows:
        claim = _claim_from_row(row)
        for anchor in KNOWN_CORPUS:
            if _anchor_represents_claim(anchor, claim):
                key = anchor.key_resolved
                assert key in registry, (
                    f"row {row['id']!r} strong-links to anchor {anchor.name!r} "
                    f"whose key {key!r} does not resolve in the registry."
                )
                linked_keys.add(key)
    assert linked_keys, "no catalogue row strong-linked to any registry key (vacuous)."


# ---------------------------------------------------------------------------
# #486 — provenance tags: an inherited-unverified citation cannot anchor a
# decision until ground-truthed (the #480 false-erratum failure mode).
# ---------------------------------------------------------------------------


def test_new_anchors_default_to_inherited_unverified() -> None:
    """Provenance defaults conservatively: an anchor that does not explicitly
    declare ``verified-against-source`` is ``inherited-unverified``.

    A constructed anchor with no provenance set must NOT be decision-grade --
    the default cannot silently certify an ungrounded citation.
    """
    fresh = CorpusAnchor(
        name="Unverified test anchor",
        primary="Sun",
        body_set=frozenset({"E", "M"}),
        authors=("Nobody",),
        keywords=("test",),
        citation="placeholder",
        doi=None,
    )
    assert fresh.provenance == "inherited-unverified"
    assert not can_anchor_decision(fresh)


def test_ground_truthed_citations_are_verified_against_source() -> None:
    """The two citations we ground-truthed this session are decision-grade."""
    for key in ("hernandez-2017-ieg-608", "jones-2017-vem-577"):
        anchor = resolve_citation_key(key)
        assert anchor.provenance == "verified-against-source", (
            f"{key!r} should be verified-against-source (ground-truthed this session)."
        )
        assert can_anchor_decision(anchor)
        # anchor_for_key combines key resolution + the provenance gate: a
        # verified, resolvable key passes it.
        assert anchor_for_key(key) is anchor


def test_inherited_unverified_citation_cannot_anchor_a_decision() -> None:
    """THE #486 gate: an ``inherited-unverified`` citation cannot anchor a
    spec / validation-gate / catalogue-row promotion until ground-truthed.

    ``anchor_for_key`` resolves the key (#484) AND enforces the provenance gate
    (#486): an inherited-unverified work raises ``PermissionError`` even though
    its key resolves cleanly. This is the exact failure mode that produced the
    #480 false erratum -- a citation present-but-ungrounded backing a decision.
    """
    # Find a real, resolvable, but still-inherited anchor in the corpus.
    inherited = next(a for a in KNOWN_CORPUS if a.provenance == "inherited-unverified")
    assert citation_key_exists(inherited.key_resolved), "test anchor must resolve by key"
    assert not can_anchor_decision(inherited)
    with pytest.raises(PermissionError):
        anchor_for_key(inherited.key_resolved)


def test_every_anchor_provenance_is_a_known_value() -> None:
    """Provenance is one of the two defined tags -- no free-text drift."""
    allowed = {"verified-against-source", "inherited-unverified"}
    for anchor in KNOWN_CORPUS:
        assert anchor.provenance in allowed, (
            f"anchor {anchor.name!r} has unknown provenance {anchor.provenance!r}"
        )
