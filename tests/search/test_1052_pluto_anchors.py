"""#1052 (Pluto part): the Pluto anchor edits and the gate outcomes they are pre-registered to give.

Pre-registration: docs/notes/2026-10-11-1052-pluto-anchor-labels-prereg.md (secs. 5, 7, 8, 10).
Rulings: D1 (b) Howett stays unlabelled, D3 Showalter-Hamilton {ephemeris}, Stern scope "one",
tadpole citation corrected to the arXiv record. Expected outcomes are the pre-registered ones under
rule U / U2 (twobody-gen2-opus, e3fba9f6): an unlabelled anchor can only make a signature
"inconclusive".
"""

from __future__ import annotations

from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)


def _anchor(fragment: str):  # type: ignore[no-untyped-def]
    hits = [a for a in KNOWN_CORPUS if fragment in a.name]
    assert len(hits) == 1, fragment
    return hits[0]


def _status(
    topology: str | None,
    sequence: tuple[str, ...] = ("Charon", "Nix"),
    working: str | None = "one",
    primary: str = "Pluto",
) -> tuple[str, str]:
    sig = CandidateSignature(
        primary=primary,
        sequence=sequence,
        period_k=4,
        vinf_per_encounter_kms=(0.178, 0.130),
        topology_label=frozenset({topology}) if topology else frozenset(),
        working_bodies=working,
        closing_repeat=False,
    )
    r = check_literature(sig, search=offline_corpus_search)
    return r.status, r.citation or ""


def test_anchor_edits() -> None:
    assert _anchor("Howett").topology_label == frozenset()  # D1 (b): left empty on purpose
    assert _anchor("Showalter-Hamilton").topology_label == frozenset({"ephemeris"})
    assert _anchor("Stern/SwRI").working_bodies_scope == "one"
    assert _anchor("Stern/SwRI").topology_label == frozenset({"mga-tour"})
    tadpole = _anchor("tadpole/horseshoe")
    assert tadpole.authors == ("Jbara",)
    assert "Chaotic Dynamics and Zero-Velocity Structures" in tadpole.citation
    assert "2510.13479" in tadpole.citation
    assert tadpole.topology_label == frozenset({"binary-coorbital"})
    assert "UNVERIFIED" in _anchor("Brozovic et al.").citation


def test_c1_c2_c3() -> None:
    # C1: a #998-type signature is not "published" (Howett, unlabelled, caps at inconclusive)
    s1, c1 = _status("repeated-moon")
    assert s1 == "inconclusive" and "Howett" in c1
    # C2: wrong system never reaches a Pluto anchor
    s2, _ = _status("repeated-moon", primary="Jupiter")
    assert s2 == "known-architecture-new-system"
    # C3: the untopologied Pluto signature moves from "published" to "inconclusive" naming Howett
    s3, c3 = _status(None)
    assert s3 == "inconclusive" and "Howett" in c3


def test_diagnostics() -> None:
    assert _status("pump-tour")[0] == "inconclusive"
    assert _status("mga-tour")[0] == "published"  # Stern, scope "one"
    assert "Stern" in _status("mga-tour")[1]
    assert _status("mga-tour", working="two")[0] == "inconclusive"  # Stern excludes "two"
    assert _status("resonant", sequence=("Nix", "Hydra"))[0] != "published"  # Showalter: ephemeris
    assert _status("binary-coorbital", sequence=("Charon",))[0] == "published"
