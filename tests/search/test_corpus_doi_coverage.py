"""Static DOI-coverage ratchet for `search/literature_check.py`'s KNOWN_CORPUS (#578).

#577 found `KNOWN_CORPUS` was missing its single most direct Jovian prior --
Russell & Strange 2009 (DOI 10.2514/1.36610) -- despite `data/catalogue.yaml`
already carrying 30 individual member rows from that exact paper (task #491
ingestion). That gap silently produced a FALSE structural "novelty clear" for
the Io-Callisto pair in #577's literature-adjacency check: nothing in
KNOWN_CORPUS could ever flag a future Io-Callisto-shaped candidate as
known-published, even though the project's OWN catalogue already held dozens
of that same paper's members.

This is a cheap, mechanical, no-schema-mapping-risk structural fix: walk every
`source: literature` catalogue row, collect every DOI it cites (its own
`first_published.doi` AND every `corroborating_sources[].doi`), and assert
each one is either registered as a `KNOWN_CORPUS` anchor DOI or explicitly
allowlisted here with a one-line justification. A DOI landing in neither
bucket is exactly the #577 failure mode: a paper the catalogue already treats
as authoritative that the literature-novelty matcher cannot see.

This does NOT auto-generate KNOWN_CORPUS from the catalogue (explicitly
rejected -- schema-mapping + circularity risk, per the answer already given to
the user's "will KNOWN_CORPUS get auto-generated?" question) and does NOT
build a runtime warn-only cross-reference layer (out of #578's scope, sketch-
only if ever wanted). It is a static, CI-enforced accounting ratchet: add a
new anchor OR add a one-line-justified allowlist entry, but do not let a
literature DOI silently fall through both.
"""

from __future__ import annotations

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.literature_check import KNOWN_CORPUS

# ---------------------------------------------------------------------------
# Explicit allowlist: literature DOIs NOT (yet) covered by a KNOWN_CORPUS
# anchor, with a one-line reason each. Seeded 2026-07-12 (#578) from the full
# `data/catalogue.yaml` scan below -- 18 of 27 distinct literature DOIs at
# HEAD (after this task's own 4 new Russell-Strange 2009 anchors landed;
# Fable's original audit counted 19/27 BEFORE those anchors existed, which
# this test independently reproduced before adding them). Most are companion/
# follow-up papers to an author already anchored (Russell-Ocampo variants,
# McConaghy/Byrnes EM-cycler follow-ups) or mission-description/infrastructure
# citations that were never candidates for a structural cycler anchor in the
# first place. Removing an entry here (because a new anchor now covers it) is
# always safe; the test will simply pass with fewer allowlist rows.
# ---------------------------------------------------------------------------

DOI_ALLOWLIST: dict[str, str] = {
    # Mission-description / infrastructure citations -- not cycler-methods
    # papers, structurally never a KNOWN_CORPUS anchor candidate.
    "10.1007/bf00200846": (
        "Kohlhase & Penzo 1977, 'Voyager Mission Description' -- mission "
        "description reference for the Voyager Grand Tour rows, not a "
        "cycler-methods paper."
    ),
    "10.1016/0032-0633(95)00107-7": (
        "Acton 1996, NAIF SPICE ancillary-data-services paper -- ephemeris "
        "infrastructure citation, not a cycler paper."
    ),
    "10.1016/0094-5765(75)90012-0": (
        "Giberson & Cunningham 1975, 'Mariner 10 mission to Venus and "
        "Mercury' -- mission description reference, not a cycler paper."
    ),
    # Russell-Ocampo Earth-Mars SnLm family: companion/variant papers to the
    # already-anchored 2004 JGCD paper (anchor doi 10.2514/1.1909).
    "10.2514/1.1011": (
        "Russell & Ocampo, 'A Systematic Method for Constructing Earth-Mars "
        "Cyclers Using Free-Return Trajectories' -- same title as the "
        "anchored 2004 JGCD paper (doi 10.2514/1.1909); a handful of "
        "catalogue corroborating_sources rows cite this DIFFERENT doi "
        "string for what appears to be the identical paper (likely a "
        "catalogue-side DOI typo/variant, flagged here rather than silently "
        "assumed -- catalogue.yaml edits are out of #578's scope)."
    ),
    "10.2514/1.5571": (
        "Russell & Ocampo 2005, 'Geometric Analysis of Free-Return "
        "Trajectories Following a Gravity-Assisted Flyby' -- the large "
        "SnLm geometric-analysis catalogue paper (companion to the anchored "
        "2004 paper), not yet separately anchored."
    ),
    "10.2514/1.13652": (
        "Russell & Ocampo 2006, 'Optimization of a Broad Class of Ephemeris "
        "Model Earth-Mars Cyclers' -- real-ephemeris follow-up to the "
        "anchored 2004 paper, not yet separately anchored."
    ),
    # McConaghy/Byrnes/Longuski Earth-Mars cycler follow-ups: named in the
    # Aldrin / Russell-Ocampo anchors' CITATION TEXT but not DOI-matchable
    # (those anchors' registered `doi` field points at a different paper in
    # the same family) -- a real, if narrow, structural gap.
    "10.2514/6.2002-4420": (
        "McConaghy/Longuski/Byrnes, AIAA-2002-4420, 'Analysis of a Broad "
        "Class of Earth-Mars Cycler Trajectories' -- named in the Aldrin "
        "anchor's citation text but the anchor's registered doi is "
        "10.2514/3.25519 (the JSR paper), so this DOI itself doesn't match."
    ),
    "10.2514/6.2002-4423": (
        "Byrnes/McConaghy/Longuski, 'Analysis of Various Two Synodic Period "
        "Earth-Mars Cycler Trajectories' -- named in the Russell-Ocampo "
        "anchor's citation text but that anchor's registered doi is "
        "10.2514/1.1909, so this DOI itself doesn't match."
    ),
    "10.2514/1.11610": (
        "Chen/Landau/McConaghy/Okutsu/Longuski/Aldrin 2005, 'Powered "
        "Earth-Mars Cycler with Three-Synodic-Period Repeat Time' -- "
        "powered/low-thrust EM cycler variant, not yet anchored."
    ),
    "10.2514/1.15215": (
        "McConaghy/Landau/Yam/Longuski 2006, 'Notable Two-Synodic-Period "
        "Earth-Mars Cycler' -- S2L1-class EM cycler paper, not yet anchored."
    ),
    "10.2514/6.1986-2009": (
        "Friedlander/Niehoff/Byrnes/Longuski 1986, 'Circulating "
        "Transportation Orbits Between Earth and Mars' -- VISIT-cycler-"
        "class paper corroborating the Aldrin/Niehoff EM cycler lineage, "
        "not yet separately anchored."
    ),
    "10.2514/6.2012-4746": (
        "Rogers/Hughes/Longuski/Aldrin 2012, 'Preliminary Analysis of "
        "Establishing Cycler Trajectories Between Earth and Mars via "
        "V-Infinity Leveraging' -- EM cycler establishment paper, not yet "
        "anchored."
    ),
    # Other Earth-Mars / Earth-Venus / cislunar cycler papers not yet folded
    # into KNOWN_CORPUS.
    "10.2307/2373181": (
        "Arenstorf 1963, 'Periodic Solutions of the Restricted Three Body "
        "Problem...' -- classical CR3BP figure-8 periodic-orbit paper, not "
        "a cycler-corpus anchor candidate."
    ),
    "10.2514/1.A35091": (
        "Sanchez Net et al. 2022, 'Cycler Orbits and Solar System Pony "
        "Express' -- modern EM cycler constellation-logistics paper, not "
        "yet anchored."
    ),
    "10.2514/1.A35160": (
        "Spreen et al. 2020, 'Design Considerations for an Earth-Mars "
        "Cycler Spacecraft Using the S1L1 Cycler' -- S1L1 spacecraft-design "
        "paper, not yet anchored."
    ),
    "10.2514/3.30134": (
        "Hollister & Menning 1970, 'Periodic Swing-By Orbits between Earth "
        "and Venus' -- foundational Earth-Venus swing-by paper predating "
        "Aldrin, not yet anchored."
    ),
    "10.1016/j.actaastro.2011.03.011": (
        "Lynam & Longuski 2011, 'Laplace-resonant triple-cyclers for "
        "missions to Jupiter' -- a Jovian IEG-adjacent triple-cycler paper "
        "distinct from Hernandez/Jones/Jesick 2017, not yet anchored."
    ),
    "10.2514/6.2022-4345": (
        "Wittal/Miaule/Asher 2022, 'BuzzCraft: Evolution of A Sturdy "
        "Cislunar Cycler Architecture for Permanent Lunar Settlement "
        "Logistics' -- modern cislunar cycler architecture paper, not yet "
        "anchored."
    ),
}


def _literature_dois() -> set[str]:
    """Every DOI cited by a `source: literature` catalogue row.

    Collects BOTH the row's own `first_published.doi` and every DOI in its
    `corroborating_sources[]` list -- a row's corroborating sources are
    additional literature this project already treats as authoritative for
    that row, so they carry the same #577-style false-clear risk as the
    row's primary citation.
    """
    cat = load_catalog()
    dois: set[str] = set()
    for row in cat.entries:
        if row.source != "literature":
            continue
        fp = row.first_published or {}
        doi = fp.get("doi")
        if doi:
            dois.add(doi.strip())
        for corroborating in row.raw.get("corroborating_sources") or ():
            if not isinstance(corroborating, dict):
                continue
            cdoi = corroborating.get("doi")
            if cdoi:
                dois.add(cdoi.strip())
    return dois


def test_every_literature_doi_is_anchored_or_allowlisted() -> None:
    """Every literature DOI is either a KNOWN_CORPUS anchor or allowlisted.

    This is the #578 structural fix: a literature DOI that is neither
    anchored nor allowlisted is exactly the #577 failure mode (a paper the
    catalogue already treats as authoritative that the literature-novelty
    matcher cannot see) and must be closed by adding an anchor or an
    allowlist entry, not silently ignored.
    """
    lit_dois = _literature_dois()
    anchor_dois = {a.doi for a in KNOWN_CORPUS if a.doi}
    uncovered = lit_dois - anchor_dois - set(DOI_ALLOWLIST)
    assert not uncovered, (
        f"{len(uncovered)} literature DOI(s) are neither a KNOWN_CORPUS "
        f"anchor nor allowlisted (add an anchor if this paper should be "
        f"structurally recognised, or a one-line-justified DOI_ALLOWLIST "
        f"entry otherwise): {sorted(uncovered)}"
    )


def test_allowlist_entries_are_still_uncovered() -> None:
    """No stale allowlist entries: every DOI_ALLOWLIST key is a REAL literature DOI.

    Catches a DOI going stale in the allowlist (removed from the catalogue,
    or a typo) rather than silently accumulating dead entries -- every key
    must actually appear among the catalogue's literature DOIs.
    """
    lit_dois = _literature_dois()
    stale = set(DOI_ALLOWLIST) - lit_dois
    assert not stale, (
        f"DOI_ALLOWLIST has {len(stale)} entr(y/ies) that are not (or no "
        f"longer) cited by any `source: literature` catalogue row -- remove "
        f"the stale entry: {sorted(stale)}"
    )


def test_allowlist_entries_are_not_already_anchored() -> None:
    """No redundant allowlist entries: nothing here should already be anchored.

    If a KNOWN_CORPUS anchor is later added for one of these DOIs, its
    allowlist entry becomes redundant (and mildly misleading -- it would read
    as "not yet anchored" when it now is) and should be deleted, not left
    behind.
    """
    anchor_dois = {a.doi for a in KNOWN_CORPUS if a.doi}
    redundant = set(DOI_ALLOWLIST) & anchor_dois
    assert not redundant, (
        f"{len(redundant)} DOI_ALLOWLIST entr(y/ies) are now covered by a "
        f"KNOWN_CORPUS anchor and should be deleted from the allowlist: "
        f"{sorted(redundant)}"
    )
