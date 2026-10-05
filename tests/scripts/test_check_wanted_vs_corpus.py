"""#964 -- scripts/check_wanted_vs_corpus.py flags wanted-list rows the corpus already holds.

Offline: a synthetic four-file corpus; no PDFs, no network. The cases are the
three misses of the #960 lists (Bruno filed as "brjuno", Henon 2001, and Rall 1969
filed as "hollister-rall-1970"), plus a genuinely missing paper and a title-only hit.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "check_wanted_vs_corpus.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("check_wanted_vs_corpus", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses resolve their module via sys.modules
    spec.loader.exec_module(mod)
    return mod


cw = _load()

HEADER = (
    "| Rank | Citation | DOI / identifier | Free PDF source | Unlocks |\n|---|---|---|---|---|\n"
)
CITATIONS = [
    (
        'Bruno, A. D. (1978a, b), "Researches on the restricted three-body problem II, III", '
        "Celest. Mech. 18:9-50",
        "10.1007/BF01233089 CONFIRMED",
    ),
    (
        'Henon, M. (2001), "Generating Families of the Restricted Three-Body Problem. II. '
        'Quantitative Study of Bifurcations", LNP 65',
        "not checked",
    ),
    ('Turner, A. (2007), "Low Road to Mars: The Venus-Mars Cycler", AAS 07-175', "none"),
    (
        'Rall, C. S. (1969), "Freefall Periodic Orbits Connecting Earth and Mars", PhD thesis, MIT',
        "none",
    ),
    (
        "**[Content held in another form; acquire only for attribution.]** "
        'Russell, R. P. & Ocampo, C. A. (2005), "Global Search for Idealized Free-Return '
        'Earth-Mars Cyclers", JGCD 28(2):194-208',
        "10.2514/1.8696 CONFIRMED",
    ),
]
WANTED = HEADER + "".join(
    f"| {n} | {cit} | {doi} | n/a | x |\n" for n, (cit, doi) in enumerate(CITATIONS, 1)
)

FILENAMES = [
    "brjuno-1978-researches-restricted-three-body-problem-II-celest-mech-18-9-doi-10.1007-BF01233089.pdf",
    "henon-2001-generating-families-restricted-three-body-problem-II-lnp-m65.pdf",
    "vaquero-2013-purdue-phd.txt",
    "hollister-rall-1970-periodic-orbits-NASA-CR.pdf",
]
TEXTS = {
    FILENAMES[
        0
    ]: "Celestial Mechanics 18 (1978) 9-50. Researches on the restricted three-body problem",
    FILENAMES[1]: "Lecture Notes in Physics m65. Michel Hénon. Generating Families",
    FILENAMES[2]: (
        "Chapter 3 ... Russell and Ocampo, Global Search for Idealized "
        "Free-Return Earth-Mars Cyclers, submitted"
    ),
    FILENAMES[3]: "FREE-FALL PERIODIC ORBITS CONNECTING EARTH AND MARS  CHARLES SHERMAN RALL",
}


@pytest.fixture
def hits() -> list:  # type: ignore[type-arg]
    rows = cw.parse_rows(WANTED)
    return [cw.check_row(r, cw.Corpus.build(TEXTS, FILENAMES, "")) for r in rows]


def test_parse_rows_reads_rank_surname_year_title_and_doi() -> None:
    rows = cw.parse_rows(WANTED)
    assert [r.rank for r in rows] == [1, 2, 3, 4, 5]
    assert (rows[0].surname, rows[0].year, rows[0].dois) == (
        "Bruno",
        "1978",
        ["10.1007/BF01233089"],
    )
    assert rows[4].surname == "Russell"  # the bold editorial note is stripped first
    assert rows[2].title == "Low Road to Mars: The Venus-Mars Cycler"


def test_transliteration_alias_catches_brjuno(hits: list) -> None:  # type: ignore[type-arg]
    h = hits[0]
    assert h.needs_look
    assert FILENAMES[0] in h.doi  # DOI in the filename
    assert FILENAMES[0] in h.name_year  # Bruno == Brjuno


def test_name_year_catches_henon_2001(hits: list) -> None:  # type: ignore[type-arg]
    assert hits[1].name_year == [FILENAMES[1]]


def test_missing_paper_is_clean(hits: list) -> None:  # type: ignore[type-arg]
    assert not hits[2].needs_look and not hits[2].title


def test_title_only_hit_does_not_fail(hits: list) -> None:  # type: ignore[type-arg]
    h = hits[4]
    assert h.title == [FILENAMES[2]] and not h.needs_look


def test_advisor_first_filename_and_report_year_catch_rall(hits: list) -> None:  # type: ignore[type-arg]
    # Rall's 1969 Sc.D. thesis is held as "hollister-rall-1970-..."; the first version missed it.
    h = hits[3]
    assert h.name_year == [FILENAMES[3]]
    assert h.title == [FILENAMES[3]]  # "Freefall" matches "FREE-FALL"


def test_accent_folding_and_aliases() -> None:
    assert cw.fold("Hénon Ollé") == "henon olle"
    assert "brjuno" in cw.aliases("Bruno") and "bruno" in cw.aliases("Brjuno")


def test_index_identity_note_catches_misleading_filename() -> None:
    row = cw.parse_rows(WANTED)[3]  # Rall 1969
    index = (
        "| smith-1970-report.pdf | note.md | [identity: Rall, C. S. (1969), Sc.D. thesis, MIT] |"
        " mined | text-layer |"
    )
    h = cw.check_row(row, cw.Corpus.build({}, ["smith-1970-report.pdf"], index))
    assert h.name_year == ["smith-1970-report.pdf"] and h.needs_look


def test_index_doi_hit() -> None:
    row = cw.parse_rows(WANTED)[4]
    h = cw.check_row(row, cw.Corpus.build({}, [], "| x.pdf | ... DOI 10.2514/1.8696 ... |"))
    assert h.doi == ["CORPUS_INDEX.md"] and h.needs_look
