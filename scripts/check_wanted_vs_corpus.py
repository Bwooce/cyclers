"""#964 -- check every row of a "wanted papers" list against the whole held corpus.

Before an acquisition list goes to the owner, no row may be a paper the corpus
already holds (docs/notes/corpus-document-policy.md, "Wanted lists"). The
first #960 list carried two held papers: Bruno 1978a,b, filed under the
transliteration "brjuno-1978", and Henon 2001. A plain filename grep missed
them, and the first version of this script missed Rall 1969 (see check 2).

For each numbered table row ``| N | citation | DOI ... |`` of the list this
reports:

1. DOI hits: a DOI from the row found in a corpus filename, in the first
   pages of a corpus PDF (or its .txt companion), or in CORPUS_INDEX.md.
2. Name+year hits: the first author's surname (with transliteration aliases
   such as Bruno/Brjuno, Henon/Hénon, Olle/Ollé) among the author tokens of a
   corpus filename (the tokens before its year), and a filename year within
   one year of the row's. Any author position counts: Rall's 1969 thesis was
   filed as "hollister-rall-1970-..." (advisor first, report year). The same
   test runs on every CORPUS_INDEX "[identity: Author (Year), Title]" note,
   which records the title page when a filename does not show it.
3. Title hits: the first 32 letters and digits of the row's quoted title,
   spaces removed ("Freefall" = "Free-Fall"), found in a corpus text. These
   are mostly citations inside other papers, so they are listed but do not
   fail the run.

Name+year and DOI hits need a human look. A held paper is removed from the
list; a paper held in another form (report, preprint, thesis) is marked
"acquire only for attribution".

Usage::

    python scripts/check_wanted_vs_corpus.py docs/notes/2026-10-05-960-wanted-papers.md
    python scripts/check_wanted_vs_corpus.py LIST --papers ../cyclers_pdf/papers --pages 3

Read-only. Extracted text is cached under ``--cache`` (default: the system
temp directory). Exit status 1 when any row has a DOI or name+year hit.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import unicodedata
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PAPERS = REPO_ROOT.parent / "cyclers_pdf" / "papers"
DEFAULT_INDEX = REPO_ROOT / "docs" / "notes" / "CORPUS_INDEX.md"

#: Groups of surname spellings that name the same author (compared after accent folding).
ALIAS_GROUPS: tuple[frozenset[str], ...] = (
    frozenset({"bruno", "brjuno", "brunno"}),
    frozenset({"henon"}),  # Hénon folds to henon
    frozenset({"olle"}),  # Ollé folds to olle
    frozenset({"stromgren", "stroemgren"}),
    frozenset({"poincare"}),
    frozenset({"szebehely"}),
    frozenset({"lyapunov", "liapunov", "ljapunov"}),
    frozenset({"chebyshev", "tchebychev", "chebyshov"}),
)

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s|;,()\]]+(?:\([^\s|;,)]*\)[^\s|;,()\]]*)*", re.IGNORECASE)
YEAR_RE = re.compile(r"(?<!\d)(1[89]\d\d|20\d\d)(?!\d)")  # also "1978a"
ROW_RE = re.compile(r"^\|\s*(\d+)\s*\|(.*)$")


def fold(s: str) -> str:
    """Lower-case and strip accents (Hénon -> henon, Ollé -> olle)."""
    nfkd = unicodedata.normalize("NFKD", s)
    return "".join(c for c in nfkd if not unicodedata.combining(c)).lower()


def squash(s: str) -> str:
    """Letters and digits only, so hyphenation and OCR letter-spacing do not matter."""
    return re.sub(r"[^a-z0-9]+", "", fold(s))


def filename_authors_year(fn: str) -> tuple[list[str], str | None]:
    """Author tokens (those before the first year token) and that year, from a corpus filename."""
    tokens = re.split(r"[^a-z0-9]+", fold(fn))
    for i, t in enumerate(tokens):
        m = re.fullmatch(r"(1[89]\d\d|20\d\d)[a-z]?", t)
        if m:
            return tokens[:i], m.group(1)
    return tokens, None


def aliases(surname: str) -> set[str]:
    f = fold(surname)
    for group in ALIAS_GROUPS:
        if f in group:
            return set(group)
    return {f}


@dataclass
class Row:
    rank: int
    citation: str
    doi_field: str
    raw: str
    dois: list[str] = field(default_factory=list)
    surname: str | None = None
    year: str | None = None
    title: str | None = None


def _clean(cell: str) -> str:
    cell = re.sub(r"\*\*\[[^\]]*\]\*\*", " ", cell)  # bold editorial notes like **[...]**
    return cell.replace("**", "").strip()


def parse_rows(markdown: str) -> list[Row]:
    rows: list[Row] = []
    for line in markdown.splitlines():
        m = ROW_RE.match(line)
        if not m:
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        citation = _clean(cells[1])
        doi_field = cells[2]
        row = Row(int(m.group(1)), citation, doi_field, line)
        row.dois = sorted({d.rstrip(".") for d in DOI_RE.findall(citation + " " + doi_field)})
        sm = re.match(r"\s*([A-Z][A-Za-zÀ-ÿ'\-]+)", citation)
        row.surname = sm.group(1) if sm else None
        ym = YEAR_RE.search(citation)
        row.year = ym.group(1) if ym else None
        tm = re.search(r'"([^"]{8,})"', citation)
        row.title = tm.group(1) if tm else None
        rows.append(row)
    return rows


@dataclass
class Hits:
    row: Row
    doi: list[str] = field(default_factory=list)
    name_year: list[str] = field(default_factory=list)
    title: list[str] = field(default_factory=list)

    @property
    def needs_look(self) -> bool:
        return bool(self.doi or self.name_year)


@dataclass(frozen=True)
class Corpus:
    """Corpus texts pre-folded once; folding per row made a 60-row list take minutes."""

    folded: dict[str, str]
    squashed: dict[str, str]
    filenames: list[str]
    folded_index: str
    identities: list[tuple[str, str]] = field(default_factory=list)

    @classmethod
    def build(cls, texts: Mapping[str, str], filenames: Iterable[str], index_text: str) -> Corpus:
        return cls(
            {fn: fold(t) for fn, t in texts.items()},
            {fn: squash(t) for fn, t in texts.items()},
            list(filenames),
            fold(index_text),
            index_identities(index_text),
        )


IDENTITY_RE = re.compile(r"\[identity:\s*([^\]]+)\]")


def index_identities(index_text: str) -> list[tuple[str, str]]:
    """(filename, folded identity) for each CORPUS_INDEX row with an [identity: ...] note.

    The note records the title-page author, year and title when the filename does not show them
    (a thesis filed under its supervisor, a report year, a transliteration).
    """
    out = []
    for line in index_text.splitlines():
        m = IDENTITY_RE.search(line)
        if not m or "|" not in line:
            continue
        fn = line.split("|")[1 if line.lstrip().startswith("|") else 0]
        fn = fn.strip().lstrip("- ").strip("*` ")
        out.append((fn, fold(m.group(1))))
    return out


def check_row(row: Row, corpus: Corpus) -> Hits:
    """Match one row against corpus texts, filenames and the index."""
    hits = Hits(row)
    fnames = corpus.filenames
    folded_index = corpus.folded_index
    for doi in row.dois:
        d = fold(doi)
        d_file = d.replace("/", "-")
        for fn in fnames:
            if d_file in fold(fn) and fn not in hits.doi:
                hits.doi.append(fn)
        for fn, txt in corpus.folded.items():
            if d in txt and fn not in hits.doi:
                hits.doi.append(fn)
        if d in folded_index and "CORPUS_INDEX.md" not in hits.doi:
            hits.doi.append("CORPUS_INDEX.md")
    if row.surname and row.year:
        names = aliases(row.surname)
        years = {str(int(row.year) + d) for d in (-1, 0, 1)}
        for fn in fnames:
            authors, year = filename_authors_year(fn)
            if names & set(authors) and year in years:
                hits.name_year.append(fn)
        for fn, ident in corpus.identities:
            words = set(re.split(r"[^a-z0-9]+", ident))
            if names & words and years & words and fn not in hits.name_year:
                hits.name_year.append(fn)
    if row.title:
        key = squash(row.title)[:32]
        if len(key) >= 15:
            for fn, txt in corpus.squashed.items():
                if key in txt:
                    hits.title.append(fn)
    return hits


def corpus_texts(papers: Path, cache: Path, pages: int) -> dict[str, str]:
    """First ``pages`` pages of every PDF (via pdftotext, cached), plus .txt companions."""
    cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, str] = {}
    for pdf in sorted(papers.glob("*.pdf")):
        cached = cache / (pdf.name + f".p{pages}.txt")
        if not cached.exists() or cached.stat().st_mtime < pdf.stat().st_mtime:
            res = subprocess.run(
                ["pdftotext", "-l", str(pages), str(pdf), "-"],
                capture_output=True,
                text=True,
                check=False,
            )  # stderr (thousands of "Badly formatted number" lines on AIAA PDFs) is discarded
            cached.write_text(res.stdout)
        out[pdf.name] = cached.read_text()
    for txt in sorted(papers.glob("*.txt")):
        out[txt.name] = txt.read_text(errors="replace")
    return out


def report(all_hits: Sequence[Hits]) -> str:
    lines = []
    for h in all_hits:
        r = h.row
        tag = "LOOK" if h.needs_look else ("title" if h.title else "ok  ")
        head = f"{tag} row {r.rank}: {r.surname or '?'} {r.year or '?'} | {r.citation[:70]}"
        lines.append(head)
        for kind, names in (("DOI", h.doi), ("name+year", h.name_year), ("title", h.title)):
            if names:
                shown = ", ".join(names[:4]) + (f" (+{len(names) - 4})" if len(names) > 4 else "")
                lines.append(f"      {kind}: {shown}")
    n_look = sum(h.needs_look for h in all_hits)
    lines.append(f"{len(all_hits)} rows; {n_look} need a human look (DOI or name+year hit)")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "wanted", type=Path, help="markdown wanted list with '| N | citation | DOI |' rows"
    )
    ap.add_argument("--papers", type=Path, default=DEFAULT_PAPERS)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--pages", type=int, default=3)
    ap.add_argument(
        "--cache", type=Path, default=Path(tempfile.gettempdir()) / "cyclers-corpus-txt"
    )
    args = ap.parse_args(argv)
    rows = parse_rows(args.wanted.read_text())
    texts = corpus_texts(args.papers, args.cache, args.pages)
    filenames = sorted(p.name for p in args.papers.iterdir())
    index_text = args.index.read_text() if args.index.exists() else ""
    corpus = Corpus.build(texts, filenames, index_text)
    all_hits = [check_row(r, corpus) for r in rows]
    print(report(all_hits))
    return 1 if any(h.needs_look for h in all_hits) else 0


if __name__ == "__main__":
    sys.exit(main())
