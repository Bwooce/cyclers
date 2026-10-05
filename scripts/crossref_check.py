"""#964 -- confirm a citation's DOI against Crossref, requiring a journal or volume match.

Crossref's bibliographic search often ranks the AIAA *conference* record
(10.2514/6.YYYY-NNNN, a ``proceedings-article``) above the journal version of
the same paper, and a bare DOI lookup only proves that the DOI exists. This
helper prints ``CONFIRMED`` only when the record's container title or volume
matches what the citation says (and, when given, the title overlaps and the
year agrees). Otherwise it prints ``MISMATCH`` with the fields that differ, or
``NOT_FOUND``.

Usage::

    # check a DOI against the citation's journal and volume
    python scripts/crossref_check.py doi 10.2514/3.29664 \\
        --container "Journal of Spacecraft and Rockets" --volume 6 --year 1969

    # search by citation text; reports every candidate and picks the first
    # record whose container/volume match (journal before conference)
    python scripts/crossref_check.py search "Gillespie Ross Venus swingby mission mode" \\
        --container "Journal of Spacecraft" --volume 4

Read-only; stdlib only (urllib). The pure functions ``assess`` and
``pick_match`` are unit-tested offline against recorded responses in
``tests/scripts/fixtures/``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

API = "https://api.crossref.org/works"
USER_AGENT = "cyclers-research crossref_check.py (https://github.com/Bwooce/cyclers)"

Record = dict[str, Any]
Fetch = Callable[[str], dict[str, Any]]


@dataclass(frozen=True)
class Expect:
    """What the citation claims; ``None`` means "not claimed, not checked"."""

    container: str | None = None
    volume: str | None = None
    title: str | None = None
    year: int | None = None


@dataclass
class Verdict:
    status: str  # CONFIRMED | MISMATCH | NOT_FOUND | UNCHECKED
    doi: str | None
    record: Record | None = None
    reasons: list[str] = field(default_factory=list)

    def line(self) -> str:
        rec = self.record or {}
        container = "; ".join(rec.get("container-title") or []) or "-"
        title = (rec.get("title") or ["-"])[0]
        vol = rec.get("volume") or "-"
        page = rec.get("page") or "-"
        year = _year(rec)
        why = f" [{'; '.join(self.reasons)}]" if self.reasons else ""
        return (
            f"{self.status} {self.doi} | {title[:90]} | {container} {vol}:{page} "
            f"({year}) | {rec.get('type', '-')}{why}"
        )


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def _year(rec: Record) -> int | None:
    for key in ("published-print", "issued", "published-online"):
        parts = (rec.get(key) or {}).get("date-parts") or []
        if parts and parts[0] and parts[0][0]:
            return int(parts[0][0])
    return None


def _container_matches(rec: Record, want: str) -> bool:
    w = _norm(want)
    return any(w in _norm(c) or _norm(c) in w for c in rec.get("container-title") or [] if c)


def _title_overlap(rec: Record, want: str) -> float:
    got = set(_norm(" ".join(rec.get("title") or [])).split())
    words = {t for t in _norm(want).split() if len(t) > 3}
    return len(words & got) / len(words) if words else 0.0


def assess(rec: Record | None, expect: Expect) -> Verdict:
    """Judge one Crossref record against the citation's claims.

    CONFIRMED needs at least one of container/volume claimed AND every claimed
    container/volume to match; a claimed title must overlap >= 60 percent of its
    long words and a claimed year must agree within one year (print vs online).
    """
    if rec is None:
        return Verdict("NOT_FOUND", None)
    doi = rec.get("DOI")
    reasons: list[str] = []
    if expect.container is None and expect.volume is None:
        return Verdict("UNCHECKED", doi, rec, ["no container or volume claimed"])
    if expect.container is not None and not _container_matches(rec, expect.container):
        reasons.append(f"container {rec.get('container-title')} != {expect.container!r}")
    if expect.volume is not None and str(rec.get("volume") or "") != str(expect.volume):
        reasons.append(f"volume {rec.get('volume')!r} != {expect.volume!r}")
    if expect.title is not None and _title_overlap(rec, expect.title) < 0.6:
        reasons.append("title words do not overlap")
    year = _year(rec)
    if expect.year is not None and (year is None or abs(year - expect.year) > 1):
        reasons.append(f"year {year} != {expect.year}")
    if rec.get("type") == "proceedings-article" and expect.container is not None:
        reasons.append("record is a conference paper")
    return Verdict("MISMATCH" if reasons else "CONFIRMED", doi, rec, reasons)


def pick_match(items: Sequence[Record], expect: Expect) -> tuple[Verdict, list[Verdict]]:
    """Assess every search candidate; return the first CONFIRMED one (else the best)."""
    verdicts = [assess(it, expect) for it in items]
    for v in verdicts:
        if v.status == "CONFIRMED":
            return v, verdicts
    if not verdicts:
        return Verdict("NOT_FOUND", None), verdicts
    return min(verdicts, key=lambda v: len(v.reasons)), verdicts


def _http_get(url: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data: dict[str, Any] = json.load(resp)
    return data


def fetch_doi(doi: str, fetch: Fetch = _http_get) -> Record | None:
    try:
        msg: Record = fetch(f"{API}/{urllib.parse.quote(doi)}")["message"]
    except Exception:  # 404 or network error: reported as NOT_FOUND
        return None
    return msg


def search(query: str, rows: int = 5, fetch: Fetch = _http_get) -> list[Record]:
    params = urllib.parse.urlencode({"query.bibliographic": query, "rows": rows})
    items: list[Record] = fetch(f"{API}?{params}")["message"]["items"]
    return items


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("mode", choices=("doi", "search"))
    ap.add_argument("value", help="a DOI (mode doi) or citation text (mode search)")
    ap.add_argument("--container", help="journal or proceedings name the citation gives")
    ap.add_argument("--volume", help="volume the citation gives")
    ap.add_argument("--title", help="title the citation gives (doi mode)")
    ap.add_argument("--year", type=int, help="year the citation gives")
    ap.add_argument("--rows", type=int, default=5)
    args = ap.parse_args(argv)
    expect = Expect(args.container, args.volume, args.title, args.year)
    if args.mode == "doi":
        verdict = assess(fetch_doi(args.value), expect)
        print(verdict.line())
    else:
        verdict, all_v = pick_match(search(args.value, args.rows), expect)
        for v in all_v:
            print("  candidate:", v.line())
        print(verdict.line())
    return 0 if verdict.status == "CONFIRMED" else 1


if __name__ == "__main__":
    sys.exit(main())
