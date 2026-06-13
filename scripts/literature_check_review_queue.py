"""Backfill the review-queue ``literature_check`` field (#261, thin driver).

A SILVER candidate from the #253 discovery daemon is only deduplicated against
our sourced catalogue + the negative registry -- a *subset* of the published
cycler literature. This driver closes that gap: it walks
``data/review_queue.jsonl``, builds each SILVER entry's structural SIGNATURE, runs
the :mod:`cyclerfinder.search.literature_check` checker against the published
record, and rewrites the entry with its ``literature_check`` block populated.

NO candidate may be claimed novel until this has run on it:
:func:`cyclerfinder.search.literature_check.is_novelty_claimable` treats an
unpopulated or ``published`` block as NOT novelty-claimable.

SEARCH BACKEND. WebSearch is a harness tool, not a Python import, so this driver
takes a pluggable ``SearchFn``. Two backends ship:

* ``--backend offline`` (default): matches each signature against the curated
  :data:`~cyclerfinder.search.literature_check.KNOWN_CORPUS` only. This is a
  conservative, deterministic FILTER -- it flags candidates whose tour falls
  inside a known published family, and reports ``inconclusive`` (NOT
  ``not-found``) for everything else, because an offline run has not actually
  searched the web and must never emit a false "not-found" clean bill.
* ``--backend results-json PATH``: replays pre-collected WebSearch results
  (JSON: ``{query: [{"title","url","snippet"}, ...]}``) -- the path used when a
  human / agent has run the live WebSearch tool and saved the hits. This is how
  a trustworthy ``not-found`` is produced (a real search returned nothing).

This separation keeps the daemon quota-proof (offline filter inline) while the
authoritative live search is run deliberately and recorded.

NO catalogue writeback: only ``data/review_queue.jsonl`` is rewritten in place.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from cyclerfinder.data.review_queue import (
    DEFAULT_REVIEW_QUEUE_PATH,
    ReviewQueueEntry,
    append_review_entry,
    load_review_queue,
)
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    SearchFn,
    SearchResult,
    check_literature,
    signature_from_review_entry,
)


def offline_corpus_search(query: str) -> Sequence[SearchResult]:
    """A deterministic offline backend: synthesise hits from the known corpus.

    For each curated anchor whose author/keyword appears in the query, emit a
    result whose title carries the anchor's family name + a tour body + the word
    "cycler" so the structural matcher can score it. This is NOT a web search --
    it only re-finds families ALREADY in the curated corpus, so a candidate
    inside a known family is flagged ``published`` while everything else falls
    through to ``inconclusive`` (never a false offline ``not-found``).
    """
    q = query.lower()
    out: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        hit = any(a.lower() in q for a in anchor.authors) or any(
            kw.lower() in q for kw in anchor.keywords
        )
        # Also fire when the query names the anchor's bodies + "cycler".
        bodies_named = sum(1 for b in anchor.body_set if b.lower() in q)
        if hit or (bodies_named >= 2 and "cycler" in q):
            bodies = " ".join(sorted(anchor.body_set))
            out.append(
                SearchResult(
                    title=f"{anchor.name} ({bodies} cycler)",
                    url=(f"https://doi.org/{anchor.doi}" if anchor.doi else anchor.citation),
                    snippet=f"{anchor.citation}. {' '.join(anchor.keywords)}. "
                    f"Authors: {', '.join(anchor.authors)}.",
                )
            )
    return out


def results_json_backend(path: Path) -> SearchFn:
    """Replay pre-collected live-WebSearch results from a JSON file.

    JSON shape: ``{query_string: [{"title","url","snippet"}, ...]}``. An exact
    query match returns its recorded hits; an unrecorded query returns ``[]``
    (so the checker sees an empty -- but real -- search, which contributes to a
    trustworthy not-found only if every query was actually recorded).
    """
    data = json.loads(path.read_text(encoding="utf-8"))

    def _search(query: str) -> Sequence[SearchResult]:
        rows = data.get(query, [])
        return [
            SearchResult(
                title=r.get("title", ""),
                url=r.get("url", ""),
                snippet=r.get("snippet", ""),
            )
            for r in rows
        ]

    return _search


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--queue",
        default=str(DEFAULT_REVIEW_QUEUE_PATH),
        help="path to data/review_queue.jsonl",
    )
    ap.add_argument(
        "--backend",
        choices=("offline", "results-json"),
        default="offline",
        help="search backend (offline corpus filter, or replay live WebSearch json)",
    )
    ap.add_argument(
        "--results-json",
        default=None,
        help="path to recorded WebSearch results (for --backend results-json)",
    )
    ap.add_argument(
        "--only-unchecked",
        action="store_true",
        help="skip entries whose literature_check is already populated",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print verdicts; do not rewrite the queue",
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    if args.backend == "results-json":
        if not args.results_json:
            raise SystemExit("--backend results-json requires --results-json PATH")
        search: SearchFn = results_json_backend(Path(args.results_json))
    else:
        search = offline_corpus_search

    queue_path = Path(args.queue)
    entries = list(load_review_queue(queue_path))
    if not entries:
        print(f"[lit-check] no entries in {queue_path}")
        return

    updated: list[ReviewQueueEntry] = []
    n_pub = n_nf = n_inc = n_skip = 0
    for entry in entries:
        if args.only_unchecked and entry.literature_check:
            updated.append(entry)
            n_skip += 1
            continue
        if entry.verdict_tier != "silver":
            updated.append(entry)
            n_skip += 1
            continue
        sig = signature_from_review_entry(entry)
        result = check_literature(sig, search=search)
        block = result.to_review_block()
        new_entry = ReviewQueueEntry(**{**entry.__dict__, "literature_check": block})
        updated.append(new_entry)
        if result.status == "published":
            n_pub += 1
        elif result.status == "not-found":
            n_nf += 1
        else:
            n_inc += 1
        print(
            f"[lit-check] {entry.candidate_id}: {result.status} "
            f"conf={result.confidence} citation={result.citation!r}"
        )

    print(
        f"[lit-check] published={n_pub} not-found={n_nf} inconclusive={n_inc} "
        f"skipped={n_skip} total={len(entries)}"
    )
    if args.dry_run:
        print("[lit-check] dry-run: queue NOT rewritten")
        return

    # Rewrite the queue in place (validation re-runs per entry).
    tmp = queue_path.with_suffix(".jsonl.tmp")
    if tmp.exists():
        tmp.unlink()
    for e in updated:
        append_review_entry(tmp, e)
    tmp.replace(queue_path)
    print(f"[lit-check] rewrote {queue_path} ({len(updated)} entries)")


if __name__ == "__main__":
    main()
