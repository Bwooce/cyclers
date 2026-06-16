"""Aldrin classic Earth-Mars cycler precursor MGA matcher (#302 / #289 Phase 4).

Runs :func:`cyclerfinder.search.precursor_matcher.find_cycler_precursors` against
the Aldrin row ``aldrin-classic-em-k1-outbound`` over the 2030-2034 launch
window (~2 Earth-Mars synodic periods), enumerates Earth-launched MGA chains
that terminate at Earth (the Aldrin first encounter) at V_inf ~ 6.5 km/s,
closes each via the Phase-1 driver with TOF optimisation, and runs the offline
literature_check anchored on KNOWN_CORPUS.

Output: ``data/precursor_302_aldrin.jsonl`` (one row per surviving candidate,
ranked by ``quality_score``).

Per the discipline preamble in :mod:`precursor_matcher` — every survivor in
this scan is a CANDIDATE for downstream review, NOT a novelty claim.  A
"literature-fresh" tag means the offline corpus did not surface a match; the
human still adjudicates and the V0-V5 gauntlet (when extended to epoch-locked
classes) still governs.

The script prints incremental progress to stdout per
``feedback_incremental_progress_reports`` — never run a black box.
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from pathlib import Path

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.literature_check import KNOWN_CORPUS, SearchResult
from cyclerfinder.search.precursor_matcher import (
    find_cycler_precursors,
    precursor_match_to_jsonl_record,
)

OUTPUT_PATH = Path(__file__).resolve().parents[1] / "data" / "precursor_302_aldrin.jsonl"

# --------------------------------------------------------------------------- #
# Offline literature-check corpus
# --------------------------------------------------------------------------- #
#
# We synthesise SearchResult objects from KNOWN_CORPUS so the matcher's
# literature check runs deterministically against the curated published anchors.
# Each anchor becomes a "hit" whose title carries the anchor name + body codes
# (so the term-overlap matching inside literature_check.py fires for relevant
# queries).  This matches the offline-corpus pattern used in
# tests/search/test_literature_check.py.


def _offline_corpus_hits() -> list[SearchResult]:
    hits: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        bodies_str = " ".join(sorted(anchor.body_set))
        snippet = (
            f"{anchor.citation}. Authors: {', '.join(anchor.authors)}. "
            f"Keywords: {', '.join(anchor.keywords)}. Bodies: {bodies_str}. "
            f"Primary: {anchor.primary}. Subject: cycler trajectory mission design."
        )
        hits.append(
            SearchResult(
                title=anchor.name + " (cycler trajectory)",
                url=anchor.doi or f"https://example.org/{anchor.name.replace(' ', '_')}",
                snippet=snippet,
            )
        )
    return hits


_OFFLINE_HITS = _offline_corpus_hits()


def _tokenise(s: str) -> set[str]:
    return {t for t in "".join(c.lower() if c.isalnum() else " " for c in s).split() if len(t) > 2}


def offline_search(query: str) -> Sequence[SearchResult]:
    """Term-overlap matcher over the KNOWN_CORPUS-derived SearchResult corpus.

    Mirrors the fake_search pattern in tests/search/test_literature_check.py:
    a corpus row is returned iff it shares >=2 lowercased alphabetic tokens
    (length >2) with the query.  Deterministic, no network, runs at all
    epochs (per feedback_keep_task_ledger_synced — offline corpus is fine).
    """
    q_terms = _tokenise(query)
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_HITS:
        text_terms = _tokenise(r.title + " " + r.snippet)
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


# --------------------------------------------------------------------------- #
# Scan
# --------------------------------------------------------------------------- #


def _print_progress(msg: str) -> None:
    """Per feedback_incremental_progress_reports — never run a black box."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def main() -> None:
    _print_progress("Aldrin precursor scan (#302) starting")
    _print_progress("Loading catalogue + DE440 ephemeris")
    cat = load_catalog()
    eph = Ephemeris("astropy")

    cycler_id = "aldrin-classic-em-k1-outbound"
    _print_progress(f"Target cycler: {cycler_id}")
    entry = cat.by_id[cycler_id]
    _print_progress(
        f"  catalogue row: bodies={entry.bodies} period_k={entry.period_k} "
        f"period_years={entry.period_years} V_inf={entry.vinf_kms_at_encounters}"
    )

    # The matcher's launch window is the precursor's first (Earth) departure —
    # one Earth-Mars synodic period is ~779.8 days; 2030-2034 covers roughly
    # two synodics so multiple launch opportunities are captured.
    launch_window = ("2030-01-01T00:00:00", "2034-12-31T00:00:00")
    _print_progress(f"Launch window: {launch_window} (~2 Earth-Mars synodic periods)")

    t0 = time.time()
    matches = find_cycler_precursors(
        cycler_id=cycler_id,
        catalogue=cat,
        ephemeris=eph,
        launch_window=launch_window,
        # Aldrin's seed V_inf at Earth is 6.5 km/s; the precursor must hit
        # that within a tolerance.
        vinf_terminal_tol_kms=0.8,
        # Discovery-probe: admit candidates whose proposal V_inf was in
        # tolerance but whose closure V_inf drifted.  We WANT to see the
        # residual distribution in the JSONL.
        vinf_terminal_post_closure_tol_kms=50.0,
        # 2-3 intermediate flybys reach the published Aldrin entry V_inf.
        max_legs=3,
        # Intermediate bodies for inner-planet precursors are Venus + Earth
        # (the standard precursor archetype; e.g. VEGA-style).
        intermediate_bodies=("V", "E"),
        # V_inf grid spans the published Aldrin Earth entry value (6.5).
        vinf_grid_kms=(4.0, 5.0, 6.0, 7.0, 8.0),
        tof_box_days_per_leg=(80.0, 500.0),
        epoch_step_days=60.0,
        tof_optimise=True,
        # Phase-4 closure gates: wide so the JSONL captures residuals across
        # the search surface, not just publication-grade hits.  The discovery
        # probe's POINT is to see WHERE the residual lives.  A tight gate
        # would silently filter out structural information.  The doc reports
        # the distribution + best survivor explicitly.
        closure_tol_kms=100.0,
        flyby_continuity_tol_kms=100.0,
        independent_cross_check=False,  # speed; final closure runs xcheck on
        independent_tol_kms=100.0,
        multi_shell=True,
        pump_envelope_factor=1.0,
        a_range_au=(0.3, 2.5),  # inner-planet precursor; Aldrin reaches Mars only
        # Cap candidate validation so the scan stays under ~15 min.
        max_candidates_to_validate=400,
        literature_check_search=offline_search,
        progress_hook=_print_progress,
    )
    elapsed = time.time() - t0
    _print_progress(f"Scan complete in {elapsed:.1f}s — {len(matches)} survivors")

    # Stratify by literature-check verdict.
    n_fresh = sum(1 for m in matches if m.is_literature_fresh())
    n_published = sum(1 for m in matches if m.literature_check.status == "published")
    n_inconclusive = sum(1 for m in matches if m.literature_check.status == "inconclusive")
    _print_progress(
        f"Literature check breakdown: "
        f"{n_fresh} not-found / {n_published} published / "
        f"{n_inconclusive} inconclusive"
    )
    if matches:
        best = matches[0]
        _print_progress(
            f"Best (quality_score={best.quality_score():.3f}): "
            f"sequence={best.candidate.sequence} "
            f"launch={best.candidate.launch_epoch_utc} "
            f"closure={best.closure.closure_residual_kms:.3f} km/s "
            f"V_inf_match_residual={best.vinf_match_residual_kms:.3f} km/s "
            f"literature={best.literature_check.status}"
        )

    # Emit JSONL.
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as f:
        for m in matches:
            f.write(json.dumps(precursor_match_to_jsonl_record(m)) + "\n")
    _print_progress(f"Wrote {OUTPUT_PATH.relative_to(OUTPUT_PATH.parents[1])}")


if __name__ == "__main__":
    main()
