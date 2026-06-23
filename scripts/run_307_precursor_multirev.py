"""#307 Task 4 — re-run the Aldrin / S1L1 precursor scans with multi-rev Lambert.

The #302 single-rev scans (``data/precursor_302_{aldrin,s1l1}.jsonl``) found
NO publication-grade precursor into either classic Earth-Mars cycler: the
binding wall is the *flyby-continuity* discontinuity (~4 km/s ballistic
mismatch at the gravity-assist node, far above the 0.10 km/s ballistic gate),
not the terminal V_inf match (which already lands ~0.06-0.21 km/s).

Task 4 wires the multi-rev Lambert branch-selection knob (``max_revs``, #307
Task 1) through the matcher and asks the empirical question: does enumerating
the multi-rev Lambert branches per leg open a low-continuity (near-ballistic)
geometry that single-rev missed?  A multi-rev branch can match a leg's
endpoints with a different heliocentric arc whose departure/arrival velocities
align better with the flyby-continuity constraint.

This script re-runs the *identical* search box as ``run_302_*`` (same launch
window, V_inf grid, TOF box, gates) with ``max_revs=2`` so the only changed
variable is the Lambert branch set.  Output:
``data/precursor_307_{target}_multirev.jsonl``.

Per the discipline preamble in :mod:`precursor_matcher` every survivor is a
CANDIDATE for review, never a novelty claim.  The expected outcome is an
honest negative (precursor-MGA space into the published E-M cyclers is itself
published — #302 Phase-4 note); this run characterises WHETHER multi-rev
changes the wall, and the doc reports the residual distribution either way.

Usage::

    uv run python scripts/run_307_precursor_multirev.py aldrin
    uv run python scripts/run_307_precursor_multirev.py s1l1
"""

from __future__ import annotations

import json
import sys
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

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Per-target config — mirrors run_302_aldrin_precursor.py / run_302_s1l1_precursor.py
# exactly, so the ONLY changed variable vs the #302 baseline is max_revs.
_TARGETS: dict[str, dict[str, object]] = {
    "aldrin": {
        "cycler_id": "aldrin-classic-em-k1-outbound",
        "vinf_grid_kms": (4.0, 5.0, 6.0, 7.0, 8.0),
        "baseline": "precursor_302_aldrin.jsonl",
    },
    "s1l1": {
        "cycler_id": "s1l1-2syn-em-cpom",
        "vinf_grid_kms": (3.0, 4.0, 5.0, 6.0, 7.0),
        "baseline": "precursor_302_s1l1.jsonl",
    },
}

# Multi-rev branch budget. max_revs=2 enumerates 0,1,2-rev Lambert branches per
# leg; the matcher picks the closure-residual-minimising combination. 2 is the
# standard inner-solar-system precursor ceiling (a 3rd rev rarely fits the TOF
# box and multiplies the per-leg Lambert solve count).
_MAX_REVS = 2


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
    q_terms = _tokenise(query)
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_HITS:
        text_terms = _tokenise(r.title + " " + r.snippet)
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def _print_progress(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _baseline_summary(path: Path) -> str:
    if not path.exists():
        return "baseline missing"
    rows = [json.loads(line) for line in path.open()]
    if not rows:
        return "baseline empty"
    cr = [r["closure"]["closure_residual_kms"] for r in rows if r.get("closure")]
    fc = [r["closure"]["flyby_continuity_max_dv_kms"] for r in rows if r.get("closure")]
    fresh = sum(1 for r in rows if r.get("is_literature_fresh"))
    return (
        f"{len(rows)} rows, {fresh} fresh, min_closure={min(cr):.4f} km/s, "
        f"min_flyby_cont={min(fc):.4f} km/s"
    )


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in _TARGETS:
        raise SystemExit(f"usage: run_307_precursor_multirev.py {{{'|'.join(_TARGETS)}}}")
    target = sys.argv[1]
    cfg = _TARGETS[target]
    cycler_id = str(cfg["cycler_id"])
    output_path = _DATA_DIR / f"precursor_307_{target}_multirev.jsonl"
    baseline_path = _DATA_DIR / str(cfg["baseline"])

    _print_progress(f"#307 Task 4 multi-rev precursor scan: {target} (max_revs={_MAX_REVS})")
    _print_progress(f"  #302 single-rev baseline: {_baseline_summary(baseline_path)}")
    _print_progress("Loading catalogue + DE440 ephemeris")
    cat = load_catalog()
    eph = Ephemeris("astropy")

    entry = cat.by_id[cycler_id]
    _print_progress(
        f"  target row: bodies={entry.bodies} period_k={entry.period_k} "
        f"V_inf={entry.vinf_kms_at_encounters}"
    )

    launch_window = ("2030-01-01T00:00:00", "2034-12-31T00:00:00")
    _print_progress(f"Launch window: {launch_window} (~2 Earth-Mars synodic periods)")

    t0 = time.time()
    matches = find_cycler_precursors(
        cycler_id=cycler_id,
        catalogue=cat,
        ephemeris=eph,
        launch_window=launch_window,
        vinf_terminal_tol_kms=0.8,
        vinf_terminal_post_closure_tol_kms=50.0,
        max_legs=3,
        intermediate_bodies=("V", "E"),
        vinf_grid_kms=tuple(cfg["vinf_grid_kms"]),  # type: ignore[arg-type]
        tof_box_days_per_leg=(80.0, 500.0),
        epoch_step_days=60.0,
        tof_optimise=True,
        closure_tol_kms=100.0,
        flyby_continuity_tol_kms=100.0,
        independent_cross_check=False,
        independent_tol_kms=100.0,
        multi_shell=True,
        pump_envelope_factor=1.0,
        a_range_au=(0.3, 2.5),
        max_candidates_to_validate=400,
        literature_check_search=offline_search,
        progress_hook=_print_progress,
        max_revs=_MAX_REVS,
    )
    elapsed = time.time() - t0
    _print_progress(f"Scan complete in {elapsed:.1f}s — {len(matches)} survivors")

    n_fresh = sum(1 for m in matches if m.is_literature_fresh())
    n_published = sum(1 for m in matches if m.literature_check.status == "published")
    n_inconclusive = sum(1 for m in matches if m.literature_check.status == "inconclusive")
    _print_progress(
        f"Literature check breakdown: {n_fresh} not-found / "
        f"{n_published} published / {n_inconclusive} inconclusive"
    )
    if matches:
        cr = [m.closure.closure_residual_kms for m in matches]
        fc = [m.closure.flyby_continuity_max_dv_kms for m in matches]
        _print_progress(
            f"Residuals: min_closure={min(cr):.4f} km/s "
            f"min_flyby_cont={min(fc):.4f} km/s (gate=0.10)"
        )
        best = matches[0]
        _print_progress(
            f"Best (quality_score={best.quality_score():.3f}): "
            f"sequence={best.candidate.sequence} launch={best.candidate.launch_epoch_utc} "
            f"closure={best.closure.closure_residual_kms:.4f} km/s "
            f"flyby_cont={best.closure.flyby_continuity_max_dv_kms:.4f} km/s "
            f"V_inf_match={best.vinf_match_residual_kms:.4f} km/s "
            f"literature={best.literature_check.status}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for m in matches:
            f.write(json.dumps(precursor_match_to_jsonl_record(m)) + "\n")
    _print_progress(f"Wrote {output_path.relative_to(_DATA_DIR.parent)}")


if __name__ == "__main__":
    main()
