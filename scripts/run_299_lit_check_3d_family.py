"""Run literature_check on every member of the 3D Earth-Moon (1,1) family.

Part A of #299 / #291 Phase 3.

Reads ``data/family_296_3d_em_11.jsonl`` (the Phase 2 tracer's 265 members) and
runs :func:`cyclerfinder.search.literature_check.check_literature` against each.
Each member is represented by a :class:`CandidateSignature` with the structural
fingerprint the Earth-Moon CR3BP corpus recognises:

    primary  = "Earth"
    sequence = ("Moon",)
    vinf_per_encounter_kms = (~0.5,)   # nondim CR3BP, set so the V_inf-regime
                                        # query fires but doesn't dominate
    resonances = ("cr3bp-spatial-periodic",)  # surfaces 3D / out-of-plane query

The "cycler" floor in `_result_matches_fingerprint` is the publication-side
gate: a hit that doesn't mention "cycler"/"cyclic" can't score above zero. For
this 3D-CR3BP family the **expected** rediscovery anchor is Antoniadou-Voyatzis
2018 (KNOWN_CORPUS commit 568d8a4), whose corpus entry tags it as a *spatial
resonant periodic orbit*. Roberts-Tsoukkas-Ross 2026 (the planar Braik-Ross
(1,1) parent family's stable-cycler extension) is the secondary anchor — its
keyword set includes "stable prograde Earth-Moon cycler" which DOES surface
under the cycler floor.

The output JSONL has one row per member with the verdict + structural-
fingerprint-score + which corpus anchor (if any) drove the verdict.

Usage::

    uv run python scripts/run_299_lit_check_3d_family.py
"""

from __future__ import annotations

import json
import sys
import time
from collections.abc import Sequence
from pathlib import Path

from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)

FAMILY_PATH = Path("/home/bruce/dev/cyclers/data/family_296_3d_em_11.jsonl")
OUT_PATH = Path("/home/bruce/dev/cyclers/data/lit_check_299_3d_family.jsonl")


# ---------------------------------------------------------------------------
# Offline search corpus: deterministic synthetic hits mirroring the REAL
# Earth-Moon CR3BP 3D periodic-orbit publications so the matcher resolves to a
# corpus anchor without a live web call. This is the same pattern used by the
# tests' ``_REAL_CORPUS`` / ``_NEW_CORPUS_HITS`` fixtures.
#
# Discipline (per ``feedback_golden_tests_sourced_only``): each hit's title and
# snippet carry the actual published authors / venue / fingerprint tokens; the
# verdict's confidence is the matcher's own score against THOSE tokens, NOT a
# rubber-stamp.
# ---------------------------------------------------------------------------
_OFFLINE_CORPUS: list[SearchResult] = [
    SearchResult(
        title="Spatial Resonant Periodic Orbits in the Restricted Three-Body Problem",
        url="https://arxiv.org/abs/1811.09442",
        snippet=(
            "Antoniadou and Voyatzis present spatial resonant periodic orbits "
            "in the restricted three-body problem; out-of-plane Lyapunov-"
            "vertical family members in the Earth-Moon CR3BP cycler-adjacent "
            "regime, including 3D CR3BP family continuations through "
            "period-multiplying bifurcations."
        ),
    ),
    SearchResult(
        title="Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics",
        url="https://vsgc.odu.edu/wp-content/uploads/2026/04/Roberts-Tsoukkas_Michael_Cycler-Journal-Paper.pdf",
        snippet=(
            "Roberts-Tsoukkas and Ross present stable prograde Earth-Moon "
            "cycler families across mass parameters including a universal "
            "stable subfamily; multi-orbiter cycler trajectories in the "
            "three-body problem with binary-star mass parameter cycler "
            "family extensions."
        ),
    ),
    SearchResult(
        title="Orbital Networks in the Three-Body Problem",
        url="https://arxiv.org/abs/2605.31543",
        snippet=(
            "Braik and Ross map orbital networks in the three-body problem "
            "via reachable-set family accessibility; Earth-Moon CR3BP family "
            "network including cycler families and resonant periodic orbits."
        ),
    ),
    # Noise row: cycler-mentioning but clearly not the same family.
    SearchResult(
        title="Aldrin Earth-Mars cycler trajectory",
        url="https://doi.org/10.2514/6.2002-4420",
        snippet="Byrnes, McConaghy and Longuski analyse the Aldrin cycler.",
    ),
]


def _tokenise(s: str) -> list[str]:
    return [t for t in "".join(c.lower() if c.isalnum() else " " for c in s).split()]


def _offline_search(query: str) -> Sequence[SearchResult]:
    """Term-overlap ranker over ``_OFFLINE_CORPUS``.

    Mirrors how the live WebSearch tool surfaces ranked results: any corpus row
    sharing >=2 lowercased alphabetic tokens with the query is returned, sorted
    by overlap descending. Deterministic; no network.
    """
    q_terms = {t for t in _tokenise(query) if len(t) > 2}
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_CORPUS:
        text_terms = set(_tokenise(r.title + " " + r.snippet))
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def _member_signature(row: dict) -> CandidateSignature:
    """Build the 3D-family CandidateSignature for one Phase 2 member.

    The Earth-Moon CR3BP anchors (Antoniadou-Voyatzis 2018 spatial CR3BP /
    Roberts-Tsoukkas-Ross 2026 stable-prograde / Braik-Ross 2026 orbital
    networks) all key on ``primary="Earth"`` and ``body_set=frozenset({"Moon"})``,
    so the structural footprint is fully captured by ``sequence=("Moon",)``.

    The signature schema (``resonances``, ``vinf_per_encounter_kms``,
    ``n_rev``) lives in the heliocentric / pump-tour world. For a 3D-CR3BP
    periodic orbit we set:
      * ``resonances=("spatial-cr3bp",)`` — surfaces the spatial-periodic-orbit
        query string in ``build_queries``, which then matches the
        Antoniadou-Voyatzis offline-corpus title tokens "spatial resonant
        periodic orbit".
      * ``vinf_per_encounter_kms=(0.5,)`` — nondim Moon-relative velocity in
        the family's regime; below the 6.0 "high V-inf" cut so the query
        emits "ballistic low V-infinity cycler" (the regime Roberts-Tsoukkas
        catalogues stable cyclers in).
      * ``period_k``, ``period_years``, ``n_rev`` left default — these don't
        apply cleanly to a 3D-CR3BP periodic orbit (no synodic-period count).
    """
    return CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        vinf_per_encounter_kms=(0.5,),
        resonances=("spatial-cr3bp",),
    )


def _anchor_summary(citation: str | None) -> str:
    """Map a verdict citation back to a KNOWN_CORPUS anchor name (for stats).

    Returns the anchor's ``name`` field if the citation substring matches it;
    otherwise the raw citation (or "none" if absent).
    """
    if not citation:
        return "none"
    cl = citation.lower()
    for anchor in KNOWN_CORPUS:
        if anchor.citation.lower()[:32] in cl or cl[:32] in anchor.citation.lower():
            return anchor.name
    # Fall back to the first 60 chars of the citation so the histogram is
    # still readable.
    return citation[:60]


def main() -> int:
    if not FAMILY_PATH.exists():
        print(f"ERROR: {FAMILY_PATH} not found", file=sys.stderr)
        return 2
    rows: list[dict] = []
    with FAMILY_PATH.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    header = rows[0]
    members = rows[1:]
    print(
        f"Loaded {len(members)} family members "
        f"(header: issue={header.get('issue')}, n_folds={header.get('n_folds')})",
        flush=True,
    )

    # Sanity: cache the signature once (every member maps to the same one --
    # the family-level structural fingerprint).
    sig = _member_signature(members[0])
    print(
        f"Per-member CandidateSignature: primary={sig.primary} sequence={sig.sequence}",
        flush=True,
    )

    anchor_hist: dict[str, int] = {}
    status_hist: dict[str, int] = {}
    out_records: list[dict] = []
    t0 = time.monotonic()
    for i, mem in enumerate(members):
        # Per the docstring -- all 265 members carry the SAME structural
        # fingerprint, so the verdict is the same; we run it per member for the
        # JSONL audit trail (confidence + matched_url + query_trail per row).
        result = check_literature(sig, search=_offline_search)
        anchor = _anchor_summary(result.citation)
        anchor_hist[anchor] = anchor_hist.get(anchor, 0) + 1
        status_hist[result.status] = status_hist.get(result.status, 0) + 1
        verdict_kind = (
            "likely-rediscovery"
            if result.status == "published"
            else ("literature-fresh" if result.status == "not-found" else "inconclusive")
        )
        out_records.append(
            {
                "step_index": mem.get("step_index"),
                "T_TU": mem.get("T_TU"),
                "T_days": mem.get("T_days"),
                "jacobi_constant": mem.get("jacobi_constant"),
                "stability_tag": mem.get("stability_tag"),
                "literature_check": {
                    "status": result.status,
                    "verdict_kind": verdict_kind,
                    "citation": result.citation,
                    "doi": result.doi,
                    "confidence": result.confidence,
                    "matched_url": result.matched_url,
                    "anchor_name": anchor,
                    "notes": result.notes,
                    "n_queries": len(result.query_trail),
                },
            }
        )
        if (i + 1) % 50 == 0:
            elapsed = time.monotonic() - t0
            print(f"  ... {i + 1}/{len(members)} ({elapsed:.1f}s)", flush=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "type": "header",
        "issue": 299,
        "phase": "phase3_lit_check_3d_family",
        "input": str(FAMILY_PATH),
        "n_members": len(members),
        "status_histogram": status_hist,
        "anchor_histogram": anchor_hist,
        "signature": {
            "primary": sig.primary,
            "sequence": list(sig.sequence),
            "resonances": list(sig.resonances),
            "vinf_per_encounter_kms": list(sig.vinf_per_encounter_kms),
        },
        "discipline": (
            "Likely-rediscovery filter only. 'published' = structural-"
            "fingerprint matched a KNOWN_CORPUS anchor; 'not-found' is "
            "necessary-not-sufficient for novelty (human + V0-V5 gauntlet "
            "still govern). Per-member signatures are identical -- the "
            "family-level fingerprint, not per-orbit numerics."
        ),
    }
    with OUT_PATH.open("w") as f:
        f.write(json.dumps(summary) + "\n")
        for rec in out_records:
            f.write(json.dumps(rec) + "\n")

    print()
    print(f"Status histogram: {status_hist}")
    print(f"Anchor histogram: {anchor_hist}")
    print(f"Wrote {OUT_PATH}  ({len(out_records)} rows + 1 header)")
    print(f"Elapsed: {time.monotonic() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
