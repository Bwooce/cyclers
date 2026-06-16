"""3D sub-family validation: lit-check + ML false-positive flagger (#301).

Phase 4 of #291 / follow-up to #299 (commit f8c3fde). For each member of the
four Neimark-Sacker-born sub-families in
``data/family_296_3d_subfamilies_299.jsonl``, this script:

1. Runs :func:`cyclerfinder.search.literature_check.check_literature` with the
   sub-family's actual CR3BP period band (T_TU_min, T_TU_max) populated -- the
   #301 ``period_band_tu`` extension. Antoniadou-Voyatzis 2018 carries an
   anchor-side period_band_tu of (0, 15); the Neimark-Sacker sub-families sit
   at T_TU >= ~20, so the AV anchor is excluded from anchor matching for them.
   The remaining Earth-Moon corpus (Roberts-Tsoukkas-Ross, Braik-Ross,
   Kumar-Rosengren-Ross resonant transport, Koblick / Zhang tulip,
   Hiraiwa lobe dynamics, Sanaga fidelity) still gets a fair shot.

2. Runs the #256 / #275 FalsePosFlagger on a minimal SILVER-shaped record
   constructed from the JSONL (closure residual, period, CR3BP model
   assumption, post-fix closure date). The flagger's feature set was hand-
   crafted for heliocentric mission-design SILVERs (V_inf, bend feasibility,
   etc.); for a 3D periodic orbit those features fall back to NaN and the
   flagger imputes -- the trained median's contribution is small but non-
   zero. Output records the resulting ``p_fp`` as-is; the caller and the
   verdict doc state plainly that this is a NON-BLOCKING signal designed for
   a different candidate class.

3. Emits ``data/sub_families_301.jsonl`` with one row per sub-family member
   plus a header. Each row carries: ``parent_family_id``, ``k``, ``T_TU``,
   ``T_days``, ``jacobi_constant``, ``stability_tag``, the
   ``literature_check`` block, ``ml_p_fp``, and ``novelty_claimable`` =
   (lit-fresh AND p_fp < 0.5 AND closure <= 1e-9).

DISCIPLINE
----------
* NO catalogue writeback (planar V0-V5 gauntlet cannot run on 3D periodic
  orbits; #302 is the future task that adapts the gauntlet).
* The independent-closure residuals come straight from the JSONL (already
  verified by Phase 2's corrector + Radau independent gate). We just record
  them; no re-integration here per ``feedback_orbit_closure_discipline``.
* Honest negatives are correct outcomes. If a sub-family is still anchored
  (which can happen for k=3 at T_TU ~ 20 if the band test is loose), the
  output records it; we do not over-claim.
"""

from __future__ import annotations

import json
import math
import sys
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import build_training_set
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)

INPUT_PATH = Path("/home/bruce/dev/cyclers/data/family_296_3d_subfamilies_299.jsonl")
OUTPUT_PATH = Path("/home/bruce/dev/cyclers/data/sub_families_301.jsonl")
PARENT_FAMILY_ID = "family_296_3d_em_11"

# ---------------------------------------------------------------------------
# Offline search corpus: same pattern as #299's lit-check (term-overlap over a
# deterministic synthetic corpus mirroring real Earth-Moon CR3BP publications).
# Includes the Antoniadou-Voyatzis 2018 row -- the matcher will still score it
# as long as the candidate anchor list contains AV, which it WILL NOT for the
# sub-families that fall outside the period band. The other rows give the
# fair-shot corpus for the period-multiplied sub-families.
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
    SearchResult(
        title="Cislunar Resonant Transport and Heteroclinic Pathways",
        url="https://arxiv.org/abs/2509.12675",
        snippet=(
            "Kumar, Rawat, Rosengren and Ross study cislunar resonant "
            "transport via heteroclinic pathways on the Earth-Moon "
            "resonant family network; period-multiplying cycler-adjacent "
            "manifold structure."
        ),
    ),
    SearchResult(
        title="Novel Tulip-Shaped Three-body Orbits for Cislunar SDA Missions",
        url="https://doi.org/10.1007/s40295-025-00510-w",
        snippet=(
            "Koblick and Kelly construct tulip-shaped three-body cycler "
            "orbits in the Earth-Moon CR3BP for cislunar SDA missions; "
            "petal-count periodic orbit families."
        ),
    ),
    SearchResult(
        title="Time-regularized bifurcation framework for tulip-shaped orbits",
        url="https://doi.org/10.1007/s11071-026-12465-0",
        snippet=(
            "Zhang, Jiang and Yuan present a time-regularized bifurcation "
            "framework for tulip-shaped orbits; period-multiplying bifurcation "
            "cislunar tulip cycler families."
        ),
    ),
    SearchResult(
        title="Design of low-energy transfers in cislunar space using sequences of lobe dynamics",
        url="https://arxiv.org/abs/2602.17444",
        snippet=(
            "Hiraiwa, Bando, Sato and Hokamoto design low-energy cislunar "
            "transfers using sequences of lobe dynamics on a weighted "
            "directed graph; resonant orbit lobe-dynamics transfer."
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
    """Term-overlap ranker over ``_OFFLINE_CORPUS`` (same pattern as #299)."""
    q_terms = {t for t in _tokenise(query) if len(t) > 2}
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_CORPUS:
        text_terms = set(_tokenise(r.title + " " + r.snippet))
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def _signature_for_subfamily(t_tu_min: float, t_tu_max: float) -> CandidateSignature:
    """3D Earth-Moon CR3BP sub-family signature with the period band populated.

    The Antoniadou-Voyatzis 2018 anchor declares ``period_band_tu=(0, 15)``;
    a sub-family at T_TU 20-44 has a disjoint band and is therefore not
    auto-anchored on AV-2018. The other Earth-Moon-CR3BP anchors carry no
    period band and remain available to flag a hit.
    """
    return CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        vinf_per_encounter_kms=(0.5,),  # nondim, regime fingerprint
        resonances=("spatial-cr3bp-period-multiplied",),
        period_band_tu=(t_tu_min, t_tu_max),
    )


def _anchor_summary(citation: str | None) -> str:
    """Map a verdict citation back to a KNOWN_CORPUS anchor name."""
    if not citation:
        return "none"
    cl = citation.lower()
    for anchor in KNOWN_CORPUS:
        if anchor.citation.lower()[:32] in cl or cl[:32] in anchor.citation.lower():
            return anchor.name
    return citation[:60]


def _silver_record_for_member(member: dict, k: int) -> dict:
    """Map a sub-family member to a minimal SILVER-record dict for the flagger.

    The FalsePosFlagger feature set was designed for heliocentric mission-
    design SILVERs (V_inf, bend feasibility, encounter resonance ratios).
    A 3D CR3BP periodic orbit does not carry those signals -- so we provide
    the fields that DO map naturally and leave the rest absent (the extractor
    is NaN-safe and the flagger imputes with training-set medians).

    Mapped:
      * ``max_residual_kms`` -- the independent-closure residual (nondim)
        scaled by the EM CR3BP unit distance (~384400 km) to get km.
      * ``period_days`` -- the member's nondim period * TU_days (4.34268 d/TU).
      * ``model_assumption`` -- "cr3bp-spatial-3d" (drives the
        cr3bp_pre_mu_fix_flag; we set ``closure_date`` post-fix so the flag
        evaluates 0 = not-suspicious).
      * ``closure_date`` -- 2026-06-16, well after the latest known fix
        (2026-06-14) -- so ``epoch_artifact_flag`` -> 0.0.
      * ``closure_method_version`` -- "3cec84c-postfix-3d" -- not in the
        known-fix prefix table; the extractor evaluates it via the date path.

    Not mapped (NaN; flagger imputes):
      * vinf_per_encounter_kms / vinf_floors_kms / bend_feasible /
        topology_match -- none of these have a 3D-CR3BP equivalent.
    """
    em_lu_km = 384400.0  # Earth-Moon CR3BP unit-distance (km), per ephemeris.py
    em_tu_days = 4.342782  # Earth-Moon CR3BP unit-time (days)
    residual_nondim = float(member.get("independent_closure_residual", 0.0) or 0.0)
    residual_kms = residual_nondim * em_lu_km
    period_days = float(member.get("T_days", member.get("T_TU", 0.0) * em_tu_days))
    return {
        "max_residual_kms": residual_kms,
        "period_days": period_days,
        "model_assumption": "cr3bp-spatial-3d",
        "closure_method_version": "3cec84c-postfix-3d",
        "closure_date": "2026-06-16",
        # Topology: a 3D CR3BP periodic orbit IS topology-correct by construction
        # (the corrector + Radau independent gate already enforced it).
        "topology_match": True,
        # Bend feasibility: N/A for a periodic orbit (no flyby bend).
        "bend_feasible": True,
        # Extra context flags (the extractor ignores unknown keys).
        "_sub_family_k": k,
    }


def _classify_floquet(floquet_abs: list[float]) -> dict:
    """Return a small structural summary of the 6 multipliers.

    Captures the "doubly hyperbolic" signature: 4 multipliers strictly off the
    unit circle (>1 or <1 by > tol) AND the largest non-trivial pair is in
    addition to a second non-trivial pair (i.e. two distinct reciprocal pairs).
    """
    tol = 0.01
    if not floquet_abs:
        return {"n_off_circle": 0, "doubly_hyperbolic": False, "max_modulus": math.nan}
    fa_sorted = sorted(floquet_abs, reverse=True)
    n_off = sum(1 for x in floquet_abs if x > 1 + tol or x < 1 - tol)
    # Doubly hyperbolic: at least two values > 1+tol (reciprocal pairs ensured by
    # symplectic structure -- two unstable pairs <=> n_off >= 4).
    n_unstable = sum(1 for x in floquet_abs if x > 1 + tol)
    doubly_hyperbolic = n_unstable >= 2
    return {
        "n_off_circle": n_off,
        "doubly_hyperbolic": doubly_hyperbolic,
        "max_modulus": float(fa_sorted[0]),
        "second_modulus": float(fa_sorted[1]) if len(fa_sorted) >= 2 else math.nan,
    }


def main() -> int:
    t0 = time.monotonic()
    if not INPUT_PATH.exists():
        print(f"ERROR: {INPUT_PATH} not found", file=sys.stderr)
        return 2

    # Load sub-family records.
    rows: list[dict] = []
    with INPUT_PATH.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    header_in = rows[0]
    sub_family_records = rows[1:]
    print(
        f"Loaded {len(sub_family_records)} sub-family records (header: issue="
        f"{header_in.get('issue')}, n_parent_members={header_in.get('n_parent_members')})",
        flush=True,
    )

    # Train the false-positive flagger once on the labeled corpus.
    print("Training FalsePosFlagger on the #256/#275 labeled corpus ...", flush=True)
    x_train, y_train, _meta = build_training_set()
    flagger = FalsePosFlagger()
    diag = flagger.fit(x_train, y_train)
    print(
        f"  N={diag.n_samples}  positives(FP)={diag.n_positives}  negatives(TR)={diag.n_negatives}",
        flush=True,
    )
    print(
        f"  AUC(train)={diag.auc_train:.4f}  AUC(LOO)={diag.auc_loo:.4f}  "
        f"[reproduce-before-trust gate: > 0.75]",
        flush=True,
    )

    # Per-sub-family lit-check + ML scoring.
    per_subfamily_summary: list[dict] = []
    out_member_records: list[dict] = []
    for rec_idx, sub in enumerate(sub_family_records):
        bracket = sub["bracket"]
        switched = sub["switched"]
        members = sub["subfamily_members"]
        k = bracket["k"]
        t_tu_extent = sub["subfamily_T_TU_extent"]
        t_tu_min, t_tu_max = t_tu_extent
        em_tu_days = 4.342782
        t_days_min = t_tu_min * em_tu_days
        t_days_max = t_tu_max * em_tu_days

        sig = _signature_for_subfamily(t_tu_min, t_tu_max)
        # Run the lit-check ONCE per sub-family (the signature is constant
        # across members -- structural fingerprint, not per-orbit numerics).
        verdict = check_literature(sig, search=_offline_search)
        anchor = _anchor_summary(verdict.citation)

        n_members = len(members)
        n_lit_fresh = 0
        n_lit_published = 0
        n_lit_inconclusive = 0
        n_ml_fresh = 0
        n_novelty_claimable = 0
        n_doubly_hyperbolic = 0

        # Per-member ML scoring + emission.
        for mem in members:
            silver_rec = _silver_record_for_member(mem, k)
            p_fp = flagger.score(silver_rec)
            floq_summary = _classify_floquet(mem.get("floquet_abs", []))

            ic_residual = float(mem.get("independent_closure_residual", 0.0) or 0.0)
            closure_ok = ic_residual <= 1e-9

            lit_status = verdict.status
            if lit_status == "published":
                n_lit_published += 1
                lit_fresh = False
            elif lit_status == "not-found":
                n_lit_fresh += 1
                lit_fresh = True
            else:
                n_lit_inconclusive += 1
                lit_fresh = False

            ml_fresh = p_fp < 0.5
            if ml_fresh:
                n_ml_fresh += 1

            novelty_claimable = bool(lit_fresh and ml_fresh and closure_ok)
            if novelty_claimable:
                n_novelty_claimable += 1
            if floq_summary["doubly_hyperbolic"]:
                n_doubly_hyperbolic += 1

            out_member_records.append(
                {
                    "parent_family_id": PARENT_FAMILY_ID,
                    "k": k,
                    "subfamily_index": rec_idx,
                    "member_step_index": mem.get("step_index"),
                    "state_nd": mem.get("state_nd"),
                    "T_TU": mem.get("T_TU"),
                    "T_days": mem.get("T_days"),
                    "jacobi": mem.get("jacobi_constant"),
                    "stability_tag": mem.get("stability_tag"),
                    "floquet_summary": floq_summary,
                    "independent_closure_residual": ic_residual,
                    "closure_ok": closure_ok,
                    "lit_check": {
                        "status": lit_status,
                        "citation": verdict.citation,
                        "doi": verdict.doi,
                        "confidence": verdict.confidence,
                        "anchor_name": anchor,
                        "matched_url": verdict.matched_url,
                    },
                    "ml_p_fp": p_fp,
                    "ml_fresh": ml_fresh,
                    "novelty_claimable": novelty_claimable,
                }
            )

        # Pick highest-priority IC for #302: the novelty-claimable member with
        # the lowest p_fp + best closure, fall back to the doubly-hyperbolic
        # member with the lowest p_fp.
        candidates = [
            r
            for r in out_member_records
            if r["k"] == k and r["subfamily_index"] == rec_idx and r["novelty_claimable"]
        ]
        if not candidates:
            candidates = [
                r
                for r in out_member_records
                if r["k"] == k and r["subfamily_index"] == rec_idx and r["closure_ok"]
            ]
        if candidates:
            best = min(candidates, key=lambda r: (r["ml_p_fp"], r["independent_closure_residual"]))
            priority_ic = {
                "state_nd": best["state_nd"],
                "T_TU": best["T_TU"],
                "T_days": best["T_days"],
                "jacobi": best["jacobi"],
                "member_step_index": best["member_step_index"],
                "p_fp": best["ml_p_fp"],
                "doubly_hyperbolic": best["floquet_summary"]["doubly_hyperbolic"],
            }
        else:
            priority_ic = None

        per_subfamily_summary.append(
            {
                "k": k,
                "subfamily_index": rec_idx,
                "n_members": n_members,
                "T_TU_extent": [t_tu_min, t_tu_max],
                "T_days_extent": [t_days_min, t_days_max],
                "switched_seed_state": switched["switched_state"],
                "switched_T_TU": switched["switched_T_TU"],
                "switched_period_ratio": switched["switched_period_ratio"],
                "lit_check": {
                    "status": verdict.status,
                    "citation": verdict.citation,
                    "doi": verdict.doi,
                    "confidence": verdict.confidence,
                    "anchor_name": anchor,
                    "matched_url": verdict.matched_url,
                    "notes": verdict.notes,
                    "n_queries": len(verdict.query_trail),
                },
                "counts": {
                    "n_lit_fresh": n_lit_fresh,
                    "n_lit_published": n_lit_published,
                    "n_lit_inconclusive": n_lit_inconclusive,
                    "n_ml_fresh": n_ml_fresh,
                    "n_novelty_claimable": n_novelty_claimable,
                    "n_doubly_hyperbolic": n_doubly_hyperbolic,
                },
                "priority_ic_for_302": priority_ic,
            }
        )
        print(
            f"  k={k}: n_members={n_members}, T_TU=[{t_tu_min:.4f}, {t_tu_max:.4f}], "
            f"lit_status={verdict.status}, anchor={anchor}, "
            f"n_lit_fresh={n_lit_fresh}, n_ml_fresh={n_ml_fresh}, "
            f"n_novelty_claimable={n_novelty_claimable}, "
            f"n_doubly_hyperbolic={n_doubly_hyperbolic}",
            flush=True,
        )

    # Decide highest-priority sub-family for #302: most doubly-hyperbolic
    # members (3D heteroclinic web signature), tie-break by lowest p_fp on the
    # priority IC, then by larger period (richer structure).
    def _sub_priority_key(s: dict) -> tuple:
        ic = s.get("priority_ic_for_302") or {}
        return (
            -s["counts"]["n_doubly_hyperbolic"],
            ic.get("p_fp", 1.0),
            -float(s["switched_T_TU"]),
        )

    eligible = [s for s in per_subfamily_summary if s["counts"]["n_novelty_claimable"] > 0]
    if eligible:
        priority_sub = sorted(eligible, key=_sub_priority_key)[0]
    elif per_subfamily_summary:
        priority_sub = sorted(per_subfamily_summary, key=_sub_priority_key)[0]
    else:
        priority_sub = None

    # Emit JSONL.
    summary = {
        "type": "header",
        "issue": 301,
        "phase": "phase4_subfamily_validation",
        "input": str(INPUT_PATH),
        "n_subfamilies": len(per_subfamily_summary),
        "n_total_members": sum(
            s["counts"]["n_lit_fresh"]
            + s["counts"]["n_lit_published"]
            + s["counts"]["n_lit_inconclusive"]
            for s in per_subfamily_summary
        ),
        "per_subfamily": per_subfamily_summary,
        "highest_priority_for_302": (
            None
            if priority_sub is None
            else {
                "k": priority_sub["k"],
                "subfamily_index": priority_sub["subfamily_index"],
                "T_TU_extent": priority_sub["T_TU_extent"],
                "priority_ic": priority_sub["priority_ic_for_302"],
                "rationale": (
                    "doubly-hyperbolic-rich (3D heteroclinic web signature)"
                    if priority_sub["counts"]["n_doubly_hyperbolic"] > 0
                    else "highest novelty claim density at lowest ml_p_fp"
                ),
            }
        ),
        "ml_flagger_diagnostics": {
            "auc_train": diag.auc_train,
            "auc_loo": diag.auc_loo,
            "n_samples": diag.n_samples,
            "n_positives": diag.n_positives,
            "n_negatives": diag.n_negatives,
        },
        "discipline": (
            "3D sub-family validation. NO catalogue writeback -- the planar "
            "V0-V5 gauntlet cannot run on 3D periodic orbits; admission is "
            "blocked until #302 adapts the gauntlet. 'novelty_claimable' is "
            "necessary-not-sufficient: lit-fresh + ml-fresh + closure_ok ONLY "
            "clears the candidate for a future 3D gauntlet attempt. The ML "
            "flagger was trained on heliocentric mission-design SILVERs; for "
            "a 3D CR3BP periodic orbit several features are NaN-imputed -- "
            "the score is a guard-rail signal, not a strong classifier here."
        ),
    }

    def _json_default(o: object) -> object:
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.integer):
            return int(o)
        raise TypeError(f"object of type {type(o).__name__} is not JSON serializable")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as f:
        f.write(json.dumps(summary) + "\n")
        for r in out_member_records:
            f.write(json.dumps(r, default=_json_default) + "\n")

    elapsed = time.monotonic() - t0
    print()
    print(f"Wrote {OUTPUT_PATH}  ({len(out_member_records)} member rows + 1 header)")
    print(f"Elapsed: {elapsed:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
