"""Prioritized repeated-moon campaign at Saturn and Uranus (#285).

The 5-tier prioritizer stack from #282 is operational on patched-conic legs
through Tier 0 (Zhang-Topputo NN); Tiers 1-5 require CR3BP representatives
which #264-style repeated-moon enumeration does not emit (#282 architectural
seam, documented in :mod:`five_tier_prioritizer`). This module wires the
applicable scorer into the discovery-campaign engine for Saturn / Uranus moon
systems where the Tisserand-graph existence priors are strongest.

Pipeline composition::

    RepeatedMoonTarget          -> enumerate candidates (sequence x n_rev)
      .close()                  -> Lambert-leg closure, residual_kms
      [residual_kms < gate]     -> SILVER survivor
      legs_from_repeated_moon_candidate -> reconstruct leg geometry (SI)
      FiveTierPrioritizer.score_candidate_legs -> Tier-0 NN per-leg dV
      _dop853_cross_check       -> independent two-body propagation per leg
      check_literature(offline) -> rediscovery filter (necessary-not-sufficient)
      FalsePosFlagger           -> p_fp (trained on labelled SILVER corpus)
      verdict in {SILVER, BRONZE, REJECT}

DISCIPLINE (golden, non-negotiable):

* **NO catalogue writeback.** SILVER survivors emit to a JSONL file ONLY.
  Promotion to the catalogue is the gauntlet's job (#274), not this scan's.
* **NO novelty claims.** A "literature-fresh" survivor here is a guard-chain
  survivor, not a discovery. The V0-V5 gauntlet is the gate.
* **Independent cross-check is mandatory** (orbit-closure discipline,
  ``feedback_orbit_closure_discipline``). The Lambert solver IS the closure
  solver, so re-running Lambert is NOT independent. We re-propagate every
  leg's departure state with :mod:`scipy.integrate.solve_ivp` DOP853 at
  rtol=atol=1e-12 (different solver, same two-body physics) and confirm the
  arrival state matches the Lambert prediction within ``cross_check_tol_km``.

VERDICT POLICY (per row):

* ``SILVER`` -- Lambert closure < gate, NN admits all legs (or NN unavailable
  AND falls back), DOP853 cross-check passes, lit-check is ``not-found`` or
  ``inconclusive``, ``p_fp <= 0.75``. Routes to the gauntlet.
* ``BRONZE`` -- Lambert closure < gate but one downstream guard caveats
  (NN rejects, lit-check ``published``, or ``p_fp > 0.75``). Recorded but
  NOT a gauntlet entry; surfaces for human triage.
* ``REJECT`` -- DOP853 cross-check FAILS (numerical disagreement >
  tolerance), i.e. the closure is suspect. Recorded as a negative for audit.
"""

from __future__ import annotations

import json
import math
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import build_training_set
from cyclerfinder.search.discovery_campaign import (
    Candidate,
    ClosureResult,
    RepeatedMoonTarget,
    moon_cycler_signature_hash,
)
from cyclerfinder.search.five_tier_prioritizer import (
    FiveTierPrioritizer,
    PatchedConicLeg,
    legs_from_repeated_moon_candidate,
)
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)


def offline_corpus_search(query: str) -> Sequence[SearchResult]:
    """Deterministic offline literature backend (inlined from #261 driver).

    Mirrors ``scripts.literature_check_review_queue.offline_corpus_search``
    verbatim: for each curated :data:`KNOWN_CORPUS` anchor whose author/
    keyword/body appears in the query, emit a synthetic hit so the structural
    matcher in :func:`check_literature` can score it. NOT a web search; this
    only re-finds families already in the curated corpus, so a candidate in a
    known family lands ``published`` while novel candidates fall through to
    ``inconclusive``/``not-found`` (never a false offline ``not-found``).
    """
    q = query.lower()
    out: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        hit = any(a.lower() in q for a in anchor.authors) or any(
            kw.lower() in q for kw in anchor.keywords
        )
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


# ---------------------------------------------------------------------------
# Verdict policy thresholds (the #274 gauntlet's V0 gate uses 0.75 for p_fp;
# we match that here so a SILVER from this scan is plumbed for the gauntlet
# without re-thresholding).
# ---------------------------------------------------------------------------

P_FP_SILVER_MAX: float = 0.75
"""``ml_flagger_p_fp <= 0.75`` for SILVER (matches gauntlet V0 gate, #274)."""

CROSS_CHECK_TOL_KM: float = 1.0
"""Per-leg DOP853 vs Lambert arrival-state agreement gate (km).

The Lambert solver is double-precision and integrates a closed-form ellipse;
DOP853 at rtol=atol=1e-12 should agree to numerical noise. We pick a generous
1 km gate (vs SMA ~10^5-10^6 km) so a "REJECT" really means a numerical
disagreement, not a tolerance grind. The actual residual is recorded.
"""


# ---------------------------------------------------------------------------
# Independent DOP853 cross-check (different solver, same two-body physics)
# ---------------------------------------------------------------------------


def _two_body_rhs(_t: float, y: np.ndarray, mu_km3_s2: float) -> np.ndarray:
    """Two-body equations of motion in km, km/s, s, with mu in km^3/s^2."""
    r = y[:3]
    v = y[3:]
    r_norm = float(np.linalg.norm(r))
    if r_norm <= 0.0:
        return np.concatenate([v, np.zeros(3)])
    a = -mu_km3_s2 * r / (r_norm**3)
    return np.concatenate([v, a])


def dop853_cross_check_leg(
    leg: PatchedConicLeg,
    *,
    rtol: float = 1.0e-12,
    atol: float = 1.0e-12,
) -> dict[str, float | bool]:
    """Re-integrate one leg's departure state with DOP853 and compare to Lambert.

    The Lambert solver produces ``(r1, v1)`` (departure) and ``(r2, v2)``
    (arrival) over time-of-flight ``dt_s``. Closure = ballistic two-body
    propagation under the central-body mu. We re-integrate ``(r1, v1)``
    forward by ``dt_s`` with :mod:`scipy.integrate.solve_ivp` DOP853 at
    ``rtol=atol=1e-12`` (independent solver, identical physics) and report
    the position/velocity disagreement with the Lambert-predicted ``(r2, v2)``.

    Units: positions in km, velocities in km/s, mu in km^3/s^2. The leg
    stores SI; we convert internally.

    Returns
    -------
    dict
        ``dr_arrival_km`` -- |r_DOP853 - r2_Lambert| (km),
        ``dv_arrival_km_s`` -- |v_DOP853 - v2_Lambert| (km/s),
        ``converged`` -- bool, ``True`` iff solver succeeded,
        ``passed`` -- bool, ``True`` iff ``dr_arrival_km <= CROSS_CHECK_TOL_KM``.
    """
    km_m = 1000.0
    r1_km = np.asarray(leg.r1_m, dtype=np.float64) / km_m
    v1_km_s = np.asarray(leg.v1_m_s, dtype=np.float64) / km_m
    r2_km = np.asarray(leg.r2_m, dtype=np.float64) / km_m
    v2_km_s = np.asarray(leg.v2_m_s, dtype=np.float64) / km_m
    mu_km3_s2 = leg.mu_m3_s2 / (km_m**3)
    y0 = np.concatenate([r1_km, v1_km_s])
    sol = solve_ivp(
        _two_body_rhs,
        (0.0, float(leg.dt_s)),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        args=(mu_km3_s2,),
        dense_output=False,
    )
    if not sol.success:
        return {
            "dr_arrival_km": float("inf"),
            "dv_arrival_km_s": float("inf"),
            "converged": False,
            "passed": False,
        }
    r_end = sol.y[:3, -1]
    v_end = sol.y[3:, -1]
    dr = float(np.linalg.norm(r_end - r2_km))
    dv = float(np.linalg.norm(v_end - v2_km_s))
    return {
        "dr_arrival_km": dr,
        "dv_arrival_km_s": dv,
        "converged": True,
        "passed": dr <= CROSS_CHECK_TOL_KM,
    }


# ---------------------------------------------------------------------------
# ML flagger: train once on the labelled corpus, score every SILVER candidate
# ---------------------------------------------------------------------------


def trained_flagger() -> FalsePosFlagger:
    """Train a fresh FalsePosFlagger on the labelled-corpus seed (#256 retrain).

    The flagger is a thin logistic regression over hand-crafted features; the
    "trained" state lives entirely in the betas, so a fresh `.fit(...)` at the
    start of each campaign gives the same retrained-as-of-`4ea5992` flagger
    the gauntlet uses. ``score()`` is non-blocking by contract (returns 0.5
    on any feature-extraction failure).
    """
    flagger = FalsePosFlagger()
    x_train, y, _meta = build_training_set()
    flagger.fit(x_train, y)
    return flagger


# ---------------------------------------------------------------------------
# Scorer chain: per-candidate plumbing
# ---------------------------------------------------------------------------


@dataclass
class ScoredCandidate:
    """One scored row (the JSONL row's structured form before json.dumps)."""

    candidate_id: str
    primary: str
    sequence: tuple[str, ...]
    n_rev: tuple[int, ...]
    # Closure
    residual_kms: float
    vinf_per_encounter_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    max_vinf_kms: float
    # Tier-0 NN
    tier0: dict[str, Any]
    # DOP853 independent cross-check
    cross_check: dict[str, Any]
    # Literature novelty
    literature: dict[str, Any]
    # ML false-positive flagger
    ml_flagger_p_fp: float
    # Final verdict
    verdict: str  # "SILVER" | "BRONZE" | "REJECT"
    verdict_reasons: list[str]
    dedup_signature: str

    def as_row(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "primary": self.primary,
            "sequence": list(self.sequence),
            "n_rev": list(self.n_rev),
            "residual_kms": self.residual_kms,
            "vinf_per_encounter_kms": list(self.vinf_per_encounter_kms),
            "tof_days": list(self.tof_days),
            "max_vinf_kms": self.max_vinf_kms,
            "tier0": self.tier0,
            "cross_check": self.cross_check,
            "literature": self.literature,
            "ml_flagger_p_fp": self.ml_flagger_p_fp,
            "verdict": self.verdict,
            "verdict_reasons": list(self.verdict_reasons),
            "dedup_signature": self.dedup_signature,
        }


def _score_one_silver(
    *,
    target: RepeatedMoonTarget,
    candidate: Candidate,
    closure: ClosureResult,
    prioritizer: FiveTierPrioritizer,
    flagger: FalsePosFlagger,
    phase_samples: int,
    cross_check_tol_km: float = CROSS_CHECK_TOL_KM,
) -> ScoredCandidate:
    """Run Tier-0 + DOP853 + literature + ML flagger for one SILVER candidate."""
    n_rev = tuple(candidate.payload["n_rev"])
    legs = legs_from_repeated_moon_candidate(
        primary=target.primary,
        sequence=candidate.sequence,
        n_rev=list(n_rev),
        phase_samples=phase_samples,
    )

    # Tier-0 NN (#282 architecture: only tier that works on patched-conic legs).
    tier0: dict[str, Any]
    if legs is None:
        tier0 = {
            "status": "no-feasible-phasing",
            "per_leg": [],
            "tier0_max_dv_kms": None,
            "tier0_sum_dv_kms": None,
            "tier0_all_admitted": False,
            "tier0_any_inference_failed": False,
            "n_legs": 0,
        }
    else:
        t0_stats = prioritizer.score_candidate_legs(legs)
        tier0 = {
            "status": "ok",
            "per_leg": [
                {
                    "label_from": p["label_from"],
                    "label_to": p["label_to"],
                    "dv_kms": p["tier0_predicted_dv_kms"],
                    "tof_days": p["tier0_predicted_tof_days"],
                    "admitted": p["tier0_admitted"],
                    "model_available": p["tier0_model_available"],
                    "fallback_used": p["tier0_fallback_used"],
                }
                for p in t0_stats["per_leg"]
            ],
            "tier0_max_dv_kms": t0_stats["tier0_max_dv_kms"],
            "tier0_sum_dv_kms": t0_stats["tier0_sum_dv_kms"],
            "tier0_all_admitted": t0_stats["tier0_all_admitted"],
            "tier0_any_inference_failed": t0_stats["tier0_any_inference_failed"],
            "n_legs": t0_stats["n_legs"],
        }

    # Independent DOP853 cross-check (the orbit-closure discipline gate).
    cross_check: dict[str, Any]
    if legs is None:
        cross_check = {
            "status": "skipped: no-feasible-phasing",
            "max_dr_arrival_km": None,
            "max_dv_arrival_km_s": None,
            "all_passed": False,
            "per_leg": [],
        }
    else:
        per_leg_xc: list[dict[str, Any]] = []
        max_dr = 0.0
        max_dv = 0.0
        all_passed = True
        for leg in legs:
            r = dop853_cross_check_leg(leg, rtol=1e-12, atol=1e-12)
            per_leg_xc.append(
                {
                    "label_from": leg.label_from,
                    "label_to": leg.label_to,
                    **r,
                }
            )
            max_dr = max(max_dr, float(r["dr_arrival_km"]))
            max_dv = max(max_dv, float(r["dv_arrival_km_s"]))
            if not r["passed"]:
                all_passed = False
        cross_check = {
            "status": "ok",
            "max_dr_arrival_km": max_dr,
            "max_dv_arrival_km_s": max_dv,
            "all_passed": all_passed,
            "per_leg": per_leg_xc,
            "tol_km": cross_check_tol_km,
        }

    # Literature novelty (offline backend).
    sig = CandidateSignature(
        primary=target.primary,
        sequence=tuple(candidate.sequence),
        period_k=len(candidate.sequence) - 1,
        vinf_per_encounter_kms=tuple(closure.vinf_per_encounter_kms),
        n_rev=n_rev,
    )
    try:
        lit_result = check_literature(sig, search=offline_corpus_search)
        literature = {
            "status": lit_result.status,
            "citation": lit_result.citation,
            "doi": lit_result.doi,
            "confidence": lit_result.confidence,
            "matched_url": lit_result.matched_url,
            "backend": "offline_corpus_search",
        }
    except Exception as exc:
        literature = {
            "status": "inconclusive",
            "citation": None,
            "doi": None,
            "confidence": 0.0,
            "matched_url": None,
            "backend": "offline_corpus_search",
            "error": repr(exc),
        }

    # ML flagger (non-blocking; 0.5 on any failure).
    silver_record_like = {
        "primary": target.primary,
        "sequence": list(candidate.sequence),
        "n_rev": list(n_rev),
        "vinf_per_encounter_kms": list(closure.vinf_per_encounter_kms),
        "tof_days": list(closure.tof_days),
        "verdict_audit": {"residual_kms": closure.residual_kms, "primary": target.primary},
        "max_vinf_kms": (
            max(closure.vinf_per_encounter_kms) if closure.vinf_per_encounter_kms else 0.0
        ),
        "bend_feasible": True,
    }
    p_fp = flagger.score(silver_record_like)

    # Verdict policy.
    reasons: list[str] = []
    verdict = "SILVER"
    if not cross_check.get("all_passed"):
        verdict = "REJECT"
        reasons.append(
            f"DOP853 cross-check failed: max_dr_arrival_km="
            f"{cross_check.get('max_dr_arrival_km')} > {cross_check_tol_km}"
        )
    else:
        # SILVER-eligible: now decide SILVER vs BRONZE on downstream guards.
        if tier0.get("status") == "ok" and not tier0.get("tier0_all_admitted"):
            verdict = "BRONZE"
            reasons.append(
                f"Tier-0 NN rejected a leg (max_dv={tier0.get('tier0_max_dv_kms')} km/s)"
            )
        if literature.get("status") == "published":
            verdict = "BRONZE"
            reasons.append(f"literature_check published: {literature.get('citation')}")
        if p_fp > P_FP_SILVER_MAX:
            verdict = "BRONZE"
            reasons.append(f"ml_flagger_p_fp={p_fp:.3f} > {P_FP_SILVER_MAX}")
        if verdict == "SILVER":
            reasons.append("all guards passed (closure + cross-check + NN + lit + ml)")

    dedup = moon_cycler_signature_hash(
        primary=target.primary,
        sequence=tuple(candidate.sequence),
        vinf_per_encounter_kms=tuple(closure.vinf_per_encounter_kms),
    )

    return ScoredCandidate(
        candidate_id=f"{target.target_id}-{candidate.index:08d}",
        primary=target.primary,
        sequence=tuple(candidate.sequence),
        n_rev=n_rev,
        residual_kms=closure.residual_kms,
        vinf_per_encounter_kms=tuple(closure.vinf_per_encounter_kms),
        tof_days=tuple(closure.tof_days),
        max_vinf_kms=(
            max(closure.vinf_per_encounter_kms) if closure.vinf_per_encounter_kms else 0.0
        ),
        tier0=tier0,
        cross_check=cross_check,
        literature=literature,
        ml_flagger_p_fp=p_fp,
        verdict=verdict,
        verdict_reasons=reasons,
        dedup_signature=dedup,
    )


# ---------------------------------------------------------------------------
# Top-level campaign runner (enumerate -> close -> SILVER -> score)
# ---------------------------------------------------------------------------


@dataclass
class CampaignSummary:
    """Aggregate counts for one (primary, moons, seq_lengths, n_rev_grid) run."""

    primary: str
    moons: tuple[str, ...]
    seq_lengths: tuple[int, ...]
    n_rev_grid: tuple[int, ...]
    enumerated: int = 0
    evaluated: int = 0
    closed: int = 0
    failed_close: int = 0
    silver_pre_guards: int = 0
    near_miss_count: int = 0  # residual in [gate, near_miss_kms)
    rows_scored: int = 0
    verdict_counts: dict[str, int] = field(
        default_factory=lambda: {"SILVER": 0, "BRONZE": 0, "REJECT": 0}
    )
    top5_near_misses: list[dict[str, Any]] = field(default_factory=list)
    elapsed_s: float = 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "primary": self.primary,
            "moons": list(self.moons),
            "seq_lengths": list(self.seq_lengths),
            "n_rev_grid": list(self.n_rev_grid),
            "enumerated": self.enumerated,
            "evaluated": self.evaluated,
            "closed": self.closed,
            "failed_close": self.failed_close,
            "silver_pre_guards": self.silver_pre_guards,
            "near_miss_count": self.near_miss_count,
            "rows_scored": self.rows_scored,
            "verdict_counts": dict(self.verdict_counts),
            "top5_near_misses": list(self.top5_near_misses),
            "elapsed_s": self.elapsed_s,
        }


def run_prioritized_scan(
    *,
    primary: str,
    moons: Sequence[str],
    seq_lengths: Sequence[int],
    n_rev_grid: Sequence[int],
    out_path: Path,
    gate_residual_kms: float = 0.05,
    near_miss_kms: float = 1.0,
    phase_samples: int = 12,
    tof_resonance_grid: Sequence[float] = (0.5, 1.0, 1.5, 2.0),
    max_candidates: int | None = None,
    progress_every: int = 50,
    git_sha: str = "uncommitted",
) -> CampaignSummary:
    """Run the prioritized repeated-moon scan over one (primary, moons) system.

    Enumerates candidates with :class:`RepeatedMoonTarget`, closes each, and
    for SILVERs (residual < gate) runs the Tier-0 + DOP853 + literature + ML
    flagger chain. Every SILVER candidate writes one JSONL row to
    ``out_path``; BRONZE / REJECT rows are written too so the negative side
    of the comparison vs the unscored #264 daemon is auditable.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    target = RepeatedMoonTarget(
        primary=primary,
        moons=tuple(moons),
        seq_lengths=tuple(seq_lengths),
        n_rev_grid=tuple(n_rev_grid),
        n_phase_samples=phase_samples,
        tof_resonance_grid=tuple(tof_resonance_grid),
        git_sha=git_sha,
    )

    # Build scorers once; the flagger is trained at the start of the run so
    # `p_fp` is the retrained model (#256, commit 4ea5992).
    prioritizer = FiveTierPrioritizer()
    flagger = trained_flagger()

    summary = CampaignSummary(
        primary=primary,
        moons=tuple(moons),
        seq_lengths=tuple(seq_lengths),
        n_rev_grid=tuple(n_rev_grid),
    )

    t0 = time.time()
    with out_path.open("w", encoding="utf-8") as fh:
        # Emit a leading metadata row so the JSONL is self-describing.
        meta = {
            "_meta": True,
            "task": "#285 prioritized repeated-moon scan",
            "primary": primary,
            "moons": list(moons),
            "seq_lengths": list(seq_lengths),
            "n_rev_grid": list(n_rev_grid),
            "gate_residual_kms": gate_residual_kms,
            "phase_samples": phase_samples,
            "tof_resonance_grid": list(tof_resonance_grid),
            "git_sha": git_sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        for cand in target.enumerate_candidates():
            summary.enumerated += 1
            if max_candidates is not None and summary.evaluated >= max_candidates:
                break
            summary.evaluated += 1
            closure = target.close(cand)
            if not closure.converged:
                summary.failed_close += 1
                if summary.evaluated % progress_every == 0:
                    elapsed = time.time() - t0
                    print(
                        f"[285] {primary} evaluated={summary.evaluated} "
                        f"closed={summary.closed} failed={summary.failed_close} "
                        f"silver_pre={summary.silver_pre_guards} "
                        f"verdicts={summary.verdict_counts} t={elapsed:.0f}s",
                        flush=True,
                    )
                continue
            summary.closed += 1
            if not math.isfinite(closure.residual_kms):
                continue
            if closure.residual_kms >= gate_residual_kms:
                # Track near-misses for the writeup: the existence-prior bands
                # may surface "closed but above the gate" candidates that the
                # multi-arc / DSM continuation work could revisit.
                if closure.residual_kms < near_miss_kms:
                    summary.near_miss_count += 1
                    nm = {
                        "candidate_index": cand.index,
                        "sequence": list(cand.sequence),
                        "n_rev": list(cand.payload.get("n_rev", [])),
                        "residual_kms": closure.residual_kms,
                        "max_vinf_kms": (
                            max(closure.vinf_per_encounter_kms)
                            if closure.vinf_per_encounter_kms
                            else 0.0
                        ),
                        "vinf_per_encounter_kms": list(closure.vinf_per_encounter_kms),
                        "tof_days": list(closure.tof_days),
                    }
                    summary.top5_near_misses.append(nm)
                    summary.top5_near_misses.sort(key=lambda d: d["residual_kms"])
                    summary.top5_near_misses = summary.top5_near_misses[:5]
                if summary.evaluated % progress_every == 0:
                    elapsed = time.time() - t0
                    print(
                        f"[285] {primary} evaluated={summary.evaluated} "
                        f"closed={summary.closed} failed={summary.failed_close} "
                        f"silver_pre={summary.silver_pre_guards} "
                        f"verdicts={summary.verdict_counts} t={elapsed:.0f}s",
                        flush=True,
                    )
                continue
            summary.silver_pre_guards += 1
            scored = _score_one_silver(
                target=target,
                candidate=cand,
                closure=closure,
                prioritizer=prioritizer,
                flagger=flagger,
                phase_samples=phase_samples,
            )
            fh.write(json.dumps(scored.as_row()) + "\n")
            fh.flush()
            summary.rows_scored += 1
            summary.verdict_counts[scored.verdict] = (
                summary.verdict_counts.get(scored.verdict, 0) + 1
            )
            print(
                f"[285] {primary} {scored.candidate_id} seq={scored.sequence} "
                f"n_rev={scored.n_rev} res={scored.residual_kms:.4f} km/s "
                f"verdict={scored.verdict} reasons={scored.verdict_reasons}",
                flush=True,
            )

    summary.elapsed_s = time.time() - t0
    return summary


__all__ = [
    "CROSS_CHECK_TOL_KM",
    "P_FP_SILVER_MAX",
    "CampaignSummary",
    "ScoredCandidate",
    "dop853_cross_check_leg",
    "run_prioritized_scan",
    "trained_flagger",
]
