"""SILVER validation rungs — n-body propagation of held candidates (design §4).

Consumer 1 of the harness. Ingests the USER-HELD SILVER candidates from the
review queue (``data/review_queue.py``), propagates one full period in the rails
n-body (Sun + E + M + J, design §2), measures the terminal closure degradation
and the node-impulse correction ΔV (:mod:`cyclerfinder.nbody.correction_dv`), and
**records the verdict to the review-queue audit trail — never auto-promotes**
(``review_queue.is_catalogue_source()`` is ``False`` by contract; the human
decides).

The rung runs an *independent integrator* (REBOUND) over the *same* DE440 BSP
astropy caches: a cross-check rung, NOT a V4 stamp (design §4, Q6). The candidate
numerics are OUR computation (SILVER = unsourced by tier definition), so the rung
asserts a closure / correction-ΔV *regime* (ROBUST / MARGINAL / ARTIFACT), never a
sourced value (golden discipline).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.review_queue import ReviewQueueEntry, load_review_queue
from cyclerfinder.nbody.correction_dv import CorrectionDV, node_impulse_correction_dv
from cyclerfinder.nbody.propagator import RestrictedNBody
from cyclerfinder.verify.gauntlet import VerdictTier

Vec3 = NDArray[np.float64]

# Human-declared rung verdict thresholds (plan Phase B, the Jones <200 m/s
# analogue). RECORDED, not gating: the rung writes the tier; the human decides.
_ROBUST_MAX_KMS = 0.200  # < 200 m/s total correction ΔV
_MARGINAL_MAX_KMS = 1.000  # 200-1000 m/s


def load_silver_candidates(path: Path | str) -> list[ReviewQueueEntry]:
    """Return the SILVER-tier entries from a review-queue file (read-only).

    Filters :func:`load_review_queue` to ``VerdictTier.SILVER`` via the
    enum-round-trip form ``VerdictTier(entry.verdict_tier) is VerdictTier.SILVER``
    rather than a raw string compare against ``.value``: this makes the comparison
    enum-typed (a renamed/aliased value tracks the enum, not a stale literal) and a
    **malformed** ``verdict_tier`` raises ``ValueError`` from the enum constructor
    rather than being silently dropped as "not SILVER". These are the held
    candidates the rung propagates; nothing here mutates the queue or promotes
    anything.
    """
    return [
        entry
        for entry in load_review_queue(path)
        if VerdictTier(entry.verdict_tier) is VerdictTier.SILVER
    ]


@dataclass(frozen=True)
class RungArc:
    """Frozen result of a single-period rung propagation (design §5 discipline)."""

    candidate_id: str
    terminal_closure_km: float
    terminal_closure_kms: float
    period_sec: float
    bodies: tuple[str, ...]
    integrator_accuracy: float
    converged: bool
    seed_node: dict[str, Vec3] = field(default_factory=dict)
    wrap_node: dict[str, Vec3] = field(default_factory=dict)


def _seed_home_earth_state(
    entry: ReviewQueueEntry, ephem: Ephemeris, t0_sec: float
) -> tuple[Vec3, Vec3]:
    """Reconstruct the spacecraft state at the home-Earth (b0) node.

    ``v_sc = v_planet + vinf_out`` (the ``verify/propagate.py:434`` reconstruction
    the design §3 cites). The recorded candidate carries only V∞ *magnitudes*
    (the seed *vectors* were not on disk; see the rung fixtures), so the outgoing
    V∞ is reconstructed as the recorded home-Earth magnitude applied along the
    Earth velocity direction (a prograde representative). NON-GOLDEN: the rung
    measures the *closure regime* of this reconstructed seed, never a value.
    """
    r_e, v_e = ephem.state(entry.sequence[0], t0_sec)
    r_e = np.asarray(r_e, dtype=np.float64)
    v_e = np.asarray(v_e, dtype=np.float64)
    vinf_mag = float(entry.vinf_per_encounter_kms[0])
    v_hat = v_e / float(np.linalg.norm(v_e))
    v_sc = v_e + vinf_mag * v_hat
    return r_e, v_sc


def propagate_one_period(
    entry: ReviewQueueEntry,
    ephem: Ephemeris,
    *,
    bodies: tuple[str, ...] = ("E", "M", "J"),
    accuracy: float = 1e-10,
    t0_sec: float = 0.0,
) -> RungArc:
    """Propagate the candidate one full period in the rails n-body (design §2).

    Seeds the spacecraft at the home-Earth node, propagates one repeat period
    (``sum(tof_days)``) in Sun + E + M + J, and measures the **terminal closure
    error**: the heliocentric gap between the propagated wrap node and the seeded
    start node. The candidate is unsourced (SILVER), so the closure error is
    *finite and recorded*, not asserted against a value.
    """
    period_days = float(sum(entry.tof_days))
    period_sec = period_days * SECONDS_PER_DAY

    r0, v0 = _seed_home_earth_state(entry, ephem, t0_sec)

    prop = RestrictedNBody("rebound")
    arc = prop.propagate(
        r0,
        v0,
        t0_sec=t0_sec,
        t1_sec=t0_sec + period_sec,
        bodies=bodies,
        accuracy=accuracy,
        ephem=ephem,
    )

    # A non-finite wrap state means the integration blew up (e.g. a sub-safe-
    # periapsis dive the softened force could not tame): an honest DIVERGENT rung
    # outcome, surfaced as converged=False with a large finite sentinel closure —
    # never a NaN/exception (mirror correct.py's honest non-converged record).
    finite = bool(np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s)))
    if finite:
        closure_km = float(np.linalg.norm(arc.r_km - r0))
        closure_kms = float(np.linalg.norm(arc.v_km_s - v0))
        converged = arc.converged
    else:
        closure_km = float("inf")
        closure_kms = float("inf")
        converged = False

    return RungArc(
        candidate_id=entry.candidate_id,
        terminal_closure_km=closure_km,
        terminal_closure_kms=closure_kms,
        period_sec=period_sec,
        bodies=bodies,
        integrator_accuracy=accuracy,
        converged=converged,
        seed_node={"r_km": r0, "v_km_s": v0},
        wrap_node={"r_km": arc.r_km, "v_km_s": arc.v_km_s},
    )


def _closure_correction_dv(arc: RungArc) -> CorrectionDV:
    """Node-impulse correction ΔV implied by the terminal closure (Phase B).

    Phase B runs a single forward propagation, not the full multiple-shooter
    (that is Phase C). The minimal node-impulse to restore periodicity from this
    seed is the velocity discontinuity at the wrap node — the burn that maps the
    propagated wrap velocity back onto the seeded start velocity. Expressed
    through :func:`node_impulse_correction_dv` so the convention (design §3, Q3) is
    identical to the shooter's. A divergent arc carries a finite (sentinel) gap
    here; the verdict path grades not-converged as ARTIFACT regardless.
    """
    seed = {"wrap": np.asarray(arc.wrap_node["v_km_s"], dtype=np.float64)}
    corr = {"wrap": np.asarray(arc.seed_node["v_km_s"], dtype=np.float64)}
    return node_impulse_correction_dv(seed, corr)


@dataclass(frozen=True)
class RungVerdict:
    """Frozen rung verdict (the recorded tier + the numbers behind it)."""

    tier: str  # "ROBUST" | "MARGINAL" | "ARTIFACT"
    correction_dv_kms: float
    terminal_closure_km: float
    converged: bool


def rung_verdict(
    correction_dv_kms: float,
    *,
    terminal_closure_km: float,
    converged: bool,
) -> RungVerdict:
    """Apply the human-declared thresholds (plan Phase B).

    ROBUST < 200 m/s; MARGINAL 200-1000 m/s; ARTIFACT >= 1000 m/s OR diverged.
    These are RECORDED, not gating — the human makes the promotion call.
    """
    if not converged or correction_dv_kms >= _MARGINAL_MAX_KMS:
        tier = "ARTIFACT"
    elif correction_dv_kms < _ROBUST_MAX_KMS:
        tier = "ROBUST"
    else:
        tier = "MARGINAL"
    return RungVerdict(
        tier=tier,
        correction_dv_kms=float(correction_dv_kms),
        terminal_closure_km=float(terminal_closure_km),
        converged=bool(converged),
    )


def record_rung_result(
    entry: ReviewQueueEntry,
    verdict: RungVerdict,
    path: Path | str,
    *,
    bodies: tuple[str, ...] = ("E", "M", "J"),
) -> dict[str, object]:
    """Append a NON-PROMOTING rung-audit record to the audit trail (design §4, Q6).

    Writes one JSON line capturing the rung verdict + numbers for the candidate.
    It MUST NOT change ``verdict_tier``, MUST NOT promote, and MUST NOT write a
    catalogue row (``review_queue.is_catalogue_source()`` is ``False``). The
    record carries ``promoted=False`` explicitly; the human decides. Returns the
    record dict that was written.
    """
    record: dict[str, object] = {
        "candidate_id": entry.candidate_id,
        "signature_hash": entry.signature_hash,
        "queue_tier": entry.verdict_tier,  # unchanged; recorded for traceability
        "rung_verdict": verdict.tier,
        "correction_dv_kms": verdict.correction_dv_kms,
        "terminal_closure_km": verdict.terminal_closure_km,
        "converged": verdict.converged,
        "bodies": list(bodies),
        "promoted": False,  # golden discipline: the rung never promotes
        "rung_kind": "nbody-silver-rung",
        "independence": "shared-DE440 cross-check (NOT a V4 stamp; design Q6)",
        "t_recorded": datetime.now(UTC).isoformat(),
    }
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, ensure_ascii=True) + "\n")
    return record


def run_rung(
    entry: ReviewQueueEntry,
    ephem: Ephemeris,
    *,
    audit_path: Path | str,
    bodies: tuple[str, ...] = ("E", "M", "J"),
    accuracy: float = 1e-9,
) -> dict[str, object]:
    """Run the full SILVER rung end-to-end and RECORD it (never promote).

    Propagate one period (design §2 Sun+E+M+J), compute the node-impulse
    correction ΔV from the terminal closure, grade the verdict, and append a
    non-promoting audit record. Returns the record dict.
    """
    arc = propagate_one_period(entry, ephem, bodies=bodies, accuracy=accuracy)
    dv = _closure_correction_dv(arc)
    verdict = rung_verdict(
        dv.total_kms,
        terminal_closure_km=arc.terminal_closure_km,
        converged=arc.converged,
    )
    return record_rung_result(entry, verdict, audit_path, bodies=bodies)


def jupiter_sensitivity(
    entry: ReviewQueueEntry,
    ephem: Ephemeris,
    *,
    accuracy: float = 1e-9,
) -> dict[str, object]:
    """Gate-4 body-inclusion arm: rerun with/without Jupiter (design §2, §5.3).

    Records the correction ΔV with the full Sun+E+M+J set and with Sun+E+M only,
    and their difference, so Jupiter's contribution at this candidate's baseline
    is *evidence*, not assertion (the standing §2 rule). Returns both verdicts +
    the delta; the verdict tier is unchanged (recorded, not gating).
    """
    arc_with = propagate_one_period(entry, ephem, bodies=("E", "M", "J"), accuracy=accuracy)
    arc_without = propagate_one_period(entry, ephem, bodies=("E", "M"), accuracy=accuracy)
    dv_with = _closure_correction_dv(arc_with).total_kms
    dv_without = _closure_correction_dv(arc_without).total_kms
    return {
        "with_jupiter": {
            "correction_dv_kms": dv_with,
            "terminal_closure_km": arc_with.terminal_closure_km,
            "converged": arc_with.converged,
        },
        "without_jupiter": {
            "correction_dv_kms": dv_without,
            "terminal_closure_km": arc_without.terminal_closure_km,
            "converged": arc_without.converged,
        },
        "delta_correction_dv_kms": dv_with - dv_without,
    }


__all__ = [
    "CorrectionDV",
    "RungArc",
    "RungVerdict",
    "jupiter_sensitivity",
    "load_silver_candidates",
    "propagate_one_period",
    "record_rung_result",
    "run_rung",
    "rung_verdict",
]
