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

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.review_queue import ReviewQueueEntry, load_review_queue
from cyclerfinder.nbody.propagator import RestrictedNBody
from cyclerfinder.verify.gauntlet import VerdictTier

Vec3 = NDArray[np.float64]


def load_silver_candidates(path: Path | str) -> list[ReviewQueueEntry]:
    """Return the SILVER-tier entries from a review-queue file (read-only).

    Filters :func:`load_review_queue` to ``verdict_tier == VerdictTier.SILVER``.
    These are the held candidates the rung propagates; nothing here mutates the
    queue or promotes anything.
    """
    return [
        entry for entry in load_review_queue(path) if entry.verdict_tier == VerdictTier.SILVER.value
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


__all__ = ["RungArc", "load_silver_candidates", "propagate_one_period"]
