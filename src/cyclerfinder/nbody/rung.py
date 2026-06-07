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

from pathlib import Path

from cyclerfinder.data.review_queue import ReviewQueueEntry, load_review_queue
from cyclerfinder.verify.gauntlet import VerdictTier


def load_silver_candidates(path: Path | str) -> list[ReviewQueueEntry]:
    """Return the SILVER-tier entries from a review-queue file (read-only).

    Filters :func:`load_review_queue` to ``verdict_tier == VerdictTier.SILVER``.
    These are the held candidates the rung propagates; nothing here mutates the
    queue or promotes anything.
    """
    return [
        entry for entry in load_review_queue(path) if entry.verdict_tier == VerdictTier.SILVER.value
    ]


__all__ = ["load_silver_candidates"]
