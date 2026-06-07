"""Phase 5 — the human-review queue (the Forge novelty gate).

SILVER candidates (machine-confirmed, unsourced) land here, each with a full
audit trail and the adversarial-panel result, **pending human approval**. This
artifact is the golden-discipline boundary: a novel candidate is *never*
auto-promoted to a catalogue row. :func:`is_catalogue_source` returns ``False``
so any loader that walks data sources treats the queue as non-catalogue — a
human must hand-promote an approved entry into ``data/catalogue.yaml``.

Default on-disk location: ``data/review_queue.jsonl`` (JSONL, one entry per
line, git-diffable + grep-able, matching the ledger convention).
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.verify.gauntlet import VerdictTier
from cyclerfinder.verify.plausibility import QuantityKind, check_publishable

# Tiers that may legitimately sit in the human-review queue. SILVER is the novel
# holding tier; GOLD is admissible only for a sourced rediscovery (the loop never
# fabricates a source, so GOLD here would come from a caller-supplied match).
# Compared via the enum (round-trip ``VerdictTier(value) in _QUEUEABLE_TIERS``)
# rather than raw strings, so a malformed tier string raises ValueError on
# coercion instead of silently failing an ``in`` test against bare strings.
_QUEUEABLE_TIERS: frozenset[VerdictTier] = frozenset({VerdictTier.SILVER, VerdictTier.GOLD})

DEFAULT_REVIEW_QUEUE_PATH: Path = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "review_queue.jsonl"
)


@dataclass(frozen=True)
class ReviewQueueEntry:
    """One queued candidate awaiting human review (frozen, audit-complete).

    Carries enough to reconstruct the finding *and* the adversarial-panel
    outcome, so a human reviewer (or a re-validation pass) sees the full chain
    without re-running the search.
    """

    candidate_id: str
    signature_hash: str
    verdict_tier: str
    match_outcome: str
    known_id: str | None
    superseded_by: tuple[str, ...]
    vinf_per_encounter_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    bend_feasible: bool
    max_vinf_kms: float
    sequence: tuple[str, ...]
    period_k: int
    model_assumption: str
    verdict_audit: dict[str, Any]
    panel: dict[str, Any]
    t_added: str


def is_catalogue_source() -> bool:
    """The review queue is NON-catalogue (golden discipline).

    Returns ``False`` always. A loader that enumerates catalogue sources must
    skip the queue: a queued candidate becomes a catalogue row only by explicit
    human action, never automatically.
    """
    return False


def validate_review_entry(entry: ReviewQueueEntry) -> None:
    """Refuse an entry that does not belong in the human-review queue.

    Raises :class:`ValueError` if:

    * the tier is not queueable (only SILVER / GOLD may await human review — a
      REJECTED or BRONZE candidate is never promoted), or
    * the adversarial panel majority-refuted the candidate (the falsification
      gate has teeth: a refuted candidate must not be queued), or
    * any recorded per-encounter V_inf breaks the elliptic-periodicity physics
      ceiling (task #127): a periodic heliocentric cycler physically cannot have
      V_inf above ``v_esc_sun(r_B) + v_B`` at a body, so such an entry is a
      degenerate / off-family artifact and must not even reach the human queue.
    """
    try:
        tier = VerdictTier(entry.verdict_tier)
    except ValueError as exc:
        raise ValueError(
            f"review-queue entry {entry.candidate_id!r} has unknown verdict tier "
            f"{entry.verdict_tier!r}: not a {VerdictTier.__name__} value."
        ) from exc
    if tier not in _QUEUEABLE_TIERS:
        queueable = sorted(t.value for t in _QUEUEABLE_TIERS)
        raise ValueError(
            f"review-queue entry {entry.candidate_id!r} has tier "
            f"{entry.verdict_tier!r}; only {queueable} may be queued "
            f"for human review (golden discipline: no auto-promotion)."
        )
    panel = entry.panel or {}
    if panel.get("majority_refute"):
        raise ValueError(
            f"review-queue entry {entry.candidate_id!r} was majority-refuted by the "
            f"adversarial panel ({panel.get('n_refuted')}/{panel.get('n_verifiers')}); "
            f"refusing to queue a refuted candidate."
        )
    # Physics-ceiling guard (task #127): each recorded encounter V_inf must be
    # publishable against its body's elliptic-periodicity ceiling.
    for body, vinf in zip(entry.sequence, entry.vinf_per_encounter_kms, strict=False):
        verdict = check_publishable(QuantityKind.VINF_KMS, vinf, {"body": body})
        if not verdict.ok:
            raise ValueError(
                f"review-queue entry {entry.candidate_id!r} records an implausible "
                f"V_inf at {body}: {verdict.reason}"
            )


def _normalise(payload: dict[str, Any]) -> dict[str, Any]:
    """Coerce JSON lists back into the tuple-typed fields for round-trip equality."""
    for key in ("superseded_by", "vinf_per_encounter_kms", "tof_days", "sequence"):
        if key in payload and isinstance(payload[key], list):
            payload[key] = tuple(payload[key])
    return payload


def append_review_entry(
    path: Path | str,
    entry: ReviewQueueEntry,
    *,
    validate: bool = True,
) -> None:
    """Append one validated entry to the review queue (creating parent dirs)."""
    if validate:
        validate_review_entry(entry)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(asdict(entry), ensure_ascii=True, default=list)
    with p.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_review_queue(path: Path | str) -> Iterator[ReviewQueueEntry]:
    """Yield every :class:`ReviewQueueEntry` from the queue file (read-only)."""
    p = Path(path)
    if not p.exists():
        return
    for raw in p.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        payload = _normalise(json.loads(line))
        yield ReviewQueueEntry(**payload)


def _now_iso() -> str:
    return datetime.now(UTC).isoformat()


__all__ = [
    "DEFAULT_REVIEW_QUEUE_PATH",
    "ReviewQueueEntry",
    "append_review_entry",
    "is_catalogue_source",
    "load_review_queue",
    "validate_review_entry",
]
