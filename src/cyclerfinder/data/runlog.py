"""Run artefact persistence — (seed -> basin) labels as committed evidence.

Origin: ``docs/notes/2026-06-07-ml-surrogate-investigation.md`` Build B
precondition. Campaigns and MBH runs already COMPUTE per-row (seed -> basin)
outcome labels (``CLOSE-AND-MATCH`` / ``CLOSE-OFF-ANCHOR`` / ``NO-CLOSE`` /
``CLOSE-MATCH-SYMMETRIC-ONLY``) with achieved-vs-sourced V∞ per body, but print
them to stdout and discard. This module is a small, schema'd JSON-lines writer /
reader that persists those labels so a corpus can accumulate. It has standalone
audit value (a versioned, grep-able trail of every campaign closure) regardless
of any future classifier.

Design (mirrors :mod:`cyclerfinder.data.ledger`)
------------------------------------------------
* **JSONL backend**, one file per campaign invocation under ``data/runs/``
  (e.g. ``data/runs/russell12-<timestamp>.jsonl``). Git-diffable, grep-able,
  review-friendly — consistent with this repo's audit culture and the committed
  ``data/gauntlet_ledger.jsonl``.
* **Frozen, typed records.** :class:`RunRecord` is the schema; the loader returns
  typed records and rejects lines missing a required field
  (:class:`RunlogError`).
* **Append semantics.** :meth:`RunLog.append` adds one record; re-opening a path
  and appending more records is supported (the file accumulates).

Why NOT the ledger's atomic-append path: ledger lines are ~300 bytes and bounded
below ``PIPE_BUF`` for lock-free parallel writes. Run records carry V∞ lists,
per-check strings and a solver-audit object, so a single line can exceed 4 KiB.
The campaign writer is single-process (one file per invocation), so a plain
buffered append is correct and simpler; we do not claim multi-writer atomicity.

Commit policy (DECISION)
------------------------
``data/runs/`` IS COMMITTED. The records are small, they are evidence artefacts,
and the repo already commits ``data/gauntlet_ledger.jsonl`` — there is no
existing convention saying run outputs are gitignored (``.gitignore`` only
excludes ``out/`` per spec §7, which is a different, throwaway sink). Committing
the runlog is what lets the (seed -> basin) corpus accumulate across runs, which
is the entire point of the Build B precondition. See
``docs/notes/2026-06-07-runlog-persistence.md``.

The ``code_version`` (``git rev-parse --short HEAD``) is captured BY THE CALLER
and passed in — this module never shells out (keeps it pure / testable and avoids
a subprocess dependency in the data layer).
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class RunlogError(Exception):
    """Raised on malformed JSON or a record missing a required field."""


def _now_iso() -> str:
    """UTC ISO-8601 timestamp with second precision."""
    return datetime.now(UTC).replace(microsecond=0).isoformat()


@dataclass(frozen=True)
class RunRecord:
    """One JSON-lines run artefact: a single (seed -> basin) label.

    Required fields (a line missing any of these is rejected on read):

    row_id:
        Catalogue id of the row solved (e.g. ``russell-ch4-4.991gG2``).
    genome:
        Genome path used (``lambert`` / ``free-return``) — the seed
        representation that produced this label.
    outcome:
        The basin label the campaign / MBH run computed:
        ``CLOSE-AND-MATCH`` / ``CLOSE-OFF-ANCHOR`` / ``NO-CLOSE`` /
        ``CLOSE-MATCH-SYMMETRIC-ONLY`` / ``NO-SEED``.
    model:
        Ephemeris model (``circular`` / ``astropy`` / ``inclined-circular``).
    code_version:
        ``git rev-parse --short HEAD`` captured by the CALLER (this module
        never shells out). The empty string is accepted (e.g. running outside
        a git checkout) but should be filled by campaign callers.

    Optional (additive) fields — default to empty / ``None`` so older lines and
    leaner callers (e.g. MBH hop logs) round-trip cleanly:

    achieved_vinf_kms:
        ``{body: V∞}`` that EMERGED from the converged solution.
    sourced_vinf_kms:
        ``{body: V∞}`` SOURCED anchors the outcome was judged against.
    sourced_anchors:
        Other sourced anchors used in the match (period, ratios, transits).
    seed:
        Seed parameters (the genome's free variables / ToF or (a, e) seed).
    residual_kms:
        Best closure residual (km/s).
    solver_audit:
        Solver bookkeeping — ``nfev`` / ``hops`` (attempted/accepted) / ``rng_seed``
        / per-mode close counts — whatever the call site carries.
    t_written:
        ISO-8601 UTC timestamp; defaults to write time.
    """

    row_id: str
    genome: str
    outcome: str
    model: str
    code_version: str
    achieved_vinf_kms: Mapping[str, float] = field(default_factory=dict)
    sourced_vinf_kms: Mapping[str, float] = field(default_factory=dict)
    sourced_anchors: Mapping[str, Any] = field(default_factory=dict)
    seed: Mapping[str, Any] = field(default_factory=dict)
    residual_kms: float | None = None
    solver_audit: Mapping[str, Any] = field(default_factory=dict)
    t_written: str = ""


_REQUIRED_FIELDS = ("row_id", "genome", "outcome", "model", "code_version")


def _record_to_jsonl(rec: RunRecord) -> str:
    """Serialise a :class:`RunRecord` to one JSON-lines line (no newline)."""
    payload: dict[str, Any] = {
        "row_id": rec.row_id,
        "genome": rec.genome,
        "outcome": rec.outcome,
        "model": rec.model,
        "code_version": rec.code_version,
        "achieved_vinf_kms": dict(rec.achieved_vinf_kms),
        "sourced_vinf_kms": dict(rec.sourced_vinf_kms),
        "sourced_anchors": dict(rec.sourced_anchors),
        "seed": dict(rec.seed),
        "residual_kms": rec.residual_kms,
        "solver_audit": dict(rec.solver_audit),
        "t_written": rec.t_written or _now_iso(),
    }
    return json.dumps(payload, ensure_ascii=True)


def _record_from_jsonl(line: str) -> RunRecord:
    """Parse one JSON-lines line into a :class:`RunRecord`.

    Raises :class:`RunlogError` on malformed JSON or a missing required field.
    """
    try:
        payload = json.loads(line)
    except json.JSONDecodeError as exc:
        raise RunlogError(f"malformed JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise RunlogError(f"expected a JSON object, got {type(payload).__name__}")
    for required in _REQUIRED_FIELDS:
        if required not in payload:
            raise RunlogError(f"missing required field: {required!r}")
    residual = payload.get("residual_kms")
    return RunRecord(
        row_id=str(payload["row_id"]),
        genome=str(payload["genome"]),
        outcome=str(payload["outcome"]),
        model=str(payload["model"]),
        code_version=str(payload["code_version"]),
        achieved_vinf_kms=dict(payload.get("achieved_vinf_kms") or {}),
        sourced_vinf_kms=dict(payload.get("sourced_vinf_kms") or {}),
        sourced_anchors=dict(payload.get("sourced_anchors") or {}),
        seed=dict(payload.get("seed") or {}),
        residual_kms=(float(residual) if residual is not None else None),
        solver_audit=dict(payload.get("solver_audit") or {}),
        t_written=str(payload.get("t_written") or ""),
    )


class RunLog:
    """Append-only JSON-lines run artefact, one file per campaign invocation.

    The file is created on construction if absent (so a fresh ``--runlog-dir``
    works without a separate mkdir). Records are appended one per line;
    :meth:`read` re-parses the whole file into typed records.
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        if not self.path.exists():
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.touch()

    def append(self, record: RunRecord) -> None:
        """Append one :class:`RunRecord` as a JSON-lines line."""
        line = _record_to_jsonl(record)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    def extend(self, records: Sequence[RunRecord]) -> None:
        """Append several records in one open — convenience for a campaign run."""
        with self.path.open("a", encoding="utf-8") as fh:
            for record in records:
                fh.write(_record_to_jsonl(record) + "\n")

    def read(self) -> list[RunRecord]:
        """Parse every line into a typed :class:`RunRecord` (in file order)."""
        return list(self.iter_records())

    def iter_records(self) -> Iterator[RunRecord]:
        """Yield each :class:`RunRecord`; raises :class:`RunlogError` on a bad line."""
        with self.path.open(encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                yield _record_from_jsonl(line)

    def __len__(self) -> int:
        return sum(1 for _ in self.iter_records())


def default_runlog_path(runlog_dir: Path | str, tag: str, timestamp: str | None = None) -> Path:
    """Compose ``<runlog_dir>/<tag>-<timestamp>.jsonl``.

    ``timestamp`` defaults to a filesystem-safe UTC stamp. Callers that want a
    reproducible / passed-in stamp (the task's ``<timestamp-arg>``) supply it.
    """
    stamp = timestamp or datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return Path(runlog_dir) / f"{tag}-{stamp}.jsonl"
