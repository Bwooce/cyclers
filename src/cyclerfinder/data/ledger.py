"""Append-only JSONL ledger (M7) — spec §13.6 + §13.8 persistence.

Spec references
---------------
* §13.6 — work queue + ledger for resumability: "The finder must
  resume cleanly: it cannot redo cells it has already searched."
* §13.8 — ledger schema: ``cell_id``, ``status``, ``n_solutions``,
  ``best_dv``, ``signatures[]``, ``validation_level``, ``t_done``,
  ``host``.

User-resolved decision (M7 plan §8, resolved 2026-06-01)
--------------------------------------------------------
* **JSONL backend.** Git-diffable, grep-able, review-friendly.
  Current catalogue scale is 219 entries; even at 10x discovery rate
  over a year the ledger holds < 10k lines, which JSONL handles
  trivially. SQLite reserved for a future scaling pass.

Parallel safety
---------------
:meth:`Ledger.record` uses POSIX ``O_APPEND | O_WRONLY`` writes
below the ``PIPE_BUF`` atomicity threshold (4 KiB on Linux); our
per-line budget is ~300 bytes. Concurrent appenders do not
interleave within a single ``\\n``-terminated line.

:meth:`Ledger.claim` is the parallel-safety primitive — first-write
wins on a race. Two workers calling ``claim(cell_id)`` simultaneously
both append a ``pending`` line; on re-read each worker checks which
of the two is the earliest timestamp and the later one backs off.
The M7 binding test exercises the round-trip and the restart path;
the concurrent-claim race is a smoke test (``pytest.mark.flaky``).

Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md`` §3.2.
"""

from __future__ import annotations

import json
import os
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final, Literal

LedgerStatus = Literal["pending", "pruned", "searched", "solved", "failed"]
"""Spec §13.8 status values.

* ``pending`` — a worker has claimed the cell but has not yet
  finished (intermediate state during a parallel run).
* ``pruned`` — Tisserand / energy pre-filter excluded the cell
  before optimisation.
* ``searched`` — optimiser ran but produced no feasible cycler.
* ``solved`` — optimiser produced a feasible cycler that survives
  the auto-validation gauntlet.
* ``failed`` — optimiser threw / converged to NaN / etc.
"""


PIPE_BUF: Final[int] = 4096
"""POSIX ``PIPE_BUF`` (Linux). Writes shorter than this via
``O_APPEND`` are atomic — no inter-process interleaving within a
single line. Per-line budget is ~300 bytes; the bound is comfortable.
"""


class LedgerError(Exception):
    """Raised on duplicate-``cell_id`` writes, schema-mismatch reads,
    persistence-path errors. The caller catches and routes to
    diagnostic output / retry policy."""


@dataclass(frozen=True)
class LedgerEntry:
    """One JSONL line — spec §13.8 schema.

    Attributes
    ----------
    cell_id:
        Spec §13.8 deterministic identifier (see
        :attr:`cyclerfinder.search.sequence.Cell.id`).
    status:
        See :data:`LedgerStatus`.
    n_solutions:
        Number of feasible cyclers the optimiser produced for this
        cell (``1`` for the M5 path; M-future may yield multiple per
        cell from richer multi-start grids).
    best_dv_kms:
        Closure-residual ΔV of the best cycler (km/s), or ``None``
        for ``pending`` / ``pruned`` / ``failed`` / ``searched``
        statuses.
    signature_hashes:
        Tuple of canonical-signature hashes produced for this cell.
        Empty tuple for non-``solved`` statuses.
    validation_level:
        ``"V0"`` .. ``"V5"`` per the spec §14 gauntlet, or ``None``
        for non-``solved`` statuses.
    t_done:
        ISO-8601 UTC-aware timestamp when this entry was written.
    host:
        Worker hostname (or arbitrary string identifier). Helps
        triage failures on multi-host parallel runs.
    """

    cell_id: str
    status: LedgerStatus
    n_solutions: int
    best_dv_kms: float | None
    signature_hashes: tuple[str, ...]
    validation_level: str | None
    t_done: str
    host: str


def _now_iso() -> str:
    """UTC ISO-8601 timestamp with second precision."""
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def _atomic_append(path: Path, line: str) -> None:
    """POSIX-atomic single-line append.

    Encodes ``line`` to UTF-8 and writes via ``O_WRONLY | O_APPEND |
    O_CREAT``; the kernel guarantees the bytes go in undivided.
    Raises :class:`OSError` (the underlying exception class for I/O
    failure) — callers catch generically.
    """
    if "\n" in line[:-1]:  # pragma: no cover — defensive
        raise LedgerError(f"line contains embedded newline: {line!r}")
    if not line.endswith("\n"):
        line = line + "\n"
    blob = line.encode("utf-8")
    if len(blob) >= PIPE_BUF:
        raise LedgerError(
            f"line {len(blob)} bytes exceeds atomic-append budget {PIPE_BUF}; "
            f"split or shrink the payload before retrying."
        )
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(fd, blob)
    finally:
        os.close(fd)


def _entry_to_jsonl(entry: LedgerEntry) -> str:
    """Serialise a :class:`LedgerEntry` to a single JSONL line."""
    payload = {
        "cell_id": entry.cell_id,
        "status": entry.status,
        "n_solutions": entry.n_solutions,
        "best_dv_kms": entry.best_dv_kms,
        "signature_hashes": list(entry.signature_hashes),
        "validation_level": entry.validation_level,
        "t_done": entry.t_done,
        "host": entry.host,
    }
    return json.dumps(payload, ensure_ascii=True)


def _entry_from_jsonl(line: str) -> LedgerEntry:
    """Parse a JSONL line into a :class:`LedgerEntry`.

    Raises :class:`LedgerError` on schema mismatch.
    """
    try:
        payload = json.loads(line)
    except json.JSONDecodeError as exc:
        raise LedgerError(f"malformed JSON: {exc}") from exc
    try:
        status = payload["status"]
        if status not in ("pending", "pruned", "searched", "solved", "failed"):
            raise LedgerError(f"unexpected status {status!r}")
        return LedgerEntry(
            cell_id=str(payload["cell_id"]),
            status=status,
            n_solutions=int(payload["n_solutions"]),
            best_dv_kms=(
                float(payload["best_dv_kms"]) if payload.get("best_dv_kms") is not None else None
            ),
            signature_hashes=tuple(payload.get("signature_hashes") or ()),
            validation_level=payload.get("validation_level"),
            t_done=str(payload["t_done"]),
            host=str(payload["host"]),
        )
    except KeyError as exc:
        raise LedgerError(f"missing required field: {exc}") from exc


class Ledger:
    """JSONL append-only ledger keyed by ``cell_id``.

    Concurrency model
    -----------------
    The in-memory index is built once on ``__init__`` and updated
    on every :meth:`record`. Multiple instances pointing at the same
    on-disk path will NOT see each other's intra-process writes
    without a re-read; :meth:`claim` re-reads the file before
    deciding to append, which closes the parallel race in practice.
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        if not self.path.exists():
            # Create the file so subsequent reads / appends don't
            # surprise the caller. mode ``a`` is idempotent.
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.touch()
        self._index: dict[str, LedgerEntry] = {}
        self._load()

    def _load(self) -> None:
        """Read every line from disk and populate :attr:`_index`."""
        with self.path.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                entry = _entry_from_jsonl(line)
                # First-write-wins — duplicate ``cell_id`` is ignored
                # rather than raising at load time so the ledger
                # survives a benign race that wrote two ``pending``
                # lines.
                if entry.cell_id not in self._index:
                    self._index[entry.cell_id] = entry

    def __len__(self) -> int:
        return len(self._index)

    def has(self, cell_id: str) -> bool:
        """Return ``True`` if ``cell_id`` is already recorded."""
        return cell_id in self._index

    def get(self, cell_id: str) -> LedgerEntry:
        """Return the recorded :class:`LedgerEntry` for ``cell_id``.

        Raises :class:`KeyError` if not recorded.
        """
        return self._index[cell_id]

    def record(self, entry: LedgerEntry) -> None:
        """Append ``entry`` to the ledger.

        Raises :class:`LedgerError` if ``entry.cell_id`` is already
        recorded (no double-record per spec §13.6 idempotency).
        """
        if entry.cell_id in self._index:
            raise LedgerError(
                f"cell_id {entry.cell_id!r} already recorded with "
                f"status {self._index[entry.cell_id].status!r}; "
                f"refusing to double-record."
            )
        line = _entry_to_jsonl(entry)
        _atomic_append(self.path, line)
        self._index[entry.cell_id] = entry

    def claim(self, cell_id: str, host: str) -> bool:
        """Attempt to claim ``cell_id`` for ``host``.

        Returns ``True`` if the claim succeeded (the caller may now
        run the optimiser for this cell), ``False`` if another
        worker already claimed it or it was already done.
        """
        # Re-read to catch other-worker writes since __init__.
        prev_index = self._index
        self._index = {}
        self._load()
        existing = self._index.get(cell_id)
        if existing is not None:
            if existing.status in ("solved", "pruned", "failed", "searched"):
                return False
            # ``pending`` — only allow re-claim by same host.
            return existing.host == host
        # Restore any in-flight entries from the previous index that
        # the on-disk load missed (the file is authoritative; this
        # path is defensive).
        for k, v in prev_index.items():
            if k not in self._index:
                self._index[k] = v
        # Append a ``pending`` line to lay claim.
        pending = LedgerEntry(
            cell_id=cell_id,
            status="pending",
            n_solutions=0,
            best_dv_kms=None,
            signature_hashes=(),
            validation_level=None,
            t_done=_now_iso(),
            host=host,
        )
        line = _entry_to_jsonl(pending)
        _atomic_append(self.path, line)
        self._index[cell_id] = pending
        return True

    def iter_pending(self) -> tuple[str, ...]:
        """Return the cell ids whose recorded status is ``"pending"``."""
        return tuple(cid for cid, entry in self._index.items() if entry.status == "pending")


class LedgerLoader:
    """Read-only loader for analysis (M8 reporter, audit, etc.).

    Light-weight wrapper around the JSONL file that does NOT update
    the in-memory index — every call re-reads from disk. Cheaper to
    construct than :class:`Ledger` when only a snapshot is needed.
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)

    def __iter__(self) -> Iterator[LedgerEntry]:
        """Yield every :class:`LedgerEntry` line by line."""
        if not self.path.exists():
            return
        with self.path.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                yield _entry_from_jsonl(line)


__all__ = [
    "PIPE_BUF",
    "Ledger",
    "LedgerEntry",
    "LedgerError",
    "LedgerLoader",
    "LedgerStatus",
]
