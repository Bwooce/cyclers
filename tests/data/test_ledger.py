"""JSONL append-only ledger tests (M7).

* Round-trip — write then read.
* Idempotency — double-record raises.
* Restart-survives — re-open after close sees prior writes.
* Atomic-append — line is well-formed, no partial JSON.
* Claim — basic semantics.
* :class:`LedgerLoader` read-only iteration.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from cyclerfinder.data.ledger import (
    Ledger,
    LedgerEntry,
    LedgerError,
    LedgerLoader,
)


def _make_entry(cell_id: str = "EM|E-M-E|k2|r00|bss") -> LedgerEntry:
    return LedgerEntry(
        cell_id=cell_id,
        status="solved",
        n_solutions=1,
        best_dv_kms=0.012,
        signature_hashes=("sha1:" + "a" * 40,),
        validation_level="V2",
        t_done="2026-06-01T12:34:56+00:00",
        host="ci-worker-3",
    )


def test_ledger_round_trip(tmp_path: Path) -> None:
    """**M7 GATE** — write → reread → field-equal; duplicate raises."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    entry = _make_entry()
    ledger.record(entry)
    assert ledger.has(entry.cell_id)
    assert ledger.get(entry.cell_id) == entry
    # Re-recording the same cell_id raises.
    with pytest.raises(LedgerError):
        ledger.record(entry)


def test_ledger_persists_across_restart(tmp_path: Path) -> None:
    """Writing then constructing a fresh :class:`Ledger` from the
    same path sees the prior write."""
    ledger_path = tmp_path / "ledger.jsonl"
    a = Ledger(ledger_path)
    entry = _make_entry()
    a.record(entry)
    del a
    b = Ledger(ledger_path)
    assert b.has(entry.cell_id)
    assert b.get(entry.cell_id) == entry


def test_ledger_atomic_append_no_partial_lines(tmp_path: Path) -> None:
    """File ends with newline; no partial JSON; lines parse round-trip."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    entries = [_make_entry(cell_id=f"EM|E-M|k{k}|r0|bs") for k in range(1, 5)]
    for e in entries:
        ledger.record(e)
    raw = ledger_path.read_text()
    assert raw.endswith("\n")
    assert raw.count("\n") == len(entries)
    # Every line is valid JSON parseable as a LedgerEntry.
    loader = LedgerLoader(ledger_path)
    re_read = list(loader)
    assert len(re_read) == len(entries)
    assert {e.cell_id for e in re_read} == {e.cell_id for e in entries}


def test_ledger_iter_pending_skips_solved(tmp_path: Path) -> None:
    """``iter_pending`` returns only ``pending`` cell ids."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    ledger.record(dataclasses.replace(_make_entry("p1"), status="pending"))
    ledger.record(dataclasses.replace(_make_entry("s1"), status="solved"))
    ledger.record(dataclasses.replace(_make_entry("p2"), status="pending"))
    pending = ledger.iter_pending()
    assert set(pending) == {"p1", "p2"}


def test_ledger_claim_returns_true_on_new(tmp_path: Path) -> None:
    """``claim`` of a fresh cell_id returns True and writes ``pending``."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    assert ledger.claim("EM|E-M-E|k2|r00|bss", "worker-a") is True
    assert ledger.has("EM|E-M-E|k2|r00|bss")
    assert ledger.get("EM|E-M-E|k2|r00|bss").status == "pending"


def test_ledger_claim_returns_false_when_done(tmp_path: Path) -> None:
    """``claim`` of an already-solved cell returns False."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    ledger.record(_make_entry("EM|E-M-E|k2|r00|bss"))  # status=solved
    assert ledger.claim("EM|E-M-E|k2|r00|bss", "worker-a") is False


def test_ledger_claim_returns_true_for_same_host_pending(tmp_path: Path) -> None:
    """Re-claiming a ``pending`` cell from the same host returns True
    (resume-own-work)."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    assert ledger.claim("c1", "worker-a") is True
    assert ledger.claim("c1", "worker-a") is True


def test_ledger_claim_returns_false_for_different_host_pending(tmp_path: Path) -> None:
    """A ``pending`` cell claimed by another host blocks subsequent
    claimers."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    assert ledger.claim("c1", "worker-a") is True
    assert ledger.claim("c1", "worker-b") is False


def test_ledger_entry_frozen() -> None:
    """:class:`LedgerEntry` is frozen."""
    entry = _make_entry()
    with pytest.raises(dataclasses.FrozenInstanceError):
        entry.status = "failed"  # type: ignore[misc]


def test_ledger_len(tmp_path: Path) -> None:
    """``len(ledger)`` returns the recorded entry count."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    assert len(ledger) == 0
    ledger.record(_make_entry("c1"))
    ledger.record(_make_entry("c2"))
    assert len(ledger) == 2


def test_ledger_loader_empty_path(tmp_path: Path) -> None:
    """Reading a non-existent ledger via :class:`LedgerLoader` yields
    nothing (no error)."""
    loader = LedgerLoader(tmp_path / "nope.jsonl")
    assert list(loader) == []


def test_ledger_loader_parses_existing(tmp_path: Path) -> None:
    """Read an existing ledger via :class:`LedgerLoader`."""
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(ledger_path)
    e1 = _make_entry("c1")
    e2 = _make_entry("c2")
    ledger.record(e1)
    ledger.record(e2)
    loader = LedgerLoader(ledger_path)
    cells = {entry.cell_id for entry in loader}
    assert cells == {"c1", "c2"}
