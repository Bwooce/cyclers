"""M8-UX Phase 4: `cyclerfinder report` writes md+json from a ledger."""

from __future__ import annotations

from pathlib import Path

import pytest

from cyclerfinder.cli import main
from cyclerfinder.data.ledger import Ledger, LedgerEntry


def _seed(path: Path) -> None:
    led = Ledger(path)
    led.record(
        LedgerEntry(
            cell_id="EM|E-M-E|k2|r00|bss",
            status="solved",
            n_solutions=1,
            best_dv_kms=0.012,
            signature_hashes=("sha1:deadbeef",),
            validation_level="V2",
            t_done="2026-06-06T00:00:00Z",
            host="test",
        )
    )


def test_report_writes_both_files(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    ledger = tmp_path / "c.jsonl"
    _seed(ledger)
    stem = tmp_path / "r"
    code = main(
        ["report", "--ledger", str(ledger), "--out", str(stem), "--format", "both"]
    )
    assert code == 0
    assert stem.with_suffix(".md").exists()
    assert stem.with_suffix(".json").exists()
    assert stem.with_suffix(".md").read_text().startswith("# Cycler campaign report")


def test_report_with_verdicts(tmp_path: Path) -> None:
    ledger = tmp_path / "c.jsonl"
    _seed(ledger)
    stem = tmp_path / "r"
    code = main(
        [
            "report",
            "--ledger",
            str(ledger),
            "--out",
            str(stem),
            "--format",
            "md",
            "--with-verdicts",
        ]
    )
    assert code == 0
    md = stem.with_suffix(".md").read_text()
    assert "verdict_tier" in md


def test_report_empty_ledger_exits_five(tmp_path: Path) -> None:
    ledger = tmp_path / "empty.jsonl"
    ledger.touch()
    stem = tmp_path / "r"
    code = main(["report", "--ledger", str(ledger), "--out", str(stem), "--format", "both"])
    assert code == 5
