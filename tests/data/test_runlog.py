"""Tests for :mod:`cyclerfinder.data.runlog` (#152 — runlog persistence).

Covers: round-trip, schema validation (reject missing required fields), append
semantics across re-opens, and a slow smoke that
``scripts/campaign_russell12.py`` with ``--runlog-dir tmpdir`` writes >=1 valid
record (driven on the cheapest single row via the free-return genome — pure
residual evals, no Lambert sweep).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from cyclerfinder.data.runlog import (
    RunLog,
    RunlogError,
    RunRecord,
    default_runlog_path,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _sample_record() -> RunRecord:
    return RunRecord(
        row_id="russell-ch4-4.991gG2",
        genome="free-return",
        outcome="CLOSE-AND-MATCH",
        model="circular",
        code_version="abc1234",
        achieved_vinf_kms={"E": 5.55, "M": 3.05},
        sourced_vinf_kms={"E": 5.6, "M": 3.0},
        sourced_anchors={"sourced_aphelion_ratio": 1.2},
        seed={"a_au": 1.3, "e": 0.25},
        residual_kms=0.004,
        solver_audit={"checks": ["vinf[E]: OK"], "n_closed_magnitude": 3},
    )


def test_round_trip(tmp_path: Path) -> None:
    log = RunLog(tmp_path / "rt.jsonl")
    rec = _sample_record()
    log.append(rec)

    loaded = log.read()
    assert len(loaded) == 1
    got = loaded[0]
    assert got.row_id == rec.row_id
    assert got.genome == rec.genome
    assert got.outcome == rec.outcome
    assert got.model == rec.model
    assert got.code_version == rec.code_version
    assert got.achieved_vinf_kms == rec.achieved_vinf_kms
    assert got.sourced_vinf_kms == rec.sourced_vinf_kms
    assert got.sourced_anchors == rec.sourced_anchors
    assert got.seed == rec.seed
    assert got.residual_kms == rec.residual_kms
    assert got.solver_audit == rec.solver_audit
    # A write time is stamped even when the record left it blank.
    assert got.t_written


def test_minimal_record_round_trips(tmp_path: Path) -> None:
    """Only the required fields set; optionals default cleanly on read."""
    log = RunLog(tmp_path / "min.jsonl")
    log.append(
        RunRecord(
            row_id="r1",
            genome="lambert",
            outcome="NO-CLOSE",
            model="astropy",
            code_version="",
        )
    )
    got = log.read()[0]
    assert got.row_id == "r1"
    assert got.achieved_vinf_kms == {}
    assert got.sourced_anchors == {}
    assert got.residual_kms is None


@pytest.mark.parametrize("missing", ["row_id", "genome", "outcome", "model", "code_version"])
def test_reject_missing_required_field(tmp_path: Path, missing: str) -> None:
    payload = {
        "row_id": "r1",
        "genome": "lambert",
        "outcome": "NO-CLOSE",
        "model": "astropy",
        "code_version": "abc1234",
    }
    del payload[missing]
    path = tmp_path / "bad.jsonl"
    path.write_text(json.dumps(payload) + "\n")

    log = RunLog(path)
    with pytest.raises(RunlogError, match=f"missing required field: {missing!r}"):
        log.read()


def test_reject_malformed_json(tmp_path: Path) -> None:
    path = tmp_path / "garbage.jsonl"
    path.write_text("{not valid json\n")
    log = RunLog(path)
    with pytest.raises(RunlogError, match="malformed JSON"):
        log.read()


def test_reject_non_object_line(tmp_path: Path) -> None:
    path = tmp_path / "list.jsonl"
    path.write_text("[1, 2, 3]\n")
    log = RunLog(path)
    with pytest.raises(RunlogError, match="expected a JSON object"):
        log.read()


def test_append_semantics_across_reopen(tmp_path: Path) -> None:
    """Re-opening the same path and appending accumulates, never truncates."""
    path = tmp_path / "append.jsonl"
    RunLog(path).append(_sample_record())

    second = RunLog(path)  # fresh handle on an existing file
    second.append(
        RunRecord(
            row_id="r2",
            genome="lambert",
            outcome="CLOSE-OFF-ANCHOR",
            model="circular",
            code_version="def5678",
        )
    )
    records = RunLog(path).read()
    assert [r.row_id for r in records] == ["russell-ch4-4.991gG2", "r2"]
    assert len(RunLog(path)) == 2


def test_extend_writes_in_order(tmp_path: Path) -> None:
    path = tmp_path / "ext.jsonl"
    log = RunLog(path)
    recs = [
        RunRecord(
            row_id=f"r{i}", genome="lambert", outcome="NO-CLOSE", model="astropy", code_version="v"
        )
        for i in range(3)
    ]
    log.extend(recs)
    assert [r.row_id for r in log.read()] == ["r0", "r1", "r2"]


def test_blank_lines_skipped(tmp_path: Path) -> None:
    path = tmp_path / "blanks.jsonl"
    log = RunLog(path)
    log.append(_sample_record())
    with path.open("a") as fh:
        fh.write("\n   \n")
    assert len(log.read()) == 1


def test_default_runlog_path() -> None:
    p = default_runlog_path("data/runs", "russell12-circular", timestamp="20260607T000000Z")
    assert p == Path("data/runs/russell12-circular-20260607T000000Z.jsonl")


@pytest.mark.slow
def test_campaign_writes_runlog(tmp_path: Path) -> None:
    """Smoke: the campaign with --runlog-dir writes >=1 schema-valid record.

    Driven via the free-return genome (pure residual evals, no Lambert sweep) on
    the circular model — the cheapest like-for-like path. We do not assert any
    scientific outcome, only that a valid runlog lands.
    """
    runlog_dir = tmp_path / "runs"
    proc = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "campaign_russell12.py"),
            "--model",
            "circular",
            "--genome",
            "free-return",
            "--phase-epochs",
            "64",
            "--runlog-dir",
            str(runlog_dir),
            "--runlog-timestamp",
            "testsmoke",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 0, f"campaign failed:\n{proc.stdout}\n{proc.stderr}"

    written = list(runlog_dir.glob("*.jsonl"))
    assert len(written) == 1, f"expected one runlog, got {written}"

    records = RunLog(written[0]).read()
    assert len(records) >= 1
    for rec in records:
        assert rec.row_id
        assert rec.genome == "free-return"
        assert rec.model == "circular"
        assert rec.outcome
