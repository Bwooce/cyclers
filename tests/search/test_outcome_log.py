"""Passive solver-outcome logger (#210).

Covers: opt-in via ``CYCLERFINDER_OUTCOME_LOG`` (records appended, well-formed,
documented fields present), the disabled NO-OP (no file, no error), numpy
coercion, and one real solver drive (``correct_periodic`` on the Arenstorf seed)
that produces a well-formed record without changing the solver's behaviour.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.outcome_log as outcome_log

# Sourced Arenstorf golden (Hairer, Nørsett, Wanner, "Solving ODEs I", p. 129).
MU = 0.012277471
X0, VY0, PERIOD = 0.994, -2.0015851063790825, 17.0652165601579625


def _read_records(path: Path) -> list[dict[str, Any]]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def test_disabled_is_noop(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With the env var unset, log_outcome writes nothing and never raises."""
    monkeypatch.delenv(outcome_log.ENV_VAR, raising=False)
    target = tmp_path / "should_not_exist.jsonl"
    # Several calls — all must be no-ops.
    for _ in range(3):
        outcome_log.log_outcome(
            solver="test.solver",
            inputs={"a": 1},
            outcome={"converged": True},
        )
    assert not target.exists()
    # No file should have been created anywhere in tmp_path.
    assert list(tmp_path.iterdir()) == []


def test_enabled_appends_well_formed_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    log_path = tmp_path / "nested" / "outcomes.jsonl"
    monkeypatch.setenv(outcome_log.ENV_VAR, str(log_path))

    outcome_log.log_outcome(
        solver="solverA",
        inputs={"x": 1.0, "vec": [1, 2, 3]},
        outcome={"converged": True, "cost": 0.5},
        meta={"note": "first"},
    )
    outcome_log.log_outcome(
        solver="solverB",
        inputs={"y": 2.0},
        outcome={"converged": False, "cost": 9.9},
    )

    assert log_path.exists()  # parent dir created on demand
    records = _read_records(log_path)
    assert len(records) == 2

    for rec in records:
        # Documented envelope fields.
        assert rec["schema_version"] == outcome_log.SCHEMA_VERSION
        assert isinstance(rec["counter"], int)
        assert "wall_time" in rec
        assert "solver" in rec
        assert "inputs" in rec
        assert "outcome" in rec
        assert "meta" in rec

    # Monotonic counter strictly increases in append order.
    assert records[0]["counter"] < records[1]["counter"]
    assert records[0]["solver"] == "solverA"
    assert records[0]["meta"] == {"note": "first"}
    assert records[1]["meta"] is None
    assert records[1]["outcome"]["converged"] is False


def test_numpy_coercion(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    log_path = tmp_path / "np.jsonl"
    monkeypatch.setenv(outcome_log.ENV_VAR, str(log_path))

    outcome_log.log_outcome(
        solver="np.solver",
        inputs={
            "scalar": np.float64(1.5),
            "intval": np.int64(7),
            "flag": np.bool_(True),
            "arr": np.array([1.0, 2.0, 3.0]),
        },
        outcome={"residual": np.float64(1e-12)},
    )
    records = _read_records(log_path)
    assert len(records) == 1
    inp = records[0]["inputs"]
    assert inp["scalar"] == 1.5
    assert inp["intval"] == 7
    assert inp["flag"] is True
    assert inp["arr"] == [1.0, 2.0, 3.0]
    assert records[0]["outcome"]["residual"] == 1e-12


def test_to_jsonable_handles_complex_and_fallback() -> None:
    out = outcome_log.to_jsonable(complex(1.0, -2.0))
    assert out == {"real": 1.0, "imag": -2.0}
    # An arbitrary object stringifies rather than raising.
    assert isinstance(outcome_log.to_jsonable(object()), str)


def test_real_solver_drive_produces_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Drive correct_periodic on the Arenstorf seed; assert the wired call-site
    emits one well-formed record AND the solver result is unchanged."""
    log_path = tmp_path / "solver.jsonl"
    monkeypatch.setenv(outcome_log.ENV_VAR, str(log_path))

    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0])
    res = cp.correct_periodic(sysm, s0, PERIOD)

    # Behaviour unchanged: the corrector still converges as in the golden test.
    assert res.converged
    assert res.closure_residual < 1e-6

    records = _read_records(log_path)
    assert len(records) == 1
    rec = records[0]
    assert rec["solver"] == "cr3bp.correct_periodic"
    assert rec["inputs"]["mu"] == MU
    assert len(rec["inputs"]["state0_guess"]) == 6
    assert rec["inputs"]["period_guess"] == PERIOD
    assert rec["outcome"]["converged"] is True
    assert rec["outcome"]["residual"] < 1e-6
    assert "period" in rec["outcome"]


def test_real_solver_disabled_no_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Driving the solver with logging OFF leaves no file and still converges."""
    monkeypatch.delenv(outcome_log.ENV_VAR, raising=False)
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0])
    res = cp.correct_periodic(sysm, s0, PERIOD)
    assert res.converged
    assert list(tmp_path.iterdir()) == []
