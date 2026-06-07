"""N-body Phase B: rung verdict thresholds + non-promoting audit record (plan Phase B).

The rung RECORDS, the human DECIDES (review_queue golden discipline: no auto-
promotion). This test pins both the thresholds and the no-promotion invariant.
"""

from __future__ import annotations

from pathlib import Path

from cyclerfinder.data.review_queue import ReviewQueueEntry
from cyclerfinder.nbody.rung import RungVerdict, record_rung_result, rung_verdict


def test_thresholds() -> None:
    assert rung_verdict(0.15, terminal_closure_km=10.0, converged=True).tier == "ROBUST"
    assert rung_verdict(0.5, terminal_closure_km=10.0, converged=True).tier == "MARGINAL"
    assert rung_verdict(2.0, terminal_closure_km=10.0, converged=True).tier == "ARTIFACT"
    assert rung_verdict(0.05, terminal_closure_km=1e9, converged=False).tier == "ARTIFACT"


def test_record_does_not_promote(tmp_path: Path, silver_fixture: ReviewQueueEntry) -> None:
    v = rung_verdict(2.0, terminal_closure_km=1e6, converged=True)
    out = record_rung_result(silver_fixture, v, tmp_path / "audit.jsonl")
    # Tier on the queue entry is untouched; the rung wrote an AUDIT note only.
    assert silver_fixture.verdict_tier == "silver"
    assert out["rung_verdict"] == "ARTIFACT"
    assert out["promoted"] is False


def test_isinstance_rung_verdict() -> None:
    v = rung_verdict(0.1, terminal_closure_km=5.0, converged=True)
    assert isinstance(v, RungVerdict)
