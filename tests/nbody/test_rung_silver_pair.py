"""N-body Phase B: SILVER rung run on the two held candidates (plan Phase B Task B.4).

End-to-end: propagate each held candidate one full period (Sun+E+M+J), compute the
node-impulse correction ΔV, grade the verdict, and RECORD it (promoted=False) to a
review-queue audit trail. NON-GOLDEN: the candidates are unsourced (SILVER); the
rung asserts the result is RECORDED with a regime verdict, never a sourced value.

Honesty boundary (plan): both candidates float in the high-V∞ basin (E∞ ~9.7 /
M∞ ~12-13 km/s), so a MARGINAL/ARTIFACT verdict is the EXPECTED, valid outcome —
the rung's job is to record, not to make them pass.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.data.review_queue import ReviewQueueEntry  # noqa: E402
from cyclerfinder.nbody.rung import run_rung  # noqa: E402


@pytest.mark.slow
def test_rung_runs_and_records_both_candidates(
    tmp_path: Path,
    silver_fixture: ReviewQueueEntry,
    silver_candidate_2: ReviewQueueEntry,
) -> None:
    audit = tmp_path / "rung_audit.jsonl"
    ephem = Ephemeris("astropy")
    records = []
    for entry in (silver_fixture, silver_candidate_2):
        rec = run_rung(entry, ephem, audit_path=audit, accuracy=1e-9)
        records.append(rec)

    # Both recorded, neither promoted, queue tiers untouched.
    assert len(records) == 2
    for rec in records:
        assert rec["promoted"] is False
        assert rec["rung_verdict"] in {"ROBUST", "MARGINAL", "ARTIFACT"}
    assert silver_fixture.verdict_tier == "silver"
    assert silver_candidate_2.verdict_tier == "silver"

    # The audit trail has one line per candidate.
    lines = [ln for ln in audit.read_text().splitlines() if ln.strip()]
    assert len(lines) == 2
    ids = {json.loads(ln)["candidate_id"] for ln in lines}
    assert ids == {"forge-silver-1", "forge-silver-2"}


@pytest.mark.slow
def test_jupiter_sensitivity_recorded(
    silver_fixture: ReviewQueueEntry,
) -> None:
    """Gate-4 body-inclusion arm against a real candidate baseline (design §2)."""
    from cyclerfinder.nbody.rung import jupiter_sensitivity

    sens = jupiter_sensitivity(silver_fixture, Ephemeris("astropy"), accuracy=1e-9)
    assert "with_jupiter" in sens and "without_jupiter" in sens
    assert "delta_correction_dv_kms" in sens
