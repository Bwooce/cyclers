"""N-body Phase B: ingest the two held SILVER candidates (plan Phase B).

NON-GOLDEN: SILVER candidate numerics are OUR computation (unsourced by tier
definition). The rung asserts closure/correction-dV regime, never a sourced V_inf.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cyclerfinder.nbody.rung import load_silver_candidates


def _row(candidate_id: str, verdict_tier: str) -> dict[str, object]:
    return {
        "candidate_id": candidate_id,
        "verdict_tier": verdict_tier,
        "sequence": ["E", "M", "E", "E"],
        "period_k": 2,
        "tof_days": [154.0, 379.0, 1027.0],
        "vinf_per_encounter_kms": [9.75, 13.01, 9.76, 9.75],
        "bend_feasible": True,
        "signature_hash": "h",
        "match_outcome": "novel",
        "known_id": None,
        "superseded_by": [],
        "max_vinf_kms": 13.01,
        "model_assumption": "circular",
        "verdict_audit": {},
        "panel": {},
        "t_added": "2026-06-06T00:00:00Z",
    }


def test_loads_only_silver_entries(tmp_path: Path) -> None:
    qfile = tmp_path / "q.jsonl"
    rows = [
        {
            "candidate_id": "c1",
            "verdict_tier": "silver",
            "sequence": ["E", "M", "E", "E"],
            "period_k": 2,
            "tof_days": [154.0, 379.0, 1027.0],
            "vinf_per_encounter_kms": [9.75, 13.01, 9.76, 9.75],
            "bend_feasible": True,
            "signature_hash": "h",
            "match_outcome": "novel",
            "known_id": None,
            "superseded_by": [],
            "max_vinf_kms": 13.01,
            "model_assumption": "circular",
            "verdict_audit": {},
            "panel": {},
            "t_added": "2026-06-06T00:00:00Z",
        },
        {
            "candidate_id": "c2",
            "verdict_tier": "bronze",
            "sequence": ["E", "M", "E"],
            "period_k": 1,
            "tof_days": [200.0],
            "vinf_per_encounter_kms": [5.0, 5.0],
            "bend_feasible": True,
            "signature_hash": "h2",
            "match_outcome": "novel",
            "known_id": None,
            "superseded_by": [],
            "max_vinf_kms": 5.0,
            "model_assumption": "circular",
            "verdict_audit": {},
            "panel": {},
            "t_added": "2026-06-06T00:00:00Z",
        },
    ]
    qfile.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    cands = load_silver_candidates(qfile)
    assert [c.candidate_id for c in cands] == ["c1"]  # bronze filtered out


def test_gold_filtered_via_enum_round_trip(tmp_path: Path) -> None:
    """GOLD is filtered out: the selection is the enum-round-trip form
    (VerdictTier(...) is VerdictTier.SILVER), not a stale string literal compare.
    """
    qfile = tmp_path / "q.jsonl"
    rows = [_row("g1", "gold"), _row("s1", "silver")]
    qfile.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    cands = load_silver_candidates(qfile)
    assert [c.candidate_id for c in cands] == ["s1"]


def test_malformed_tier_raises_not_silently_dropped(tmp_path: Path) -> None:
    """A malformed verdict_tier raises ValueError from the enum constructor — it is
    NOT silently treated as "not SILVER" and dropped (the I3 round-trip fix: the
    old `== VerdictTier.SILVER.value` compare swallowed corrupt tiers).
    """
    qfile = tmp_path / "q.jsonl"
    qfile.write_text(json.dumps(_row("bad", "platinum")) + "\n")
    with pytest.raises(ValueError, match="platinum"):
        load_silver_candidates(qfile)
