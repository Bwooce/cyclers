"""N-body Phase B: ingest the two held SILVER candidates (plan Phase B).

NON-GOLDEN: SILVER candidate numerics are OUR computation (unsourced by tier
definition). The rung asserts closure/correction-dV regime, never a sourced V_inf.
"""

from __future__ import annotations

import json
from pathlib import Path

from cyclerfinder.nbody.rung import load_silver_candidates


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
