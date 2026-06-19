"""§14 V4 annual-sweep frozen-gate for branch_C32_b0 — STRUCTURAL_FAIL_ALL_EPOCHS.

The #338-pattern annual launch-epoch sweep (2000-2099) on branch_C32_b0's V4.
Where the #338 SILVER sweep found an interior PASS run (EFFECTIVELY_CYCLIC),
branch_C32_b0 FAILS V4 at every single 21st-century launch epoch: the real
DE440 solar tide destabilizes the far-amplitude orbit into escape (~1e9 km
drift) regardless of launch phase. This upgrades the P389.5 single-epoch HALT to
a complete structural negative — the candidate cannot be a real-ephemeris cycler
at any 21st-century launch epoch and stays unadmitted.

This frozen-gate test ties the recorded sweep verdict to that conclusion: any
change that produces even one PASS epoch trips a failure and forces a re-review.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from P389.6.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

VERDICT_PATH = Path("data/branch_c32_b0_v4_annual_sweep_389.jsonl")
CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
EXPECTED_VERDICT_LABEL = "STRUCTURAL_FAIL_ALL_EPOCHS"
EXPECTED_N_EPOCHS = 100
V4_ESCAPE_FLOOR_KMS = 1.0e8


def _load_records() -> dict[str, Any]:
    assert VERDICT_PATH.exists(), f"frozen sweep file missing: {VERDICT_PATH}"
    records: dict[str, Any] = {"rows": []}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            kind = row.get("kind")
            if kind == "annual_sweep_row":
                records["rows"].append(row)
            elif kind:
                records[kind] = row
    return records


def test_sweep_is_structural_fail_at_every_epoch() -> None:
    """0/100 epochs PASS; the boundary verdict is STRUCTURAL_FAIL_ALL_EPOCHS."""
    records = _load_records()
    rows = records["rows"]
    assert len(rows) == EXPECTED_N_EPOCHS, (
        f"expected {EXPECTED_N_EPOCHS} sweep epochs, got {len(rows)}"
    )
    n_pass = sum(1 for r in rows if r["passes_v4"])
    assert n_pass == 0, (
        f"{n_pass} epoch(s) now PASS V4 — the structural-fail HALT has flipped; "
        "re-review before any admission."
    )
    assert "boundary_verdict" in records, f"boundary_verdict missing from {VERDICT_PATH}"
    bv = records["boundary_verdict"]
    assert bv["verdict_label"] == EXPECTED_VERDICT_LABEL, (
        f"sweep verdict {bv['verdict_label']!r} != {EXPECTED_VERDICT_LABEL!r}"
    )
    assert bv["epoch_dependent"] is False, (
        "sweep now reports epoch_dependent=True — the failure is no longer "
        "uniformly structural; re-review."
    )


def test_sweep_every_epoch_is_an_escape() -> None:
    """Every epoch's drift is an escape (>= 1e8 km), not a marginal miss."""
    records = _load_records()
    for r in records["rows"]:
        assert float(r["max_v4_drift_kms"]) >= V4_ESCAPE_FLOOR_KMS, (
            f"epoch {r['year']} drift {r['max_v4_drift_kms']:.3e} km "
            f"< {V4_ESCAPE_FLOOR_KMS:.0e} km — no longer a clear escape."
        )


def test_sweep_footer_consistent() -> None:
    """The footer mirrors the boundary verdict (0 PASS / 100 FAIL)."""
    records = _load_records()
    assert "footer" in records, f"footer missing from {VERDICT_PATH}"
    footer = records["footer"]
    assert footer["candidate_id"] == CANDIDATE_ID
    assert int(footer["n_pass"]) == 0
    assert int(footer["n_fail"]) == EXPECTED_N_EPOCHS
    assert footer["verdict_label"] == EXPECTED_VERDICT_LABEL
