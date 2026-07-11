"""#569 §14 V4 frozen-gate evidence for the 5 #563 symmetric-closure representatives.

This module is a FROZEN-GATE wrapper asserting that the catalogue V4 claim on
the five sibling representatives of the #312 Uranian symmetric-closure
quasi-cycler family says exactly what the recorded project output on disk says.
Each of the five rows (Titania-Oberon / Ariel-Oberon / Umbriel-Titania /
Ariel-Titania / Ariel-Umbriel) is registered at #312-equivalent V4 (windowed)
in ``src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE``; this gate ties that
claim to two frozen files:

* ``data/gauntlet_566_five_representatives.jsonl`` (#566): the
  V2->V3->V4-scipy->V4-strict gauntlet headline. Each candidate must report
  ``chain_verdict == "PASS_AS_QUASI_CYCLER"`` with ``v4_strict_all_pass`` true
  and ``v2_status == "FAIL_QUASI_BOUNDED"`` — the same admission path #312
  itself used (strict V2 fails the 50,000 km strict-periodic floor by design;
  promotion is on the structural V3/V4/V4-strict verdicts).
* ``data/scan_567_epoch_robustness.jsonl`` (#567): the daily-2000 + daily-2030
  epoch scan. The feasible synodic DUTY CYCLE
  ``n_PASS / (n_PASS + n_planet_crossing_infeasible)`` over the ``daily`` sweep
  must equal the ``synodic_duty_cycle_pct`` recorded in each catalogue row's
  ``validity_window`` (per the #568 verdict Sec 2). This is the recomputed,
  non-circular check: the EXPECTED side is the frozen scan, the row merely
  records the number the scan produces.

Sourced-golden discipline: the EXPECTED side is the frozen project output
(#566 gauntlet + #567 scan). This test is intentionally NOT ``@pytest.mark.slow``
— it only reads JSONL (no propagation), runs in well under a second, and must
stay in the default suite so the V4 claim is verified in CI, not skipped.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

GAUNTLET_PATH = Path("data/gauntlet_566_five_representatives.jsonl")
SCAN_PATH = Path("data/scan_567_epoch_robustness.jsonl")

# candidate_id (in the frozen JSONLs) -> (catalogue row id, expected duty cycle %).
# Duty cycles are the #568 verdict Sec 2 table, independently recomputed here.
REPRESENTATIVES: dict[str, tuple[str, float]] = {
    "enum563-line57-titania-oberon-titania": (
        "titania-oberon-1-1-uranian-quasi-cycler-2026",
        74.4,
    ),
    "enum563-line18-ariel-oberon-ariel": (
        "ariel-oberon-1-1-uranian-quasi-cycler-2026",
        71.0,
    ),
    "enum563-line26-umbriel-titania-umbriel": (
        "umbriel-titania-1-1-uranian-quasi-cycler-2026",
        68.5,
    ),
    "enum563-line12-ariel-titania-ariel": (
        "ariel-titania-1-1-uranian-quasi-cycler-2026",
        66.5,
    ),
    "enum563-line2-ariel-umbriel-ariel": (
        "ariel-umbriel-1-1-uranian-quasi-cycler-2026",
        61.7,
    ),
}


def _load_gauntlet_headline() -> dict[str, dict[str, Any]]:
    """Return ``{candidate_id: candidate_summary}`` from the #566 headline record."""
    assert GAUNTLET_PATH.exists(), f"frozen #566 gauntlet file missing: {GAUNTLET_PATH}"
    headline: dict[str, Any] | None = None
    with GAUNTLET_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "headline":
                headline = row
    assert headline is not None, f"#566 headline record missing from {GAUNTLET_PATH}"
    return {c["candidate_id"]: c for c in headline["candidates"]}


def _daily_duty_cycles() -> dict[str, float]:
    """Recompute each candidate's feasible synodic duty cycle from the #567 scan.

    duty = n_PASS / (n_PASS + n_planet_crossing_infeasible) over ``sweep_type ==
    'daily'`` (the combined daily-2000 + daily-2030 windows), matching the #568
    verdict Sec 2 (drift-floor fails excluded from the denominator).
    """
    assert SCAN_PATH.exists(), f"frozen #567 scan file missing: {SCAN_PATH}"
    n_pass: dict[str, int] = {}
    n_pc: dict[str, int] = {}
    with SCAN_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("_meta") or row.get("sweep_type") != "daily":
                continue
            cid = row.get("candidate_id")
            if cid not in REPRESENTATIVES:
                continue
            if row.get("passes_v4_strict"):
                n_pass[cid] = n_pass.get(cid, 0) + 1
            elif row.get("epoch_failure_mode") == "planet_crossing_infeasible":
                n_pc[cid] = n_pc.get(cid, 0) + 1
    duty: dict[str, float] = {}
    for cid in REPRESENTATIVES:
        p = n_pass.get(cid, 0)
        f = n_pc.get(cid, 0)
        assert p + f > 0, f"no daily rows for {cid} in {SCAN_PATH}"
        duty[cid] = 100.0 * p / (p + f)
    return duty


@pytest.mark.parametrize("candidate_id", list(REPRESENTATIVES))
def test_566_representative_passes_v4_strict_as_quasi_cycler(candidate_id: str) -> None:
    """Each representative reaches PASS_AS_QUASI_CYCLER through V4-strict (#566)."""
    summaries = _load_gauntlet_headline()
    assert candidate_id in summaries, f"#566 headline missing candidate {candidate_id}"
    s = summaries[candidate_id]
    assert s["chain_verdict"] == "PASS_AS_QUASI_CYCLER", (
        f"{candidate_id}: #566 chain_verdict {s['chain_verdict']!r} != "
        "'PASS_AS_QUASI_CYCLER' — the catalogue V4 claim is no longer backed."
    )
    assert s["v4_strict_all_pass"] is True, f"{candidate_id}: #566 v4_strict_all_pass is not True."
    assert s["v2_status"] == "FAIL_QUASI_BOUNDED", (
        f"{candidate_id}: #566 v2_status {s['v2_status']!r} != 'FAIL_QUASI_BOUNDED' "
        "(the bounded-drift quasi_cycler admission path #312 itself used)."
    )


@pytest.mark.parametrize("candidate_id", list(REPRESENTATIVES))
def test_566_representative_duty_cycle_matches_catalogue(candidate_id: str) -> None:
    """The recomputed #567 daily duty cycle equals the catalogue row's value."""
    duty = _daily_duty_cycles()
    _, expected_pct = REPRESENTATIVES[candidate_id]
    actual = duty[candidate_id]
    assert abs(actual - expected_pct) < 0.1, (
        f"{candidate_id}: recomputed #567 daily duty cycle {actual:.1f}% != "
        f"catalogue validity_window.synodic_duty_cycle_pct {expected_pct}% "
        "(within 0.1) — the #568 characterization has drifted from the row."
    )
