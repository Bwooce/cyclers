"""§14 V4 frozen-gate evidence for the #327 SILVER (URA111 real-eph + epoch sweep).

This module is a FROZEN-GATE wrapper asserting that the #338 V4-strict annual
epoch sweep boundary verdict on disk says what the catalogue row claims. The
V4 question is whether the V3 bounded-drift signature survives real-ephemeris
Uranian dynamics (JPL/NAIF URA111 SPICE kernel with real Umbriel/Oberon
eccentricities and inclinations). The annual sweep across 2000-2099 returned
94/100 PASS, with 6/6 failures clustered in the last 15 yr of URA111 coverage
(a kernel-edge / extrapolation artifact). The interior PASS rate is 85/85 =
100% — verdict ``EFFECTIVELY_CYCLIC``, recommended launch in the 84-yr PASS
run 2000-2083.

What this gate asserts:

* The #338 Part B boundary verdict at
  ``data/silver_327_v4_strict_boundary_338.jsonl`` reports
  ``verdict_label = EFFECTIVELY_CYCLIC``.
* All 6 failures sit in the kernel-edge window (2084-2099) — interior PASS
  count is 85/85 = 100% (``n_fails_interior == 0`` in the kernel-edge breakdown
  record).
* The longest PASS run is at least 84 yr (``longest_pass_run_length_yr >= 84``)
  — the basis for the catalogue row's 2000-2083 validity_window.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from #338. This test ties the catalogue row's V4 claim to the
recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

VERDICT_PATH = Path("data/silver_327_v4_strict_boundary_338.jsonl")
EXPECTED_VERDICT_LABEL = "EFFECTIVELY_CYCLIC"
EXPECTED_LONGEST_PASS_RUN_YR = 84
EXPECTED_LONGEST_PASS_START_YEAR = 2000
EXPECTED_LONGEST_PASS_END_YEAR = 2083


def _load_v4_records() -> dict[str, dict[str, Any]]:
    """Return ``{kind: record}`` for the #338 boundary JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    records: dict[str, dict[str, Any]] = {}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            kind = row.get("kind")
            if kind:
                records[kind] = row
            elif row.get("_meta"):
                records["_meta"] = row
    return records


def test_silver_327_v4_strict_boundary_verdict_is_effectively_cyclic() -> None:
    """The #338 boundary verdict_label is EFFECTIVELY_CYCLIC."""
    records = _load_v4_records()
    assert "boundary_verdict" in records, (
        f"#338 boundary_verdict record missing from {VERDICT_PATH}"
    )
    verdict = records["boundary_verdict"]
    assert verdict["verdict_label"] == EXPECTED_VERDICT_LABEL, (
        f"#338 V4-strict boundary verdict {verdict['verdict_label']!r} != "
        f"{EXPECTED_VERDICT_LABEL!r} — registry entry must be re-validated."
    )
    assert verdict["kernel_edge_dominated"] is True, (
        "#338 verdict no longer reports kernel_edge_dominated=True — the "
        "fail-mode argument has changed."
    )


def test_silver_327_v4_strict_interior_pass_rate_is_100pct() -> None:
    """All 6 failures cluster at the URA111 kernel edge; interior is 85/85 = 100%."""
    records = _load_v4_records()
    assert "kernel_edge_breakdown" in records, (
        f"#338 kernel_edge_breakdown record missing from {VERDICT_PATH}"
    )
    breakdown = records["kernel_edge_breakdown"]
    assert int(breakdown["n_fails_interior"]) == 0, (
        f"#338 V4-strict interior failures: {breakdown['n_fails_interior']} "
        "(expected 0 — interior PASS rate is the basis of the EFFECTIVELY_CYCLIC "
        "verdict)."
    )
    n_fails_total = int(breakdown["n_fails_total"])
    n_fails_edge = int(breakdown["n_fails_near_edge"])
    assert n_fails_total == n_fails_edge, (
        f"#338 V4-strict failures: {n_fails_total} total != {n_fails_edge} edge — "
        "the kernel-edge artifact attribution is broken."
    )


def test_silver_327_v4_strict_longest_pass_run_is_at_least_84_years() -> None:
    """The longest PASS run is at least 84 yr (basis for the validity_window)."""
    records = _load_v4_records()
    assert "longest_pass_run" in records, (
        f"#338 longest_pass_run record missing from {VERDICT_PATH}"
    )
    run = records["longest_pass_run"]
    length_yr = int(run["longest_pass_run_length_yr"])
    assert length_yr >= EXPECTED_LONGEST_PASS_RUN_YR, (
        f"#338 longest PASS run is {length_yr} yr < {EXPECTED_LONGEST_PASS_RUN_YR} yr — "
        "the catalogue row's 84-yr validity_window claim is no longer backed."
    )
    assert int(run["longest_pass_run_start_year"]) == EXPECTED_LONGEST_PASS_START_YEAR, (
        f"#338 longest PASS run start {run['longest_pass_run_start_year']} != "
        f"{EXPECTED_LONGEST_PASS_START_YEAR}"
    )
    assert int(run["longest_pass_run_end_year"]) == EXPECTED_LONGEST_PASS_END_YEAR, (
        f"#338 longest PASS run end {run['longest_pass_run_end_year']} != "
        f"{EXPECTED_LONGEST_PASS_END_YEAR}"
    )
