"""§14 V1 frozen-gate evidence for the #327 SILVER (Umbriel-Oberon-Umbriel).

This module is a FROZEN-GATE wrapper asserting that the V1-aspect of the #306
Phase 1 Part D / #331 V3 gauntlet verdict JSONL on disk says what the catalogue
row claims. It does NOT recompute anything: the SILVER's full provenance was
built by #306 / #330 / #331 / #332 / #335 / #338 and the per-cycle, per-leg
residuals are frozen into project-output JSONLs.

What this gate asserts:

* The #306 Phase 1 Part D verdict at
  ``data/silver_327_v1_v2_verdicts.jsonl`` says ``passes_v1 = True`` for the
  SILVER candidate id, AND the independent cross-check arrival residual
  (DOP853 re-propagation vs. Lambert leg endpoints) is < 1e-3 m/s, the spec
  §14 V1 floor (``v1_floor_kms`` in the verdict, the moontour's per-leg
  independent-cross-check headline).

If this gate breaks, EITHER the JSONL was regenerated with a different answer
(in which case the SILVER row's V1 evidence registry entry must be revisited)
OR the JSONL file was moved/deleted (in which case the row's evidence pointer
is stale). Either way, the V1 registry entry for
``umbriel-oberon-1-1-uranian-quasi-cycler-2026`` is what is at stake.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, which
is sourced by the #306 / #331 verdict chain (project output, frozen). This
test is the wrapper that ties the catalogue claim to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "repeated-moon-uranus-00000041"
VERDICT_PATH = Path("data/silver_327_v1_v2_verdicts.jsonl")
V1_FLOOR_KMS = 1.0e-3  # spec §14 V1: agreement < 1 m/s


def _load_v1_verdict() -> dict[str, Any]:
    """Return the V1-verdict row for the SILVER from the #306 Phase 1 Part D JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v1_verdict_3d_moontour"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                return row
    raise AssertionError(
        f"V1 verdict row for {CANDIDATE_ID!r} not found in {VERDICT_PATH}",
    )


def test_silver_327_v1_passes_independent_crosscheck() -> None:
    """The SILVER's V1 verdict reports passes_v1=True with agreement < 1 m/s."""
    row = _load_v1_verdict()
    assert row["passes_v1"] is True, (
        f"SILVER V1 verdict reports passes_v1={row['passes_v1']!r} — "
        "registry entry must be re-validated."
    )
    # The moontour's V1 headline is the independent cross-check arrival residual
    # in km/s (DOP853 re-propagation vs the Lambert-derived velocity at next
    # encounter). The spec §14 V1 floor is 1 m/s = 1e-3 km/s.
    headline_kms = float(row["headline_kms"])
    assert headline_kms < V1_FLOOR_KMS, (
        f"SILVER V1 headline {headline_kms:.3e} km/s is NOT below the spec §14 "
        f"V1 floor {V1_FLOOR_KMS:.0e} km/s — registry entry must be re-validated."
    )
    # Cross-check the recorded v1_floor_kms matches our local constant — if the
    # JSONL was regenerated under a different spec floor, that must be reviewed.
    assert float(row["v1_floor_kms"]) == V1_FLOOR_KMS, (
        f"JSONL v1_floor_kms ({row['v1_floor_kms']}) drifted from spec "
        f"§14 V1 floor ({V1_FLOOR_KMS}) — registry entry must be re-validated."
    )
