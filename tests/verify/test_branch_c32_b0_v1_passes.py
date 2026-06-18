"""§14 V1 frozen-gate evidence for the #389 branch_C32_b0 candidate.

This module is a FROZEN-GATE wrapper asserting that the #389 P389.1 V1 verdict
JSONL on disk says what the catalogue admission claims. It does NOT recompute
anything: the branched orbit's V1 verdict was built by
``scripts/branch_c32_b0_v1_verify.py`` (P389.1 commit) and the corrector +
independent closure residuals are frozen into the project-output JSONL.

What this gate asserts:

* The #389 P389.1 verdict at ``data/branch_c32_b0_v1_verdict.jsonl`` reports
  ``passes_v1 = True`` for the branch_C32_b0 candidate id.
* The corrector residual is below the V1 floor (1e-10), and the independent
  Radau closure is below the V1 closure floor (1e-6).
* The corrected period matches the Phase 2 sweep period to 8 decimals (the
  IC's structural identity must be preserved by the V1 re-correction).

If this gate breaks, EITHER the JSONL was regenerated with a different answer
(in which case the candidate's V1 registry entry must be revisited) OR the
JSONL file was moved/deleted (in which case the row's evidence pointer is
stale). Either way, the V1 registry entry for
``branch-c32-b0-em-3-3-quasi-cycler-2026`` is what is at stake.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, which
is sourced by the #389 P389.1 closure run (project output, frozen). This test
is the wrapper that ties the catalogue claim to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
VERDICT_PATH = Path("data/branch_c32_b0_v1_verdict.jsonl")
V1_CORRECTOR_FLOOR = 1.0e-10
V1_INDEPENDENT_CLOSURE_FLOOR = 1.0e-6
PHASE2_PERIOD_TU = 23.355184434547017


def _load_v1_verdict() -> dict[str, Any]:
    """Return the V1-verdict row for branch_C32_b0."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v1_verdict_cr3bp_periodic"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                return row
    raise AssertionError(
        f"V1 verdict row for {CANDIDATE_ID!r} not found in {VERDICT_PATH}",
    )


def test_branch_c32_b0_v1_passes_compound_gate() -> None:
    """The candidate's V1 verdict reports passes_v1=True with corrector + closure floors held."""
    row = _load_v1_verdict()
    assert row["passes_v1"] is True, (
        f"branch_C32_b0 V1 verdict reports passes_v1={row['passes_v1']!r} — "
        "registry entry must be re-validated."
    )
    corrector = float(row["corrector_residual"])
    assert corrector < V1_CORRECTOR_FLOOR, (
        f"branch_C32_b0 V1 corrector residual {corrector:.3e} is NOT below the "
        f"V1 floor {V1_CORRECTOR_FLOOR:.0e} — registry entry must be re-validated."
    )
    independent = float(row["independent_closure_residual"])
    assert independent < V1_INDEPENDENT_CLOSURE_FLOOR, (
        f"branch_C32_b0 V1 independent closure {independent:.3e} is NOT below "
        f"the V1 closure floor {V1_INDEPENDENT_CLOSURE_FLOOR:.0e} — registry "
        "entry must be re-validated."
    )
    # Cross-check the recorded floors match this gate's constants — if the
    # JSONL was regenerated under different floors, that must be reviewed.
    assert float(row["v1_corrector_residual_floor"]) == V1_CORRECTOR_FLOOR, (
        f"JSONL v1_corrector_residual_floor ({row['v1_corrector_residual_floor']}) "
        f"drifted from gate floor ({V1_CORRECTOR_FLOOR}) — registry entry must be re-validated."
    )
    assert float(row["v1_independent_closure_floor"]) == V1_INDEPENDENT_CLOSURE_FLOOR, (
        f"JSONL v1_independent_closure_floor ({row['v1_independent_closure_floor']}) "
        f"drifted from gate floor ({V1_INDEPENDENT_CLOSURE_FLOOR}) — registry entry "
        "must be re-validated."
    )


def test_branch_c32_b0_v1_preserves_phase2_period() -> None:
    """The corrected period matches the Phase 2 sweep period to 8 decimals.

    The V1 re-correction must NOT migrate the orbit to a different period —
    the IC's structural identity (T=101.56d / (3, 3) topology) is what the
    catalogue admission rests on. A period drift past 8 decimals would indicate
    the V1 corrector landed on a different member of the family or a different
    family entirely; that would invalidate the topology / jacobi claims.
    """
    row = _load_v1_verdict()
    corrected = float(row["period_corrected_TU"])
    assert abs(corrected - PHASE2_PERIOD_TU) < 1e-8, (
        f"branch_C32_b0 V1 corrected period {corrected:.15f} TU diverged from the "
        f"Phase 2 period {PHASE2_PERIOD_TU:.15f} TU by "
        f"{abs(corrected - PHASE2_PERIOD_TU):.3e} TU — structural identity broken."
    )


def test_branch_c32_b0_v1_orbit_is_planar() -> None:
    """The corrected orbit retains the Phase 2 degenerate-planar character.

    branch_C32_b0 was classified ``degenerate_planar=True`` in Phase 2
    (z0~1.7e-22, zdot0~8.3e-24). The V1 re-correction must NOT spawn an out-of-
    plane component — that would change the topology label from planar-(3,3) to
    a 3D orbit and break the catalogue row's published characterisation.
    """
    row = _load_v1_verdict()
    assert row["degenerate_planar"] is True, (
        f"branch_C32_b0 V1 corrected orbit reports degenerate_planar="
        f"{row['degenerate_planar']!r} — planar (3, 3) topology claim broken."
    )
