"""§14 P389.2 frozen-gate evidence for branch_C32_b0 — physical / lit / ML.

This module is a FROZEN-GATE wrapper asserting that the #389 P389.2 verdict
JSONL on disk says what the catalogue admission claims. Three sub-gates are
asserted:

1. Physical-sanity (#324): every of the 3 lunar close-approaches in the
   (3, 3) cycle reports a patched-conic max_bend >= 5° at the Moon's safe
   periapsis. NOTE: the actual close-approach distances are far above the
   Moon's Hill sphere — branch_C32_b0 is a far-amplitude bound Earth-system
   orbit, NOT a lunar-flyby tour. The structural-feasibility check is asserted
   in case a future mission design DID place a spacecraft near the Moon at
   the indicated V_∞.

2. Lit-fresh (#346/#349): the offline check_literature returned a non-
   ``published`` status against the present KNOWN_CORPUS (no Earth-Moon CR3BP
   anchor publishes a (3, 3) planar cycler at jacobi=3.797). Status
   ``inconclusive`` is acceptable here per the SILVER precedent — the project
   has no live web search lane.

3. ML flagger (#256): the false-positive flagger scores branch_C32_b0 below
   the 0.75 spec §16.5 routing threshold.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, which
is sourced by the #389 P389.2 combined run (project output, frozen). This test
ties the catalogue claim to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
VERDICT_PATH = Path("data/branch_c32_b0_p389_2_verdict.jsonl")
N_LUNAR_ENCOUNTERS_EXPECTED = 3  # (3, 3) topology -> 3 lunar windings per cycle
MIN_BEND_DEG = 5.0
ML_FP_THRESHOLD = 0.75


def _load_verdict() -> dict[str, Any]:
    """Return the ``p389_2_verdict`` row from the JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "p389_2_verdict" and row.get("candidate_id") == CANDIDATE_ID:
                return row
    raise AssertionError(
        f"P389.2 verdict row for {CANDIDATE_ID!r} not found in {VERDICT_PATH}",
    )


def test_branch_c32_b0_p389_2_physical_sanity_passes_all_lunar_encounters() -> None:
    """All 3 lunar close-approaches clear the 5° max-bend feasibility floor."""
    row = _load_verdict()
    assert row["physical_sanity_pass"] is True, (
        f"branch_C32_b0 physical-sanity gate reports "
        f"physical_sanity_pass={row['physical_sanity_pass']!r}; expected True."
    )
    encounters = row["physical_per_encounter"]
    assert len(encounters) == N_LUNAR_ENCOUNTERS_EXPECTED, (
        f"Expected {N_LUNAR_ENCOUNTERS_EXPECTED} lunar encounters (one per "
        f"(3, 3) winding); got {len(encounters)}."
    )
    for enc in encounters:
        bend = float(enc["max_bend_deg"])
        assert bend >= MIN_BEND_DEG, (
            f"Encounter {enc['encounter_index']} max_bend {bend:.3f} deg is below "
            f"the {MIN_BEND_DEG:.2f} deg structural-feasibility floor."
        )
        assert enc["is_useful"] is True, (
            f"Encounter {enc['encounter_index']} reports is_useful="
            f"{enc['is_useful']!r}; expected True."
        )


def test_branch_c32_b0_p389_2_lit_fresh_not_published() -> None:
    """The offline check_literature returned a non-``published`` status.

    Status ``not-found`` (clean offline miss) and ``inconclusive`` (no results
    returned) both clear here — the project's web-search lane is not wired in
    this run, mirroring the SILVER's #328 offline-corpus precedent.
    """
    row = _load_verdict()
    lit = row["lit_fresh"]
    assert lit["lit_check_status"] in ("not-found", "inconclusive"), (
        f"branch_C32_b0 lit-fresh check reports status="
        f"{lit['lit_check_status']!r}; expected not-found or inconclusive."
    )
    sig = lit["candidate_signature"]
    assert sig["primary"] == "Earth", f"signature primary={sig['primary']!r}; expected Earth"
    assert sig["sequence"] == ["Moon"], f"signature sequence={sig['sequence']!r}; expected ['Moon']"
    assert "repeated-moon" in sig["topology_label"], (
        f"signature topology_label={sig['topology_label']!r}; expected to include 'repeated-moon'"
    )


def test_branch_c32_b0_p389_2_ml_flagger_below_threshold() -> None:
    """The ML false-positive flagger scores below the 0.75 routing threshold."""
    row = _load_verdict()
    ml = row["ml_flagger"]
    p_fp = float(ml["p_false_positive"])
    assert p_fp <= ML_FP_THRESHOLD, (
        f"branch_C32_b0 ML flagger p_fp={p_fp:.4f} exceeds spec §16.5 routing "
        f"threshold {ML_FP_THRESHOLD:.2f}."
    )
    assert ml["ml_passes"] is True, (
        f"ML flagger ml_passes={ml['ml_passes']!r}; expected True given p_fp under threshold."
    )


def test_branch_c32_b0_p389_2_combined_verdict_passes() -> None:
    """The combined P389.2 verdict reports all three gates passed."""
    row = _load_verdict()
    assert row["p389_2_passes"] is True, (
        f"branch_C32_b0 P389.2 combined verdict reports p389_2_passes="
        f"{row['p389_2_passes']!r}; expected True. Inspect physical/lit/ML "
        "sub-verdicts for the failing gate."
    )
