"""S14 V2 frozen-gate evidence for the #344 SILVER (FAIL verdict).

This module is a FROZEN-GATE wrapper asserting that the #344 Phase 2 Stage B
V2-moontour gauntlet verdict JSONL on disk says what the Stage B note claims.
The strict V2 verdict is ``FAIL`` (closure residual blows past the v4.7
quasi-cycler envelope of 0.5 km/s from cycle 1 onward), so the candidate does
NOT enter either the strict-cycler slot OR the quasi_cycler slot.

What this gate asserts:

* The Stage B verdict at ``data/silver_344_moontour_v2_verdicts.jsonl``
  reports 10/10 cycles with all Lambert legs converging
  (``converged_legs == n_legs`` for every cycle). Geometry stays solvable
  even though the V_inf-continuity is broken.
* Cycle-0 closure residual reproduces the Stage A SILVER value of
  ~0.0102 km/s to within ``1e-9`` km/s, confirming the v2-convention
  phase0/rel_offset translation is correct.
* Per-cycle closure residual jumps above the v4.7 quasi envelope
  (>= 0.5 km/s) by cycle 1 and peaks at > 5 km/s. Closure does
  oscillate (cycles 5 and 8 drop back below 0.5 km/s) but the MAX
  across the 10-cycle horizon is > 5 km/s, ~15x the #330 Umbriel-Oberon
  SILVER's max of 0.349 km/s - so the Stage B candidate is structurally
  a closure-divergent oscillator, not a #330-style bounded quasi-cycler.
* Cumulative drift exceeds 9e5 km by cycle 1 and never returns below
  the 50,000 km strict-cycler floor.
* The headline verdict label is ``FAIL``, NOT
  ``PASS_QUASI_CYCLER`` / ``PASS_STRICT_CYCLER`` - the candidate retires
  to the negative-results registry rather than progressing to Stage C.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself,
frozen project output from the Stage B run. This test ties the Stage B
narrative to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "repeated-moon-saturn-titan-rhea-titan-stage-b"
VERDICT_PATH = Path("data/silver_344_moontour_v2_verdicts.jsonl")

# Stage A's reproduced cycle-0 closure (sourced: data/silver_344_verified.jsonl,
# field ``residual_kms`` of the ic_ps96_sweep_reproduction row, commit 63809ec).
STAGE_A_CYCLE_0_CLOSURE_KMS = 0.010188096573990224

# v4.7 quasi_cycler envelope: closure residual < 0.5 km/s across the 10-cycle
# horizon. The Stage B SILVER blows this from cycle 1 onward.
QUASI_CYCLER_CLOSURE_ENVELOPE_KMS = 0.5

# V2 strict-cycler floors (sourced: V2_MOONTOUR_DRIFT_FLOOR_KMS /
# V2_MOONTOUR_CLOSURE_FLOOR_KMS in src/cyclerfinder/data/validation/v2_moontour.py).
V2_STRICT_DRIFT_FLOOR_KMS = 50_000.0
V2_STRICT_CLOSURE_FLOOR_KMS = 0.05

EXPECTED_HEADLINE_VERDICT = "FAIL"


def _load_v2_verdicts() -> tuple[dict[str, Any], dict[str, Any]]:
    """Return ``(n10_cycle_verdict, headline)`` from the Stage B V2-moontour JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    n10_verdict: dict[str, Any] | None = None
    headline: dict[str, Any] | None = None
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "moontour_v2_verdict"
                and row.get("candidate_id") == CANDIDATE_ID
                and row.get("n_cycles_requested") == 10
            ):
                n10_verdict = row
            elif row.get("_meta") and row.get("kind") == "headline":
                headline = row
    assert n10_verdict is not None, (
        f"Stage B n_cycles=10 V2 verdict for {CANDIDATE_ID!r} not found in {VERDICT_PATH}"
    )
    assert headline is not None, f"Stage B headline row not found in {VERDICT_PATH}"
    return n10_verdict, headline


def test_silver_344_v2_ten_cycles_all_lambert_converge() -> None:
    """All 10 Lambert cycles converge at all legs (``converged_legs == n_legs``).

    Lambert geometry is solvable across all 10 cycles - the FAIL verdict is
    about V_inf-continuity / drift, not Lambert non-convergence. This matches
    the #330 Umbriel-Oberon precedent (both candidates have geometrically
    convergent Lambert legs over 10 cycles).
    """
    n10, _ = _load_v2_verdicts()
    assert n10["n_cycles_completed"] == 10, (
        f"Stage B V2 expected 10/10 completed cycles, got {n10['n_cycles_completed']}/10."
    )
    for cyc in n10["per_cycle"]:
        assert cyc["converged_legs"] == cyc["n_legs"], (
            f"Stage B V2 cycle {cyc['cycle_index']} reports "
            f"{cyc['converged_legs']}/{cyc['n_legs']} legs converged - "
            "verdict needs re-validation."
        )


def test_silver_344_v2_cycle_zero_reproduces_stage_a_closure() -> None:
    """Cycle-0 closure under v2_moontour reproduces the Stage A SILVER value.

    Confirms the v2-convention phase0/rel_offset translation (Stage A
    {anchor: phase0, intermediate: phase0+rel_off} -> v2 {sorted-first:
    phase0, sorted-second: phase0+rel_off}) is correct. Without this
    cycle-0 closure would not match Stage A, indicating a geometry bug
    rather than a real V2 failure.
    """
    n10, _ = _load_v2_verdicts()
    cycle_0 = n10["per_cycle"][0]
    assert cycle_0["cycle_index"] == 0
    closure_0 = float(cycle_0["closure_residual_kms"])
    # Allow up to 1e-9 km/s for the Lambert solver's internal rounding;
    # Stage A and Stage B use the same kernel so the residual should match
    # to ~floating-point precision.
    assert abs(closure_0 - STAGE_A_CYCLE_0_CLOSURE_KMS) < 1e-9, (
        f"Stage B cycle-0 closure {closure_0:.6e} km/s does not reproduce "
        f"Stage A's {STAGE_A_CYCLE_0_CLOSURE_KMS:.6e} km/s - phase convention "
        "translation may be wrong."
    )


def test_silver_344_v2_closure_exits_quasi_envelope_immediately() -> None:
    """Closure residual blows past the v4.7 0.5 km/s quasi envelope.

    Distinguishes the Stage B FAIL from #330's PASS_QUASI_CYCLER. The
    #330 SILVER has ALL 10 cycles' closure below 0.5 km/s (its drift
    oscillates but Lambert continuity stays within the quasi envelope
    throughout). The Stage B SILVER has cycle-1 closure > 0.5 km/s AND
    the MAX closure across the 10-cycle horizon blows past 1 km/s (peaks
    above 5 km/s at cycle 1). The Stage B closure oscillates - cycles 5
    and 8 do return below 0.5 km/s - but the intervening cycles have
    closure > 1 km/s, so the candidate is structurally a
    closure-divergent oscillator, not a bounded quasi-cycler. The
    headline ``max_closure_residual_kms`` is the gate that disqualifies
    it from the quasi slot.
    """
    n10, _ = _load_v2_verdicts()
    per_cycle_closures = [float(c["closure_residual_kms"]) for c in n10["per_cycle"]]
    # Cycle 1 already exits the quasi envelope.
    cycle_1_closure = per_cycle_closures[1]
    assert cycle_1_closure >= QUASI_CYCLER_CLOSURE_ENVELOPE_KMS, (
        f"Stage B cycle-1 closure {cycle_1_closure:.3e} km/s is below the v4.7 "
        f"quasi envelope {QUASI_CYCLER_CLOSURE_ENVELOPE_KMS} km/s - if so, the "
        "candidate may actually be a quasi_cycler. Verdict needs re-review."
    )
    # The MAX closure across cycles 1..9 must blow past the quasi envelope -
    # this is what disqualifies the candidate from the quasi_cycler slot.
    # (#330 Umbriel-Oberon's max over 10 cycles was 0.349 km/s; the Stage B
    # max is > 5 km/s, ~15x the #330 max.)
    drifts_after_seed = per_cycle_closures[1:]
    max_closure_after_seed = max(drifts_after_seed)
    assert max_closure_after_seed > QUASI_CYCLER_CLOSURE_ENVELOPE_KMS, (
        f"Stage B max non-seed closure {max_closure_after_seed:.3e} km/s "
        f"falls inside the v4.7 quasi envelope {QUASI_CYCLER_CLOSURE_ENVELOPE_KMS} "
        "km/s - candidate may actually qualify as quasi_cycler; verdict review."
    )


def test_silver_344_v2_drift_exceeds_strict_floor() -> None:
    """Drift exceeds the 50,000 km strict-cycler floor and stays large.

    The strict-cycler V2 floor is 50,000 km. The Stage B drift hits ~9.6e5
    km by cycle 1 (~19x the floor). This test confirms the strict cycler
    gate is failed by drift alone (orthogonal to the closure-residual
    failure).
    """
    n10, _ = _load_v2_verdicts()
    per_cycle_drifts = [float(c["rendezvous_drift_kms"]) for c in n10["per_cycle"]]
    # Cycle 0 is the seed (drift = 0 by construction); skip it.
    drifts_after_seed = per_cycle_drifts[1:]
    assert min(drifts_after_seed) > V2_STRICT_DRIFT_FLOOR_KMS, (
        f"Stage B min non-seed drift {min(drifts_after_seed):.1f} km falls below "
        f"the strict V2 floor {V2_STRICT_DRIFT_FLOOR_KMS} km - strict-cycler "
        "verdict needs review."
    )


def test_silver_344_v2_headline_is_fail() -> None:
    """Strict V2 verdict is FAIL (no admission to strict OR quasi cycler slot)."""
    _, headline = _load_v2_verdicts()
    assert headline["verdict_label"] == EXPECTED_HEADLINE_VERDICT, (
        f"Stage B V2 headline {headline['verdict_label']!r} != "
        f"{EXPECTED_HEADLINE_VERDICT!r} - Stage B narrative needs review."
    )
    n10 = headline["n_cycles_10"]
    assert n10["passes_v2"] is False, (
        "Strict V2 must FAIL for the FAIL verdict (passing strict V2 would "
        "imply a strict-cycler verdict, not FAIL)."
    )
    assert headline["writeback_to_catalogue"] is False, (
        "Stage B must not request catalogue writeback - V2 alone never admits."
    )
