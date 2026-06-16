"""§14 V2 frozen-gate evidence for the #327 SILVER (quasi-cycler verdict).

This module is a FROZEN-GATE wrapper asserting that the #330 V2-moontour
gauntlet verdict JSONL on disk says what the catalogue row claims. The strict
V2 verdict is ``FAIL_QUASI_BOUNDED`` (drift exceeds the 50,000 km strict-cycler
floor), BUT the bounded-oscillation pattern matches the v4.7 ``quasi_cycler``
structural definition exactly. THAT is the V2 evidence the row stands on.

What this gate asserts:

* The #330 verdict at ``data/silver_327_moontour_v2_verdicts.jsonl`` reports
  10/10 cycles convergent at all 10 Lambert legs (``converged_legs == n_legs``
  for every cycle).
* Per-cycle drift oscillates between ~86,000 km and ~530,000 km — peaks at
  cycle 1, returns near 86,000 km by cycle 5, then oscillates within that
  band — the bounded-drift signature of a near-5:1 synodic-resonance
  quasi-cycler.
* The strict V2 verdict label is ``FAIL_QUASI_BOUNDED``, NOT
  ``FAIL_DIVERGING`` (so the row enters the v4.7 ``quasi_cycler`` slot, not a
  rejected slot).

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from #330. This test ties the catalogue row's quasi_cycler
claim to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "repeated-moon-uranus-00000041"
VERDICT_PATH = Path("data/silver_327_moontour_v2_verdicts.jsonl")

# The V2 strict floor is 50,000 km cycler-class drift; the SILVER FAILS_V2
# precisely because it BOUNDS at ~86k-530k km rather than diverging. The
# catalogue's quasi_cycler slot encodes "bounded oscillation, NOT monotonic
# divergence" — these are the bounds we assert.
EXPECTED_DRIFT_LOWER_BOUND_KMS = 80_000.0  # cycle 5 oscillation floor ~86k
EXPECTED_DRIFT_UPPER_BOUND_KMS = 600_000.0  # cycle 7 oscillation peak ~530k
EXPECTED_HEADLINE_VERDICT = "FAIL_QUASI_BOUNDED"


def _load_v2_verdicts() -> tuple[dict[str, Any], dict[str, Any]]:
    """Return ``(n10_cycle_verdict, headline)`` from the #330 V2-moontour JSONL."""
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
        f"#330 n_cycles=10 V2 verdict for {CANDIDATE_ID!r} not found in {VERDICT_PATH}"
    )
    assert headline is not None, f"#330 headline row not found in {VERDICT_PATH}"
    return n10_verdict, headline


def test_silver_327_v2_ten_cycles_all_converge() -> None:
    """All 10 Lambert cycles converge at all legs (``converged_legs == n_legs``)."""
    n10, _ = _load_v2_verdicts()
    assert n10["n_cycles_completed"] == 10, (
        f"#330 V2 expected 10/10 completed cycles, got {n10['n_cycles_completed']}/10."
    )
    for cyc in n10["per_cycle"]:
        assert cyc["converged_legs"] == cyc["n_legs"], (
            f"#330 V2 cycle {cyc['cycle_index']} reports "
            f"{cyc['converged_legs']}/{cyc['n_legs']} legs converged — "
            "registry entry must be re-validated."
        )


def test_silver_327_v2_drift_is_bounded_oscillation() -> None:
    """Per-cycle drift oscillates inside the ~86k-530k km bounded band."""
    n10, _ = _load_v2_verdicts()
    per_cycle_drifts = [float(c["rendezvous_drift_kms"]) for c in n10["per_cycle"]]
    # Cycle 0 is the seed (drift = 0 by construction); skip it for the
    # bounded-oscillation check.
    drifts_after_seed = per_cycle_drifts[1:]
    max_drift = max(drifts_after_seed)
    assert max_drift < EXPECTED_DRIFT_UPPER_BOUND_KMS, (
        f"#330 V2 max drift {max_drift:.1f} km exceeds bounded-oscillation upper "
        f"bound {EXPECTED_DRIFT_UPPER_BOUND_KMS:.0f} km — the quasi_cycler "
        "characterization needs re-validation."
    )
    # Within the 10-cycle band, the trajectory must return CLOSE to the seed
    # (not monotonically diverge). The minimum non-seed drift is ~86,000 km at
    # cycle 5; we require it to be < EXPECTED_DRIFT_LOWER_BOUND_KMS upper-side
    # check to confirm at least one return to that band.
    min_drift = min(drifts_after_seed)
    assert min_drift < EXPECTED_DRIFT_LOWER_BOUND_KMS + 10_000.0, (
        f"#330 V2 min non-seed drift {min_drift:.1f} km does not return to the "
        f"~{EXPECTED_DRIFT_LOWER_BOUND_KMS:.0f} km band — quasi_cycler bounded "
        "oscillation signature unverified."
    )


def test_silver_327_v2_headline_is_quasi_bounded_fail() -> None:
    """Strict V2 verdict is FAIL_QUASI_BOUNDED (quasi_cycler slot, not rejected)."""
    _, headline = _load_v2_verdicts()
    assert headline["verdict_label"] == EXPECTED_HEADLINE_VERDICT, (
        f"#330 V2 headline {headline['verdict_label']!r} != "
        f"{EXPECTED_HEADLINE_VERDICT!r} — quasi_cycler admission needs review."
    )
    n10 = headline["n_cycles_10"]
    assert n10["passes_v2"] is False, (
        "Strict V2 must FAIL for the quasi_cycler slot (passing strict V2 would "
        "imply strict-periodic cycler, not quasi)."
    )
