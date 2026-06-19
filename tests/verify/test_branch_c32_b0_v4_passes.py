"""§14 V4 frozen-gate evidence for branch_C32_b0 — an HONEST-NEGATIVE verdict.

Unlike the SILVER's V4 (which PASSED and earned a catalogue row), branch_C32_b0
FAILS V4: under a real DE440 n-body model (Earth + Moon + Sun + Mars + Jupiter)
the far-amplitude (3,3) orbit is destabilized by real solar tides and escapes
the Earth-Moon neighborhood within a few periods. The candidate cleared V1-V3
(all of which live in or near the idealized, Sun-free CR3BP) but cannot clear V4.

This is a deliberate HALT verdict per `feedback_orbit_closure_discipline`: a
clean negative at V4 is a legitimate outcome. branch_C32_b0 is NOT admitted to
the catalogue. This frozen-gate test ties the recorded V4 evidence to that
HALT decision, so any future change that would FLIP the verdict (e.g. a model
fix that makes it pass) trips a loud test failure and forces a re-review before
admission.

What this gate asserts (against ``data/branch_c32_b0_v4_verdict.jsonl``):

* every V4 n_cycles record reports ``passes_v4 == False`` and the footer reports
  ``all_pass == False`` (the honest-negative);
* the full-model V4 drift is catastrophically large (>= 1e8 km — an escape, not
  a bounded cycle), confirming the failure mode is escape rather than a marginal
  miss;
* the Earth+Moon-only control drift stays bounded (< 1e6 km), so the failure is
  attributable to the real solar/planetary tides (the perturbers), NOT to a
  seeding/registration artifact;
* the structural diagnostic confirms the cause: the orbit amplitude is a large
  fraction of the Earth-Sun Hill radius (>= 0.5) and the solar tidal
  acceleration is a non-negligible fraction of Earth's gravity (>= 0.1) — the
  regime the CR3BP cannot model.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from P389.5. This test is the catalogue-protection ratchet for
the HALT decision.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

VERDICT_PATH = Path("data/branch_c32_b0_v4_verdict.jsonl")
CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"

# Failure-mode floors (the verdict must remain a clear escape, not a marginal miss).
V4_ESCAPE_FLOOR_KMS = 1.0e8  # full-model drift at/above this = escape
EM_ONLY_BOUNDED_CEIL_KMS = 1.0e6  # EM-only control must stay below this (bounded)
HILL_FRACTION_FLOOR = 0.5  # amplitude / Earth-Sun Hill radius
SOLAR_TIDE_RATIO_FLOOR = 0.1  # solar tide / Earth gravity


def _load_records() -> dict[str, Any]:
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    records: dict[str, Any] = {"v4": []}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            kind = row.get("kind")
            if kind == "v4_verdict_cr3bp_periodic":
                records["v4"].append(row)
            elif kind:
                records[kind] = row
    return records


def test_branch_c32_b0_v4_is_honest_negative() -> None:
    """Every V4 n_cycles record FAILS and the footer reports all_pass=False."""
    records = _load_records()
    assert records["v4"], f"no V4 verdict records in {VERDICT_PATH}"
    for row in records["v4"]:
        assert row["candidate_id"] == CANDIDATE_ID
        assert row["passes_v4"] is False, (
            f"V4 n_cycles={row['n_cycles_requested']} unexpectedly PASSES — "
            "the HALT verdict has flipped; re-review before any admission."
        )
    assert "footer" in records, f"footer missing from {VERDICT_PATH}"
    assert records["footer"]["all_pass"] is False, (
        "V4 footer all_pass is no longer False — the candidate would now admit; "
        "re-review the model + the HALT decision."
    )


def test_branch_c32_b0_v4_failure_mode_is_solar_tide_escape() -> None:
    """Full-model V4 escapes (huge drift) while the EM-only control stays bounded."""
    records = _load_records()
    max_full_drift = max(float(r["max_v4_drift_kms"]) for r in records["v4"])
    assert max_full_drift >= V4_ESCAPE_FLOOR_KMS, (
        f"full-model V4 max drift {max_full_drift:.3e} km < {V4_ESCAPE_FLOOR_KMS:.0e} km — "
        "no longer a clear escape; the failure mode has changed."
    )
    assert "control_em_only" in records, (
        f"Earth+Moon-only control record missing from {VERDICT_PATH}"
    )
    ctrl = records["control_em_only"]
    ctrl_max = float(ctrl["max_drift_kms"])
    assert ctrl_max < EM_ONLY_BOUNDED_CEIL_KMS, (
        f"EM-only control max drift {ctrl_max:.3e} km >= {EM_ONLY_BOUNDED_CEIL_KMS:.0e} km — "
        "the control no longer isolates the perturbers as the cause; the failure "
        "may be a seeding artifact and needs investigation."
    )
    # The full-model drift must dwarf the control to attribute the failure to tides.
    assert max_full_drift > 100.0 * ctrl_max, (
        "full-model V4 drift no longer dwarfs the EM-only control — the "
        "solar-tide attribution is broken."
    )


def test_branch_c32_b0_v4_structural_cause_is_near_hill_amplitude() -> None:
    """The orbit amplitude is a large fraction of the Earth-Sun Hill radius."""
    records = _load_records()
    assert "structural_diagnostic" in records, (
        f"structural_diagnostic record missing from {VERDICT_PATH}"
    )
    diag = records["structural_diagnostic"]
    hill_frac = float(diag["amplitude_fraction_of_hill"])
    tide_ratio = float(diag["solar_tide_to_earth_gravity_ratio"])
    assert hill_frac >= HILL_FRACTION_FLOOR, (
        f"orbit amplitude is only {hill_frac:.3f} of the Earth-Sun Hill radius "
        f"(< {HILL_FRACTION_FLOOR}) — the near-Hill structural cause no longer holds."
    )
    assert tide_ratio >= SOLAR_TIDE_RATIO_FLOOR, (
        f"solar tide / Earth gravity is only {tide_ratio:.3f} (< {SOLAR_TIDE_RATIO_FLOOR}) — "
        "the solar-tide-dominance structural cause no longer holds."
    )
