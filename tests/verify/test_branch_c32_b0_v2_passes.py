"""§14 V2 frozen-gate evidence for branch_C32_b0 — bounded-cycle CR3BP gate.

This module is a FROZEN-GATE wrapper asserting that the #389 P389.3 V2 verdict
JSONL on disk says what the catalogue admission claims. The V2 question for
branch_C32_b0 — a closed CR3BP periodic orbit, NOT a Lambert-tour — is
whether the orbit remains bounded under DOP853 propagation for n consecutive
cycles WITHOUT recorrecting. The spec §14 same-model 50,000 km floor applies.

What this gate asserts (per #389 P389.3 across n_cycles in {3, 5, 10}):

* ``passes_v2 = True`` at every n_cycles
* ``max_drift_kms`` is far below the 50,000 km floor (the orbit is
  essentially stable per its Phase 2 max_floquet_mag = 1.000000000000617,
  sigma_d/day = 6.08e-15 — the per-cycle drift is integrator round-off, NOT
  dynamical instability).
* ``converged_at_each_return = True`` (propagator never failed mid-flight).

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from #389 P389.3. This test ties the catalogue claim to the
recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
VERDICT_PATH = Path("data/branch_c32_b0_v2_verdict.jsonl")
V2_DRIFT_FLOOR_KMS = 50_000.0  # spec §14 same-model floor
EXPECTED_N_CYCLES = (3, 5, 10)
# branch_C32_b0 is essentially-stable (max_floquet_mag = 1.000000000000617).
# The expected drift across 10 cycles is dominated by DOP853 round-off, not
# dynamics — we assert it stays under 1 km (~5 orders of magnitude below the
# spec floor) to catch any future integrator regression or model change.
TIGHT_INTEGRATOR_DRIFT_FLOOR_KMS = 1.0


def _load_v2_verdicts() -> dict[int, dict[str, Any]]:
    """Return ``{n_cycles: verdict_row}`` from the V2 JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    verdicts: dict[int, dict[str, Any]] = {}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v2_verdict_cr3bp_periodic"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                verdicts[int(row["n_cycles_requested"])] = row
    for n in EXPECTED_N_CYCLES:
        assert n in verdicts, (
            f"#389 P389.3 V2 verdict at n_cycles={n} for {CANDIDATE_ID!r} "
            f"not found in {VERDICT_PATH}"
        )
    return verdicts


def test_branch_c32_b0_v2_passes_at_all_three_n_cycles() -> None:
    """V2 passes at n_cycles in {3, 5, 10} — bounded drift below the spec floor."""
    verdicts = _load_v2_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["passes_v2"] is True, (
            f"#389 V2 verdict at n_cycles={n} reports passes_v2={row['passes_v2']!r} — "
            "registry entry must be re-validated."
        )
        max_drift = float(row["max_drift_kms"])
        assert max_drift < V2_DRIFT_FLOOR_KMS, (
            f"#389 V2 max_drift at n_cycles={n} is {max_drift:.3e} km — "
            f"NOT below the {V2_DRIFT_FLOOR_KMS:.0f} km spec §14 floor."
        )
        assert row["converged_at_each_return"] is True, (
            f"#389 V2 verdict at n_cycles={n} reports converged_at_each_return="
            f"{row['converged_at_each_return']!r}; expected True."
        )


def test_branch_c32_b0_v2_drift_is_integrator_floor_not_dynamics() -> None:
    """The drift across 10 cycles is dominated by integrator round-off, not dynamics.

    branch_C32_b0's Phase 2 max_floquet_mag is 1.000000000000617 — essentially
    unity to 13 decimals. The per-cycle drift over 10 cycles must stay at
    integrator round-off levels (we assert < 1 km, ~5 orders below the spec
    floor). A drift growing into the spec-floor band would indicate either an
    integrator regression or a discovery that the orbit is not as stable as
    Phase 2 reported.
    """
    verdicts = _load_v2_verdicts()
    row = verdicts[10]
    max_drift = float(row["max_drift_kms"])
    assert max_drift < TIGHT_INTEGRATOR_DRIFT_FLOOR_KMS, (
        f"#389 V2 max_drift at n_cycles=10 is {max_drift:.3e} km — exceeds the "
        f"{TIGHT_INTEGRATOR_DRIFT_FLOOR_KMS:.1e} km integrator-floor band the "
        "essentially-stable Floquet character predicts. Either the Phase 2 "
        "Floquet was wrong or the integrator regressed."
    )


def test_branch_c32_b0_v2_drift_is_below_spec_floor_by_many_orders() -> None:
    """The V2 drift sits ~9 orders of magnitude below the 50,000 km spec floor."""
    verdicts = _load_v2_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        max_drift = float(row["max_drift_kms"])
        # ~5 orders of magnitude below the floor is the structural margin we
        # expect for an essentially-stable orbit; ~9 is what we actually see.
        assert max_drift < V2_DRIFT_FLOOR_KMS / 1.0e5, (
            f"#389 V2 max_drift at n_cycles={n} is {max_drift:.3e} km — "
            f"loss of margin against the {V2_DRIFT_FLOOR_KMS:.0f} km spec floor."
        )


def test_branch_c32_b0_v2_drift_floor_matches_spec() -> None:
    """The recorded drift_floor_kms matches the spec §14 same-model floor."""
    verdicts = _load_v2_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert float(row["drift_floor_kms"]) == V2_DRIFT_FLOOR_KMS, (
            f"#389 V2 verdict at n_cycles={n} drift_floor_kms="
            f"{row['drift_floor_kms']} != {V2_DRIFT_FLOOR_KMS} — "
            "registry entry must be re-validated."
        )
