"""§14 V3 frozen-gate evidence for branch_C32_b0 — REBOUND IAS15 n-body cross-check.

This module is a FROZEN-GATE wrapper asserting that the #389 P389.4 V3 verdict
JSONL on disk says what the catalogue admission claims. The V3 question is
whether the V2 round-off-floor bounded-drift signature is a REAL dynamical
property of the Earth+Moon 2-body system or a DOP853 artifact. REBOUND IAS15
(an independent integrator, symplectic Gauss-Radau, epsilon=1e-12) propagates
the spacecraft in the inertial frame with Earth and Moon as point-mass
primaries on Kepler-consistent circular orbits matching the CR3BP setup,
the IAS15 result is converted back to the rotating frame and compared per
cycle vs the V2 driver's drift series. The agreement floor is spec §14's
100 km.

What this gate asserts (per #389 P389.4 across n_cycles in {3, 5, 10}):

* ``passes_v3 = True``: V3-vs-V2 drift-agreement headline (max |V3-V2| in km)
  is below the 100 km agreement floor.
* The integrator is ``REBOUND IAS15`` (the named independent cross-check; if
  the JSONL records a different integrator, the registry pointer is stale).
* All cycles converge at every n_cycles {3, 5, 10}.
* The agreement headline sits in the integrator-floor band (< 1 km),
  matching the essentially-stable Floquet character.

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from #389 P389.4. This test ties the catalogue row's V3 claim
to the recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "branch-c32-b0-em-3-3-quasi-cycler-2026"
VERDICT_PATH = Path("data/branch_c32_b0_v3_verdict.jsonl")
V3_AGREEMENT_FLOOR_KMS = 100.0
EXPECTED_INTEGRATOR = "REBOUND IAS15"
EXPECTED_N_CYCLES = (3, 5, 10)
TIGHT_INTEGRATOR_AGREEMENT_FLOOR_KMS = 1.0


def _load_v3_verdicts() -> dict[int, dict[str, Any]]:
    """Return ``{n_cycles: verdict_row}`` from the V3 JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    verdicts: dict[int, dict[str, Any]] = {}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if (
                row.get("kind") == "v3_verdict_cr3bp_periodic"
                and row.get("candidate_id") == CANDIDATE_ID
            ):
                verdicts[int(row["n_cycles_requested"])] = row
    for n in EXPECTED_N_CYCLES:
        assert n in verdicts, (
            f"#389 P389.4 V3 verdict at n_cycles={n} for {CANDIDATE_ID!r} "
            f"not found in {VERDICT_PATH}"
        )
    return verdicts


def test_branch_c32_b0_v3_passes_at_all_three_n_cycles() -> None:
    """V3 passes at n_cycles in {3, 5, 10} — IAS15 reproduces V2 drift below 100 km."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["passes_v3"] is True, (
            f"#389 V3 verdict at n_cycles={n} reports passes_v3={row['passes_v3']!r} — "
            "registry entry must be re-validated."
        )
        agreement = float(row["drift_agreement_kms"])
        assert agreement < V3_AGREEMENT_FLOOR_KMS, (
            f"#389 V3 drift_agreement at n_cycles={n} is {agreement:.3e} km — "
            f"NOT below the {V3_AGREEMENT_FLOOR_KMS:.0f} km spec §14 V3 floor."
        )


def test_branch_c32_b0_v3_uses_rebound_ias15_integrator() -> None:
    """The recorded integrator is REBOUND IAS15 (the named V3 cross-check)."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["integrator"] == EXPECTED_INTEGRATOR, (
            f"#389 V3 verdict at n_cycles={n} integrator={row['integrator']!r}, "
            f"expected {EXPECTED_INTEGRATOR!r} — registry pointer is stale."
        )


def test_branch_c32_b0_v3_all_cycles_converge() -> None:
    """Every cycle converges at every n_cycles {3, 5, 10}."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["converged_at_each_cycle"] is True, (
            f"#389 V3 verdict at n_cycles={n} reports "
            f"converged_at_each_cycle={row['converged_at_each_cycle']!r}; expected True."
        )


def test_branch_c32_b0_v3_agreement_in_integrator_floor_band() -> None:
    """The V3-vs-V2 agreement sits in the integrator-round-off band (< 1 km).

    branch_C32_b0 is essentially stable; both V2 (DOP853) and V3 (IAS15) drifts
    are integrator round-off, so their difference is too. An agreement creeping
    into the spec-floor band would indicate either a model regression or a
    discovery that the orbit is structurally different than Phase 2 reported.
    """
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        agreement = float(row["drift_agreement_kms"])
        assert agreement < TIGHT_INTEGRATOR_AGREEMENT_FLOOR_KMS, (
            f"#389 V3 drift_agreement at n_cycles={n} is {agreement:.3e} km — "
            f"exceeds the {TIGHT_INTEGRATOR_AGREEMENT_FLOOR_KMS:.1e} km "
            "integrator-floor band the essentially-stable Floquet character predicts."
        )
