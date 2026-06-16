"""§14 V3 frozen-gate evidence for the #327 SILVER (REBOUND IAS15 n-body).

This module is a FROZEN-GATE wrapper asserting that the #331 V3 n-body gauntlet
verdict JSONL on disk says what the catalogue row claims. The V3 question is
whether the V2 bounded-drift signature is a REAL dynamical property or merely
a Lambert/DOP853 integrator artifact. The independent integrator (REBOUND
IAS15, symplectic Gauss-Radau, eps=1e-12) reproduces the V2 driver's per-cycle
drift series to nanometer precision — the signature is REAL.

What this gate asserts (per #331 across n_cycles in {3, 5, 10}):

* ``passes_v3 = True``: the V3-vs-V2 drift-agreement headline (max |V3-V2| in
  km) is below the 100 km agreement floor.
* The integrator is ``REBOUND IAS15`` (the named independent cross-check; if
  the JSONL records a different integrator, the registry pointer is stale).
* All legs converge across all cycles (``converged_legs == n_legs``).

Sourced-golden discipline: the EXPECTED side is the JSONL file itself, frozen
project output from #331. This test ties the catalogue row's V3 claim to the
recorded evidence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "repeated-moon-uranus-00000041"
VERDICT_PATH = Path("data/silver_327_v3_verdicts.jsonl")
V3_AGREEMENT_FLOOR_KMS = 100.0  # spec §14 V3 V3-vs-V2 agreement floor
EXPECTED_INTEGRATOR = "REBOUND IAS15"
EXPECTED_N_CYCLES = (3, 5, 10)


def _load_v3_verdicts() -> dict[int, dict[str, Any]]:
    """Return ``{n_cycles: verdict_row}`` from the #331 V3 n-body JSONL."""
    assert VERDICT_PATH.exists(), f"frozen verdict file missing: {VERDICT_PATH}"
    verdicts: dict[int, dict[str, Any]] = {}
    with VERDICT_PATH.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "moontour_v3_verdict" and row.get("candidate_id") == CANDIDATE_ID:
                verdicts[int(row["n_cycles_propagated"])] = row
    for n in EXPECTED_N_CYCLES:
        assert n in verdicts, (
            f"#331 V3 verdict at n_cycles={n} for {CANDIDATE_ID!r} not found in {VERDICT_PATH}"
        )
    return verdicts


def test_silver_327_v3_passes_at_all_three_n_cycles() -> None:
    """V3 passes at n_cycles in {3, 5, 10} — bounded drift survives REBOUND IAS15."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["passes_v3"] is True, (
            f"#331 V3 verdict at n_cycles={n} reports passes_v3={row['passes_v3']!r} — "
            "registry entry must be re-validated."
        )
        agreement = float(row["drift_agreement_kms"])
        assert agreement < V3_AGREEMENT_FLOOR_KMS, (
            f"#331 V3 drift_agreement at n_cycles={n} is {agreement:.3e} km — "
            f"NOT below the {V3_AGREEMENT_FLOOR_KMS} km agreement floor."
        )


def test_silver_327_v3_uses_named_independent_integrator() -> None:
    """The recorded integrator is REBOUND IAS15 (the named V3 cross-check)."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        assert row["integrator"] == EXPECTED_INTEGRATOR, (
            f"#331 V3 verdict at n_cycles={n} integrator={row['integrator']!r}, "
            f"expected {EXPECTED_INTEGRATOR!r} — registry pointer is stale."
        )


def test_silver_327_v3_all_legs_converge_across_all_cycles() -> None:
    """Every leg converges in every cycle at every n_cycles {3, 5, 10}."""
    verdicts = _load_v3_verdicts()
    for n in EXPECTED_N_CYCLES:
        row = verdicts[n]
        for cyc in row["per_cycle"]:
            assert cyc["converged_legs"] == cyc["n_legs"], (
                f"#331 V3 n_cycles={n} cycle {cyc['cycle_index']} reports "
                f"{cyc['converged_legs']}/{cyc['n_legs']} legs converged."
            )
