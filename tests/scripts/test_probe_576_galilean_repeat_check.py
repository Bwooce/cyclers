"""#576 -- regression coverage for the Galilean multi-cycle repeat-instrumentation
control, mirroring ``tests/scripts/test_probe_575_repeat_check.py``'s convention.

Pins the two load-bearing behaviors the #576 dispatch relies on:

1. Every one of the 36 gate-passing symmetric closures in
   ``data/enumerate_576_jupiter_galilean_symmetric_closures.jsonl`` (step 3's
   already-committed output) repeats to machine precision under
   :func:`run_v2_moontour` -- proving the construction is genuinely periodic,
   not a one-off closure.
2. A deliberately non-symmetric point (same pair/n_rev/commensurate-tof as a
   genuine Ganymede-Callisto survivor, but rel_offset=90 deg -- outside the
   {0,180} symmetric set the construction enumerates over) does NOT repeat,
   proving the instrumentation discriminates.

Both sides trace to this dispatch's own already-committed construction output
or its own construction machinery -- not a value this test's own code invented.
"""

from __future__ import annotations

import scripts.probe_576_galilean_repeat_check as probe576


def test_all_36_survivors_repeat_to_machine_precision() -> None:
    survivors = probe576.load_576_survivors()
    assert len(survivors) == 36  # pinned: 36 gate-passing symmetric closures, 6 pairs

    for rec in survivors:
        result = probe576.repeat_check_survivor(rec)
        assert result["n_cycles_completed"] == probe576.N_CYCLES, rec
        assert result["max_closure_residual_kms"] < 1e-6, rec
        assert result["repeats_to_machine_precision"] is True, rec


def test_non_symmetric_rel90_negative_control_does_not_repeat() -> None:
    survivors = probe576.load_576_survivors()
    neg_pt, negative_result = probe576.build_negative_control(survivors)

    # This exact point is Lambert-feasible but does not close (residual gate fails
    # at cycle 0 already -- it is not a symmetric construction).
    assert neg_pt is not None
    assert neg_pt["residual_kms"] >= 0.05  # GATE_RESIDUAL_KMS
    assert negative_result["repeats_to_machine_precision"] is False
