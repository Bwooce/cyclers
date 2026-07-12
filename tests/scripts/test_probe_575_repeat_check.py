"""#575 C2 -- regression coverage for the two-sided repeat-instrumentation control.

Pins the two load-bearing behaviors the #575 dispatch report relies on:

1. The known-bad #571 branch-1 negative control (documented 117.6/78.0/38.3 deg
   per-cycle transfer-angle drift in ``run_574_stageB_saturn_gauntlet.py``) must NOT
   repeat under :func:`run_v2_saturn_3d` even at e=0 -- proving the instrumentation
   actually discriminates a one-off closure from a genuine periodic one.
2. The #571 cross-check scans exactly the 187 candidates the project's own
   ``data/scan_571_saturn_titan_pairs_index.jsonl`` records as
   ``n_all_gates_passed`` for the two Titan-Iapetus directions (69 + 118), per
   ``data/OUTSTANDING.md``'s own "187 Titan-Iapetus candidates" figure -- a count
   sourced from already-committed data, not invented by this dispatch.

Both expected values trace to already-committed data files (``data/probe_574_*``,
``data/scan_571_saturn_titan_pairs_index.jsonl``), not values this test's own code
computed -- not a circular golden.
"""

from __future__ import annotations

import scripts.probe_575_titan_iapetus_repeat_check as probe575


def test_branch1_negative_control_does_not_repeat() -> None:
    params = probe575.load_branch1_negative_control()
    assert params.n_rev == (1, 1)
    assert params.e_titan == 0.0
    assert params.e_iapetus == 0.0
    assert params.inclination_deg == 0.0

    result = probe575.repeat_check("branch1_e0_test", params)
    # Branch 1's own documented failure mode: the same-n_rev Lambert transfer
    # physically ceases to exist past cycle 0 (confirmed at e=0 too, per the
    # #574 Stage B write-up) -- it must not complete all N_CYCLES.
    assert result["n_cycles_completed"] < probe575.N_CYCLES
    assert result["repeats_to_machine_precision"] is False


def test_575_symmetric_survivors_all_repeat() -> None:
    _t_syn, sqrt_papb = probe575._t_syn_and_sqrt_papb()
    survivors = probe575.load_575_titan_anchored_survivors()
    assert len(survivors) == 9  # pinned: 9 Titan-anchored gate-passing symmetric closures

    for rec in survivors:
        params = probe575.params_from_575_record(rec, sqrt_papb)
        result = probe575.repeat_check(f"n{rec['n_commensurate_int']}_test", params)
        assert result["n_cycles_completed"] == probe575.N_CYCLES, rec
        assert result["max_closure_residual_kms"] < 1e-6, rec
        assert result["repeats_to_machine_precision"] is True, rec


def test_571_cross_check_count_matches_committed_187() -> None:
    t_syn, _sqrt_papb = probe575._t_syn_and_sqrt_papb()
    cross = probe575.cross_check_571_against_symmetric_condition(t_syn)
    # data/scan_571_saturn_titan_pairs_index.jsonl: n_all_gates_passed 69 (Titan
    # anchored) + 118 (Iapetus anchored) = 187, matching data/OUTSTANDING.md's own
    # "187 Titan-Iapetus candidates" figure (#552/#571).
    assert cross["n_571_records_scanned"] == 187
