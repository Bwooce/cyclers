"""#655 -- regression coverage for the Rhea-Titan inclination-extension probe.

Pins the load-bearing facts the #655 dispatch report relies on:

1. The smoke test (``iapetus_state_3d`` reduces to ``_moon_state`` at inc=0,
   reused verbatim from #572) still passes when the module-level globals are
   overridden to Rhea/Titan -- this is the ONLY guard that a rotation-algebra
   bug in the reused #572 machinery hasn't silently broken.
2. Exactly 3 Rhea-anchored candidates are loaded from the committed #655
   enumeration output (never hand-transcribed).
3. All 3 candidates find at least one gate-passing closing basin at
   Rhea-Titan's real (conservative-upper-bound 0.7 deg) mutual inclination --
   this is the numeric result the #655 OUTSTANDING.md bullet reports.

Expected values trace to the already-committed
``data/enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl`` (produced by
the genericized, golden-validated #563 script), not values this test's own
code invented -- not a circular golden.
"""

from __future__ import annotations

import scripts.probe_655_rhea_titan_3d_closure as probe655


def test_smoke_test_reduction_passes_for_rhea_titan_globals() -> None:
    assert probe655.p572.ANCHOR == "Rhea"
    assert probe655.p572.FLYBY == "Titan"
    assert probe655.p572._smoke_test_reduction() is True


def test_loads_exactly_3_rhea_anchored_candidates() -> None:
    candidates = probe655.load_655_candidates()
    assert len(candidates) == 3
    labels = {c["label"] for c in candidates}
    assert labels == {
        "n3_nrev[0, 0]_rel0",
        "n5_nrev[1, 1]_rel180",
        "n8_nrev[2, 2]_rel0",
    }
    for c in candidates:
        assert c["coplanar_residual_kms"] < 1e-9  # machine-precision coplanar closure


def test_all_3_candidates_find_a_closing_basin_at_real_inclination() -> None:
    """Pins the #655 dispatch's key intermediate result: unlike Titan-Iapetus's
    15.5deg inclination, Rhea-Titan's ~0.7deg conservative-bound inclination is
    small enough that a single-cycle closing basin (residual + physical-bend
    gate) survives for every candidate. This does NOT by itself mean these are
    genuine repeating cyclers -- see test_probe_655_rhea_titan_repeat_check.py
    for the multi-cycle repeat check, which is where they fail."""
    candidates = probe655.load_655_candidates()
    for cand in candidates:
        sweep = probe655.p572.sweep_node_alignment(cand, n_omega=360)
        closing = [
            b
            for b in sweep["basins"]
            if b["residual_kms"] < probe655.GATE_RESIDUAL_KMS
            and probe655.candidate_passes_physical_gate(
                (probe655.p572.ANCHOR, probe655.p572.FLYBY, probe655.p572.ANCHOR),
                tuple(b["vinf_kms"]),
                min_useful_bend_deg=probe655.DEFAULT_MIN_USEFUL_BEND_DEG,
            )[0]
        ]
        assert closing, f"expected >=1 closing basin for {cand['label']}"
