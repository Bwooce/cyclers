"""#655 -- regression coverage for the Rhea-Titan multi-cycle repeat check.

Pins the #655 dispatch's headline negative result: the 3 single-cycle closing
basins found by ``probe_655_rhea_titan_3d_closure.py`` at Rhea-Titan's real
(conservative-bound 0.7 deg) mutual inclination do NOT survive as genuine
repeating 3-cycle cyclers -- the same qualitative failure #575 found for
Titan-Iapetus at 15.5 deg, now reproduced at a mutual inclination ~22x
smaller. Loaded directly from the already-committed
``data/probe_655_rhea_titan_3d_closure.jsonl`` (never hand-transcribed).
"""

from __future__ import annotations

import scripts.probe_655_rhea_titan_repeat_check as repeat655


def test_loads_exactly_3_closing_basins() -> None:
    basins = repeat655.load_closing_basins()
    assert len(basins) == 3
    labels = {b["label"] for b in basins}
    assert labels == {
        "n3_nrev[0, 0]_rel0",
        "n5_nrev[1, 1]_rel180",
        "n8_nrev[2, 2]_rel0",
    }


def test_none_of_the_3_basins_repeat_as_genuine_cycles() -> None:
    """The headline #655 negative: 0/3 -- Lambert converges every cycle (no
    infeasibility), but the V_inf-continuity residual blows past the project
    gate (0.05 km/s) to 0.58-1.16 km/s and inter-cycle drift reaches
    3e4-1e6 km, matching #575's Titan-Iapetus failure signature even though
    the mutual inclination here is far smaller."""
    import json

    enum_by_label = {}
    enum_path = repeat655.DATA_DIR / "enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl"
    with enum_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "pass" and d.get("anchor") == "Rhea":
                lbl = f"n{d['n_commensurate_int']}_nrev{d['n_rev']}_rel{d['rel_offset_deg']:.0f}"
                enum_by_label[lbl] = d

    basins = repeat655.load_closing_basins()
    assert len(basins) == 3
    import re

    n_repeat = 0
    for b in basins:
        m = re.match(r"n(\d+)_nrev\[(\d+), (\d+)\]_rel(\d+)", b["label"])
        assert m is not None
        n_rev = (int(m.group(2)), int(m.group(3)))
        rel_offset_deg = enum_by_label[b["label"]]["rel_offset_deg"]
        res = repeat655.repeat_check(b, rel_offset_deg, n_rev)
        assert res["n_cycles_completed"] == 3  # Lambert converges every cycle
        assert res["max_residual_kms"] > repeat655.GATE_RESIDUAL_KMS  # but blows the gate
        if res["repeats_as_genuine_cycle"]:
            n_repeat += 1
    assert n_repeat == 0
