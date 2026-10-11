# ruff: noqa: E501  (the row template is YAML with long lines)
"""#1055: DRAFT catalogue row for gc4-19 (held by the owner until GanCal#1 on jup365 reports).

Reuses the #1025 writeback (``scripts/writeback_1025_gc.py``, unchanged) for the spec §14 V1
evidence and the row template, with the full-mass values from the #1055 run
(``data/1055_gc4-19/gc4-19.json``, reached from the corrected seed) and the #1055 caveats. Writes
``data/1055_gc4-19/v1_evidence.json`` and the draft YAML to --rows-out. No catalogue edit.

Usage: uv run python scripts/draft_1055_gc4_19_row.py --rows-out data/1055_gc4-19/gc4-19_row_draft.yaml
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.lambert import lambert_crosscheck
from cyclerfinder.search.two_working_body import LambertLeg, ResonantLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
M = "gc4-19"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


W = _load("writeback_1025_gc", REPO / "scripts" / "writeback_1025_gc.py")


def evidence() -> dict[str, Any]:
    """As writeback_1025_gc.evidence, with the full-mass record from #1055."""
    _, key, vg, vc = next(m for m in W.BATCH.MEMBERS if m[0] == M)
    _, c = W.BATCH.gauntlet_candidate(key, vg, vc)
    circ, a, b = W.ENUM.cell_system("gc")
    _, cyc = W.ENUM.parse_cycle_key(c["key"], circ, a, b)
    x = np.asarray(c["x_days"]) * DAY
    legs = eval_lambert_legs(circ, cyc, x)
    assert legs is not None
    diffs = []
    for ev, li in zip(legs, cyc.lambert_index, strict=True):
        leg = cyc.legs[li]
        assert isinstance(leg, LambertLeg)
        r1, _ = circ.state(leg.frm, ev.t_dep)
        r2, _ = circ.state(leg.to, ev.t_arr)
        res = lambert_crosscheck(
            r1, r2, ev.t_arr - ev.t_dep, mu=circ.mu, n_revs=leg.nrev, branch=leg.branch
        )
        diffs.append(float(res["max_diff_mps"]))
    tt, li = [], 0
    for leg in cyc.legs:
        if isinstance(leg, LambertLeg):
            tt.append(round((legs[li].t_arr - legs[li].t_dep) / DAY, 3))
            li += 1
        else:
            assert isinstance(leg, ResonantLeg)
            tt.append(round(leg.body_revs * circ.period_s(leg.body) / DAY, 3))
    fm = json.loads((REPO / "data" / "1055_gc4-19" / f"{M}.json").read_text())
    assert fm["outcome"] == "EXISTS", fm["outcome"]
    by_body: dict[str, dict[str, float]] = {}
    for f in c["flyby_table"]:
        d = by_body.setdefault(f["body"], {"turn": 0.0, "alt": float("inf"), "ratio": 0.0})
        d["turn"] = max(d["turn"], f["turn_deg"])
        d["alt"] = min(d["alt"], f["required_alt_km"])
        d["ratio"] = max(d["ratio"], f["ratio"])
    return {
        "member": M,
        "key": c["key"],
        "x_days": c["x_days"],
        "period_days": c["period_days"],
        "k": int(c["key"].split("|")[0][1:]),
        "vinf_kms": c["vinf_kms"],
        "by_body": by_body,
        "worst_ratio": c["worst_ratio"],
        "sequence": [f["body"] for f in sorted(c["flyby_table"], key=lambda f: f["t_days"])],
        "r_min_km": c["r_min_km"],
        "r_max_km": c["r_max_km"],
        "transit_times_days": tt,
        "lamberthub_max_diff_mps_per_leg": diffs,
        "dop853_refly_max_miss_km": c["cross_check"]["max_arrival_miss_km"],
        "dop853_max_vinf_vector_error_kms": c["cross_check"]["max_vinf_vector_error_kms"],
        "gate_on_integrated_vectors": c["cross_check"]["gate_status_integrated"],
        "full_mass_nodes": [
            {k: n[k] for k in ("moon", "t_days", "vinf_kms", "rp_km", "turn_deg")}
            for n in fm["full_mass"]
        ],
        "ias15_max_km": fm["ias15"]["max_dr_km"],
    }


REPLACE = [
    (
        "(k = 4, both moons bend; #973/#1025)",
        "(k = 4, both moons bend; #973/#1025/#1055)",
    ),
    (
        "  our_status: candidate-novel   # OWNER DECISION 2026-10-08 (#1025): a member of gc-1's two-working-body Ganymede-Callisto class, written candidate-novel with gc-1's attribution; the gate result is recorded in notes. docs/notes/2026-10-08-1025-owner-summary.md.",
        "  our_status: candidate-novel   # DRAFT, OWNER DECISION PENDING (held 2026-10-11 until GanCal#1 on jup365 reports): proposed as for the 15 #1025 members (a member of gc-1's class, gc-1's attribution).",
    ),
    (
        "continuous-gravity ideal-model rung (#1025)",
        "continuous-gravity ideal-model rung (#1025 rules, reached from the #1055 corrected seed)",
    ),
    (
        '    todo_ref: "#1025"\n',
        '    todo_ref: "#1025"\n'
        '  - path: "notes"\n'
        '    kind: "method"\n'
        '    note: "Continuous-gravity closure (ideal model) reached from the #1055 degenerate-family-corrected seed: the raw seed does not start (#1025 NUMERICAL STOP at sigma 0.02). The method passed its pre-registered controls A (gc-1 ideal) and B (GanEur#316 jup365), and the negative control holds (4/4 off-root starts fail). The branch\'s sigma -> 0 limit is this #973 member (Richardson, within 14 s; #973 gate re-check, fd5ee948)."\n'
        '    source_hint: "docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md secs. 8, 10, 10.1; docs/notes/2026-10-11-1055-gc4-19-973-gate.md"\n'
        '    todo_ref: "#1055"\n',
    ),
    ('    task_chain: "#943 -> #973 -> #1025"', '    task_chain: "#943 -> #973 -> #1025 -> #1055"'),
    (
        '["data/973_gc/k4_gauntlet.json", "data/1025_sigma/gc4-19.json", "data/1025_writeback_v1.json", "data/1025_litcheck_973_members.json"]',
        '["data/973_gc/k4_gauntlet.json", "data/1055_gc4-19/gc4-19.json", "data/1055_gc4-19/fix.json", "data/1055_gc4-19/neg.json", "data/1055_gc4-19/limit.json", "data/1055_gc4-19/v1_evidence.json", "data/1025_litcheck_973_members.json"]',
    ),
    (
        "    V1 evidence (data/1025_writeback_v1.json)",
        "    V1 evidence (data/1055_gc4-19/v1_evidence.json)",
    ),
    (
        "    Continuous gravity (#1025, docs/notes/2026-10-08-1025-sigma-batch.md sec. 10; ideal model with\n"
        "    point-mass moons at real GMs, not the ephemeris): EXISTS.",
        "    Continuous gravity (#1025 rules, docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md\n"
        "    sec. 10; ideal model with point-mass moons at real GMs, not the ephemeris): EXISTS, reached\n"
        "    from the #1055 degenerate-family-corrected seed (the raw seed does not start, #1025).",
    ),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows-out", type=Path, required=True)
    args = ap.parse_args()
    lit = {
        m["name"]: m
        for m in json.loads((REPO / "data" / "1025_litcheck_973_members.json").read_text())[
            "members"
        ]
    }
    e = evidence()
    e["literature_status"] = lit[M]["status"]
    assert max(e["lamberthub_max_diff_mps_per_leg"]) < 1e-3
    assert e["dop853_refly_max_miss_km"] < 1.0 and e["gate_on_integrated_vectors"] == "pass"
    assert e["literature_status"] == "inconclusive", e["literature_status"]
    (REPO / "data" / "1055_gc4-19" / "v1_evidence.json").write_text(json.dumps(e, indent=1))
    text = W.row_yaml(e)
    for a, b in REPLACE:
        assert text.count(a) == 1, a[:60]
        text = text.replace(a, b)
    args.rows_out.write_text(text)
    print(
        M,
        f"lamberthub max {max(e['lamberthub_max_diff_mps_per_leg']):.1e} m/s",
        f"re-fly {e['dop853_refly_max_miss_km']:.1e} km",
        e["gate_on_integrated_vectors"],
        e["literature_status"],
        "->",
        args.rows_out,
    )


if __name__ == "__main__":
    main()
