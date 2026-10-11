# ruff: noqa: E501  (the row template is YAML with long lines)
"""#1025: catalogue rows for the 10 Earth-Venus members with Standish >= 4/5 (owner decision
2026-10-11, option "(a) 10 rows, Standish >= 4/5": candidate-novel V1, one row each).

For each member: the spec §14 V1 evidence (lamberthub izzo2015 + gooding1990 on every Lambert leg of
the #973 ideal-model zero; the #973 gauntlet's independent DOP853 re-fly), the ideal-model values,
the #1025 rung (d) results (Standish per epoch, the DE440 blend per epoch, the #415 band, the DE440
shoot) and the literature-gate result. Writes data/1025_ev_writeback_v1.json and the draft rows
(YAML text) to --rows-out; the catalogue edit itself is made in the catalogue window.

Usage: uv run python scripts/writeback_1025_ev.py --rows-out <path>
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.lambert import lambert_crosscheck
from cyclerfinder.search.two_working_body import LambertLeg, ResonantLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
AU_KM = 149_597_870.7
RE = REPO / "data" / "1025_ev_realeph"
#: Owner decision 2026-10-11: the members with Standish >= 4/5 (docs/notes/2026-10-08-1025-owner-summary.md sec. 7).
ROWS = ["ev4-1", "ev4-4", "ev4-6", "ev4-7", "ev5-4", "ev5-5", "ev5-6", "ev5-7", "ev5-8", "ev5-9"]
#: Members sharing one transfer and Venus block (identical real-ephemeris dates at every epoch;
#: docs/notes/2026-10-11-1025-ev-realeph-rung-d.md sec. 8).
SIBLINGS = {
    "ev4-4": ["ev4-6", "ev4-7"],
    "ev4-6": ["ev4-4", "ev4-7"],
    "ev4-7": ["ev4-4", "ev4-6"],
    "ev5-4": ["ev5-5"],
    "ev5-5": ["ev5-4"],
    "ev5-6": ["ev5-7"],
    "ev5-7": ["ev5-6"],
}
NAME = {"E": "Earth", "V": "Venus"}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")


def gauntlet_candidate(member: str) -> dict[str, Any]:
    k, i = member[2], int(member.split("-")[1])
    g = json.loads((REPO / "data" / "973_ev" / f"k{k}_gauntlet.json").read_text())
    c: dict[str, Any] = g["candidates"][i]
    assert c["screen_status"] == "pass", (member, c["screen_status"])
    return c


def shoot_best_residual(member: str) -> float | None:
    best = math.inf
    for p in (RE / "de440_shoot" / member).glob("*/e*/realeph_chain.json"):
        for r in json.loads(p.read_text())[0].get("shoot_restarts") or []:
            v = r.get("max_residual")
            if v is not None and math.isfinite(v):
                best = min(best, float(v))
    return None if best == math.inf else best


def evidence(
    member: str, summary: dict[str, Any], shoot: dict[str, Any], lit: dict[str, Any]
) -> dict[str, Any]:
    c = gauntlet_candidate(member)
    circ, a, b = ENUM.cell_system("ev")
    _, cyc = ENUM.parse_cycle_key(c["key"], circ, a, b)
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
    by_body: dict[str, dict[str, Any]] = {}
    for f in c["flyby_table"]:
        d = by_body.setdefault(f["body"], {"turns": [], "alt": math.inf, "ratio": 0.0})
        d["turns"].append(round(f["turn_deg"], 2))
        d["alt"] = min(d["alt"], f["required_alt_km"])
        d["ratio"] = max(d["ratio"], f["ratio"])
    s = summary[member]
    std = []
    for e in s["standish"]:
        ok = [r for r in e["routes"] if r["outcome"] == e["outcome"]]
        r = ok[0] if ok else (e["routes"][-1] if e["routes"] else {})
        std.append(
            {
                "epoch": e["epoch"],
                "epoch_jd": r.get("epoch_jd"),
                "outcome": e["outcome"],
                "route": r.get("route"),
                "worst_ratio": r.get("worst_ratio"),
                "refly_miss_km": r.get("refly_miss_km"),
            }
        )
    deb = [
        {
            "epoch": e["epoch"],
            "epoch_jd": e.get("epoch_jd"),
            "reached": e.get("reached"),
            "minimax_gate": e.get("minimax_gate"),
            "per7_total_mps": e.get("per7_total_mps"),
            "per7_venus_plus_deficit_mps": e.get("per7_venus_plus_deficit_mps"),
            "closure_dv_earth_mps": e.get("closure_dv_earth_mps"),
            "closure_dv_venus_mps": e.get("closure_dv_venus_mps"),
        }
        for e in s["de440_blend"]["per_epoch"]
    ]
    return {
        "member": member,
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
        "screen_status": c["screen_status"],
        "standish": std,
        "standish_counts": s["standish_counts"],
        "standish_band_best": s["standish_band"]["best"],
        "de440_blend": deb,
        "de440_blend_reached": s["de440_blend"]["reached"],
        "de440_blend_gate_pass": s["de440_blend"]["minimax_gate_pass"],
        "de440_blend_best": s["de440_blend"]["best"],
        "de440_blend_range_per7_mps": s["de440_blend"]["range_per7_total_mps"],
        "de440_shoot": [{"epoch": e["epoch"], "outcome": e["outcome"]} for e in shoot[member]],
        "de440_shoot_best_residual": shoot_best_residual(member),
        "literature_status": lit[member]["status"],
        "siblings": SIBLINGS.get(member, []),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows-out", type=Path, required=True)
    args = ap.parse_args()
    summary = json.loads((RE / "summary.json").read_text())
    shoot = json.loads((RE / "de440_shoot" / "summary.json").read_text())
    lit = {
        m["name"]: m
        for m in json.loads((REPO / "data" / "1025_litcheck_973_members.json").read_text())[
            "members"
        ]
    }
    ev = {m: evidence(m, summary, shoot, lit) for m in ROWS}
    for m, e in ev.items():
        c = e["standish_counts"]
        print(
            m,
            e["key"],
            f"lamberthub max {max(e['lamberthub_max_diff_mps_per_leg']):.1e} m/s",
            f"re-fly {e['dop853_refly_max_miss_km']:.1e} km",
            e["gate_on_integrated_vectors"],
            f"Standish {c['PASS']}/5",
            e["literature_status"],
            flush=True,
        )
        assert max(e["lamberthub_max_diff_mps_per_leg"]) < 1e-3
        assert e["dop853_refly_max_miss_km"] < 1.0 and e["gate_on_integrated_vectors"] == "pass"
        assert c["PASS"] >= 4, (m, c)
        assert all(x["outcome"] != "PASS" for x in e["de440_shoot"]), m
    (REPO / "data" / "1025_ev_writeback_v1.json").write_text(json.dumps(ev, indent=1))
    args.rows_out.write_text("".join(row_yaml(e) for e in ev.values()))
    print("wrote data/1025_ev_writeback_v1.json and", args.rows_out, flush=True)


def _std_text(e: dict[str, Any]) -> str:
    out = []
    for x in e["standish"]:
        if x["outcome"] == "PASS":
            out.append(
                f"JD {x['epoch_jd']:.1f} PASS (route {x['route']}, worst {x['worst_ratio']:.3f}, re-fly {x['refly_miss_km']:.1e} km)"
            )
        else:
            out.append(f"epoch {x['epoch']} {x['outcome']}")
    return "; ".join(out)


def _deb_text(e: dict[str, Any]) -> str:
    out = []
    for x in e["de440_blend"]:
        if not x["reached"]:
            out.append(f"epoch {x['epoch']} not reached")
            continue
        g = x["minimax_gate"]
        out.append(
            f"JD {x['epoch_jd']:.1f} gate {g['status']} {g['worst_ratio']:.3f}, closure E/V {x['closure_dv_earth_mps']:.1f}/{x['closure_dv_venus_mps']:.1f} m/s per chain"
        )
    return "; ".join(out)


def row_yaml(e: dict[str, Any]) -> str:
    m = e["member"]
    k = e["k"]
    n_cyc = 3 if k == 4 else 2
    ve, vv = e["vinf_kms"]["E"], e["vinf_kms"]["V"]
    be, bv = e["by_body"]["E"], e["by_body"]["V"]
    seq = "-".join(e["sequence"])
    years = e["period_days"] / 365.25
    lam = max(e["lamberthub_max_diff_mps_per_leg"])
    c = e["standish_counts"]
    sb = e["standish_band_best"]
    db = e["de440_blend_best"]
    rng = e["de440_blend_range_per7_mps"]
    shoot_counts = {
        o: sum(x["outcome"] == o for x in e["de440_shoot"])
        for o in ("PASS", "FAIL", "NOT REACHED", "NOT RUN")
    }
    res = e["de440_shoot_best_residual"]
    sib = (
        f"Shares its transfer and Venus block with {' and '.join(e['siblings'])}: identical ideal-model V_inf and identical real-ephemeris dates at every epoch; the members differ only in their Earth full-rev returns (docs/notes/2026-10-11-1025-ev-realeph-rung-d.md sec. 8)."
        if e["siblings"]
        else "No sibling among the 22 #973 ev members shares its transfer and Venus block."
    )
    vdv = (
        [
            """  - authors: ["VanderVeen, A. A."]
    year: 1969
    title: "The 1978 Venus-Swingby-To-Mercury Mission"
    venue: "Bellcomm Technical Memorandum TM-69-1013-2 (NASA-CR-103637, NTRS 19690023106)"
    doi:
    note: "PERIOD lineage only: names the Earth-Venus 8-yr (five-synodic-period) near repeat; not a cycler."
"""
        ]
        if k == 5
        else []
    )
    return f'''
- id: earth-venus-two-working-body-cycler-{m}-2026
  name: "Computed Earth-Venus two-working-body ballistic cycler {m} (k = {k}, both planets bend; #973/#1025)"
  source: this-project
  our_status: candidate-novel   # OWNER DECISION 2026-10-11 (#1025), option "(a) 10 rows, Standish >= 4/5": candidate-novel V1; literature_check v3-A1 returns not-found. docs/notes/2026-10-08-1025-owner-summary.md secs. 5-7.
  trajectory_regime: ballistic
  model_assumption: circular-coplanar
  cycler_class: multi-arc
  invariants:   # schema v4 (2026-06-03); cycle-level identity for multi-arc cyclers; see docs/spec.md §16.7.4
    aphelion_ratio: null   # not applicable: no Earth-Mars outbound/inbound arc pair
    transit_times_days: {e["transit_times_days"]}   # DERIVED (this project): per-leg durations in cycle order (key {e["key"]})
    turn_ratio: null   # not applicable; the demanded/available turn ratios are in vinf_kms_at_encounters notes
  orbit_class: cycler
  epoch_locked: false
  n_returns: infinite
  validation_level: V1   # spec §14 V1 (#1025 writeback): lamberthub Izzo+Gooding agreement and the #973 independent DOP853 re-fly; see validate._LEVEL_EVIDENCE. NOT V2.
  orbit_fidelity: circular-coplanar
  vinf_fidelity: circular-coplanar
  delta_v_kms:
  v_infinity_leveraging_dv_kms:
  fleet_size:
  flyby_mechanics:
  - body: "V"
    turning_angle_deg: {max(bv["turns"]):.2f}
    min_altitude_km: {bv["alt"]:.0f}
  - body: "E"
    turning_angle_deg: {max(be["turns"]):.2f}
    min_altitude_km: {be["alt"]:.0f}
  primary: "Sun"
  bodies: ["E", "V"]
  sequence_canonical: "{seq}"
  sense: "n/a"
  period:
    pair: "E-V"
    k: {k}
    years: {years:.3f}
    note: "DERIVED (this project): {k} Earth-Venus synodic periods = {e["period_days"]:.3f} d in the ideal model (Venus period 0.61520 yr, Earth 1 yr)."
  vinf_kms_at_encounters:
  - body: "E"
    vinf_kms: {ve:.3f}
    note: "DERIVED (this project, #973 ideal-model zero; data/973_ev/k{k}_gauntlet.json). Earth turns {" / ".join(f"{t:.2f}" for t in be["turns"])} deg; largest demanded/available ratio {be["ratio"]:.3f} at the 200-km floor."
  - body: "V"
    vinf_kms: {vv:.3f}
    note: "DERIVED (this project). Venus turns {" / ".join(f"{t:.2f}" for t in bv["turns"])} deg; largest ratio {bv["ratio"]:.3f} at the 300-km floor (required altitude {bv["alt"]:,.0f} km)."
  orbit_elements:
    a_au:
    e:
    perihelion_au:
    aphelion_au:
    inclination_deg: 0.0
    note: "Multi-arc. Distance from the Sun {e["r_min_km"] / AU_KM:.3f}-{e["r_max_km"] / AU_KM:.3f} AU (#973 extents, all legs, the full-revs included)."
  source_ephemeris: "ideal circular-coplanar model; real-ephemeris patched-conic rung (d) on Standish J2000 mean elements and on DE440 (blend to lambda 1, closure dV; the DE440 shoot)"
  data_gaps:
  - path: "trajectory_regime"
    kind: "uncertain"
    note: "BAND AT THE BEST EPOCH, WINDOW-DEPENDENT. Standish: ballistic and gate-passing at {c["PASS"]} of 5 epochs ({c["FAIL"]} fail, {c["NOT REACHED"]} not reached), turn deficit 0 per 7 cycles at the passing epochs (strictly ballistic by construction). DE440 blend (no shoot; #942 6.18/6.24 reading): reaches lambda 1 at {e["de440_blend_reached"]}/5 with the minimax gate passing at {e["de440_blend_gate_pass"]}/5; mid-course closure dV per 7 cycles at the best epoch {db["per7_total_mps"]:.1f} m/s ({db["band_total_upper_bound"]}; the Earth legs are an UPPER BOUND, DE440's Earth carries the lunar reflex), of which the Venus legs plus turn deficit {db["per7_venus_plus_deficit_mps"]:.1f} m/s ({db["band_venus_plus_deficit"]}); range over the reached epochs {rng[0]:.1f}-{rng[1]:.1f} m/s."
    source_hint: "docs/notes/2026-10-11-1025-ev-realeph-rung-d.md secs. 6, 5a; data/1025_ev_realeph/summary.json; data/1025_ev_writeback_v1.json"
    todo_ref: "#1025"
  - path: "validation_level"
    kind: "unknown"
    note: "DE440 ballistic closure untested in effect: the DE440 shoot found no gate-passing closure for this member ({shoot_counts["PASS"]} pass, {shoot_counts["FAIL"]} fail, {shoot_counts["NOT REACHED"]} no closure of 5 epochs; best shoot residual {res:.4f}, the worst junction component, where 0.001 means a 1-m/s V_inf mismatch or a 1-km miss). The shoot route cannot judge a full-rev Earth-Venus chain on DE440: its control ev-A fails with both the free and the anchored shoot (#1057 steps 1 and 1b). The negative is conditional on the formulation (restart defect #1058, date guard #1059; backlog, low priority, owner 2026-10-11)."
    source_hint: "docs/notes/2026-10-11-1025-ev-realeph-rung-d.md secs. 8, 9.3, 9.5; data/1025_ev_realeph/de440_shoot/"
    todo_ref: "#1057"
  legs: []
  first_published:
  corroborating_sources:
  - authors: ["Hollister, W. M."]
    year: 1969
    title: "Periodic Orbits for Interplanetary Flight"
    venue: "Journal of Spacecraft and Rockets 6(4):366-369"
    doi: "10.2514/3.29664"
    note: "CLASS source: the Earth-Venus two-working-body periodic orbits (orbits I-III)."
  - authors: ["Hollister, W. M.", "Menning, M. D."]
    year: 1970
    title: "Periodic Swing-By Orbits between Earth and Venus"
    venue: "Journal of Spacecraft and Rockets 7(10):1193-1199"
    doi: "10.2514/3.30134"
    note: "CLASS and METHOD source: orbits built from full-revolution (FR, 1:1) and symmetric (SY) returns only (Menning 1968 App. A-1). {m} carries n:m full-revolution returns other than 1:1, outside those itineraries."
{"".join(vdv)}  discovery_run:
    finder: "#973 two-working-body k-extension (scripts/run_973_enumerate.py; the #942/#943 generator)"
    discovery_date: "2026-10-08"
    description: "Exact zero of the Hollister-Menning date residual in the Earth-Venus cell at k = {k} (Earth returns full-rev only), gate-passing at both planets, clean under the #973 amendment-2.6/2.7 screens"
    task_chain: "#942 -> #971 -> #973 -> #1025"
    evidence_files: ["data/973_ev/k{k}_gauntlet.json", "data/1025_ev_writeback_v1.json", "data/1025_ev_realeph/summary.json", "data/1025_ev_realeph/de440_shoot/summary.json", "data/1025_litcheck_973_members.json"]
    finder_version:
    run_id: "#973 k{k} ev cell (lead launch)"
    cell_id: "ev {e["key"]}"
  priority_date: "2026-10-08"
  notes: |
    Member {m} of the #973 Earth-Venus two-working-body set: both planets bend; k = {k}
    ({e["period_days"]:.3f} d, {years:.3f} yr); structure {e["key"]}. Ideal-model Lambert-leg
    starts {", ".join(f"{v:.6f}" for v in e["x_days"])} d (model epoch: both planets at angle 0 at
    t = 0). Worst demanded/available turn ratio {e["worst_ratio"]:.3f}.

    Family: the same two-working-body Earth-Venus class as ev-A
    (earth-venus-two-working-body-cycler-eva-2026) and ev-B (earth-venus-two-working-body-cycler-evb-2026):
    both planets bend, with returns outside the Hollister-Menning FR/SY itineraries. ev-C
    (earth-venus-venus-hosted-cycler-evc-2026) is a different architecture (Earth massless). {sib}
    #1027 (direction-dependent tilted-circle flag) does not apply: the #973 screen status is "pass".

    V1 evidence (data/1025_ev_writeback_v1.json): lamberthub <= {lam:.1e} m/s on every Lambert leg;
    #973 DOP853 re-fly {e["dop853_refly_max_miss_km"]:.1e} km, gate passes on the integrated vectors.

    Real-ephemeris standing (#1025 rung (d), docs/notes/2026-10-11-1025-ev-realeph-rung-d.md; chains of
    {n_cyc} cycles = {n_cyc * years:.1f} yr; route controls Hollister 1H, ev-A and ev-B reproduce):
    - Standish J2000 mean elements: {c["PASS"]}/5 PASS. {_std_text(e)}.
      Band at the best epoch: {sb["band"]} (turn deficit {sb["per7_beyond_bend_mps"]:.1e} m/s per 7 cycles, numerical zero).
    - DE440 blend (lambda 1, minimax directions, no shoot): {_deb_text(e)}.
    - DE440 shoot (#1025 sec. 8): no gate-passing closure; the route cannot judge (#1057), see data_gaps.

    Literature: literature_check v3-A1 (#1025) returns "{e["literature_status"]}" (necessary, not
    sufficient; data/1025_litcheck_973_members.json). The Hollister/H&M anchor is excluded on return
    types. Human prior-art reading: docs/notes/2026-10-08-1025-literature-step.md sec. 5 (the ev-A/ev-B
    search, docs/notes/2026-10-07-942-evAB-prior-art-search.md, applies; no new search was run).
    OWNER DECISION 2026-10-11: candidate-novel V1 under option "(a) 10 rows, Standish >= 4/5",
    attributed to Hollister 1969 and Hollister & Menning 1970. Wording: "first computed", never
    "new class".
'''


if __name__ == "__main__":
    main()
