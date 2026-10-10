# ruff: noqa: E501  (the row template is YAML with long lines)
"""#1025: catalogue rows for the 15 Ganymede-Callisto members that EXIST at full mass in the
continuous ideal model (owner decision 2026-10-08 evening: candidate-novel V1, one row each).

For each member: the spec §14 V1 evidence (lamberthub izzo2015 + gooding1990 on every Lambert leg of
the ideal-model zero; the #973 gauntlet's independent DOP853 re-fly), the patched-conic values, the
#1025 full-mass values and the literature-gate result. Writes data/1025_writeback_v1.json and the
draft rows (YAML text) to --rows-out; the catalogue edit itself is made in the catalogue window.

Usage: uv run python scripts/writeback_1025_gc.py --rows-out <path>
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
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.two_working_body import LambertLeg, ResonantLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
EXISTS = [
    "gc4-1", "gc4-2", "gc4-6", "gc4-8", "gc4-10", "gc4-25", "gc4-28", "gc4-38", "gc4-41",
    "gc5-3", "gc5-5", "gc5-6",
    "gc6-0", "gc6-1", "gc6-2",
]  # fmt: skip
SHORT = {"Ganymede": "G", "Callisto": "C"}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


BATCH = _load("run_1025_sigma_batch", REPO / "scripts" / "run_1025_sigma_batch.py")
ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")


def evidence(member: str) -> dict[str, Any]:
    _, key, vg, vc = next(m for m in BATCH.MEMBERS if m[0] == member)
    _, c = BATCH.gauntlet_candidate(key, vg, vc)
    circ, a, b = ENUM.cell_system("gc")
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
    # per-leg durations in cycle order (Lambert from the solve; full-rev = body_revs body periods)
    tt, li = [], 0
    for leg in cyc.legs:
        if isinstance(leg, LambertLeg):
            tt.append(round((legs[li].t_arr - legs[li].t_dep) / DAY, 3))
            li += 1
        else:
            assert isinstance(leg, ResonantLeg)
            tt.append(round(leg.body_revs * circ.period_s(leg.body) / DAY, 3))
    fm = json.loads((REPO / "data" / "1025_sigma" / f"{member}.json").read_text())
    assert fm["outcome"] == "EXISTS", (member, fm["outcome"])
    by_body: dict[str, dict[str, float]] = {}
    for f in c["flyby_table"]:
        d = by_body.setdefault(f["body"], {"turn": 0.0, "alt": float("inf"), "ratio": 0.0})
        d["turn"] = max(d["turn"], f["turn_deg"])
        d["alt"] = min(d["alt"], f["required_alt_km"])
        d["ratio"] = max(d["ratio"], f["ratio"])
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
        "full_mass_nodes": [
            {k: n[k] for k in ("moon", "t_days", "vinf_kms", "rp_km", "turn_deg")}
            for n in fm["full_mass"]
        ],
        "ias15_max_km": fm["ias15"]["max_dr_km"],
    }


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
    ev = {m: evidence(m) for m in EXISTS}
    for m, e in ev.items():
        e["literature_status"] = lit[m]["status"]
        print(
            m,
            e["key"],
            f"lamberthub max {max(e['lamberthub_max_diff_mps_per_leg']):.1e} m/s",
            f"re-fly {e['dop853_refly_max_miss_km']:.1e} km",
            e["gate_on_integrated_vectors"],
            e["literature_status"],
            flush=True,
        )
        assert max(e["lamberthub_max_diff_mps_per_leg"]) < 1e-3
        assert e["dop853_refly_max_miss_km"] < 1.0 and e["gate_on_integrated_vectors"] == "pass"
    (REPO / "data" / "1025_writeback_v1.json").write_text(json.dumps(ev, indent=1))
    args.rows_out.write_text("".join(row_yaml(e) for e in ev.values()))
    print("wrote data/1025_writeback_v1.json and", args.rows_out, flush=True)


def _fm_text(e: dict[str, Any]) -> str:
    parts = []
    for n in e["full_mass_nodes"]:
        r = SATELLITES[n["moon"]].radius_eq_km
        parts.append(
            f"{SHORT[n['moon']]} {n['vinf_kms']:.3f} km/s, altitude {n['rp_km'] - r:,.0f} km, "
            f"turn {n['turn_deg']:.1f} deg"
        )
    return "; ".join(parts)


def row_yaml(e: dict[str, Any]) -> str:
    m = e["member"]
    k = e["k"]
    g, c = e["vinf_kms"]["Ganymede"], e["vinf_kms"]["Callisto"]
    bg, bc = e["by_body"]["Ganymede"], e["by_body"]["Callisto"]
    seq = "-".join(e["sequence"])
    if m == "gc6-0":
        gate = (
            'literature_check v3-A1 (#1025): "published" via Campagnola et al. 2019 GCGC AS RUN; '
            "a signature-convention artefact (the cycle is G-G-C, not alternating; started at "
            'Callisto it reads "inconclusive" naming Russell-Strange 2009 Ganymede-Callisto, '
            "data/1025_litcheck_gc6-0_diagnostic.json; #1056)"
        )
    else:
        gate = (
            'literature_check v3-A1 (#1025): "inconclusive", naming the Russell-Strange 2009 '
            "Ganymede-Callisto ideal-model moon cycler (one working body; rule F7)"
        )
    lam = max(e["lamberthub_max_diff_mps_per_leg"])
    return f'''
- id: ganymede-callisto-two-working-body-cycler-{m}-2026
  name: "Computed Ganymede-Callisto two-working-body ballistic cycler {m} (k = {k}, both moons bend; #973/#1025)"
  source: this-project
  our_status: candidate-novel   # OWNER DECISION 2026-10-08 (#1025): a member of gc-1's two-working-body Ganymede-Callisto class, written candidate-novel with gc-1's attribution; the gate result is recorded in notes. docs/notes/2026-10-08-1025-owner-summary.md.
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
  validation_level: V1   # spec §14 V1 (#1025 writeback): lamberthub Izzo+Gooding agreement and the #973 independent DOP853 re-fly; see validate._LEVEL_EVIDENCE. NOT V2 (no real-ephemeris chain).
  orbit_fidelity: circular-coplanar
  vinf_fidelity: circular-coplanar
  delta_v_kms:
  v_infinity_leveraging_dv_kms:
  fleet_size:
  flyby_mechanics:
  - body: "Ganymede"
    turning_angle_deg: {bg["turn"]:.2f}
    min_altitude_km: {bg["alt"]:.0f}
  - body: "Callisto"
    turning_angle_deg: {bc["turn"]:.2f}
    min_altitude_km: {bc["alt"]:.0f}
  primary: "Jupiter"
  bodies: ["Ganymede", "Callisto"]
  sequence_canonical: "{seq}"
  sense: "n/a"
  period:
    pair: "Ganymede-Callisto"
    k: {k}
    years:
    note: "DERIVED (this project): {k} Ganymede-Callisto synodic periods = {e["period_days"]:.3f} d in the Russell & Strange 2009 Table 2 ideal model."
  vinf_kms_at_encounters:
  - body: "Ganymede"
    vinf_kms: {g:.3f}
    note: "DERIVED (this project, #973 ideal-model patched-conic zero; data/973_gc/k{k}_gauntlet.json). Largest demanded/available turn ratio {bg["ratio"]:.3f} at the 100-km floor."
  - body: "Callisto"
    vinf_kms: {c:.3f}
    note: "DERIVED (this project). Largest demanded/available turn ratio {bc["ratio"]:.3f} at the 200-km floor."
  orbit_elements:
    a_au:
    e:
    perihelion_au:
    aphelion_au:
    inclination_deg: 0.0
    note: "Multi-arc (no single ellipse). Distance from Jupiter {e["r_min_km"]:,.0f}-{e["r_max_km"]:,.0f} km (#973 extents, all legs)."
  source_ephemeris: "ideal circular-coplanar Galilean model (Russell & Strange 2009 Table 2 constants); continuous-gravity ideal-model rung (#1025)"
  data_gaps:
  - path: "vinf_kms_at_encounters"
    kind: "unknown"
    note: "Real-ephemeris standing untested. The patched-conic real-ephemeris chain (rung (d)) was not run for this member (#973 note), and real-ephemeris continuous gravity is untested (lane limit: the open-chain n-body formulation does not converge across a full-revolution return, #1044/#1046)."
    source_hint: "docs/notes/2026-10-07-973-two-working-body-k-extension.md; docs/notes/2026-10-08-1025-sigma-batch.md sec. 10"
    todo_ref: "#1025"
  legs: []
  first_published:
  corroborating_sources:
  - authors: ["Hollister, W. M.", "Menning, M. D."]
    year: 1970
    title: "Periodic Swing-By Orbits between Earth and Venus"
    venue: "Journal of Spacecraft and Rockets 7(10):1193-1199"
    doi: "10.2514/3.30134"
    note: "METHOD source: the two-working-body date-residual corrector, generalised here to the Jovian moons (as for gc-1)."
  - authors: ["Russell, R. P.", "Strange, N. J."]
    year: 2009
    title: "Cycler Trajectories in Planetary Moon Systems"
    venue: "Journal of Guidance, Control, and Dynamics 32(1):143-157"
    doi: "10.2514/1.36610"
    note: "MODEL source: the ideal circular-coplanar Galilean model; their cyclers have one working body and a massless target (AAS 07-118 p.2), and AAS 07-118 p.18 names removing that assumption as future work."
  - authors: ["Campagnola, S.", "Buffington, B. B.", "Lam, T.", "Petropoulos, A. E.", "Pellegrini, E."]
    year: 2019
    title: "Tour Design Techniques for the Europa Clipper Mission"
    venue: "Journal of Guidance, Control, and Dynamics 42(12):2615-2626"
    doi: "10.2514/1.G004309"
    note: "PRIOR ART: the one published two-working-body Ganymede-Callisto cycler (GCGC, a single trajectory, not a searched class); as for gc-1."
  discovery_run:
    finder: "#973 two-working-body k-extension (scripts/run_973_enumerate.py; the #942/#943 generator)"
    discovery_date: "2026-10-07"
    description: "Exact zero of the Hollister-Menning date residual in the Ganymede-Callisto cell at k = {k}, gate-passing at both moons, clean under the #973 amendment-2.6 screens"
    task_chain: "#943 -> #973 -> #1025"
    evidence_files: ["data/973_gc/k{k}_gauntlet.json", "data/1025_sigma/{m}.json", "data/1025_writeback_v1.json", "data/1025_litcheck_973_members.json"]
    finder_version:
    run_id: "#973 k{k} gc cell"
    cell_id: "gc {e["key"]}"
  priority_date: "2026-10-07"
  notes: |
    Member {m} of the two-working-body Ganymede-Callisto class of gc-1
    (ganymede-callisto-two-working-body-cycler-gc1-2026): both moons bend; k = {k}
    ({e["period_days"]:.3f} d); structure {e["key"]}. Ideal-model Lambert-leg starts
    {", ".join(f"{v:.6f}" for v in e["x_days"])} d (model epoch: all bodies at angle 0 at t = 0).
    Worst demanded/available turn ratio {e["worst_ratio"]:.3f}.

    V1 evidence (data/1025_writeback_v1.json): lamberthub <= {lam:.1e} m/s on every Lambert leg;
    #973 DOP853 re-fly {e["dop853_refly_max_miss_km"]:.1e} km, gate passes on the integrated vectors.

    Continuous gravity (#1025, docs/notes/2026-10-08-1025-sigma-batch.md sec. 10; ideal model with
    point-mass moons at real GMs, not the ephemeris): EXISTS. The joint-mass branch reaches the full
    masses (identity by the limit, IAS15 re-fly {e["ias15_max_km"]:.1e} km, no unscheduled minimum
    inside 3 Hill radii). Full-mass nodes: {_fm_text(e)}. Patched-conic values above.

    Real ephemeris: untested. Rung (d) (the patched-conic chain tool) was not run for this member,
    and real-ephemeris continuous gravity is untested (lane limit, #1044/#1046).

    Literature: {gate}. Resolving notes: as for gc-1 (docs/notes/2026-10-06-943-gc-prior-art-search.md)
    and docs/notes/2026-10-08-1025-literature-step.md. OWNER DECISION 2026-10-08: written
    candidate-novel as a member of gc-1's class, with gc-1's attribution (Hollister & Menning 1970
    method; Russell & Strange 2007/2009 model; Campagnola et al. 2019 nearest prior art). Wording:
    "first computed", never "new class".
'''


if __name__ == "__main__":
    main()
