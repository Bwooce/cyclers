"""#942/#943 gauntlet: adjudicate the gate-passing cyclers of one enumeration cell.

Pre-registered in docs/notes/2026-10-05-942-943-two-working-body-generator.md sec. 6.1.
For every physically distinct gate-passing cycler (analyse_942_enumeration.collect):

1. Independent cross-check: every leg is re-flown with scipy DOP853 (rtol 1e-13,
   atol 1e-8 km) on the two-body equations, not with the Lambert solver or the Kepler
   step. Lambert legs start
   from the solved departure V_inf; fixed (full-rev, half-rev) legs start from the chosen
   flyby direction. Reported: the largest arrival miss (km) and V_inf vector error (km/s),
   the largest junction magnitude mismatch built from the integrated arrivals, and the turn
   gate re-run on the integrated vectors.
2. SOI self-consistency: largest miss / smallest Laplace sphere of influence.
3. Literal collisions: the Russell & Strange 2007/2009 Table 3 rows (same body pair, V_inf at
   both bodies within 0.05 km/s and period within 0.3 d = LITERAL; within 0.3 km/s = NEAR);
   Campagnola 2019 GCGC (Ganymede-Callisto, V_inf 3.5/4.5 km/s, within 0.5 = NEAR); Hollister's
   circular-coplanar orbits I-III (Earth-Venus, k = 2, same topology = PUBLISHED FAMILY);
   catalogue rows on the same body pair (id, our_status, V_inf listed for comparison);
   empty_regions.jsonl entries naming both bodies.
4. literature_check.check_literature with the offline corpus backend (KNOWN_CORPUS, which
   now holds the Rall 1969/1971, Pisarevsky 2008, R-S 2007 and Campagnola 2019 anchors).
   This is not a web search.

Nothing here calls a cycler novel. Output: <out>.json and a summary table on stdout.

Usage: uv run python scripts/gauntlet_942.py --cell vm DIR [DIR ...] --out FILE
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]
from scipy.integrate import solve_ivp

from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)
from cyclerfinder.search.two_working_body import (
    Flyby,
    _blocks,
    cycle_flybys,
    eval_lambert_legs,
    fixed_duration_s,
    gate_cycle,
    sphere_of_influence_km,
)

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
ANALYSE = _load("analyse_942_enumeration", REPO / "scripts" / "analyse_942_enumeration.py")

#: Russell & Strange 2007 (AAS 07-118) Table 3 (p.9), via the #960 digest:
#: (id, flyby body A, target B, V_inf A, V_inf B, period d, number of legs). The number of
#: legs (Table 3) is the number of flyby-body encounters per cycle; the target is met once.
RS_ROWS = [
    ("VenMar#45", "V", "M", 8.22, 12.96, 667.8, 1),
    ("EurGan#93", "Europa", "Ganymede", 2.37, 4.10, 28.2, 3),
    ("EurGan#131", "Europa", "Ganymede", 2.40, 4.10, 21.2, 2),
    ("EurGan#159", "Europa", "Ganymede", 2.45, 4.11, 28.2, 3),
    ("GanCal#1", "Ganymede", "Callisto", 3.18, 3.26, 37.6, 3),
    ("GanCal#5", "Ganymede", "Callisto", 3.24, 3.34, 37.6, 2),
    ("GanEur#5", "Ganymede", "Europa", 1.66, 2.57, 35.3, 1),
    ("GanEur#43", "Ganymede", "Europa", 1.87, 3.89, 14.1, 1),
    ("GanEur#316", "Ganymede", "Europa", 3.20, 3.81, 49.4, 4),
]
#: Hollister 1969 p.367 circular-coplanar orbits I-III = Menning 1H-3H (k = 2, E-V).
HOLLISTER_TOPOLOGIES = {
    "orbit I (1H)": ("RE/1:1", "LE>V/0s", "RV/1:1", "RV/1:1", "LV>E/0s"),
    "orbit II (2H)": ("RE/1:1", "LE>V/0s", "RV/1:1", "LV>V/1*", "LV>E/0s"),
    "orbit III (3H)": ("LE>E/1*", "LE>V/0s", "RV/1:1", "RV/1:1", "LV>E/0s"),
}
CATALOGUE_CODE = {"Ganymede": "Ganymede", "Callisto": "Callisto", "Europa": "Europa"}


def _twobody(_t: float, y: np.ndarray, mu: float) -> np.ndarray:
    r = y[:3]
    return np.concatenate([y[3:], -mu * r / float(np.linalg.norm(r)) ** 3])


def fly(mu: float, r0: np.ndarray, v0: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    sol = solve_ivp(
        _twobody,
        (0.0, dt),
        np.concatenate([r0, v0]),
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-8,
    )
    return sol.y[:3, -1], sol.y[3:, -1]


def cross_check(system: Any, cycle: Any, x: np.ndarray, flybys: list[Flyby]) -> dict[str, float]:
    """Re-fly every leg with DOP853 and rebuild the junctions and the gate."""
    legs = eval_lambert_legs(system, cycle, x)
    assert legs is not None
    worst_miss = 0.0
    worst_dv = 0.0
    integrated_in: dict[float, np.ndarray] = {}
    for ev, li in zip(legs, cycle.lambert_index, strict=True):
        leg = cycle.legs[li]
        r0, w0 = system.state(leg.frm, ev.t_dep)
        r1, v1 = fly(system.mu, r0, w0 + ev.vinf_dep, ev.t_arr - ev.t_dep)
        rb, wb = system.state(leg.to, ev.t_arr)
        worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
        worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - ev.vinf_arr)))
        integrated_in[round(ev.t_arr, 3)] = v1 - wb
    # fixed legs: start from each chosen direction, compare with the next flyby's inbound
    for blk in _blocks(system, cycle, legs):
        for leg, t0 in blk.fixed:
            fb_out = next(f for f in flybys if abs(f.t_s - t0) < 1.0 and f.body == blk.body)
            dt = fixed_duration_s(system, leg)
            r0, w0 = system.state(blk.body, t0)
            r1, v1 = fly(system.mu, r0, w0 + fb_out.vinf_out, dt)
            rb, wb = system.state(blk.body, t0 + dt)
            worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
            nxt = next(f for f in flybys if abs(f.t_s - (t0 + dt)) < 1.0 and f.body == blk.body)
            worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - nxt.vinf_in)))
            integrated_in[round(t0 + dt, 3)] = v1 - wb
    # gate on integrated inbound vectors where available
    rebuilt = []
    for f in flybys:
        vin = integrated_in.get(round(f.t_s, 3))
        if vin is None:
            # the wrap arrival is stored at t + period
            vin = integrated_in.get(round(f.t_s + cycle.period_s, 3), f.vinf_in)
        rebuilt.append(Flyby(f.body, f.t_s, np.asarray(vin), f.vinf_out))
    rep = gate_cycle(system, rebuilt)
    mism = max(
        (
            abs(float(np.linalg.norm(f.vinf_in)) - float(np.linalg.norm(f.vinf_out)))
            for f in rebuilt
        ),
        default=0.0,
    )
    return {
        "max_arrival_miss_km": worst_miss,
        "max_vinf_vector_error_kms": worst_dv,
        "max_junction_mismatch_kms": mism,
        "gate_status_integrated": rep.status,
        "worst_ratio_integrated": rep.gate.worst_ratio,
    }


def collisions(cell: str, line: dict, system: Any, a: str, b: str) -> list[str]:
    out = []
    v = line["vinf_kms"]
    period_d = line["k"] * system.synodic_s(a, b) / DAY
    n_enc = {c: sum(1 for f in line["flybys"] if f[0] == c) for c in (a, b)}
    for rid, fa, fb, va, vb, per, n_legs in RS_ROWS:
        if {fa, fb} != {a, b} or fa not in v or fb not in v:
            continue
        dv = max(abs(v[fa] - va), abs(v[fb] - vb))
        # V_inf and period alone do not separate the members of one R-S family (the EurGan
        # rows share their V_inf), so LITERAL also needs R-S's encounter counts.
        # (a massless target is not in the flyby sequence)
        same_count = n_enc[fa] == n_legs and (n_enc[fb] == 1 or system.body(fb).massless)
        if dv <= 0.05 and abs(period_d - per) <= 0.3 and same_count:
            out.append(f"LITERAL {rid} (dV_inf {dv:.3f}, period {period_d:.1f} vs {per})")
        elif dv <= 0.3:
            out.append(
                f"NEAR {rid} (dV_inf {dv:.3f}, period {period_d:.1f} vs {per}, "
                f"encounters {n_enc[fa]}+{n_enc[fb]} vs {n_legs}+1)"
            )
    if {a, b} == {"Ganymede", "Callisto"}:
        dv = max(abs(v.get("Ganymede", 0) - 3.5), abs(v.get("Callisto", 0) - 4.5))
        dv2 = max(abs(v.get("Ganymede", 0) - 4.5), abs(v.get("Callisto", 0) - 3.5))
        if min(dv, dv2) <= 0.5:
            out.append(f"NEAR Campagnola 2019 GCGC (V_inf 3.5/4.5; dV {min(dv, dv2):.2f})")
    if cell == "ev" and line["k"] == 2:
        for key in line["keys"]:
            toks = tuple(t.replace("/1l", "/1*").replace("/1h", "/1*") for t in key.split("|")[1:])
            for name, topo in HOLLISTER_TOPOLOGIES.items():
                rots = {toks[i:] + toks[:i] for i in range(len(toks))}
                if topo in rots:
                    out.append(f"PUBLISHED FAMILY Hollister 1969 {name}")
    return sorted(set(out))


def catalogue_rows(a: str, b: str) -> list[dict]:
    cat = yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())
    rows = cat if isinstance(cat, list) else cat.get("cyclers", [])
    out = []
    for r in rows:
        bodies = set(r.get("bodies") or [])
        if {a, b} <= bodies or {CATALOGUE_CODE.get(a, a), CATALOGUE_CODE.get(b, b)} <= bodies:
            out.append(
                {
                    "id": r.get("id"),
                    "our_status": r.get("our_status"),
                    "bodies": sorted(bodies),
                    "vinf": [
                        (e.get("body"), e.get("vinf_kms"))
                        for e in (r.get("vinf_kms_at_encounters") or [])
                    ],
                }
            )
    return out


def empty_regions(a: str, b: str) -> list[str]:
    out = []
    for line in (REPO / "data" / "empty_regions.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        text = json.dumps(rec)
        names = {"V": "Venus", "M": "Mars", "E": "Earth"}
        if names.get(a, a) in text and names.get(b, b) in text:
            out.append(rec.get("region_id", "?"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    system, a, b = ENUM.cell_system(args.cell)
    n_struct, recs, groups, passing = ANALYSE.collect(args.dirs)
    n_struct_err = sum(
        1
        for d in args.dirs
        if (d / "structures.jsonl").exists()
        for line in (d / "structures.jsonl").open()
        if json.loads(line).get("error")
    )
    n_zero_err = sum(1 for r in recs if r.get("status") == "error")
    print(
        f"cell {args.cell}: structures {n_struct}, zeros {len(recs)}, physical {len(groups)}, "
        f"gate-passing {len(passing)}"
    )
    cat = catalogue_rows(a, b)
    reg = empty_regions(a, b)
    results = []
    for (k, seq), g in sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        best = min(
            (m for m in g["members"] if ANALYSE.is_pass(m)),
            key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9,
        )
        kk, cycle = ENUM.parse_cycle_key(best["key"], system, a, b)
        x = np.asarray(best["x_days"]) * DAY
        fl = cycle_flybys(system, cycle, x)
        assert fl is not None
        xc = cross_check(system, cycle, x, fl)
        soi = min(sphere_of_influence_km(system, c) for c in (a, b))
        line = {
            "k": k,
            "flybys": seq,
            "keys": sorted({m["key"] for m in g["members"]}),
            "key": best["key"],
            "mirror_pair": g["mirror_pair"],
            "vinf_kms": best["vinf_kms"],
        }
        sig = CandidateSignature(
            primary="Jupiter" if a in ("Ganymede", "Callisto", "Europa") else "Sun",
            sequence=tuple(f.body for f in fl),
            period_k=k,
            vinf_per_encounter_kms=tuple(round(f.vinf_kms, 3) for f in fl),
            topology_label=frozenset({"repeated-moon"}),
        )
        lit = check_literature(sig, search=offline_corpus_search)
        rec = line | {
            "x_days": best["x_days"],
            "period_days": kk * system.synodic_s(a, b) / DAY,
            "flyby_table": best["flybys"],
            "r_min_km": best["r_min_km"],
            "r_max_km": best["r_max_km"],
            "min_required_alt_km": best["min_required_alt_km"],
            "max_turn_deg": best["max_turn_deg"],
            "cross_check": xc,
            "soi_fraction": xc["max_arrival_miss_km"] / soi,
            "collisions": collisions(args.cell, line, system, a, b),
            "literature_offline": {
                "status": lit.status,
                "citation": getattr(lit, "citation", None),
                "notes": getattr(lit, "notes", None),
            },
        }
        results.append(rec)
        print(
            f"k={k} vinf={ {kk2: round(vv, 3) for kk2, vv in best['vinf_kms'].items()} } "
            f"xc_miss={xc['max_arrival_miss_km']:.2e}km "
            f"xc_dv={xc['max_vinf_vector_error_kms']:.1e} "
            f"gate_int={xc['gate_status_integrated']} lit={lit.status} "
            f"collisions={rec['collisions']}"
        )
    out = {
        "cell": args.cell,
        "n_structures": n_struct,
        "n_zeros": len(recs),
        "n_physical": len(groups),
        "n_gate_passing": len(passing),
        "n_structure_errors": n_struct_err,
        "n_zero_assessment_errors": n_zero_err,
        "catalogue_rows_same_pair": cat,
        "empty_regions_same_pair": reg,
        "candidates": results,
    }
    args.out.write_text(json.dumps(out, indent=1, default=float))
    print(f"catalogue rows on this pair: {len(cat)}; registry entries: {len(reg)}")
    print(
        f"ERRORS: structures {n_struct_err} (not searched, NOT empty), "
        f"zero assessments {n_zero_err}"
    )


if __name__ == "__main__":
    main()
