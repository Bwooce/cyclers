"""#973: the two-working-body generator at longer cycle periods, and a Europa-Callisto cell.

Routes R11-R13 of the #971 corpus review (docs/notes/2026-10-07-971-fable-corpus-review-2.md
sec. 3), run with the EXISTING #942/#943 generator and driver (scripts/run_942_enumerate.py):

  R11  gc  Ganymede-Callisto, k = 4, 5, 6 (cell-3 settings of the #943 run)
  R13  ev  Earth-Venus, k = 4, 5 (cell-5 settings of the #942 run)
  R12  ec  Europa-Callisto, both massive, k = 1-4 (NEW cell, defined here)

Subcommands:

  enumerate ARGS    the run_942 driver unchanged (same flags, same outputs), with this task's
                    preflight and the ec cell added to its cell table
  recall --out F    the in-run controls: published / earlier structures re-solved with the
                    production seeds (gc-1, GanCal#5 at gc k = 3; ev-A, ev-C, Hollister 1H at
                    ev k = 2; the #576 Europa-Callisto symmetric closures at ec k = 3, 4, 5)
  liang --out F     Liang et al. 2024 Table 3 C-G and C-E halves as OPEN two-leg segments in
                    the gc / ec cells (not closures)
  gauntlet --cell C DIR [DIR ...] --out F
                    scripts/gauntlet_942.py's DOP853 re-fly, SOI check and literal-collision
                    checks, WITHOUT its literature step (deferred until #972 lands), plus the
                    Liang 2024 Tables 3/5/7 and #576 closure comparisons

The ec cell: Europa and Callisto with Russell & Strange 2009 Table 2 constants (p.148; the
same ``rs_moon_system`` as the gc and ge cells) and the registry floors (Europa 100 km,
Callisto 200 km). Europa is body A (inner), as Ganymede is in gc.

Usage:
  uv run python scripts/run_973_enumerate.py enumerate --cell gc --k 4 --out DIR \
      --max-returns 1,1 --transfer-revs 0,1,2 --generic-revs 1,2 \
      --n-phase 36 --n-split 12 --n-refine 40 [--timing-pilot-s S]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
import types
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import CircularSystem, cycle_flybys, moon_circular
from cyclerfinder.search.two_working_body_enum import (
    assess,
    flyby_table,
    leg_extent,
    solve_structure,
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

#: Cells this task adds to run_942's table (additive; run_942's cells are unchanged).
EXTRA_CELLS = frozenset({"ec"})
_RUN_942_CELL_SYSTEM = ENUM.cell_system


def cell_system(cell: str) -> tuple[CircularSystem, str, str]:
    if cell == "ec":
        return ENUM.rs_moon_system(["Europa", "Callisto"], []), "Europa", "Callisto"
    if cell == "ec576":
        # the #563/#576 symmetric-closure model (registry sma, registry Jupiter GM), controls only
        return (
            moon_circular(PRIMARIES["Jupiter"], ["Europa", "Callisto"]),
            "Europa",
            "Callisto",
        )
    result: tuple[CircularSystem, str, str] = _RUN_942_CELL_SYSTEM(cell)
    return result


def _preflight_973(**kwargs: Any) -> None:
    preflight_search(task_no=973, script_path=Path(__file__), **kwargs)


def cmd_enumerate(argv: list[str]) -> None:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--cell", required=True)
    cell = pre.parse_known_args(argv)[0].cell
    ENUM.cell_system = cell_system  # main() looks the cell up through the module global
    sys.argv = [sys.argv[0], *argv]
    try:
        ENUM.main(preflight=_preflight_973, x1=cell in ENUM.X1_CELLS)
    finally:
        ENUM.cell_system = _RUN_942_CELL_SYSTEM


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------

#: Production seeds of the #942/#943 cells (generator note sec. 6.1).
SEEDS = {"n_phase": 36, "n_split": 12, "n_refine": 40}

#: Recall controls: (name, cell, key, source of the expected values, tolerance km/s).
#: The #942/#943 rows are code-path regressions: the expected V_inf are the full-precision values
#: those runs stored (gauntlet JSON), at 1e-3 km/s. GanCal#5 is also held to its published
#: R-S 2009 Table 3 values (3.24 / 3.34 km/s) at the #942 LITERAL tolerance, 0.05 km/s.
RECALL = [
    (
        "gc-1 (#943 cell gc, k = 3)",
        "gc",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        "data/943_cell_gc_gauntlet.json",
        1e-3,
    ),
    (
        "GanCal#5 as found by the #943 gc run (k = 3)",
        "gc",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s",
        "data/943_cell_gc_gauntlet.json",
        1e-3,
    ),
    (
        "GanCal#5 against R-S 2009 Table 3 (published 3.24 / 3.34 km/s)",
        "gc",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s",
        "published",
        0.05,
    ),
    (
        "ev-A (#942 cell ev, k = 2)",
        "ev",
        "k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s",
        "data/942_cell_ev_gauntlet.json",
        1e-3,
    ),
    (
        "ev-C (#942 cell ev, k = 2)",
        "ev",
        "k2|LE>V/0s|LV>V/1h|LV>E/0s",
        "data/942_cell_ev_gauntlet.json",
        1e-3,
    ),
    (
        "Hollister 1H, orbit I topology (#942 in-run control, k = 2, the 2.99 / 3.19 member)",
        "ev",
        "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s",
        "data/942_cell_ev_gauntlet.json",
        1e-3,
    ),
]
#: Picked by hand where one key holds more than one stored gate-passer (ev-C's key also holds
#: the 4.40 / 8.82 cycler; 1H's key also holds 6.14 / 5.28). Full-precision stored values.
RECALL_PICK = {
    "ev-C (#942 cell ev, k = 2)": {"E": 9.074921385818076, "V": 13.166434198896294},
    "Hollister 1H, orbit I topology (#942 in-run control, k = 2, the 2.99 / 3.19 member)": {
        "E": 2.994071820948365,
        "V": 3.190271920792486,
    },
}
PUBLISHED = {
    "GanCal#5 against R-S 2009 Table 3 (published 3.24 / 3.34 km/s)": {
        "Ganymede": 3.24,
        "Callisto": 3.34,
    }
}

#: #576 control tolerance in its own model (registry moons, cell ``ec576``): 1e-3 km/s. In the ec
#: cell (R-S Table 2 constants) the closest zero is reported, not judged (#973 note sec. 2).
TOL_576 = 1e-3


def _expected(name: str, key: str, src: str) -> dict[str, float]:
    if src == "published":
        return PUBLISHED[name]
    if name in RECALL_PICK:
        return RECALL_PICK[name]
    d = json.loads((REPO / src).read_text())
    hits = [c for c in d["candidates"] if key in c["keys"]]
    if len(hits) != 1:
        raise SystemExit(f"{name}: {len(hits)} stored candidates hold {key}")
    return {c: float(v) for c, v in hits[0]["vinf_kms"].items()}


def _576_ec_closures() -> list[dict[str, Any]]:
    """The #576 Europa-Callisto symmetric closures (anchor Europa; the anchor-Callisto rows are
    the same cyclers started at the other moon)."""
    out = []
    path = REPO / "data" / "enumerate_576_jupiter_galilean_symmetric_closures.jsonl"
    for line in path.read_text().splitlines():
        d = json.loads(line)
        if d.get("kind") == "pass" and {d["anchor"], d["flyby"]} == {"Europa", "Callisto"}:
            out.append(d)
    return out


def _solve_and_assess(cell: str, key: str) -> list[dict[str, Any]]:
    system, a, b = cell_system(cell)
    k, cyc = ENUM.parse_cycle_key(key, system, a, b)
    zeros = solve_structure(system, cyc, phase_period_s=system.synodic_s(a, b), **SEEDS)
    out = []
    for z in zeros:
        ass = assess(system, z)
        out.append(
            {
                "key": key,
                "k": k,
                "x_days": (z.x / DAY).tolist(),
                "residual_kms": z.residual_kms,
                "status": ass.status,
                "worst_ratio": ass.report.gate.worst_ratio if ass.report else None,
                "max_encounter_miss_km": ass.max_encounter_miss_km,
                "vinf_kms": ass.vinf_kms,
                "flybys": flyby_table(system, ass),
            }
        )
    return out


def _match(zeros: list[dict[str, Any]], want: dict[str, float], tol: float) -> list[dict]:
    return [
        z
        for z in zeros
        if z["vinf_kms"]
        and all(abs(z["vinf_kms"].get(c, math.inf) - v) <= tol for c, v in want.items())
    ]


def cmd_recall(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="recall")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", type=str, default="", help="comma list of cells (gc, ev, ec)")
    args = ap.parse_args(argv)
    cells = {c for c in args.only.split(",") if c}
    controls = [
        {"name": n, "cell": c, "key": k, "expect": (_expected(n, k, src), tol, src)}
        for n, c, k, src, tol in RECALL
    ]
    # #576 closures: E->C->E, both legs 0-rev, period n synodic (tof per leg n T_syn / 2). Each
    # physical cycler is listed twice in #576 (anchor Europa, anchor Callisto); both cyclic
    # orders are solved, in #576's own model (ec576) and in the ec cell.
    for d in _576_ec_closures():
        n = d["n_commensurate_int"]
        v = d["vinf_per_encounter_kms"]
        an, fl = d["anchor"], d["flyby"]
        key = f"k{n}|L{an}>{fl}/0s|L{fl}>{an}/0s"
        for cell in ("ec576", "ec"):
            controls.append(
                {
                    "name": f"#576 closure anchor {an} n = {n} in {cell}",
                    "cell": cell,
                    "key": key,
                    "expect": ({an: v[0], fl: v[1]}, TOL_576, "#576 symmetric-closure jsonl"),
                }
            )
    results = []
    t0 = time.time()
    for c in controls:
        if cells and c["cell"] not in cells:
            continue
        zeros = _solve_and_assess(c["cell"], c["key"])
        want, tol, src = c["expect"]
        hits = _match(zeros, want, tol)
        rec = {
            "name": c["name"],
            "cell": c["cell"],
            "key": c["key"],
            "expected_vinf_kms": want,
            "tolerance_kms": tol,
            "expected_source": src,
            "n_zeros": len(zeros),
            "recalled": bool(hits),
            "hits": hits,
            "all_zeros": zeros,
        }
        results.append(rec)
        best = min(
            (
                max(abs(z["vinf_kms"].get(b, math.inf) - v) for b, v in want.items())
                for z in zeros
                if z["vinf_kms"]
            ),
            default=math.inf,
        )
        print(
            f"{time.time() - t0:7.1f}s {c['name']}: {len(zeros)} zeros, recalled={bool(hits)} "
            f"(closest max dV_inf {best:.4f} km/s; tol {tol}); hit statuses "
            f"{[h['status'] for h in hits]}",
            flush=True,
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=1, default=float))


# ---------------------------------------------------------------------------
# Liang 2024 open segments
# ---------------------------------------------------------------------------


def _liang_segment(cell: str, member: str, legs: tuple[int, int]) -> dict[str, Any]:
    """One two-leg half of Liang's CGCEC cycle (Table 3/5/7 legs ``legs``) in our cell.

    Liang's moons are placed by their Table 1 mean motions and phases (cge_scaffold, the #222
    reproduction). Our system is rotated so that the first moon is at Liang's angle at the
    segment start and the second at Liang's angle at the middle flyby; the end moon's angle
    then differs from Liang's only by the period difference of the two models (reported).
    The two Lambert legs are solved at Liang's printed times of flight, on every revolution
    count 0-4 and branch; the reported solution is the one closest to the printed V_inf."""
    from cyclerfinder.search import cge_scaffold as cg

    spec = cg.LIANG_MEMBERS[member]
    rep = cg.reproduce_member(member)
    seq = cg.CGCEC_SEQUENCE
    i0, i1 = legs
    m0, m1, m2 = seq[i0], seq[i1], seq[i1 + 1]
    t0, t1, t2 = (rep.flybys[j].epoch_days for j in (i0, i1, i1 + 1))
    liang_angle = {
        j: spec.phases_rad[seq[j]]
        + cg.LIANG_MEAN_MOTIONS_RAD_DAY[seq[j]] * rep.flybys[j].epoch_days
        for j in (i0, i1, i1 + 1)
    }
    system, _a, _b = cell_system(cell)
    # theta(t) = theta0 + 2 pi t / P; t measured from the segment start
    bodies = dict(system.bodies)
    for moon, j, tj in ((m0, i0, 0.0), (m1, i1, (t1 - t0) * DAY)):
        sma, per, _ = bodies[moon]
        bodies[moon] = (sma, per, liang_angle[j] - 2.0 * math.pi * tj / per)
    sys2 = replace(system, bodies=bodies, _fb={})
    _sma2, per2, th2 = bodies[m2]
    ours_end = th2 + 2.0 * math.pi * (t2 - t0) * DAY / per2
    end_angle_diff = math.remainder(ours_end - liang_angle[i1 + 1], 2.0 * math.pi)
    printed = spec.vinf_printed_kms
    out_legs = []
    vin: dict[int, np.ndarray] = {}
    vout: dict[int, np.ndarray] = {}
    for fa, fb, ja, jb in ((m0, m1, i0, i1), (m1, m2, i1, i1 + 1)):
        ta = (rep.flybys[ja].epoch_days - t0) * DAY
        tb = (rep.flybys[jb].epoch_days - t0) * DAY
        ra, wa = sys2.state(fa, ta)
        rb, wb = sys2.state(fb, tb)
        sols = lambert(ra, rb, tb - ta, mu=sys2.mu, max_revs=4)

        def cost(s: Any, ja: int = ja, jb: int = jb, wa: Any = wa, wb: Any = wb) -> float:
            return abs(float(np.linalg.norm(s.v1 - wa)) - printed[ja]) + abs(
                float(np.linalg.norm(s.v2 - wb)) - printed[jb]
            )

        ranked = sorted(sols, key=cost)
        best = ranked[0]
        vout[ja] = best.v1 - wa
        vin[jb] = best.v2 - wb
        out_legs.append(
            {
                "from": fa,
                "to": fb,
                "tof_days": (tb - ta) / DAY,
                "n_revs": best.n_revs,
                "branch": best.branch,
                "vinf_dep_kms": float(np.linalg.norm(vout[ja])),
                "vinf_arr_kms": float(np.linalg.norm(vin[jb])),
                "printed_dep_arr_kms": (printed[ja], printed[jb]),
                "selection_margin_kms": cost(ranked[1]) - cost(best) if len(ranked) > 1 else None,
                "n_lambert_solutions": len(sols),
            }
        )
    mid_in = float(np.linalg.norm(vin[i1]))
    mid_out = float(np.linalg.norm(vout[i1]))
    # pass bar (#973 note sec. 2): the #222 print tolerance at Liang's epoch of each flyby
    tol = {
        j: cg.vinf_print_tolerance_kms(rep.flybys[j].epoch_days, seq[j], rep.radii_km)
        for j in (i0, i1, i1 + 1)
    }
    vals = [
        (m0, out_legs[0]["vinf_dep_kms"], printed[i0], tol[i0]),
        (m1, mid_in, printed[i1], tol[i1]),
        (m1, mid_out, printed[i1], tol[i1]),
        (m2, out_legs[1]["vinf_arr_kms"], printed[i1 + 1], tol[i1 + 1]),
    ]
    return {
        "member": member,
        "cell": cell,
        "sequence": [m0, m1, m2],
        "segment_days": t2 - t0,
        "legs": out_legs,
        "middle_flyby_vinf_in_out_kms": (mid_in, mid_out),
        "middle_flyby_mismatch_kms": abs(mid_in - mid_out),
        "max_abs_vinf_minus_printed_kms": max(abs(v - p) for _, v, p, _t in vals),
        "per_flyby": [
            {"moon": m, "ours_kms": v, "printed_kms": p, "tolerance_kms": t} for m, v, p, t in vals
        ],
        "reproduced": all(abs(v - p) <= t for _, v, p, t in vals)
        and abs(mid_in - mid_out) <= tol[i1],
        "end_moon_angle_offset_deg": math.degrees(end_angle_diff),
        "liang_member_reproduction_worst_kms": rep.max_vinf_residual_kms,
    }


def cmd_liang(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="liang")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    res = []
    for member in ("A", "B", "C"):
        for cell, legs in (("gc", (0, 1)), ("ec", (2, 3))):
            r = _liang_segment(cell, member, legs)
            res.append(r)
            print(
                f"member {member} {'-'.join(r['sequence'])} ({r['segment_days']:.4f} d): legs "
                f"{[(lg['n_revs'], lg['branch']) for lg in r['legs']]} "
                f"max |V_inf - printed| {r['max_abs_vinf_minus_printed_kms']:.4f} km/s, "
                f"middle |in|-|out| {r['middle_flyby_mismatch_kms']:.4f}, end-moon angle "
                f"offset {r['end_moon_angle_offset_deg']:.4f} deg; reproduced={r['reproduced']}",
                flush=True,
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(res, indent=1, default=float))


# ---------------------------------------------------------------------------
# Gauntlet (literature step deferred)
# ---------------------------------------------------------------------------


def _load_gauntlet_without_literature() -> Any:
    """scripts/gauntlet_942.py with its literature_check import replaced by a stub.

    The #942/#943 gauntlet's re-fly, collision and catalogue functions are reused as they are;
    the scripted literature step is deferred for #973 (literature_check.py is being edited for
    #972), so the stub raises if anything calls it."""

    def _deferred(*_a: Any, **_k: Any) -> Any:
        raise RuntimeError("#973: the literature step is deferred until #972 lands")

    stub = types.ModuleType("cyclerfinder.search.literature_check")
    stub.CandidateSignature = _deferred  # type: ignore[attr-defined]
    stub.check_literature = _deferred  # type: ignore[attr-defined]
    stub.offline_corpus_search = _deferred  # type: ignore[attr-defined]
    saved = sys.modules.get("cyclerfinder.search.literature_check")
    sys.modules["cyclerfinder.search.literature_check"] = stub
    try:
        return _load("gauntlet_942", REPO / "scripts" / "gauntlet_942.py")
    finally:
        if saved is None:
            del sys.modules["cyclerfinder.search.literature_check"]
        else:
            sys.modules["cyclerfinder.search.literature_check"] = saved


def _liang_rows() -> list[tuple[str, dict[str, list[float]]]]:
    """Liang et al. 2024 Tables 3, 5, 7 printed V_inf per moon (transcribed in cge_scaffold)."""
    from cyclerfinder.search import cge_scaffold as cg

    rows = []
    for m, spec in cg.LIANG_MEMBERS.items():
        per: dict[str, list[float]] = {}
        for moon, v in zip(cg.CGCEC_SEQUENCE, spec.vinf_printed_kms, strict=True):
            per.setdefault(moon, []).append(v)
        rows.append((f"Liang 2024 member {m}", per))
    return rows


def liang_576_collisions(line: dict[str, Any]) -> list[str]:
    out = []
    v = line["vinf_kms"]
    for name, per in _liang_rows():
        if not set(v) <= set(per):
            continue
        dv = max(min(abs(v[c] - x) for x in per[c]) for c in v)
        if dv <= 0.3:
            out.append(f"NEAR {name} (per-moon dV_inf {dv:.3f}; a three-moon CGCEC cycle)")
    for d in _576_ec_closures() + _576_gc_closures():
        pair = {
            d["anchor"]: d["vinf_per_encounter_kms"][0],
            d["flyby"]: d["vinf_per_encounter_kms"][1],
        }
        if set(pair) != set(v):
            continue
        dv = max(abs(v[c] - pair[c]) for c in v)
        if dv <= 0.3:
            out.append(
                f"NEAR #576 symmetric closure {d['anchor']}-{d['flyby']} n = "
                f"{d['n_commensurate_int']} (dV_inf {dv:.3f}, k {d['n_commensurate_int']} vs "
                f"{line['k']})"
            )
    return out


def _576_gc_closures() -> list[dict[str, Any]]:
    path = REPO / "data" / "enumerate_576_jupiter_galilean_symmetric_closures.jsonl"
    return [
        d
        for d in map(json.loads, path.read_text().splitlines())
        if d.get("kind") == "pass" and {d["anchor"], d["flyby"]} == {"Ganymede", "Callisto"}
    ]


def cmd_gauntlet(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="gauntlet")
    ap.add_argument("--cell", required=True)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    g = _load_gauntlet_without_literature()
    system, a, b = cell_system(args.cell)
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
    cat = g.catalogue_rows(a, b)
    reg = g.empty_regions(a, b)
    results = []
    for (k, seq), grp in sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        best = min(
            (m for m in grp["members"] if ANALYSE.is_pass(m)),
            key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9,
        )
        kk, cycle = ENUM.parse_cycle_key(best["key"], system, a, b)
        x = np.asarray(best["x_days"]) * DAY
        fl = cycle_flybys(system, cycle, x)
        assert fl is not None
        ext = leg_extent(system, cycle, x, flybys=fl)
        xc = g.cross_check(system, cycle, x, fl)
        soi = min(g.sphere_of_influence_km(system, c) for c in (a, b))
        line = {
            "k": k,
            "flybys": seq,
            "keys": sorted({m["key"] for m in grp["members"]}),
            "key": best["key"],
            "mirror_pair": grp["mirror_pair"],
            "vinf_kms": best["vinf_kms"],
        }
        rec = line | {
            "x_days": best["x_days"],
            "period_days": kk * system.synodic_s(a, b) / DAY,
            "worst_ratio": best["worst_ratio"],
            "status_hm_floor": best["status_hm_floor"],
            "flyby_table": best["flybys"],
            "r_min_km": ext[0],
            "r_max_km": ext[1],
            "min_required_alt_km": best["min_required_alt_km"],
            "max_turn_deg": best["max_turn_deg"],
            "cross_check": xc,
            "soi_fraction": xc["max_arrival_miss_km"] / soi,
            "collisions": g.collisions(args.cell, line, system, a, b)
            + g.catalogue_matches(line, k, cat)
            + liang_576_collisions(line),
            "literature": "DEFERRED (#973 note; runs after #972 lands)",
        }
        results.append(rec)
        print(
            f"k={k} vinf={ {c: round(v, 3) for c, v in best['vinf_kms'].items()} } "
            f"worst={best['worst_ratio']:.3f} xc_miss={xc['max_arrival_miss_km']:.2e}km "
            f"xc_dv={xc['max_vinf_vector_error_kms']:.1e} gate_int={xc['gate_status_integrated']} "
            f"collisions={rec['collisions']}",
            flush=True,
        )
    out = {
        "cell": args.cell,
        "dirs": [str(d) for d in args.dirs],
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


def main() -> None:
    cmds = {"enumerate": cmd_enumerate, "recall": cmd_recall, "liang": cmd_liang}
    cmds["gauntlet"] = cmd_gauntlet
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        raise SystemExit(f"usage: run_973_enumerate.py {{{','.join(cmds)}}} ...")
    cmds[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
