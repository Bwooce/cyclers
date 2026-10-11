"""#1024: the barycentric circular-model option for binaries, applied to the #998 Pluto cells.

Pre-registration: docs/notes/2026-10-11-1024-barycentric-binary-model.md (secs. 2-7). The option
is src/cyclerfinder/search/two_working_body_binary.py; the #998 cells, gauntlet rules and plu060
check are reused (scripts/run_998_enumerate.py, scripts/gauntlet_998.py logic,
scripts/realeph_998.py).

Subcommands:
  controls --out F                the sec. 5 controls 1-2 (off path bitwise; VenMar#45 and gc-1 with
                                  the option on at their low mass fractions)
  continue --cell pn [--members i/n] --out F
                                  sec. 4: continue each #998 gate-passing cyclers' representative
                                  zero from mass fraction 0 to mu_Charon, then re-assess (resumable:
                                  members already in F are skipped)
  realeph --out F                 sec. 5 control 3: the #998 sec. 6 defect for the 12 best strong
                                  cyclers, seeded from their continued barycentric zeros

Usage: uv run python scripts/run_1024_barycentric.py continue --cell pn --members 0/4 --out F
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    correct_dates,
    date_residual,
    eval_lambert_legs,
    kepler_step,
    sphere_of_influence_km,
)
from cyclerfinder.search.two_working_body_binary import (
    BarycentricCircularSystem,
    barycentric_binary,
    primary_min_distance_km,
    with_mass_fraction,
)
from cyclerfinder.search.two_working_body_enum import Zero, assess, flyby_table, leg_extent
from cyclerfinder.search.two_working_body_refly import cross_check

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
OUT_DIR = REPO / "data" / "1024_barycentric"
PLUTO_RADIUS_KM = 1188.3  # as scripts/gauntlet_998.py (Nimmo et al. 2017)
CHARON_RADIUS_KM = 606.0
CHARON_SOI_998_KM = 8070.0  # #998 note sec. 3 (registry sma 19,600 km)
GM_CHARON = 106.1  # registry


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
DRV998 = _load("run_998_enumerate", REPO / "scripts" / "run_998_enumerate.py")


def log(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')} {msg}", flush=True)


def _preflight(region: str, n: int) -> None:
    preflight_search(
        task_no=1024,
        region_id=region,
        method=MethodCapability(
            genome="#998 one-working-node cycles, barycentric binary option (continuation)",
            corrector="two_working_body.correct_dates, mass-fraction ramp",
            capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=n,
        timing_pilot_seconds_per_point=5.0,
    )


# ---------------------------------------------------------------------------------------------
# continuation
# ---------------------------------------------------------------------------------------------


def _lambert_vinf(system: CircularSystem, cyc: Cycle, x: np.ndarray) -> np.ndarray | None:
    legs = eval_lambert_legs(system, cyc, x)
    if legs is None:
        return None
    return np.array(
        [[float(np.linalg.norm(ev.vinf_dep)), float(np.linalg.norm(ev.vinf_arr))] for ev in legs]
    ).ravel()


def ramp(
    sys0: BarycentricCircularSystem,
    secondary: str,
    a_rel: float,
    mass_fraction: float,
    cyc: Cycle,
    x0: np.ndarray,
) -> dict[str, Any]:
    """Continue the zero x0 from mass fraction 0 to ``mass_fraction`` (sec. 4)."""
    lam, dlam, x = 0.0, 0.1, np.asarray(x0, dtype=float)
    v_prev = _lambert_vinf(sys0, cyc, x)
    steps = 0
    rejected = {"unconverged": 0, "jump": 0}
    while lam < 1.0:
        nxt = min(1.0, lam + dlam)
        if 1.0 - nxt < 1e-9:  # amendment A2: ten steps of 0.1 sum to 0.9999999999999999
            nxt = 1.0
        sysl = with_mass_fraction(sys0, secondary, a_rel, nxt * mass_fraction)
        sol = correct_dates(sysl, cyc, x, tol_kms=1e-8, max_nfev=80)
        ok = bool(sol.converged)
        reason = "unconverged"
        if ok:
            # the guard is per 0.1 of lambda, scaled to the step; the smallest legal step is 1/640
            scale = max((nxt - lam) / 0.1, (1.0 / 640.0) / 0.1)
            v_new = _lambert_vinf(sysl, cyc, sol.x)
            date_jump = float(np.max(np.abs(sol.x - x))) / DAY
            v_jump = (
                math.inf
                if v_new is None or v_prev is None
                else float(np.max(np.abs(v_new - v_prev)))
            )
            if date_jump > 0.5 * scale or v_jump > 0.02 * scale:
                ok, reason = False, "jump"
        steps += 1
        if ok:
            lam, x, v_prev = nxt, np.asarray(sol.x), v_new
            dlam = min(0.1, 2.0 * dlam)
        else:
            rejected[reason] += 1
            dlam /= 2.0
            if dlam < 1.0 / 640.0:
                return {
                    "outcome": "FOLDS/STOPS",
                    "last_lambda": lam,
                    "x": x,
                    "steps": steps,
                    "rejected": rejected,
                }
    return {
        "outcome": "CONTINUES",
        "last_lambda": 1.0,
        "x": x,
        "steps": steps,
        "rejected": rejected,
    }


def inside_fraction(
    system: BarycentricCircularSystem, cyc: Cycle, x: np.ndarray, radius_km: float, n: int = 400
) -> float:
    """Fraction of the Lambert-leg time inside ``radius_km`` of the barycentre (diagnostic)."""
    legs = eval_lambert_legs(system, cyc, x)
    assert legs is not None
    t_in = t_tot = 0.0
    for ev, li in zip(legs, cyc.lambert_index, strict=True):
        leg = cyc.legs[li]
        r0, w0 = system.state(leg.frm, ev.t_dep)  # type: ignore[union-attr]
        dt = ev.t_arr - ev.t_dep
        ts = np.linspace(0.0, dt, n)
        inside = sum(
            float(np.linalg.norm(kepler_step(r0, w0 + ev.vinf_dep, float(t), system.mu)[0]))
            < radius_km
            for t in ts
        )
        t_in += dt * inside / n
        t_tot += dt
    return t_in / t_tot if t_tot > 0 else 0.0


def reassess(system: BarycentricCircularSystem, cyc: Cycle, x: np.ndarray) -> dict[str, Any]:
    """The #998 assessment and strong filters at the barycentric zero (sec. 4)."""
    z = Zero(cyc, x, float(np.max(np.abs(date_residual(system, cyc, x)))))  # type: ignore[arg-type]
    ass = assess(system, z)
    out: dict[str, Any] = {
        "residual_kms": z.residual_kms,
        "status": ass.status,
        "status_hm_floor": ass.report_hm_floor.status if ass.report_hm_floor else None,
        "worst_ratio": ass.report.gate.worst_ratio if ass.report else None,
        "max_turn_deg": ass.report.max_turn_deg if ass.report else None,
        "max_encounter_miss_km": ass.max_encounter_miss_km,
        "vinf_kms": ass.vinf_kms,
    }
    if ass.flybys is None or ass.report is None:
        return out
    fl = ass.flybys
    table = flyby_table(system, ass)
    charon = [f for f in table if f["body"] == "Charon" and "ratio" in f]
    r_min, r_max = leg_extent(system, cyc, x, flybys=fl)
    pa = primary_min_distance_km(system, cyc, x, fl)
    soi_run = sphere_of_influence_km(system, "Charon")
    rp_max = max(float(f["required_alt_km"]) for f in charon) + CHARON_RADIUS_KM
    a_c = system.bodies["Charon"][0]
    xc = cross_check(system, cyc, x, fl)
    out |= {
        "flybys": table,
        "r_min_km": r_min,
        "r_max_km": r_max,
        "pluto_min_distance_km": pa.min_km,
        "pluto_surface_clear": bool(pa.min_km > PLUTO_RADIUS_KM),
        "max_required_rp_charon_km": rp_max,
        "charon_soi_run_km": soi_run,
        "rp_within_charon_soi": bool(rp_max <= soi_run),
        "rp_within_charon_soi_8070": bool(rp_max <= CHARON_SOI_998_KM),
        "r_min_over_a_charon": r_min / a_c,
        "lambert_time_fraction_inside_2a": inside_fraction(system, cyc, x, 2.0 * a_c),
        "cross_check": xc,
    }
    out["pass"] = bool(out["status"] == "pass" and (out["max_encounter_miss_km"] or 9.0) < 1.0)
    out["strong"] = bool(out["pass"] and out["pluto_surface_clear"] and out["rp_within_charon_soi"])
    out["strong_8070"] = bool(
        out["pass"] and out["pluto_surface_clear"] and out["rp_within_charon_soi_8070"]
    )
    return out


def cmd_continue(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="continue")
    ap.add_argument("--cell", required=True, choices=sorted(DRV998.CELLS))
    ap.add_argument("--members", default="0/1")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    moon = DRV998.CELLS[args.cell]
    g = json.loads(
        (REPO / "data" / "998_pluto_smallmoons" / f"{args.cell}_gauntlet.json").read_text()
    )
    i_sh, n_sh = (int(v) for v in args.members.split("/"))
    # gauntlet keys are not unique (two merged groups can share a representative key): index ids
    cands = [(j, c) for j, c in enumerate(g["candidates"]) if j % n_sh == i_sh]
    _preflight(f"pluto-{args.cell}-barycentric-continuation-1024", len(cands))
    done: dict[str, Any] = {}
    if args.out.exists():
        done = {r["index"]: r for r in json.loads(args.out.read_text())}
    off = DRV998.pluto_system(moon)
    a_rel = off.bodies["Charon"][0]
    mu_c = GM_CHARON / off.mu
    sys0 = barycentric_binary(off, "Charon", 0.0)
    sys1 = barycentric_binary(off, "Charon", mu_c)
    t_run = time.time()
    for jj, (j, c) in enumerate(cands):
        if j in done:
            continue
        t0 = time.time()
        rec: dict[str, Any] = {
            "index": j,
            "key": c["key"],
            "strong_998": c["strong"],
            "worst_ratio_998": c["worst_ratio"],
            "vinf_charon_998": c["vinf_charon_kms"],
            "vinf_moon_998": c["vinf_moon_kms"],
            "x_days_998": c["x_days"],
        }
        _, cyc = ENUM.parse_cycle_key(c["key"], off, "Charon", moon)
        x0 = np.asarray(c["x_days"]) * DAY
        r_off = date_residual(off, cyc, x0)
        r_zero = date_residual(sys0, cyc, x0)
        rec["lambda0_residual_off_kms"] = None if r_off is None else float(np.max(np.abs(r_off)))
        rec["lambda0_residual_kms"] = None if r_zero is None else float(np.max(np.abs(r_zero)))
        if "|H" not in c["key"]:
            a0 = assess(sys0, Zero(cyc, x0, 0.0))
            rec["lambda0_max_dvinf_charon_kms"] = abs(
                a0.vinf_kms.get("Charon", math.inf) - float(c["vinf_charon_kms"])
            )
        if "|H" in c["key"]:
            rec["outcome"] = "EXCLUDED (half-rev leg)"
        else:
            rr = ramp(sys0, "Charon", a_rel, mu_c, cyc, x0)
            rec |= {k: v for k, v in rr.items() if k != "x"}
            rec["x_days"] = (np.asarray(rr["x"]) / DAY).tolist()
            if rr["outcome"] == "CONTINUES":
                rec |= reassess(sys1, cyc, np.asarray(rr["x"]))
        rec["secs"] = time.time() - t0
        done[j] = rec
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(list(done.values()), indent=1, default=float))
        log(
            f"{args.cell} [{jj + 1}/{len(cands)}] #{j} {rec.get('outcome')} pass={rec.get('pass')} "
            f"strong={rec.get('strong')} worst={rec.get('worst_ratio')} "
            f"({rec['secs']:.1f}s; {time.time() - t_run:.0f}s total)"
        )


# ---------------------------------------------------------------------------------------------
# controls (sec. 5, controls 1-2)
# ---------------------------------------------------------------------------------------------

#: Low-mass-fraction on-path controls: (name, cell, key, gauntlet, V_inf pick, secondary,
#: secondary GM, primary GM of the cell's ``mu``, tolerance km/s = 10 mu v_S, published V_inf).
ON_PATH = [
    (
        "VenMar#45",
        "vm",
        "k2|LV>M/0s|LM>V/0s",
        "data/942_cell_vm_gauntlet.json",
        "V",
        324_860.0,
        8.6e-4,
        {"V": 8.22, "M": 12.96},
    ),
    (
        "gc-1",
        "gc",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        "data/943_cell_gc_gauntlet.json",
        "Ganymede",
        9_887.83,
        8.5e-3,
        None,
    ),
]
#: Off-path bitwise controls: (name, cell, key, gauntlet).
OFF_PATH = [
    ON_PATH[0][:4],
    ON_PATH[1][:4],
    (
        "GanCal#5",
        "gc",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s",
        "data/943_cell_gc_gauntlet.json",
    ),
]


def _stored(path: str, key: str) -> dict[str, Any]:
    g = json.loads((REPO / path).read_text())
    hits = [c for c in g["candidates"] if c["key"] == key]
    assert len(hits) == 1, (path, key, len(hits))
    return hits[0]


def cmd_controls(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="controls")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    out: dict[str, Any] = {"off_path": [], "on_path": []}
    for name, cell, key, path in OFF_PATH:
        system, a, b = ENUM.cell_system(cell)
        _, cyc = ENUM.parse_cycle_key(key, system, a, b)
        c = _stored(path, key)
        x = np.asarray(c["x_days"]) * DAY
        ass = assess(system, Zero(cyc, x, 0.0))
        dv = max(abs(ass.vinf_kms[q] - float(v)) for q, v in c["vinf_kms"].items())
        rec = {"name": name, "max_abs_dvinf_kms": dv, "status": ass.status, "ok": dv <= 1e-12}
        out["off_path"].append(rec)
        log(f"off path {name}: max |dV_inf| {dv:.1e} km/s, status {ass.status}, ok={rec['ok']}")
    for name, cell, key, path, sec, gm_s, tol, published in ON_PATH:
        system, a, b = ENUM.cell_system(cell)
        _, cyc = ENUM.parse_cycle_key(key, system, a, b)
        c = _stored(path, key)
        x0 = np.asarray(c["x_days"]) * DAY
        mu_s = gm_s / (system.mu + gm_s)
        a_rel = system.bodies[sec][0]
        sys0 = barycentric_binary(system, sec, 0.0)
        sys1 = barycentric_binary(system, sec, mu_s)
        rr = ramp(sys0, sec, a_rel, mu_s, cyc, x0)
        rec = {"name": name, "secondary": sec, "mass_fraction": mu_s, "tolerance_kms": tol,
               "outcome": rr["outcome"], "last_lambda": rr["last_lambda"]}  # fmt: skip
        if rr["outcome"] == "CONTINUES":
            ass = assess(sys1, Zero(cyc, np.asarray(rr["x"]), 0.0))
            off = assess(system, Zero(cyc, x0, 0.0))
            dv = max(abs(ass.vinf_kms[q] - off.vinf_kms[q]) for q in off.vinf_kms)
            rec |= {
                "vinf_on_kms": ass.vinf_kms,
                "vinf_off_kms": off.vinf_kms,
                "max_abs_dvinf_kms": dv,
                "status_on": ass.status,
                "status_off": off.status,
                "date_shift_s_max": float(np.max(np.abs(np.asarray(rr["x"]) - x0))),
            }
            ok = dv <= tol and ass.status == off.status
            if published:
                rec["published_dev_kms"] = max(
                    abs(ass.vinf_kms[q] - v) for q, v in published.items()
                )
                ok = ok and rec["published_dev_kms"] <= 0.05
            rec["ok"] = ok
        else:
            rec["ok"] = False
        out["on_path"].append(rec)
        log(f"on path {name}: {json.dumps(rec, default=float)}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1, default=float))
    bad = [r["name"] for r in out["off_path"] + out["on_path"] if not r["ok"]]
    log(f"CONTROLS NOT OK: {bad or 'none'}")


# ---------------------------------------------------------------------------------------------
# real-world test (sec. 5, control 3)
# ---------------------------------------------------------------------------------------------


def cmd_realeph(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="realeph")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    import spiceypy as sp

    from cyclerfinder.verify.spice_kernels import ensure_pluto_kernel

    reph = _load("realeph_998", REPO / "scripts" / "realeph_998.py")
    sp.furnsh(str(reph.LSK))
    sp.furnsh(ensure_pluto_kernel())
    et_a = sp.str2et("2030-01-01T00:00:00")
    frame = reph.mean_plane_frame(et_a, et_a + 4 * 365.25 * DAY)
    out = []
    for cell, moon in DRV998.CELLS.items():
        old = json.loads(
            (REPO / "data" / "998_pluto_smallmoons" / f"{cell}_realeph_n1.json").read_text()
        )
        recs = json.loads((OUT_DIR / f"{cell}_continue.json").read_text())
        off = DRV998.pluto_system(moon)
        for cand in old["candidates"]:
            # keys are not unique: match the #998 member by key and its ideal worst ratio
            hits = [
                x
                for x in recs
                if x["key"] == cand["key"]
                and abs(x["worst_ratio_998"] - cand["ideal"]["worst_ratio"]) < 1e-9
            ]
            assert len(hits) == 1, (cell, cand["key"], len(hits))
            r = hits[0]
            rec: dict[str, Any] = {"cell": cell, "key": cand["key"]}
            old_def = [e["residual_before_kms"] for e in cand["epochs"]]
            rec["defect_998_mps"] = [1e3 * v for v in old_def]
            if r is None or r.get("outcome") != "CONTINUES":
                rec["outcome"] = None if r is None else r.get("outcome")
                out.append(rec)
                continue
            _, cyc = ENUM.parse_cycle_key(cand["key"], off, "Charon", moon)
            x = np.asarray(r["x_days"]) * DAY
            new_def = []
            for ie in range(len(cand["epochs"])):
                et0 = reph.next_conjunction(frame, moon, et_a + ie * 400.0 * DAY)
                eph = reph.EphSystem(off, et0, frame)
                rb = date_residual(eph, cyc, x)
                new_def.append(None if rb is None else 1e3 * float(np.max(np.abs(rb))))
            rec["defect_barycentric_mps"] = new_def
            rec["max_998_mps"] = max(rec["defect_998_mps"])
            rec["max_barycentric_mps"] = max(v for v in new_def if v is not None)
            out.append(rec)
            log(
                f"{cell} {cand['key']}: defect max {rec['max_998_mps']:.1f} -> "
                f"{rec['max_barycentric_mps']:.1f} m/s"
            )
    args.out.write_text(json.dumps(out, indent=1, default=float))
    pairs = [
        (r["max_998_mps"], r["max_barycentric_mps"]) for r in out if "max_barycentric_mps" in r
    ]
    fell = sum(b < a for a, b in pairs)
    med_old = float(np.median([a for a, _ in pairs])) if pairs else math.nan
    med_new = float(np.median([b for _, b in pairs])) if pairs else math.nan
    log(f"defect fell for {fell} of {len(pairs)}; median {med_old:.1f} -> {med_new:.1f} m/s")


def main() -> None:
    cmds = {"continue": cmd_continue, "controls": cmd_controls, "realeph": cmd_realeph}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        raise SystemExit(f"usage: run_1024_barycentric.py {{{','.join(cmds)}}} ...")
    cmds[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
