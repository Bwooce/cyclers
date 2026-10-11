"""#1054: the #1000 seed-design fix (perigee at least 1.3 times the floor, omega from the target
lunar-pass distance, a perigee step limit). Pre-registration:
docs/notes/2026-10-11-1054-seed-design-fix.md.

    uv run python scripts/run_1054_seed_design.py predict   # old records' rates per seed subset
    uv run python scripts/run_1054_seed_design.py arm B|C|D --max-seconds 420   # pilot ablation
    uv run python scripts/run_1054_seed_design.py control --max-seconds 420     # 7-3b control
    uv run python scripts/run_1054_seed_design.py report

`predict` reads the #1000 grid and #1048 gaps records and tabulates, per seed subset (all; r_p at
least 1.3 x floor; skeleton lunar pass in the target band; both), the impact, converged and
cycler-class-candidate rates and the families found. These are the predictions for the two
seed-selection parts of the fix; only the step limit adds dynamics the old records cannot show.
Output: data/1054_seed_design/predict.json.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_shooting as sh

Arr = NDArray[np.float64]
REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "data" / "1000_complement" / "shooting"
OUT = REPO / "data" / "1054_seed_design"
L_KM = sh.base.L_KM
FLOOR_KM = sh.base.FLOOR_E
RP_MIN_KM = 1.3 * FLOOR_KM
BAND_KM = (20000.0, 60000.0)
# pre-registration sec. 2
RP_NEW = [r / L_KM for r in np.geomspace(RP_MIN_KM, sh.base.GEO_KM, 8)]
TARGETS_KM = (5000.0, 10000.0, 20000.0, 40000.0, 80000.0)
OM_SCAN = np.radians(np.arange(360.0))
DROP = 0.2
N_ARM = 240


def skeleton_min_moon_km(p: int, q: int, sense: int, rp: float, omegas: Arr) -> Arr:
    """Closest approach to the Moon (km) of the two-body p:q skeleton over q lunar periods, for
    each omega (vectorised over omega; the #1048 per-seed function, same sampling)."""
    a = (q / p) ** (2.0 / 3.0) * sh.GM_E ** (1.0 / 3.0)
    e = 1.0 - rp / a
    n = math.sqrt(sh.GM_E / a**3)
    t = np.linspace(0.0, 2 * math.pi * q, 6000)
    m_anom = n * t
    ecc = m_anom.copy()
    for _ in range(30):
        ecc = ecc - (ecc - e * np.sin(ecc) - m_anom) / (1 - e * np.cos(ecc))
    nu = 2 * np.arctan2(math.sqrt(1 + e) * np.sin(ecc / 2), math.sqrt(1 - e) * np.cos(ecc / 2))
    r = a * (1 - e * np.cos(ecc))
    ang = np.asarray(omegas)[:, None] + sense * nu[None, :] - t[None, :]
    x = r[None, :] * np.cos(ang) - sh.MU
    y = r[None, :] * np.sin(ang)
    return np.asarray(np.min(np.hypot(x - 1 + sh.MU, y), axis=1) * L_KM)


def target_omegas(p: int, q: int, sense: int, rp: float) -> list[tuple[float, float]]:
    """F2: every omega where the skeleton's closest lunar approach crosses a target distance
    (1-degree scan, linear interpolation); returns (omega, target_km) pairs."""
    d = skeleton_min_moon_km(p, q, sense, rp, OM_SCAN)
    out = []
    for tgt in TARGETS_KM:
        g = d - tgt
        for i in range(len(OM_SCAN)):
            j = (i + 1) % len(OM_SCAN)
            if g[i] == 0.0 or g[i] * g[j] < 0.0:
                w = g[i] / (g[i] - g[j]) if g[i] != g[j] else 0.0
                out.append((float(OM_SCAN[i] + w * math.radians(1.0)) % (2 * math.pi), tgt))
    return out


def new_seeds() -> tuple[list[dict[str, Any]], int]:
    """F1 + F2 seeds for the #1000 PQ set (lunar-orbit crossing only); count of skipped combos."""
    seeds, skipped = [], 0
    for p, q in sh.PQ:
        a = (q / p) ** (2.0 / 3.0) * sh.GM_E ** (1.0 / 3.0)
        for sense in (1, -1):
            for irn, rp in enumerate(RP_NEW):
                if not (rp < 1.0 < 2 * a - rp):
                    continue
                oms = target_omegas(p, q, sense, rp)
                if not oms:
                    skipped += 1
                for k, (om, tgt) in enumerate(oms):
                    seeds.append(
                        {"p": p, "q": q, "sense": sense, "irn": irn, "k": k, "om": om, "tgt": tgt}
                    )
    return seeds, skipped


def old_seeds() -> list[dict[str, Any]]:
    seeds = []
    for p, q in sh.PQ:
        a = (q / p) ** (2.0 / 3.0) * sh.GM_E ** (1.0 / 3.0)
        for sense in (1, -1):
            for ir, rp in enumerate(sh.RP):
                if not (rp < 1.0 < 2 * a - rp):
                    continue
                for iw in range(len(sh.OMEGA)):
                    seeds.append({"p": p, "q": q, "sense": sense, "ir": ir, "iw": iw})
    return seeds


def seed_values(sd: dict[str, Any]) -> tuple[float, float]:
    if "ir" in sd:
        return sh.RP[sd["ir"]], sh.OMEGA[sd["iw"]]
    return RP_NEW[sd["irn"]], sd["om"]


def skeleton_arc_impact(sd: dict[str, Any]) -> bool:
    """True when any arc of the uncorrected skeleton hits the Earth or the Moon (the test's
    first-arc impact: the mechanism of the old 93-94 % bins)."""
    rp, om = seed_values(sd)
    nodes, tp = sh.skeleton_nodes(sd["p"], sd["q"], sd["sense"], rp, om)
    for s0 in nodes:
        sol = solve_ivp(
            sh.base.eom,
            (0.0, tp),
            s0,
            args=(sh.MU,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
            events=[sh.base._hit_e, sh.base._hit_m],
        )
        if sol.status == 1:
            return True
    return False


def arms() -> dict[str, list[dict[str, Any]]]:
    """Pilot samples (pre-registration sec. 4): 240 old seeds (seed 1054), and 240 new seeds with
    the same (p, q, sense) counts; arms C and D share the new sample."""
    rng = random.Random(1054)
    old = rng.sample(old_seeds(), N_ARM)
    new, _ = new_seeds()
    by: dict[tuple[int, int, int], list[dict[str, Any]]] = {}
    for sd in new:
        by.setdefault((sd["p"], sd["q"], sd["sense"]), []).append(sd)
    want = Counter((sd["p"], sd["q"], sd["sense"]) for sd in old)
    new_s = []
    for key, n in sorted(want.items()):
        new_s += rng.sample(by[key], min(n, len(by[key])))
    # E (F1 alone, pre-registration sec. 6): A's seeds paired, same (p, q, sense, row index,
    # omega), with the perigee row moved to the F1 value RP_NEW[ir]; no F3
    f1 = [
        {"p": x["p"], "q": x["q"], "sense": x["sense"], "irn": x["ir"], "k": x["iw"]}
        | {"om": sh.OMEGA[x["iw"]]}
        for x in old
    ]
    return {"B": old, "C": new_s, "D": new_s, "E": f1}


def solve(job: tuple[dict[str, Any], float | None, float | None]) -> dict[str, Any]:
    sd, drop, c_target = job
    rp, om = seed_values(sd)
    rec = dict(sd)
    rec["drop"] = drop
    rec["C_target"] = c_target
    t0 = time.time()
    rec = sh.solve_values(rec, sd["p"], sd["q"], sd["sense"], rp, om, c_target, drop)
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def seed_id(sd: dict[str, Any]) -> str:
    keys = ("p", "q", "sense", "ir", "iw", "irn", "k")
    return "|".join(str(sd.get(k)) for k in keys)


def run_jobs(
    tag: str, jobs: list[tuple[dict[str, Any], float | None, float | None]], max_s: float
) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{tag}.jsonl"
    done = set()
    if f.exists():
        done = {seed_id(json.loads(x)) for x in f.read_text().splitlines()}
    todo = [j for j in jobs if seed_id(j[0]) not in done]
    print(f"{time.strftime('%H:%M:%S')} {tag}: {len(done)} done, {len(todo)} to go", flush=True)
    t0, n = time.time(), 0
    with Pool(2) as pool, f.open("a") as fh:
        for rec in pool.imap_unordered(solve, todo):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n += 1
            if n % 25 == 0:
                print(
                    f"{time.strftime('%H:%M:%S')} {n} ({(time.time() - t0) / n:.2f} s/seed)",
                    flush=True,
                )
            if time.time() - t0 > max_s:
                print(f"time budget reached after {n}", flush=True)
                pool.terminate()
                return
    print(f"{time.strftime('%H:%M:%S')} {tag}: complete ({n} this call)", flush=True)


def stats(recs: list[dict[str, Any]]) -> dict[str, Any]:
    recs = [r for r in recs if r["status"] != "no lunar-orbit crossing"]
    n = len(recs)
    st: Counter[str] = Counter()
    for r in recs:
        s = r["status"] if r["status"] in ("converged", "impact") else "other"
        if s == "impact":
            s = "impact@1" if r["it"] == 1 else "impact@>1"
        st[s] += 1
    cand = [r for r in recs if r.get("cycler_class_candidate")]
    return {
        "n": n,
        "impact": round((st["impact@1"] + st["impact@>1"]) / n, 3),
        "impact@1": round(st["impact@1"] / n, 3),
        "impact@>1": round(st["impact@>1"] / n, 3),
        "converged": round(st["converged"] / n, 3),
        "other": round(st["other"] / n, 3),
        "candidates_per_solve": round(len(cand) / n, 4),
        "families": sorted({family_key(r) for r in cand}),
        "steps_halved_mean": round(float(np.mean([r.get("n_halved", 0) for r in recs])), 2),
    }


def family_key(r: dict[str, Any]) -> str:
    sense = "pro" if r["sense"] > 0 else "ret"
    sym = "sym" if r["symmetric"] else "asym"
    return f"{r['p']}:{r['q']} {sense} wE{round(r['wind_E'])} wM{round(r['wind_M'])} {sym}"


def predict() -> None:
    omegas = np.asarray(sh.OMEGA)
    out: dict[str, Any] = {"rp_min_km": RP_MIN_KM, "band_km": BAND_KM}
    for name in ("grid", "gaps"):
        recs = [json.loads(x) for x in (SRC / f"{name}.jsonl").read_text().splitlines()]
        recs = [r for r in recs if r["status"] != "no lunar-orbit crossing"]
        cache: dict[tuple[int, int, int, int], Arr] = {}
        for r in recs:
            k = (r["p"], r["q"], r["sense"], r["ir"])
            if k not in cache:
                cache[k] = skeleton_min_moon_km(r["p"], r["q"], r["sense"], sh.RP[r["ir"]], omegas)
            r["d_moon_km"] = float(cache[k][r["iw"]])
        subsets = {
            "all": lambda r: True,
            "rp>=1.3floor": lambda r: sh.RP[r["ir"]] * L_KM >= RP_MIN_KM - 1.0,
            "band": lambda r: BAND_KM[0] <= r["d_moon_km"] <= BAND_KM[1],
            "rp>=1.3floor & band": lambda r: (
                sh.RP[r["ir"]] * L_KM >= RP_MIN_KM - 1.0
                and BAND_KM[0] <= r["d_moon_km"] <= BAND_KM[1]
            ),
            "rp>=1.3floor & d>=5000": lambda r: (
                sh.RP[r["ir"]] * L_KM >= RP_MIN_KM - 1.0 and r["d_moon_km"] >= 5000.0
            ),
            "rp>=1.3floor & 5000<=d<=60000": lambda r: (
                sh.RP[r["ir"]] * L_KM >= RP_MIN_KM - 1.0 and 5000.0 <= r["d_moon_km"] <= 60000.0
            ),
            "rp>=1.3floor & 5000<=d<=20000": lambda r: (
                sh.RP[r["ir"]] * L_KM >= RP_MIN_KM - 1.0 and 5000.0 <= r["d_moon_km"] <= 20000.0
            ),
            "rp_floor_row (old 94%)": lambda r: r["ir"] == 0,
            "d<5000 (old 93%)": lambda r: r["d_moon_km"] < 5000.0,
        }
        res: dict[str, Any] = {}
        all_fams = {family_key(r) for r in recs if r.get("cycler_class_candidate")}
        for sname, fn in subsets.items():
            sel = [r for r in recs if fn(r)]
            st: Counter[str] = Counter()
            for r in sel:
                s = r["status"] if r["status"] in ("converged", "impact") else "other"
                if s == "impact":
                    s = "impact@1" if r["it"] == 1 else "impact@>1"
                st[s] += 1
            n = len(sel)
            cand = [r for r in sel if r.get("cycler_class_candidate")]
            fams = {family_key(r) for r in cand}
            res[sname] = {
                "n": n,
                "impact": round((st["impact@1"] + st["impact@>1"]) / n, 3) if n else None,
                "impact@1": round(st["impact@1"] / n, 3) if n else None,
                "impact@>1": round(st["impact@>1"] / n, 3) if n else None,
                "converged": round(st["converged"] / n, 3) if n else None,
                "candidates_per_solve": round(len(cand) / n, 4) if n else None,
                "n_candidates": len(cand),
                "n_families": len(fams),
                "families_lost_vs_all": sorted(all_fams - fams),
            }
        out[name] = res
        # band coverage: (p, q, sense, ir) combinations with no omega in the band
        nocover = [k for k, d in cache.items() if not np.any((d >= BAND_KM[0]) & (d <= BAND_KM[1]))]
        out[f"{name}_combos_without_band_omega"] = len(nocover)
        out[f"{name}_combos"] = len(cache)
        print(name, json.dumps(res, indent=1), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "predict.json").write_text(json.dumps(out, indent=1) + "\n")


def report() -> None:
    out: dict[str, Any] = {}
    for arm in ("B", "C", "D", "E"):
        f = OUT / f"arm_{arm}.jsonl"
        if f.exists():
            out[arm] = stats([json.loads(x) for x in f.read_text().splitlines()])
    # arm A: the old records for the same old sample
    old = {seed_id(sd) for sd in arms()["B"]}
    recs = [json.loads(x) for x in (SRC / "grid.jsonl").read_text().splitlines()]
    out["A (records, B's seeds)"] = stats([r for r in recs if seed_id(r) in old])
    for tag in ("control", "control_f1"):
        f = OUT / f"{tag}.jsonl"
        if f.exists():
            out[f"{tag}_7-3b"] = control_hits([json.loads(x) for x in f.read_text().splitlines()])
    (OUT / "report.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def control_hits(recs: list[dict[str, Any]]) -> dict[str, Any]:
    rid = "casoliva-7-3b-em-cycler-2010"
    per, t_row = sh.row_perigees(rid)
    hits = []
    for r in recs:
        if r["status"] != "converged":
            continue
        z = r["z0"]
        d = min(
            max(abs(z[0] - a), abs((z[1] - b + math.pi) % (2 * math.pi) - math.pi)) for a, b in per
        )
        if d < 1e-6 and abs(r["T"] - t_row) / t_row < 1e-6:
            hits.append(seed_id(r))
    return {
        "n": len(recs),
        "converged": sum(r["status"] == "converged" for r in recs),
        "hits": len(hits),
        "old_control": "22 of 1,152 seeds",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("predict")
    sub.add_parser("report")
    sp = sub.add_parser("arm")
    sp.add_argument("arm", choices=["B", "C", "D", "E"])
    sp.add_argument("--max-seconds", type=float, default=420.0)
    sp = sub.add_parser("control")
    sp.add_argument("--max-seconds", type=float, default=420.0)
    sp.add_argument("--design", choices=["f1f2f3", "f1"], default="f1f2f3")
    args = ap.parse_args()
    preflight_search(
        task_no=1054,
        region_id=f"1000-seed-design-{args.cmd}",
        method=MethodCapability(
            genome="#1000 p:q skeleton seeds, seed-design fix",
            corrector="none (record analysis)" if args.cmd == "predict" else "multiple shooting",
            capability_tags=frozenset({"cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    if args.cmd == "predict":
        predict()
    elif args.cmd == "report":
        report()
    elif args.cmd == "arm":
        drop = None if args.arm in ("C", "E") else DROP
        run_jobs(f"arm_{args.arm}", [(sd, drop, None) for sd in arms()[args.arm]], args.max_seconds)
    elif args.cmd == "control" and args.design == "f1":
        # F1 alone: F1 rows x the uniform 10-degree omega grid, both senses, no F3
        c73b = sh.CONTROLS["casoliva-7-3b-em-cycler-2010"]
        seeds = [
            {"p": 7, "q": 3, "sense": sense, "irn": irn, "k": iw, "om": sh.OMEGA[iw]}
            for sense in (-1, 1)
            for irn in range(len(RP_NEW))
            for iw in range(len(sh.OMEGA))
        ]
        run_jobs("control_f1", [(sd, None, c73b) for sd in seeds], args.max_seconds)
    elif args.cmd == "control":
        c73b = sh.CONTROLS["casoliva-7-3b-em-cycler-2010"]
        seeds = [sd for sd in new_seeds()[0] if (sd["p"], sd["q"]) == (7, 3)]
        # retrograde first: the 7-3b row is retrograde (all 22 old control hits were); order only
        seeds.sort(key=lambda sd: sd["sense"])
        run_jobs("control", [(sd, DROP, c73b) for sd in seeds], args.max_seconds)


if __name__ == "__main__":
    main()
