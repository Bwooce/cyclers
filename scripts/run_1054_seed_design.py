"""#1054: the #1000 seed-design fix (perigee at least 1.3 times the floor, omega from the target
lunar-pass distance, a perigee step limit). Pre-registration:
docs/notes/2026-10-11-1054-seed-design-fix.md.

    uv run python scripts/run_1054_seed_design.py predict   # old records' rates per seed subset

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
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

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


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("predict")
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


if __name__ == "__main__":
    main()
