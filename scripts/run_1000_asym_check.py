"""#1000: check that families A, B, D (and C) are really asymmetric, arc by arc from converged
multiple-shooting nodes (single shooting over a whole period drifts on these orbits), and give
b from the product of the per-arc STMs for every asymmetric grid member.
Output: data/1000_complement/asym_check.json.
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base
import run_1000_continue as cn
import run_1000_shooting as sh

REPO = Path(__file__).resolve().parent.parent
D = REPO / "data" / "1000_complement"


def arc_min_xdot(nodes: list[list[float]], taus: list[float]) -> float:
    def ey(t: float, s: np.ndarray, mu: float) -> float:
        return float(s[1])

    best = float("inf")
    for s, t in zip(nodes, taus, strict=True):
        sol = solve_ivp(
            base.eom,
            (0.0, t),
            np.array(s),
            args=(base.MU,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
            events=ey,
        )
        for y in sol.y_events[0]:
            best = min(best, abs(float(y[2])))
        if abs(s[1]) < 1e-12:
            best = min(best, abs(float(s[2])))
    return best


def check(s_per: list[float], period: float, c: float, p: int) -> dict[str, Any]:
    nodes, taus = cn.nodes_from(np.array(s_per), period, p)
    r = sh.shoot(nodes, taus, c, step_max=0.01, max_it=60)
    if not r["converged"]:
        return {"resolved": False}
    b, lam = cn.ms_monodromy(r["nodes"], r["taus"])
    return {
        "resolved": True,
        "min_xdot_y0": arc_min_xdot(r["nodes"], r["taus"]),
        "b_ms": b,
        "lambda_max_ms": lam,
        "C": c,
        "T": r["T"],
    }


def main() -> None:
    preflight_search(
        task_no=1000,
        region_id="em-complement-asymmetry-check",
        method=MethodCapability(
            genome="asymmetric family members re-solved by multiple shooting",
            corrector="run_1000_shooting.shoot at fixed C",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar", "asymmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=150,
    )
    out: dict[str, Any] = {"grid": [], "continuation": {}}
    grid = [json.loads(line) for line in (D / "shooting" / "grid.jsonl").read_text().splitlines()]
    for r in grid:
        if r["status"] == "converged" and r["cycler_class_candidate"] and not r["symmetric"]:
            sense = "pro" if r["sense"] > 0 else "ret"
            key = f"{r['p']}_{r['q']}_{sense}_wE{round(r['wind_E'])}_wM{round(r['wind_M'])}"
            res = check(r["s_perigee"], r["T"], r["C"], r["p"])
            res.update({"family": key, "C_grid": r["C"], "b_single_shooting": r["b_h"]})
            out["grid"].append(res)
            print(key, res, flush=True)
    for f in sorted(glob.glob(str(D / "continuation" / "*.jsonl"))):
        rows = [json.loads(line) for line in Path(f).read_text().splitlines()]
        rows = [x for x in rows if ("s_perigee" in x and "stop" not in x) or "s_perigee" in x]
        p = int(Path(f).name.split("_")[0])
        mins = []
        for x in rows[::5]:
            res = check(x["s_perigee"], x["T"], x["C"], p)
            if res["resolved"]:
                mins.append(res["min_xdot_y0"])
        out["continuation"][Path(f).name] = {
            "n_checked": len(mins),
            "min_xdot_y0_min": min(mins) if mins else None,
            "min_xdot_y0_max": max(mins) if mins else None,
        }
        print(Path(f).name, out["continuation"][Path(f).name], flush=True)
    (D / "asym_check.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
