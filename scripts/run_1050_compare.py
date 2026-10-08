"""#1050: compare the upper ends of #1000 families A and D with Liang, Xu & Xu 2017's 5:2 and 7:3
polygonal-like orbits (match rule in docs/notes/2026-10-08-1050-continue-a-d-toward-liang-2017.md
sec. 0). Osculating Earth-centred elements at each perigee of a member: a, e, and omega = the
rotating-frame angle of the perigee. Output: data/1000_complement/liang2017_compare.json.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base

REPO = Path(__file__).resolve().parent.parent
CONT = REPO / "data" / "1000_complement" / "continuation"
GM = 1.0 - base.MU
LIANG = {
    "5:2": {"a": 0.5620, "e": 0.4660, "omega": 1.3634, "C": 3.0996, "p": 5, "q": 2},
    "7:3": {"a": 0.5604, "e": 0.3430, "omega": 4.1493, "C": 3.1858, "p": 7, "q": 3},
}


def osc(s: np.ndarray) -> tuple[float, float, float]:
    """Earth-centred osculating a, e and the perigee angle (rotating frame) at state s."""
    r = np.array([s[0] + base.MU, s[1]])
    v_in = np.array([s[2] - s[1], s[3] + s[0] + base.MU])  # + omega x (r - r_E)
    rn = float(np.linalg.norm(r))
    a = 1.0 / (2.0 / rn - float(v_in @ v_in) / GM)
    h = r[0] * v_in[1] - r[1] * v_in[0]
    e = math.sqrt(max(0.0, 1.0 - h * h / (GM * a)))
    return a, e, math.atan2(r[1], r[0]) % (2 * math.pi)


def perigees(s0: np.ndarray, period: float) -> list[tuple[float, float, float]]:
    sol = solve_ivp(
        base.eom,
        (0.0, period * 1.0001),
        s0,
        args=(base.MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=base._dr1,
    )
    return [osc(np.asarray(y)) for y in sol.y_events[0]]


def liang_state(el: dict[str, float]) -> np.ndarray:
    a, e, om = el["a"], el["e"], el["omega"]
    rp = a * (1 - e)
    vp = math.sqrt(GM * (1 + e) / rp)  # periapsis speed (prograde)
    x, y = rp * math.cos(om), rp * math.sin(om)
    vx_in, vy_in = -vp * math.sin(om), vp * math.cos(om)
    return np.array([x - base.MU, y, vx_in + y, vy_in - x])  # minus omega x (r - r_E)


def main() -> None:
    preflight_search(
        task_no=1050,
        region_id="liang2017-plpo-vs-1000-families-A-D",
        method=MethodCapability(
            genome="#1000 families A and D upper ends vs Liang 2017 PLPOs",
            corrector="none (comparison)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=2,
    )
    out: dict[str, Any] = {}
    for fam, key, ext in (
        ("A", "5:2", "5_2_pro_wE3_wM0_p_ext.jsonl"),
        ("D", "7:3", "7_3_pro_wE4_wM-1_p_ext.jsonl"),
    ):
        rows = [json.loads(line) for line in (CONT / ext).read_text().splitlines()]
        members = [r for r in rows if "s_perigee" in r]
        top = max(members, key=lambda r: r["C"])
        per = perigees(np.array(top["s_perigee"]), top["T"])
        el = LIANG[key]
        ls = liang_state(el)
        c_l = base.omega_eff(ls[0], ls[1]) - ls[2] ** 2 - ls[3] ** 2
        lp = perigees(ls, 2 * math.pi * el["q"])
        best = min(
            per,
            key=lambda x: abs(x[0] / el["a"] - 1) + abs(x[1] / el["e"] - 1),
        )
        match = (
            abs(top["C"] - el["C"]) < 0.01
            and abs(best[0] / el["a"] - 1) < 0.02
            and abs(best[1] / el["e"] - 1) < 0.02
        )
        out[fam] = {
            "family_top_C": top["C"],
            "family_top_T": top["T"],
            "family_stop": rows[-1].get("stop"),
            "family_top_perigees_a_e_omega": per,
            "liang": el,
            "liang_state_C_check": c_l,
            "liang_state_perigees_a_e_omega_over_2pi_q": lp,
            "closest_family_perigee": best,
            "match": bool(match),
        }
        print(fam, json.dumps({k: v for k, v in out[fam].items() if "perigees" not in k}))
    (REPO / "data" / "1000_complement" / "liang2017_compare.json").write_text(
        json.dumps(out, indent=1) + "\n"
    )


if __name__ == "__main__":
    main()
