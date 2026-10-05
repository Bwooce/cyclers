"""#942 ladder step 3: Hollister 1969 Table 1 V_inf directions control (orbit I).

Pre-registration: docs/notes/2026-10-05-942-943-two-working-body-generator.md sec. 3.7.
Expected values: Hollister, JSR 6(4) 1969, p.368, Table 1 (read from the page image).

Usage: uv run python scripts/run_942_hollister1969_directions.py --out <dir>
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.hollister_menning_1970 import build_cycle, load_table3, pick_branches
from cyclerfinder.search.two_working_body import (
    EMOS_KMS,
    MeanElementSystem,
    correct_dates,
    eval_lambert_legs,
)

DAY = SECONDS_PER_DAY
#: (event, body, JD-2440000, V EMOS, Ang deg, Elev deg), Hollister 1969 Table 1.
TABLE1 = [
    ("LV", "E", 806, 0.154, 163, -52),
    ("AR", "V", 971, 0.178, 38, 55),
    ("LV", "V", 1421, 0.178, 335, -60),
    ("AR", "E", 1592, 0.154, 201, 50),
    ("LV", "E", 1957, 0.154, 143, -44),
    ("AR", "V", 2125, 0.203, 8, 67),
    ("LV", "V", 2575, 0.203, 294, 14),
    ("AR", "E", 2797, 0.191, 240, 12),
    ("LV", "E", 3163, 0.191, 186, 58),
    ("AR", "V", 3316, 0.194, 36, -62),
    ("LV", "V", 3765, 0.194, 343, 65),
    ("AR", "E", 3935, 0.156, 215, 46),
    ("LV", "E", 4300, 0.156, 154, 50),
    ("AR", "V", 4471, 0.192, 17, -64),
    ("LV", "V", 4921, 0.192, 319, 58),
    ("AR", "E", 5077, 0.174, 183, -56),
    ("LV", "E", 5442, 0.174, 123, -12),
    ("AR", "V", 5664, 0.223, 69, -15),
    ("LV", "V", 6114, 0.223, 8, -69),
    ("AR", "E", 6285, 0.154, 215, 47),
    ("LV", "E", 6650, 0.154, 163, -52),
]
TOL = {"v": 0.005, "ang": 10.0, "elev": 10.0}


def ang_elev(
    system: MeanElementSystem, body: str, t_s: float, v: np.ndarray, sign: int
) -> tuple[float, float]:
    r, w = system.state(body, t_s)
    n = np.cross(r, w)
    n /= np.linalg.norm(n)
    r_hat = r / np.linalg.norm(r)
    c_hat = np.cross(n, r_hat)
    vin = v - (v @ n) * n
    ang = math.degrees(math.atan2(sign * float(vin @ r_hat), float(vin @ c_hat))) % 360.0
    elev = math.degrees(math.asin(float(v @ n) / float(np.linalg.norm(v))))
    return ang, elev


def evaluate(system: MeanElementSystem, x: np.ndarray, cyc, sign: int) -> list[dict]:
    legs = eval_lambert_legs(system, cyc, x)
    assert legs is not None
    ends = []
    for ev, li in zip(legs, cyc.lambert_index, strict=True):
        leg = cyc.legs[li]
        ends.append(("LV", leg.frm, ev.t_dep, ev.vinf_dep))
        ends.append(("AR", leg.to, ev.t_arr, ev.vinf_arr))
    out = []
    for kind, body, date, v_pr, a_pr, e_pr in TABLE1:
        cands = [e for e in ends if e[0] == kind and e[1] == body]
        best = min(
            cands, key=lambda e: min(abs(e[2] / DAY - date + k * 5844.0) for k in (-1, 0, 1))
        )
        dd = min(abs(best[2] / DAY - date + k * 5844.0) for k in (-1, 0, 1))
        a, el = ang_elev(system, body, best[2], best[3], sign)
        da = (a - a_pr + 180.0) % 360.0 - 180.0
        rec = {
            "event": kind,
            "body": body,
            "date": date,
            "ddate": dd,
            "v_pr": v_pr,
            "v": float(np.linalg.norm(best[3])) / EMOS_KMS,
            "ang_pr": a_pr,
            "ang": a,
            "dang": da,
            "elev_pr": e_pr,
            "elev": el,
            "delev": el - e_pr,
        }
        rec["ok"] = bool(
            abs(rec["v"] - v_pr) <= TOL["v"]
            and abs(da) <= TOL["ang"]
            and abs(rec["delev"]) <= TOL["elev"]
        )
        out.append(rec)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    preflight_search(
        task_no=942,
        region_id="hollister-1969-table1-direction-control",
        method=MethodCapability(
            genome="Hollister periodic orbit I (Table 3 orbit 1 structure)",
            corrector="two_working_body.correct_dates",
            capability_tags=frozenset({"ballistic", "patched-conic", "3d", "inclined-elliptic"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=4,
    )
    rows = load_table3()[1]
    result = {}
    for name, system in (
        (
            "exact_periodicity_SOURCED",
            MeanElementSystem(
                periods_days={"E": 5844.0 / 16, "V": 5844.0 / 26}, anchor_jd=2443363.0
            ),
        ),
        ("real_periods", MeanElementSystem()),
    ):
        br, _ = pick_branches(system, 1, rows)
        cyc, x0 = build_cycle(system, 1, rows, br)
        sol = correct_dates(system, cyc, x0)
        for where, x in (("converged", sol.x), ("printed_dates", x0)):
            for conv, sign in (("A_toward_outward_radial", 1), ("B_toward_inward_radial", -1)):
                ev = evaluate(system, x, cyc, sign)
                frac = float(np.mean([e["ok"] for e in ev]))
                key = f"{name}|{where}|{conv}"
                result[key] = {
                    "frac_ok": frac,
                    "pass": frac >= 0.9,
                    "converged": sol.converged,
                    "events": ev,
                }
                med_a = float(np.median([abs(e["dang"]) for e in ev]))
                med_e = float(np.median([abs(e["delev"]) for e in ev]))
                print(
                    f"{key}: frac_ok={frac:.2f} pass={frac >= 0.9} "
                    f"median|dAng|={med_a:.1f} median|dElev|={med_e:.1f}",
                    flush=True,
                )
    (args.out / "hollister1969_directions.json").write_text(
        json.dumps(result, indent=1, default=float)
    )


if __name__ == "__main__":
    main()
