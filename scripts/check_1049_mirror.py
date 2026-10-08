"""#1049: are the catalogue's casoliva-7-3b and casoliva-7-3c rows one orbit and its mirror?
Each row is tested, directly and mirrored, against the #1000 family-C continuation branch solved at
the row's own C, and the 7-3c row is compared with the mirror of the 7-3b row.

Solves the #1000 family C (continuation members in
data/1000_complement/continuation/7_3_ret_wE-10_wM-1_p.jsonl) at the 7-3c row's C, then compares
that member's perigee section points with the 7-3c row and with its mirror image
(x, y, xdot, ydot) -> (x, -y, -xdot, ydot) under time reversal.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base
import run_1000_continue as cn
import run_1000_shooting as sh

REPO = Path(__file__).resolve().parent.parent


def perigees(s0: np.ndarray, period: float) -> list[tuple[float, float]]:
    sol = solve_ivp(
        base.eom,
        (0.0, period),
        s0,
        args=(base.MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=base._dr1,
    )
    return [base.to_section(y) for y in sol.y_events[0]]


def distance(a_set: list[tuple[float, float]], b_set: list[tuple[float, float]]) -> float:
    return min(
        max(abs(math.log(a[0] / b[0])), abs((a[1] - b[1] + math.pi) % (2 * math.pi) - math.pi))
        for a in a_set
        for b in b_set
    )


def main() -> None:
    with (REPO / "data" / "catalogue.yaml").open() as fh:
        rows = {r["id"]: r for r in yaml.safe_load(fh)}
    fam_file = REPO / "data" / "1000_complement" / "continuation" / "7_3_ret_wE-10_wM-1_p.jsonl"
    with fam_file.open() as fh:
        fam = [json.loads(line) for line in fh if "s_perigee" in line]
    sets = {}
    for rid in ("casoliva-7-3b-em-cycler-2010", "casoliva-7-3c-em-cycler-2010"):
        c = rows[rid]["orbit_elements"]["cr3bp"]
        s6 = np.array(c["state_nd"])
        s = np.array([s6[0], s6[1], s6[3], s6[4]])
        mirror = np.array([s[0], -s[1], -s[2], s[3]])
        near = min(fam, key=lambda m: abs(m["C"] - c["jacobi_constant"]))
        nodes, taus = cn.nodes_from(np.array(near["s_perigee"]), near["T"], 7)
        res = sh.shoot(nodes, taus, c["jacobi_constant"], step_max=0.01, max_it=60)
        member = perigees(np.array(res["nodes"][0]), res["T"])
        direct = distance(perigees(s, c["period_nd"]), member)
        mirrored = distance(perigees(mirror, c["period_nd"]), member)
        print(
            f"{rid}: C {c['jacobi_constant']:.9f}; vs the family-C branch member at this C: "
            f"direct {direct:.2e}, mirrored {mirrored:.2e}"
        )
        sets[rid] = (s, mirror, c["period_nd"])
    b, c7 = sets["casoliva-7-3b-em-cycler-2010"], sets["casoliva-7-3c-em-cycler-2010"]
    c_set = perigees(c7[0], c7[2])
    d_mir = distance(c_set, perigees(b[1], b[2]))
    d_dir = distance(c_set, perigees(b[0], b[2]))
    print(f"7-3c row vs mirror(7-3b row): {d_mir:.2e}; 7-3c row vs 7-3b row: {d_dir:.2e}")


if __name__ == "__main__":
    main()
