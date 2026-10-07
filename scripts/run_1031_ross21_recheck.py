"""#1031: re-check the stability claim of ross-rt-em-cycler-21-2025 (Ross & Roberts-Tsoukkas
2025, AAS 25-621, Table 3: (2,1) C^stable = 3.129389531088256, C^max = 3.129389531092325,
T^stable = 19.44043166795154 TU, stable-subfamily perilune width Delta_p_m = 4.23 km).

C^max - C^stable = 4.1e-12: the published stable window sits against a fold of the family in C.
A fixed-C corrector is ill-conditioned there and can land on either branch. So the family is
parametrised by x0 (the perpendicular crossing at the row's start, x0 > 0) with ydot0 solved by
Newton on xdot at the half-period crossing (fixed crossing index), as in
scripts/run_997_lineage.py. Along an x0 scan around the row's state this records C, T, the
perilune and perigee distances, the Barden index nu = k_par / 2 (in-plane) and k_perp (vertical)
from the full-period monodromy. It then locates nu = 0 and the |nu| = 1 window edges by
bisection in x0.

Model: the row's own mu = 0.012150584270572 (Ross p.3), L = 384,400 km, R_Moon = 1737.4 km.
Output: data/1031_ross21/scan.json.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_997_lineage as lin

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1031_ross21" / "scan.json"
MU = 1.2150584270572e-2
X0_ROW, YD0_ROW = 0.7237335857, 0.4137707374
C_STABLE, C_MAX, T_STABLE = 3.129389531088256, 3.129389531092325, 19.44043166795154


def member(x0: float, yd_guess: float, n: int, t_guess: float) -> dict[str, Any] | None:
    res = lin.correct_fixed_x0(MU, x0, yd_guess, n, t_guess)
    if res is None:
        return None
    yd0, t_half = res
    m = lin.member_report(MU, x0, yd0, t_half)
    m["b_v"] = lin.vertical_index(MU, x0, yd0, m["T"])
    m["nu"] = m["b_h"] / 2.0
    m["dC_from_C_stable"] = m["C"] - C_STABLE
    m["dC_from_C_max"] = m["C"] - C_MAX
    return m


def main() -> None:
    preflight_search(
        task_no=1031,
        region_id="ross-rt-em-cycler-21-stability-recheck",
        method=MethodCapability(
            genome="planar x-axis-symmetric (2,1) prograde EM cycler family, x0 scan",
            corrector="fixed-x0 perpendicular-crossing Newton (run_997_lineage)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=80,
    )
    # crossing index of the half period on the row's own state
    times = lin.crossing_times(MU, X0_ROW, YD0_ROW, 0.6 * T_STABLE)
    n = int(np.argmin(np.abs(np.array(times) - T_STABLE / 2))) + 1
    print(f"half-period crossing index n = {n} (t = {times[n - 1]:.5f})", flush=True)
    rows: list[dict[str, Any]] = []
    yd, th = YD0_ROW, times[n - 1]
    # scan outward from the row's x0 in both directions (continuation in x0)
    for sign in (+1.0, -1.0):
        yd, th = YD0_ROW, times[n - 1]
        for k in range(0 if sign > 0 else 1, 41):
            x0 = X0_ROW + sign * k * 5e-7
            m = member(x0, yd, n, th)
            if m is None:
                print(f"  x0={x0:.10f}: no closure; stop this direction", flush=True)
                break
            yd, th = m["ydot0"], m["t_half"]
            rows.append(m)
            print(
                f"  x0={x0:.10f} C={m['C']:.15f} T={m['T']:.8f} nu={m['nu']:+.4f} "
                f"k_perp={m['b_v']:+.4f} perilune={m['periselene_km']:.3f} km",
                flush=True,
            )
    rows.sort(key=lambda r: r["x0"])

    # bisection for nu = 0, nu = +1, nu = -1 crossings between scan neighbours
    def bisect(target: float) -> list[dict[str, Any]]:
        found = []
        for a, b in itertools.pairwise(rows):
            fa, fb = a["nu"] - target, b["nu"] - target
            if fa * fb > 0:
                continue
            lo, hi, mlo = a["x0"], b["x0"], a
            mm = a
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                mm = member(mid, mlo["ydot0"], n, mlo["t_half"]) or mm
                if (mm["nu"] - target) * fa > 0:
                    lo, mlo = mid, mm
                else:
                    hi = mid
                if hi - lo < 1e-12:
                    break
            found.append(
                {k: mm[k] for k in ("x0", "ydot0", "C", "T", "nu", "b_v", "periselene_km")}
            )
        return found

    edges = {f"nu={t:+.0f}": bisect(t) for t in (0.0, 1.0, -1.0)}
    for k, v in edges.items():
        for e in v:
            print(
                f"{k}: x0={e['x0']:.12f} C={e['C']:.15f} T={e['T']:.10f} nu={e['nu']:+.2e} "
                f"k_perp={e['b_v']:+.4f} perilune={e['periselene_km']:.4f} km",
                flush=True,
            )
    keep = (
        "x0",
        "ydot0",
        "C",
        "T",
        "nu",
        "b_h",
        "b_v",
        "periselene_km",
        "perigee_km",
        "dC_from_C_stable",
        "dC_from_C_max",
        "closure",
        "det_M4",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps(
            {
                "mu": MU,
                "row_state": [X0_ROW, YD0_ROW],
                "paper": {
                    "C_stable": C_STABLE,
                    "C_max": C_MAX,
                    "T_stable": T_STABLE,
                    "delta_p_m_km": 4.23,
                },
                "n_half_crossing": n,
                "scan": [{k: r[k] for k in keep} for r in rows],
                "edges": edges,
            },
            indent=1,
        )
        + "\n"
    )
    print(f"wrote {OUT}; perilune unit check: 1 nd = {lin.L_KM} km, R_M = {lin.R_M} km")


if __name__ == "__main__":
    main()
