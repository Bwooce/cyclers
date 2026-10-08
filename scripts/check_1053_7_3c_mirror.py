"""#1053: re-derive casoliva-7-3c as the EXACT mirror image (x, -y, -t) of the casoliva-7-3b row
(lead ruling: option (a)), and record its match to Casoliva 2010 Table 3's printed 7-3c values.

The mirror of a planar CR3BP state (x, y, xdot, ydot) under (x, y, t) -> (x, -y, -t) is
(x, -y, -xdot, ydot); the mirrored orbit has the same C and period. Output:
data/1053_casoliva_7-3c_rederive/rederive.json.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, cr3bp_system, jacobi_constant
from cyclerfinder.search.earth_moon_resonant_families import (
    TABLE3_ROWS,
    planar_stability_index,
    table3_seed_state,
)

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1053_casoliva_7-3c_rederive" / "rederive.json"


def main() -> None:
    sysm = cr3bp_system("Earth", "Moon")
    with (REPO / "data" / "catalogue.yaml").open() as fh:
        rows = {r["id"]: r for r in yaml.safe_load(fh)}
    b = rows["casoliva-7-3b-em-cycler-2010"]["orbit_elements"]["cr3bp"]
    c_old = rows["casoliva-7-3c-em-cycler-2010"]["orbit_elements"]["cr3bp"]
    sb = np.array(b["state_nd"], dtype=float)
    s_new = np.array([sb[0], -sb[1], sb[2], -sb[3], sb[4], -sb[5]])
    s_new[2] = 0.0
    s_new[5] = 0.0
    period = float(b["period_nd"])
    mu = sysm.mu
    jc = jacobi_constant(s_new, mu)

    def ey(t: float, y: np.ndarray, m: float) -> float:
        return float(y[1])

    def d1(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] + m) * y[3] + y[1] * y[4])

    def d2(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] - 1 + m) * y[3] + y[1] * y[4])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, period),
        s_new,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=[ey, d1, d2],
    )
    closure = float(np.linalg.norm(sol.y[:, -1] - s_new))
    r1 = [math.hypot(y[0] + mu, y[1]) for y in sol.y_events[1]]
    r2 = [math.hypot(y[0] - 1 + mu, y[1]) for y in sol.y_events[2]]
    printed = next(r for r in TABLE3_ROWS if r.designation == "7-3c")
    ps = table3_seed_state(printed)  # project frame, Casoliva's printed 7-3c crossing (y = 0)
    crossings = [np.asarray(y) for y in sol.y_events[0]]
    best = min(crossings, key=lambda y: float(np.linalg.norm(y[[0, 3, 4]] - ps[[0, 3, 4]])))
    # store the row at the printed section crossing of the mirrored orbit (same orbit, re-phased)
    s_row = best.copy()
    s_row[1] = 0.0
    s_row[2] = 0.0
    s_row[5] = 0.0
    sol_r = solve_ivp(
        cr3bp_eom, (0.0, period), s_row, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
    )
    closure_row = float(np.linalg.norm(sol_r.y[:, -1] - s_row))
    stab = planar_stability_index(sysm, s_row, period, designation="7-3c")
    stab_b = planar_stability_index(sysm, sb, period, designation="7-3b")
    out = {
        "state_nd_mirror_of_7_3b_start": s_new.tolist(),
        "state_nd_new_row": s_row.tolist(),
        "closure_row_state": closure_row,
        "jacobi_row_state": jacobi_constant(s_row, mu),
        "jacobi_constant_new": jc,
        "jacobi_7_3b_row": b["jacobi_constant"],
        "period_nd_new": period,
        "tof_days_new": period * sysm.t_s / 86400.0,
        "closure_full_period": closure,
        "periselene_km": min(r2) * sysm.l_km,
        "perigee_km": min(r1) * sysm.l_km,
        "apogee_km": max(r1) * sysm.l_km,
        "printed_7_3c": {
            "C": printed.c_j,
            "T": printed.period,
            "x_i": printed.x_i,
            "u_i": printed.u_i,
            "v_i": printed.v_i,
            "k": printed.k,
            "r_pm": printed.r_pm,
            "r_pe": printed.r_pe,
            "r_ae": printed.r_ae,
        },
        "nearest_crossing_to_printed": best.tolist(),
        "crossing_minus_printed_x_vx_vy": (best[[0, 3, 4]] - ps[[0, 3, 4]]).tolist(),
        "rel_dC_printed": (jc - printed.c_j) / printed.c_j,
        "rel_dT_printed": (period - printed.period) / printed.period,
        "rel_dx0_printed": float((best[0] - ps[0]) / ps[0]),
        "k_par": stab.k_par,
        "k_perp": stab.k_perp,
        "k_signed_casoliva": stab.k_signed,
        "k_par_7_3b_row": stab_b.k_par,
        "k_perp_7_3b_row": stab_b.k_perp,
        "rel_dk_printed": (stab.k_signed - printed.k) / printed.k,
        "old_row": {
            "state_nd": c_old["state_nd"],
            "jacobi_constant": c_old["jacobi_constant"],
            "period_nd": c_old["period_nd"],
            "stability_index": c_old["stability_index"],
        },
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    for k, v in out.items():
        print(k, v)


if __name__ == "__main__":
    main()
