"""#1007 EXPLORATORY pre-check of the DRAFT same-phase drift metric (sec. 8.4), awaiting the owner.

Not a spec or catalogue change. For each chain: D'_i = max_j |t_(i+L)j - t_ij - L P_mean| and
W'_i = max_j | |V_inf,dep|_(i+L)j - |V_inf,dep|_ij |, i = 0 .. n-L-1; band D' <= 5 % of the chain's
shortest Lambert leg, W' <= 10 % of its smallest |V_inf,dep|; no growth (second-half max <= 2 x
first-half max, floors 0.5 d / 0.05 km/s). Junction-only continuity is reported beside it.

Chains: gc-1 and gc-2 (existing 10-cycle jup365 chains, L = 4); H&M 1970 orbit 1 from its printed
dates over two H&M periods (10 cycles, Standish, direct LM; L = 5); the broken-date negative (gc-1
epoch 0, 2nd Lambert start + 0.05 d).

Usage: uv run python scripts/v2_samephase_1007.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.search.hollister_menning_1970 import load_table3
from cyclerfinder.search.two_working_body import Cycle, LambertLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
OUT = REPO / "data" / "1007_v2" / "samephase_exploratory.json"
LABEL = "exploratory, under the draft metric, awaiting the owner"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


V2 = _load("v2_chain_1007", REPO / "scripts" / "v2_chain_1007.py")
CHAIN = V2.CHAIN


def legs_data(
    cell: str, key: str, n: int, mode: str, x0: float, y: list[float], real: Any
) -> tuple[np.ndarray, list[float], list[float], float]:
    """(Lambert start dates in days, |V_inf,dep| per Lambert leg, flight times in days,
    junction-only continuity in km/s)."""
    circ, a, b = CHAIN.ENUM.cell_system(cell)
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    legs_chain = one.legs * n
    nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
    yy = np.asarray(y, dtype=float)
    x = np.concatenate([[x0], yy[: n * nl - 1]])
    if mode == "ramp":
        cyc = Cycle(legs_chain, float(yy[n * nl - 1]) * DAY)
        lam = eval_lambert_legs(real, cyc, x * DAY)
        assert lam is not None
        res = CHAIN.date_chain_residual(real, legs_chain, x0 * DAY, yy)
        return (
            x,
            [float(np.linalg.norm(lg.vinf_dep)) for lg in lam],
            [(lg.t_arr - lg.t_dep) / DAY for lg in lam],
            float(np.max(np.abs(res))),
        )
    sysm = CHAIN.Blend(circ, real, np.eye(3), 0.0, 0.0, 1.0)
    ev = CHAIN.chain_eval(sysm, legs_chain, x0 * DAY, yy)
    assert ev is not None
    mags, q = [], 0
    for blk in ev.block_flybys:
        q += 3 * (len(blk) - 1)
        mags.append(abs(float(ev.residual[q])))
        q += 1
    vdep, tof = [], []
    for frm, t0, v0, _to, t1 in ev.segments[: n * nl]:
        _, w = real.state(frm, t0)
        vdep.append(float(np.linalg.norm(np.asarray(v0) - w)))
        tof.append((t1 - t0) / DAY)
    return x, vdep, tof, max(mags)


def same_phase(label: str, x: np.ndarray, vdep: list[float], tof: list[float], cont: float,
               n: int, nl: int, period_d: float, big_l: int) -> dict[str, Any]:  # fmt: skip
    p_mean = period_d / n
    m = n - big_l
    d = [
        max(abs(x[(i + big_l) * nl + j] - x[i * nl + j] - big_l * p_mean) for j in range(nl))
        for i in range(m)
    ]
    w = [
        max(abs(vdep[(i + big_l) * nl + j] - vdep[i * nl + j]) for j in range(nl)) for i in range(m)
    ]
    band_d, band_w = 0.05 * min(tof), 0.10 * min(vdep)
    h = max(1, m // 2)

    def growth(v: list[float], floor: float) -> tuple[bool, float]:
        first, second = max(v[:h]), max(v[h:]) if v[h:] else 0.0
        r = second / first if first > 0 else (float("inf") if second > 0 else 0.0)
        return second <= floor or second <= 2.0 * first, r

    gd_ok, gd = growth(d, 0.5)
    gw_ok, gw = growth(w, 0.05)
    ok = max(d) <= band_d and max(w) <= band_w and gd_ok and gw_ok
    out = {
        "label": label, "status": LABEL, "L": big_l, "cycles": n, "comparisons": m,
        "junction_continuity_kms": cont, "Dprime_d": d, "Wprime_kms": w,
        "band_D_d": band_d, "band_W_kms": band_w, "growth_D": gd, "growth_W": gw,
        "drift_pass": bool(ok), "continuity_pass": bool(cont < 1e-6),
    }  # fmt: skip
    print(
        f"{label}: L={big_l} n={n} cont={cont:.1e} D'max={max(d):.3f}/{band_d:.3f} d "
        f"W'max={max(w):.4f}/{band_w:.3f} km/s growth D {gd:.2f} W {gw:.2f} -> drift "
        f"{'PASS' if ok else 'FAIL'}, continuity {'PASS' if cont < 1e-6 else 'FAIL'} [{LABEL}]",
        flush=True,
    )
    return out


def main() -> None:
    res: dict[str, Any] = {"status": LABEL}
    real_gc = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    rows = []
    g1 = wb["gc-1"]
    for i in range(5):
        e = json.loads(
            (REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json").read_text()
        )[0]
        x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], e["final_y"], real_gc)
        rows.append(same_phase(f"gc-1 JD {e['epoch_jd']:.1f}", x, v, t, c, 10, 3,
                               e["final_y"][29], 4))  # fmt: skip
    res["gc-1"] = rows
    g2 = wb["gc-2"]
    rows = []
    ladder = json.loads((REPO / "data" / "942_943_ladder_realeph_chains.json").read_text())["gc2"]
    for e in ladder:
        xs = e["final_x_days"]
        y = list(xs[1:])
        x, v, t, c = legs_data("gc", g2["key"], 10, "blend", xs[0], y, real_gc)
        rows.append(same_phase(f"gc-2 JD {e['epoch_jd']:.1f}", x, v, t, c, 10, 4, y[39], 4))
    res["gc-2"] = rows
    # broken-date negative
    e = json.loads((REPO / "data" / "943_gc1_realeph" / "e0" / "realeph_chain.json").read_text())[0]
    y_bad = list(e["final_y"])
    y_bad[0] += 0.05
    x, v, t, c = legs_data("gc", g1["key"], 10, "blend", e["x0_days"], y_bad, real_gc)
    res["negative"] = same_phase("NEGATIVE gc-1 e0 broken", x, v, t, c, 10, 3, y_bad[29], 4)
    # H&M orbit 1, printed dates, two periods
    real_ev = CHAIN.real_ephemeris("ev", "mean")
    key = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
    t3 = load_table3()[1]
    one_period = []
    for k in range(5):
        one_period += [t3[5 * k + 1].date, t3[5 * k + 4].date]
    starts = one_period + [s + 5844.0 for s in one_period]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], 2 * 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    legs = one.legs * 10
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real_ev, legs, x0 * DAY, yy),
        y_seed, method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200 * len(y_seed),
    )  # fmt: skip
    x, v, t, c = legs_data("ev", key, 10, "ramp", x0, list(sol.x), real_ev)
    hm = same_phase("POSITIVE H&M orbit 1 (printed, 10 cycles)", x, v, t, c, 10, 2, sol.x[-1], 5)
    hm["max_shift_from_print_d"] = float(np.max(np.abs(sol.x[:-1] - y_seed[:-1])))
    hm["chain_period_d"] = float(sol.x[-1])
    print(f"H&M: max shift from print {hm['max_shift_from_print_d']:.3f} d, chain period "
          f"{hm['chain_period_d']:.3f} d (2 x 5844 = 11688)", flush=True)  # fmt: skip
    res["hm_orbit1"] = hm
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print("wrote", OUT.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
