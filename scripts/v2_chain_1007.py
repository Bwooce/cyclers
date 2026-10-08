"""#1007: V2-ballistic re-judgement of real-ephemeris chains (spec §14 amendment 2026-10-08).

Pre-registration: docs/notes/2026-10-08-1007-v2-ballistic-chain-campaign.md. For each stored chain
(no re-solve) at one epoch, the four criteria:
1. continuity: the chain residual at lambda = 1 is < 1e-6 (km/s; shot fixed legs also carry their
   arrival miss / 1000 km);
2. gate: the #888/#937 demanded-turn gate passes at every interior flyby;
3. re-fly: DOP853 arrival miss < 1 km on every leg (ramp: also V_inf vector error < 1e-6 km/s);
4. bounded drift: D_i (Lambert-start dates against a strictly periodic chain at the chain's own
   mean period) and W_i (Lambert departure |V_inf| against the ideal template) inside the band, with
   no growth (second-half max <= 2 x first-half max, floors 0.5 d and 0.05 km/s).

Usage: uv run python scripts/v2_chain_1007.py --lane gc   (or --lane ev-hm, see the note)
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.search.two_working_body import (
    Cycle,
    LambertLeg,
    cycle_flybys,
    eval_lambert_legs,
)

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
OUT = REPO / "data" / "1007_v2"

CONT_TOL = 1e-6
REFLY_KM = 1.0
VINF_VEC_TOL = 1e-6
BAND_D_FRAC = 0.05
BAND_W_FRAC = 0.10
GROWTH = 2.0
FLOOR_D_DAYS = 0.5
FLOOR_W_KMS = 0.05


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CHAIN = _load("run_942_realeph_chain", REPO / "scripts" / "run_942_realeph_chain.py")
GAUNTLET = _load("gauntlet_942", REPO / "scripts" / "gauntlet_942.py")


def fly(mu: float, r0: np.ndarray, v0: np.ndarray, dt: float) -> np.ndarray:
    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        return np.concatenate([y[3:], -mu * y[:3] / np.linalg.norm(y[:3]) ** 3])

    sol = solve_ivp(
        rhs, (0.0, dt), np.concatenate([r0, v0]), method="DOP853", rtol=1e-13, atol=1e-8
    )
    return np.asarray(sol.y[:3, -1])


def template(cell: str, key: str, x_days: list[float]) -> dict[str, Any]:
    circ, a, b = CHAIN.ENUM.cell_system(cell)
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    legs = eval_lambert_legs(circ, one, np.asarray(x_days) * DAY)
    assert legs is not None
    return {
        "circ": circ,
        "one": one,
        "vinf_dep": [float(np.linalg.norm(lg.vinf_dep)) for lg in legs],
        "min_tof_d": min((lg.t_arr - lg.t_dep) / DAY for lg in legs),
    }


def _growth_ok(vals: list[float], floor: float) -> tuple[bool, float]:
    """Cycles 1..n-1: second-half max <= GROWTH x first-half max (growth below floor ignored)."""
    v = vals[1:]
    h = max(1, len(v) // 2)
    first, second = max(v[:h]), max(v[h:]) if v[h:] else 0.0
    ratio = second / first if first > 0 else float("inf") if second > 0 else 0.0
    return (second <= floor or second <= GROWTH * first), ratio


def judge(
    cell: str,
    key: str,
    x_days: list[float],
    n: int,
    mode: str,
    x0_days: float,
    y: list[float],
    real: Any,
    label: str,
) -> dict[str, Any]:
    tpl = template(cell, key, x_days)
    one = tpl["one"]
    legs_chain = one.legs * n
    nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
    yy = np.asarray(y, dtype=float)
    x = np.concatenate([[x0_days], yy[: n * nl - 1]]) * DAY
    period = float(yy[n * nl - 1]) * DAY
    out: dict[str, Any] = {"label": label, "mode": mode, "cycles": n}
    if mode == "ramp":
        cyc = Cycle(legs_chain, period)
        res = CHAIN.date_chain_residual(real, legs_chain, x0_days * DAY, yy)
        gate = CHAIN.interior_gate(real, cyc, x)
        lam = eval_lambert_legs(real, cyc, x)
        flybys = cycle_flybys(real, cyc, x)
        if lam is None or flybys is None:
            out["error"] = "lambert-fail"
            out["pass"] = False
            return out
        cc = GAUNTLET.cross_check(real, cyc, x, flybys)
        refly_km = float(cc["max_arrival_miss_km"])
        vec_err = float(cc["max_vinf_vector_error_kms"])
        vdep = [float(np.linalg.norm(lg.vinf_dep)) for lg in lam]
    else:
        sysm = CHAIN.Blend(tpl["circ"], real, np.eye(3), 0.0, 0.0, 1.0)
        ev = CHAIN.chain_eval(sysm, legs_chain, x0_days * DAY, yy)
        if ev is None:
            out["error"] = "chain-eval-fail"
            out["pass"] = False
            return out
        res = ev.residual
        # the |V_inf| magnitude residual of each junction (the last entry of each block; the
        # shot fixed legs' 3 position residuals per leg come first), reported beside the
        # pre-registered mixed residual
        mags, q = [], 0
        for blk in ev.block_flybys:
            q += 3 * (len(blk) - 1)
            mags.append(abs(float(res[q])))
            q += 1
        out["junction_vinf_max_kms"] = max(mags)
        gate = CHAIN.gate_eval(sysm, ev)
        refly_km = 0.0
        for frm, t0, v0, to, t1 in ev.segments:
            r0, _ = real.state(frm, t0)
            rb, _ = real.state(to, t1)
            refly_km = max(refly_km, float(np.linalg.norm(fly(real.mu, r0, v0, t1 - t0) - rb)))
        vec_err = 0.0
        vdep = []
        for frm, t0, v0, _to, _t1 in ev.segments[: n * nl]:
            _, w = real.state(frm, t0)
            vdep.append(float(np.linalg.norm(np.asarray(v0) - w)))
    cont = float(np.max(np.abs(res)))
    p_mean = period / n
    d_i = [max(abs(x[i * nl + j] - x[j] - i * p_mean) / DAY for j in range(nl)) for i in range(n)]
    w_i = [max(abs(vdep[i * nl + j] - tpl["vinf_dep"][j]) for j in range(nl)) for i in range(n)]
    band_d = BAND_D_FRAC * tpl["min_tof_d"]
    band_w = BAND_W_FRAC * min(tpl["vinf_dep"])
    gd_ok, gd = _growth_ok(d_i, FLOOR_D_DAYS)
    gw_ok, gw = _growth_ok(w_i, FLOOR_W_KMS)
    c1 = cont < CONT_TOL
    c2 = gate.get("status") == "pass"
    c3 = refly_km < REFLY_KM and vec_err < VINF_VEC_TOL
    c4 = max(d_i) <= band_d and max(w_i) <= band_w and gd_ok and gw_ok
    out |= {
        "continuity_max": cont,
        "gate": {k: gate.get(k) for k in ("status", "worst_ratio")},
        "refly_max_km": refly_km,
        "vinf_vector_error_kms": vec_err,
        "D_max_d": max(d_i),
        "W_max_kms": max(w_i),
        "band_D_d": band_d,
        "band_W_kms": band_w,
        "growth_D": gd,
        "growth_W": gw,
        "D_by_cycle": d_i,
        "W_by_cycle": w_i,
        "c1_continuity": c1,
        "c2_gate": c2,
        "c3_refly": c3,
        "c4_drift": c4,
        "pass": bool(c1 and c2 and c3 and c4),
    }
    print(
        f"{label}: n={n} cont={cont:.1e} gate={gate.get('status')} "
        f"({gate.get('worst_ratio', float('nan')):.3f}) refly={refly_km:.2e} km "
        f"D={max(d_i):.3f}/{band_d:.3f} d W={max(w_i):.3f}/{band_w:.3f} km/s "
        f"growth D {gd:.2f} W {gw:.2f} -> {'PASS' if out['pass'] else 'FAIL'} "
        f"[{int(c1)}{int(c2)}{int(c3)}{int(c4)}]",
        flush=True,
    )
    return out


def _blend_records(path: Path) -> list[dict[str, Any]]:
    return [e for e in json.loads(path.read_text()) if e.get("final_y")]


def lane_gc() -> dict[str, Any]:
    real = CHAIN.real_ephemeris("gc", "spice")
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())
    r316 = next(
        r
        for r in json.loads((REPO / "data" / "943_ganeur316_recall.json").read_text())
        if r.get("status") == "pass"
    )
    out: dict[str, Any] = {}

    def run_blend(
        name: str, cell: str, key: str, xd: list[float], files: list[Path], n: int
    ) -> None:
        rows = []
        for f in files:
            for e in json.loads(f.read_text()):
                if not e.get("final_y"):
                    rows.append(
                        {
                            "label": f"{name} JD {e['epoch_jd']:.1f}",
                            "pass": False,
                            "error": "no chain stored (did not converge)",
                        }
                    )
                    print(f"{name} JD {e['epoch_jd']:.1f}: no chain stored -> FAIL", flush=True)
                    continue
                rows.append(
                    judge(
                        cell,
                        key,
                        xd,
                        n,
                        "blend",
                        e["x0_days"],
                        e["final_y"],
                        real,
                        f"{name} JD {e['epoch_jd']:.1f}",
                    )
                )
        out[name] = rows

    # positive control first
    d316 = REPO / "data" / "943_ganeur316_realeph"
    run_blend(
        "CONTROL+ GanEur#316",
        "ge",
        r316["key"],
        r316["x_days"],
        [d316 / "n10_std" / "realeph_chain.json"],
        10,
    )
    run_blend(
        "REPORT GanEur#316@2019",
        "ge",
        r316["key"],
        r316["x_days"],
        [d316 / "n10_rs2019_rel" / "realeph_chain.json"],
        10,
    )
    # GanCal#1@2013: reported only (gate indeterminate by construction)
    s_gc, a, b = CHAIN.ENUM.cell_system("gc")
    key1 = "k3|LGanymede>Ganymede/1l|RGanymede/2:1|LGanymede>Callisto/0s|LCallisto>Ganymede/0s"
    from cyclerfinder.search.two_working_body import correct_dates

    _, cyc1 = CHAIN.ENUM.parse_cycle_key(key1, s_gc, a, b)
    sol = correct_dates(s_gc, cyc1, np.array([0.991426, 26.049946, 35.954310]) * DAY, tol_kms=1e-8)
    run_blend(
        "REPORT GanCal#1@2013",
        "gc",
        key1,
        [float(v) / DAY for v in sol.x],
        [REPO / "data" / "943_c4_rs2013" / "n10_grow" / "realeph_chain.json"],
        10,
    )
    # negative control: gc-1 epoch 0, 2nd Lambert start + 0.05 d
    g1 = wb["gc-1"]
    e0 = json.loads((REPO / "data" / "943_gc1_realeph" / "e0" / "realeph_chain.json").read_text())[
        0
    ]
    y_bad = list(e0["final_y"])
    y_bad[0] += 0.05
    neg = judge(
        "gc",
        g1["key"],
        g1["x_days"],
        10,
        "blend",
        e0["x0_days"],
        y_bad,
        real,
        "CONTROL- gc-1 e0 broken",
    )
    out["CONTROL- gc-1 e0 broken"] = [neg]
    pos_ok = sum(r["pass"] for r in out["CONTROL+ GanEur#316"]) >= 3
    neg_ok = not neg["pass"]
    out["controls_ok"] = {"positive": pos_ok, "negative": neg_ok}
    print(
        f"CONTROLS: positive {'PASS' if pos_ok else 'FAIL'}, negative "
        f"{'FAILS as required' if neg_ok else 'PASSED (BAD)'}",
        flush=True,
    )
    if not (pos_ok and neg_ok):
        print("STOP: a control did not behave; rows not judged", flush=True)
        return out
    # rows
    run_blend(
        "gc-1",
        "gc",
        g1["key"],
        g1["x_days"],
        [REPO / "data" / "943_gc1_realeph" / f"e{i}" / "realeph_chain.json" for i in range(5)],
        10,
    )
    g2 = wb["gc-2"]
    ladder = json.loads((REPO / "data" / "942_943_ladder_realeph_chains.json").read_text())["gc2"]
    rows = []
    for e in ladder:
        # ladder layout (note 6.40): the 4 n Lambert-leg start dates, then the chain period
        xs = e["final_x_days"]
        y = list(xs[1:])
        n = (len(xs) - 1) // 4
        rows.append(
            judge(
                "gc",
                g2["key"],
                g2["x_days"],
                n,
                "blend",
                xs[0],
                y,
                real,
                f"gc-2 JD {e['epoch_jd']:.1f}",
            )
        )
    out["gc-2"] = rows
    for name in ("gc-1", "gc-2"):
        k = sum(r["pass"] for r in out[name])
        print(f"{name}: V2 {'PASS' if k >= 3 else 'FAIL'} ({k}/5 epochs)", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lane", required=True, choices=["gc"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    res = lane_gc()
    (OUT / f"v2_chain_{args.lane}.json").write_text(json.dumps(res, indent=1, default=float))
    print("wrote", (OUT / f"v2_chain_{args.lane}.json").relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
