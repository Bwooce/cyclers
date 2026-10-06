"""#942: closure Delta-V of full-rev legs on a non-Keplerian ephemeris (note sec. 6.24).

Input: a blend-mode chain written by scripts/run_942_realeph_chain.py (the lambda = 1 date solution,
fixed legs timed at the model period from the minimax directions). For each full-rev leg:
- D: the miss between the spacecraft (DOP853 two-body re-fly, real GM) and the planet at
  arrival, km;
- dv_mid: the exact mid-course correction at half the flight time that puts the arc on the planet at
  the same arrival time (two-body Newton on the 3 arrival-position equations), m/s;
- dvinf: the change of the arrival |V_inf| that this correction causes, m/s.
Descriptive only, no pass/fail.

Usage: uv run python scripts/analyse_942_fullrev_closure_dv.py --cell ev --key KEY --n-cycles N
        --real de440 --chain FILE [FILE ...] --out data/....json
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.search.two_working_body import (
    _blocks,
    cycle_flybys,
    eval_lambert_legs,
    fixed_duration_s,
    kepler_step,
)

DAY = SECONDS_PER_DAY
REPO = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CHECK = _load("check_942_realeph_chain", REPO / "scripts" / "check_942_realeph_chain.py")
CHAIN = CHECK.CHAIN


def mid_course(
    r_m: np.ndarray, v_m: np.ndarray, dt: float, r_target: np.ndarray, mu: float
) -> np.ndarray:
    """Exact impulse at (r_m, v_m) that reaches r_target after dt (two-body Newton)."""
    dv = np.zeros(3)
    for _ in range(30):
        r_end, _ = kepler_step(r_m, v_m + dv, dt, mu)
        f = r_end - r_target
        if float(np.linalg.norm(f)) < 1e-6:
            return dv
        jac = np.empty((3, 3))
        h = 1e-7
        for k in range(3):
            e = np.zeros(3)
            e[k] = h
            jac[:, k] = (kepler_step(r_m, v_m + dv + e, dt, mu)[0] - r_end) / h
        dv = dv - np.linalg.solve(jac, f)
    raise RuntimeError("mid-course Newton did not converge")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--n-cycles", type=int, required=True)
    ap.add_argument("--real", default="de440", choices=["de440", "spice", "mean"])
    ap.add_argument("--chain", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    circ, a, b = CHAIN.ENUM.cell_system(args.cell)
    _, one = CHAIN.ENUM.parse_cycle_key(args.key, circ, a, b)
    legs_chain = one.legs * args.n_cycles
    real = CHAIN.real_ephemeris(args.cell, args.real)
    mu = float(real.mu)
    out: list[dict[str, Any]] = []
    for path in args.chain:
        for ep in json.loads(path.read_text()):
            if ep.get("last_converged_lambda") != 1.0:
                print(f"epoch JD {ep['epoch_jd']:.1f}: blend did not reach lambda = 1, skipped")
                continue
            sysm = CHAIN.Blend(circ, real, np.eye(3), 0.0, 0.0, 1.0)  # lambda = 1: real only
            y = np.asarray(ep["final_y_dates"])
            x = np.concatenate([[ep["x0_days"]], y[:-1]]) * DAY
            cyc = CHAIN.Cycle(legs_chain, float(y[-1]) * DAY)
            legs = eval_lambert_legs(sysm, cyc, x)
            flybys = cycle_flybys(sysm, cyc, x)
            assert legs is not None and flybys is not None
            rows = []
            for blk in _blocks(sysm, cyc, legs):
                for leg, t0 in blk.fixed:
                    fb = next(f for f in flybys if abs(f.t_s - t0) < 1.0 and f.body == blk.body)
                    dt = fixed_duration_s(sysm, leg)
                    r0, w0 = real.state(blk.body, t0)
                    r1, v1 = CHECK.fly(mu, r0, w0 + fb.vinf_out, dt)
                    rb, wb = real.state(blk.body, t0 + dt)
                    # Correction half a spacecraft revolution before arrival (amendment 6.24a):
                    # a whole revolution left would make the position map singular.
                    t_rem = 0.5 * dt / getattr(leg, "sc_revs", 1)
                    r_m, v_m = CHECK.fly(mu, r0, w0 + fb.vinf_out, dt - t_rem)
                    dv = mid_course(r_m, v_m, t_rem, rb, mu)
                    r2, v2 = CHECK.fly(mu, r_m, v_m + dv, t_rem)  # independent re-fly
                    rows.append(
                        {
                            "body": blk.body,
                            "t0_days": t0 / DAY,
                            "miss_km": float(np.linalg.norm(r1 - rb)),
                            "dv_mid_ms": 1e3 * float(np.linalg.norm(dv)),
                            "refly_miss_after_km": float(np.linalg.norm(r2 - rb)),
                            "dvinf_ms": 1e3
                            * (float(np.linalg.norm(v2 - wb)) - float(np.linalg.norm(v1 - wb))),
                        }
                    )
            steps = [s for s in ep["steps"] if s.get("converged")]
            rec = {
                "epoch_jd": ep["epoch_jd"],
                "gate_at_lambda1_minimax": {
                    "status": steps[-1].get("status"),
                    "worst_ratio": steps[-1].get("worst_ratio"),
                }
                if steps
                else None,
                "shoot_best": ep.get("shoot_best"),
                "legs": rows,
                "sum_dv_mid_ms": sum(r["dv_mid_ms"] for r in rows),
                "max_dv_mid_ms": max((r["dv_mid_ms"] for r in rows), default=0.0),
                "max_abs_dvinf_ms": max((abs(r["dvinf_ms"]) for r in rows), default=0.0),
                "max_miss_km": max((r["miss_km"] for r in rows), default=0.0),
            }
            out.append(rec)
            print(
                f"epoch JD {ep['epoch_jd']:.1f}: {len(rows)} full-rev legs, miss <= "
                f"{rec['max_miss_km']:.0f} km, mid-course dv sum {rec['sum_dv_mid_ms']:.2f} m/s "
                f"(max {rec['max_dv_mid_ms']:.2f}), |V_inf| change <= "
                f"{rec['max_abs_dvinf_ms']:.2f} m/s; ballistic shoot: "
                f"{(ep.get('shoot_best') or {}).get('status', 'none')} "
                f"{(ep.get('shoot_best') or {}).get('worst_ratio', float('nan')):.3f}",
                flush=True,
            )
    args.out.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
