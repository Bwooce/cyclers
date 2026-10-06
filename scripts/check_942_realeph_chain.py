"""#942/#943: independent re-fly of a converged real-ephemeris chain (rung d).

For each epoch's final (lambda = 1) chain written by scripts/run_942_realeph_chain.py, every
Lambert leg is re-flown with scipy DOP853 (rtol 1e-13, atol 1e-8 km) from the real-ephemeris
body state plus the solved departure V_inf, and compared with the arrival body's real position
(miss, km) and with the solved arrival V_inf (km/s). The interior-flyby V_inf magnitude
mismatch is rebuilt from the integrated arrivals. Fixed (full-rev) legs are re-flown from the
minimax-chosen direction. Not the Lambert solver, not the Kepler step.

Usage: uv run python scripts/check_942_realeph_chain.py --cell gc --key KEY --chain FILE
        --n-cycles N
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

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.search.two_working_body import (
    Cycle,
    _blocks,
    eval_lambert_legs,
    fixed_duration_s,
    optimise_block,
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


CHAIN = _load("run_942_realeph_chain", REPO / "scripts" / "run_942_realeph_chain.py")


def fly(mu: float, r0: np.ndarray, v0: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        return np.concatenate([y[3:], -mu * y[:3] / np.linalg.norm(y[:3]) ** 3])

    sol = solve_ivp(
        rhs, (0.0, dt), np.concatenate([r0, v0]), method="DOP853", rtol=1e-13, atol=1e-8
    )
    return sol.y[:3, -1], sol.y[3:, -1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--chain", type=Path, nargs="+", required=True)
    ap.add_argument("--n-cycles", type=int, required=True)
    ap.add_argument("--real", default="auto", choices=["auto", "mean", "de440", "spice"])
    args = ap.parse_args()
    circ, a, b = CHAIN.ENUM.cell_system(args.cell)
    _, one = CHAIN.ENUM.parse_cycle_key(args.key, circ, a, b)
    legs_chain = one.legs * args.n_cycles
    real = CHAIN.real_ephemeris(args.cell, args.real)
    for path in args.chain:
        for ep in json.loads(path.read_text()):
            if not ep["rung_pass"]:
                print(f"epoch JD {ep['epoch_jd']:.1f}: rung not passed, skipped")
                continue
            xs = ep["final_x_days"]
            x0, y = xs[0], np.asarray(xs[1:])
            cyc = Cycle(legs_chain, float(y[-1]) * DAY)
            x = np.concatenate([[x0], y[:-1]]) * DAY
            sysm = CHAIN.Blend(circ, real, np.eye(3), 0.0, 0.0, 1.0)  # lambda = 1: real only
            legs = eval_lambert_legs(sysm, cyc, x)
            assert legs is not None
            worst_miss = worst_dv = 0.0
            for ev, li in zip(legs, cyc.lambert_index, strict=True):
                leg = cyc.legs[li]
                r0, w0 = real.state(leg.frm, ev.t_dep)
                r1, v1 = fly(sysm.mu, r0, w0 + ev.vinf_dep, ev.t_arr - ev.t_dep)
                rb, wb = real.state(leg.to, ev.t_arr)
                worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
                worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - ev.vinf_arr)))
            for blk in _blocks(sysm, cyc, legs)[:-1]:
                res = optimise_block(sysm, blk)
                assert res is not None
                fl = res[0]
                for i, (leg, t0) in enumerate(blk.fixed):
                    dt = fixed_duration_s(sysm, leg)
                    r0, w0 = real.state(blk.body, t0)
                    r1, v1 = fly(sysm.mu, r0, w0 + fl[i].vinf_out, dt)
                    rb, wb = real.state(blk.body, t0 + dt)
                    worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
                    worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - fl[i + 1].vinf_in)))
            print(
                f"epoch JD {ep['epoch_jd']:.1f}: legs {len(legs)}, DOP853 max arrival miss "
                f"{worst_miss:.3e} km, max V_inf vector error {worst_dv:.3e} km/s"
            )


if __name__ == "__main__":
    main()
