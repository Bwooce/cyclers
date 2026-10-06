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
from cyclerfinder.search.two_working_body import cycle_flybys

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
GAUNTLET = _load("gauntlet_942", REPO / "scripts" / "gauntlet_942.py")


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
    ap.add_argument(
        "--include-failed",
        action="store_true",
        help="also re-fly a stored closure whose rung did not pass (e.g. a gate fail)",
    )
    args = ap.parse_args()
    circ, a, b = CHAIN.ENUM.cell_system(args.cell)
    _, one = CHAIN.ENUM.parse_cycle_key(args.key, circ, a, b)
    legs_chain = one.legs * args.n_cycles
    real = CHAIN.real_ephemeris(args.cell, args.real)
    n_lam = sum(isinstance(lg, CHAIN.LambertLeg) for lg in legs_chain)
    for path in args.chain:
        for ep in json.loads(path.read_text()):
            if not ep["rung_pass"] and not (args.include_failed and "final_y" in ep):
                print(f"epoch JD {ep['epoch_jd']:.1f}: rung not passed, skipped")
                continue
            if ep.get("mode") == "ramp":
                # Ramp mode stores dates only; at lambda = 1 the ramped model IS the
                # Standish mean-element system (checked to 1e-5 km), so re-fly against it:
                # Lambert legs from the dates, full-rev legs from the minimax directions.
                y = np.asarray(ep["final_y_dates"])
                x = np.concatenate([[ep["x0_days"]], y[:-1]]) * DAY
                cyc = CHAIN.Cycle(legs_chain, float(y[-1]) * DAY)
                flybys = cycle_flybys(real, cyc, x)
                assert flybys is not None
                cc = GAUNTLET.cross_check(real, cyc, x, flybys)
                print(
                    f"epoch JD {ep['epoch_jd']:.1f}: ramp mode, DOP853 max arrival miss "
                    f"{cc['max_arrival_miss_km']:.3e} km, max V_inf vector error "
                    f"{cc['max_vinf_vector_error_kms']:.3e} km/s (all legs, full-revs included), "
                    f"max junction |V_inf| mismatch {cc['max_junction_mismatch_kms']:.3e} km/s "
                    f"(includes the chain wrap, closed in magnitude only)"
                )
                continue
            sysm = CHAIN.Blend(circ, real, np.eye(3), 0.0, 0.0, 1.0)  # lambda = 1: real only
            ev = CHAIN.chain_eval(sysm, legs_chain, ep["x0_days"] * DAY, np.asarray(ep["final_y"]))
            assert ev is not None
            worst_miss = worst_dv = 0.0
            for frm, t0, v0, to, t1 in ev.segments:
                r0, _ = real.state(frm, t0)
                r1, v1 = fly(real.mu, r0, v0, t1 - t0)  # the real GM, not the solver's
                rb, _ = real.state(to, t1)
                _, vk = CHAIN.kepler_step(r0, v0, t1 - t0, real.mu)
                worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
                worst_dv = max(worst_dv, float(np.linalg.norm(v1 - vk)))
            print(
                f"epoch JD {ep['epoch_jd']:.1f}: segments {len(ev.segments)} (Lambert {n_lam}), "
                f"DOP853 max arrival miss {worst_miss:.3e} km, "
                f"max velocity difference vs the solver's conic {worst_dv:.3e} km/s"
            )


if __name__ == "__main__":
    main()
