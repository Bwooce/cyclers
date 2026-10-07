"""#998 real-ephemeris check of the 'strong' Pluto small-moon cyclers (plu060.bsp).

The ideal-model cycle is placed on the real Charon and small-moon states (NAIF 901..905 relative
to the Pluto-system barycentre, 9) as a patched-conic ephemeris model: the spacecraft moves on
conics about the barycentre with the system GM 975.5 km^3/s^2 (the Russell-Strange 'patched-conic
ephemeris model'), the bodies move as in the kernel. The frame is the mean orbit plane of Charon
over 2030-2034 (z along its orbit normal), rotated so that Charon is at angle 0 at the epoch.

Epoch alignment: the ideal cycle starts with both bodies at angle 0, so the epoch T0 is the first
real conjunction (moon and Charon at the same in-plane angle) after a nominal date. The ideal start
dates are then used as the first guess, and the Lambert-leg dates are corrected on the real states
(two_working_body.correct_dates, least squares on the H&M magnitude/vector residual). Reported per
epoch: the residual before and after correction (km/s), the date shifts, and the demanded-turn
gate with the real V_inf vectors (registry floors).

Usage: uv run python scripts/realeph_998.py --cell pn --top 3 --n-cycles 1 --out F
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import spiceypy as sp
from scipy.optimize import brentq

from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    FlybyBody,
    correct_dates,
    cycle_flybys,
    date_residual,
    encounter_self_consistency,
    gate_cycle,
)
from cyclerfinder.verify.spice_kernels import ensure_pluto_kernel

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
NAIF = {"Charon": 901, "Styx": 905, "Nix": 902, "Kerberos": 904, "Hydra": 903}
LSK = REPO / "src" / "cyclerfinder" / "verify" / "kernels" / "naif0012.tls"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


class EphSystem:
    """plu060 states in the Charon mean-plane frame, Charon at angle 0 at ``et0``."""

    def __init__(self, ideal: CircularSystem, et0: float, frame: np.ndarray) -> None:
        self.mu = ideal.mu
        self.massless = ideal.massless
        self._ideal = ideal
        self.et0 = et0
        self.frame = frame
        st, _ = sp.spkezr("901", et0, "J2000", "NONE", "9")
        p = frame @ st[:3]
        self.rot = _rz(-math.atan2(p[1], p[0]))

    def state(self, code: str, t_s: float) -> tuple[np.ndarray, np.ndarray]:
        st, _ = sp.spkezr(str(NAIF[code]), self.et0 + t_s, "J2000", "NONE", "9")
        return self.rot @ self.frame @ st[:3], self.rot @ self.frame @ st[3:]

    def period_s(self, code: str) -> float:
        return self._ideal.period_s(code)

    def body(self, code: str) -> FlybyBody:
        return self._ideal.body(code)

    def wrap_rotation(self, code: str, dt_s: float) -> np.ndarray:
        return self._ideal.wrap_rotation(code, dt_s)

    def synodic_s(self, a: str, b: str) -> float:
        return self._ideal.synodic_s(a, b)


def mean_plane_frame(et_a: float, et_b: float) -> np.ndarray:
    ts = np.linspace(et_a, et_b, 600)
    pos, vel = [], []
    for t in ts:
        st, _ = sp.spkezr("901", t, "J2000", "NONE", "9")
        pos.append(st[:3])
        vel.append(st[3:])
    h = np.cross(np.array(pos), np.array(vel)).mean(axis=0)
    n = h / np.linalg.norm(h)
    e1 = np.array(pos[0]) - np.dot(pos[0], n) * n
    e1 /= np.linalg.norm(e1)
    return np.vstack([e1, np.cross(n, e1), n])


def phase_diff(frame: np.ndarray, moon: str, et: float) -> float:
    ang = []
    for nid in (NAIF[moon], 901):
        st, _ = sp.spkezr(str(nid), et, "J2000", "NONE", "9")
        p = frame @ st[:3]
        ang.append(math.atan2(p[1], p[0]))
    d = ang[0] - ang[1]
    return (d + math.pi) % (2.0 * math.pi) - math.pi


def next_conjunction(frame: np.ndarray, moon: str, et_start: float) -> float:
    """First time after ``et_start`` the moon and Charon share the in-plane angle (moon passes
    Charon from behind: the phase difference crosses zero going from negative to positive is not
    required; either sign accepted, the first zero crossing)."""
    step = 0.05 * DAY
    t = et_start
    f0 = phase_diff(frame, moon, t)
    for _ in range(int(80 * DAY / step)):
        f1 = phase_diff(frame, moon, t + step)
        if f0 * f1 < 0 and abs(f0) < 1.0 and abs(f1) < 1.0:
            return float(brentq(lambda s: phase_diff(frame, moon, s), t, t + step, xtol=1e-3))
        t += step
        f0 = f1
    raise RuntimeError("no conjunction found")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--n-cycles", type=int, default=1)
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    drv = _load("run_998_enumerate", REPO / "scripts" / "run_998_enumerate.py")
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    moon = drv.CELLS[args.cell]
    ideal = drv.pluto_system(moon)
    g = json.loads(
        (REPO / "data" / "998_pluto_smallmoons" / f"{args.cell}_gauntlet.json").read_text()
    )
    strong = sorted((c for c in g["candidates"] if c["strong"]), key=lambda c: c["worst_ratio"])
    picks = strong[: args.top]

    sp.furnsh(str(LSK))
    sp.furnsh(ensure_pluto_kernel())
    et_a = sp.str2et("2030-01-01T00:00:00")
    frame = mean_plane_frame(et_a, et_a + 4 * 365.25 * DAY)
    out = []
    for cand in picks:
        _, cyc1 = enum.parse_cycle_key(cand["key"], ideal, "Charon", moon)
        n = args.n_cycles
        cyc = Cycle(cyc1.legs * n, cyc1.period_s * n)
        x0 = np.tile(np.asarray(cand["x_days"]) * DAY, n)
        x0 = x0 + np.repeat(np.arange(n) * cyc1.period_s, len(cand["x_days"]))
        rows = []
        for ie in range(args.epochs):
            et_nom = et_a + ie * 400.0 * DAY
            et0 = next_conjunction(frame, moon, et_nom)
            eph = EphSystem(ideal, et0, frame)
            r_ideal = date_residual(ideal, cyc, x0)
            r_before = date_residual(eph, cyc, x0)
            rec: dict[str, Any] = {
                "epoch_utc": sp.et2utc(et0, "ISOC", 3),
                "residual_ideal_model_kms": None
                if r_ideal is None
                else float(np.max(np.abs(r_ideal))),
                "residual_before_kms": None
                if r_before is None
                else float(np.max(np.abs(r_before))),
            }
            sol = correct_dates(eph, cyc, x0, tol_kms=1e-6)
            rec["residual_after_kms"] = sol.max_abs_residual_kms
            rec["converged_1e-6"] = bool(sol.converged)
            rec["date_shift_hours_max"] = float(np.max(np.abs(sol.x - x0)) / 3600.0)
            fl = cycle_flybys(eph, cyc, sol.x)
            if fl is None:
                rec["gate_status"] = "no-directions"
            else:
                rep = gate_cycle(eph, fl)
                rec["gate_status"] = rep.status
                rec["worst_ratio"] = rep.gate.worst_ratio
                rec["max_turn_deg"] = rep.max_turn_deg
                rec["vinf_charon_kms"] = [
                    float(np.linalg.norm(f.vinf_in)) for f in fl if f.body == "Charon"
                ]
            rec["encounter_miss_km"] = encounter_self_consistency(eph, cyc, sol.x)
            rows.append(rec)
            print(cand["key"][:60], rec, flush=True)
        out.append(
            {
                "key": cand["key"],
                "k": cand["k"],
                "ideal": {
                    "worst_ratio": cand["worst_ratio"],
                    "vinf_charon_kms": cand["vinf_charon_kms"],
                    "vinf_moon_kms": cand["vinf_moon_kms"],
                    "period_days": cand["period_days"],
                },
                "n_cycles": n,
                "epochs": rows,
            }
        )
    args.out.write_text(json.dumps({"cell": args.cell, "moon": moon, "candidates": out}, indent=1))


if __name__ == "__main__":
    main()
