"""#942 robustness (d) for vm2-1: patched-conic real-ephemeris continuation (R-S 2007 style).

Pre-registered in docs/notes/2026-10-05-942-943-two-working-body-generator.md sec. 6.5 (d).

A chain of n consecutive vm2-1 cycles (3 Lambert legs each: V->V 2-rev low, V->M, M->V) has its
3n + 1 encounter dates d_0 .. d_3n; d_0 (the epoch) is fixed and d_1 .. d_3n are free. Residuals:
the V_inf magnitude difference at the 3n - 1 interior encounters, plus the final arrival magnitude
against the first departure magnitude (H&M closure "B"). That makes the system square.

The planets are blended from the R-S circular model, rotated and time-shifted to match the real
mean longitudes at the epoch, to fixed Standish & Williams J2000 mean elements (real periods, e, i)
by linear interpolation of position and velocity with weight lambda (R-S 2007 p.6 homotopy),
lambda = 0, 0.1, ..., 1. At each step the chain is re-solved from the previous one, and every
interior encounter is judged by the #888/#937 gate with R-S 2007 Table 2 constants and registry
floors.

Epochs: 5, spread over one 32-yr cycle from 2030-01-01 (every 6.4 yr). At each, the nearest date
whose Venus-Mars relative mean longitude equals vm2-1's starting configuration is used.

Usage: uv run python scripts/run_942_vm2_realeph.py --out DIR [--n-cycles 7] [--epochs 5]
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Flyby,
    FlybyBody,
    MeanElementSystem,
    Vec,
    gate_cycle,
)

DAY = SECONDS_PER_DAY
#: vm2-1 (results note sec. 6.4): Lambert-leg starts in the R-S circular model (d).
VM2_1_STARTS_D = (52.916704638791344, 614.929974502681, 834.8083489268107)
VM2_1_PERIOD_D = 1001.7700187121935
LEGS = (("V", "V", 2, "low"), ("V", "M", 0, "single"), ("M", "V", 0, "single"))


def rs_circular() -> CircularSystem:
    mu = 1.3271244e11
    bodies = {}
    for c, per in (("V", 19_414_153.0), ("M", 59_354_429.0)):
        bodies[c] = ((mu * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0), per, 0.0)
    over = {
        "V": FlybyBody("V", 324_860.0, 6052.0, 300.0),
        "M": FlybyBody("M", 42_828.3, 3399.0, 200.0),
    }
    return CircularSystem(mu, bodies, frozenset(), flyby_overrides=over)


def mean_longitude(eph: MeanElementSystem, code: str, t_s: float) -> float:
    r, _ = eph.state(code, t_s)
    return math.atan2(float(r[1]), float(r[0]))


@dataclass
class Blend:
    """lambda-blend of a rotated, time-shifted circular model and the mean-element ephemeris.

    Times are seconds past JD 2440000.0 (MeanElementSystem's convention).
    """

    circ: CircularSystem
    eph: MeanElementSystem
    t_shift: float  # circular-model time = t - t_shift
    rot: float  # rotation of the circular model about z
    lam: float = 0.0
    mu: float = field(init=False)

    def __post_init__(self) -> None:
        self.mu = self.circ.mu

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        rc, vc = self.circ.state(code, t_s - self.t_shift)
        c, s = math.cos(self.rot), math.sin(self.rot)
        m = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        rc, vc = m @ rc, m @ vc
        if self.lam == 0.0:
            return rc, vc
        re, ve = self.eph.state(code, t_s)
        return (1 - self.lam) * rc + self.lam * re, (1 - self.lam) * vc + self.lam * ve

    def body(self, code: str) -> FlybyBody:
        return self.circ.body(code)


def chain_legs(n: int) -> list[tuple[str, str, int, str]]:
    return [LEGS[i % 3] for i in range(3 * n)]


def leg_vinf(
    sysm: Blend, leg: tuple[str, str, int, str], t0: float, t1: float
) -> tuple[Vec, Vec] | None:
    frm, to, nrev, br = leg
    if t1 <= t0:
        return None
    r1, w1 = sysm.state(frm, t0)
    r2, w2 = sysm.state(to, t1)
    try:
        sols = lambert(r1, r2, t1 - t0, mu=sysm.mu, max_revs=nrev)
    except (LambertError, ValueError):
        return None
    sol = next((s for s in sols if s.n_revs == nrev and s.branch == br), None)
    if sol is None:
        return None
    return sol.v1 - w1, sol.v2 - w2


def residual(sysm: Blend, legs: list, d0: float, y: np.ndarray) -> np.ndarray:
    dates = np.concatenate([[d0], y]) * DAY
    vs = []
    for i, leg in enumerate(legs):
        v = leg_vinf(sysm, leg, dates[i], dates[i + 1])
        if v is None:
            return np.full(len(legs), 1e3)
        vs.append(v)
    res = [
        float(np.linalg.norm(vs[i][1]) - np.linalg.norm(vs[i + 1][0])) for i in range(len(legs) - 1)
    ]
    res.append(float(np.linalg.norm(vs[-1][1]) - np.linalg.norm(vs[0][0])))
    return np.asarray(res)


def gate_chain(sysm: Blend, legs: list, dates_d: np.ndarray) -> dict:
    dates = dates_d * DAY
    vs = [leg_vinf(sysm, leg, dates[i], dates[i + 1]) for i, leg in enumerate(legs)]
    if any(v is None for v in vs):
        return {"status": "lambert-fail"}
    fl = [Flyby(legs[i][1], dates[i + 1], vs[i][1], vs[i + 1][0]) for i in range(len(legs) - 1)]
    rep = gate_cycle(sysm, fl)  # type: ignore[arg-type]
    ratios = {"V": [], "M": []}
    for f, e in zip(fl, rep.gate.encounters, strict=True):
        ratios[f.body].append(e.ratio)
    return {
        "status": rep.status,
        "worst_ratio": rep.gate.worst_ratio,
        "max_venus_ratio": max(ratios["V"]),
        "max_mars_ratio": max(ratios["M"]),
        "min_required_alt_km": rep.gate.min_required_alt_km,
        "vinf_range_kms": [
            min(float(np.linalg.norm(v[0])) for v in vs),
            max(float(np.linalg.norm(v[0])) for v in vs),
        ],
    }


def epoch_for(
    eph: MeanElementSystem, circ: CircularSystem, near_jd: float
) -> tuple[float, float, float]:
    """(epoch t_s past JD 2440000, t_shift, rot): the real date near ``near_jd`` whose
    Venus-Mars relative mean longitude equals vm2-1's start configuration."""
    t_start_c = VM2_1_STARTS_D[0] * DAY
    th = {c: 2 * math.pi * t_start_c / circ.period_s(c) for c in ("V", "M")}
    target = (th["M"] - th["V"]) % (2 * math.pi)

    def g(t_s: float) -> float:
        d = (mean_longitude(eph, "M", t_s) - mean_longitude(eph, "V", t_s) - target) % (2 * math.pi)
        return (d + math.pi) % (2 * math.pi) - math.pi

    t0 = (near_jd - 2440000.0) * DAY
    syn = circ.synodic_s("V", "M")
    grid = np.linspace(t0, t0 + syn, 200)
    vals = [g(float(t)) for t in grid]
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0 and abs(vals[i] - vals[i + 1]) < math.pi:
            te = brentq(g, float(grid[i]), float(grid[i + 1]))
            rot = mean_longitude(eph, "V", te) - th["V"]
            return te, te - t_start_c, rot
    raise RuntimeError("no phase match")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-cycles", type=int, default=7)
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--steps", type=int, default=10)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    preflight_search(
        task_no=942,
        region_id=f"vm2-1-realeph-chain-n{args.n_cycles}",
        method=MethodCapability(
            genome="vm2-1 chain of n cycles, open dates, H&M closure B",
            corrector="least squares on dates, circular-to-mean-element homotopy",
            capability_tags=frozenset({"ballistic", "patched-conic", "3d", "inclined-elliptic"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=args.epochs,
    )
    circ = rs_circular()
    eph = MeanElementSystem()
    legs = chain_legs(args.n_cycles)
    base = [VM2_1_STARTS_D[i % 3] + (i // 3) * VM2_1_PERIOD_D for i in range(3 * args.n_cycles)]
    base.append(VM2_1_STARTS_D[0] + args.n_cycles * VM2_1_PERIOD_D)
    out = []
    t_run = time.time()
    for ie in range(args.epochs):
        near = 2462502.5 + ie * 32.0 * 365.25 / args.epochs  # 2030-01-01 + i * 6.4 yr
        te, tshift, rot = epoch_for(eph, circ, near)
        sysm = Blend(circ, eph, tshift, rot, 0.0)
        y = np.array(base[1:]) + tshift / DAY
        d0 = base[0] + tshift / DAY
        steps = []
        for k in range(args.steps + 1):
            sysm.lam = k / args.steps
            sol = least_squares(
                lambda yy, s_=sysm, d_=d0: residual(s_, legs, d_, yy),
                y,
                method="lm",
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
                max_nfev=200 * len(y),
            )
            res = residual(sysm, legs, d0, sol.x)
            ok = bool(np.max(np.abs(res)) < 1e-6)
            g = gate_chain(sysm, legs, np.concatenate([[d0], sol.x]))
            rec = {
                "lambda": sysm.lam,
                "max_residual_kms": float(np.max(np.abs(res))),
                "converged": ok,
                "max_date_shift_d": float(np.max(np.abs(sol.x - y))),
            } | g
            steps.append(rec)
            print(
                f"{time.strftime('%H:%M:%S')} epoch {ie} (JD {te / DAY + 2440000:.1f}) "
                f"lam={sysm.lam:.1f} conv={ok} res={rec['max_residual_kms']:.1e} "
                f"gate={g.get('status')} V={g.get('max_venus_ratio', float('nan')):.3f} "
                f"M={g.get('max_mars_ratio', float('nan')):.3f} "
                f"shift={rec['max_date_shift_d']:.2f}d "
                f"[elapsed {time.time() - t_run:.0f}s]",
                flush=True,
            )
            if not ok:
                break
            y = sol.x
        out.append(
            {"epoch_jd": te / DAY + 2440000.0, "steps": steps, "final_dates_d": [d0, *y.tolist()]}
        )
        (args.out / "realeph.json").write_text(json.dumps(out, indent=1, default=float))
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
