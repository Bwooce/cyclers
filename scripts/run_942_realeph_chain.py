"""#942/#943 ladder rung (d), generic: real-ephemeris chain continuation of one ideal-model cycler.

Generalises scripts/run_942_vm2_realeph.py (vm2-1 only, Lambert legs only) to any cell and any
cycle key, including full-revolution and half-revolution legs, and to the Jovian moons.

A chain of n consecutive cycles is one Cycle of n x (legs). The first Lambert-leg start (the
epoch) is fixed; every other start date and the chain period (hence the final arrival) are free,
as in scripts/run_942_vm2_realeph.py. Fixed legs (full-rev, half-rev) are SHOT (2026-10-06 fix):
direction (2 angles) and flight time are unknowns and the arrival must be at the body (3
equations), so they close on the real ephemeris (the earlier version timed them at the ideal
period and they did not close). The wrap junction is the H&M closure "B"
magnitude match; it is a closure condition, NOT a flyby, so the wrap block's flybys are not
gated. The planets/moons are a lambda-blend (position and velocity) of the ideal circular model,
rotated into the bodies' mean orbital plane and phase-matched at the epoch, with the real
ephemeris:
  heliocentric cells: Standish & Williams J2000 fixed mean elements (MeanElementSystem);
  Jovian cells: NAIF jup365 via Ephemeris('spice', center='Jupiter') (J2000 equatorial).
lambda advances adaptively (start 0.1, halved on failure down to 1/640, as in the vm2-1
amendment). At every converged lambda the interior flybys are gated (#888/#937, registry
floors, the cell's flyby constants), with free directions by minimax.

Pass at this rung (pre-registration, results note sec. 6.13): lambda = 1 reached AND every
interior flyby gate-passes, at >= 1 epoch.

Usage: uv run python scripts/run_942_realeph_chain.py --cell gc --key KEY --x-days a,b,c
        --n-cycles N --epochs 5 --only-epoch I --out DIR
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    Flyby,
    FlybyBody,
    LambertLeg,
    MeanElementSystem,
    Vec,
    _blocks,
    date_residual,
    eval_lambert_legs,
    fixed_duration_s,
    gate_cycle,
    kepler_step,
    optimise_block,
)

DAY = SECONDS_PER_DAY
REPO = Path(__file__).resolve().parents[1]
JD_2440000_S_FROM_J2000 = (2440000.0 - 2451545.0) * DAY


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")


class DE440Planets:
    """Heliocentric DE440 states (astropy backend, J2000 ecliptic); ``t_s`` past JD 2440000.0."""

    def __init__(self) -> None:
        self.eph = Ephemeris("astropy")

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        r, v = self.eph.state(code, t_s + JD_2440000_S_FROM_J2000)
        return np.asarray(r, dtype=float), np.asarray(v, dtype=float)


def real_ephemeris(cell: str, which: str) -> Any:
    """``which``: "auto" (SPICE for Jovian cells, Standish mean elements otherwise),
    "mean", "de440" or "spice"."""
    if which == "auto":
        which = "spice" if cell in ENUM.X1_CELLS else "mean"
    return {"spice": SpiceMoons, "de440": DE440Planets, "mean": MeanElementSystem}[which]()


class SpiceMoons:
    """Jupiter-centred real moon states; ``t_s`` is seconds past JD 2440000.0."""

    def __init__(self) -> None:
        self.eph = Ephemeris("spice", center="Jupiter")

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        r, v = self.eph.state(code, t_s + JD_2440000_S_FROM_J2000)
        return np.asarray(r, dtype=float), np.asarray(v, dtype=float)


def plane_basis(real: Any, codes: tuple[str, str], t_s: float) -> np.ndarray:
    """Rotation whose columns are (x, y, z) of the mean orbital plane of the two bodies
    (z = the average orbit normal at ``t_s``; x = the J2000 x-axis projected into the plane)."""
    n = np.zeros(3)
    for c in codes:
        r, v = real.state(c, t_s)
        h = np.cross(r, v)
        n += h / np.linalg.norm(h)
    z = n / np.linalg.norm(n)
    x = np.array([1.0, 0.0, 0.0]) - z[0] * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.column_stack([x, y, z])


def in_plane_longitude(real: Any, basis: np.ndarray, code: str, t_s: float) -> float:
    r, _ = real.state(code, t_s)
    p = basis.T @ r
    return math.atan2(float(p[1]), float(p[0]))


@dataclass
class Blend:
    circ: CircularSystem
    real: Any
    basis: np.ndarray
    t_shift: float
    rot: float
    lam: float = 0.0
    mu: float = field(init=False)

    def __post_init__(self) -> None:
        self.mu = self.circ.mu

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        rc, vc = self.circ.state(code, t_s - self.t_shift)
        c, s = math.cos(self.rot), math.sin(self.rot)
        rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        m = self.basis @ rz
        rc, vc = m @ rc, m @ vc
        if self.lam == 0.0:
            return rc, vc
        re, ve = self.real.state(code, t_s)
        return (1 - self.lam) * rc + self.lam * re, (1 - self.lam) * vc + self.lam * ve

    def period_s(self, code: str) -> float:
        return self.circ.period_s(code)

    def body(self, code: str) -> FlybyBody:
        return self.circ.body(code)

    def wrap_rotation(self, code: str, dt_s: float) -> Vec:
        return np.eye(3)


@dataclass
class RampedKepler:
    """Hollister 1969 p.368 homotopy: "increase the eccentricity and inclination to their actual
    values". Keplerian at EVERY lambda, so a full-rev return timed at the body's current period
    with |v_sc| = |V_P| is exact and the free-direction minimax stays valid.

    lambda = 0 is the ideal circular body (rotated, phase-matched at the epoch t_e). lambda = 1 is
    the fixed Standish & Williams J2000 mean-element orbit (MeanElementSystem) exactly. In
    between, a, mu and the epoch mean longitude are interpolated linearly, e and i are ramped
    from 0, and Omega and varpi are the Standish values. The period is Kepler's at (a, mu).
    """

    circ: CircularSystem
    real: MeanElementSystem
    t_e: float
    t_shift: float
    rot: float
    lam: float = 0.0
    mu: float = field(init=False)

    def __post_init__(self) -> None:
        self.mu = self.circ.mu

    def _params(self, code: str) -> tuple[float, float, float, float, float, float, float]:
        lam = self.lam
        a_c = self.circ.bodies[code][0]
        a_s, e_s, i_s, lan, varpi, l0 = self.real._elements(code)
        mu = (1 - lam) * self.circ.mu + lam * self.real.mu
        a = (1 - lam) * a_c + lam * a_s
        n = math.sqrt(mu / a**3)
        th_c = 2 * math.pi * (self.t_e - self.t_shift) / self.circ.period_s(code) + self.rot
        n_s = 2 * math.pi / self.real.period_s(code)
        l_s = l0 + n_s * (self.t_e + JD_2440000_S_FROM_J2000)
        l_s = th_c + ((l_s - th_c + math.pi) % (2 * math.pi) - math.pi)
        l_e = (1 - lam) * th_c + lam * l_s
        return a, lam * e_s, lam * i_s, lan, varpi, l_e, n

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        a, e, inc, lan, varpi, l_e, n = self._params(code)
        m_anom = (l_e - varpi + n * (t_s - self.t_e)) % (2 * math.pi)
        ecc_anom = m_anom
        for _ in range(60):
            d = (ecc_anom - e * math.sin(ecc_anom) - m_anom) / (1 - e * math.cos(ecc_anom))
            ecc_anom -= d
            if abs(d) < 1e-15:
                break
        ce, se = math.cos(ecc_anom), math.sin(ecc_anom)
        b = math.sqrt(1 - e * e)
        x, y = a * (ce - e), a * b * se
        rdot = n * a / (1 - e * ce)
        vx, vy = -rdot * se, rdot * b * ce
        cz, sz = math.cos(lan), math.sin(lan)
        ci, si = math.cos(inc), math.sin(inc)
        w = varpi - lan
        cw, sw = math.cos(w), math.sin(w)
        rz1 = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1.0]])
        rx = np.array([[1.0, 0, 0], [0, ci, -si], [0, si, ci]])
        rz2 = np.array([[cw, -sw, 0], [sw, cw, 0], [0, 0, 1.0]])
        m = rz1 @ rx @ rz2
        return m @ np.array([x, y, 0.0]), m @ np.array([vx, vy, 0.0])

    def period_s(self, code: str) -> float:
        return 2 * math.pi / self._params(code)[6]

    def body(self, code: str) -> FlybyBody:
        return self.circ.body(code)

    def wrap_rotation(self, code: str, dt_s: float) -> Vec:
        return np.eye(3)


def date_chain_residual(sysm: Any, legs: tuple, x0: float, y: np.ndarray) -> np.ndarray:
    """Ramp mode: dates [x0, y[:-1]] (days), chain period y[-1] (days); fixed legs timed at
    the model's current body period (exact in a Keplerian model)."""
    cyc = Cycle(legs, float(y[-1]) * DAY)
    r = date_residual(sysm, cyc, np.concatenate([[x0], np.asarray(y[:-1]) * DAY]))
    n = sum(isinstance(lg, LambertLeg) for lg in legs)
    return np.full(n, 1e3) if r is None else r


def interior_gate(sysm: Any, cycle: Cycle, x: np.ndarray) -> dict:
    """Ramp mode gate: minimax free directions; the last (closure) block is not a flyby."""
    legs = eval_lambert_legs(sysm, cycle, x)
    if legs is None:
        return {"status": "lambert-fail"}
    fl = []
    for blk in _blocks(sysm, cycle, legs)[:-1]:
        res = optimise_block(sysm, blk)
        if res is None:
            return {"status": "no-directions"}
        fl.extend(res[0])
    rep = gate_cycle(sysm, fl)
    by_body: dict[str, float] = {}
    for f, e in zip(fl, rep.gate.encounters, strict=True):
        by_body[f.body] = max(by_body.get(f.body, 0.0), e.ratio)
    return {
        "status": rep.status,
        "worst_ratio": rep.gate.worst_ratio,
        "worst_ratio_by_body": by_body,
        "min_required_alt_km": rep.gate.min_required_alt_km,
        "n_flybys": len(fl),
    }


def _fixed_basis(w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    w_hat = w / np.linalg.norm(w)
    ref = np.array([0.0, 0.0, 1.0])
    e1 = ref - (ref @ w_hat) * w_hat
    e1 /= np.linalg.norm(e1)
    return w_hat, e1, np.cross(w_hat, e1)


def _dir(theta: float, phi: float, w: np.ndarray) -> np.ndarray:
    w_hat, e1, e2 = _fixed_basis(w)
    return (
        math.sin(theta) * math.cos(phi) * e1
        + math.sin(theta) * math.sin(phi) * e2
        + math.cos(theta) * w_hat
    )


def _angles(u: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    w_hat, e1, e2 = _fixed_basis(w)
    uh = u / np.linalg.norm(u)
    return math.acos(max(-1.0, min(1.0, float(uh @ w_hat)))), math.atan2(
        float(uh @ e2), float(uh @ e1)
    )


#: Position residual scale: 1 km of arrival miss counts as 1e-3 km/s of V_inf mismatch.
POS_SCALE_KM = 1000.0


@dataclass
class ChainEval:
    residual: np.ndarray
    block_flybys: list[list[Flyby]]
    segments: list[tuple[str, float, np.ndarray, str, float]]  # (from, t0, v0_sc, to, t1)


def _res(sysm: Blend, legs: tuple, x0: float, y: np.ndarray) -> np.ndarray:
    ev = chain_eval(sysm, legs, x0, y)
    n_fix = sum(not isinstance(lg, LambertLeg) for lg in legs)
    n_lam = len(legs) - n_fix
    return np.full(n_lam + 3 * n_fix, 1e3) if ev is None else ev.residual


def chain_eval(sysm: Blend, legs: tuple, x0: float, y: np.ndarray) -> ChainEval | None:
    """Evaluate the open chain with SHOT fixed legs (full-rev / half-rev).

    ``y`` (days, rad): starts of Lambert legs 2..N (days), the chain period (days), then
    (theta, phi, tau_days) for each fixed leg in chain order. Each fixed leg leaves its body
    with the junction's |V_inf| in direction (theta, phi) about the body velocity and is
    propagated for tau; its arrival must be at the body (3 residuals). Each block ends with the
    magnitude match against the next Lambert leg (1 residual). Square system.
    """
    lam_idx = [i for i, lg in enumerate(legs) if isinstance(lg, LambertLeg)]
    fix_idx = [i for i, lg in enumerate(legs) if not isinstance(lg, LambertLeg)]
    nl = len(lam_idx)
    dates = [x0, *(np.asarray(y[: nl - 1]) * DAY).tolist()]
    period = float(y[nl - 1]) * DAY
    fp = {i: y[nl + 3 * q : nl + 3 * q + 3] for q, i in enumerate(fix_idx)}
    res: list[float] = []
    block_flybys: list[list[Flyby]] = []
    segments: list[tuple[str, float, np.ndarray, str, float]] = []
    lam_v: list[tuple[np.ndarray, np.ndarray, float, float]] = []
    for j in range(nl):
        i0, i1 = lam_idx[j], lam_idx[(j + 1) % nl]
        between = []
        q = (i0 + 1) % len(legs)
        while q != i1:
            between.append(q)
            q = (q + 1) % len(legs)
        t_next = dates[j + 1] if j < nl - 1 else dates[0] + period
        t_arr = t_next - sum(float(fp[q][2]) * DAY for q in between)
        leg = legs[i0]
        if t_arr <= dates[j]:
            return None
        r1, w1 = sysm.state(leg.frm, dates[j])
        r2, w2 = sysm.state(leg.to, t_arr)
        try:
            sols = lambert(r1, r2, t_arr - dates[j], mu=sysm.mu, max_revs=leg.nrev)
        except (LambertError, ValueError):
            return None
        sol = next((s_ for s_ in sols if s_.n_revs == leg.nrev and s_.branch == leg.branch), None)
        if sol is None:
            return None
        lam_v.append((sol.v1 - w1, sol.v2 - w2, dates[j], t_arr))
        segments.append((leg.frm, dates[j], sol.v1, leg.to, t_arr))
    for j in range(nl):
        i0, i1 = lam_idx[j], lam_idx[(j + 1) % nl]
        body = legs[i0].to
        prev = lam_v[j][1]
        t = lam_v[j][3]
        fl: list[Flyby] = []
        q = (i0 + 1) % len(legs)
        while q != i1:
            theta, phi, tau = float(fp[q][0]), float(fp[q][1]), float(fp[q][2]) * DAY
            if tau <= 0.0:
                return None
            r, w = sysm.state(body, t)
            u = float(np.linalg.norm(prev)) * _dir(theta, phi, w)
            fl.append(Flyby(body, t, prev, u))
            r_end, v_end = kepler_step(r, w + u, tau, sysm.mu)
            rb, wb = sysm.state(body, t + tau)
            res.extend(((r_end - rb) / POS_SCALE_KM).tolist())
            segments.append((body, t, w + u, body, t + tau))
            prev = v_end - wb
            t += tau
            q = (q + 1) % len(legs)
        v_out = lam_v[(j + 1) % nl][0]
        fl.append(Flyby(body, t, prev, v_out))
        res.append(float(np.linalg.norm(prev) - np.linalg.norm(v_out)))
        block_flybys.append(fl)
    return ChainEval(np.asarray(res), block_flybys, segments)


def initial_fixed_params(
    sysm: Blend, cycle_ideal: Cycle, x_ideal: np.ndarray
) -> list[tuple[float, float, float]]:
    """(theta, phi, tau_days) per fixed leg IN ASCENDING LEG INDEX (chain_eval's order), from
    the ideal (lambda = 0) minimax directions. Walks the legs exactly as chain_eval does."""
    legs_all = cycle_ideal.legs
    lam_idx = cycle_ideal.lambert_index
    legs = eval_lambert_legs(sysm, cycle_ideal, x_ideal)  # type: ignore[arg-type]
    assert legs is not None
    blocks = _blocks(sysm, cycle_ideal, legs)  # type: ignore[arg-type]
    by_index: dict[int, tuple[float, float, float]] = {}
    for j, blk in enumerate(blocks):
        if not blk.fixed:
            continue
        res = optimise_block(sysm, blk)  # type: ignore[arg-type]
        assert res is not None
        fl = res[0]
        q = (lam_idx[j] + 1) % len(legs_all)
        for i, (leg, t0) in enumerate(blk.fixed):
            assert legs_all[q] is leg or legs_all[q] == leg
            _, w = sysm.state(blk.body, t0)
            th, ph = _angles(fl[i].vinf_out, w)
            by_index[q] = (th, ph, fixed_duration_s(sysm, leg) / DAY)  # type: ignore[arg-type]
            q = (q + 1) % len(legs_all)
    return [by_index[i] for i in sorted(by_index)]


def gate_eval(sysm: Blend, ev: ChainEval) -> dict:
    fl = [f for blk in ev.block_flybys[:-1] for f in blk]  # last block is the closure
    rep = gate_cycle(sysm, fl)  # type: ignore[arg-type]
    by_body: dict[str, float] = {}
    for f, e in zip(fl, rep.gate.encounters, strict=True):
        by_body[f.body] = max(by_body.get(f.body, 0.0), e.ratio)
    return {
        "status": rep.status,
        "worst_ratio": rep.gate.worst_ratio,
        "worst_ratio_by_body": by_body,
        "min_required_alt_km": rep.gate.min_required_alt_km,
        "n_flybys": len(fl),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--x-days", required=True)
    ap.add_argument("--n-cycles", type=int, required=True)
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--only-epoch", type=int, default=None)
    ap.add_argument("--first-epoch-jd", type=float, default=2462502.5)  # 2030-01-01
    ap.add_argument("--epoch-span-yr", type=float, default=32.0)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--real", default="auto", choices=["auto", "mean", "de440", "spice"])
    ap.add_argument("--shoot-restarts", type=int, default=20)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    preflight_search(
        task_no=942,
        region_id=f"realeph-chain-{args.cell}-{args.key}-n{args.n_cycles}",
        method=MethodCapability(
            genome="n-cycle chain of one ideal-model cycler, dates free, closure B",
            corrector="correct_dates; circular-to-real-ephemeris homotopy, adaptive step",
            capability_tags=frozenset({"ballistic", "patched-conic", "3d", "real-ephemeris"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=args.epochs,
    )
    circ, a, b = ENUM.cell_system(args.cell)
    if circ.massless:
        # A massless target in 3-D adds 2 equations per passage that dates alone cannot
        # meet (plane and magnitude); R-S's ephemeris model has every body massive (p.2
        # restricts masslessness to the ideal model). Use the both-massive cell.
        raise SystemExit(f"cell {args.cell!r} has a massless body; use its both-massive cell")
    _, one = ENUM.parse_cycle_key(args.key, circ, a, b)
    t_cyc = one.period_s
    legs_chain = one.legs * args.n_cycles
    x1 = np.array([float(v) for v in args.x_days.split(",")]) * DAY
    xs = np.concatenate([x1 + i * t_cyc for i in range(args.n_cycles)])
    real = real_ephemeris(args.cell, args.real)
    out = []
    t_run = time.time()
    for ie in range(args.epochs):
        if args.only_epoch is not None and ie != args.only_epoch:
            continue
        near_s = (
            args.first_epoch_jd + ie * args.epoch_span_yr * 365.25 / args.epochs - 2440000.0
        ) * DAY
        basis = plane_basis(real, (a, b), near_s)
        th = {c: 2 * math.pi * x1[0] / circ.period_s(c) for c in (a, b)}
        target = (th[b] - th[a]) % (2 * math.pi)

        def g(t_s: float, basis: np.ndarray = basis, target: float = target) -> float:
            d = in_plane_longitude(real, basis, b, t_s) - in_plane_longitude(real, basis, a, t_s)
            return ((d - target) % (2 * math.pi) + math.pi) % (2 * math.pi) - math.pi

        syn = circ.synodic_s(a, b)
        grid = np.linspace(near_s, near_s + 1.1 * syn, 440)  # > one period: a crossing exists
        vals = [g(float(t)) for t in grid]
        te = None
        for i in range(len(grid) - 1):
            if vals[i] * vals[i + 1] < 0 and abs(vals[i] - vals[i + 1]) < math.pi:
                lo, hi = float(grid[i]), float(grid[i + 1])
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if g(lo) * g(mid) <= 0:
                        hi = mid
                    else:
                        lo = mid
                te = 0.5 * (lo + hi)
                break
        assert te is not None, "no phase match"
        rot = in_plane_longitude(real, basis, a, te) - th[a]
        t_shift = te - x1[0]
        x0 = xs[0] + t_shift
        ramp = isinstance(real, MeanElementSystem)
        sysm: Any = (
            RampedKepler(circ, real, te, t_shift, rot, 0.0)
            if ramp
            else Blend(circ, real, basis, t_shift, rot, 0.0)
        )
        # Phase 1 (both modes): dates + chain period, fixed legs timed at the model's body
        # period with minimax directions. Exact at every lambda in ramp mode; a continuation
        # device only in blend mode.
        y = np.concatenate([(xs[1:] + t_shift) / DAY, [args.n_cycles * t_cyc / DAY]])
        r0 = date_chain_residual(sysm, legs_chain, x0, y)
        assert float(np.max(np.abs(r0))) < 1e-6, f"lambda=0 residual {np.max(np.abs(r0))}"
        lam, dlam, lam_done = 0.0, 0.1, 0.0
        steps = []
        while True:
            sysm.lam = lam
            sol = least_squares(
                lambda yy, s_=sysm, x0_=x0: date_chain_residual(s_, legs_chain, x0_, yy),
                y,
                method="lm",
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
                max_nfev=50 * len(y),
            )
            res = date_chain_residual(sysm, legs_chain, x0, sol.x)
            conv = bool(np.max(np.abs(res)) < 1e-6)
            cyc_now = Cycle(legs_chain, float(sol.x[-1]) * DAY)
            x_now = np.concatenate([[x0], sol.x[:-1] * DAY])
            gate = interior_gate(sysm, cyc_now, x_now) if conv else {"status": "unconverged"}
            rec = {
                "phase": "ramp" if ramp else "blend",
                "lambda": lam,
                "converged": conv,
                "max_residual_kms": float(np.max(np.abs(res))),
                "max_date_shift_d": float(np.max(np.abs(sol.x - y))),
                "chain_period_d": float(sol.x[-1]),
            } | gate
            steps.append(rec)
            print(
                f"{time.strftime('%H:%M:%S')} {args.cell} epoch {ie} (JD {te / DAY + 2440000:.1f}) "
                f"{rec['phase']} lam={lam:.4f} conv={conv} res={rec['max_residual_kms']:.1e} "
                f"gate={gate.get('status')} worst={gate.get('worst_ratio', float('nan')):.3f} "
                f"[elapsed {time.time() - t_run:.0f}s]",
                flush=True,
            )
            if conv:
                y = sol.x
                if lam >= 1.0:
                    break
                lam_done = lam
                lam = min(1.0, lam + dlam)
            else:
                dlam /= 2.0
                if dlam < 1.0 / 640:
                    break
                lam = lam_done + dlam
        conv_steps = [s_ for s_ in steps if s_["converged"]]
        reached = bool(conv_steps) and conv_steps[-1]["lambda"] >= 1.0
        rec_out: dict[str, Any] = {
            "epoch_jd": te / DAY + 2440000.0,
            "mode": "ramp" if ramp else "blend+shoot",
            "last_converged_lambda": max((s_["lambda"] for s_ in conv_steps), default=None),
            "steps": steps,
            "x0_days": x0 / DAY,
            "final_y_dates": y.tolist(),
        }
        if ramp:
            rec_out["rung_pass"] = reached and conv_steps[-1].get("status") == "pass"
        elif reached and any(not isinstance(lg, LambertLeg) for lg in legs_chain):
            # Phase 2 (blend mode, lambda = 1 = the real ephemeris): shoot every fixed leg,
            # dates as start, restarts over the fixed-leg directions; keep the closing solution
            # with the smallest worst gate ratio.
            sysm.lam = 1.0
            cyc_now = Cycle(legs_chain, float(y[-1]) * DAY)
            x_now = np.concatenate([[x0], y[:-1] * DAY])
            fixed0 = np.ravel(initial_fixed_params(sysm, cyc_now, x_now))
            rng = np.random.default_rng(943 + ie)
            best = None
            tried = []
            for r in range(args.shoot_restarts):
                f0 = fixed0.copy()
                if r > 0:
                    f0[0::3] += rng.normal(0.0, 0.3, f0[0::3].size)
                    f0[1::3] += rng.normal(0.0, 0.6, f0[1::3].size)
                ys = np.concatenate([y, f0])
                sol = least_squares(
                    lambda yy, s_=sysm, x0_=x0: _res(s_, legs_chain, x0_, yy),
                    ys,
                    method="lm",
                    xtol=1e-14,
                    ftol=1e-14,
                    gtol=1e-14,
                    max_nfev=30 * len(ys),
                )
                ev = chain_eval(sysm, legs_chain, x0, sol.x)
                if ev is None or float(np.max(np.abs(ev.residual))) >= 1e-6:
                    tried.append({"restart": r, "converged": False})
                    continue
                g = gate_eval(sysm, ev)
                tried.append({"restart": r, "converged": True} | g)
                if best is None or g["worst_ratio"] < best[0]["worst_ratio"]:
                    best = (g, sol.x)
                print(
                    f"{time.strftime('%H:%M:%S')}   shoot restart {r}: closes, "
                    f"gate={g['status']} worst={g['worst_ratio']:.3f}",
                    flush=True,
                )
            rec_out["shoot_restarts"] = tried
            rec_out["rung_pass"] = best is not None and best[0]["status"] == "pass"
            if best is not None:
                rec_out["shoot_best"] = best[0]
                rec_out["final_y"] = best[1].tolist()
        else:
            rec_out["rung_pass"] = reached and conv_steps[-1].get("status") == "pass"
            if reached:
                rec_out["final_y"] = y.tolist()
        out.append(rec_out)
        (args.out / "realeph_chain.json").write_text(json.dumps(out, indent=1, default=float))
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
