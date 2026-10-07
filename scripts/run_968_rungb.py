"""#968 rung (b): GanEur#316@2019 on jup365 rails in continuous gravity (note sec. 9, amendment 8).

Open chain of N cycles. Node 0 is ANCHORED: the patched-conic state one day after the first
departure from Ganymede (on the first Lambert leg; independent of sigma). Nodes 1..M are the
chain's flybys (periapsis gauge relative to each node's moon). Forward-backward mid-time matches.
Ganymede and Europa are point masses with GM and softening radius scaled by sigma; Io and Callisto
are off (step 8 switches them on at sigma = 1). Arcs: ``jovian_stm.propagate_with_stm`` (DOP853 +
analytic STM) on a cubic-spline wrapper of jup365 (validated against spkezr in --stage ephem).

    uv run python scripts/run_968_rungb.py --stage ephem
    uv run python scripts/run_968_rungb.py --stage sigma --n-cycles 1
    uv run python scripts/run_968_rungb.py --stage ias15 --n-cycles 1
    uv run python scripts/run_968_rungb.py --stage fourmoon --n-cycles 1

Checkpoints in data/968_rungb/; live log in data/968_rungb/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import (
    MU_JUPITER_KM3_S2,
    JovianEphemeris,
    JovianRailsCache,
    JovianRestrictedNBody,
    periapsis_node,
)
from cyclerfinder.nbody.jovian_stm import _accel_and_gradient, propagate_with_stm
from cyclerfinder.search.two_working_body import kepler_step

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "968_rungb"
DAY = 86400.0
OFF = (2440000.0 - 2451545.0) * DAY  # chain time (s past JD 2440000) -> lane TDB s past J2000
KERNEL = Path.home() / "dev" / "references" / "kernels" / "jup365.bsp"
GAN, EUR = "Ganymede", "Europa"
ALL4 = ("Io", "Europa", "Ganymede", "Callisto")
W = np.array([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
W_GAUGE = 1e4
RTOL, ATOL = 1e-13, 1e-12  # amendment 10
FLOOR_KM = {GAN: 100.0, EUR: 100.0}
SOI_LIMIT = 0.5  # rung (b) acceptance (amendment 8 item 6); callers may override
Arr = NDArray[np.float64]


def _log(msg: str) -> None:
    (OUT / "live").mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "live" / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


class SplineEphem:
    """Cubic-spline wrapper of jup365 with ``state()`` (velocity from the spline derivative)."""

    def __init__(self, t0: float, t1: float, step_days: float = 0.004) -> None:
        self.j = JovianEphemeris(str(KERNEL))
        n = int((t1 - t0) / (step_days * DAY)) + 1
        grid = np.linspace(t0, t1, n)
        self.sp: dict[str, CubicSpline] = {}
        for m in ALL4:
            pos = np.array([self.j.state(m, float(t))[0] for t in grid])
            self.sp[m] = CubicSpline(grid, pos, axis=0)

    def state(self, moon: str, t_sec: float) -> tuple[Arr, Arr]:
        s = self.sp[moon]
        return np.asarray(s(t_sec), dtype=np.float64), np.asarray(s(t_sec, 1), dtype=np.float64)

    def position(self, moon: str, t_sec: float) -> Arr:
        return np.asarray(self.sp[moon](t_sec), dtype=np.float64)


@dataclass
class Chain:
    eph: SplineEphem
    anchor: Arr  # 6-state of node 0
    t_anchor: float
    moons: list[str]  # moon of nodes 1..M
    seed: Arr  # z for sigma = 1 periapsis nodes (unscaled offsets)
    chain_vinf: list[float]  # the reconstructed chain's V_inf at nodes 1..M
    vin_first: Arr  # amendment 12: the chain's inbound V_inf vector at node 1
    vout_last: Arr  # and its outbound V_inf vector at node M
    force: tuple[str, ...] = (GAN, EUR)
    mus: dict[str, float] = field(default_factory=dict)
    surf: dict[str, float] = field(default_factory=dict)

    @property
    def m(self) -> int:
        return len(self.moons)

    def set_sigma(self, s: float) -> None:
        self.mus = {k: s * SATELLITES[k].mu_km3_s2 for k in self.force}
        self.surf = {k: s * SATELLITES[k].radius_eq_km for k in self.force}
        for k in self.force:
            if k not in (GAN, EUR):  # Io / Callisto (step 8) always at full mass
                self.mus[k] = SATELLITES[k].mu_km3_s2
                self.surf[k] = SATELLITES[k].radius_eq_km


def build(n_cycles: int) -> Chain:
    d = json.loads((OUT / "seed_chain.json").read_text())
    fl = d["flybys"][: 5 * n_cycles]
    x0 = d["x0_s_past_jd2440000"] + OFF
    t_end = fl[-1]["t_s_past_jd2440000"] + OFF
    eph = SplineEphem(x0 - 2 * DAY, t_end + 2 * DAY)
    rg, vg = eph.j.state(GAN, x0)
    v0 = vg + np.asarray(d["first_departure"]["vinf_out"])
    ra, va = kepler_step(rg, v0, DAY, MU_JUPITER_KM3_S2)
    z = []
    for f in fl:
        t = f["t_s_past_jd2440000"] + OFF
        r, v, _ = periapsis_node(
            f["body"], t, np.asarray(f["vinf_in"]), np.asarray(f["vinf_out"]), eph.j
        )
        rm, vm = eph.state(f["body"], t)
        z.extend([*(r - rm), *(v - vm), t])  # moon-relative node coordinates
    return Chain(
        eph,
        np.concatenate([ra, va]),
        x0 + DAY,
        [f["body"] for f in fl],
        np.asarray(z, dtype=np.float64),
        [f["vinf_kms"] for f in fl],
        np.asarray(fl[0]["vinf_in"], dtype=np.float64),
        np.asarray(fl[-1]["vinf_out"], dtype=np.float64),
    )


def scale_offsets(c: Chain, z: Arr, f: float) -> Arr:
    """Scale every node's moon-relative position by f (patched-conic r_p is linear in GM)."""
    z = z.copy()
    for k in range(c.m):
        z[7 * k : 7 * k + 3] *= f
    return z


def nodes(c: Chain, z: Arr) -> list[tuple[Arr, Arr, float, Arr]]:
    n = 7 * c.m
    out: list[tuple[Arr, Arr, float, Arr]] = []  # amendment 9: no anchor
    for k in range(c.m):
        # Unknowns are moon-relative (r_rel, v_rel) and the epoch: x = moon(t) + rel, so an
        # epoch change carries the node with its moon (the absolute form is stiff at small sigma).
        t = float(z[7 * k + 6])
        rm, vm = c.eph.state(c.moons[k], t)
        h = 1.0
        am = (c.eph.state(c.moons[k], t + h)[1] - c.eph.state(c.moons[k], t - h)[1]) / (2 * h)
        dx = np.zeros((6, n))
        dx[:, 7 * k : 7 * k + 6] = np.eye(6)
        dx[0:3, 7 * k + 6] = vm
        dx[3:6, 7 * k + 6] = am
        dt = np.zeros(n)
        dt[7 * k + 6] = 1.0
        out.append(
            (np.concatenate([rm + z[7 * k : 7 * k + 3], vm + z[7 * k + 3 : 7 * k + 6]]), dx, t, dt)
        )
    return out


def accel(c: Chain, x: Arr, t: float) -> Arr:
    a, _ = _accel_and_gradient(x[:3], t, ephem=c.eph, moons=c.force, mus=c.mus, surf=c.surf)  # type: ignore[arg-type]
    return np.concatenate([x[3:], a])


def arc(c: Chain, x: Arr, t0: float, t1: float) -> tuple[Arr, Arr]:
    rf, vf, phi = propagate_with_stm(
        x[:3],
        x[3:],
        t0,
        t1,
        ephem=c.eph,  # type: ignore[arg-type]
        moons=c.force,
        rtol=RTOL,
        atol=ATOL,
        mu_overrides=c.mus,
        radius_overrides=c.surf,
    )
    return np.concatenate([rf, vf]), phi


def gauge_and_grad(c: Chain, moon: str, x: Arr, t: float) -> tuple[float, Arr, float]:
    rm, vm = c.eph.state(moon, t)
    h = 1.0
    am = (c.eph.state(moon, t + h)[1] - c.eph.state(moon, t - h)[1]) / (2 * h)
    dr, dv = x[:3] - rm, x[3:] - vm
    a, b = float(np.linalg.norm(dr)), float(np.linalg.norm(dv))
    g = float(dr @ dv) / (a * b)
    g_dr = dv / (a * b) - g * dr / a**2
    g_dv = dr / (a * b) - g * dv / b**2
    return g, np.concatenate([g_dr, g_dv]), float(-(g_dr @ vm) - (g_dv @ am))


def asymptote(mu: float, r: Arr, v: Arr, which: str) -> Arr:
    """Inbound or outbound V_inf vector of the moon-relative two-body hyperbola through (r, v)."""
    h = np.cross(r, v)
    ev = np.cross(v, h) / mu - r / float(np.linalg.norm(r))
    e = float(np.linalg.norm(ev))
    vinf = math.sqrt(max(float(v @ v) - 2.0 * mu / float(np.linalg.norm(r)), 0.0))
    p_hat = ev / e
    q_hat = np.cross(h / float(np.linalg.norm(h)), p_hat)
    nu = math.acos(-1.0 / e)
    if which == "out":
        return vinf * (math.cos(nu) * p_hat + math.sin(nu) * q_hat)
    return -vinf * (math.cos(nu) * p_hat - math.sin(nu) * q_hat)


END_MODE = "vector"  # "vector": #968 amendment 12; "direction": #1039 formulation (b)


def _perp(u: Arr) -> tuple[Arr, Arr]:
    a = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, a)
    e1 /= float(np.linalg.norm(e1))
    return e1, np.cross(u, e1)


def end_rows_direction(c: Chain, z: Arr) -> tuple[Arr, Arr]:
    """#1039 (b): end V_inf directions pinned (2 + 2 rows), |V_inf in, first| = |V_inf out, last|
    (1 row), and the chain span t_last - t_first equal to the seed's (1 row). Rows are scaled so
    that a row value of 1e-3 is 1e-6 rad, 1e-6 km/s, or 1e-3 s."""
    n = 7 * c.m
    k_last = c.m - 1
    span0 = float(c.seed[7 * k_last + 6] - c.seed[6])
    u_in = c.vin_first / float(np.linalg.norm(c.vin_first))
    u_out = c.vout_last / float(np.linalg.norm(c.vout_last))
    p_in, p_out = _perp(u_in), _perp(u_out)
    mu0, mu1 = c.mus[c.moons[0]], c.mus[c.moons[k_last]]
    i0, i1 = 7 * k_last, 7 * k_last + 3

    def f(zz: Arr) -> Arr:
        vi = asymptote(mu0, zz[0:3], zz[3:6], "in")
        vo = asymptote(mu1, zz[i0:i1], zz[i1 : i1 + 3], "out")
        ui, uo = vi / float(np.linalg.norm(vi)), vo / float(np.linalg.norm(vo))
        return 1e3 * np.array(
            [
                float(p_in[0] @ ui),
                float(p_in[1] @ ui),
                float(p_out[0] @ uo),
                float(p_out[1] @ uo),
                float(np.linalg.norm(vi) - np.linalg.norm(vo)),
                (float(zz[7 * k_last + 6] - zz[6]) - span0) * 1e-3,
            ]
        )

    res = f(z)
    jac = np.zeros((6, n))
    for j in [*range(0, 7), *range(7 * k_last, 7 * k_last + 7)]:
        jj = j - (0 if j < 7 else 7 * k_last)
        hstep = 1e-4 * max(1.0, abs(float(z[j]))) if jj < 3 else (1e-8 if jj < 6 else 1e-2)
        zp, zm = z.copy(), z.copy()
        zp[j] += hstep
        zm[j] -= hstep
        jac[:, j] = (f(zp) - f(zm)) / (2 * hstep)
    return res, jac


def end_rows(c: Chain, z: Arr, nd: list[tuple[Arr, Arr, float, Arr]]) -> tuple[Arr, Arr]:
    """Amendment 12: inbound asymptote at node 1, outbound at node M (6 rows, weight 1e3)."""
    if END_MODE == "direction":
        return end_rows_direction(c, z)
    n = 7 * c.m
    res = np.zeros(6)
    jac = np.zeros((6, n))
    for row0, k, which, target in ((0, 0, "in", c.vin_first), (3, c.m - 1, "out", c.vout_last)):
        mu = c.mus[c.moons[k]]

        def f(zz: Arr, k: int = k, which: str = which, mu: float = mu) -> Arr:
            return asymptote(mu, zz[7 * k : 7 * k + 3], zz[7 * k + 3 : 7 * k + 6], which)

        res[row0 : row0 + 3] = 1e3 * (f(z) - target)
        for j in range(7 * k, 7 * k + 6):
            hstep = 1e-4 * max(1.0, abs(float(z[j]))) if (j - 7 * k) < 3 else 1e-8
            zp, zm = z.copy(), z.copy()
            zp[j] += hstep
            zm[j] -= hstep
            jac[row0 : row0 + 3, j] = 1e3 * (f(zp) - f(zm)) / (2 * hstep)
    return res, jac


def residual_and_jac(c: Chain, z: Arr, want_jac: bool = True) -> tuple[Arr, Arr, dict[str, Any]]:
    nd = nodes(c, z)
    res: list[float] = []
    rows: list[Arr] = []
    info: dict[str, Any] = {"dr": [], "dv": [], "gauge": []}
    for k in range(c.m - 1):
        xp, dxp, tp, dtp = nd[k]
        xn, dxn, tn, dtn = nd[k + 1]
        tm = 0.5 * (tp + tn)
        try:
            xf, pf = arc(c, xp, tp, tm)
            xb, pb = arc(c, xn, tn, tm)
        except RuntimeError:
            if want_jac:
                raise
            # a trial step into a moon's core: an honest failed evaluation, rejected by the caller
            return (
                np.full(7 * c.m, np.inf),
                np.zeros((0, 7 * c.m)),
                {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf]},
            )
        d = xf - xb
        info["dr"].append(float(np.linalg.norm(d[:3])))
        info["dv"].append(float(np.linalg.norm(d[3:])))
        res.extend(W * d)
        if want_jac:
            ff, fb = accel(c, xf, tm), accel(c, xb, tm)
            dxf = (
                pf @ dxp
                + np.outer(-pf @ accel(c, xp, tp) + 0.5 * ff, dtp)
                + np.outer(0.5 * ff, dtn)
            )
            dxb = (
                pb @ dxn
                + np.outer(-pb @ accel(c, xn, tn) + 0.5 * fb, dtn)
                + np.outer(0.5 * fb, dtp)
            )
            rows.append(W[:, None] * (dxf - dxb))
    for k in range(c.m):
        x, dx, t, dt = nd[k]
        g, gx, gt = gauge_and_grad(c, c.moons[k], x, t)
        info["gauge"].append(g)
        res.append(W_GAUGE * g)
        if want_jac:
            rows.append(W_GAUGE * (gx @ dx + gt * dt)[None, :])
    er, ej = end_rows(c, z, nd)
    res.extend(er.tolist())
    info["end_asymptote_kms"] = [
        float(np.max(np.abs(er[:3]))) / 1e3,
        float(np.max(np.abs(er[3:]))) / 1e3,
    ]
    if want_jac:
        rows.append(ej)
    return np.asarray(res), (np.vstack(rows) if want_jac else np.zeros((0, 7 * c.m))), info


def converged(info: dict[str, Any]) -> bool:
    return (
        max(info["dr"]) < 1e-3
        and max(info["dv"]) < 1e-6
        and max(abs(g) for g in info["gauge"]) < 1e-9
    )


def newton(c: Chain, z: Arr, iters: int, tag: str) -> tuple[Arr, dict[str, Any]]:
    info: dict[str, Any] = {}
    for _ in range(iters):
        try:
            r, j, info = residual_and_jac(c, z)
        except RuntimeError:  # a predictor or step into a moon's core: a failed solve
            return z, {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf]}
        _log(f"{tag} |r| {np.linalg.norm(r):.3e} dr {max(info['dr']):.2e} dv {max(info['dv']):.2e}")
        if converged(info):
            return z, info
        d = np.linalg.norm(j, axis=0)
        d[d == 0.0] = 1.0
        step = np.linalg.lstsq(j / d, -r, rcond=None)[0] / d
        # Step cap: no node's moon-relative position moves by more than half its length.
        cap = 1.0
        for k in range(c.m):
            rr = float(np.linalg.norm(z[7 * k : 7 * k + 3]))
            dd = float(np.linalg.norm(step[7 * k : 7 * k + 3]))
            if dd > 0.5 * rr:
                cap = min(cap, 0.5 * rr / dd)
        step = cap * step
        f0 = float(np.linalg.norm(r))
        for alpha in [0.5**i for i in range(11)]:  # backtracking down to 1e-3
            rn = residual_and_jac(c, z + alpha * step, False)[0]
            if np.all(np.isfinite(rn)) and float(np.linalg.norm(rn)) < f0:
                z = z + alpha * step
                break
        else:
            break
    _, _, info = residual_and_jac(c, z, False)
    return z, info


def describe(c: Chain, z: Arr, sigma: float) -> dict[str, Any]:
    nd = nodes(c, z)
    out = []
    for k in range(c.m):
        moon = c.moons[k]
        x, _, t, _ = nd[k]
        rm, vm = c.eph.state(moon, t)
        d = float(np.linalg.norm(x[:3] - rm))
        mu = c.mus[moon]
        vinf2 = float(np.linalg.norm(x[3:] - vm)) ** 2 - 2.0 * mu / d
        soi = SATELLITES[moon].sma_km * (mu / MU_JUPITER_KM3_S2) ** 0.4 if mu > 0 else 0.0
        out.append(
            {
                "moon": moon,
                "t_days_from_x0": (t - c.t_anchor + DAY) / DAY,
                "vinf_kms": math.sqrt(max(vinf2, 0.0)),
                "rp_km": d,
                "alt_km": d - sigma * SATELLITES[moon].radius_eq_km,
                "rp_over_soi": d / soi if soi > 0 else float("nan"),
                "floor_ok": d - sigma * SATELLITES[moon].radius_eq_km >= sigma * FLOOR_KM[moon],
            }
        )
    return {"nodes": out}


def unscheduled(c: Chain, z: Arr, sigma: float) -> list[dict[str, Any]]:
    """Distance minima to all four moons over the chain, other than the node encounters."""
    nd = nodes(c, z)
    found = []
    for k in range(c.m - 1):
        x0, _, t0, _ = nd[k]
        t1 = nd[k + 1][2]
        sol = solve_ivp(
            lambda t, y: accel(c, y, t),
            (t0, t1),
            x0,
            method="DOP853",
            rtol=RTOL,
            atol=ATOL,
            dense_output=True,
        )
        ts = np.linspace(t0, t1, int((t1 - t0) / (0.002 * DAY)) + 2)
        ys = sol.sol(ts)
        for moon in ALL4:
            mu = SATELLITES[moon].mu_km3_s2 * (sigma if moon in (GAN, EUR) else 1.0)
            hill = SATELLITES[moon].sma_km * (mu / (3 * MU_JUPITER_KM3_S2)) ** (1 / 3)
            dist = np.array(
                [
                    np.linalg.norm(ys[:3, i] - c.eph.position(moon, float(t)))
                    for i, t in enumerate(ts)
                ]
            )
            for i in range(1, len(ts) - 1):
                if dist[i] <= dist[i - 1] and dist[i] <= dist[i + 1] and dist[i] < hill:
                    found.append(
                        {
                            "leg": k,
                            "moon": moon,
                            "t_days": float(ts[i] / DAY),
                            "dist_km": float(dist[i]),
                            "hill_km": hill,
                        }
                    )
    # the scheduled encounters sit at the leg ends; keep only interior minima
    return found


def stage_ephem() -> None:
    c = build(1)
    worst_r = worst_v = 0.0
    for t in np.linspace(c.t_anchor, c.t_anchor + 45 * DAY, 20):
        for m in ALL4:
            r1, v1 = c.eph.j.state(m, float(t))
            r2, v2 = c.eph.state(m, float(t))
            worst_r = max(worst_r, float(np.linalg.norm(r1 - r2)))
            worst_v = max(worst_v, float(np.linalg.norm(v1 - v2)))
    ok = worst_r < 0.1 and worst_v < 1e-6
    _log(f"ephem: spline vs spkezr worst |dr| {worst_r:.3e} km |dv| {worst_v:.3e} km/s pass={ok}")
    (OUT / "ephem_check.json").write_text(
        json.dumps({"worst_dr_km": worst_r, "worst_dv_kms": worst_v, "pass": ok})
    )


def stage_sigma(n_cycles: int) -> None:
    c = build(n_cycles)
    path = OUT / f"sigma_n{n_cycles}.json"
    rec = json.loads(path.read_text()) if path.exists() else {"points": []}
    pts = rec["points"]
    if not pts:
        sg = 0.02  # amendment 10
        c.set_sigma(sg)
        z, info = newton(c, scale_offsets(c, c.seed, sg), 40, f"n{n_cycles} sigma={sg}")
        if not converged(info):
            _log("sigma: no convergence at 0.02; stop")
            return
        pts.append({"sigma": sg, "z": z.tolist(), **describe(c, z, sg)})
        path.write_text(json.dumps(rec))
    fac = rec.get("fac", 1.2)
    while pts[-1]["sigma"] < 1.0:
        s1, z1 = pts[-1]["sigma"], np.asarray(pts[-1]["z"])
        sg = min(1.0, s1 * fac)
        if sg - s1 < 1e-5:
            _log(f"sigma: step below 1e-5 at {s1:.6f}; stop")
            rec["stopped_at"] = s1
            break
        if len(pts) >= 2:
            s0, z0 = pts[-2]["sigma"], np.asarray(pts[-2]["z"])
            zs = z1 + (z1 - z0) * (sg - s1) / (s1 - s0)  # linear in sigma (r_rel ~ sigma)
        else:
            zs = scale_offsets(c, z1, sg / s1)
        c.set_sigma(sg)
        z, info = newton(c, zs, 25, f"n{n_cycles} sigma={sg:.6g}")
        floor_ok = bool(converged(info))
        noise_ok = sg < 1.0 and max(info["dr"]) < 1e-2 and max(info["dv"]) < 1e-6  # amendment 11
        if floor_ok or noise_ok:
            desc = describe(c, z, sg)
            bad = [n for n in desc["nodes"] if not n["floor_ok"] or n["rp_over_soi"] > SOI_LIMIT]
            if bad:
                _log(f"sigma={sg:.6g}: acceptance FAILED at nodes {bad}; stop")
                rec["acceptance_fail"] = {"sigma": sg, "nodes": bad}
                break
            pts.append(
                {
                    "sigma": sg,
                    "z": z.tolist(),
                    **desc,
                    "noise_floor_limited": not floor_ok,
                    "dr": max(info["dr"]),
                    "dv": max(info["dv"]),
                }
            )
            fac = min(1.2, 1.0 + 2.0 * (fac - 1.0))
            _log(
                f"sigma={sg:.6g} ok (noise-floor-limited={not floor_ok}) "
                f"vinf {[round(n['vinf_kms'], 4) for n in desc['nodes']]}"
            )
        else:
            fac = 1.0 + 0.5 * (fac - 1.0)
            _log(f"sigma={sg:.6g} failed; factor -> {fac:.6g}")
        rec["fac"] = fac
        path.write_text(json.dumps(rec))
    if pts[-1]["sigma"] >= 1.0:
        z = np.asarray(pts[-1]["z"])
        c.set_sigma(1.0)
        rec["unscheduled_at_1"] = unscheduled(c, z, 1.0)
        _log(f"sigma reached 1; unscheduled minima inside a Hill radius: {rec['unscheduled_at_1']}")
    path.write_text(json.dumps(rec))


def stage_down(n_cycles: int) -> None:
    """Amendment 10, criterion 7 (ii): from sigma = 0.02 down to 0.01 and 0.005."""
    c = build(n_cycles)
    path = OUT / f"sigma_n{n_cycles}.json"
    rec = json.loads(path.read_text())
    q = rec["points"][0]
    assert abs(q["sigma"] - 0.02) < 1e-12
    z, s_prev = np.asarray(q["z"]), 0.02
    down = []
    for sg in (0.01, 0.005):
        c.set_sigma(sg)
        zn, info = newton(c, scale_offsets(c, z, sg / s_prev), 40, f"n{n_cycles} down sigma={sg}")
        ok = bool(converged(info))
        row: dict[str, Any] = {
            "sigma": sg,
            "converged": ok,
            "dr": max(info["dr"]),
            "dv": max(info["dv"]),
        }
        if ok:
            desc = describe(c, zn, sg)
            row.update(desc)
            row["z"] = zn.tolist()
            row["max_dvinf_vs_chain"] = max(
                abs(n["vinf_kms"] - v) for n, v in zip(desc["nodes"], c.chain_vinf, strict=True)
            )
            z, s_prev = zn, sg
        down.append(row)
        _log(f"down sigma={sg}: converged={ok} max |dV_inf| {row.get('max_dvinf_vs_chain')}")
        if not ok:
            break
    good = [r for r in down if r["converged"]]
    rec["down"] = down
    rec["criterion_ii"] = bool(good) and good[-1]["max_dvinf_vs_chain"] < 0.01
    path.write_text(json.dumps(rec))
    _log(f"criterion (ii): {rec['criterion_ii']}")


def stage_ias15(n_cycles: int) -> None:
    c = build(n_cycles)
    rec = json.loads((OUT / f"sigma_n{n_cycles}.json").read_text())
    q = rec["points"][-1]
    assert q["sigma"] >= 1.0
    z = np.asarray(q["z"])
    c.set_sigma(1.0)
    nd = nodes(c, z)
    lo = min(t for _, _, t, _ in nd) - DAY
    hi = max(t for _, _, t, _ in nd) + DAY
    cache = JovianRailsCache((GAN, EUR), c.eph.j, lo, hi)
    prop = JovianRestrictedNBody()
    rows = []
    for k in range(c.m - 1):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for lab, x0, t0 in (("fwd", xp, tp), ("bwd", xn, tn)):
            xd, _ = arc(c, x0, t0, tm)
            a = prop.propagate(
                x0[:3], x0[3:], t0, tm, moons=(GAN, EUR), cache=cache, max_wall_sec=3000.0
            )
            if not a.converged or abs(a.t1_sec - tm) > 1e-6:
                rows.append({"leg": k, "dir": lab, "status": "timeout"})
                continue
            rows.append(
                {
                    "leg": k,
                    "dir": lab,
                    "status": "ok",
                    "dr_km": float(np.linalg.norm(a.r_km - xd[:3])),
                    "dv_kms": float(np.linalg.norm(a.v_km_s - xd[3:])),
                }
            )
    ok = all(r["status"] == "ok" and r["dr_km"] < 1e-2 and r["dv_kms"] < 1e-7 for r in rows)
    rec["ias15"] = {"rows": rows, "pass": ok}
    (OUT / f"sigma_n{n_cycles}.json").write_text(json.dumps(rec))
    _log(
        f"ias15 n{n_cycles}: pass={ok} max dr {max(r.get('dr_km', float('nan')) for r in rows):.2e}"
    )


def stage_fourmoon(n_cycles: int) -> None:
    c = build(n_cycles)
    rec = json.loads((OUT / f"sigma_n{n_cycles}.json").read_text())
    z = np.asarray(rec["points"][-1]["z"])
    c.force = ALL4
    c.set_sigma(1.0)
    zn, info = newton(c, z, 25, f"n{n_cycles} four-moon")
    out: dict[str, Any] = {
        "converged": bool(converged(info)),
        "dr": max(info["dr"]),
        "dv": max(info["dv"]),
    }
    if converged(info):
        out.update(describe(c, zn, 1.0))
        out["z"] = zn.tolist()
    rec["fourmoon"] = out
    (OUT / f"sigma_n{n_cycles}.json").write_text(json.dumps(rec))
    _log(f"four-moon n{n_cycles}: converged={out['converged']} dv {out['dv']:.2e}")


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="ganeur316-rs2019-jup365-continuous-rungb",
        method=MethodCapability(
            genome="GanEur#316 open chain on jup365 (one structure)",
            corrector="forward-backward shooting, DOP853 + STM; sigma continuation",
            capability_tags=frozenset({"ballistic", "n-body", "real-ephemeris"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage", choices=["ephem", "sigma", "down", "ias15", "fourmoon"], required=True
    )
    ap.add_argument("--n-cycles", type=int, default=1)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {
        "ephem": lambda: stage_ephem(),
        "sigma": lambda: stage_sigma(args.n_cycles),
        "down": lambda: stage_down(args.n_cycles),
        "ias15": lambda: stage_ias15(args.n_cycles),
        "fourmoon": lambda: stage_fourmoon(args.n_cycles),
    }[args.stage]()


if __name__ == "__main__":
    main()
