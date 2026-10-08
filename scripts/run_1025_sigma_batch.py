"""#1025: ideal-model continuous-gravity batch for the clean #973 Ganymede-Callisto members.

Pre-registration: docs/notes/2026-10-08-1025-sigma-batch.md. Per member: the #1034 joint-sigma
continuation, generalised to any node count. Closed chain in the R-S 2009 circular model (cell gc,
both moons massive), periapsis nodes at every flyby (moon-relative coordinates, epoch, periapsis
gauge), forward-backward mid-leg matches, rotation-aware wrap (Ganymede's advance over T = k S_GC),
DOP853 + analytic STM (rtol 1e-13), Ganymede and Callisto GMs and radii scaled by sigma.

Per member: converge at sigma = 0.02; identity by the limit going DOWN (0.01, 0.005, 0.0025; gate);
natural continuation UP (factor 1.2, secant in log sigma, halving to a 1e-5 step); on a stop,
monitor-variable continuation (no SVD tangent: the enforced tangent rule), then the fold checks
(two solutions at one sigma, direct Newton at sigma = 1) or a resume if sigma keeps rising; impact
when a periapsis is inside sigma x R; at sigma = 1 the IAS15 re-fly of every half-arc. Outcome per
member: EXISTS / FOLDS at sigma_f / IMPACT at sigma_i / IDENTITY FAIL / NUMERICAL STOP.

    uv run python scripts/run_1025_sigma_batch.py --members gc-1          # control (#1034)
    uv run python scripts/run_1025_sigma_batch.py --members gc-2          # control (#1004 fold)
    uv run python scripts/run_1025_sigma_batch.py --members batch --worker 0/2   # the 20 members
    uv run python scripts/run_1025_sigma_batch.py --summary

Resumable per member and per point (data/1025_sigma/<id>.json); a call stops cleanly after
--max-wall seconds (default 420). Live log data/1025_sigma/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2, periapsis_node
from cyclerfinder.nbody.jovian_stm import _accel_and_gradient, propagate_with_stm
from cyclerfinder.search.two_working_body import CircularSystem, cycle_flybys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1025_sigma"
DAY = 86400.0
GAN, CAL = "Ganymede", "Callisto"
MOONS = (GAN, CAL)
W = np.array([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
W_GAUGE = 1e4
RTOL, ATOL = 1e-13, 1e-12
FLOOR_KM = {GAN: 100.0, CAL: 200.0}
Arr = NDArray[np.float64]

# (id, k, key, V_inf G, V_inf C) from the #973 note secs. 10-12 (clean members) and the controls.
MEMBERS: list[tuple[str, str, float, float]] = [
    (
        "gc-1",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        2.397,
        1.807,
    ),
    (
        "gc-2",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s",
        3.617,
        3.039,
    ),
    (
        "gc4-1",
        "k4|RGanymede/3:2|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        1.4719,
        1.5814,
    ),
    (
        "gc4-2",
        "k4|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        1.7380,
        1.8572,
    ),
    ("gc4-6", "k4|RGanymede/3:2|LGanymede>Callisto/2l|LCallisto>Ganymede/0s", 5.0073, 2.1348),
    (
        "gc4-8",
        "k4|LGanymede>Ganymede/1h|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s",
        2.4191,
        2.4855,
    ),
    ("gc4-10", "k4|RGanymede/1:1|LGanymede>Callisto/2l|LCallisto>Ganymede/2h", 6.4468, 2.7512),
    (
        "gc4-19",
        "k4|RGanymede/2:1|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        2.3565,
        4.1090,
    ),
    (
        "gc4-25",
        "k4|LGanymede>Ganymede/1h|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s",
        8.1525,
        5.2279,
    ),
    (
        "gc4-28",
        "k4|LGanymede>Ganymede/2h|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s",
        8.8570,
        5.5509,
    ),
    (
        "gc4-38",
        "k4|LGanymede>Ganymede/1h|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/1l",
        11.4600,
        7.2168,
    ),
    (
        "gc4-41",
        "k4|LGanymede>Ganymede/2l|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        11.8283,
        7.5126,
    ),
    (
        "gc5-0",
        "k5|LGanymede>Ganymede/2l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/1h",
        1.490,
        1.255,
    ),
    (
        "gc5-3",
        "k5|LGanymede>Ganymede/2l|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        1.536,
        1.831,
    ),
    (
        "gc5-4",
        "k5|RGanymede/3:2|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        1.745,
        1.902,
    ),
    ("gc5-5", "k5|LGanymede>Callisto/0s|LCallisto>Callisto/2l|LCallisto>Ganymede/0s", 1.936, 2.760),
    (
        "gc5-6",
        "k5|RGanymede/2:1|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        2.191,
        3.206,
    ),
    (
        "gc5-7",
        "k5|LGanymede>Ganymede/1h|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/1l",
        4.701,
        4.597,
    ),
    ("gc6-0", "k6|LGanymede>Ganymede/2h|LGanymede>Callisto/2h|LCallisto>Ganymede/2h", 1.673, 1.257),
    (
        "gc6-1",
        "k6|LGanymede>Ganymede/2h|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        1.607,
        1.837,
    ),
    (
        "gc6-2",
        "k6|LGanymede>Ganymede/2l|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/0s",
        2.468,
        2.012,
    ),
    (
        "gc6-3",
        "k6|LGanymede>Ganymede/1h|LGanymede>Callisto/0s|LCallisto>Callisto/1l|LCallisto>Ganymede/1h",
        1.664,
        2.192,
    ),
]
CONTROLS = ("gc-1", "gc-2")


def _log(msg: str) -> None:
    (OUT / "live").mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "live" / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


def _enum() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_942_enumerate", ROOT / "scripts" / "run_942_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_942_enumerate"] = mod
    spec.loader.exec_module(mod)
    return mod


def rot_z(a: float) -> Arr:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


@dataclass
class P:
    circ: CircularSystem
    period: float
    q6: Arr
    moons: list[str]  # node moons in time order
    seed: Arr  # moon-relative node coordinates + epochs (sigma = 1 offsets)
    pc_vinf: list[float]  # patched-conic V_inf per node
    mus: dict[str, float] = field(default_factory=dict)
    surf: dict[str, float] = field(default_factory=dict)

    @property
    def m(self) -> int:
        return len(self.moons)

    def set_sigma(self, s: float) -> None:
        self.mus = {k: s * SATELLITES[k].mu_km3_s2 for k in MOONS}
        self.surf = {k: s * SATELLITES[k].radius_eq_km for k in MOONS}


def gauntlet_candidate(key: str, vg: float, vc: float) -> tuple[int, dict[str, Any]]:
    k = int(key.split("|")[0][1:])
    path = ROOT / "data" / ("943_cell_gc_gauntlet.json" if k == 3 else f"973_gc/k{k}_gauntlet.json")
    cands = json.loads(path.read_text())["candidates"]
    for c in cands:
        if (
            (key == c["key"] or key in c.get("keys", []))
            and abs(c["vinf_kms"][GAN] - vg) < 2e-3
            and abs(c["vinf_kms"][CAL] - vc) < 2e-3
        ):
            return k, c
    raise KeyError(f"no gauntlet candidate for {key} at V_inf {vg}/{vc}")


def build(member: str) -> P:
    _, key, vg, vc = next(m for m in MEMBERS if m[0] == member)
    _, cand = gauntlet_candidate(key, vg, vc)
    enum = _enum()
    circ, a, b = enum.cell_system("gc")
    _, cyc = enum.parse_cycle_key(cand["key"], circ, a, b)
    fl = cycle_flybys(circ, cyc, np.asarray(cand["x_days"]) * DAY)
    assert fl is not None
    fl = sorted(fl, key=lambda f: f.t_s)
    q = rot_z(2.0 * math.pi * cyc.period_s / circ.period_s(GAN))
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    z = []
    for f in fl:
        r, v, _ = periapsis_node(f.body, f.t_s, f.vinf_in, f.vinf_out, circ, max_offset_km=1.0e12)  # type: ignore[arg-type]
        rm, vm = circ.state(f.body, f.t_s)
        z.extend([*(r - rm), *(v - vm), f.t_s])
    p = P(
        circ,
        cyc.period_s,
        q6,
        [f.body for f in fl],
        np.asarray(z, dtype=np.float64),
        [f.vinf_kms for f in fl],
    )
    p.set_sigma(1.0)
    return p


def nodes(p: P, z: Arr) -> list[tuple[Arr, Arr, float, Arr]]:
    n = 7 * p.m
    out = []
    for k in range(p.m):
        t = float(z[7 * k + 6])
        rm, vm = p.circ.state(p.moons[k], t)
        nm = 2.0 * math.pi / p.circ.period_s(p.moons[k])
        am = -(nm**2) * rm
        dx = np.zeros((6, n))
        dx[:, 7 * k : 7 * k + 6] = np.eye(6)
        dx[0:3, 7 * k + 6] = vm
        dx[3:6, 7 * k + 6] = am
        dt = np.zeros(n)
        dt[7 * k + 6] = 1.0
        out.append(
            (np.concatenate([rm + z[7 * k : 7 * k + 3], vm + z[7 * k + 3 : 7 * k + 6]]), dx, t, dt)
        )
    x0, dx0, t0, dt0 = out[0]
    out.append((p.q6 @ x0, p.q6 @ dx0, t0 + p.period, dt0.copy()))
    return out


def accel(p: P, x: Arr, t: float) -> Arr:
    a, _ = _accel_and_gradient(x[:3], t, ephem=p.circ, moons=MOONS, mus=p.mus, surf=p.surf)  # type: ignore[arg-type]
    return np.concatenate([x[3:], a])


def arc(p: P, x: Arr, t0: float, t1: float) -> tuple[Arr, Arr]:
    rf, vf, phi = propagate_with_stm(
        x[:3],
        x[3:],
        t0,
        t1,
        ephem=p.circ,
        moons=MOONS,
        rtol=RTOL,
        atol=ATOL,
        mu_overrides=p.mus,
        radius_overrides=p.surf,  # type: ignore[arg-type]
    )
    return np.concatenate([rf, vf]), phi


def gauge_grad(p: P, moon: str, x: Arr, t: float) -> tuple[float, Arr, float]:
    rm, vm = p.circ.state(moon, t)
    nm = 2.0 * math.pi / p.circ.period_s(moon)
    am = -(nm**2) * rm
    dr, dv = x[:3] - rm, x[3:] - vm
    a, b = float(np.linalg.norm(dr)), float(np.linalg.norm(dv))
    g = float(dr @ dv) / (a * b)
    g_dr = dv / (a * b) - g * dr / a**2
    g_dv = dr / (a * b) - g * dv / b**2
    return g, np.concatenate([g_dr, g_dv]), float(-(g_dr @ vm) - (g_dv @ am))


def residual_and_jac(p: P, z: Arr, want_jac: bool = True) -> tuple[Arr, Arr, dict[str, Any]]:
    nd = nodes(p, z)
    res: list[float] = []
    rows: list[Arr] = []
    info: dict[str, Any] = {"dr": [], "dv": [], "gauge": []}
    for k in range(p.m):
        xp, dxp, tp, dtp = nd[k]
        xn, dxn, tn, dtn = nd[k + 1]
        tm = 0.5 * (tp + tn)
        try:
            xf, pf = arc(p, xp, tp, tm)
            xb, pb = arc(p, xn, tn, tm)
        except RuntimeError:
            if want_jac:
                raise
            return (
                np.full(7 * p.m, np.inf),
                np.zeros((0, 7 * p.m)),
                {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf]},
            )
        d = xf - xb
        info["dr"].append(float(np.linalg.norm(d[:3])))
        info["dv"].append(float(np.linalg.norm(d[3:])))
        res.extend(W * d)
        if want_jac:
            ff, fb = accel(p, xf, tm), accel(p, xb, tm)
            dxf = (
                pf @ dxp
                + np.outer(-pf @ accel(p, xp, tp) + 0.5 * ff, dtp)
                + np.outer(0.5 * ff, dtn)
            )
            dxb = (
                pb @ dxn
                + np.outer(-pb @ accel(p, xn, tn) + 0.5 * fb, dtn)
                + np.outer(0.5 * fb, dtp)
            )
            rows.append(W[:, None] * (dxf - dxb))
    for k in range(p.m):
        x, dx, t, dt = nd[k]
        g, gx, gt = gauge_grad(p, p.moons[k], x, t)
        info["gauge"].append(g)
        res.append(W_GAUGE * g)
        if want_jac:
            rows.append(W_GAUGE * (gx @ dx + gt * dt)[None, :])
    return np.asarray(res), (np.vstack(rows) if want_jac else np.zeros((0, 7 * p.m))), info


def converged(info: dict[str, Any]) -> bool:
    return (
        max(info["dr"]) < 1e-3
        and max(info["dv"]) < 1e-6
        and max(abs(g) for g in info["gauge"]) < 1e-9
    )


def noise_ok(info: dict[str, Any]) -> bool:
    """The noise-floor rule (#968): intermediate points (sigma < 1) accepted at dr < 1e-2 km."""
    return (
        max(info["dr"]) < 1e-2
        and max(info["dv"]) < 1e-6
        and max(abs(g) for g in info["gauge"]) < 1e-7
    )


class _DoneError(Exception):
    """Raised inside the trf residual once the floors are met."""


def solve(
    p: P, z0: Arr, newton_iters: int, trf_nfev: int, accept_noise: bool
) -> tuple[Arr, dict[str, Any]]:
    """#1004 solve: damped Newton (column-scaled lstsq, step cap 0.5 r_p, backtracking), then
    trf."""

    def ok(info: dict[str, Any]) -> bool:
        return converged(info) or (accept_noise and noise_ok(info))

    bad: dict[str, Any] = {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf]}
    z = np.asarray(z0, dtype=np.float64).copy()
    info = bad
    for _ in range(newton_iters):
        try:
            r, j, info = residual_and_jac(p, z)
        except RuntimeError:
            return z, bad
        if ok(info):
            return z, info
        d = np.linalg.norm(j, axis=0)
        d[d == 0.0] = 1.0
        step = np.linalg.lstsq(j / d, -r, rcond=None)[0] / d
        cap = 1.0
        for k in range(p.m):
            rr = float(np.linalg.norm(z[7 * k : 7 * k + 3]))
            dd = float(np.linalg.norm(step[7 * k : 7 * k + 3]))
            if dd > 0.5 * rr:
                cap = min(cap, 0.5 * rr / dd)
        step = cap * step
        f0 = float(np.linalg.norm(r))
        for alpha in [0.5**i for i in range(8)]:
            rn = residual_and_jac(p, z + alpha * step, False)[0]
            if np.all(np.isfinite(rn)) and float(np.linalg.norm(rn)) < f0:
                z = z + alpha * step
                break
        else:
            break
    if trf_nfev > 0:
        cache: dict[str, Any] = {}

        def fun(y: Arr) -> Arr:
            try:
                r, j, inf = residual_and_jac(p, y)
            except RuntimeError:
                return np.full(7 * p.m, 1e12)
            cache["z"], cache["j"] = y.copy(), j
            if ok(inf):
                cache["done"] = y.copy()
                raise _DoneError
            return r

        def jac(y: Arr) -> Arr:
            if "z" in cache and np.array_equal(cache["z"], y):
                return np.asarray(cache["j"])
            return residual_and_jac(p, y)[1]

        try:
            sol = least_squares(
                fun,
                z,
                jac=jac,
                method="trf",
                x_scale="jac",
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
                max_nfev=trf_nfev,
            )
            z = np.asarray(sol.x)
        except _DoneError:
            z = np.asarray(cache["done"])
        except RuntimeError:
            return z, bad
    try:
        _, _, info = residual_and_jac(p, z, False)
    except RuntimeError:
        info = bad
    return z, info


def describe(p: P, z: Arr, sigma: float) -> list[dict[str, Any]]:
    nd = nodes(p, z)
    out = []
    for k in range(p.m):
        moon = p.moons[k]
        x, _, t, _ = nd[k]
        rm, vm = p.circ.state(moon, t)
        d = float(np.linalg.norm(x[:3] - rm))
        mu = p.mus[moon]
        vinf2 = float(np.linalg.norm(x[3:] - vm)) ** 2 - 2.0 * mu / d
        e = 1.0 + d * vinf2 / mu
        out.append(
            {
                "moon": moon,
                "t_days": t / DAY,
                "vinf_kms": math.sqrt(max(vinf2, 0.0)),
                "rp_km": d,
                "rp_over_sigma": d / sigma,
                "turn_deg": math.degrees(2.0 * math.asin(1.0 / e)) if e > 1.0 else float("nan"),
                "impact": d <= sigma * SATELLITES[moon].radius_eq_km,
                "floor_ok": d - sigma * SATELLITES[moon].radius_eq_km >= sigma * FLOOR_KM[moon],
            }
        )
    return out


def scale_offsets(p: P, z: Arr, f: float) -> Arr:
    z = z.copy()
    for k in range(p.m):
        z[7 * k : 7 * k + 3] *= f
    return z


def _js(p: P, z: Arr, s: float) -> Arr:
    out = []
    for f in (1.0 + 1e-6, 1.0 - 1e-6):
        p.set_sigma(s * f)
        out.append(residual_and_jac(p, z, False)[0])
    p.set_sigma(s)
    return (out[0] - out[1]) / (2e-6 * s)


def _point(p: P, z: Arr, sg: float, info: dict[str, Any]) -> dict[str, Any]:
    return {
        "sigma": sg,
        "z": z.tolist(),
        "nodes": describe(p, z, sg),
        "floor": bool(converged(info)),
        "dr": max(info["dr"]),
    }


def monitor_step(p: P, rec: dict[str, Any]) -> bool:
    """One monitor-variable step (#1034 sec. 2): sigma is an unknown and the state component that
    changed most over the last secant (column-scaled) is stepped and fixed by an extra row. No SVD
    tangent is ever used (the enforced tangent rule). The step grows x1.5 after a fast convergence
    (at most 4 secants) and halves on failure (at most six halvings)."""
    pts = rec["upper"]
    a, b = pts[-2], pts[-1]
    za, zb = np.asarray(a["z"]), np.asarray(b["z"])
    p.set_sigma(b["sigma"])
    _, jz, _ = residual_and_jac(p, zb)
    d = np.linalg.norm(jz, axis=0)
    d[d == 0.0] = 1.0
    j = int(np.argmax(np.abs((zb - za) * d)))
    ratio = float(rec.get("mon_ratio", 1.0))
    for _half in range(7):
        target = zb[j] + ratio * (zb[j] - za[j])
        y = np.concatenate(
            [zb + ratio * (zb - za), [b["sigma"] + ratio * (b["sigma"] - a["sigma"])]]
        )
        for it in range(15):
            sg = float(y[-1])
            p.set_sigma(sg)
            try:
                r, jzz, info = residual_and_jac(p, y[:-1])
            except RuntimeError:
                break
            if (converged(info) or noise_ok(info)) and abs(y[j] - target) < 1e-9 * max(
                1.0, abs(target)
            ):
                pts.append(_point(p, y[:-1], sg, info))
                rec["mon_ratio"] = min(4.0, 1.5 * ratio) if it <= 4 else ratio
                _log(f"{rec['member']}: monitor sigma={sg:.6f} ratio {ratio:.3g}")
                return True
            js = _js(p, y[:-1], sg)
            row = np.zeros(len(y))
            row[j] = 1.0
            jac = np.vstack([np.hstack([jzz, js[:, None]]), row[None, :]])
            f = np.concatenate([r, [y[j] - target]])
            if not (np.all(np.isfinite(jac)) and np.all(np.isfinite(f))):
                break  # a sigma-derivative probe hit a failed propagation: a failed attempt
            dc = np.linalg.norm(jac, axis=0)
            dc[dc == 0.0] = 1.0
            try:
                step = np.linalg.lstsq(jac / dc, -f, rcond=None)[0] / dc
            except np.linalg.LinAlgError:
                break
            f0 = float(np.linalg.norm(f))
            for alpha in [0.5**i for i in range(8)]:
                yt = y + alpha * step
                p.set_sigma(float(yt[-1]))
                rt = residual_and_jac(p, yt[:-1], False)[0]
                if (
                    np.all(np.isfinite(rt))
                    and float(np.linalg.norm(np.concatenate([rt, [yt[j] - target]]))) < f0
                ):
                    y = yt
                    break
            else:
                break
        ratio *= 0.5
    rec["monitor_stop"] = "six halvings without convergence"
    return False


def two_solutions(p: P, rec: dict[str, Any]) -> None:
    """#1004 amendment-1 check 1, generalised: at 0.915 and 0.976 sigma_f (gc-2: 0.150, 0.160),
    fixed-sigma Newton from the nearest lower- and upper-branch points; supported if both converge
    and some node differs by > 0.02 km/s in V_inf and > 100 km in r_p / sigma."""
    lower = rec["points"]
    upper = rec["upper"][2:]
    out: dict[str, Any] = {}
    for sg in (round(0.915 * rec["sigma_f"], 5), round(0.976 * rec["sigma_f"], 5)):
        res: dict[str, Any] = {}
        for lab, pts in (("lower", lower), ("upper", upper)):
            q = min(pts, key=lambda t: abs(t["sigma"] - sg))
            p.set_sigma(sg)
            z, info = solve(p, scale_offsets(p, np.asarray(q["z"]), sg / q["sigma"]), 25, 30, False)
            res[lab] = describe(p, z, sg) if converged(info) else None
        row: dict[str, Any] = {
            "lower": res["lower"],
            "upper": res["upper"],
            "fold_supported": False,
        }
        if res["lower"] and res["upper"]:
            dif = [
                (abs(x["vinf_kms"] - y["vinf_kms"]), abs(x["rp_over_sigma"] - y["rp_over_sigma"]))
                for x, y in zip(res["lower"], res["upper"], strict=True)
            ]
            row["d_vinf_d_rp_over_sigma"] = dif
            row["fold_supported"] = any(dv > 0.02 and drp > 100.0 for dv, drp in dif)
        out[str(sg)] = row
        _log(f"{rec['member']}: two solutions at {sg}: supported={row['fold_supported']}")
    rec["two_solutions"] = out


def direct_sigma1(p: P, rec: dict[str, Any]) -> None:
    """#1034 sec. 3 check 3: damped Newton at sigma = 1 from the patched-conic seed and from each
    branch's last point (offsets scaled). A converged orbit is recorded, not identified with the
    member."""
    p.set_sigma(1.0)
    out = {}
    starts = {
        "seed": (p.seed, 1.0),
        "lower": (np.asarray(rec["points"][-1]["z"]), rec["points"][-1]["sigma"]),
    }
    if rec.get("upper"):
        starts["upper"] = (np.asarray(rec["upper"][-1]["z"]), rec["upper"][-1]["sigma"])
    for lab, (z0, s0) in starts.items():
        z, info = solve(p, scale_offsets(p, z0, 1.0 / s0), 25, 0, False)
        out[lab] = {
            "converged": bool(converged(info)),
            "dr": max(info["dr"]),
            "nodes": describe(p, z, 1.0) if converged(info) else None,
        }
    rec["direct_sigma1"] = out


def ias15_scaled(p: P, x0: Arr, t0: float, t1: float) -> Arr:
    """REBOUND IAS15 with the scaled moon GMs and softening radii (#1004 verify check 3)."""
    import rebound

    sim = rebound.Simulation()
    sim.G = MU_JUPITER_KM3_S2
    sim.integrator = "ias15"
    sim.integrator.epsilon = 1e-11
    sim.add(m=1.0)
    sim.add(m=0.0, x=x0[0], y=x0[1], z=x0[2], vx=x0[3], vy=x0[4], vz=x0[5])
    sim.t = float(t0)
    mus, surf = dict(p.mus), dict(p.surf)

    def forces(ptr: Any) -> None:
        s = ptr.contents
        sc = s.particles[1]
        r = np.array([sc.x, sc.y, sc.z])
        acc = np.zeros(3)
        for m in MOONS:
            rm = p.circ.state(m, float(s.t))[0]
            dd = rm - r
            dn = max(float(np.linalg.norm(dd)), surf[m])
            acc += mus[m] * (dd / dn**3 - rm / float(np.linalg.norm(rm)) ** 3)
        sc.ax += acc[0]
        sc.ay += acc[1]
        sc.az += acc[2]

    sim.additional_forces = forces
    sim.force_is_velocity_dependent = 0
    sim.integrate(float(t1))
    sc = sim.particles[1]
    return np.array([sc.x, sc.y, sc.z, sc.vx, sc.vy, sc.vz])


def ias15(p: P, z: Arr, rec: dict[str, Any], deadline: float) -> bool:
    """Every half-arc at the end point; resumable per half-arc. Pass: < 1e-2 km and < 1e-7 km/s."""
    nd = nodes(p, z)
    arcs = rec.setdefault("ias15_arcs", [])
    for k in range(len(arcs) // 2, p.m):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for x0, t0 in ((xp, tp), (xn, tn)):
            if time.monotonic() > deadline:
                return False
            xd, _ = arc(p, x0, t0, tm)
            xi = ias15_scaled(p, x0, t0, tm)
            arcs.append(
                [float(np.linalg.norm(xi[:3] - xd[:3])), float(np.linalg.norm(xi[3:] - xd[3:]))]
            )
    wr, wv = max(a[0] for a in arcs), max(a[1] for a in arcs)
    rec["ias15"] = {
        "max_dr_km": wr,
        "max_dv_kms": wv,
        "pass": wr < 1e-2 and wv < 1e-7,
        "mu_jupiter": MU_JUPITER_KM3_S2,
        "mus": p.mus,
    }
    return True


def run_member(member: str, deadline: float) -> dict[str, Any]:
    path = OUT / f"{member}.json"
    rec: dict[str, Any] = (
        json.loads(path.read_text())
        if path.exists()
        else {"member": member, "points": [], "stage": "start"}
    )
    if rec.get("stage") == "done":
        return rec
    p = build(member)

    def save() -> None:
        path.write_text(json.dumps(rec))

    def finish(outcome: str, **kw: Any) -> None:
        rec.update(stage="done", outcome=outcome, **kw)
        _log(f"{member}: {outcome} {kw}")

    while time.monotonic() < deadline and rec["stage"] != "done":
        st = rec["stage"]
        if st == "start":
            sg = 0.02
            p.set_sigma(sg)
            z, info = solve(p, scale_offsets(p, p.seed, sg), 40, 60, True)
            if not (converged(info) or noise_ok(info)):
                finish("NUMERICAL STOP", reason="no convergence at sigma 0.02")
            else:
                rec["points"] = [_point(p, z, sg, info)]
                rec["stage"] = "identity"
                _log(f"{member}: sigma=0.02 converged")
        elif st == "identity":
            z, s_prev = np.asarray(rec["points"][0]["z"]), 0.02
            rows = []
            for sg in (0.01, 0.005, 0.0025):
                p.set_sigma(sg)
                zn, info = solve(p, scale_offsets(p, z, sg / s_prev), 40, 60, True)
                row: dict[str, Any] = {
                    "sigma": sg,
                    "ok": bool(converged(info) or noise_ok(info)),
                    "dr": max(info["dr"]),
                }
                if row["ok"]:
                    d = describe(p, zn, sg)
                    row["shift"] = max(
                        abs(n["vinf_kms"] - v) for n, v in zip(d, p.pc_vinf, strict=True)
                    )
                    z, s_prev = zn, sg
                rows.append(row)
                if not row["ok"]:
                    break
            sh = [
                max(
                    abs(n["vinf_kms"] - v)
                    for n, v in zip(rec["points"][0]["nodes"], p.pc_vinf, strict=True)
                )
            ]
            sh += [r["shift"] for r in rows if r["ok"]]
            small = [r for r in rows if r["ok"] and r["sigma"] <= 0.005]
            mono = all(b <= a + 1e-4 for a, b in itertools.pairwise(sh))
            rec["identity"] = {"rows": rows, "shift_002": sh[0], "monotone": mono}
            rec["identity_pass"] = bool(small) and small[-1]["shift"] < 0.005 and mono
            _log(
                f"{member}: identity shifts {[round(x, 4) for x in sh]} pass={rec['identity_pass']}"
            )
            if rec["identity_pass"]:
                rec["stage"], rec["fac"] = "up", 1.2
            else:
                finish("IDENTITY FAIL")
        elif st == "up":
            pts = rec["points"]
            s1, z1 = pts[-1]["sigma"], np.asarray(pts[-1]["z"])
            if s1 >= 1.0:
                rec["stage"] = "ias15"
                save()
                continue
            fac = rec["fac"]
            sg = min(1.0, s1 * fac)
            if sg - s1 < 1e-5:
                rec["sigma_f"] = s1
                rec["stage"] = "fold"
                rec["upper"] = [pts[-2], pts[-1]]
                _log(f"{member}: natural step below 1e-5 at sigma {s1:.6f}; monitor")
                save()
                continue
            if len(pts) >= 2:
                s0, z0 = pts[-2]["sigma"], np.asarray(pts[-2]["z"])
                zs = z1 + (z1 - z0) * math.log(sg / s1) / math.log(s1 / s0)
            else:
                zs = scale_offsets(p, z1, sg / s1)
            p.set_sigma(sg)
            z, info = solve(p, zs, 25, 15, sg < 1.0)
            if converged(info) or (sg < 1.0 and noise_ok(info)):
                q = _point(p, z, sg, info)
                if any(n["impact"] for n in q["nodes"]):
                    finish("IMPACT", sigma_i=sg, impact_point=q)
                else:
                    pts.append(q)
                    rec["fac"] = min(1.2, 1.0 + 2.0 * (fac - 1.0))
                    vs = [round(n["vinf_kms"], 4) for n in q["nodes"]]
                    _log(f"{member}: sigma={sg:.6g} ok vinf {vs}")
            else:
                rec["fac"] = 1.0 + 0.5 * (fac - 1.0)
                _log(f"{member}: sigma={sg:.6g} failed; factor {rec['fac']:.6g}")
        elif st == "fold":
            up = rec["upper"]
            sf = rec["sigma_f"]
            if up[-1]["sigma"] > sf + 1e-4:
                # The monitor passed the stop with sigma still rising: not a turning point. Resume.
                rec["points"].append(up[-1])
                rec.setdefault("passed_stops", []).append(sf)
                rec.pop("upper")
                rec.pop("mon_ratio", None)
                rec["stage"], rec["fac"] = "up", 1.01
                _log(f"{member}: monitor passed the stop at {sf:.6f}, sigma rising; resume")
            elif up[-1]["sigma"] < 0.91 * sf:
                rec["stage"] = "fold_checks"
            elif not monitor_step(p, rec):
                finish(
                    "NUMERICAL STOP",
                    reason=f"natural stop at {sf:.6g}; monitor: {rec['monitor_stop']}",
                )
        elif st == "fold_checks":
            two_solutions(p, rec)
            direct_sigma1(p, rec)
            sup = all(v["fold_supported"] for v in rec["two_solutions"].values())
            if sup:
                finish("FOLDS", sigma_f=rec["sigma_f"])
            else:
                finish("NUMERICAL STOP", reason="turn passed but two-solution check not supported")
        elif st == "ias15":
            q = rec["points"][-1]
            p.set_sigma(1.0)
            if ias15(p, np.asarray(q["z"]), rec, deadline):
                rec["full_mass"] = q["nodes"]
                finish("EXISTS", ias15_pass=rec["ias15"]["pass"])
        save()
    save()
    return rec


def summary() -> None:
    rows = []
    for mid, _key, _vg, _vc in MEMBERS:
        path = OUT / f"{mid}.json"
        if not path.exists():
            rows.append({"member": mid, "outcome": "not run"})
            continue
        r = json.loads(path.read_text())
        row: dict[str, Any] = {
            "member": mid,
            "outcome": r.get("outcome", f"in progress ({r.get('stage')})"),
            "identity_pass": r.get("identity_pass"),
            "sigma_f": r.get("sigma_f"),
            "sigma_i": r.get("sigma_i"),
            "last_sigma": r["points"][-1]["sigma"] if r.get("points") else None,
            "n_points": len(r.get("points", [])),
            "passed_stops": r.get("passed_stops"),
            "ias15_pass": (r.get("ias15") or {}).get("pass"),
            "floor_ok_all": all(n["floor_ok"] for q in r.get("points", []) for n in q["nodes"])
            if r.get("points")
            else None,
        }
        if r.get("full_mass"):
            row["full_mass_vinf"] = [round(n["vinf_kms"], 4) for n in r["full_mass"]]
            row["full_mass_alt_km"] = [
                round(n["rp_km"] - SATELLITES[n["moon"]].radius_eq_km, 1) for n in r["full_mass"]
            ]
        rows.append(row)
    (OUT / "summary.json").write_text(json.dumps(rows, indent=1))
    for row in rows:
        print(" | ".join(f"{k}={v}" for k, v in row.items() if v is not None))


def main() -> None:
    preflight_search(
        task_no=1025,
        region_id="gc-k4-k6-clean-members-ideal-continuous",
        method=MethodCapability(
            genome="clean #973 gc members (k = 4-6), R-S circular coplanar, both moons massive",
            corrector="forward-backward shooting, DOP853 + STM; joint sigma; fold checks",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=len(MEMBERS),
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", default="batch", help="'batch', or comma list of ids")
    ap.add_argument("--worker", default="0/1")
    ap.add_argument("--max-wall", type=float, default=420.0)
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.summary:
        summary()
        return
    ids = (
        [m[0] for m in MEMBERS if m[0] not in CONTROLS]
        if args.members == "batch"
        else args.members.split(",")
    )
    w, nw = (int(v) for v in args.worker.split("/"))
    ids = [mid for i, mid in enumerate(ids) if i % nw == w]
    deadline = time.monotonic() + args.max_wall
    for mid in ids:
        if time.monotonic() > deadline:
            break
        rec = run_member(mid, deadline)
        _log(f"{mid}: stage {rec.get('stage')} outcome {rec.get('outcome')}")


if __name__ == "__main__":
    main()
