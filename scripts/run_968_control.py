"""#968: Jovian n-body positive control, GanCal#5 in the R-S circular model made continuous.

Pre-registration: docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 3 (committed
before this script was first run). Corrector A (this file): forward-backward multiple shooting
with mid-leg match points, Ganymede flyby nodes at periapsis (gauge), a pinned massless-Callisto
hit and a rotation-aware periodic wrap; arcs and Jacobian from ``jovian_stm.propagate_with_stm``
(DOP853 + analytic STM). Corrector B: the lane's own ``jovian_defect_residual`` with the fixed
wrap. Cross-check: REBOUND IAS15 re-fly of every half-arc.

Usage::

    uv run python scripts/run_968_control.py --stage a      # corrector A, checkpointed
    uv run python scripts/run_968_control.py --stage check  # criteria + IAS15 cross-check
    uv run python scripts/run_968_control.py --stage b      # the lane, fixed wrap

Each call stays under about 8 minutes; the LM loop resumes from ``data/968_control/a_state.json``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import JovianRestrictedNBody, periapsis_node
from cyclerfinder.nbody.jovian_stm import _accel_and_gradient, propagate_with_stm
from cyclerfinder.search.two_working_body import CircularSystem, cycle_flybys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "968_control"
DAY = 86400.0
KEY = "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s"
W = np.array([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
W_GAUGE = 1e4
W_JACOBI = 1e3
RTOL, ATOL = 1e-12, 1e-10
R_GAN_RS_KM = 2634.0  # R-S 2009 Table 2 radius, for altitudes compared with the paper
GAN = "Ganymede"
CAL = "Callisto"
FORCE_MOONS: tuple[str, ...] = (GAN,)  # Callisto massless, as in R-S

Arr = NDArray[np.float64]


def _load_enum() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_942_enumerate", ROOT / "scripts" / "run_942_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _log(msg: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


def rot_z(a: float) -> Arr:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


@dataclass
class Problem:
    circ: CircularSystem
    period_s: float
    q6: Arr  # 6x6 blockdiag rotation by Ganymede's advance over the period
    seed: Arr
    mus: dict[str, float]
    surf: dict[str, float]


def build_problem() -> Problem:
    enum = _load_enum()
    circ, a, b = enum.cell_system("gc1")  # R-S Table 2, Callisto massless
    cands = {
        c["key"]: c
        for c in json.loads((ROOT / "data/943_cell_gc_gauntlet.json").read_text())["candidates"]
    }
    _, cyc = enum.parse_cycle_key(KEY, circ, a, b)
    x = np.asarray(cands[KEY]["x_days"]) * DAY
    fl = cycle_flybys(circ, cyc, x)
    assert fl is not None and [f.body for f in fl] == [GAN, CAL, GAN]
    period = cyc.period_s
    adv = 2.0 * math.pi * period / circ.period_s(GAN)
    q = rot_z(adv)
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    fb, fc, fa = fl
    rb, vb, _ = periapsis_node(GAN, fb.t_s, fb.vinf_in, fb.vinf_out, circ)  # type: ignore[arg-type]
    ra, va, _ = periapsis_node(GAN, fa.t_s, fa.vinf_in, fa.vinf_out, circ)  # type: ignore[arg-type]
    _, vcal = circ.state(CAL, fc.t_s)
    seed = np.concatenate([rb, vb, [fb.t_s], vcal + fc.vinf_in, [fc.t_s], ra, va, [fa.t_s]]).astype(
        np.float64
    )
    mus = {m: SATELLITES[m].mu_km3_s2 for m in FORCE_MOONS}
    surf = {m: SATELLITES[m].radius_eq_km for m in FORCE_MOONS}
    return Problem(circ, period, q6, seed, mus, surf)


# Layouts. "pin" (pre-registered, 18 vars): B state 0:6, tB 6, C velocity 7:10, tC 10,
# A state 11:17, tA 17; Callisto hit pinned, T fixed. "free" (diagnostic, 20 vars): B state
# 0:6 (tB fixed), C state 6:12 (tC fixed), A state 12:18, tA 18, T 19; no Callisto pin, plus a
# Jacobi-constant target at B. Nodes in time order B -> C -> A -> B' = (Q(T) x_B, tB + T).

K_ROT = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])


@dataclass
class Layout:
    mode: str  # "pin" | "free"
    tb: float = 0.0  # fixed epochs for "free"
    tc: float = 0.0
    j_target: float = 0.0
    jacobi_row: bool = True
    t_target_s: float | None = None  # "free": a period row instead of the Jacobi row

    @property
    def n(self) -> int:
        return 18 if self.mode == "pin" else 20


def q_of(p: Problem, period: float) -> tuple[Arr, Arr]:
    n_g = 2.0 * math.pi / p.circ.period_s(GAN)
    q = rot_z(n_g * period)
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    dq6 = np.zeros((6, 6))
    dq6[:3, :3] = n_g * q @ K_ROT
    dq6[3:, 3:] = n_g * q @ K_ROT
    return q6, dq6


def nodes(p: Problem, z: Arr, lay: Layout | None = None) -> list[tuple[Arr, Arr, float, Arr]]:
    """Each node: (state, d state/dz, epoch, d epoch/dz)."""
    lay = lay or Layout("pin")
    n = lay.n
    out = []
    if lay.mode == "pin":
        dxb = np.zeros((6, n))
        dxb[:, 0:6] = np.eye(6)
        dtb = np.zeros(n)
        dtb[6] = 1.0
        out.append((z[0:6].copy(), dxb, float(z[6]), dtb))
        tc = float(z[10])
        rc, vc = p.circ.state(CAL, tc)
        dxc = np.zeros((6, n))
        dxc[3:6, 7:10] = np.eye(3)
        dxc[0:3, 10] = vc  # the pinned position moves with Callisto
        dtc = np.zeros(n)
        dtc[10] = 1.0
        out.append((np.concatenate([rc, z[7:10]]), dxc, tc, dtc))
        dxa = np.zeros((6, n))
        dxa[:, 11:17] = np.eye(6)
        dta = np.zeros(n)
        dta[17] = 1.0
        out.append((z[11:17].copy(), dxa, float(z[17]), dta))
        out.append((p.q6 @ z[0:6], p.q6 @ dxb, float(z[6]) + p.period_s, dtb.copy()))
        return out
    dxb = np.zeros((6, n))
    dxb[:, 0:6] = np.eye(6)
    out.append((z[0:6].copy(), dxb, lay.tb, np.zeros(n)))
    dxc = np.zeros((6, n))
    dxc[:, 6:12] = np.eye(6)
    out.append((z[6:12].copy(), dxc, lay.tc, np.zeros(n)))
    dxa = np.zeros((6, n))
    dxa[:, 12:18] = np.eye(6)
    dta = np.zeros(n)
    dta[18] = 1.0
    out.append((z[12:18].copy(), dxa, float(z[18]), dta))
    period = float(z[19])
    q6, dq6 = q_of(p, period)
    dxw = q6 @ dxb
    dxw[:, 19] = dq6 @ z[0:6]
    dtw = np.zeros(n)
    dtw[19] = 1.0
    out.append((q6 @ z[0:6], dxw, lay.tb + period, dtw))
    return out


def accel(p: Problem, x: Arr, t: float) -> Arr:
    a, _ = _accel_and_gradient(x[:3], t, ephem=p.circ, moons=FORCE_MOONS, mus=p.mus, surf=p.surf)  # type: ignore[arg-type]
    return np.concatenate([x[3:], a])


def arc(p: Problem, x: Arr, t0: float, t1: float) -> tuple[Arr, Arr]:
    rf, vf, phi = propagate_with_stm(
        x[:3],
        x[3:],
        t0,
        t1,
        ephem=p.circ,  # type: ignore[arg-type]
        moons=FORCE_MOONS,
        rtol=RTOL,
        atol=ATOL,
        mu_overrides=p.mus,
        radius_overrides=p.surf,
    )
    return np.concatenate([rf, vf]), phi


def gauge(p: Problem, x: Arr, t: float) -> float:
    rg, vg = p.circ.state(GAN, t)
    dr, dv = x[:3] - rg, x[3:] - vg
    return float(dr @ dv) / (float(np.linalg.norm(dr)) * float(np.linalg.norm(dv)))


def gauge_grad(p: Problem, x: Arr, t: float) -> tuple[Arr, float]:
    """Analytic d g/d x (6) and d g/d t for the periapsis gauge."""
    rg, vg = p.circ.state(GAN, t)
    n_g = 2.0 * math.pi / p.circ.period_s(GAN)
    ag = -(n_g**2) * rg
    dr, dv = x[:3] - rg, x[3:] - vg
    a, b = float(np.linalg.norm(dr)), float(np.linalg.norm(dv))
    g = float(dr @ dv) / (a * b)
    g_dr = dv / (a * b) - g * dr / a**2
    g_dv = dr / (a * b) - g * dv / b**2
    return np.concatenate([g_dr, g_dv]), float(-(g_dr @ vg) - (g_dv @ ag))


def jacobi(p: Problem, x: Arr, t: float) -> float:
    """Jacobi-like integral (energy in Ganymede's rotating frame, lane force model)."""
    from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2

    n_g = 2.0 * math.pi / p.circ.period_s(GAN)
    w = np.array([0.0, 0.0, n_g])
    r, v = x[:3], x[3:]
    rg, _ = p.circ.state(GAN, t)
    u = v - np.cross(w, r)
    mu_g = p.mus[GAN]
    return float(
        0.5 * u @ u
        - MU_JUPITER_KM3_S2 / np.linalg.norm(r)
        - mu_g / np.linalg.norm(r - rg)
        + mu_g * float(r @ rg) / float(np.linalg.norm(rg)) ** 3
        - 0.5 * float(np.linalg.norm(np.cross(w, r))) ** 2
    )


def jacobi_grad(p: Problem, x: Arr, t: float) -> tuple[Arr, float]:
    hs = np.array([1e-2, 1e-2, 1e-2, 1e-7, 1e-7, 1e-7])
    gx = np.zeros(6)
    for i in range(6):
        e = np.zeros(6)
        e[i] = hs[i]
        gx[i] = (jacobi(p, x + e, t) - jacobi(p, x - e, t)) / (2 * hs[i])
    gt = (jacobi(p, x, t + 1.0) - jacobi(p, x, t - 1.0)) / 2.0
    return gx, gt


def residual_and_jac(
    p: Problem, z: Arr, *, want_jac: bool = True, lay: Layout | None = None
) -> tuple[Arr, Arr, dict[str, Any]]:
    lay = lay or Layout("pin")
    nd = nodes(p, z, lay)
    res: list[float] = []
    jac_rows: list[Arr] = []
    info: dict[str, Any] = {"match_dr_km": [], "match_dv_kms": [], "gauge": []}
    for k in range(3):
        xp, dxp, tp, dtp = nd[k]
        xn, dxn, tn, dtn = nd[k + 1]
        tm = 0.5 * (tp + tn)
        xf, pf = arc(p, xp, tp, tm)
        xb, pb = arc(p, xn, tn, tm)
        d = xf - xb
        info["match_dr_km"].append(float(np.linalg.norm(d[:3])))
        info["match_dv_kms"].append(float(np.linalg.norm(d[3:])))
        res.extend(W * d)
        if want_jac:
            ff = accel(p, xf, tm)
            fb = accel(p, xb, tm)
            dxf = (
                pf @ dxp
                + np.outer(-pf @ accel(p, xp, tp) + 0.5 * ff, dtp)
                + np.outer(0.5 * ff, dtn)
            )
            dxb_ = (
                pb @ dxn
                + np.outer(-pb @ accel(p, xn, tn) + 0.5 * fb, dtn)
                + np.outer(0.5 * fb, dtp)
            )
            jac_rows.append(W[:, None] * (dxf - dxb_))
    for idx in (0, 2):
        x, dx, t, dt = nd[idx]
        g = gauge(p, x, t)
        res.append(W_GAUGE * g)
        info["gauge"].append(g)
        if want_jac:
            gx, gt = gauge_grad(p, x, t)
            jac_rows.append(W_GAUGE * (gx @ dx + gt * dt)[None, :])
    if lay.mode == "free":
        info["jacobi"] = jacobi(p, nd[0][0], nd[0][2])
    if lay.mode == "free" and lay.jacobi_row:
        x, dx, t, dt = nd[0]
        jv = jacobi(p, x, t)
        res.append(W_JACOBI * (jv - lay.j_target))
        info["jacobi"] = jv
        if want_jac:
            gx, gt = jacobi_grad(p, x, t)
            jac_rows.append(W_JACOBI * (gx @ dx + gt * dt)[None, :])
    if lay.mode == "free" and lay.t_target_s is not None:
        res.append((float(z[19]) - lay.t_target_s) / 60.0)
        if want_jac:
            row = np.zeros(lay.n)
            row[19] = 1.0 / 60.0
            jac_rows.append(row[None, :])
    jac = np.vstack(jac_rows) if want_jac else np.zeros((0, lay.n))
    return np.asarray(res), jac, info


def check_jacobian(p: Problem, z: Arr, lay: Layout | None = None) -> dict[str, float]:
    """FD check of the analytic Jacobian, every column (instrument check)."""
    lay = lay or Layout("pin")
    _, jac, _ = residual_and_jac(p, z, lay=lay)
    out = {}
    for j in range(lay.n):
        h = 1e-6 * max(1.0, abs(float(z[j])))
        zp = z.copy()
        zp[j] += h
        zm = z.copy()
        zm[j] -= h
        fd = (
            residual_and_jac(p, zp, want_jac=False, lay=lay)[0]
            - residual_and_jac(p, zm, want_jac=False, lay=lay)[0]
        ) / (2 * h)
        out[f"col{j}"] = float(np.linalg.norm(fd - jac[:, j]) / max(np.linalg.norm(fd), 1e-30))
    return out


def solve_trf(
    p: Problem, z0: Arr, lay: Layout, max_nfev: int, tag: str
) -> tuple[Arr, dict[str, Any]]:
    """scipy trf with the analytic Jacobian; logs every evaluation."""
    from scipy.optimize import least_squares

    cache: dict[str, Any] = {}

    def fun(z: Arr) -> Arr:
        r, j, info = residual_and_jac(p, z, lay=lay)
        cache["z"], cache["j"] = z.copy(), j
        _log(
            f"{tag} |r| {np.linalg.norm(r):.4e} dr {max(info['match_dr_km']):.3e} km "
            f"dv {max(info['match_dv_kms']):.3e} km/s "
            f"gauge {max(abs(g) for g in info['gauge']):.2e}"
            + (f" J {info['jacobi']:.6f} T {z[19] / DAY:.6f} d" if lay.mode == "free" else "")
        )
        return r

    def jac(z: Arr) -> Arr:
        if "z" in cache and np.array_equal(cache["z"], z):
            return np.asarray(cache["j"])
        return residual_and_jac(p, z, lay=lay)[1]

    sol = least_squares(
        fun,
        z0,
        jac=jac,
        method="trf",
        x_scale="jac",
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=max_nfev,
    )
    r, _, info = residual_and_jac(p, sol.x, want_jac=False, lay=lay)
    info["norm"] = float(np.linalg.norm(r))
    info["nfev"] = int(sol.nfev)
    info["status"] = int(sol.status)
    return np.asarray(sol.x), info


def converged(info: dict[str, Any]) -> bool:
    return (
        max(info["match_dr_km"]) < 1e-3
        and max(info["match_dv_kms"]) < 1e-6
        and max(abs(g) for g in info["gauge"]) < 1e-9
    )


def stage_a(max_nfev: int) -> None:
    """Pre-registered corrector A (amendment 1: scipy trf replaces the hand LM loop)."""
    p = build_problem()
    lay = Layout("pin")
    state_path = OUT / "a_state.json"
    if state_path.exists():
        z = np.asarray(json.loads(state_path.read_text())["z"])
        _log("stage a: resume from a_state.json")
    else:
        z = p.seed.copy()
        jc = check_jacobian(p, z, lay)
        _log(f"stage a: Jacobian FD check {jc}")
        (OUT / "a_jacobian_check.json").write_text(json.dumps(jc, indent=1))
    z, info = solve_trf(p, z, lay, max_nfev, "stage a")
    state_path.write_text(json.dumps({"z": z.tolist(), "info": info}))
    _log(f"stage a: end {info['norm']:.4e} converged={converged(info)} nfev {info['nfev']}")


def stage_free(max_nfev: int, d_j: float, start: str) -> None:
    """Diagnostic (note amendment 1): T free, Callisto unpinned, Jacobi at B fixed."""
    p = build_problem()
    zp = np.asarray(json.loads((OUT / start).read_text())["z"])
    path = OUT / "free_state.json"
    if start == "free_state.json":
        st = json.loads(path.read_text())
        lay = Layout("free", st["tb"], st["tc"], st["j_target"] + d_j)
        z = np.asarray(st["z"])
    else:
        nd = nodes(p, zp)
        tb, tc = nd[0][2], nd[1][2]
        lay = Layout("free", tb, tc, jacobi(p, nd[0][0], tb) + d_j)
        z = np.concatenate([nd[0][0], nd[1][0], nd[2][0], [nd[2][2]], [p.period_s]])
        jc = check_jacobian(p, z, lay)
        _log(f"stage free: Jacobian FD check {jc}")
    z, info = solve_trf(p, z, lay, max_nfev, f"free dJ={d_j:+.4f}")
    path.write_text(
        json.dumps(
            {"z": z.tolist(), "tb": lay.tb, "tc": lay.tc, "j_target": lay.j_target, "info": info}
        )
    )
    _log(
        f"stage free: end {info['norm']:.4e} converged={converged(info)} T {z[19] / DAY:.6f} d "
        f"(3 S = {p.period_s / DAY:.6f} d)"
    )


def sample_orbit(p: Problem, z: Arr) -> dict[str, Any]:
    """Dense DOP853 over one period from B: Ganymede approaches, r range, Callisto miss."""
    x0 = z[0:6]
    t0 = float(z[6])

    def rhs(t: float, y: Arr) -> Arr:
        return accel(p, y, t)

    sol = solve_ivp(
        rhs, (t0, t0 + p.period_s), x0, method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True
    )
    ts = np.arange(t0, t0 + p.period_s, 0.01 * DAY)
    ys = sol.sol(ts)
    rr = np.linalg.norm(ys[:3], axis=0)
    dg = np.array(
        [np.linalg.norm(ys[:3, i] - p.circ.state(GAN, float(t))[0]) for i, t in enumerate(ts)]
    )
    hill = 31_715.0
    mins = []
    for i in range(1, len(ts) - 1):
        if dg[i] <= dg[i - 1] and dg[i] <= dg[i + 1] and dg[i] < hill:
            res = minimize_scalar(
                lambda t: float(np.linalg.norm(sol.sol(t)[:3] - p.circ.state(GAN, t)[0])),
                bounds=(ts[i - 1], ts[i + 1]),
                method="bounded",
                options={"xatol": 1e-3},
            )
            mins.append(
                {
                    "t_days": float(res.x) / DAY,
                    "dist_km": float(res.fun),
                    "alt_rs_km": float(res.fun) - R_GAN_RS_KM,
                }
            )
    # endpoints of the window (a flyby at the start/end node)
    for t in (t0, t0 + p.period_s):
        d0 = float(np.linalg.norm(sol.sol(t)[:3] - p.circ.state(GAN, t)[0]))
        if d0 < hill:
            mins.append(
                {
                    "t_days": t / DAY,
                    "dist_km": d0,
                    "alt_rs_km": d0 - R_GAN_RS_KM,
                    "window_edge": True,
                }
            )
    rmin = minimize_scalar(
        lambda t: float(np.linalg.norm(sol.sol(t)[:3])),
        bounds=(ts[np.argmin(rr)] - 0.02 * DAY, ts[np.argmin(rr)] + 0.02 * DAY),
        method="bounded",
    )
    rmax = minimize_scalar(
        lambda t: -float(np.linalg.norm(sol.sol(t)[:3])),
        bounds=(ts[np.argmax(rr)] - 0.02 * DAY, ts[np.argmax(rr)] + 0.02 * DAY),
        method="bounded",
    )
    xend = sol.sol(t0 + p.period_s)
    return {
        "ganymede_close_approaches": mins,
        "r_min_km": float(rmin.fun),
        "r_max_km": float(-rmax.fun),
        "single_shot_return_dr_km": float(np.linalg.norm(xend[:3] - (p.q6 @ x0)[:3])),
        "single_shot_return_dv_kms": float(np.linalg.norm(xend[3:] - (p.q6 @ x0)[3:])),
    }


def vinf_gan(p: Problem, x: Arr, t: float) -> tuple[float, float]:
    rg, vg = p.circ.state(GAN, t)
    d = float(np.linalg.norm(x[:3] - rg))
    v2 = float(np.linalg.norm(x[3:] - vg)) ** 2
    return math.sqrt(v2 - 2.0 * p.mus[GAN] / d), d


class _ExactCache:
    def __init__(self, circ: CircularSystem) -> None:
        self.circ = circ

    def position(self, moon: str, t_sec: float) -> Arr:
        return self.circ.state(moon, t_sec)[0]


def ias15(p: Problem, x: Arr, t0: float, t1: float) -> tuple[Arr | None, str]:
    prop = JovianRestrictedNBody()
    arc_ = prop.propagate(
        x[:3],
        x[3:],
        t0,
        t1,
        moons=FORCE_MOONS,
        cache=_ExactCache(p.circ),
        accuracy=1e-11,
        max_wall_sec=3000.0,  # type: ignore[arg-type]
    )
    if not arc_.converged:
        return None, "timeout_or_nonfinite"
    if abs(arc_.t1_sec - t1) > 1e-6:
        return None, f"stopped_at_{arc_.t1_sec}"
    return np.concatenate([arc_.r_km, arc_.v_km_s]), "ok"


def stage_check() -> None:
    p = build_problem()
    z = np.asarray(json.loads((OUT / "a_state.json").read_text())["z"])
    _, _, info = residual_and_jac(p, z, want_jac=False)
    nd = nodes(p, z)  # check reads the pinned layout
    out: dict[str, Any] = {
        "note": "docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 3.4"
    }
    out["match_dr_km"] = info["match_dr_km"]
    out["match_dv_kms"] = info["match_dv_kms"]
    out["gauge"] = info["gauge"]
    rc, vc = p.circ.state(CAL, nd[1][2])
    out["callisto_hit_km"] = float(np.linalg.norm(nd[1][0][:3] - rc))
    out["callisto_rel_speed_kms"] = float(np.linalg.norm(nd[1][0][3:] - vc))
    vb, db = vinf_gan(p, nd[0][0], nd[0][2])
    va, da = vinf_gan(p, nd[2][0], nd[2][2])
    out["ganymede"] = {
        "B": {"t_days": nd[0][2] / DAY, "vinf_kms": vb, "alt_rs_km": db - R_GAN_RS_KM},
        "A": {"t_days": nd[2][2] / DAY, "vinf_kms": va, "alt_rs_km": da - R_GAN_RS_KM},
    }
    out["epochs_days"] = [n[2] / DAY for n in nd]
    out.update(sample_orbit(p, z))
    # IAS15 cross-check of every half-arc
    xcheck = []
    for k in range(3):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for lab, x0, t0 in (("fwd", xp, tp), ("bwd", xn, tn)):
            xd, _ = arc(p, x0, t0, tm)
            t_w = time.monotonic()
            xi, status = ias15(p, x0, t0, tm)
            row: dict[str, Any] = {
                "leg": k,
                "dir": lab,
                "status": status,
                "wall_s": time.monotonic() - t_w,
            }
            if xi is not None:
                row["dr_km"] = float(np.linalg.norm(xi[:3] - xd[:3]))
                row["dv_kms"] = float(np.linalg.norm(xi[3:] - xd[3:]))
            xcheck.append(row)
            _log(f"IAS15 leg {k} {lab}: {row}")
    out["ias15_halfarc"] = xcheck
    xi, status = ias15(p, nd[0][0], nd[0][2], nd[0][2] + p.period_s)
    if xi is not None:
        tgt = p.q6 @ nd[0][0]
        out["ias15_full_period"] = {
            "dr_km": float(np.linalg.norm(xi[:3] - tgt[:3])),
            "dv_kms": float(np.linalg.norm(xi[3:] - tgt[3:])),
        }
    else:
        out["ias15_full_period"] = {"status": status}
    (OUT / "a_check.json").write_text(json.dumps(out, indent=1))
    _log(f"stage check: wrote a_check.json: {json.dumps(out)[:1500]}")


def family_point(p: Problem, z: Arr, lay: Layout) -> dict[str, Any]:
    nd = nodes(p, z, lay)
    vb, db = vinf_gan(p, nd[0][0], nd[0][2])
    va, da = vinf_gan(p, nd[2][0], nd[2][2])
    return {
        "T_days": float(z[19]) / DAY,
        "jacobi": jacobi(p, nd[0][0], nd[0][2]),
        "vinf_B": vb,
        "alt_B_rs_km": db - R_GAN_RS_KM,
        "vinf_A": va,
        "alt_A_rs_km": da - R_GAN_RS_KM,
        "tA_days": float(z[18]) / DAY,
    }


def stage_family(n_steps: int, ds0: float, direction: int) -> None:
    """Diagnostic pseudo-arclength continuation of the T-free family (note amendment 2)."""
    from scipy.optimize import least_squares

    p = build_problem()
    path = OUT / f"family_{'up' if direction > 0 else 'down'}.json"
    st0 = json.loads((OUT / "free_state_J0.json").read_text())
    lay = Layout("free", st0["tb"], st0["tc"], 0.0, jacobi_row=False)
    if path.exists():
        fam = json.loads(path.read_text())
    else:
        z0 = np.asarray(st0["z"])
        _, j0, _ = residual_and_jac(p, z0, lay=lay)
        d = np.linalg.norm(j0, axis=0)
        d[d == 0.0] = 1.0
        fam = {
            "d": d.tolist(),
            "ds": ds0,
            "points": [{"z": z0.tolist(), **family_point(p, z0, lay)}],
        }
    d = np.asarray(fam["d"])
    t_prev: Arr | None = np.asarray(fam["tangent"]) if "tangent" in fam else None
    for _ in range(n_steps):
        z = np.asarray(fam["points"][-1]["z"])
        _, jz, _ = residual_and_jac(p, z, lay=lay)
        _, sv, vt = np.linalg.svd(jz / d, full_matrices=True)
        tan = vt[-1]
        if t_prev is None:
            tan = tan * (direction * np.sign(tan[19]) or 1.0)
        elif float(tan @ t_prev) < 0.0:
            tan = -tan
        ds = float(fam["ds"])
        u0 = z * d
        ok = False
        for _attempt in range(5):
            u_pred = u0 + ds * tan

            def fun(u: Arr, u0: Arr = u0, tan: Arr = tan, ds: float = ds) -> Arr:
                r, _, _ = residual_and_jac(p, u / d, want_jac=False, lay=lay)
                return np.concatenate([r, [1e3 * (float(tan @ (u - u0)) - ds)]])

            def jac(u: Arr, tan: Arr = tan) -> Arr:
                _, j, _ = residual_and_jac(p, u / d, lay=lay)
                return np.vstack([j / d, 1e3 * tan[None, :]])

            sol = least_squares(
                fun, u_pred, jac=jac, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=25
            )
            zn = np.asarray(sol.x) / d
            _, _, info = residual_and_jac(p, zn, want_jac=False, lay=lay)
            if converged(info):
                ok = True
                break
            ds *= 0.5
            _log(f"family: step failed (|r| dv {max(info['match_dv_kms']):.2e}); ds -> {ds:.3e}")
        if not ok:
            _log("family: step size underflow; stop")
            break
        pt = {"z": zn.tolist(), **family_point(p, zn, lay), "sv_min3": sv[-3:].tolist()}
        fam["points"].append(pt)
        fam["ds"] = min(ds * 1.5, 50.0 * ds0)
        fam["tangent"] = tan.tolist()
        t_prev = tan
        path.write_text(json.dumps(fam))
        _log(
            f"family {len(fam['points']) - 1}: T {pt['T_days']:.6f} d J {pt['jacobi']:.6f} "
            f"vinf B/A {pt['vinf_B']:.4f}/{pt['vinf_A']:.4f} alt B/A {pt['alt_B_rs_km']:.1f}/"
            f"{pt['alt_A_rs_km']:.1f} sv {sv[-3:]}"
        )
        if pt["T_days"] * DAY > p.period_s and direction > 0:
            _log("family: T crossed 3 S_GC; stop (bracket found)")
            break


def stage_tcont(t_list_days: list[float], max_nfev: int) -> None:
    """Diagnostic natural-parameter continuation in T (note amendment 2): J free."""
    p = build_problem()
    path = OUT / "tcont.json"
    st0 = json.loads((OUT / "free_state_J0.json").read_text())
    rec = json.loads(path.read_text()) if path.exists() else {"points": []}
    good = [q for q in rec["points"] if q["converged"]]
    z = np.asarray(good[-1]["z"]) if good else np.asarray(st0["z"])
    for t_d in t_list_days:
        lay = Layout("free", st0["tb"], st0["tc"], 0.0, jacobi_row=False, t_target_s=t_d * DAY)
        z_try = z.copy()
        z_try[19] = t_d * DAY
        zn, info = solve_trf(p, z_try, lay, max_nfev, f"tcont T={t_d:.4f}")
        ok = converged(info)
        pt = {"T_target_days": t_d, "converged": ok, "norm": info["norm"], "z": zn.tolist()}
        if ok:
            pt.update(family_point(p, zn, lay))
            z = zn
        rec["points"].append(pt)
        path.write_text(json.dumps(rec))
        _log(
            f"tcont T={t_d:.4f}: converged={ok} |r| {info['norm']:.3e} "
            + (
                f"J {pt['jacobi']:.6f} vinf B/A {pt['vinf_B']:.4f}/{pt['vinf_A']:.4f} "
                f"alt B/A {pt['alt_B_rs_km']:.1f}/{pt['alt_A_rs_km']:.1f}"
                if ok
                else ""
            )
        )
        if not ok:
            break


def stage_pinseed() -> None:
    """Amendment 3: pinned-layout seed from the T = 3 S_GC family member, time-shifted so that
    Callisto sits at the orbit's inbound crossing of Callisto's radius nearest the old C node."""
    from scipy.optimize import brentq

    p = build_problem()
    st0 = json.loads((OUT / "free_state_J0.json").read_text())
    rec = json.loads((OUT / "tcont.json").read_text())
    last = [q for q in rec["points"] if q["converged"]][-1]
    z = np.asarray(last["z"])
    lay = Layout("free", st0["tb"], st0["tc"], 0.0, jacobi_row=False)
    nd = nodes(p, z, lay)
    xc, tc = nd[1][0], nd[1][2]
    a_c = p.circ.bodies[CAL][0]

    def rhs(t: float, y: Arr) -> Arr:
        return accel(p, y, t)

    sol = solve_ivp(
        rhs, (tc, tc - 3 * DAY), xc, method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True
    )
    sol2 = solve_ivp(
        rhs, (tc, tc + 2 * DAY), xc, method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True
    )

    def dens(t: float) -> Arr:
        return np.asarray(sol.sol(t) if t <= tc else sol2.sol(t))

    ts = np.linspace(tc - 3 * DAY, tc + 2 * DAY, 2001)
    f = [float(np.linalg.norm(dens(t)[:3])) - a_c for t in ts]
    roots = []
    for i in range(len(ts) - 1):
        if f[i] * f[i + 1] < 0.0:
            tr = brentq(
                lambda t: float(np.linalg.norm(dens(t)[:3])) - a_c, ts[i], ts[i + 1], xtol=1e-6
            )
            y = dens(tr)
            if float(y[:3] @ y[3:]) < 0.0:  # inbound
                roots.append(tr)
    t_p = min(roots, key=lambda t: abs(t - tc))
    xp = dens(t_p)
    n_g = 2.0 * math.pi / p.circ.period_s(GAN)
    n_c = 2.0 * math.pi / p.circ.period_s(CAL)
    syn = 2.0 * math.pi / abs(n_c - n_g)
    th_c = 2.0 * math.pi * t_p / p.circ.period_s(CAL) + p.circ.bodies[CAL][2]
    phi = math.atan2(float(xp[1]), float(xp[0]))
    tau = ((phi - th_c) / (n_c - n_g)) % syn
    if tau > 0.5 * syn:
        tau -= syn
    q, _ = q_of(p, tau)
    zp = np.concatenate(
        [q @ nd[0][0], [nd[0][2] + tau], (q @ xp)[3:], [t_p + tau], q @ nd[2][0], [nd[2][2] + tau]]
    )
    rc, _ = p.circ.state(CAL, t_p + tau)
    miss = float(np.linalg.norm((q @ xp)[:3] - rc))
    _log(f"pinseed: t_x {t_p / DAY:.4f} C {tc / DAY:.4f} tau {tau / DAY:.5f} d miss {miss:.2e}")
    (OUT / "a_state.json").write_text(
        json.dumps({"z": zp.tolist(), "from": "pinseed (amendment 3)"})
    )


def stage_gm(s_list: list[float], max_nfev: int) -> None:
    """Amendment 3.2: continuation of the pinned control in Ganymede's GM, s_G from 1 down.

    Softening radius scaled linearly with s (amendment 5). Predictor: secant in log s, or at the
    start the Ganymede-relative node offsets scaled with s (patched-conic r_p is proportional to
    GM at fixed V_inf and turn), velocities kept."""
    p = build_problem()
    mu0, r0 = p.mus[GAN], p.surf[GAN]
    path = OUT / "gm.json"
    rec = json.loads(path.read_text()) if path.exists() else {"points": []}
    good = [q for q in rec["points"] if q["converged"]]
    if good:
        z, s_prev = np.asarray(good[-1]["z"]), float(good[-1]["s"])
    else:
        z, s_prev = np.asarray(json.loads((OUT / "a_state.json").read_text())["z"]), 1.0
    lay = Layout("pin")
    hist = [(float(q["s"]), np.asarray(q["z"])) for q in good[-2:]]
    for s_g in s_list:
        if len(hist) == 2:  # secant predictor in log s
            (s0, z0), (s1, z1) = hist
            zs = z1 + (z1 - z0) * (math.log(s_g) - math.log(s1)) / (math.log(s1) - math.log(s0))
        else:
            zs = z.copy()
            for i0, it in ((0, 6), (11, 17)):
                rg, _ = p.circ.state(GAN, float(zs[it]))
                zs[i0 : i0 + 3] = rg + (zs[i0 : i0 + 3] - rg) * (s_g / s_prev)
        p.mus[GAN] = s_g * mu0
        p.surf[GAN] = s_g * r0  # amendment 5: linear, like r_p
        zn, info = solve_trf(p, zs, lay, max_nfev, f"gm s={s_g:.4g}")
        ok = converged(info)
        pt: dict[str, Any] = {"s": s_g, "converged": ok, "norm": info["norm"], "z": zn.tolist()}
        if ok:
            nd = nodes(p, zn, lay)
            gan = []
            for idx in (0, 2):
                v, d = vinf_gan(p, nd[idx][0], nd[idx][2])
                gan.append({"vinf_kms": v, "rp_km": d, "rp_over_s_km": d / s_g})
            _, vc = p.circ.state(CAL, nd[1][2])
            pt["ganymede"] = gan
            pt["callisto_speed_kms"] = float(np.linalg.norm(nd[1][0][3:] - vc))
            z, s_prev = zn, s_g
            hist = [*hist[-1:], (s_g, zn)]
        rec["points"].append(pt)
        path.write_text(json.dumps(rec))
        msg = json.dumps(pt.get("ganymede")) if ok else ""
        _log(f"gm s={s_g:.4g}: converged={ok} |r| {info['norm']:.3e} {msg}")
        if not ok:
            break


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="gancal5-continuous-ideal-positive-control",
        method=MethodCapability(
            genome="R-S 2009 GanCal#5 (one structure), circular coplanar, Callisto massless",
            corrector="forward-backward multiple shooting, DOP853 + analytic STM",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage",
        choices=["a", "free", "family", "tcont", "pinseed", "gm", "check", "b"],
        required=True,
    )
    ap.add_argument("--max-nfev", type=int, default=40)
    ap.add_argument("--dj", type=float, default=0.0)
    ap.add_argument("--start", default="a_state.json")
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--t-list", default="")
    ap.add_argument("--ds", type=float, default=1e-3)
    ap.add_argument("--direction", type=int, default=1)
    args = ap.parse_args()
    if args.stage == "a":
        stage_a(args.max_nfev)
    elif args.stage == "gm":
        stage_gm([float(v) for v in args.t_list.split(",")], args.max_nfev)
    elif args.stage == "pinseed":
        stage_pinseed()
    elif args.stage == "tcont":
        stage_tcont([float(v) for v in args.t_list.split(",")], args.max_nfev)
    elif args.stage == "family":
        stage_family(args.steps, args.ds, args.direction)
    elif args.stage == "free":
        stage_free(args.max_nfev, args.dj, args.start)
    elif args.stage == "check":
        stage_check()
    else:
        raise SystemExit("stage b: see run_968_lane.py")


if __name__ == "__main__":
    main()
