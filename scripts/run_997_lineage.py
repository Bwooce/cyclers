"""#997: Earth-Moon both-primary lineage. Pseudo-arclength continuation of planar x-axis-symmetric
periodic orbits seeded from Schwaniger 1963 (#970), Newton 1959 Table 1, Newton's other
alpha/beta generating ellipses, and Hoelker & Winston 1968 Fig. 89/90.

Pre-registration (families, continuation, stopping rules, cycler-class test, database gate):
docs/notes/2026-10-07-997-earth-moon-both-primary-lineage.md sec. 0. This script implements it.

Usage (each call stays under --max-seconds; re-running resumes from the checkpoint):
    uv run python scripts/run_997_lineage.py seed --family F2        # build seeds -> seeds JSON
    uv run python scripts/run_997_lineage.py cont --seed F1-0 --dir +1
Output: data/997_lineage/seeds.json and data/997_lineage/<seed>_<dir>.jsonl (one member a line).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

Arr = NDArray[np.float64]
REPO = Path(__file__).resolve().parent.parent
OUTDIR = REPO / "data" / "997_lineage"
SYS = cr3bp_system("Earth", "Moon")
MU_REG = SYS.mu
L_KM = SYS.l_km
TU_S = SYS.t_s
R_E = PLANETS["E"].radius_eq_km
R_M = SATELLITES["Moon"].radius_eq_km
FLOOR_E = R_E + PLANETS["E"].safe_alt_km  # 6578.137 km
FLOOR_M = R_M + SATELLITES["Moon"].safe_alt_km  # 1837.4 km
GEO_KM = 42164.0
HILL_KM = (MU_REG / 3.0) ** (1.0 / 3.0) * L_KM
RTOL = ATOL = 1e-12

# pre-registered continuation controls (note sec. 0.3-0.4)
DS0, DS_MAX, DS_MIN = 2e-3, 2e-2, 1e-6
MAX_NEWTON, EASY_ITERS = 8, 4
G_TOL = 1e-11
# Sensitive orbits (dG/dydot0 ~ 1e4-1e5, e.g. Hoelker-Winston Fig. 89) stall at |G| ~ 1e-9, the
# integration-noise floor; accept a stalled Newton step (|dz| < DZ_STAG) when |G| < G_STAG.
DZ_STAG, G_STAG = 1e-12, 1e-7
MAX_MEMBERS, T_MAX, R1_MAX = 400, 60.0, 5.0
# Collision guard, as a fraction of the physical radius. 1.0 everywhere except the F5 seed
# search, whose published Hoelker-Winston neighbours pass inside the Moon's radius (point-mass
# model); a member below the surface still stops the continuation (rule 1).
GUARD = [1.0]
ALPHA_BETA = [(2, 5), (3, 8), (4, 11), (3, 7), (2, 3), (3, 4)]


def _ts(*a: object) -> None:
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------ planar dynamics


def eom4(t: float, s: Arr, mu: float) -> Arr:
    x, y, vx, vy = s[0], s[1], s[2], s[3]
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    c1 = (1.0 - mu) / r1**3
    c2 = mu / r2**3
    ax = x - c1 * (x + mu) - c2 * (x - 1.0 + mu) + 2.0 * vy
    ay = y - c1 * y - c2 * y - 2.0 * vx
    return np.array([vx, vy, ax, ay])


def eom4_stm(t: float, z: Arr, mu: float) -> Arr:
    x, y = z[0], z[1]
    phi = z[4:].reshape(4, 4)
    xa, xb = x + mu, x - 1.0 + mu
    r1 = math.hypot(xa, y)
    r2 = math.hypot(xb, y)
    om = 1.0 - mu
    r13, r23, r15, r25 = r1**3, r2**3, r1**5, r2**5
    uxx = 1.0 - om / r13 - mu / r23 + 3 * om * xa * xa / r15 + 3 * mu * xb * xb / r25
    uyy = 1.0 - om / r13 - mu / r23 + 3 * om * y * y / r15 + 3 * mu * y * y / r25
    uxy = 3 * om * xa * y / r15 + 3 * mu * xb * y / r25
    a = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [uxx, uxy, 0, 2.0], [uxy, uyy, -2.0, 0]], dtype=float)
    out = np.empty(20)
    out[:4] = eom4(t, z[:4], mu)
    out[4:] = (a @ phi).ravel()
    return out


def jacobi(s: Arr, mu: float) -> float:
    x, y, vx, vy = s[:4]
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    return float(x * x + y * y + 2 * (1 - mu) / r1 + 2 * mu / r2 - vx * vx - vy * vy)


def _yevent(t: float, z: Arr, mu: float) -> float:
    return float(z[1])


def nth_crossing(
    mu: float, x0: float, yd0: float, n: int, t_guess: float, stm: bool
) -> tuple[float, Arr] | None:
    """The n-th y = 0 crossing (t > 0) of (x0, 0, 0, yd0), with the 4 x 4 STM if asked."""
    z0 = np.array([x0, 0.0, 0.0, yd0])
    f: Callable[..., Arr] = eom4
    if stm:
        z0 = np.concatenate([z0, np.eye(4).ravel()])
        f = eom4_stm
    count = [0]
    t_lo = 1e-6  # the start point is itself on y = 0; ignore that root

    def ev(t: float, z: Arr, m: float) -> float:
        return float(z[1])

    # Guard: a pass inside either body (its physical radius) is an impact orbit for this
    # point-mass model; stop instead of grinding the step size down (seen at F4 rho = 0.005).
    def hit_e(t: float, z: Arr, m: float) -> float:
        return math.hypot(z[0] + m, z[1]) - GUARD[0] * R_E / L_KM

    def hit_m(t: float, z: Arr, m: float) -> float:
        return math.hypot(z[0] - 1 + m, z[1]) - GUARD[0] * R_M / L_KM

    hit_e.terminal = True  # type: ignore[attr-defined]
    hit_m.terminal = True  # type: ignore[attr-defined]
    hits: list[tuple[float, Arr]] = []
    t0, zc, t_hi = 0.0, z0, max(1.5 * t_guess, 1.0)
    for _ in range(6):
        sol = solve_ivp(
            f,
            (t0, t_hi),
            zc,
            args=(mu,),
            method="DOP853",
            rtol=RTOL,
            atol=ATOL,
            events=[ev, hit_e, hit_m],
        )
        if sol.status != 0:
            return None
        for te, ye in zip(sol.t_events[0], sol.y_events[0], strict=True):
            if te > t_lo and (not hits or te - hits[-1][0] > 1e-9):
                hits.append((float(te), np.asarray(ye)))
        count[0] = len(hits)
        if len(hits) >= n:
            return hits[n - 1]
        t0, zc = float(sol.t[-1]), sol.y[:, -1]
        t_hi = t0 + max(t_guess, 1.0)
        if t0 > T_MAX:
            return None
    return None


def crossing_times(mu: float, x0: float, yd0: float, t_end: float) -> list[float]:
    """All y = 0 crossing times in (0, t_end], one integration; [] on a collision."""

    def hit_e(t: float, z: Arr, m: float) -> float:
        return math.hypot(z[0] + m, z[1]) - GUARD[0] * R_E / L_KM

    def hit_m(t: float, z: Arr, m: float) -> float:
        return math.hypot(z[0] - 1 + m, z[1]) - GUARD[0] * R_M / L_KM

    hit_e.terminal = True  # type: ignore[attr-defined]
    hit_m.terminal = True  # type: ignore[attr-defined]
    sol = solve_ivp(
        eom4,
        (0.0, t_end),
        np.array([x0, 0.0, 0.0, yd0]),
        args=(mu,),
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=[_yevent, hit_e, hit_m],
    )
    if sol.status != 0:
        return []
    return [float(t) for t in sol.t_events[0] if t > 1e-6]


def residual(mu: float, x0: float, yd0: float, n: int, t_guess: float) -> dict[str, Any] | None:
    """G = xdot at the n-th crossing and dG/d(x0, yd0) with the crossing time projected out."""
    hit = nth_crossing(mu, x0, yd0, n, t_guess, stm=True)
    if hit is None:
        return None
    t, z = hit
    s = z[:4]
    phi = z[4:].reshape(4, 4)
    acc = eom4(t, s, mu)
    if abs(s[3]) < 1e-12:
        return None
    k = acc[2] / s[3]  # xddot / ydot
    grad = np.array([phi[2, 0] - k * phi[1, 0], phi[2, 3] - k * phi[1, 3]])
    return {"G": float(s[2]), "grad": grad, "t": t, "state": s}


# ------------------------------------------------------------------ member diagnostics


def member_report(mu: float, x0: float, yd0: float, t_half: float) -> dict[str, Any]:
    s0 = np.array([x0, 0.0, 0.0, yd0])
    period = 2.0 * t_half

    def d1(t: float, z: Arr, m: float) -> float:
        return float((z[0] + m) * z[2] + z[1] * z[3])

    def d2(t: float, z: Arr, m: float) -> float:
        return float((z[0] - 1 + m) * z[2] + z[1] * z[3])

    z0 = np.concatenate([s0, np.eye(4).ravel()])
    sol = solve_ivp(
        eom4_stm,
        (0, period),
        z0,
        args=(mu,),
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=[d1, d2, _yevent],
        dense_output=True,
    )
    m4 = sol.y[4:, -1].reshape(4, 4)
    sf = sol.y[:4, -1]
    r1e = [math.hypot(z[0] + mu, z[1]) for z in sol.y_events[0]]
    r2e = [math.hypot(z[0] - 1 + mu, z[1]) for z in sol.y_events[1]]
    r1e.append(math.hypot(x0 + mu, 0.0))
    r2e.append(math.hypot(x0 - 1 + mu, 0.0))
    # winding numbers about each primary over one period (dense, refined near passes)
    tt = np.unique(np.concatenate([sol.t, np.linspace(0, period, 4001)]))
    fine = [tt[:1]]
    for a, b in itertools.pairwise(tt):
        fine.append(np.linspace(a, b, 9)[1:])
    tg = np.concatenate(fine)
    zz = sol.sol(tg)
    th1 = np.unwrap(np.arctan2(zz[1], zz[0] + mu))
    th2 = np.unwrap(np.arctan2(zz[1], zz[0] - 1 + mu))
    w1 = (th1[-1] - th1[0]) / (2 * math.pi)
    w2 = (th2[-1] - th2[0]) / (2 * math.pi)
    eig = np.linalg.eigvals(m4)
    n_half = int(sum(1 for t in sol.t_events[2] if 1e-6 < t < t_half + 1e-7))
    j_half = int(np.argmin(np.abs(sol.t_events[2] - t_half)))
    x_half = float(sol.y_events[2][j_half][0])
    r1min, r2min = min(r1e) * L_KM, min(r2e) * L_KM
    c = jacobi(s0, mu)
    return {
        "x0": x0,
        "ydot0": yd0,
        "mu": mu,
        "C": c,
        "C_with_mu1mu": c + mu * (1 - mu),
        "t_half": t_half,
        "T": period,
        "T_days": period * TU_S / 86400.0,
        "perigee_km": r1min,
        "perigee_alt_km": r1min - R_E,
        "apogee_km": max(r1e) * L_KM,
        "periselene_km": r2min,
        "periselene_alt_km": r2min - R_M,
        "n_r1_extrema_below_geo": len([r for r in r1e if r * L_KM < GEO_KM]),
        "b_h": float(np.trace(m4) - 2.0),
        "lambda_max": float(np.max(np.abs(eig))),
        "det_M4": float(np.linalg.det(m4)),
        "closure": float(np.linalg.norm(sf - s0)),
        "crossings_to_half": n_half,
        "wind_E": float(w1),
        "wind_M": float(w2),
        "r1max_L": max(r1e),
        "r2max_km": max(r2e) * L_KM,
        "x_half": x_half,
    }


def vertical_index(mu: float, x0: float, yd0: float, period: float) -> float:
    """b_v = tr of the z/zdot block of the full monodromy (vertical variational equation)."""

    def f(t: float, z: Arr, m: float) -> Arr:
        s = z[:4]
        r1 = math.hypot(s[0] + m, s[1])
        r2 = math.hypot(s[0] - 1 + m, s[1])
        uzz = -(1 - m) / r1**3 - m / r2**3
        p = z[4:].reshape(2, 2)
        a = np.array([[0.0, 1.0], [uzz, 0.0]])
        return np.concatenate([eom4(t, s, m), (a @ p).ravel()])

    z0 = np.concatenate([[x0, 0.0, 0.0, yd0], np.eye(2).ravel()])
    sol = solve_ivp(f, (0, period), z0, args=(mu,), method="DOP853", rtol=RTOL, atol=ATOL)
    return float(np.trace(sol.y[4:, -1].reshape(2, 2)))


def classify(m: dict[str, Any]) -> dict[str, Any]:
    pe, ps = m["perigee_km"], m["periselene_km"]
    out = {
        "passes_earth_floor": pe >= FLOOR_E,
        "passes_moon_floor": ps >= FLOOR_M,
        "earth_pass_below_geo": pe <= GEO_KM,
        "moon_pass_inside_hill": ps <= HILL_KM,
        "below_surface": pe < R_E or ps < R_M,
    }
    out["cycler_class_candidate"] = all(
        out[k]
        for k in (
            "passes_earth_floor",
            "passes_moon_floor",
            "earth_pass_below_geo",
            "moon_pass_inside_hill",
        )
    )
    # Restrepo-Russell 2018 global-search domain (note sec. 0.6): a perpendicular crossing within
    # 5 x_L1 of the Moon, N <= 10 crossings to T/2, max r1 <= 5.
    x_l1 = 0.1509
    xm = 1 - m["mu"]
    out["rr_domain"] = bool(
        abs(m["x0"] - xm) < 5 * x_l1 and m["crossings_to_half"] <= 10 and m["r1max_L"] <= 5.0
    )
    # Franz-Russell 2022 discard any orbit that is ever more than ~350,000 km from the Moon
    # (note sec. 0.6): the test is the maximum Moon distance over the whole period.
    out["franz_russell_domain"] = bool(m["r2max_km"] <= 350000.0)
    return out


# ------------------------------------------------------------------ correctors


def correct_fixed_x0(
    mu: float, x0: float, yd0: float, n: int, t_guess: float, max_iter: int = 30
) -> tuple[float, float] | None:
    for _ in range(max_iter):
        r = residual(mu, x0, yd0, n, t_guess)
        if r is None:
            return None
        if abs(r["G"]) < G_TOL:
            return yd0, r["t"]
        d = r["grad"][1]
        if abs(d) < 1e-14:
            return None
        step = -r["G"] / d
        if abs(step) < DZ_STAG and abs(r["G"]) < G_STAG:
            return yd0, r["t"]  # stagnated at the integration-noise floor (sensitive orbit)
        step = max(-0.05, min(0.05, step))
        yd0 += step
        t_guess = r["t"]
    return None


def arclength_step(
    mu: float, z: Arr, tau: Arr, ds: float, n: int, t_guess: float
) -> tuple[Arr, float, int] | None:
    zp = z + ds * tau
    zc = zp.copy()
    for it in range(1, MAX_NEWTON + 1):
        r = residual(mu, zc[0], zc[1], n, t_guess)
        if r is None:
            return None
        f = np.array([r["G"], float(tau @ (zc - zp))])
        if abs(r["G"]) < G_TOL and it > 1:
            return zc, r["t"], it
        jac = np.vstack([r["grad"], tau])
        try:
            dz = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            return None
        if it > 1 and float(np.linalg.norm(dz)) < DZ_STAG and abs(r["G"]) < G_STAG:
            return zc, r["t"], it  # noise floor of a sensitive orbit
        zc = zc + dz
        t_guess = r["t"]
    r = residual(mu, zc[0], zc[1], n, t_guess)
    if r is not None and abs(r["G"]) < G_TOL:
        return zc, r["t"], MAX_NEWTON
    return None


def tangent(mu: float, z: Arr, n: int, t_guess: float, prev: Arr | None) -> Arr:
    r = residual(mu, z[0], z[1], n, t_guess)
    assert r is not None
    g = r["grad"]
    tau = np.array([-g[1], g[0]])
    tau /= np.linalg.norm(tau)
    if prev is not None and float(tau @ prev) < 0:
        tau = -tau
    return tau


# ------------------------------------------------------------------ seeds

NEWTON_MU = 1.0 / 82.45
NEWTON_ROWS = [  # rho_M0, -ydot0, class (Newton 1959 Table 1, image-checked in the digest)
    (0.012129, 2.3314, "retro"),
    (0.062129, 1.7739, "retro"),
    (0.112129, 1.6604, "retro"),
    (0.162129, 1.5842, "retro"),
    (0.212129, 1.498, "retro"),
    (0.012129, 1.5364, "mixed"),
    (0.03211, 1.0200, "zero-vel"),
    (0.062129, 0.8475, "direct"),
    (0.112129, 0.8303, "direct"),
    (0.162129, 0.9124, "direct"),
    (0.212129, 1.049, "direct"),
]
NEWTON_PERIODS = [
    7.8925,
    6.5227,
    6.3841,
    6.3335,
    6.302,
    5.4292,
    5.5718,
    5.7574,
    5.975,
    6.1081,
    6.192,
]
HW_MU = 1.0 / 80.0


def mu_homotopy(
    x0: float, yd0: float, n: int, t_half: float, mu_from: float, steps: int
) -> tuple[float, float] | None:
    for mu in np.linspace(mu_from, MU_REG, steps + 1)[1:]:
        res = correct_fixed_x0(float(mu), x0, yd0, n, t_half)
        if res is None:
            return None
        yd0, t_half = res
    return yd0, t_half


def build_seeds(family: str, only_ab: str | None = None) -> list[dict[str, Any]]:
    seeds: list[dict[str, Any]] = []
    if family == "F1":
        ctrl = json.loads((REPO / "data" / "970_schwaniger" / "control.json").read_text())
        reg = next(r for r in ctrl["results"] if r["label"] == "registry")
        seeds.append(
            {"id": "F1-0", "x0": reg["x0"], "ydot0": reg["ydot0"], "n": 2, "t_half": reg["t_half"]}
        )
    elif family in ("F2", "F3"):
        for i, ((rho, ydm, cls), tp) in enumerate(zip(NEWTON_ROWS, NEWTON_PERIODS, strict=True)):
            if (family == "F2") != (cls == "retro"):
                continue
            x0 = 1 - NEWTON_MU + rho
            res = correct_fixed_x0(NEWTON_MU, x0, -ydm, 3, tp / 2)
            if res is None:
                _ts(f"  {family} row {i + 1}: no closure at Newton's mu")
                continue
            res2 = mu_homotopy(x0, res[0], 3, res[1], NEWTON_MU, 4)
            if res2 is None:
                _ts(f"  {family} row {i + 1}: mu homotopy failed")
                continue
            seeds.append(
                {
                    "id": f"{family}-{i + 1}",
                    "x0": x0,
                    "ydot0": res2[0],
                    "n": 3,
                    "t_half": res2[1],
                    "src": f"Newton 1959 Table 1 row {i + 1} ({cls})",
                    "newton_mu_ydot0": res[0],
                }
            )
    elif family == "F4":
        mu = MU_REG
        for al, be in ALPHA_BETA:
            if only_ab is not None and only_ab != f"{al}/{be}":
                continue
            a = (al / be) ** (2.0 / 3.0)
            for sense in (+1.0, -1.0):
                # 0.005-0.02 added before any of their results were seen (note sec. 0.7)
                for rho in (0.005, 0.01, 0.02, 0.03, 0.06, 0.1, 0.15, 0.2):
                    ra = 1.0 + rho
                    if ra >= 2 * a:
                        continue
                    va = math.sqrt((1 - mu) * (2 / ra - 1 / a))
                    x0 = ra - mu
                    yd0 = sense * va - ra
                    t_half_guess = math.pi * al
                    # crossing index nearest pi*alpha on the mu = 0 seed path
                    # topology from the mu = 0 rotating-Kepler seed path (note sec. 0.2): the
                    # crossing nearest pi alpha on that path fixes the index n
                    hits = crossing_times(0.0, ra, yd0, 1.3 * t_half_guess)
                    if not hits:
                        continue
                    n = int(np.argmin(np.abs(np.array(hits) - t_half_guess))) + 1
                    res = correct_fixed_x0(mu, x0, yd0, n, hits[n - 1])
                    tag = f"F4-{al}/{be}-{'d' if sense > 0 else 'r'}-rho{rho}"
                    if res is None:
                        _ts(f"  {tag}: n={n} no closure")
                        continue
                    if abs(res[1] - t_half_guess) > 0.15 * t_half_guess:
                        # closed, but not as the alpha/beta type: the half period is far from
                        # pi alpha (typically a lunar capture with many loops); not a seed
                        _ts(f"  {tag}: n={n} closed at t_half={res[1]:.4f}, not the type; dropped")
                        continue
                    _ts(f"  {tag}: n={n} closed, t_half={res[1]:.4f} (pi*alpha={t_half_guess:.4f})")
                    seeds.append(
                        {
                            "id": tag.replace("/", "_"),
                            "x0": x0,
                            "ydot0": res[0],
                            "n": n,
                            "t_half": res[1],
                            "alpha_beta": [al, be],
                            "sense": "direct" if sense > 0 else "retrograde",
                            "rho": rho,
                        }
                    )
    elif family == "F5":
        GUARD[0] = 0.01
        x0 = 1.4875
        for fig, yd in (("89", -0.848358), ("90", -0.84835)):
            for n in range(1, 11):
                h = nth_crossing(HW_MU, x0, yd, n, 5.0, stm=False)
                if h is None:
                    break
                res = correct_fixed_x0(HW_MU, x0, yd, n, h[0])
                if res is None or abs(res[0] - yd) > 2e-3:
                    continue
                res2 = mu_homotopy(x0, res[0], n, res[1], HW_MU, 8)
                if res2 is None:
                    _ts(
                        f"  F5 Fig {fig} n={n}: closed at 1/80 (ydot {res[0]:.7f}), "
                        "mu homotopy failed"
                    )
                    continue
                _ts(f"  F5 Fig {fig} n={n}: ydot(1/80)={res[0]:.7f} t_half={res[1]:.4f}")
                seeds.append(
                    {
                        "id": f"F5-fig{fig}-n{n}",
                        "x0": x0,
                        "ydot0": res2[0],
                        "n": n,
                        "t_half": res2[1],
                        "hw_mu_ydot0": res[0],
                        "hw_mu_t_half": res[1],
                    }
                )
    return seeds


# ------------------------------------------------------------------ continuation driver


def topology(m: dict[str, Any]) -> tuple[int, int, int]:
    return (m["crossings_to_half"], round(m["wind_E"]), round(m["wind_M"]))


def continue_seed(seed: dict[str, Any], direction: int, max_seconds: float) -> str:
    out = OUTDIR / f"{seed['id']}_{'p' if direction > 0 else 'm'}.jsonl"
    rows = [json.loads(s) for s in out.read_text().splitlines()] if out.exists() else []
    t_start = time.time()
    mu, n = MU_REG, seed["n"]
    if rows and rows[-1].get("stop"):
        return f"already stopped: {rows[-1]['stop']}"
    if rows:
        last = rows[-1]
        z = np.array([last["x0"], last["ydot0"]])
        t_half, ds = last["t_half"], last["ds"]
        tau = np.array(last["tau"])
        topo0 = tuple(rows[0]["topology"])
    else:
        z = np.array([seed["x0"], seed["ydot0"]])
        t_half, ds = seed["t_half"], DS0
        tau = tangent(mu, z, n, t_half, None) * direction
        m = member_report(mu, z[0], z[1], t_half)
        m.update(classify(m))
        m["b_v"] = vertical_index(mu, z[0], z[1], m["T"])
        m.update({"k": 0, "ds": ds, "tau": tau.tolist(), "topology": list(topology(m))})
        topo0 = topology(m)
        rows.append(m)
        with out.open("a") as fh:
            fh.write(json.dumps(m) + "\n")
    easy = 0
    stop = None

    def write(rec: dict[str, Any]) -> None:
        with out.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()

    while time.time() - t_start < max_seconds:
        k = rows[-1]["k"] + 1
        if k > MAX_MEMBERS:
            stop = "max members"
            break
        res = arclength_step(mu, z, tau, ds, n, t_half)
        if res is None:
            ds *= 0.5
            easy = 0
            if ds < DS_MIN:
                stop = "step failure (ds < ds_min)"
                break
            continue
        z_new, t_new, iters = res
        m = member_report(mu, z_new[0], z_new[1], t_new)
        m.update(classify(m))
        topo = topology(m)
        prev = rows[-1]
        if topo != topo0:
            if ds > 4 * DS_MIN:  # refine onto the topology change before stopping
                ds *= 0.25
                continue
            stop = f"topology change {topo0} -> {topo}"
        if m["T"] > T_MAX or m["r1max_L"] > R1_MAX:
            stop = "escape/runaway (T or max r1)"
        # Floors (note sec. 0.4 rule 1): stop when a member drops from above a floor to below
        # it. A run that starts below a floor (F1: 177 km perigee) is labelled and runs on to
        # the physical surface.
        was_ok = prev["passes_earth_floor"] and prev["passes_moon_floor"]
        now_ok = m["passes_earth_floor"] and m["passes_moon_floor"]
        if m["below_surface"]:
            stop = "impact at the physical surface"
        elif was_ok and not now_ok:
            stop = "impact at floor"
        if (prev["b_h"] + 2.0) * (m["b_h"] + 2.0) < 0:
            stop = "period doubling (b_h crosses -2)"
        # The start crossing and the T/2 crossing swap: the curve has passed an orbit with
        # period T/2 (a period-doubling bifurcation of a half-period family, where b_h -> +2)
        # and now retraces the same orbits from the other symmetric crossing.
        if (prev["x0"] - prev["x_half"]) * (m["x0"] - m["x_half"]) < 0:
            stop = "start and T/2 crossings swap (period-halving end point; retrace beyond)"
        m["b_h_crossed_plus2"] = bool((prev["b_h"] - 2.0) * (m["b_h"] - 2.0) < 0)
        m["b_v"] = vertical_index(mu, z_new[0], z_new[1], m["T"])
        new_tau = tangent(mu, z_new, n, t_new, tau)
        m.update(
            {"k": k, "ds": ds, "tau": new_tau.tolist(), "topology": list(topo), "iters": iters}
        )
        if stop:
            m["stop"] = stop
        rows.append(m)
        write(m)
        if stop:
            break
        z, t_half, tau = z_new, t_new, new_tau
        easy = easy + 1 if iters <= EASY_ITERS else 0
        if easy >= 3:
            ds = min(ds * 1.5, DS_MAX)
            easy = 0
    last = rows[-1]
    return (
        f"{out.name}: {len(rows)} members, stop={stop}, last C={last['C']:.6f} "
        f"T={last['T']:.4f} perigee {last['perigee_alt_km']:.0f} km alt, "
        f"periselene {last['periselene_alt_km']:.0f} km alt"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s1 = sub.add_parser("seed")
    s1.add_argument("--family", required=True, choices=["F1", "F2", "F3", "F4", "F5"])
    s1.add_argument("--ab", default=None, help="F4 only: one alpha/beta, e.g. 2/5")
    s2 = sub.add_parser("cont")
    s2.add_argument("--seed", required=True)
    s2.add_argument("--dir", type=int, choices=[1, -1], required=True)
    s2.add_argument("--max-seconds", type=float, default=420.0)
    args = ap.parse_args()
    region = args.family if args.cmd == "seed" else args.seed
    if args.cmd == "seed" and args.ab:
        region += "-" + args.ab.replace("/", "_")
    preflight_search(
        task_no=997,
        region_id=f"em-both-primary-lineage-{args.cmd}-{region}",
        method=MethodCapability(
            genome="planar x-axis-symmetric periodic orbit, CR3BP Earth-Moon (registry mu)",
            corrector="perpendicular-crossing Newton + pseudo-arclength in (x0, ydot0)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=MAX_MEMBERS,
    )
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seeds_path = OUTDIR / "seeds.json"
    seeds = json.loads(seeds_path.read_text()) if seeds_path.exists() else {}
    _ts(f"mu={MU_REG!r} L={L_KM} km TU={TU_S} s floors E {FLOOR_E} km M {FLOOR_M} km")
    if args.cmd == "seed":
        new = build_seeds(args.family, args.ab)
        if args.ab is None:
            seeds[args.family] = new
        else:  # merge one alpha/beta into the family's seed list
            keep = [
                x
                for x in seeds.get(args.family, [])
                if x.get("alpha_beta") != [int(v) for v in args.ab.split("/")]
            ]
            seeds[args.family] = keep + new
        seeds_path.write_text(json.dumps(seeds, indent=1) + "\n")
        _ts(f"{args.family}: {len(seeds[args.family])} seeds")
        return
    seed = next(s for fam in seeds.values() for s in fam if s["id"] == args.seed)
    _ts(continue_seed(seed, args.dir, args.max_seconds))


if __name__ == "__main__":
    main()
