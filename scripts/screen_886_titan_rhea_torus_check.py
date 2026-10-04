"""#886 part 1 -- invariant circles of the Saturn-Titan-Rhea CCR4BP against Kumar (IAC-25-C1.9.6).

The published check: B. Kumar, "Analysis of Unstable Resonant Orbits for Saturn Tour Design:
Between Titan and Rhea", IAC-25-C1.9.6 (filed in the private paper corpus as
``kumar-2025-analysis-unstable-resonant-orbits-saturn-tour-design-titan-rhea-IAC-25-C1.9.6-doi-10.52202-083087-0076.pdf``),
Section 6: the Titan 3:2 low/mid-e unstable resonant family, continued in Rhea's mass from
0 to mu3 = 4.05746e-6 as invariant circles of the stroboscopic map, persists for most Jacobi
constants and fails near the secondary resonances Tp/T = 4/21, 9/47, 5/26, 6/31, 7/36, 8/41
(Tp/T over the family from 0.1898475 to 0.1952735).

This script tests the #882 invariant-circle corrector
(:mod:`cyclerfinder.search.ccr4bp_strob_connection`, used read-only) against that statement.
Pre-registered criteria: ``docs/notes/2026-10-04-886-titan-rhea-torus-published-check.md``.

Stages (each well under 8 minutes; the scan stage is chunked and checkpointed)::

    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage constants
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage family
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage members
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage scan --budget-s 420
        (repeat until it reports nothing left; appends to scan.jsonl)
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage refine --budget-s 420
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage closure --budget-s 420
    uv run python scripts/screen_886_titan_rhea_torus_check.py --stage summary

Outputs: ``data/found/886_titan_rhea_torus/``.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from fractions import Fraction  # noqa: E402
from itertools import pairwise  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq, minimize_scalar  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp as ccr4bp  # noqa: E402
import cyclerfinder.search.ccr4bp_strob_connection as sc  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "886_titan_rhea_torus"

# ---------------------------------------------------------------------------
# Constants exactly as printed by the paper (pp. 2, 12).
# ---------------------------------------------------------------------------
MU = 2.36639e-4  # Saturn-Titan mass ratio m2/(m1+m2)
MU3 = 4.05746e-6  # Rhea, m3/(m1+m2): the same definition as CCR4BPSystem.mu_gan
TP = 2.48376  # Rhea's Titan-relative synodic period (Titan time units)
R13_PRINTED = 0.4315  # "about 0.4315 normalized units"
OMEGA = 2.0 * math.pi / TP  # Rhea is inside Titan's orbit, so it advances: omega > 0
A3 = ((1.0 - MU + MU3) / (1.0 + OMEGA) ** 2) ** (1.0 / 3.0)  # Kepler radius for that rate
TITAN_RADIUS = SATELLITES["Titan"].radius_eq_km / SATELLITES["Titan"].sma_km
LISTED = ((4, 21), (9, 47), (5, 26), (6, 31), (7, 36), (8, 41))
RATIO_LO, RATIO_HI = 0.1898475, 0.1952735

# Pre-registered numbers (see the note, section 4).
TOL = 1e-10  # node invariance residual required at every continuation step
OFFNODE_TOL = 1e-8  # residual at half-node angles of the final circle
TAIL_TOL = 1e-6  # top-quarter Fourier tail of the final circle
WINDOW = 1.5e-4  # agreement window in Tp/T around each listed ratio
TWIST_BAND = 2e-5  # exclusion band in Tp/T around each extremum of T along the family
N_LEVELS = (151, 201, 251, 321, 401, 501, 601, 701)
H0 = 0.05  # first continuation step, fraction of MU3
H_MAX = 0.25
H_FLOOR = 1.0 / 1024.0
OFFSETS = (0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4)


def system(mu_gan: float) -> ccr4bp.CCR4BPSystem:
    """The paper's Saturn-Titan-Rhea system in the project's CCR4BP (perturber phase 0 at t=0)."""
    return ccr4bp.CCR4BPSystem(mu=MU, mu_gan=mu_gan, a_gan=A3, omega_gan=OMEGA, theta_gan0=0.0)


def log(msg: str, t_start: float) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{stamp} +{time.time() - t_start:7.1f}s] {msg}", flush=True)


def jacobi(s4: np.ndarray, mu: float = MU) -> float:
    x, y, vx, vy = (float(v) for v in s4[:4])
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    return x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - vx * vx - vy * vy


def farey_in_range(lo: float, hi: float, qmax: int) -> list[tuple[int, int]]:
    """Reduced fractions p/q in (lo, hi) with q <= qmax, sorted by value."""
    out = set()
    for q in range(1, qmax + 1):
        for p in range(math.ceil(lo * q), math.floor(hi * q) + 1):
            f = Fraction(p, q)
            if lo < f < hi:
                out.add((f.numerator, f.denominator))
    return sorted(out, key=lambda pq: pq[0] / pq[1])


def nearest_rational(r: float, qmax: int) -> tuple[int, int, float]:
    best = (0, 1, math.inf)
    for q in range(1, qmax + 1):
        p = round(r * q)
        d = abs(r - p / q)
        if d < best[2] - 1e-15:
            f = Fraction(p, q)
            best = (f.numerator, f.denominator, d)
    return best


# ---------------------------------------------------------------------------
# Three-body family (mu3 = 0): x-axis-symmetric, periapse on the -x axis.
# ---------------------------------------------------------------------------


def orbit(x: float, vy_guess: float, t_guess: float) -> tuple[np.ndarray, float, float]:
    s4, period, res = sc.symmetric_periodic_orbit(MU, x, vy_guess, 0.5 * t_guess)
    return s4, period, res


def monodromy(s4: np.ndarray, period: float) -> np.ndarray:
    y0 = np.concatenate([s4, np.eye(4).reshape(-1)])
    sol = solve_ivp(
        sc.planar_rhs_batch,
        (0.0, period),
        y0,
        args=(system(0.0), 1, True),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    return np.asarray(sol.y[4:, -1].reshape(4, 4))


def orbit_geometry(s4: np.ndarray, period: float, n: int = 6000) -> dict[str, float]:
    """Minimum distance to Titan, periapse/apoapse radius about Saturn, along one period."""
    sol = solve_ivp(
        sc.planar_rhs_batch,
        (0.0, period),
        s4,
        args=(system(0.0), 1, False),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    tt = np.linspace(0.0, period, n)
    yy = sol.sol(tt)
    rs = np.hypot(yy[0] + MU, yy[1])
    rt = np.hypot(yy[0] - 1.0 + MU, yy[1])
    return {
        "min_dist_titan": float(rt.min()),
        "r_min_saturn": float(rs.min()),
        "r_max_saturn": float(rs.max()),
        "min_gap_to_rhea_orbit": float(np.abs(rs - A3).min()),
    }


def family_row(x: float, s4: np.ndarray, period: float, res: float) -> dict[str, Any]:
    m = monodromy(s4, period)
    k = float(np.trace(m)) - 2.0  # = lam + 1/lam for the hyperbolic pair
    lam = float(np.max(np.abs(np.linalg.eigvals(m))))
    row = {
        "x": x,
        "vy": float(s4[3]),
        "T": period,
        "C": jacobi(s4),
        "ratio": TP / period,
        "lam_po": lam,
        "stab_index": k,
        "perp_residual": res,
    }
    row.update(orbit_geometry(s4, period))
    return row


def stage_constants(t_start: float) -> dict[str, Any]:
    gs = PRIMARIES["Saturn"]
    gt = SATELLITES["Titan"].mu_km3_s2
    gr = SATELLITES["Rhea"].mu_km3_s2
    a_t = SATELLITES["Titan"].sma_km
    a_r = SATELLITES["Rhea"].sma_km
    eph_mu = gt / gs
    # gs is the Saturn SYSTEM GM (planet plus moons), the same convention the project's
    # cr3bp systems use; m1+m2 and the system differ by Rhea and the smaller moons (~1e-5).
    eph_mu3 = gr / gs
    eph_a3 = a_r / a_t
    om_r13 = math.sqrt((1.0 - MU + MU3) / R13_PRINTED**3) - 1.0
    om_eph = math.sqrt((1.0 - MU + MU3) / eph_a3**3) - 1.0
    out = {
        "paper": {"mu": MU, "mu3": MU3, "Tp": TP, "r13": R13_PRINTED},
        "used": {"mu": MU, "mu_gan": MU3, "omega_gan": OMEGA, "a_gan": A3, "Tp_model": TP},
        "Tp_from_printed_r13": 2.0 * math.pi / om_r13,
        "Tp_from_ephemeris_a3": 2.0 * math.pi / om_eph,
        "a3_from_Tp": A3,
        "ephemeris": {
            "GM_saturn_system": gs,
            "GM_titan": gt,
            "GM_rhea": gr,
            "mu_titan_over_saturn_system": eph_mu,
            "mu3_rhea": eph_mu3,
            "a_rhea_over_a_titan": eph_a3,
            "rel_diff_mu": eph_mu / MU - 1.0,
            "rel_diff_mu3": eph_mu3 / MU3 - 1.0,
            "rel_diff_a3": eph_a3 / A3 - 1.0,
        },
        "titan_radius_titan_units": TITAN_RADIUS,
        "farey_q_lt_50_in_range": farey_in_range(RATIO_LO, RATIO_HI, 49),
        "farey_q_le_100_in_range": farey_in_range(RATIO_LO, RATIO_HI, 100),
    }
    log(
        f"Tp(model) {2 * math.pi / system(MU3).omega_gan:.6f}; Tp from r13=0.4315: "
        f"{out['Tp_from_printed_r13']:.6f}; from ephemeris a3: {out['Tp_from_ephemeris_a3']:.6f}",
        t_start,
    )
    log(f"Farey q<50 in range: {out['farey_q_lt_50_in_range']}", t_start)
    return out


def _continue_family(
    x0: float, vy0: float, t0: float, direction: float, stop: Any, t_start: float
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    x, vy, tg = x0, vy0, t0
    h = 0.0025
    prev: dict[str, Any] | None = None
    prev2: dict[str, Any] | None = None
    while h >= 2e-5:
        s4, period, res = orbit(x, vy, tg)
        ok = res < 1e-11
        if ok and prev is not None:
            ok = abs(period - prev["T"]) < 0.03 and abs(jacobi(s4) - prev["C"]) < 0.005
        if not ok:
            if prev is None:
                raise RuntimeError("family seed did not converge")
            h *= 0.5
            x = prev["x"] + direction * h
            if prev2 is not None:
                fr = h / abs(prev["x"] - prev2["x"])
                vy = prev["vy"] + fr * (prev["vy"] - prev2["vy"])
                tg = prev["T"] + fr * (prev["T"] - prev2["T"])
            else:
                vy, tg = prev["vy"], prev["T"]
            continue
        row = family_row(x, s4, period, res)
        rows.append(row)
        if len(rows) % 10 == 0:
            log(
                f"family x {x:.5f} C {row['C']:.6f} T {period:.6f} Tp/T {row['ratio']:.7f} "
                f"lam {row['lam_po']:.3f} dTitan {row['min_dist_titan']:.4f}",
                t_start,
            )
        if stop(row):
            break
        prev2, prev = prev, row
        h = min(h * 1.3, 0.0025)
        x = row["x"] + direction * h
        if prev2 is not None:
            fr = h / abs(prev["x"] - prev2["x"])
            vy = prev["vy"] + fr * (prev["vy"] - prev2["vy"])
            tg = prev["T"] + fr * (prev["T"] - prev2["T"])
        else:
            vy, tg = row["vy"], row["T"]
    return rows


def stage_family(t_start: float, reuse: dict[str, Any] | None = None) -> dict[str, Any]:
    if reuse is not None:
        rows = reuse["rows"]
        log(f"family: reusing {len(rows)} stored rows", t_start)
    else:
        a = (2.0 / 3.0) ** (2.0 / 3.0)
        rp = 0.70
        vin = math.sqrt((1.0 - MU) * (2.0 / rp - 1.0 / a))
        x0 = -rp - MU
        s4, period, res = orbit(x0, -vin - x0, 2.0 * math.pi * 2.0)
        assert res < 1e-11
        # toward the near-circular end: stop past the stability boundary (lam_po -> 1)
        low_e = _continue_family(
            x0, float(s4[3]), period, -1.0, lambda r: r["stab_index"] < 2.0 - 1e-3, t_start
        )
        # toward the high-e end: stop at Titan's surface or when continuation fails
        high_e = _continue_family(
            x0, float(s4[3]), period, +1.0, lambda r: r["min_dist_titan"] < TITAN_RADIUS, t_start
        )
        rows = sorted(low_e + high_e[1:], key=lambda r: r["x"])
    xs = np.array([r["x"] for r in rows])
    vys = np.array([r["vy"] for r in rows])
    ts = np.array([r["T"] for r in rows])

    def solve_at(x: float) -> dict[str, Any]:
        tg = float(np.interp(x, xs, ts))
        s4, period, res = orbit(x, float(np.interp(x, xs, vys)), tg)
        if res > 1e-11 or abs(period - tg) > 0.01:
            raise RuntimeError(f"family member at x={x} left the branch (T {period} vs {tg})")
        return family_row(x, s4, period, res)

    # stability boundary (near-circular end)
    stab = [(r["x"], r["stab_index"] - 2.0) for r in rows]
    boundary = None
    for (xa, fa), (xb, fb) in pairwise(stab):
        if fa * fb < 0:
            xb_ = brentq(lambda x: solve_at(x)["stab_index"] - 2.0, xa, xb, xtol=1e-10)
            boundary = solve_at(xb_)
    # maximum of C (the fold of Fig. 2)
    i_c = int(np.argmax([r["C"] for r in rows]))
    opt_c = minimize_scalar(
        lambda x: -solve_at(x)["C"],
        bounds=(xs[i_c - 1], xs[i_c + 1]),
        method="bounded",
        options={"xatol": 1e-9},
    )
    c_max = solve_at(float(opt_c.x))
    # extrema of T along the family
    extrema = []
    for i in range(1, len(rows) - 1):
        if (ts[i] - ts[i - 1]) * (ts[i + 1] - ts[i]) < 0:
            sign = 1.0 if ts[i] < ts[i - 1] else -1.0
            opt = minimize_scalar(
                lambda x, s=sign: s * solve_at(x)["T"],
                bounds=(xs[i - 1], xs[i + 1]),
                method="bounded",
                options={"xatol": 1e-9},
            )
            e_row = solve_at(float(opt.x))
            e_row["kind"] = "min_T" if sign > 0 else "max_T"
            extrema.append(e_row)
    # where the family reaches the paper's lower Tp/T bound
    paper_end = None
    for i in range(len(rows) - 1, 0, -1):
        if (rows[i]["ratio"] - RATIO_LO) * (rows[i - 1]["ratio"] - RATIO_LO) < 0:
            xe = brentq(lambda x: solve_at(x)["ratio"] - RATIO_LO, xs[i - 1], xs[i], xtol=1e-11)
            paper_end = solve_at(xe)
            paper_end["titan_altitude_km"] = (
                paper_end["min_dist_titan"] - TITAN_RADIUS
            ) * SATELLITES["Titan"].sma_km
            break
    log(
        f"stability boundary: {boundary and {k: boundary[k] for k in ('x', 'C', 'T', 'ratio')}}",
        t_start,
    )
    log(f"C maximum: x {c_max['x']:.8f} C {c_max['C']:.8f}", t_start)
    for e in extrema:
        log(
            f"extremum {e['kind']}: x {e['x']:.9f} C {e['C']:.7f} T {e['T']:.9f} "
            f"Tp/T {e['ratio']:.9f}",
            t_start,
        )
    if paper_end:
        log(
            f"paper's lower Tp/T bound at x {paper_end['x']:.7f} C {paper_end['C']:.6f} "
            f"Titan altitude {paper_end['titan_altitude_km']:.0f} km",
            t_start,
        )
    return {
        "rows": rows,
        "stability_boundary": boundary,
        "c_max": c_max,
        "extrema": extrema,
        "paper_lower_ratio_point": paper_end,
        "end_high_e": rows[-1],
    }


# ---------------------------------------------------------------------------
# Family lookup used by the later stages.
# ---------------------------------------------------------------------------


class Family:
    def __init__(self, data: dict[str, Any]) -> None:
        bx = data["stability_boundary"]["x"]
        self.rows = [r for r in data["rows"] if r["x"] >= bx]
        self.x = np.array([r["x"] for r in self.rows])
        self.extrema = data["extrema"]
        self.boundary = data["stability_boundary"]

    def at(self, x: float) -> tuple[np.ndarray, float]:
        vy = float(np.interp(x, self.x, [r["vy"] for r in self.rows]))
        tg = float(np.interp(x, self.x, [r["T"] for r in self.rows]))
        s4, period, res = orbit(x, vy, tg)
        if res > 1e-11:
            raise RuntimeError(f"orbit at x={x} did not converge ({res:.1e})")
        return s4, period

    def ratio(self, x: float) -> float:
        return TP / self.at(x)[1]

    def segments(self) -> list[tuple[float, float]]:
        """Monotone segments of Tp/T in x, between the boundary, the extrema and the end."""
        cuts = [float(self.x[0]), *sorted(e["x"] for e in self.extrema), float(self.x[-1])]
        return list(pairwise(cuts))

    def solve_ratio(self, target: float, seg: tuple[float, float]) -> float | None:
        """Member ``x`` in the monotone segment ``seg`` with ``Tp/T = target`` (None if absent).

        Bracketed from the stored family rows, then a safeguarded secant on the true family.
        """
        a, b = seg
        ratios = np.array([r["ratio"] for r in self.rows])
        inside = (self.x >= a - 1e-12) & (self.x <= b + 1e-12)
        xs, rs = self.x[inside], ratios[inside]
        xs = np.concatenate([[a], xs, [b]])
        rs = np.concatenate([[self.ratio(a)], rs, [self.ratio(b)]])
        idx = np.nonzero((rs[:-1] - target) * (rs[1:] - target) <= 0)[0]
        if idx.size == 0:
            return None
        i = int(idx[0])
        lo, hi = float(xs[i]), float(xs[i + 1])
        flo, fhi = float(rs[i] - target), float(rs[i + 1] - target)
        if flo == 0.0:
            return lo
        x_new = lo - flo * (hi - lo) / (fhi - flo)
        for _ in range(40):
            f_new = self.ratio(x_new) - target
            if abs(f_new) < 1e-13 or hi - lo < 1e-13:
                return x_new
            if f_new * flo < 0:
                hi, fhi = x_new, f_new
            else:
                lo, flo = x_new, f_new
            x_sec = lo - flo * (hi - lo) / (fhi - flo)
            # Illinois-style safeguard: keep the iterate inside the bracket
            x_new = x_sec if lo < x_sec < hi else 0.5 * (lo + hi)
        return x_new


def _solve_member(args: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    fam_data, m = args
    fam = Family(fam_data)
    m = dict(m)
    m["x"] = fam.solve_ratio(m["target_ratio"], fam.segments()[m["segment"]])
    return m


def stage_members(fam_data: dict[str, Any], workers: int, t_start: float) -> list[dict[str, Any]]:
    fam = Family(fam_data)
    members: list[dict[str, Any]] = []
    # 1) uniform grid in x across the unstable family
    xg = np.arange(fam.x[0] + 0.0005, fam.x[-1] - 1e-6, 0.004)
    for x in xg:
        members.append({"kind": "grid", "x": float(x)})
    # 2) every crossing of a listed ratio, with signed offsets in Tp/T
    segs = fam.segments()
    seg_rng = [(min(fam.ratio(a), fam.ratio(b)), max(fam.ratio(a), fam.ratio(b))) for a, b in segs]
    jobs = []
    for si, (lo, hi) in enumerate(seg_rng):
        for p, q in LISTED:
            if not lo < p / q < hi:
                continue
            for d in OFFSETS:
                for sgn in (1.0, -1.0) if d > 0 else (1.0,):
                    t = p / q + sgn * d
                    if lo < t < hi:
                        jobs.append(
                            {
                                "kind": "crossing",
                                "p": p,
                                "q": q,
                                "segment": si,
                                "offset": sgn * d,
                                "target_ratio": t,
                            }
                        )
    log(f"members: solving {len(jobs)} crossing/offset members on {workers} workers", t_start)
    with Pool(workers) as pool:
        solved = pool.map(_solve_member, [(fam_data, j) for j in jobs])
    members.extend(m for m in solved if m["x"] is not None)
    # 3) points at the twistless extrema themselves
    for e in fam.extrema:
        members.append({"kind": "extremum", "x": e["x"], "extremum": e["kind"]})
    for i, m in enumerate(members):
        m["id"] = i
    log(f"{len(members)} members", t_start)
    for m in members:
        if m["kind"] == "crossing" and m["offset"] == 0.0:
            log(f"crossing {m['p']}/{m['q']} segment {m['segment']} at x {m['x']:.9f}", t_start)
    return members


# ---------------------------------------------------------------------------
# The test proper: continuation of one member's circle in Rhea's mass.
# ---------------------------------------------------------------------------


def offnode_residual(circle: sc.InvariantCircle) -> float:
    """``max |F(u(theta)) - u(theta + rho)|`` at the half-node angles (not collocated)."""
    th = circle.thetas + math.pi / circle.n_nodes
    pts = circle.state(th)
    orb = sc.strob_iterates(circle.system, pts, n=1, t0=circle.t0)
    return float(np.max(np.abs(orb.states[1] - circle.state(th + circle.rho))))


def next_odd(n: float) -> int:
    k = math.ceil(n)
    return k if k % 2 == 1 else k + 1


def choose_nodes(
    s4: np.ndarray, period: float, levels: tuple[int, ...]
) -> tuple[sc.InvariantCircle, float, list[tuple[int, float]]]:
    """Smallest node count whose mu3 = 0 seed (the sampled periodic orbit) is invariant to TOL.

    Also returns the seed residual at every level tried (the decay used to estimate the node
    count an unresolved member would need).
    """
    c0 = None
    res = math.inf
    tried: list[tuple[int, float]] = []
    for n in levels:
        c0 = sc.seed_circle_from_periodic_orbit(system(0.0), s4, period, n)
        res, _ = sc.circle_residual(c0)
        tried.append((n, res))
        if res <= TOL:
            break
    assert c0 is not None
    return c0, res, tried


def estimate_nodes_needed(tried: list[tuple[int, float]]) -> float | None:
    """Node count at which a log-linear fit of seed residual against N reaches TOL."""
    pts = [(n, math.log10(r)) for n, r in tried if r > 0 and math.isfinite(r)]
    if len(pts) < 3:
        return None
    n_arr = np.array([p[0] for p in pts], dtype=float)
    l_arr = np.array([p[1] for p in pts])
    slope, icpt = np.polyfit(n_arr, l_arr, 1)
    if slope >= 0:
        return None
    return float((math.log10(TOL) - icpt) / slope)


def continue_member(
    x: float,
    fam: Family,
    *,
    levels: tuple[int, ...] = N_LEVELS,
    h0: float = H0,
    h_floor: float = H_FLOOR,
    harmonic_q: int | None = None,
    with_bundles: bool = True,
    target: float = MU3,
) -> dict[str, Any]:
    """Continue the invariant circle of member ``x`` from mu3 = 0 to MU3 at fixed rho.

    Adaptive natural-parameter steps with a secant predictor; a failed step is halved and
    retried from the last converged circle, down to ``h_floor`` (fraction of MU3).
    """
    t_start = time.time()
    s4, period = fam.at(x)
    ratio = TP / period
    c_seed, seed_res, tried = choose_nodes(s4, period, levels)
    out: dict[str, Any] = {
        "x": x,
        "T": period,
        "ratio": ratio,
        "C": jacobi(s4),
        "n_nodes": c_seed.n_nodes,
        "seed_residual": seed_res,
        "seed_tail": c_seed.fourier_tail(),
        "h0": h0,
        "h_floor": h_floor,
        "target_mu3": target,
        "seed_levels": tried,
    }
    if seed_res > TOL:
        out.update(
            status="seed_unresolved",
            frac_reached=0.0,
            n_needed_estimate=estimate_nodes_needed(tried),
        )
        return out
    cur = c_seed
    prev: sc.InvariantCircle | None = None
    frac, prev_frac, h = 0.0, 0.0, h0
    n_steps = n_fail = n_eval = 0
    last_fail: dict[str, Any] | None = None
    while frac < 1.0:
        f_try = min(1.0, frac + h)
        guess = cur.nodes
        if prev is not None and frac > prev_frac:
            guess = cur.nodes + (cur.nodes - prev.nodes) * (f_try - frac) / (frac - prev_frac)
        nxt = sc.correct_invariant_circle(
            system(f_try * target),
            guess,
            cur.rho,
            t0=0.0,
            tol=TOL,
            max_iter=12,
            phase_ref=cur.nodes,
        )
        n_eval += len(nxt.residual_history)
        if nxt.converged:
            prev, prev_frac = cur, frac
            cur, frac = nxt, f_try
            n_steps += 1
            h = min(h * 1.5, H_MAX)
        else:
            n_fail += 1
            last_fail = {"frac_try": f_try, "history": list(nxt.residual_history)}
            h *= 0.5
            if h < h_floor:
                break
    out.update(
        frac_reached=frac,
        n_steps=n_steps,
        n_failed_steps=n_fail,
        n_evals=n_eval,
        last_fail=last_fail,
        reached_target=frac >= 1.0,
    )
    final = cur
    out["final_residual"] = float(final.residual) if frac > 0 else float(seed_res)
    out["final_tail"] = final.fourier_tail()
    out["offnode_residual"] = offnode_residual(final)
    amp_d = sc.fourier_amplitudes(final.nodes - c_seed.nodes)
    out["deformation_max"] = float(np.max(np.abs(final.nodes - c_seed.nodes)))
    if harmonic_q is not None and harmonic_q < amp_d.size:
        out["harmonic_q"] = harmonic_q
        out["deformation_harmonic_q"] = float(amp_d[harmonic_q])
    out["deformation_peak_harmonic"] = int(np.argmax(amp_d[1:]) + 1)
    if with_bundles and frac > 0:
        try:
            b = sc.hyperbolic_bundles(final)
            out["lam_u_map"] = float(b.lam_u)
            out["lam_s_map"] = float(b.lam_s)
        except Exception as exc:
            out["bundle_error"] = repr(exc)
    persisted = (
        out["reached_target"]
        and out["final_residual"] <= TOL
        and out["offnode_residual"] <= OFFNODE_TOL
        and out["final_tail"] <= TAIL_TOL
    )
    out["status"] = "persist" if persisted else ("unresolved" if out["reached_target"] else "fail")
    out["wall_s"] = time.time() - t_start
    if persisted:
        out["final_nodes"] = final.nodes.tolist()
    return out


def _worker(args: tuple[dict[str, Any], dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    member, fam_data, kw = args
    fam = Family(fam_data)
    q = member.get("q")
    if q is None:
        q = nearest_rational(fam.ratio(member["x"]), 49)[1]
    kw2 = dict(kw)
    if "target" in member:
        kw2["target"] = member["target"]
    try:
        res = continue_member(member["x"], fam, harmonic_q=q, **kw2)
    except Exception as exc:
        res = {"x": member["x"], "status": "crash", "error": repr(exc)}
    res["id"] = member["id"]
    res["member"] = member
    return res


def _done_ids(path: Path) -> set[int]:
    if not path.exists():
        return set()
    return {json.loads(line)["id"] for line in path.read_text().splitlines() if line.strip()}


def run_pool(
    todo: list[dict[str, Any]],
    fam_data: dict[str, Any],
    kw: dict[str, Any],
    path: Path,
    budget_s: float,
    workers: int,
    t_start: float,
) -> None:
    done = 0
    with Pool(workers) as pool:
        it = pool.imap_unordered(_worker, [(m, fam_data, kw) for m in todo])
        for res in it:
            with path.open("a") as fh:
                fh.write(json.dumps(res) + "\n")
                fh.flush()
            done += 1
            log(
                f"member {res['id']:4d} x {res['x']:.6f} "
                f"ratio {res.get('ratio', float('nan')):.8f} N {res.get('n_nodes')} "
                f"-> {res['status']} frac {res.get('frac_reached', 0):.4f} "
                f"res {res.get('final_residual', float('nan')):.1e} "
                f"tail {res.get('final_tail', float('nan')):.1e} "
                f"({res.get('wall_s', 0):.0f}s) [{done}/{len(todo)}]",
                t_start,
            )
            if time.time() - t_start > budget_s:
                log("budget reached: terminating pool, rerun the stage to continue", t_start)
                pool.terminate()
                break


# ---------------------------------------------------------------------------
# Independent checks on converged circles.
# ---------------------------------------------------------------------------


def closure_check(
    circle_nodes: np.ndarray, rho: float, n_periods: int = 5, n_pick: int = 16
) -> dict[str, Any]:
    """Map nodes forward ``n_periods`` in ONE integration each with the core 6-state EOM.

    Uses :func:`cyclerfinder.core.ccr4bp.ccr4bp_eom` (not the batched planar RHS used by the
    corrector) and compares the state at each multiple of the forcing period with the
    interpolated circle rotated by ``k * rho``.
    """
    sysp = system(MU3)
    period = sc.forcing_period(sysp)
    n = circle_nodes.shape[0]
    idx = np.linspace(0, n - 1, n_pick).round().astype(int)
    th = sc.node_angles(n)
    errs = np.zeros((n_pick, n_periods))
    for a, j in enumerate(idx):
        sol = solve_ivp(
            ccr4bp.ccr4bp_eom,
            (0.0, n_periods * period),
            sc.to_state6(circle_nodes[j]),
            args=(sysp,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
            t_eval=period * np.arange(1, n_periods + 1),
        )
        for k in range(n_periods):
            pred = sc.fourier_eval(circle_nodes, th[j] + (k + 1) * rho)
            errs[a, k] = float(np.max(np.abs(sc.to_state4(sol.y[:, k]) - pred)))
    return {"max_err_per_period": errs.max(axis=0).tolist(), "nodes_checked": idx.tolist()}


def _saturn_centred_rhs(t: float, yflat: np.ndarray, m: int) -> np.ndarray:
    """Kumar Eq. [4] in velocity form: Rhea on a circle about Saturn at (-mu, 0)."""
    yy = yflat.reshape(4, m)
    x, y, vx, vy = yy
    dx1, dx2 = x + MU, x - 1.0 + MU
    r1 = np.hypot(dx1, y) ** 3
    r2 = np.hypot(dx2, y) ** 3
    th = OMEGA * t
    x3 = -MU + A3 * math.cos(th)
    y3 = A3 * math.sin(th)
    d3 = np.hypot(x - x3, y - y3) ** 3
    ax = x - (1 - MU) * dx1 / r1 - MU * dx2 / r2 - MU3 * (x - x3) / d3 - MU3 * math.cos(th) / A3**2
    ay = y - (1 - MU) * y / r1 - MU * y / r2 - MU3 * (y - y3) / d3 - MU3 * math.sin(th) / A3**2
    return np.concatenate([vx, vy, ax + 2 * vy, ay - 2 * vx])


def saturn_centred_residual(circle_nodes: np.ndarray, rho: float) -> float:
    """Invariance residual of a converged circle under the paper's Saturn-centred perturber."""
    n = circle_nodes.shape[0]
    sol = solve_ivp(
        _saturn_centred_rhs,
        (0.0, TP),
        circle_nodes.T.reshape(-1).copy(),
        args=(n,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    fu = sol.y[:, -1].reshape(4, n).T
    g = sc.shift_matrix(n, -rho) @ fu - circle_nodes
    return float(np.max(np.abs(g)))


def refine_distance(
    a_nodes: np.ndarray, b_nodes: np.ndarray, rho: float, n_pick: int = 12
) -> float:
    """Max over sampled nodes of ``b`` of the distance to the curve ``a``.

    Newton on ``<u(theta) - p, u'(theta)> = 0`` from the nearest point of a dense sample. The
    module's ``distance_to_circle`` minimises ``|u(theta) - p|`` with a bounded Brent search,
    whose relative tolerance in theta (about 1.5e-8 * theta) floors the distance at about
    1e-8, too coarse for this check.
    """
    del rho
    m = 16 * a_nodes.shape[0]
    th_dense = sc.node_angles(m)
    dense = sc.fourier_eval(a_nodes, th_dense)
    idx = np.linspace(0, b_nodes.shape[0] - 1, n_pick).round().astype(int)
    worst = 0.0
    for j in idx:
        p = b_nodes[j]
        th = float(th_dense[int(np.argmin(np.linalg.norm(dense - p[None, :], axis=1)))])
        for _ in range(30):
            d = sc.fourier_eval(a_nodes, th) - p
            d1 = sc.fourier_eval(a_nodes, th, deriv=1)
            d2 = sc.fourier_eval(a_nodes, th, deriv=2)
            step = float(np.dot(d, d1) / (np.dot(d1, d1) + np.dot(d, d2)))
            th -= step
            if abs(step) < 1e-15:
                break
        worst = max(worst, float(np.linalg.norm(sc.fourier_eval(a_nodes, th) - p)))
    return worst


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------


def _load(name: str) -> Any:
    return json.loads((OUT_DIR / name).read_text())


def _save(name: str, obj: Any) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / name).write_text(json.dumps(obj, indent=1))


def _read_jsonl(name: str) -> list[dict[str, Any]]:
    p = OUT_DIR / name
    if not p.exists():
        return []
    return [json.loads(line) for line in p.read_text().splitlines() if line.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage",
        required=True,
        choices=[
            "constants",
            "family",
            "members",
            "scan",
            "control",
            "refine",
            "posthoc",
            "closure",
            "summary",
        ],
    )
    ap.add_argument("--budget-s", type=float, default=420.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--posthoc-factor", type=float, default=3.0)
    ap.add_argument("--fresh", action="store_true", help="family: recompute the rows")
    args = ap.parse_args()
    t_start = time.time()
    workers = min(4, args.workers)
    if args.stage == "constants":
        _save("constants.json", stage_constants(t_start))
    elif args.stage == "family":
        fpath = OUT_DIR / "family.json"
        reuse = _load("family.json") if fpath.exists() and not args.fresh else None
        _save("family.json", stage_family(t_start, reuse))
    elif args.stage == "members":
        _save("members.json", stage_members(_load("family.json"), workers, t_start))
    elif args.stage == "scan":
        members = _load("members.json")
        done = _done_ids(OUT_DIR / "scan.jsonl")
        todo = [m for m in members if m["id"] not in done]
        log(f"scan: {len(todo)} of {len(members)} members left", t_start)
        run_pool(
            todo, _load("family.json"), {}, OUT_DIR / "scan.jsonl", args.budget_s, workers, t_start
        )
    elif args.stage == "control":
        stage_control(args.budget_s, workers, t_start)
    elif args.stage == "posthoc":
        stage_posthoc(args.budget_s, workers, t_start, args.posthoc_factor)
    elif args.stage == "refine":
        stage_refine(args.budget_s, workers, t_start)
    elif args.stage == "closure":
        stage_closure(args.budget_s, t_start)
    elif args.stage == "summary":
        _save("summary.json", stage_summary(t_start))


CONTROL_CROSSING = (8, 41, 1)  # p, q, segment: a low-e crossing inside Region A
CONTROL_FACTORS = (100.0, 1000.0)


def stage_control(budget_s: float, workers: int, t_start: float) -> None:
    """Method check (not a paper check): the same crossing with Rhea's mass raised 100x, 1000x.

    Resonance widths grow like sqrt(mu3), so at a raised mass the detector must report failure
    near the ratio and persistence away from it.
    """
    fam = Family(_load("family.json"))
    p, q, si = CONTROL_CROSSING
    seg = fam.segments()[si]
    path = OUT_DIR / "control.jsonl"
    done = _done_ids(path)
    todo = []
    k = 0
    for fac in CONTROL_FACTORS:
        for d in OFFSETS:
            for sgn in (1.0, -1.0) if d > 0 else (1.0,):
                k += 1
                if k in done:
                    continue
                xm = fam.solve_ratio(p / q + sgn * d, seg)
                if xm is None:
                    continue
                todo.append(
                    {
                        "id": k,
                        "kind": "control",
                        "p": p,
                        "q": q,
                        "segment": si,
                        "offset": sgn * d,
                        "x": xm,
                        "factor": fac,
                        "target": fac * MU3,
                    }
                )
    log(f"control: {len(todo)} members left", t_start)
    run_pool(todo, _load("family.json"), {}, path, budget_s, workers, t_start)


def stage_refine(budget_s: float, workers: int, t_start: float) -> None:
    """Every non-persisting member, rerun with more nodes and a quarter of the step."""
    scan = {r["id"]: r for r in _read_jsonl("scan.jsonl")}
    done = _done_ids(OUT_DIR / "refine.jsonl")
    todo = []
    for r in scan.values():
        if r["status"] == "persist" or r["id"] in done:
            continue
        if r["status"] in ("seed_unresolved", "crash"):
            continue
        levels = (next_odd(1.5 * r["n_nodes"]),)
        m = dict(r["member"])
        todo.append((m, levels))
    log(f"refine: {len(todo)} members left", t_start)
    fam_data = _load("family.json")
    by_levels: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    for m, lv in todo:
        by_levels.setdefault(lv, []).append(m)
    for lv, ms in by_levels.items():
        kw = {"levels": lv, "h0": H0 / 4, "h_floor": H_FLOOR / 4}
        run_pool(ms, fam_data, kw, OUT_DIR / "refine.jsonl", budget_s, workers, t_start)
        if time.time() - t_start > budget_s:
            return


def stage_closure(budget_s: float, t_start: float) -> None:
    """Independent closure, Saturn-centred residual and node refinement for sampled successes."""
    scan = [r for r in _read_jsonl("scan.jsonl") if r["status"] == "persist"]
    done = _done_ids(OUT_DIR / "closure.jsonl")
    scan.sort(key=lambda r: r["x"])
    pick = scan[:: max(1, len(scan) // 12)]
    fam = Family(_load("family.json"))
    for r in pick:
        if r["id"] in done:
            continue
        nodes = np.asarray(r["final_nodes"])
        rho = (2.0 * math.pi * TP / r["T"]) % (2.0 * math.pi)
        out: dict[str, Any] = {
            "id": r["id"],
            "x": r["x"],
            "ratio": r["ratio"],
            "n_nodes": r["n_nodes"],
        }
        out["closure"] = closure_check(nodes, rho)
        out["saturn_centred_residual"] = saturn_centred_residual(nodes, rho)
        n2 = next_odd(1.5 * r["n_nodes"])
        if n2 > 701:  # a dense solve at more than 701 nodes exceeds the per-command budget
            ref = {"status": "skipped (N too large for the dense solve)"}
        else:
            ref = continue_member(r["x"], fam, levels=(n2,), with_bundles=False)
        out["refined_status"] = ref["status"]
        out["refined_n"] = n2
        if ref["status"] == "persist":
            out["refined_distance"] = refine_distance(nodes, np.asarray(ref["final_nodes"]), rho)
        with (OUT_DIR / "closure.jsonl").open("a") as fh:
            fh.write(json.dumps(out) + "\n")
        log(
            f"closure id {r['id']} x {r['x']:.5f}: errs "
            f"{[f'{e:.1e}' for e in out['closure']['max_err_per_period']]} "
            f"saturn-centred res {out['saturn_centred_residual']:.1e} refined {ref['status']} "
            f"dist {out.get('refined_distance', float('nan')):.1e}",
            t_start,
        )
        if time.time() - t_start > budget_s:
            return


def stage_posthoc(budget_s: float, workers: int, t_start: float, factor: float = 3.0) -> None:
    """NOT pre-registered: members not PERSIST in scan and refine that DID reach the target mass
    (off-node or tail test failed), rerun at three times the scan's node count (cap 1001)."""
    scan = {r["id"]: r for r in _read_jsonl("scan.jsonl")}
    ref = {r["id"]: r for r in _read_jsonl("refine.jsonl")}
    done = _done_ids(OUT_DIR / "posthoc.jsonl")
    by_n: dict[int, list[dict[str, Any]]] = {}
    for i, r in ref.items():
        if i in done or r["status"] != "unresolved":
            continue
        n3 = next_odd(factor * scan[i]["n_nodes"])
        if n3 > 1001:
            continue
        by_n.setdefault(n3, []).append(dict(scan[i]["member"]))
    log(f"posthoc: {sum(len(v) for v in by_n.values())} members left", t_start)
    for n3, ms in sorted(by_n.items()):
        kw = {"levels": (n3,), "h0": H0 / 4, "h_floor": H_FLOOR / 4}
        run_pool(
            ms, _load("family.json"), kw, OUT_DIR / "posthoc.jsonl", budget_s, workers, t_start
        )
        if time.time() - t_start > budget_s:
            return


def classify_member(r_scan: dict[str, Any], r_ref: dict[str, Any] | None) -> str:
    """Pre-registered classes: persist / rescued / fail-floor / fail-unresolved / region-B."""
    if r_scan["status"] == "seed_unresolved":
        return "region-B"
    if r_scan["status"] == "persist":
        return "persist"
    if r_ref is None:
        return "crash" if r_scan["status"] == "crash" else "not-refined"
    if r_ref["status"] == "persist":
        return "rescued"
    if r_ref["status"] == "unresolved":
        return "fail-unresolved"
    return "fail-floor"


def locate(ratio: float, extrema: list[float]) -> dict[str, Any]:
    near = min(LISTED, key=lambda pq: abs(ratio - pq[0] / pq[1]))
    d = ratio - near[0] / near[1]
    d_ext = min(abs(ratio - e) for e in extrema)
    where = "window" if abs(d) <= WINDOW else ("twist-band" if d_ext <= TWIST_BAND else "outside")
    p100, q100, d100 = nearest_rational(ratio, 100)
    return {
        "nearest_listed": f"{near[0]}/{near[1]}",
        "delta_listed": d,
        "delta_extremum": d_ext,
        "where": where,
        "nearest_q100": f"{p100}/{q100}",
        "delta_q100": d100,
    }


def stage_summary(t_start: float) -> dict[str, Any]:
    fam = _load("family.json")
    extrema = [e["ratio"] for e in fam["extrema"]]
    scan = {r["id"]: r for r in _read_jsonl("scan.jsonl")}
    ref = {r["id"]: r for r in _read_jsonl("refine.jsonl")}
    post = {r["id"]: r for r in _read_jsonl("posthoc.jsonl")}
    table = []
    for i, r in sorted(scan.items(), key=lambda kv: kv[1]["x"]):
        m = r["member"]
        cls = classify_member(r, ref.get(i))
        row: dict[str, Any] = {
            "id": i,
            "kind": m["kind"],
            "segment": m.get("segment"),
            "offset": m.get("offset"),
            "x": r["x"],
            "C": r.get("C"),
            "ratio": r.get("ratio"),
            "class": cls,
            "n_scan": r.get("n_nodes"),
            "frac_scan": r.get("frac_reached"),
            "offnode_scan": r.get("offnode_residual"),
            "n_refine": ref.get(i, {}).get("n_nodes"),
            "frac_refine": ref.get(i, {}).get("frac_reached"),
            "offnode_refine": ref.get(i, {}).get("offnode_residual"),
            "tail_refine": ref.get(i, {}).get("final_tail"),
        }
        best = r if cls == "persist" else ref.get(i, r)
        row["final_residual"] = best.get("final_residual")
        row["final_tail"] = best.get("final_tail")
        row["deformation_max"] = best.get("deformation_max")
        row["deformation_harmonic_q"] = best.get("deformation_harmonic_q")
        row["harmonic_q"] = best.get("harmonic_q")
        row["lam_u_map"] = best.get("lam_u_map")
        if r.get("ratio") is not None:
            row.update(locate(r["ratio"], extrema))
        if cls == "region-B":
            row["n_needed_estimate"] = r.get("n_needed_estimate")
        if cls.startswith("fail"):
            hist = (ref.get(i) or r).get("last_fail") or {}
            h = hist.get("history") or []
            row["tolerance_edge"] = bool(h) and min(h) < 2.0 * TOL
        if i in post:
            row["posthoc_status"] = post[i]["status"]
            row["posthoc_n"] = post[i].get("n_nodes")
            row["posthoc_offnode"] = post[i].get("offnode_residual")
        table.append(row)
    fam_rows = {r["x"]: r for r in fam["rows"]}
    del fam_rows
    tested = [t for t in table if t["class"] != "region-B"]
    fails = [t for t in tested if t["class"].startswith("fail")]
    summary = {
        "counts": {
            c: sum(1 for t in table if t["class"] == c) for c in sorted({t["class"] for t in table})
        },
        "fails_by_where": {
            w: sum(1 for t in fails if t["where"] == w) for w in ("window", "twist-band", "outside")
        },
        "fails_outside": [t for t in fails if t["where"] == "outside"],
        "max_abs_delta_of_window_fail": {},
        "region_b": [
            {k: t.get(k) for k in ("id", "x", "C", "ratio", "kind", "offset", "n_needed_estimate")}
            for t in table
            if t["class"] == "region-B"
        ],
        "table": table,
    }
    for p, q in LISTED:
        tag = f"{p}/{q}"
        per_seg: dict[str, Any] = {}
        for t in tested:
            if t["nearest_listed"] != tag or t["kind"] != "crossing":
                continue
            seg = str(t["segment"])
            d = per_seg.setdefault(seg, {"fail_offsets": [], "pass_offsets": []})
            key = "fail_offsets" if t["class"].startswith("fail") else "pass_offsets"
            d[key].append(t["offset"])
        for d in per_seg.values():
            d["fail_offsets"].sort()
            d["pass_offsets"].sort()
            d["max_abs_fail_offset"] = max((abs(o) for o in d["fail_offsets"]), default=None)
        summary["max_abs_delta_of_window_fail"][tag] = per_seg
    control = _read_jsonl("control.jsonl")
    summary["control"] = sorted(
        (
            {
                "factor": c["member"]["factor"],
                "offset": c["member"]["offset"],
                "status": c["status"],
                "frac": c.get("frac_reached"),
                "offnode": c.get("offnode_residual"),
                "tail": c.get("final_tail"),
                "deformation_harmonic_q": c.get("deformation_harmonic_q"),
            }
            for c in control
        ),
        key=lambda c: (c["factor"], c["offset"]),
    )
    summary["closure"] = _read_jsonl("closure.jsonl")
    log(f"classes: {summary['counts']}", t_start)
    log(f"fails by location: {summary['fails_by_where']}", t_start)
    for tag, segs in summary["max_abs_delta_of_window_fail"].items():
        log(f"{tag}: {segs}", t_start)
    for t in summary["fails_outside"]:
        log(
            f"FAIL OUTSIDE: id {t['id']} x {t['x']:.5f} ratio {t['ratio']:.7f} {t['class']} "
            f"nearest q<=100 {t['nearest_q100']} ({t['delta_q100']:.1e})",
            t_start,
        )
    return summary


if __name__ == "__main__":
    main()
