"""#1000: Earth-Moon complement search. Periodic orbits as fixed points of the perigee return map
at fixed Jacobi constant, symmetric or not; the asymmetric ones are outside the Restrepo-Russell
and Franz-Russell databases by construction.

Pre-registration: docs/notes/2026-10-08-1000-earth-moon-complement-search.md sec. 0 (and its
amendments). Every call runs in the foreground under --max-seconds and resumes from checkpoints.

    uv run python scripts/run_1000_complement.py pilot            # 50-seed timing pilot
    uv run python scripts/run_1000_complement.py scan             # grid scan, resumable
    uv run python scripts/run_1000_complement.py refine           # Newton on candidates
"""

from __future__ import annotations

import argparse
import json
import math
import time
from multiprocessing import Pool
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
OUT = REPO / "data" / "1000_complement"
SYS = cr3bp_system("Earth", "Moon")
MU = SYS.mu
L_KM = SYS.l_km
R_E = PLANETS["E"].radius_eq_km
R_M = SATELLITES["Moon"].radius_eq_km
FLOOR_E = R_E + PLANETS["E"].safe_alt_km
FLOOR_M = R_M + SATELLITES["Moon"].safe_alt_km
GEO_KM = 42164.0
HILL_KM = (MU / 3.0) ** (1.0 / 3.0) * L_KM
RTOL = ATOL = 1e-12

# pre-registered grid (note sec. 0.3) + the two control levels (sec. 0.6)
C_LEVELS = [round(0.5 + 0.1 * i, 4) for i in range(27)] + [1.0687, 1.0672]
RP_KM = list(np.geomspace(FLOOR_E, GEO_KM, 8))
OMEGA = [math.radians(10.0 * i) for i in range(36)]
SENSES = (+1, -1)
K_MAX = 7
T_MAX = 30.0  # TU; 7 returns of the 7-3 class take 18.85 TU
RES_MAX = 0.3


def eom(t: float, s: Arr, mu: float) -> Arr:
    x, y, vx, vy = s
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    c1 = (1.0 - mu) / r1**3
    c2 = mu / r2**3
    return np.array(
        [
            vx,
            vy,
            x - c1 * (x + mu) - c2 * (x - 1.0 + mu) + 2.0 * vy,
            y - c1 * y - c2 * y - 2.0 * vx,
        ]
    )


def omega_eff(x: float, y: float) -> float:
    r1 = math.hypot(x + MU, y)
    r2 = math.hypot(x - 1.0 + MU, y)
    return x * x + y * y + 2 * (1 - MU) / r1 + 2 * MU / r2


def section_state(c: float, rp: float, om: float, sense: int) -> Arr | None:
    """Rotating-frame state at a perigee (r_p nondimensional, omega from +x about the Earth)."""
    x = -MU + rp * math.cos(om)
    y = rp * math.sin(om)
    v2 = omega_eff(x, y) - c
    if v2 <= 0:
        return None
    v = math.sqrt(v2)
    return np.array([x, y, -sense * v * math.sin(om), sense * v * math.cos(om)])


def to_section(s: Arr) -> tuple[float, float]:
    dx, dy = s[0] + MU, s[1]
    return math.hypot(dx, dy), math.atan2(dy, dx) % (2 * math.pi)


def _dr1(t: float, s: Arr, mu: float) -> float:
    return float((s[0] + mu) * s[2] + s[1] * s[3])


def _dr2(t: float, s: Arr, mu: float) -> float:
    return float((s[0] - 1 + mu) * s[2] + s[1] * s[3])


def _hit_e(t: float, s: Arr, mu: float) -> float:
    return math.hypot(s[0] + mu, s[1]) - R_E / L_KM


def _hit_m(t: float, s: Arr, mu: float) -> float:
    return math.hypot(s[0] - 1 + mu, s[1]) - R_M / L_KM


_dr1.direction = 1.0  # type: ignore[attr-defined]  (perigee: dr/dt from - to +)
_dr2.direction = 1.0  # type: ignore[attr-defined]
_hit_e.terminal = True  # type: ignore[attr-defined]
_hit_m.terminal = True  # type: ignore[attr-defined]


def returns(s0: Arr, k_max: int = K_MAX, t_max: float = T_MAX) -> dict[str, Any]:
    """Integrate to the first k_max perigees; per return (t, r_p, omega, min r2 so far)."""
    sol = solve_ivp(
        eom,
        (0.0, t_max),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=[_dr1, _dr2, _hit_e, _hit_m],
    )
    out: list[list[float]] = []
    t_r2 = sol.t_events[1]
    r2s = [math.hypot(y[0] - 1 + MU, y[1]) for y in sol.y_events[1]]
    for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True):
        if t < 1e-6:
            continue
        rp, om = to_section(y)
        prior = [r for tt, r in zip(t_r2, r2s, strict=True) if tt <= t]
        out.append([float(t), rp, om, min(prior) if prior else float("inf")])
        if len(out) >= k_max:
            break
    status = "ok"
    if sol.status == 1 and (len(sol.t_events[2]) or len(sol.t_events[3])):
        status = "impact"
    return {"returns": out, "status": status}


def residual(rp0: float, om0: float, rp1: float, om1: float) -> float:
    d_om = (om1 - om0 + math.pi) % (2 * math.pi) - math.pi
    return math.hypot(math.log(rp1 / rp0), d_om)


def scan_one(args: tuple[float, int, int, int]) -> dict[str, Any]:
    c, sense, ir, iw = args
    rp = RP_KM[ir] / L_KM
    om = OMEGA[iw]
    rec: dict[str, Any] = {"C": c, "sense": sense, "ir": ir, "iw": iw}
    s0 = section_state(c, rp, om, sense)
    if s0 is None:
        rec["status"] = "forbidden"
        return rec
    # require a perigee (r1'' > 0) at the start
    a = eom(0.0, s0, MU)
    rdd = (s0[2] ** 2 + s0[3] ** 2) + (s0[0] + MU) * a[2] + s0[1] * a[3]
    if rdd <= 0:
        rec["status"] = "apogee"
        return rec
    r = returns(s0)
    rec["status"] = r["status"]
    rec["res"] = [residual(rp, om, x[1], x[2]) for x in r["returns"]]
    rec["t"] = [x[0] for x in r["returns"]]
    rec["minr2_km"] = [x[3] * L_KM for x in r["returns"]]
    return rec


def scan_file(c: float, sense: int) -> Path:
    return OUT / "scan" / f"C{c:.4f}_s{'p' if sense > 0 else 'm'}.jsonl"


def cmd_scan(max_seconds: float, levels: list[float] | None, pilot: int | None) -> None:
    (OUT / "scan").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_done = 0
    with Pool(2) as pool:
        for c in levels or C_LEVELS:
            for sense in SENSES:
                f = scan_file(c, sense)
                done = set()
                if f.exists():
                    for line in f.read_text().splitlines():
                        d = json.loads(line)
                        done.add((d["ir"], d["iw"]))
                todo = [
                    (c, sense, ir, iw)
                    for ir in range(len(RP_KM))
                    for iw in range(len(OMEGA))
                    if (ir, iw) not in done
                ]
                if pilot is not None:
                    todo = todo[:pilot]
                if not todo:
                    continue
                with f.open("a") as fh:
                    for rec in pool.imap_unordered(scan_one, todo, chunksize=2):
                        fh.write(json.dumps(rec) + "\n")
                        fh.flush()
                        n_done += 1
                        if time.time() - t0 > max_seconds:
                            print(f"time budget reached after {n_done} seeds", flush=True)
                            pool.terminate()
                            return
                el = time.time() - t0
                print(
                    f"{time.strftime('%H:%M:%S')} C={c:.4f} sense={sense:+d}: {len(todo)} seeds; "
                    f"total {n_done} in {el:.0f} s ({el / max(n_done, 1):.2f} s/seed)",
                    flush=True,
                )
                if pilot is not None:
                    return
    print(f"scan complete: {n_done} seeds this call in {time.time() - t0:.0f} s")


# ------------------------------------------------------------------ candidates and Newton

F_TOL, CLOSE_TOL, PERP_TOL = 1e-8, 1e-7, 1e-6  # note sec. 0.9 amendment 2


def candidates(f: Path) -> list[dict[str, Any]]:
    rs = [json.loads(line) for line in f.read_text().splitlines()]
    g = {(r["ir"], r["iw"]): r for r in rs}
    out = []
    for k in range(1, K_MAX + 1):
        grid = np.full((len(RP_KM), len(OMEGA)), np.inf)
        for (ir, iw), r in g.items():
            if len(r.get("res", [])) >= k and r["minr2_km"][k - 1] <= HILL_KM:
                grid[ir, iw] = r["res"][k - 1]
        for ir in range(len(RP_KM)):
            for iw in range(len(OMEGA)):
                v = grid[ir, iw]
                if not v < RES_MAX:
                    continue
                nb = [
                    grid[i, (iw + j) % len(OMEGA)]
                    for i in (ir - 1, ir, ir + 1)
                    for j in (-1, 0, 1)
                    if 0 <= i < len(RP_KM) and not (i == ir and j == 0)
                ]
                if all(v <= x for x in nb):
                    r = g[(ir, iw)]
                    out.append(
                        {
                            **{q: r[q] for q in ("C", "sense", "ir", "iw")},
                            "k": k,
                            "res0": v,
                            "t_k": r["t"][k - 1],
                        }
                    )
    return out


def pk(c: float, sense: int, z: Arr, k: int, t_max: float) -> tuple[Arr, float] | None:
    rp, om = math.exp(z[0]), z[1]
    s0 = section_state(c, rp, om, sense)
    if s0 is None:
        return None
    r = returns(s0, k_max=k, t_max=t_max)
    if len(r["returns"]) < k:
        return None
    t, rp1, om1, _ = r["returns"][k - 1]
    d_om = (om1 - om + math.pi) % (2 * math.pi) - math.pi
    return np.array([math.log(rp1 / rp), d_om]), t


def newton(c: float, sense: int, z0: Arr, k: int, t_k: float) -> dict[str, Any]:
    z = z0.copy()
    t_max = 1.6 * t_k + 1.0
    h = 1e-7
    for it in range(1, 16):
        f0 = pk(c, sense, z, k, t_max)
        if f0 is None:
            return {"converged": False, "why": "no k-th return", "it": it}
        fz, t = f0
        if float(np.linalg.norm(fz)) < F_TOL:
            return {
                "converged": True,
                "z": z.tolist(),
                "T": t,
                "it": it,
                "F": float(np.linalg.norm(fz)),
            }
        jac = np.zeros((2, 2))
        for j in range(2):
            e = np.zeros(2)
            e[j] = h
            fp, fm = pk(c, sense, z + e, k, t_max), pk(c, sense, z - e, k, t_max)
            if fp is None or fm is None:
                return {"converged": False, "why": "FD eval failed", "it": it}
            jac[:, j] = (fp[0] - fm[0]) / (2 * h)
        a = jac - np.eye(2)
        try:
            dz = np.linalg.solve(a, -fz)
        except np.linalg.LinAlgError:
            return {"converged": False, "why": "singular", "it": it}
        nrm = float(np.linalg.norm(dz))
        if nrm > 0.2:
            dz *= 0.2 / nrm
        z = z + dz
        z[1] %= 2 * math.pi
    return {"converged": False, "why": "max iterations", "it": 15}


def refine_one(cand: dict[str, Any]) -> dict[str, Any]:
    c, sense, k = cand["C"], cand["sense"], cand["k"]
    z0 = np.array([math.log(RP_KM[cand["ir"]] / L_KM), OMEGA[cand["iw"]]])
    res = newton(c, sense, z0, k, cand["t_k"])
    out = {**cand, **res}
    if not res["converged"]:
        return out
    z = np.array(res["z"])
    # minimal period: smallest divisor d of k with P^d(z) = z
    k_min = k
    for d in range(1, k):
        if k % d:
            continue
        f = pk(c, sense, z, d, 1.6 * res["T"] + 1.0)
        if f is not None and float(np.linalg.norm(f[0])) < 1e-6:
            k_min = d
            out["T"] = f[1]
            break
    out["k_min"] = k_min
    s0 = section_state(c, math.exp(z[0]), z[1], sense)
    assert s0 is not None
    out["state0"] = s0.tolist()
    out.update(classify(s0, out["T"]))
    return out


def classify(s0: Arr, period: float) -> dict[str, Any]:
    def ey(t: float, s: Arr, mu: float) -> float:
        return float(s[1])

    sol = solve_ivp(
        eom,
        (0.0, period),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=[ey, _dr1, _dr2],
        dense_output=True,
    )
    perp = [
        [float(t), float(y[0]), float(y[3])]
        for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True)
        if abs(y[2]) < PERP_TOL
    ]

    # all r1 and r2 extrema (both directions) for minima and maxima
    def extrema(fun: Any) -> list[float]:
        fun.direction = 0.0
        s2 = solve_ivp(
            eom, (0.0, period), s0, args=(MU,), method="DOP853", rtol=RTOL, atol=ATOL, events=fun
        )
        fun.direction = 1.0
        return [float(t) for t in s2.t_events[0]], s2.y_events[0]

    _, y1 = extrema(_dr1)
    _, y2 = extrema(_dr2)
    r1 = [math.hypot(y[0] + MU, y[1]) * L_KM for y in y1] + [math.hypot(s0[0] + MU, s0[1]) * L_KM]
    r2 = [math.hypot(y[0] - 1 + MU, y[1]) * L_KM for y in y2]
    tg = np.linspace(0.0, period, 20001)
    zz = sol.sol(tg)
    w1 = (np.unwrap(np.arctan2(zz[1], zz[0] + MU))[-1] - math.atan2(s0[1], s0[0] + MU)) / (
        2 * math.pi
    )
    w2 = (np.unwrap(np.arctan2(zz[1], zz[0] - 1 + MU))[-1] - math.atan2(s0[1], s0[0] - 1 + MU)) / (
        2 * math.pi
    )
    from cyclerfinder.core.cr3bp import propagate

    s6 = np.array([s0[0], s0[1], 0.0, s0[2], s0[3], 0.0])
    arc = propagate(SYS, s6, period, with_stm=True, stm_mode="fixed_path")
    m = np.asarray(arc.stm)
    m4 = m[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])]
    mz = m[np.ix_([2, 5], [2, 5])]
    eig = np.linalg.eigvals(m4)
    order = np.argsort(np.abs(eig - 1.0))
    lam = eig[order[2:]]
    b = float(np.real(lam[0] + 1.0 / lam[0]))
    c = float(omega_eff(s0[0], s0[1]) - s0[2] ** 2 - s0[3] ** 2)
    pe, ps = min(r1), min(r2) if r2 else float("inf")
    return {
        "C_check": c,
        "closure": float(np.linalg.norm(sol.y[:, -1] - s0)),
        "symmetric": bool(perp),
        "perp_crossings": perp,
        "perigee_km": pe,
        "perigee_alt_km": pe - R_E,
        "apogee_km": max(r1),
        "periselene_km": ps,
        "periselene_alt_km": ps - R_M,
        "max_moon_km": max(r2) if r2 else float("nan"),
        "b_h": b,
        "lambda_max": float(np.max(np.abs(eig))),
        "k_perp": float(np.trace(mz)),
        "det_M4": float(np.linalg.det(m4)),
        "wind_E": float(w1),
        "wind_M": float(w2),
        "T_days": period * SYS.t_s / 86400.0,
        "cycler_class_candidate": bool(FLOOR_E <= pe <= GEO_KM and FLOOR_M <= ps <= HILL_KM),
        "passes_with_50km_moon_floor": bool(FLOOR_E <= pe <= GEO_KM and R_M + 50 <= ps <= HILL_KM),
    }


def cmd_refine(max_seconds: float) -> None:
    out_f = OUT / "refined.jsonl"
    done = set()
    if out_f.exists():
        for line in out_f.read_text().splitlines():
            d = json.loads(line)
            done.add((d["C"], d["sense"], d["ir"], d["iw"], d["k"]))
    todo = []
    # control levels (note sec. 0.6) first, so a failed control stops the run early
    files = sorted((OUT / "scan").glob("*.jsonl"), key=lambda f: ("1.06" not in f.name, f.name))
    for f in files:
        for cd in candidates(f):
            if (cd["C"], cd["sense"], cd["ir"], cd["iw"], cd["k"]) not in done:
                todo.append(cd)
    print(f"{len(done)} refined, {len(todo)} to go", flush=True)
    t0 = time.time()
    n = 0
    with Pool(2) as pool, out_f.open("a") as fh:
        for rec in pool.imap_unordered(refine_one, todo):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n += 1
            if n % 10 == 0:
                print(
                    f"{time.strftime('%H:%M:%S')} {n} refined ({time.time() - t0:.0f} s)",
                    flush=True,
                )
            if time.time() - t0 > max_seconds:
                print(f"time budget reached after {n}", flush=True)
                pool.terminate()
                return
    print(f"refine complete: {n} this call in {time.time() - t0:.0f} s")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("pilot")
    p1.add_argument("--n", type=int, default=50)
    p2 = sub.add_parser("scan")
    p2.add_argument("--max-seconds", type=float, default=420.0)
    p2.add_argument("--levels", type=float, nargs="*", default=None)
    p2.add_argument("--pilot-s", type=float, required=True, help="measured s/seed from pilot")
    p3 = sub.add_parser("refine")
    p3.add_argument("--max-seconds", type=float, default=420.0)
    p3.add_argument("--pilot-s", type=float, required=True, help="measured s/candidate")
    args = ap.parse_args()
    preflight_search(
        task_no=1000,
        region_id=f"em-complement-perigee-map-{args.cmd}",
        method=MethodCapability(
            genome="planar CR3BP perigee return map fixed points (symmetric or not), Earth-Moon",
            corrector="Newton on P^k(z) - z, finite-difference Jacobian",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar", "asymmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=args.n
        if args.cmd == "pilot"
        else len(C_LEVELS) * len(RP_KM) * len(OMEGA) * len(SENSES),
        timing_pilot_seconds_per_point=getattr(args, "pilot_s", None),
    )
    print(
        f"mu={MU!r} floors E {FLOOR_E:.1f} km M {FLOOR_M:.1f} km Hill {HILL_KM:.0f} km", flush=True
    )
    if args.cmd == "pilot":
        cmd_scan(420.0, [1.1], args.n)
    elif args.cmd == "scan":
        cmd_scan(args.max_seconds, args.levels, None)
    elif args.cmd == "refine":
        cmd_refine(args.max_seconds)


if __name__ == "__main__":
    main()
