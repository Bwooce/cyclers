"""#1000 method (1): multiple shooting from near-Keplerian p:q skeletons (amendment A,
docs/notes/2026-10-08-1000-earth-moon-complement-search.md sec. 2.1).

    uv run python scripts/run_1000_shooting.py control   # 7:3 seeds at the 7-3b/7-3c C values
    uv run python scripts/run_1000_shooting.py grid      # full seed grid, C free

Every call runs in the foreground under --max-seconds and resumes from its checkpoint.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base

Arr = NDArray[np.float64]
REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1000_complement" / "shooting"
MU = base.MU
GM_E = 1.0 - MU
PQ = [(2, 1), (1, 2), (3, 2), (5, 2), (1, 3), (2, 3), (4, 3), (5, 3), (7, 3), (8, 3)]
RP = [r / base.L_KM for r in base.RP_KM]
OMEGA = base.OMEGA
CONTROLS = {
    "casoliva-7-3b-em-cycler-2010": 1.068655371747616,
    "casoliva-7-3c-em-cycler-2010": 1.0671969118897233,
}
MAX_IT, STEP_MAX, CONT_TOL, CLOSE_TOL = 25, 0.1, 1e-10, 1e-6


def eom_stm(t: float, z: Arr, mu: float) -> Arr:
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
    a = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [uxx, uxy, 0, 2.0], [uxy, uyy, -2.0, 0]])
    out = np.empty(20)
    out[:4] = base.eom(t, z[:4], mu)
    out[4:] = (a @ phi).ravel()
    return out


def arc(s: Arr, tau: float) -> tuple[Arr, Arr] | None:
    z0 = np.concatenate([s, np.eye(4).ravel()])
    sol = solve_ivp(
        eom_stm,
        (0.0, tau),
        z0,
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=[base._hit_e, base._hit_m],
    )
    if sol.status != 0:
        return None
    return sol.y[:4, -1], sol.y[4:, -1].reshape(4, 4)


def skeleton_nodes(p: int, q: int, sense: int, rp: float, om: float) -> tuple[list[Arr], float]:
    a = (q / p) ** (2.0 / 3.0) * GM_E ** (1.0 / 3.0)
    tp = 2.0 * math.pi * q / p
    vp = math.sqrt(GM_E * (2.0 / rp - 1.0 / a))
    ra = 2.0 * a - rp
    va = vp * rp / ra
    # apocentre of the inertial ellipse whose perigee points along omega at t = 0
    x_in = -ra * np.array([math.cos(om), math.sin(om)])
    v_in = -sense * va * np.array([-math.sin(om), math.cos(om)])
    nodes = []
    for i in range(p):
        # nodes at the skeleton's APOCENTRES (amendment A2): slow dynamics there, so the
        # period shift from the lunar perturbation costs little in the node mismatch
        t = (i + 0.5) * tp
        c, s = math.cos(-t), math.sin(-t)
        rot = np.array([[c, -s], [s, c]])
        xr = rot @ x_in
        vr = rot @ v_in + np.array([xr[1], -xr[0]])  # minus omega x r
        nodes.append(np.array([xr[0] - MU, xr[1], vr[0], vr[1]]))
    return nodes, tp


def jac_c(s: Arr) -> tuple[float, Arr]:
    x, y, vx, vy = s
    r1 = math.hypot(x + MU, y)
    r2 = math.hypot(x - 1 + MU, y)
    c = base.omega_eff(x, y) - vx * vx - vy * vy
    dcx = 2 * x - 2 * (1 - MU) * (x + MU) / r1**3 - 2 * MU * (x - 1 + MU) / r2**3
    dcy = 2 * y - 2 * (1 - MU) * y / r1**3 - 2 * MU * y / r2**3
    return c, np.array([dcx, dcy, -2 * vx, -2 * vy])


def shoot(nodes: list[Arr], taus: list[float], c_target: float | None) -> dict[str, Any]:
    n = len(nodes)
    xv = np.concatenate([np.concatenate(nodes), np.array(taus)])
    for it in range(1, MAX_IT + 1):
        s = [xv[4 * i : 4 * i + 4] for i in range(n)]
        tau = xv[4 * n :]
        if np.any(tau <= 0):
            return {"converged": False, "why": "negative duration", "it": it}
        m_rows = 4 * n + 1 + (1 if c_target is not None else 0)
        f = np.zeros(m_rows)
        jm = np.zeros((m_rows, 5 * n))
        for i in range(n):
            res = arc(s[i], float(tau[i]))
            if res is None:
                return {"converged": False, "why": "impact", "it": it}
            sf, phi = res
            j = (i + 1) % n
            f[4 * i : 4 * i + 4] = sf - s[j]
            jm[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] += phi
            jm[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] -= np.eye(4)
            jm[4 * i : 4 * i + 4, 4 * n + i] = base.eom(0.0, sf, MU)
        x0, y0, vx0, vy0 = s[0]
        f[4 * n] = (x0 + MU) * vx0 + y0 * vy0
        jm[4 * n, 0:4] = [vx0, vy0, x0 + MU, y0]
        if c_target is not None:
            c, g = jac_c(s[0])
            f[4 * n + 1] = c - c_target
            jm[4 * n + 1, 0:4] = g
        cont = float(np.max(np.abs(f[: 4 * n])))
        if (
            cont < CONT_TOL
            and abs(f[4 * n]) < CONT_TOL
            and (c_target is None or abs(f[4 * n + 1]) < CONT_TOL)
        ):
            return {
                "converged": True,
                "it": it,
                "s0": s[0].tolist(),
                "T": float(np.sum(tau)),
                "taus": tau.tolist(),
                "cont_res": cont,
            }
        dx = np.linalg.lstsq(jm, -f, rcond=None)[0]
        nrm = float(np.linalg.norm(dx[: 4 * n], ord=np.inf))
        if nrm > STEP_MAX:
            dx *= STEP_MAX / nrm
        xv = xv + dx
    return {"converged": False, "why": "max iterations", "it": MAX_IT}


def solve_seed(seed: tuple[int, int, int, int, int, float | None]) -> dict[str, Any]:
    p, q, sense, ir, iw, c_target = seed
    rec: dict[str, Any] = {"p": p, "q": q, "sense": sense, "ir": ir, "iw": iw, "C_target": c_target}
    rp = RP[ir]
    a = (q / p) ** (2.0 / 3.0) * GM_E ** (1.0 / 3.0)
    if not (rp < 1.0 < 2 * a - rp):
        rec["status"] = "no lunar-orbit crossing"
        return rec
    nodes, tp = skeleton_nodes(p, q, sense, rp, OMEGA[iw])
    res = shoot(nodes, [tp] * p, c_target)
    rec.update(res)
    rec["status"] = "converged" if res["converged"] else res["why"]
    if not res["converged"]:
        return rec
    s0 = np.array(res["s0"])
    cl = solve_ivp(
        base.eom, (0.0, res["T"]), s0, args=(MU,), method="DOP853", rtol=1e-12, atol=1e-12
    )
    rec["closure"] = float(np.linalg.norm(cl.y[:, -1] - s0))
    if rec["closure"] > CLOSE_TOL:
        rec["status"] = "closure fail"
        return rec
    # k-minimal: node 0 is an apogee; move to the first perigee and check P^d, d | p
    sp = solve_ivp(
        base.eom,
        (0.0, res["T"]),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=base._dr1,
    )
    if not len(sp.t_events[0]):
        rec["status"] = "no perigee"
        return rec
    s0 = np.asarray(sp.y_events[0][0])
    rec["s_perigee"] = s0.tolist()
    rpp, omp = base.to_section(s0)
    h = (s0[0] + MU) * s0[3] - s0[1] * s0[2]
    sheet = 1 if h > 0 else -1
    c_now = base.omega_eff(s0[0], s0[1]) - s0[2] ** 2 - s0[3] ** 2
    rec["C"] = c_now
    rec["sheet"] = sheet
    rec["z0"] = [math.log(rpp), omp]
    k_min = p
    for d in range(1, p):
        if p % d:
            continue
        f = base.pk(c_now, sheet, np.array([math.log(rpp), omp]), d, 1.6 * res["T"] + 1.0)
        if f is not None and float(np.linalg.norm(f[0])) < 1e-6:
            k_min = d
            rec["T"] = f[1]
            break
    rec["k_min"] = k_min
    rec.update({k: v for k, v in base.classify(s0, rec["T"]).items() if k != "perp_crossings"})
    return rec


def row_perigees(rid: str) -> tuple[list[tuple[float, float]], float]:
    rows = {r["id"]: r for r in yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())}
    c = rows[rid]["orbit_elements"]["cr3bp"]
    s6 = np.array(c["state_nd"])
    s = np.array([s6[0], s6[1], s6[3], s6[4]])
    r = base.returns(s, k_max=20, t_max=c["period_nd"] * 1.001)
    return [(math.log(x[1]), x[2]) for x in r["returns"]], float(c["period_nd"])


def run(tag: str, seeds: list[tuple[Any, ...]], max_seconds: float) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{tag}.jsonl"
    done = set()
    if f.exists():
        for line in f.read_text().splitlines():
            d = json.loads(line)
            done.add((d["p"], d["q"], d["sense"], d["ir"], d["iw"], d["C_target"]))
    todo = [s for s in seeds if tuple(s) not in done]
    print(f"{tag}: {len(done)} done, {len(todo)} to go", flush=True)
    t0, n = time.time(), 0
    with Pool(2) as pool, f.open("a") as fh:
        for rec in pool.imap_unordered(solve_seed, todo):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n += 1
            if n % 25 == 0:
                el = time.time() - t0
                print(f"{time.strftime('%H:%M:%S')} {n} ({el / n:.2f} s/seed)", flush=True)
            if time.time() - t0 > max_seconds:
                print(f"time budget reached after {n}", flush=True)
                pool.terminate()
                return
    print(f"{tag}: complete ({n} this call, {time.time() - t0:.0f} s)")


def control_report() -> None:
    f = OUT / "control.jsonl"
    recs = [json.loads(line) for line in f.read_text().splitlines()]
    conv = [r for r in recs if r["status"] == "converged"]
    print(f"control: {len(recs)} solved, {len(conv)} converged")
    for rid, c_t in CONTROLS.items():
        per, t_row = row_perigees(rid)
        hits = []
        for r in conv:
            if r["C_target"] != c_t:
                continue
            z = r["z0"]
            d = min(
                max(abs(z[0] - a), abs((z[1] - b + math.pi) % (2 * math.pi) - math.pi))
                for a, b in per
            )
            if d < 1e-6 and abs(r["T"] - t_row) / t_row < 1e-6:
                hits.append((r["sense"], r["ir"], r["iw"], r["symmetric"], d))
        print(f"{rid}: {len(hits)} seeds converged onto the row; first: {hits[:3]}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("control", "grid", "pilot"):
        sp = sub.add_parser(name)
        sp.add_argument("--max-seconds", type=float, default=420.0)
        sp.add_argument("--pilot-s", type=float, default=None)
    sub.add_parser("control-report")
    args = ap.parse_args()
    if args.cmd == "control-report":
        control_report()
        return
    control_seeds = [
        (7, 3, sense, ir, iw, c_t)
        for c_t in CONTROLS.values()
        for sense in (1, -1)
        for ir in range(len(RP))
        for iw in range(len(OMEGA))
    ]
    grid_seeds = [
        (p, q, sense, ir, iw, None)
        for p, q in PQ
        for sense in (1, -1)
        for ir in range(len(RP))
        for iw in range(len(OMEGA))
    ]
    seeds = {"pilot": control_seeds[::58][:20], "control": control_seeds, "grid": grid_seeds}[
        args.cmd
    ]
    preflight_search(
        task_no=1000,
        region_id=f"em-complement-multiple-shooting-{args.cmd}",
        method=MethodCapability(
            genome="near-Keplerian p:q skeleton, multiple shooting (N = p arcs), Earth-Moon CR3BP",
            corrector="minimum-norm Newton on continuity + perigee phase (+ C in controls)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar", "asymmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=len(seeds),
        timing_pilot_seconds_per_point=args.pilot_s,
    )
    run(args.cmd, seeds, args.max_seconds)


if __name__ == "__main__":
    main()
