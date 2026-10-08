"""#1000 method (2): symmetry-breaking at b = +2 pitchforks (amendment A, note sec. 2.2).

    uv run python scripts/run_1000_pitchfork.py control     # 7-3b/c branch -> parent -> branch
    uv run python scripts/run_1000_pitchfork.py parents997  # pitchforks along the #997 families
    uv run python scripts/run_1000_pitchfork.py parentsRR   # RR records with |b_h - 2| < 0.02

A parent is a symmetric periodic orbit with in-plane b = +2. Branching: perturb the parent's
start state along the monodromy eigenvector of eigenvalue +1 that is ANTISYMMETRIC under the
x-axis reflection (x, y, xdot, ydot) -> (x, -y, -xdot, ydot), and solve the multiple-shooting
problem at C_parent + dC on both sides and with both signs. A solution is a new branch if it is
asymmetric (no perpendicular crossing within 1e-6) and more than 1e-5 from the parent.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_997_lineage as lin
import run_1000_complement as base
import run_1000_continue as cn
import run_1000_shooting as sh

REPO = Path(__file__).resolve().parent.parent
GRID = REPO / "data" / "1000_complement" / "shooting" / "grid.jsonl"
OUT = REPO / "data" / "1000_complement" / "pitchfork"
MU = base.MU
REFL = np.diag([1.0, -1.0, -1.0, 1.0])


def min_perp(s0: np.ndarray, period: float) -> float:
    def ey(t: float, s: np.ndarray, mu: float) -> float:
        return float(s[1])

    sol = solve_ivp(
        base.eom,
        (0.0, period),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=ey,
    )
    return min((abs(float(y[2])) for y in sol.y_events[0]), default=float("inf"))


def ms_solution(s0: np.ndarray, period: float, n: int, c: float) -> dict[str, Any]:
    nodes, taus = cn.nodes_from(s0, period, n)
    return sh.shoot(nodes, taus, c)


def monodromy(res: dict[str, Any]) -> np.ndarray:
    m = np.eye(4)
    for s, t in zip(res["nodes"], res["taus"], strict=True):
        a = sh.arc(np.array(s), float(t))
        assert a is not None
        m = a[1] @ m
    return m


def rotate_to_perp(res: dict[str, Any]) -> tuple[list[np.ndarray], list[float]] | None:
    """Re-order a symmetric multiple-shooting solution so node 0 sits at a perpendicular
    x-axis crossing (found inside one arc; no long re-propagation)."""
    nodes = [np.array(x) for x in res["nodes"]]
    taus = list(res["taus"])

    def ey(t: float, s: np.ndarray, mu: float) -> float:
        return float(s[1])

    for i, (sn, tn) in enumerate(zip(nodes, taus, strict=True)):
        sol = solve_ivp(
            base.eom,
            (0.0, tn),
            sn,
            args=(MU,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
            events=ey,
        )
        for te, ye in zip(sol.t_events[0], sol.y_events[0], strict=True):
            if abs(ye[2]) < 1e-6 and 1e-6 < te < tn - 1e-6:
                n = len(nodes)
                new_nodes = [np.asarray(ye)] + [nodes[(i + k) % n] for k in range(1, n)]
                new_nodes.append(nodes[i])
                new_taus = [tn - te] + [taus[(i + k) % n] for k in range(1, n)] + [te]
                # the appended node i duplicates nothing: arcs are crossing -> i+1 -> ... -> i
                # -> crossing
                return new_nodes, new_taus
    return None


def shoot_eta(nodes: list[np.ndarray], taus: list[float], eta: float) -> dict[str, Any]:
    """Multiple shooting with the symmetry-breaking parameter: node 0 on y = 0 with
    xdot(node 0) = eta (0 on the symmetric family), C free."""
    n = len(nodes)
    xv = np.concatenate([np.concatenate(nodes), np.array(taus)])
    for it in range(1, 121):
        s = [xv[4 * i : 4 * i + 4] for i in range(n)]
        tau = xv[4 * n :]
        if np.any(tau <= 0):
            return {"converged": False, "why": "negative duration", "it": it}
        f = np.zeros(4 * n + 2)
        jm = np.zeros((4 * n + 2, 5 * n))
        for i in range(n):
            res = sh.arc(s[i], float(tau[i]))
            if res is None:
                return {"converged": False, "why": "impact", "it": it}
            sf, phi = res
            j = (i + 1) % n
            f[4 * i : 4 * i + 4] = sf - s[j]
            jm[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] += phi
            jm[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] -= np.eye(4)
            jm[4 * i : 4 * i + 4, 4 * n + i] = base.eom(0.0, sf, MU)
        f[4 * n] = s[0][1]
        jm[4 * n, 1] = 1.0
        f[4 * n + 1] = s[0][2] - eta
        jm[4 * n + 1, 2] = 1.0
        if float(np.max(np.abs(f))) < sh.CONT_TOL:
            s0 = s[0]
            c = base.omega_eff(s0[0], s0[1]) - s0[2] ** 2 - s0[3] ** 2
            return {
                "converged": True,
                "it": it,
                "nodes": [x.tolist() for x in s],
                "taus": tau.tolist(),
                "T": float(np.sum(tau)),
                "C": c,
            }
        dx = np.linalg.lstsq(jm, -f, rcond=None)[0]
        nrm = float(np.linalg.norm(dx[: 4 * n], ord=np.inf))
        if nrm > 0.005:
            dx *= 0.005 / nrm
        xv = xv + dx
    return {"converged": False, "why": "max iterations", "it": 120}


def branch_eta(res: dict[str, Any], label: str) -> list[dict[str, Any]]:
    """Follow the symmetry-breaking parameter eta = xdot at the parent's perpendicular crossing
    from 0 to +-eta_max; a solution with eta != 0 is asymmetric by construction."""
    rot = rotate_to_perp(res)
    if rot is None:
        return [{"label": label, "status": "no perpendicular crossing on the parent"}]
    out = []
    for sign in (+1.0, -1.0):
        nodes, taus = rot
        for eta_abs in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2):
            r = shoot_eta(nodes, taus, sign * eta_abs)
            rec: dict[str, Any] = {"label": label, "eta": sign * eta_abs, **r}
            out.append(rec)
            if not r["converged"]:
                break
            nodes = [np.array(x) for x in r["nodes"]]
            taus = list(r["taus"])
    return out


def branch(
    res: dict[str, Any], c_par: float, label: str, dcs: tuple[float, ...] = (1e-4, -1e-4)
) -> list[dict[str, Any]]:
    """Branch solves at a symmetric parent (converged multiple-shooting solution ``res``)
    whose b is close to +2. Each node is perturbed by the linear image of the start
    perturbation under the arc STMs, so the perturbed node set stays near-continuous."""
    stms = []
    for s, t in zip(res["nodes"], res["taus"], strict=True):
        a = sh.arc(np.array(s), float(t))
        assert a is not None
        stms.append(a[1])
    m = np.eye(4)
    for ph in stms:
        m = ph @ m
    b_par = float(np.trace(m) - 2.0)
    s0 = np.array(res["nodes"][0])
    f = base.eom(0.0, s0, MU)
    w, v = np.linalg.eig(m)
    dirs = []
    for i in np.argsort(np.abs(w - 1.0))[:4]:
        vec = np.real(v[:, i])
        vec = vec - (vec @ f) / (f @ f) * f
        if np.linalg.norm(vec) > 1e-8:
            dirs.append(vec / np.linalg.norm(vec))
    out = []
    for d_i, vec in enumerate(dirs[:3]):
        deltas = [vec]
        for ph in stms[:-1]:
            deltas.append(ph @ deltas[-1])
        for sign in (+1.0, -1.0):
            for eps in (1e-3, 1e-2):
                for dc in dcs:
                    nodes = [
                        np.array(n) + sign * eps * dl
                        for n, dl in zip(res["nodes"], deltas, strict=True)
                    ]
                    r = sh.shoot(nodes, list(res["taus"]), c_par + dc, step_max=0.005, max_it=80)
                    rec: dict[str, Any] = {
                        "label": label,
                        "b_parent": b_par,
                        "C_parent": c_par,
                        "dir": d_i,
                        "eps": eps,
                        "dC": dc,
                        "sign": sign,
                        "converged": r["converged"],
                    }
                    if r["converged"]:
                        sn = np.array(r["nodes"][0])
                        mp = min_perp(sn, r["T"])
                        rec.update(
                            {
                                "T": r["T"],
                                "s0": sn.tolist(),
                                "nodes": r["nodes"],
                                "taus": r["taus"],
                                "min_perp_xdot": mp,
                                "asymmetric": bool(mp > 1e-6),
                            }
                        )
                    out.append(rec)
    return out


def cmd_control() -> None:
    """Continue the 7:3 retrograde asymmetric family (the 7-3b/c class) from its lowest-C grid
    member downward until it becomes symmetric (its pitchfork parent), then branch there."""
    recs = [json.loads(line) for line in GRID.read_text().splitlines()]
    fam = [
        r
        for r in recs
        if r["status"] == "converged"
        and not r["symmetric"]
        and r["p"] == 7
        and r["sense"] == -1
        and round(r["wind_E"]) == -10
    ]
    fam.sort(key=lambda r: r["C"])
    seed = fam[0]
    s, period, c = np.array(seed["s_perigee"]), seed["T"], seed["C"]
    res = ms_solution(s, period, 7, c)
    log = []
    dc = 2e-4
    for k in range(80):
        b = float(np.trace(monodromy(res)) - 2.0)
        sn = np.array(res["nodes"][0])
        mp = min_perp(sn, res["T"])
        log.append({"k": k, "C": c, "T": res["T"], "b": b, "min_perp_xdot": mp})
        print(f"k={k} C={c:.6f} T={res['T']:.5f} b={b:.5f} min|xdot|={mp:.2e}", flush=True)
        if mp < 1e-6:
            break
        c_try = c - dc
        r2 = sh.shoot(
            [np.array(x) for x in res["nodes"]], res["taus"], c_try, step_max=0.005, max_it=80
        )
        while not r2["converged"] and dc > 1e-7:
            dc *= 0.5
            c_try = c - dc
            r2 = sh.shoot(
                [np.array(x) for x in res["nodes"]], res["taus"], c_try, step_max=0.005, max_it=80
            )
        if not r2["converged"]:
            print("lost convergence before the parent", flush=True)
            break
        res, c = r2, c_try
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "control_approach.json").write_text(json.dumps(log, indent=1) + "\n")
    last = log[-1]
    if last["min_perp_xdot"] >= 1e-6:
        print("CONTROL: no symmetric parent reached", flush=True)
        return

    # walk the SYMMETRIC family back up in C to b = +2 (the pitchfork), bisecting
    def solve_at(cc: float, base_res: dict[str, Any]) -> dict[str, Any]:
        return sh.shoot(
            [np.array(x) for x in base_res["nodes"]],
            base_res["taus"],
            cc,
            step_max=0.005,
            max_it=80,
        )

    b_now = float(np.trace(monodromy(res)) - 2.0)
    lo_c, lo_res, lo_b = c, res, b_now
    hi_c = hi_res = None
    for _ in range(60):
        cc = lo_c + 2e-5
        r = solve_at(cc, lo_res)
        if not r["converged"]:
            break
        bb = float(np.trace(monodromy(r)) - 2.0)
        if (lo_b - 2.0) * (bb - 2.0) <= 0:
            hi_c, hi_res = cc, r
            break
        lo_c, lo_res, lo_b = cc, r, bb
    if hi_res is None:
        print("CONTROL: symmetric family b does not reach +2", flush=True)
        return
    for _ in range(30):
        mid = 0.5 * (lo_c + hi_c)
        r = solve_at(mid, lo_res)
        if not r["converged"]:
            break
        bb = float(np.trace(monodromy(r)) - 2.0)
        if (lo_b - 2.0) * (bb - 2.0) > 0:
            lo_c, lo_res, lo_b = mid, r, bb
        else:
            hi_c, hi_res = mid, r
        if hi_c - lo_c < 1e-9:
            break
    c_pf = lo_c
    print(f"pitchfork parent: C={c_pf:.9f} T={lo_res['T']:.6f} b={lo_b:.6f}", flush=True)
    log.append(
        {
            "pitchfork_C": c_pf,
            "pitchfork_T": lo_res["T"],
            "pitchfork_b": lo_b,
            "nodes": lo_res["nodes"],
            "taus": lo_res["taus"],
        }
    )
    (OUT / "control_approach.json").write_text(json.dumps(log, indent=1) + "\n")
    br = branch_eta(lo_res, "7-3b/c parent")
    (OUT / "control_branch.json").write_text(json.dumps(br, indent=1) + "\n")
    # does any asymmetric branch solution lie on the 7-3b/c family? compare with the family's
    # grid members interpolated in C (T to 1e-4)
    for r in br:
        if r.get("converged"):
            print(f"branch eta={r['eta']:+.1e}: C={r['C']:.7f} T={r['T']:.6f}", flush=True)
        else:
            print(f"branch eta={r.get('eta')}: {r.get('why') or r.get('status')}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["control", "parents997", "parentsRR"])
    ap.add_argument("--max-seconds", type=float, default=420.0)
    args = ap.parse_args()
    preflight_search(
        task_no=1000,
        region_id=f"em-complement-pitchfork-{args.cmd}",
        method=MethodCapability(
            genome="symmetric Earth-Moon CR3BP orbits at b = +2, antisymmetric branching",
            corrector="multiple shooting at fixed C (run_1000_shooting.shoot)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar", "asymmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    t0 = time.time()
    if args.cmd == "control":
        cmd_control()
    print(f"done in {time.time() - t0:.0f} s")
    _ = (glob, math, lin)


if __name__ == "__main__":
    main()
