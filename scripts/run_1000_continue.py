"""#1000: continue the asymmetric cycler-class families found by method (1) in C, with the
pre-registered stop rules (note sec. 0.3 and 2.1), and locate b = +2 points along them (the
pitchfork parents of method (2), sec. 2.2).

    uv run python scripts/run_1000_continue.py --family <id> --dir +1|-1

Families are the clusters of asymmetric candidates in data/1000_complement/shooting/grid.jsonl
by (p, q, sense, winding about the Earth, winding about the Moon); the seed is the member with
the median C. One JSONL per family and direction under data/1000_complement/continuation/.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_complement as base
import run_1000_shooting as sh

REPO = Path(__file__).resolve().parent.parent
GRID = REPO / "data" / "1000_complement" / "shooting" / "grid.jsonl"
OUT = REPO / "data" / "1000_complement" / "continuation"
DC0, MAX_STEPS, MAX_HALVINGS = 0.01, 40, 3


def families() -> dict[str, dict[str, Any]]:
    recs = [json.loads(line) for line in GRID.read_text().splitlines()]
    asym = [
        r
        for r in recs
        if r["status"] == "converged" and r["cycler_class_candidate"] and not r["symmetric"]
    ]
    groups: dict[str, list[dict[str, Any]]] = {}
    for r in asym:
        sense = "pro" if r["sense"] > 0 else "ret"
        key = f"{r['p']}_{r['q']}_{sense}_wE{round(r['wind_E'])}_wM{round(r['wind_M'])}"
        groups.setdefault(key, []).append(r)
    out = {}
    for key, g in groups.items():
        g.sort(key=lambda r: r["C"])
        out[key] = {"n": len(g), "C_range": [g[0]["C"], g[-1]["C"]], "seed": g[len(g) // 2]}
    return out


def nodes_from(s0: np.ndarray, period: float, n: int) -> tuple[list[np.ndarray], list[float]]:
    ts = np.linspace(0.0, period, n + 1)
    sol = solve_ivp(
        base.eom,
        (0.0, period),
        s0,
        args=(base.MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        t_eval=ts[:-1],
    )
    return [sol.y[:, i].copy() for i in range(n)], [period / n] * n


def predict(
    nodes: list[np.ndarray],
    taus: list[float],
    older: tuple[list[np.ndarray], list[float]] | None,
    dc_prev: float,
    step: float,
) -> tuple[list[np.ndarray], list[float]]:
    """Secant predictor from the previous member (a corrector aid, not a rule change)."""
    if older is None or dc_prev == 0.0 or os.environ.get("NO_PREDICTOR"):
        return nodes, taus
    w = step / dc_prev
    pn = [a + w * (a - b) for a, b in zip(nodes, older[0], strict=True)]
    pt = [a + w * (a - b) for a, b in zip(taus, older[1], strict=True)]
    return pn, pt


def ms_monodromy(nodes: list[list[float]], taus: list[float]) -> tuple[float, float]:
    """b = tr(M) - 2 and max |lambda| from the product of the per-arc STMs of a converged
    multiple-shooting solution (single shooting over a whole period of a lambda ~ 1e2-1e3
    orbit drifts off the orbit and gives a wrong monodromy)."""
    m = np.eye(4)
    for s, t in zip(nodes, taus, strict=True):
        res = sh.arc(np.array(s), float(t))
        assert res is not None
        m = res[1] @ m
    return float(np.trace(m) - 2.0), float(np.max(np.abs(np.linalg.eigvals(m))))


def attempt(
    nodes: list[np.ndarray],
    taus: list[float],
    older: tuple[list[np.ndarray], list[float]] | None,
    dc_prev: float,
    step: float,
    c_t: float,
) -> dict[str, Any]:
    """One continuation solve: secant predictor first, then the previous nodes unchanged."""
    r = sh.shoot(*predict(nodes, taus, older, dc_prev, step), c_t)
    return r if r["converged"] else sh.shoot(nodes, taus, c_t)


def member(s0: np.ndarray, period: float, p: int, ref: np.ndarray | None = None) -> dict[str, Any]:
    # move to a perigee for the section / classification, as in method (1); along a
    # continuation, take the perigee nearest the previous member's (so nodes correspond)
    sp = solve_ivp(
        base.eom,
        (0.0, period),
        s0,
        args=(base.MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=base._dr1,
    )
    ys = [np.asarray(y) for y in sp.y_events[0]]
    sper = ys[0] if ref is None else min(ys, key=lambda y: float(np.linalg.norm(y[:2] - ref[:2])))
    m = base.classify(sper, period)
    m.pop("perp_crossings", None)
    m["s_perigee"] = sper.tolist()
    m["C"] = base.omega_eff(sper[0], sper[1]) - sper[2] ** 2 - sper[3] ** 2
    m["T"] = period
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default=None)
    ap.add_argument("--dir", type=int, choices=[1, -1], default=1)
    ap.add_argument("--max-seconds", type=float, default=420.0)
    ap.add_argument("--list", action="store_true")
    # #1050: extend an existing run past its cap: new file, start from its last converged member
    ap.add_argument("--extend", action="store_true")
    ap.add_argument("--max-steps", type=int, default=MAX_STEPS)
    ap.add_argument("--c-max", type=float, default=None)
    args = ap.parse_args()
    fams = families()
    if args.list or args.family is None:
        for k, v in fams.items():
            print(k, v["n"], [round(c, 4) for c in v["C_range"]])
        return
    preflight_search(
        task_no=1000,
        region_id=f"em-complement-continuation-{args.family}-{args.dir:+d}",
        method=MethodCapability(
            genome="asymmetric Earth-Moon CR3BP family, multiple shooting at fixed C",
            corrector="minimum-norm Newton (run_1000_shooting.shoot) with C constraint",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar", "asymmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=MAX_STEPS,
    )
    fam = fams[args.family]
    seed = fam["seed"]
    p = seed["p"]
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{args.family}_{'p' if args.dir > 0 else 'm'}.jsonl"
    if args.extend:
        src = [json.loads(line) for line in f.read_text().splitlines()]
        f = OUT / f"{args.family}_{'p' if args.dir > 0 else 'm'}_ext.jsonl"
        if not f.exists():
            last_ok = [x for x in src if "stop" not in x][-1]
            start = {**last_ok, "k": 0, "extended_from_C": last_ok["C"]}
            f.write_text(json.dumps(start) + "\n")
    rows = [json.loads(line) for line in f.read_text().splitlines()] if f.exists() else []
    if rows and rows[-1].get("stop"):
        print("already stopped:", rows[-1]["stop"])
        return
    if rows:
        last = rows[-1]
        s0, period, dc = np.array(last["s_perigee"]), last["T"], last["dC"]
        topo0 = tuple(rows[0]["topology"])
    else:
        s0, period, dc = np.array(seed["s_perigee"]), seed["T"], DC0
        m = member(s0, period, p)
        m.update({"k": 0, "dC": dc, "topology": [round(m["wind_E"]), round(m["wind_M"])]})
        n0, t0s = nodes_from(s0, period, p)
        r0 = sh.shoot(n0, t0s, m["C"])
        if r0["converged"]:
            m["b_h_single_shooting"] = m["b_h"]
            m["b_h"], m["lambda_max"] = ms_monodromy(r0["nodes"], r0["taus"])
        topo0 = tuple(m["topology"])
        rows.append(m)
        with f.open("a") as fh:
            fh.write(json.dumps(m) + "\n")
    t0 = time.time()
    while time.time() - t0 < args.max_seconds:
        prev = rows[-1]
        k = prev["k"] + 1
        if k > args.max_steps:
            stop = "max steps"
            prev["stop"] = stop
            with f.open("a") as fh:
                fh.write(json.dumps({**prev, "k": k, "stop": stop}) + "\n")
            break
        nodes, taus = nodes_from(np.array(prev["s_perigee"]), prev["T"], p)
        # secant predictor from the previous member (a corrector aid, not a rule change)
        older = None
        if len(rows) >= 2 and "s_perigee" in rows[-2]:
            older = nodes_from(np.array(rows[-2]["s_perigee"]), rows[-2]["T"], p)
        dc_prev = abs(prev["C"] - rows[-2]["C"]) if older is not None else 0.0

        c_new = prev["C"] + args.dir * dc

        res = attempt(nodes, taus, older, dc_prev, dc, c_new)
        halv = 0
        while not res["converged"] and halv < MAX_HALVINGS:
            dc *= 0.5
            halv += 1
            c_new = prev["C"] + args.dir * dc
            res = attempt(nodes, taus, older, dc_prev, dc, c_new)
        stop = None
        if not res["converged"]:
            stop = f"loss of convergence ({res.get('why')})"
            with f.open("a") as fh:
                fh.write(json.dumps({**prev, "k": k, "stop": stop}) + "\n")
            break
        m = member(np.array(res["s0"]), res["T"], p, ref=np.array(prev["s_perigee"]))
        m.update({"k": k, "dC": dc, "topology": [round(m["wind_E"]), round(m["wind_M"])]})
        m["b_h_single_shooting"] = m["b_h"]
        m["b_h"], m["lambda_max"] = ms_monodromy(res["nodes"], res["taus"])
        if tuple(m["topology"]) != topo0:
            stop = f"topology change {topo0} -> {tuple(m['topology'])}"
        if m["symmetric"]:
            stop = "became symmetric (met a symmetric orbit: candidate pitchfork parent)"
        below = m["perigee_km"] < base.FLOOR_E or m["periselene_km"] < base.FLOOR_M
        if below:
            stop = "impact at floor"
        if (prev["b_h"] + 2.0) * (m["b_h"] + 2.0) < 0:
            stop = "period doubling (b crosses -2)"
        if args.c_max is not None and args.dir * (m["C"] - args.c_max) >= 0:
            stop = f"reached C limit {args.c_max}"
        m["b_crossed_plus2"] = bool((prev["b_h"] - 2.0) * (m["b_h"] - 2.0) < 0)
        if stop:
            m["stop"] = stop
        rows.append(m)
        with f.open("a") as fh:
            fh.write(json.dumps(m) + "\n")
        print(
            f"{time.strftime('%H:%M:%S')} k={k} C={m['C']:.6f} T={m['T']:.5f} "
            f"pe={m['perigee_alt_km']:.0f} ps={m['periselene_alt_km']:.0f} b={m['b_h']:.4g} "
            f"sym={m['symmetric']} {stop or ''}",
            flush=True,
        )
        if stop:
            break
        if halv == 0 and dc < DC0:
            dc = min(DC0, dc * 1.5)
    print(f"{f.name}: {len(rows)} members")


if __name__ == "__main__":
    main()
