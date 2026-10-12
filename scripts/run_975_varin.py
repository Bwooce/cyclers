"""#975: the Varin closed-family cascade i_5-i_8 at the satellite mass ratios.

Pre-registration: docs/notes/2026-10-11-975-varin-cascade-satellite-mass-ratios.md.

    uv run python scripts/run_975_varin.py p1     # positive control P1: Table 1's 8 orbits
    uv run python scripts/run_975_varin.py curve --name N --start mu2p --y2-dir 1 --mu-stop X
    uv run python scripts/run_975_varin.py loop --name L --curve N
    uv run python scripts/run_975_varin.py controls-report

Every continuation is checkpointed (data/975_varin/curves/<name>.jsonl, loops/<name>.jsonl), resumes
from its last points, and stops at --max-seconds (default 420) with "time budget".

P1 (hard stop): each Table 1 orbit (Bruno & Varin 2009, printed to 8 digits) is converted to the
project's frame and (a) integrated at its printed mu with no correction (|y1| at the x2 = 0 crossing
nearest T/2 <= 1e-8), (b) its Jacobi constant compared with C_printed - mu (1 - mu) (1e-8), (c)
its modified traces compared with the printed values (Tr~ = 2 to 1e-4, Tr~_v to 1e-6), and (d) the
critical-point solver, started from the printed orbit with its four unknowns perturbed by 1e-4
(relative), returns to the printed (x1(0), y2(0), T, mu) to 5e-8 relative.
Output: data/975_varin/p1_control.json.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core import cr3bp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search import varin_closed_families as vcf

REPO = Path(__file__).resolve().parent.parent
TABLES = REPO / "data" / "sources" / "varin-bruno-2009-closed-families-tables.yaml"
OUT = REPO / "data" / "975_varin"


def table1() -> list[dict[str, Any]]:
    """Table 1 rows with their mu (Table 2) and the converted crossing."""
    d = yaml.safe_load(TABLES.read_text())
    mus = {r["k"]: (float(r["mu_p"]), float(r["mu_pp"])) for r in d["table2_bifurcation_mu"]}
    rows = []
    for r in d["table1_bifurcation_orbits"]:
        lab = r["mu"]
        k = lab[2]
        mu = mus[k][1] if lab.endswith("pp") else mus[k][0]
        a, e = float(r["a"]), float(r["e"])
        x1, y2 = vcf.crossing_from_varin(a, e, mu)
        rows.append(
            {
                "label": lab,
                "k": int(k),
                "kind": "elliptic" if lab.endswith("pp") else "saddle",
                "mu": mu,
                "a": a,
                "e": e,
                "x1": x1,
                "y2": y2,
                "T": 2.0 * math.pi * float(r["T_over_2pi"]),
                "C_printed": float(r["C"]),
                "Trv_printed": float(r["Trv"]),
            }
        )
    return rows


def p1() -> None:
    out = []
    ok_all = True
    for r in table1():
        t0 = time.time()
        mu, x1, y2, tt = r["mu"], r["x1"], r["y2"], r["T"]
        tc, y1 = vcf.crossing_near_half(x1, y2, tt, mu)
        cj = cr3bp.jacobi_constant(vcf.state_from_crossing(x1, y2, mu), mu)
        dc = cj - vcf.jacobi_from_printed(r["C_printed"], mu)
        tr, trv = vcf.traces(x1, y2, tt, mu)
        z = np.array([x1, y2, tt, mu])
        cp = vcf.solve_critical_point(z * (1.0 + 1e-4 * np.array([1.0, -1.0, 1.0, -1.0])))
        zc = np.array([cp.x1, cp.y2, cp.t_full, cp.mu])
        rel = float(np.max(np.abs((zc - z) / z)))
        checks = {
            "a_y1": y1 <= 1e-8,
            "b_jacobi": abs(dc) <= 1e-8,
            "c_tr": abs(vcf.modified_trace(tr) - 2.0) <= 1e-4,
            "c_trv": abs(vcf.modified_trace(trv) - r["Trv_printed"]) <= 1e-6,
            "d_solver": cp.converged and rel <= 5e-8,
        }
        ok = all(checks.values())
        ok_all &= ok
        rec = {
            **r,
            "y1_at_half": y1,
            "T_over_2pi_crossing": tc / math.pi,
            "jacobi": cj,
            "dC": dc,
            "Tr_mod": vcf.modified_trace(tr),
            "Trv_mod": vcf.modified_trace(trv),
            "solver": {
                "converged": cp.converged,
                "iterations": cp.iterations,
                "residual": cp.residual,
                "max_rel_diff": rel,
            },
            "checks": checks,
            "pass": ok,
            "seconds": round(time.time() - t0, 1),
        }
        out.append(rec)
        print(
            f"{r['label']}: pass {ok} y1 {y1:.1e} dC {dc:.1e} Tr~ {rec['Tr_mod']:.6f} "
            f"Trv~ {rec['Trv_mod']:.8f} ({r['Trv_printed']}) solver rel {rel:.1e} "
            f"({rec['seconds']} s)",
            flush=True,
        )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "p1_control.json").write_text(json.dumps({"pass": ok_all, "rows": out}, indent=1) + "\n")
    print("P1 PASS" if ok_all else "P1 FAIL: stop", flush=True)


CURVES = OUT / "curves"
LOOPS = OUT / "loops"


def _jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(x) for x in path.read_text().splitlines()] if path.exists() else []


def _log(fh: Any, rec: dict[str, Any]) -> None:
    fh.write(json.dumps(rec) + "\n")
    fh.flush()


def cmd_curve(name: str, start: str, y2_dir: float, mu_stop: float, max_s: float) -> None:
    """y2-extremum curve from a Table 1 orbit, toward mu_stop, checkpointed."""
    CURVES.mkdir(parents=True, exist_ok=True)
    f = CURVES / f"{name}.jsonl"
    old = [r for r in _jsonl(f) if "w" in r]
    ends = [r for r in _jsonl(f) if "end" in r]
    if ends:
        print(f"{name}: already ended: {ends[-1]}", flush=True)
        return
    row = {r["label"]: r for r in table1()}[start]
    if len(old) >= 2:
        w0 = np.array(old[-1]["w"])
        direction = w0 - np.array(old[-2]["w"])
    else:
        w0 = np.array([row["x1"], row["y2"], row["T"], math.log(row["mu"])])
        direction = np.array([0.0, y2_dir, 0.0, math.copysign(1.0, math.log(mu_stop) - w0[3])])
    t0 = time.time()
    with f.open("a") as fh:
        if not old:
            _log(
                fh,
                {"w": w0.tolist(), "mu": row["mu"], "start": start, "t": time.strftime("%H:%M:%S")},
            )

        def cb(p: vcf.CurvePoint) -> None:
            _log(fh, {"w": p.w, "mu": p.mu, "t": time.strftime("%H:%M:%S")})
            print(
                f"{time.strftime('%H:%M:%S')} {name}: mu {p.mu:.6e} x1 {p.w[0]:.6f} "
                f"y2 {p.w[1]:.6f} T~ {p.w[2] / (2 * math.pi):.5f}",
                flush=True,
            )

        stats: dict[str, int] = {}
        pts, why = vcf.follow_extremum_curve(
            w0,
            direction=direction,
            mu_stop=mu_stop,
            on_point=cb,
            deadline=t0 + max_s,
            stats=stats,
        )
        if why != "time budget":
            _log(fh, {"end": why, "rejected": stats, "t": time.strftime("%H:%M:%S")})
    print(f"{name}: {why} ({len(pts)} points this call; A2 rejections {stats})", flush=True)


def cmd_loop(
    name: str, curve: str, max_s: float, max_members: int = 20000, x1_dir: float = 1.0
) -> None:
    """Closed-family trace at the mu of a curve's final point (which must have reached mu_stop)."""
    LOOPS.mkdir(parents=True, exist_ok=True)
    c = _jsonl(CURVES / f"{curve}.jsonl")
    if not c or c[-1].get("end") != "reached mu_stop":
        print(f"{name}: curve {curve} has not reached its mu_stop", flush=True)
        return
    start = [r for r in c if "w" in r][-1]
    mu = float(start["mu"])
    f = LOOPS / f"{name}.jsonl"
    old = [r for r in _jsonl(f) if "p" in r]
    if [r for r in _jsonl(f) if "end" in r]:
        print(f"{name}: already ended", flush=True)
        return
    w = start["w"]
    origin = np.array([w[0], w[1], w[2] / (2 * math.pi)])
    start_tangent = None
    if len(old) >= 2:
        p_last = np.array(old[-1]["p"])
        direction = p_last - np.array(old[-2]["p"])
        start_tangent = np.array(old[1]["p"]) - np.array(old[0]["p"])
    else:
        p_last = origin
        direction = np.array([math.copysign(1.0, x1_dir), 0.0, 0.0])
    t0 = time.time()
    with f.open("a") as fh:
        if not old:
            _log(fh, {"p": origin.tolist(), "mu": mu, "curve": curve})

        count = [len(old)]
        stats: dict[str, int] = {}

        def cb(m: vcf.LoopMember) -> None:
            _log(fh, {"p": [m.x1, m.y2, m.t_full / (2 * math.pi)], "t": time.strftime("%H:%M:%S")})
            count[0] += 1
            if count[0] % 50 == 0:
                d = float(np.linalg.norm(np.array([m.x1, m.y2, m.t_full / (2 * math.pi)]) - origin))
                print(
                    f"{time.strftime('%H:%M:%S')} {name}: {count[0]} members, x1 {m.x1:.5f} "
                    f"y2 {m.y2:.5f} T~ {m.t_full / (2 * math.pi):.5f}, distance to start {d:.3e}, "
                    f"A2 rejections {stats}",
                    flush=True,
                )

        members, why, closure = vcf.trace_loop(
            float(p_last[0]),
            float(p_last[1]),
            2 * math.pi * float(p_last[2]),
            mu,
            direction=direction,
            origin=origin,
            n_done=max(len(old) - 1, 0),
            on_member=cb,
            deadline=t0 + max_s,
            max_members=max_members,
            start_tangent=start_tangent,
            stats=stats,
        )
        rec = {"end": why, "closure": closure, "rejected": stats, "t": time.strftime("%H:%M:%S")}
        if why != "time budget":
            _log(fh, rec)
    print(
        f"{name}: {why} closure {closure:.2e} ({len(members)} members this call; "
        f"A2 rejections this call {stats})",
        flush=True,
    )


def cmd_compare(loop_a: str, loop_b: str) -> None:
    """P2: does each loop's start lie on the other loop (distance < 1e-6 after refinement)?"""
    out: dict[str, Any] = {}
    for a, b in ((loop_a, loop_b), (loop_b, loop_a)):
        ra = _jsonl(LOOPS / f"{a}.jsonl")
        rb = _jsonl(LOOPS / f"{b}.jsonl")
        mu = float(ra[0]["mu"])
        assert abs(mu - float(rb[0]["mu"])) < 1e-15
        members = [np.array(r["p"]) for r in ra if "p" in r]
        d = vcf.distance_to_loop(np.array(rb[0]["p"]), members, mu)
        out[f"start of {b} on {a}"] = d
        out[f"{a} end"] = [r for r in ra if "end" in r]
        print(f"start of {b} on {a}: {d:.2e}", flush=True)
    out["coincide"] = all(v < 1e-6 for k, v in out.items() if k.startswith("start"))
    (OUT / f"compare_{loop_a}_{loop_b}.json").write_text(json.dumps(out, indent=1) + "\n")
    print("coincide" if out["coincide"] else "DO NOT coincide", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("p1")
    sp = sub.add_parser("curve")
    sp.add_argument("--name", required=True)
    sp.add_argument("--start", required=True)
    sp.add_argument("--y2-dir", type=float, default=1.0)
    sp.add_argument("--mu-stop", type=float, required=True)
    sp.add_argument("--max-seconds", type=float, default=420.0)
    sp = sub.add_parser("compare")
    sp.add_argument("--a", required=True)
    sp.add_argument("--b", required=True)
    sp = sub.add_parser("loop")
    sp.add_argument("--name", required=True)
    sp.add_argument("--curve", required=True)
    sp.add_argument("--max-seconds", type=float, default=420.0)
    sp.add_argument("--max-members", type=int, default=20000)
    sp.add_argument("--x1-dir", type=float, default=1.0, help="start direction in x1(0) (A2 v2)")
    args = ap.parse_args()
    preflight_search(
        task_no=975,
        region_id=f"varin-closed-families-{args.cmd}",
        method=MethodCapability(
            genome="Varin closed families i_k (symmetric planar CR3BP)",
            corrector="Varin 4x4 critical-point Newton; symmetric single shooting",
            capability_tags=frozenset({"cr3bp", "planar", "symmetric"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    if args.cmd == "p1":
        p1()
    elif args.cmd == "curve":
        cmd_curve(args.name, args.start, args.y2_dir, args.mu_stop, args.max_seconds)
    elif args.cmd == "loop":
        cmd_loop(args.name, args.curve, args.max_seconds, args.max_members, args.x1_dir)
    elif args.cmd == "compare":
        cmd_compare(args.a, args.b)


if __name__ == "__main__":
    main()
