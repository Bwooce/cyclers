"""#1034: gc-1 in continuous gravity, joint-sigma continuation (pre-registration:
docs/notes/2026-10-08-1034-gc1-continuous-sigma.md). Reuses the #1004 corrector
(``scripts/run_1004_gc2.py``) with gc-1's structure.

    uv run python scripts/run_1034_gc1.py --stage sigma     # natural continuation up to 1
    uv run python scripts/run_1034_gc1.py --stage monitor   # pass a fold (monitor variable)
    uv run python scripts/run_1034_gc1.py --stage two --at 0.1,0.12
    uv run python scripts/run_1034_gc1.py --stage sigma1
    uv run python scripts/run_1034_gc1.py --stage ias15 --at 1.0

Checkpoints in data/1034_gc1/; live log in data/1034_gc1/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


G = _load("run_1004_gc2", ROOT / "scripts" / "run_1004_gc2.py")
V = _load("run_1004_verify", ROOT / "scripts" / "run_1004_verify.py")
G.KEY = "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s"
G.OUT = ROOT / "data" / "1034_gc1"
OUT = G.OUT
PC_VINF = (2.397, 1.807)
RADII = {G.GAN: SATELLITES[G.GAN].radius_eq_km, G.CAL: SATELLITES[G.CAL].radius_eq_km}


def _rec(name: str) -> dict[str, Any]:
    path = OUT / name
    return json.loads(path.read_text()) if path.exists() else {"points": []}


def _put(name: str, rec: dict[str, Any]) -> None:
    (OUT / name).write_text(json.dumps(rec))


def impact(p: Any, z: Any, sigma: float) -> bool:
    d = G.describe(p, z)
    return any(n["rp_km"] <= sigma * RADII[n["moon"]] for n in d["nodes"])


def point(p: Any, z: Any, sigma: float) -> dict[str, Any]:
    d = G.describe(p, z)
    return {
        "sigma": sigma,
        "z": np.asarray(z).tolist(),
        "vinf": [n["vinf_kms"] for n in d["nodes"]],
        "rp_km": [n["rp_km"] for n in d["nodes"]],
        "rp_over_sigma": [n["rp_km"] / sigma for n in d["nodes"]],
        "turn_deg": [n["turn_deg"] for n in d["nodes"]],
        "c1c2_a_km": d["c1c2_mid_a_km"],
        "c1c2_e": d["c1c2_mid_e"],
    }


def stage_sigma() -> None:
    """Natural continuation in sigma from 0.05, step x1.2, halving to a 1e-5 absolute step."""
    p = G.build()
    rec = _rec("sigma.json")
    pts = rec["points"]
    if not pts:
        sg = 0.05
        G.set_sc(p, sg, "J")
        z, info = G.solve(p, G.scale_offsets(p, p.seed, sg), 250, f"sigma={sg:.5g}")
        if not G.converged(info):
            G._log("sigma: no convergence at 0.05; stop (pre-registered)")
            return
        pt = point(p, z, sg)
        ok_id = abs(pt["vinf"][0] - PC_VINF[0]) < 0.02 and abs(pt["vinf"][1] - PC_VINF[1]) < 0.02
        pt["identity_pass"] = ok_id
        pts.append(pt)
        _put("sigma.json", rec)
        G._log(f"sigma=0.05 identity_pass={ok_id} vinf {pt['vinf']}")
        if not ok_id:
            return
    if not pts[0].get("identity_pass") and not _rec("identity.json").get("pass"):
        G._log("sigma: identity not established (sec. 2 / amendment 1); stop")
        return
    fac = rec.get("fac", 1.2)
    while pts[-1]["sigma"] < 1.0:
        s1, z1 = pts[-1]["sigma"], np.asarray(pts[-1]["z"])
        sg = min(1.0, s1 * fac)
        if sg - s1 < 1e-5:
            G._log(f"sigma: step below 1e-5 at {s1:.6f}: stop (fold or numerical; run checks)")
            rec["stopped_at"] = s1
            _put("sigma.json", rec)
            return
        if len(pts) >= 2:
            s0, z0 = pts[-2]["sigma"], np.asarray(pts[-2]["z"])
            zs = z1 + (z1 - z0) * (math.log(sg) - math.log(s1)) / (math.log(s1) - math.log(s0))
        else:
            zs = G.scale_offsets(p, z1, sg / s1)
        G.set_sc(p, sg, "J")
        z, info = G.solve(p, zs, 40, f"sigma={sg:.6g}")
        if G.converged(info) and not impact(p, z, sg):
            pts.append(point(p, z, sg))
            fac = min(1.2, 1.0 + 2.0 * (fac - 1.0))
            G._log(
                f"sigma={sg:.6g} ok vinf {np.round(pts[-1]['vinf'], 4).tolist()} "
                f"rp/s {np.round(pts[-1]['rp_over_sigma'], 1).tolist()}"
            )
        elif G.converged(info):
            G._log(f"sigma={sg:.6g}: converged INSIDE a scaled body radius (impact end); stop")
            rec["impact_at"] = sg
            _put("sigma.json", rec)
            return
        else:
            fac = 1.0 + 0.5 * (fac - 1.0)
            G._log(f"sigma={sg:.6g} failed; factor -> {fac:.6g}")
        rec["fac"] = fac
        _put("sigma.json", rec)
    G._log("sigma: reached 1")


def stage_identity() -> None:
    """Amendment 1: identity by the sigma -> 0 limit (0.02, 0.01, 0.005)."""
    p = G.build()
    start = _rec("sigma.json")["points"][0]
    z, s_prev = np.asarray(start["z"]), float(start["sigma"])
    rows = []
    for sg in (0.02, 0.01, 0.005):
        G.set_sc(p, sg, "J")
        zn, info = G.solve(p, G.scale_offsets(p, z, sg / s_prev), 80, f"identity sigma={sg}")
        ok = bool(G.converged(info))
        row = point(p, zn, sg) if ok else {"sigma": sg}
        row["converged"] = ok
        rows.append(row)
        G._log(f"identity sigma={sg}: converged={ok} vinf {row.get('vinf')}")
        if not ok:
            break
        z, s_prev = zn, sg
    last = rows[-1]
    passed = bool(
        last.get("converged")
        and abs(last["vinf"][0] - 2.3972) < 0.005
        and abs(last["vinf"][1] - 1.8067) < 0.005
        and all(
            abs(a["vinf"][0] - 2.3972) > abs(b["vinf"][0] - 2.3972)
            for a, b in zip([start, *rows[:-1]], rows, strict=True)
        )
    )
    _put("identity.json", {"points": rows, "pass": passed})
    G._log(f"identity: pass={passed}")


def gap_of(p: Any, z: Any, sigma: float) -> float:
    _, jz, _ = G.residual_and_jac(p, z)
    js = G.j_s(p, z, sigma, "J")
    d = np.linalg.norm(jz, axis=0)
    d[d == 0.0] = 1.0
    sv = np.linalg.svd(
        np.hstack([jz / d, (js / max(np.linalg.norm(js), 1e-300))[:, None]]), compute_uv=False
    )
    return float(sv[-2] / max(sv[-1], 1e-300))


def stage_monitor(n_points: int) -> None:
    """Pass a fold: sigma unknown, the fastest-changing state component (scaled) stepped and fixed.
    (The SVD tangent is used only when the gap test passes; the gap is logged, never bypassed.)"""
    p = G.build()
    rec = _rec("monitor.json")
    if not rec["points"]:
        sig = _rec("sigma.json")["points"]
        rec["points"] = [sig[-2], sig[-1]]
    for _ in range(n_points):
        a, b = rec["points"][-2], rec["points"][-1]
        za, zb = np.asarray(a["z"]), np.asarray(b["z"])
        G.set_sc(p, b["sigma"], "J")
        _, jz, _ = G.residual_and_jac(p, zb)
        d = np.linalg.norm(jz, axis=0)
        d[d == 0.0] = 1.0
        gap = gap_of(p, zb, b["sigma"])
        j = int(np.argmax(np.abs((zb - za) * d)))
        h = float(zb[j] - za[j])
        target = zb[j] + h
        y = np.concatenate([zb + (zb - za), [b["sigma"] + (b["sigma"] - a["sigma"])]])
        ok = False
        for _half in range(6):
            for _it in range(20):
                G.set_sc(p, float(y[-1]), "J")
                r, jzz, info = G.residual_and_jac(p, y[:-1])
                if G.converged(info) and abs(y[j] - target) < 1e-9 * max(1.0, abs(target)):
                    ok = True
                    break
                js = G.j_s(p, y[:-1], float(y[-1]), "J")
                row = np.zeros(G.NZ + 1)
                row[j] = 1.0
                jac = np.vstack([np.hstack([jzz, js[:, None]]), row[None, :]])
                f = np.concatenate([r, [y[j] - target]])
                dcol = np.linalg.norm(jac, axis=0)
                dcol[dcol == 0.0] = 1.0
                step = np.linalg.lstsq(jac / dcol, -f, rcond=None)[0] / dcol
                f0 = float(np.linalg.norm(f))
                for alpha in (1.0, 0.5, 0.25, 0.125):
                    yt = y + alpha * step
                    G.set_sc(p, float(yt[-1]), "J")
                    rt = G.residual_and_jac(p, yt[:-1], False)[0]
                    if float(np.linalg.norm(np.concatenate([rt, [yt[j] - target]]))) < f0:
                        y = yt
                        break
                else:
                    break
            if ok:
                break
            h *= 0.5
            target = zb[j] + h
            y = np.concatenate(
                [zb + 0.5 * (y[:-1] - zb), [b["sigma"] + 0.5 * (float(y[-1]) - b["sigma"])]]
            )
        if not ok:
            G._log("monitor: six halvings without convergence: NUMERICAL stop")
            break
        sg = float(y[-1])
        G.set_sc(p, sg, "J")
        if impact(p, y[:-1], sg):
            G._log(f"monitor: impact at sigma {sg:.6g}; stop")
            break
        pt = point(p, y[:-1], sg)
        pt["monitor_index"] = j
        pt["gap_svd"] = gap
        rec["points"].append(pt)
        _put("monitor.json", rec)
        G._log(
            f"monitor sigma={sg:.6g} idx {j} gap {gap:.2f} "
            f"vinf {np.round(pt['vinf'], 4).tolist()} "
            f"rp/s {np.round(pt['rp_over_sigma'], 1).tolist()}"
        )
        if sg >= 1.0:
            G._log("monitor: reached sigma 1")
            break


def stage_two(at: list[float]) -> None:
    p = G.build()
    lower = _rec("sigma.json")["points"]
    upper = [q for q in _rec("monitor.json")["points"][2:]]
    out = {}
    for sg in at:
        zl = np.asarray(min(lower, key=lambda q: abs(q["sigma"] - sg))["z"])
        zu = np.asarray(min(upper, key=lambda q: abs(q["sigma"] - sg))["z"])
        res = {}
        for lab, z0 in (("lower", zl), ("upper", zu)):
            G.set_sc(p, sg, "J")
            z, info = G.solve(p, z0, 40, f"two {lab} {sg}")
            res[lab] = point(p, z, sg) if G.converged(info) else None
        row: dict[str, Any] = {"lower": res["lower"], "upper": res["upper"]}
        if res["lower"] and res["upper"]:
            row["d_vinf_C"] = abs(res["lower"]["vinf"][1] - res["upper"]["vinf"][1])
            row["d_rp_over_sigma_C"] = abs(
                res["lower"]["rp_over_sigma"][1] - res["upper"]["rp_over_sigma"][1]
            )
            row["fold_supported"] = row["d_vinf_C"] > 0.02 and row["d_rp_over_sigma_C"] > 100.0
        out[str(sg)] = row
        G._log(
            f"two sigma={sg}: {row.get('d_vinf_C')} {row.get('d_rp_over_sigma_C')} "
            f"{row.get('fold_supported')}"
        )
    rec = _rec("checks.json")
    rec["two"] = out
    _put("checks.json", rec)


def stage_sigma1() -> None:
    p = G.build()
    starts = {"patched_conic_seed": p.seed}
    lo = _rec("sigma.json")["points"]
    if lo:
        starts["lower_last"] = G.scale_offsets(p, np.asarray(lo[-1]["z"]), 1.0 / lo[-1]["sigma"])
    up = _rec("monitor.json")["points"]
    if len(up) > 2:
        starts["upper_last"] = G.scale_offsets(p, np.asarray(up[-1]["z"]), 1.0 / up[-1]["sigma"])
    out = {}
    for lab, z0 in starts.items():
        G.set_sc(p, 1.0, "J")
        z, info = G.solve(p, z0, 40, f"sigma1 {lab}")
        out[lab] = {
            "converged": bool(G.converged(info)),
            "dr": max(info["dr"]),
            "dv": max(info["dv"]),
        }
        if G.converged(info):
            out[lab].update(point(p, z, 1.0))
        G._log(f"sigma1 {lab}: {out[lab]['converged']} dv {out[lab]['dv']:.3e}")
    rec = _rec("checks.json")
    rec["sigma1"] = out
    _put("checks.json", rec)


def stage_ias15(at: float) -> None:
    p = G.build()
    pts = _rec("sigma.json")["points"] + _rec("monitor.json")["points"]
    q = min(pts, key=lambda t: abs(t["sigma"] - at))
    z, sg = np.asarray(q["z"]), float(q["sigma"])
    G.set_sc(p, sg, "J")
    nd = G.nodes(p, z)

    class Exact:
        def position(self, moon: str, t: float) -> Any:
            return p.circ.state(moon, t)[0]

    rows = []
    for k in range(4):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for lab, x0, t0 in (("fwd", xp, tp), ("bwd", xn, tn)):
            xd, _ = G.arc(p, x0, t0, tm)
            xi = V._ias15_scaled(p, x0, t0, tm, Exact())
            rows.append(
                {
                    "leg": k,
                    "dir": lab,
                    "dr_km": float(np.linalg.norm(xi[:3] - xd[:3])),
                    "dv_kms": float(np.linalg.norm(xi[3:] - xd[3:])),
                }
            )
    rec = _rec("checks.json")
    rec.setdefault("ias15", {})[f"{sg:.6g}"] = rows
    _put("checks.json", rec)
    G._log(
        f"ias15 sigma={sg:.6g}: max dr {max(r['dr_km'] for r in rows):.2e} "
        f"dv {max(r['dv_kms'] for r in rows):.2e}"
    )


def main() -> None:
    preflight_search(
        task_no=1034,
        region_id="gc1-joint-sigma-continuation-continuous",
        method=MethodCapability(
            genome="gc-1 (one structure), R-S circular coplanar, Ganymede + Callisto massive",
            corrector="forward-backward shooting, DOP853 + STM; natural and monitor continuation",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage", choices=["sigma", "identity", "monitor", "two", "sigma1", "ias15"], required=True
    )
    ap.add_argument("--at", default="1.0")
    ap.add_argument("--points", type=int, default=10)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "sigma":
        stage_sigma()
    elif args.stage == "identity":
        stage_identity()
    elif args.stage == "monitor":
        stage_monitor(args.points)
    elif args.stage == "two":
        stage_two([float(v) for v in args.at.split(",")])
    elif args.stage == "sigma1":
        stage_sigma1()
    else:
        stage_ias15(float(args.at))


if __name__ == "__main__":
    main()
