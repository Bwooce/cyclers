"""#1073: Pluto continuous-gravity lane driver.

Pre-registration: docs/notes/2026-10-12-1073-pluto-continuous-gravity-lane-prereg.md.
Code: src/cyclerfinder/nbody/pluto_cg.py. Output: data/1073_pluto_cg/ (one JSON per run, resumable;
live logs in data/1073_pluto_cg/live/, gitignored).

Subcommands:
  pc1 --row K           PC1: Bruno-Varin 2009 Table 4 (mu = 0.1) orbit K at Pluto scale (sec. 7.1)
  pc2                   PC2: GanCal#5 in the Jupiter-Ganymede binary through the lane (sec. 7.2)
  member --cell C --index I [--free-bz]
                        a #1065 strong_1065 member (sec. 6); --free-bz is the M2 diagnostic solve
  n1                    N1: the pnb gate-failing zero with the lowest worst ratio above 1
  n2                    N2: identity of phb 26 against its sibling zero (sec. 7.3)
  h1 --cell C --index I H1 test on a spatial member (sec. 4.2)

Every call stops cleanly after --max-wall seconds (default 420) and resumes from its JSON.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from numpy.typing import NDArray

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.pluto_cg import (
    BinaryCG,
    Chain,
    NodeSpec,
    closest_approaches,
    converged,
    flyby_seed_offsets,
    identity_check,
    jacobi,
    noise_ok,
    propagate,
    solve,
)

Arr = NDArray[np.float64]
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "1073_pluto_cg"
LIVE = OUT / "live"
DAY = 86400.0
GM_TOTAL = 975.5  # PRIMARIES["Pluto"], the #998 cell mu
GM_CHARON = 106.1
R_CHARON = 606.0
R_PLUTO = 1188.3  # Nimmo et al. 2017, as scripts/gauntlet_998.py
PLUTO_FLAG_KM = R_PLUTO + 100.0  # Stern et al. 2020 design floor (flag only)
LADDER_MIN = 0.00125
_LOGF: Path | None = None


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')} {msg}"
    print(line, flush=True)
    if _LOGF is not None:
        with _LOGF.open("a") as f:
            f.write(line + "\n")


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _preflight(region: str, n: int) -> None:
    preflight_search(
        task_no=1073,
        region_id=region,
        method=MethodCapability(
            genome="#1065 strong_1065 Pluto-Charon one-working-node members (and lane controls)",
            corrector="pluto_cg closed chain, DOP853 + analytic STM, sigma homotopy",
            capability_tags=frozenset({"ballistic", "continuous-gravity", "cr3bp", "circular"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=n,
        timing_pilot_seconds_per_point=600.0,
    )


def _save(path: Path, rec: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, indent=1, default=float))
    tmp.replace(path)


# ---------------------------------------------------------------------------------------------
# independent re-flies (sec. 6.1 item 6)
# ---------------------------------------------------------------------------------------------


def ias15_rails(m: BinaryCG, sigma: float, x0: Arr, t0: float, t1: float) -> Arr:
    """REBOUND IAS15 with the lane's own field on the rails (a Python force callback)."""
    import rebound

    from cyclerfinder.nbody.pluto_cg import accel_and_gradient

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.integrator.epsilon = 1e-11
    sim.add(m=0.0, x=x0[0], y=x0[1], z=x0[2], vx=x0[3], vy=x0[4], vz=x0[5])
    sim.t = float(t0)

    def forces(ptr: Any) -> None:
        s = ptr.contents
        p = s.particles[0]
        a, _ = accel_and_gradient(m, sigma, np.array([p.x, p.y, p.z]), float(s.t))
        p.ax += a[0]
        p.ay += a[1]
        p.az += a[2]

    sim.additional_forces = forces
    sim.force_is_velocity_dependent = 0
    sim.integrate(float(t1), exact_finish_time=1)
    p = sim.particles[0]
    return np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz])


def live_bodies(m: BinaryCG, x0: Arr, t0: float, t1: float) -> Arr:
    """REBOUND IAS15 with the primary and the secondary as LIVE massive particles (sigma = 1),
    started on their circular states; the craft a test particle. Checks the rails."""
    import rebound

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.integrator.epsilon = 1e-11
    rp, vp, _ = m.primary_state(t0)
    rs, vs, _ = m.secondary_state(t0)
    sim.add(m=m.gm_primary, x=rp[0], y=rp[1], z=rp[2], vx=vp[0], vy=vp[1], vz=vp[2])
    sim.add(m=m.gm_secondary, x=rs[0], y=rs[1], z=rs[2], vx=vs[0], vy=vs[1], vz=vs[2])
    sim.add(m=0.0, x=x0[0], y=x0[1], z=x0[2], vx=x0[3], vy=x0[4], vz=x0[5])
    sim.N_active = 2
    sim.t = float(t0)
    sim.integrate(float(t1), exact_finish_time=1)
    p = sim.particles[2]
    return np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz])


def refly(ch: Chain, z: Arr, rec: dict[str, Any], deadline: float, key: str) -> bool:
    """Both re-flies of every half-arc against DOP853 (resumable per half-arc). PASS < 1e-2 km and
    < 1e-7 km/s. ``live`` only at sigma = 1."""
    out = rec.setdefault(key, {"rails": [], "live": []})
    arcs = ch.arcs(z)
    for i in range(len(out["rails"]), len(arcs)):
        if time.monotonic() > deadline:
            return False
        x0, t0, t1 = arcs[i]
        xd, _ = propagate(ch.model, ch.sigma, x0, t0, t1, stm=False)
        xr = ias15_rails(ch.model, ch.sigma, x0, t0, t1)
        out["rails"].append(
            [float(np.linalg.norm(xr[:3] - xd[:3])), float(np.linalg.norm(xr[3:] - xd[3:]))]
        )
        if ch.sigma == 1.0:
            xl = live_bodies(ch.model, x0, t0, t1)
            out["live"].append(
                [float(np.linalg.norm(xl[:3] - xd[:3])), float(np.linalg.norm(xl[3:] - xd[3:]))]
            )
        log(f"refly arc {i + 1}/{len(arcs)}: rails {out['rails'][-1]} live {out['live'][-1:]}")
    for k in ("rails", "live"):
        if out[k]:
            out[f"{k}_max_dr_km"] = max(a[0] for a in out[k])
            out[f"{k}_max_dv_kms"] = max(a[1] for a in out[k])
            out[f"{k}_pass"] = bool(out[f"{k}_max_dr_km"] < 1e-2 and out[f"{k}_max_dv_kms"] < 1e-7)
    out["done"] = True
    return True


# ---------------------------------------------------------------------------------------------
# PC1: Bruno and Varin 2009, Table 4 (mu = 0.1)
# ---------------------------------------------------------------------------------------------

BV_YAML = REPO / "data" / "sources" / "bruno-varin-2009-family-h-small-mu-tables.yaml"


def pc1_model(mu: float) -> BinaryCG:
    per = 6.38722 * DAY
    a = (GM_TOTAL * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
    # the paper's point masses: no softening radius
    return BinaryCG(GM_TOTAL, mu * GM_TOTAL, a, per, 0.0, 0.0)


def bv_rotating(m: BinaryCG, x: Arr, t: float) -> tuple[float, float, float, float, float]:
    """Bruno-Varin coordinates (P1 at the origin, P2 at (+1, 0); y = canonical momenta) and the
    printed C = -2H + mu, from a barycentric inertial state."""
    mu = m.mass_fraction
    a, n = m.a_rel_km, m.n
    th = m.theta0 + n * t
    c, s = math.cos(th), math.sin(th)
    rot = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    rr = rot @ x[:3] / a
    v_p1 = np.array([0.0, -mu, 0.0])  # P1's inertial velocity in rotating axes, nondimensional
    yy = rot @ x[3:] / (a * n) - v_p1
    x1, x2 = rr[0] + mu, rr[1]
    y1, y2 = yy[0], yy[1]
    r = math.hypot(x1, x2)
    r2 = math.hypot(x1 - 1.0, x2)
    h = 0.5 * (y1 * y1 + y2 * y2) + x2 * y1 - x1 * y2 - 1.0 / r + mu * (1.0 / r + x1 - 1.0 / r2)
    return x1, x2, y1, y2, -2.0 * h + mu


def cmd_pc1(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="pc1")
    ap.add_argument("--row", required=True)
    args = ap.parse_args(argv)
    tab = yaml.safe_load(BV_YAML.read_text())["table_4"]
    row = next(r for r in tab["rows"] if r["k"] == args.row)
    mu = float(tab["mu"])
    _preflight(f"pluto-cg-pc1-bv2009-t4-row{args.row}", 1)
    m = pc1_model(mu)
    a, n = m.a_rel_km, m.n
    x1_0, y2_0 = float(row["x1_0"]), float(row["y2_0"])
    t_norm = float(row["T_norm_printed"])
    period = 2.0 * math.pi * t_norm / (1.0 - mu) / n
    x0 = np.array([a * (x1_0 - mu), 0.0, 0.0, 0.0, a * n * (y2_0 - mu), 0.0])
    ch = Chain(
        m,
        [NodeSpec("free", t_fixed=0.0)],
        period_s=period,
        sigma=1.0,
        free_period=True,
        extra=((0, a * (x1_0 - mu)), (1, 0.0)),
    )
    z0 = np.concatenate([x0, [period]])
    t0 = time.time()
    _, _, info0 = ch.residual_and_jac(z0, want_jac=False)
    z, info, stats = solve(ch, z0, newton_iters=40, trf_nfev=200)
    per = float(z[-1])
    xh, _ = propagate(m, 1.0, z[:6], 0.0, 0.5 * per, stm=False)
    bv0 = bv_rotating(m, z[:6], 0.0)
    bvh = bv_rotating(m, xh, 0.5 * per)
    j0 = jacobi(m, 1.0, z[:6], 0.0)
    jh = jacobi(m, 1.0, xh, 0.5 * per)
    t_norm_got = (1.0 - mu) * per * n / (2.0 * math.pi)
    printed = {
        "x1_0": x1_0,
        "y2_0": y2_0,
        "x1_T2": float(row["x1_T2"]),
        "y2_T2": float(row["y2_T2"]),
        "T_norm": t_norm,
        "C": float(row["C"]),
    }
    got = {
        "x1_0": bv0[0],
        "y2_0": bv0[3],
        "x1_T2": bvh[0],
        "y2_T2": bvh[3],
        "T_norm": t_norm_got,
        "C": bv0[4],
    }
    tol = {"x1_0": 1e-9, "y2_0": 2e-5, "x1_T2": 5e-5, "y2_T2": 5e-5, "T_norm": 3e-5, "C": 5e-5}
    checks = {k: {"printed": printed[k], "got": got[k], "tol": tol[k]} for k in printed}
    for v in checks.values():
        v["ok"] = bool(abs(v["got"] - v["printed"]) <= v["tol"])
    rec: dict[str, Any] = {
        "control": "PC1",
        "source": "Bruno and Varin 2009, Sol. Syst. Res. 43(1):2-25, Table 4 (mu = 0.1)",
        "row": args.row,
        "model": {"gm_total": GM_TOTAL, "mu": mu, "a_rel_km": a, "period_binary_s": m.period_s},
        "seed_dr_km": max(info0["dr"]),
        "seed_dv_kms": max(info0["dv"]),
        "converged_floors": bool(converged(info)),
        "dr_km": max(info["dr"]),
        "dv_kms": max(info["dv"]),
        "stats": stats,
        "z": z.tolist(),
        "period_days": per / DAY,
        "vx0_kms": float(z[3]),
        "y1_at_T2": bvh[2],
        "x2_at_T2": bvh[1],
        "checks": checks,
        "jacobi_rel_drift": abs(jh - j0) / abs(j0),
        "jacobi_szebehely": j0 / (a * n) ** 2,
        "C_minus_mu1mu": float(row["C"]) - mu * (1 - mu),
        "secs": time.time() - t0,
    }
    path = OUT / f"pc1_row{args.row}.json"
    refly(ch, z, rec, time.monotonic() + 1e9, "refly")
    rec["closest"] = closest_approaches(ch, z)
    rec["pass"] = bool(
        rec["converged_floors"]
        and all(v["ok"] for v in checks.values())
        and rec["refly"].get("rails_pass")
        and rec["jacobi_rel_drift"] < 1e-10
    )
    _save(path, rec)
    log(
        f"PC1 row {args.row}: floors={rec['converged_floors']} dr={rec['dr_km']:.2e} "
        f"dv={rec['dv_kms']:.2e} checks={ {k: v['ok'] for k, v in checks.items()} } "
        f"refly={rec['refly'].get('rails_pass')}/{rec['refly'].get('live_pass')} "
        f"J drift {rec['jacobi_rel_drift']:.1e} PASS={rec['pass']}"
    )


# ---------------------------------------------------------------------------------------------
# members and PC2: the per-member procedure (sec. 6)
# ---------------------------------------------------------------------------------------------


@dataclass
class Problem:
    tag: str
    model: BinaryCG
    specs: list[NodeSpec]
    period_s: float
    seed: Arr  # flyby offsets at sigma = 1
    ref: list[dict[str, Any]]
    meta: dict[str, Any] = field(default_factory=dict)

    def chain(self, sigma: float, free_bz: bool = False) -> Chain:
        return Chain(self.model, self.specs, self.period_s, sigma, free_bz=free_bz)

    def seed_z(self, sigma: float, free_bz: bool = False) -> Arr:
        ch = self.chain(1.0)
        z = ch.scale_offsets(self.seed, sigma)
        if not free_bz:
            return z
        out: list[float] = []
        for s, o in zip(ch.specs, ch.offsets, strict=True):
            w = {"flyby": 7, "pin": 4, "free": 6}[s.kind]
            out += list(z[o : o + w])
            if s.kind == "pin":
                out.append(0.0)
        return np.asarray(out)


def problem_from_flybys(
    tag: str,
    model: BinaryCG,
    secondary: str,
    flybys: list[Any],
    period_s: float,
    meta: dict[str, Any],
) -> Problem:
    fl = sorted(flybys, key=lambda f: f.t_s)
    specs: list[NodeSpec] = []
    seed: list[float] = []
    ref: list[dict[str, Any]] = []
    for f in fl:
        vin, vout = np.asarray(f.vinf_in), np.asarray(f.vinf_out)
        if f.body == secondary:
            dr, dv = flyby_seed_offsets(model.gm_secondary, vin, vout)
            specs.append(NodeSpec("flyby"))
            seed += [*dr, *dv, f.t_s]
            ref.append(
                {
                    "kind": "flyby",
                    "vinf_kms": float(np.linalg.norm(vin)),
                    "rp_km": float(np.linalg.norm(dr)),
                    "turn_deg": float(f.turn_deg),
                    "t_days": f.t_s / DAY,
                }
            )
        else:
            specs.append(NodeSpec("pin", f.body))
            seed += [*vin, f.t_s]
            ref.append(
                {"kind": "pin", "vinf_kms": float(np.linalg.norm(vin)), "t_days": f.t_s / DAY}
            )
    return Problem(tag, model, specs, period_s, np.asarray(seed), ref, meta)


def pluto_model(system: Any, moon: str) -> BinaryCG:
    """The CR3BP model of a #1065 barycentric cell (sec. 2)."""
    a_c, per, th0 = system.bodies["Charon"]
    a_rel = a_c / (1.0 - GM_CHARON / system.mu)
    return BinaryCG(
        gm_total=system.mu,
        gm_secondary=GM_CHARON,
        a_rel_km=a_rel,
        period_s=per,
        secondary_radius_km=R_CHARON,
        primary_radius_km=R_PLUTO,
        targets={moon: system.bodies[moon]},
        theta0=th0,
    )


def member_problem(cell: str, index: int) -> Problem:
    r1065 = _load(
        "run_1065_barycentric_enumerate", REPO / "scripts" / "run_1065_barycentric_enumerate.py"
    )
    from cyclerfinder.search.two_working_body import cycle_flybys

    moon = r1065.CELLS[cell]
    system = r1065.bary_system(moon)
    gaunt = json.loads((REPO / "data" / "1065_barycentric" / f"{cell}_gauntlet.json").read_text())
    c = next(c for c in gaunt["candidates"] if c["index"] == index)
    _, cyc = r1065.ENUM.parse_cycle_key(c["key"], system, "Charon", moon)
    x = np.asarray(c["x_days"]) * DAY
    fl = cycle_flybys(system, cyc, x)
    assert fl is not None
    m = pluto_model(system, moon)
    meta = {
        "cell": cell,
        "index": index,
        "key": c["key"],
        "x_days": c["x_days"],
        "worst_ratio": c["worst_ratio"],
        "strong_1065": c.get("strong_1065"),
        "max_out_of_plane": max(
            max(abs(f.vinf_in[2]), abs(f.vinf_out[2])) / float(np.linalg.norm(f.vinf_in))
            for f in fl
        ),
    }
    return problem_from_flybys(f"{cell}_{index}", m, "Charon", list(fl), cyc.period_s, meta)


# --- the per-member state machine --------------------------------------------------------------


def _point(ch: Chain, z: Arr, info: dict[str, Any], scan: bool = True) -> dict[str, Any]:
    q: dict[str, Any] = {
        "sigma": ch.sigma,
        "z": z.tolist(),
        "nodes": ch.describe(z),
        "floor": bool(converged(info)),
        "dr": max(info["dr"]),
        "dv": max(info["dv"]),
    }
    if scan:
        q["closest"] = closest_approaches(ch, z)
    return q


def _impact(q: dict[str, Any]) -> str | None:
    if any(n.get("impact") for n in q["nodes"]):
        return "Charon"
    if "closest" in q and q["closest"]["primary"] <= R_PLUTO:
        return "Pluto"
    return None


def _js(prob: Problem, z: Arr, s: float, free_bz: bool) -> Arr:
    out = []
    for f in (1.0 + 1e-6, 1.0 - 1e-6):
        out.append(prob.chain(s * f, free_bz).residual_and_jac(z, False)[0])
    return (out[0] - out[1]) / (2e-6 * s)


def monitor_step(prob: Problem, rec: dict[str, Any], free_bz: bool) -> bool:
    """#1034/#1025 monitor-variable step (no SVD tangent; the enforced tangent rule)."""
    pts = rec["upper"]
    a, b = pts[-2], pts[-1]
    za, zb = np.asarray(a["z"]), np.asarray(b["z"])
    _, jz, _ = prob.chain(b["sigma"], free_bz).residual_and_jac(zb)
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
            ch = prob.chain(sg, free_bz)
            try:
                r, jzz, info = ch.residual_and_jac(y[:-1])
            except RuntimeError:
                break
            if (converged(info) or noise_ok(info)) and abs(y[j] - target) < 1e-9 * max(
                1.0, abs(target)
            ):
                pts.append(_point(ch, y[:-1], info))
                rec["mon_ratio"] = min(4.0, 1.5 * ratio) if it <= 4 else ratio
                log(f"{prob.tag}: monitor sigma={sg:.6f} ratio {ratio:.3g}")
                return True
            js = _js(prob, y[:-1], sg, free_bz)
            row = np.zeros(len(y))
            row[j] = 1.0
            jac = np.vstack([np.hstack([jzz, js[:, None]]), row[None, :]])
            f = np.concatenate([r, [y[j] - target]])
            if not (np.all(np.isfinite(jac)) and np.all(np.isfinite(f))):
                break
            dc = np.linalg.norm(jac, axis=0)
            dc[dc == 0.0] = 1.0
            try:
                step = np.linalg.lstsq(jac / dc, -f, rcond=None)[0] / dc
            except np.linalg.LinAlgError:
                break
            f0 = float(np.linalg.norm(f))
            for alpha in [0.5**i for i in range(8)]:
                yt = y + alpha * step
                rt = prob.chain(float(yt[-1]), free_bz).residual_and_jac(yt[:-1], False)[0]
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


def two_solutions(prob: Problem, rec: dict[str, Any], free_bz: bool) -> None:
    """#1025 sec. 3 item 6: two solutions at 0.915 and 0.976 sigma_f; supported if both converge and
    some flyby node differs by > 0.02 V_inf_rel... (Pluto scale: 0.4 percent x 5 = 2 percent of
    V_inf) and > 5 percent in r_p / sigma."""
    lower = rec["points"]
    upper = rec["upper"][2:]
    out: dict[str, Any] = {}
    for sg in (round(0.915 * rec["sigma_f"], 7), round(0.976 * rec["sigma_f"], 7)):
        res: dict[str, Any] = {}
        for lab, pts in (("lower", lower), ("upper", upper)):
            if not pts:
                res[lab] = None
                continue
            q = min(pts, key=lambda t: abs(t["sigma"] - sg))
            ch = prob.chain(sg, free_bz)
            z, info, _ = solve(ch, ch.scale_offsets(np.asarray(q["z"]), sg / q["sigma"]), 25, 30)
            res[lab] = ch.describe(z) if converged(info) else None
        row: dict[str, Any] = {
            "lower": res["lower"],
            "upper": res["upper"],
            "fold_supported": False,
        }
        if res["lower"] and res["upper"]:
            dif = [
                (
                    abs(x["vinf_kms"] - y["vinf_kms"]) / x["vinf_kms"],
                    abs(x["rp_over_sigma"] - y["rp_over_sigma"]) / x["rp_over_sigma"],
                )
                for x, y in zip(res["lower"], res["upper"], strict=True)
                if x["kind"] == "flyby"
            ]
            row["rel_d_vinf_rel_d_rp"] = dif
            row["fold_supported"] = any(dv > 0.02 and drp > 0.05 for dv, drp in dif)
        out[str(sg)] = row
        log(f"{prob.tag}: two solutions at {sg}: supported={row['fold_supported']}")
    rec["two_solutions"] = out


def direct_sigma1(prob: Problem, rec: dict[str, Any], free_bz: bool) -> None:
    ch = prob.chain(1.0, free_bz)
    out = {}
    starts = {"seed": (prob.seed_z(1.0, free_bz), 1.0)}
    if rec["points"]:
        starts["lower"] = (np.asarray(rec["points"][-1]["z"]), rec["points"][-1]["sigma"])
    if rec.get("upper"):
        starts["upper"] = (np.asarray(rec["upper"][-1]["z"]), rec["upper"][-1]["sigma"])
    for lab, (z0, s0) in starts.items():
        z, info, _ = solve(ch, ch.scale_offsets(z0, 1.0 / s0), 25, 0)
        out[lab] = {
            "converged": bool(converged(info)),
            "dr": max(info["dr"]),
            "nodes": ch.describe(z) if converged(info) else None,
        }
    rec["direct_sigma1"] = out


def run_problem(
    prob: Problem, path: Path, deadline: float, free_bz: bool = False
) -> dict[str, Any]:
    rec: dict[str, Any] = (
        json.loads(path.read_text())
        if path.exists()
        else {"tag": prob.tag, "meta": prob.meta, "ref": prob.ref, "points": [], "stage": "start"}
    )
    rec["free_bz"] = free_bz

    def finish(outcome: str, **kw: Any) -> None:
        rec.update(stage="done", outcome=outcome, **kw)
        log(f"{prob.tag}: {outcome} {kw}")

    while time.monotonic() < deadline and rec["stage"] != "done":
        st = rec["stage"]
        t_stage = time.time()
        if st == "start":
            tries = rec.setdefault("start_tries", [])
            for sg in (0.02, 0.005, 0.00125)[len(tries) :]:
                ch = prob.chain(sg, free_bz)
                z, info, stats = solve(ch, prob.seed_z(sg, free_bz), 40, 200)
                tries.append(
                    {
                        "sigma": sg,
                        "converged": bool(converged(info)),
                        "dr": max(info["dr"]),
                        "dv": max(info["dv"]),
                        "gauge": max((abs(g) for g in info["gauge"]), default=0.0),
                        "stats": stats,
                        "secs": time.time() - t_stage,
                        "nodes": ch.describe(z),
                        "z": z.tolist(),
                    }
                )
                log(
                    f"{prob.tag}: start sigma={sg}: converged={converged(info)} "
                    f"dr={max(info['dr']):.3e} "
                    f"dv={max(info['dv']):.3e} ({time.time() - t_stage:.0f}s)"
                )
                if converged(info):
                    rec["points"] = [_point(ch, z, info)]
                    rec["ladder"] = [rec["points"][0]]
                    rec["stage"] = "ladder"
                    break
                _save(path, rec)
                if time.monotonic() > deadline:
                    break
            else:
                finish("NO START", reason="no convergence at sigma 0.02, 0.005, 0.00125")
        elif st == "ladder":
            lad = rec["ladder"]
            last = lad[-1]
            if last["sigma"] <= LADDER_MIN * 1.0001:
                rec["stage"] = "identity"
                continue
            sg = last["sigma"] / 2.0
            ch = prob.chain(sg, free_bz)
            z, info, _ = solve(
                ch, ch.scale_offsets(np.asarray(last["z"]), 0.5), 40, 100, accept_noise=True
            )
            if converged(info) or noise_ok(info):
                lad.append(_point(ch, z, info, scan=False))
                log(f"{prob.tag}: ladder sigma={sg:g} ok ({time.time() - t_stage:.0f}s)")
            else:
                rec["ladder_stop"] = sg
                log(f"{prob.tag}: ladder sigma={sg:g} FAILED")
                rec["stage"] = "identity"
        elif st == "identity":
            idc = identity_check(rec["ladder"], prob.ref)
            rec["identity"] = idc
            log(
                f"{prob.tag}: identity {idc['verdict']} ok={idc.get('ok')} "
                f"pair={idc.get('sigma_pair')}"
            )
            if idc["verdict"] == "FAIL":
                finish("IDENTITY FAIL")
            else:
                rec["stage"], rec["fac"] = "up", 1.2
                imp = _impact(rec["points"][0])
                if imp:
                    finish("IMPACT", sigma_i=rec["points"][0]["sigma"], body=imp)
        elif st == "up":
            pts = rec["points"]
            s1, z1 = pts[-1]["sigma"], np.asarray(pts[-1]["z"])
            if s1 >= 1.0:
                rec["stage"] = "end"
                continue
            fac = rec["fac"]
            sg = min(1.0, s1 * fac)
            if sg - s1 < 1e-5:
                rec["sigma_f"] = s1
                rec["stage"] = "fold"
                rec["upper"] = [pts[-2], pts[-1]] if len(pts) >= 2 else [pts[-1]]
                log(f"{prob.tag}: natural step below 1e-5 at sigma {s1:.6f}; monitor")
                if len(rec["upper"]) < 2:
                    finish("NUMERICAL STOP", reason="stop at the first point")
                continue
            ch = prob.chain(sg, free_bz)
            if len(pts) >= 2:
                s0, z0 = pts[-2]["sigma"], np.asarray(pts[-2]["z"])
                zs = z1 + (z1 - z0) * math.log(sg / s1) / math.log(s1 / s0)
            else:
                zs = ch.scale_offsets(z1, sg / s1)
            z, info, _ = solve(ch, zs, 25, 15, accept_noise=sg < 1.0)
            if converged(info) or (sg < 1.0 and noise_ok(info)):
                q = _point(ch, z, info)
                imp = _impact(q)
                if imp:
                    pts.append(q)
                    finish("IMPACT", sigma_i=sg, body=imp)
                else:
                    pts.append(q)
                    rec["fac"] = min(1.2, 1.0 + 2.0 * (fac - 1.0))
                    vs = [round(n["vinf_kms"], 4) for n in q["nodes"]]
                    log(
                        f"{prob.tag}: sigma={sg:.6g} ok vinf {vs} "
                        f"pluto {q['closest']['primary']:.0f} km "
                        f"({time.time() - t_stage:.0f}s)"
                    )
            else:
                rec["fac"] = 1.0 + 0.5 * (fac - 1.0)
                log(f"{prob.tag}: sigma={sg:.6g} failed; factor {rec['fac']:.6g}")
        elif st == "fold":
            up = rec["upper"]
            sf = rec["sigma_f"]
            if up[-1]["sigma"] > sf + 1e-4:
                rec["points"].append(up[-1])
                rec.setdefault("passed_stops", []).append(sf)
                rec.pop("upper")
                rec.pop("mon_ratio", None)
                rec["stage"], rec["fac"] = "up", 1.01
                log(f"{prob.tag}: monitor passed the stop at {sf:.6f}; resume")
            elif up[-1]["sigma"] < 0.91 * sf:
                rec["stage"] = "fold_checks"
            elif _impact(up[-1]):
                finish("IMPACT", sigma_i=up[-1]["sigma"], body=_impact(up[-1]), on="monitor branch")
            elif not monitor_step(prob, rec, free_bz):
                finish(
                    "NUMERICAL STOP",
                    reason=f"natural stop at {sf:.6g}; monitor: {rec['monitor_stop']}",
                )
        elif st == "fold_checks":
            two_solutions(prob, rec, free_bz)
            direct_sigma1(prob, rec, free_bz)
            if all(v["fold_supported"] for v in rec["two_solutions"].values()):
                finish("FOLDS", sigma_f=rec["sigma_f"])
            else:
                finish(
                    "NUMERICAL STOP", reason="turn passed but the two-solution check not supported"
                )
        elif st == "end":
            q = rec["points"][-1]
            ch = prob.chain(1.0, free_bz)
            if refly(ch, np.asarray(q["z"]), rec, deadline, "refly"):
                rf = rec["refly"]
                rec["full_mass"] = q["nodes"]
                rec["full_mass_closest"] = q.get("closest")
                ok = bool(rf.get("rails_pass") and rf.get("live_pass"))
                verdict = rec["identity"]["verdict"]
                finish(
                    ("EXISTS" if verdict == "PASS" else "EXISTS-SIBLING")
                    if ok
                    else "EXISTS (re-fly FAIL)",
                    refly_pass=ok,
                )
        _save(path, rec)
    _save(path, rec)
    return rec


def cmd_member(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="member")
    ap.add_argument("--cell", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--free-bz", action="store_true")
    ap.add_argument("--max-wall", type=float, default=420.0)
    args = ap.parse_args(argv)
    global _LOGF
    tag = f"{args.cell}_{args.index}" + ("_bz" if args.free_bz else "")
    LIVE.mkdir(parents=True, exist_ok=True)
    _LOGF = LIVE / f"{tag}.live.log"
    _preflight(f"pluto-cg-member-{tag}", 1)
    prob = member_problem(args.cell, args.index)
    rec = run_problem(prob, OUT / f"{tag}.json", time.monotonic() + args.max_wall, args.free_bz)
    log(f"{tag}: stage {rec['stage']} outcome {rec.get('outcome')}")


def main() -> None:
    cmds = {"pc1": cmd_pc1, "member": cmd_member}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        raise SystemExit(f"usage: {Path(__file__).name} {{{','.join(cmds)}}} ...")
    cmds[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
