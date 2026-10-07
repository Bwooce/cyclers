"""#968 second published control: GanEur#43 (Russell & Strange 2009), note sec. 5.3.

Same model and corrector as ``run_968_control.py`` (corrector A: forward-backward multiple
shooting, DOP853 + analytic STM, Ganymede the only massive moon, the massless target's centre hit
pinned, rotation-aware periodic wrap, periapsis gauge at the Ganymede node), written for a general
node list. GanEur#43: one Ganymede flyby node B and one massless Europa node E per period
(2 G-E synodic periods). Lane corrector B (``jovian_defect_residual``, fixed wrap) is seeded from
corrector A's solution. Cross-check: REBOUND IAS15 re-fly of every half-arc.

    uv run python scripts/run_968_control2.py --stage a      # resumable, a_state.json
    uv run python scripts/run_968_control2.py --stage check
    uv run python scripts/run_968_control2.py --stage lane
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import (
    _W_VEL,
    JovianRailsCache,
    jovian_defect_residual,
    periapsis_node,
)
from cyclerfinder.nbody.jovian_stm import jovian_stm_jacobian
from cyclerfinder.nbody.shooter import ShootingSeed, _seed_with_states, _states_to_x, _x_to_states
from cyclerfinder.search.two_working_body import cycle_flybys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "968_control2"
DAY = 86400.0
KEY = "k2|LGanymede>Europa/1h|LEuropa>Ganymede/1l"
TARGET = "Europa"
GAN = "Ganymede"
Arr = NDArray[np.float64]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CTL = _load("run_968_control", ROOT / "scripts" / "run_968_control.py")


def _log(msg: str) -> None:
    (OUT / "live").mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "live" / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


def build() -> tuple[Any, Arr]:
    """Problem (from run_968_control) and the patched-conic seed z = [x_B, t_B, v_E, t_E]."""
    enum = CTL._load_enum()
    circ, a, b = enum.cell_system("ge1")  # R-S Table 2, Europa massless
    cands = {
        c["key"]: c
        for c in json.loads((ROOT / "data/943_cell_ge_gauntlet.json").read_text())["candidates"]
    }
    _, cyc = enum.parse_cycle_key(KEY, circ, a, b)
    x = np.asarray(cands[KEY]["x_days"]) * DAY
    fl = cycle_flybys(circ, cyc, x)
    assert fl is not None and [f.body for f in fl] == [TARGET, GAN]
    fe, fg = fl
    period = cyc.period_s
    q = CTL.rot_z(2.0 * math.pi * period / circ.period_s(GAN))
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    # Node B: the Ganymede flyby one period earlier (so B precedes E in time).
    tg = fg.t_s - period
    qi = q.T
    rb, vb, _ = periapsis_node(GAN, tg, qi @ fg.vinf_in, qi @ fg.vinf_out, circ)  # type: ignore[arg-type]
    _, ve = circ.state(TARGET, fe.t_s)
    seed = np.concatenate([rb, vb, [tg], ve + fe.vinf_in, [fe.t_s]]).astype(np.float64)
    mus = {GAN: SATELLITES[GAN].mu_km3_s2}
    surf = {GAN: SATELLITES[GAN].radius_eq_km}
    return CTL.Problem(circ, period, q6, seed, mus, surf), seed


N = 11  # z: B state 0:6, tB 6, E velocity 7:10, tE 10


def nodes(p: Any, z: Arr) -> list[tuple[Arr, Arr, float, Arr]]:
    dxb = np.zeros((6, N))
    dxb[:, 0:6] = np.eye(6)
    dtb = np.zeros(N)
    dtb[6] = 1.0
    te = float(z[10])
    re, ve = p.circ.state(TARGET, te)
    dxe = np.zeros((6, N))
    dxe[3:6, 7:10] = np.eye(3)
    dxe[0:3, 10] = ve
    dte = np.zeros(N)
    dte[10] = 1.0
    return [
        (z[0:6].copy(), dxb, float(z[6]), dtb),
        (np.concatenate([re, z[7:10]]), dxe, te, dte),
        (p.q6 @ z[0:6], p.q6 @ dxb, float(z[6]) + p.period_s, dtb.copy()),
    ]


def residual_and_jac(p: Any, z: Arr, want_jac: bool = True) -> tuple[Arr, Arr, dict[str, Any]]:
    nd = nodes(p, z)
    res: list[float] = []
    rows: list[Arr] = []
    info: dict[str, Any] = {"match_dr_km": [], "match_dv_kms": [], "gauge": []}
    for k in range(2):
        xp, dxp, tp, dtp = nd[k]
        xn, dxn, tn, dtn = nd[k + 1]
        tm = 0.5 * (tp + tn)
        xf, pf = CTL.arc(p, xp, tp, tm)
        xb, pb = CTL.arc(p, xn, tn, tm)
        d = xf - xb
        info["match_dr_km"].append(float(np.linalg.norm(d[:3])))
        info["match_dv_kms"].append(float(np.linalg.norm(d[3:])))
        res.extend(CTL.W * d)
        if want_jac:
            ff = CTL.accel(p, xf, tm)
            fb = CTL.accel(p, xb, tm)
            dxf = (
                pf @ dxp
                + np.outer(-pf @ CTL.accel(p, xp, tp) + 0.5 * ff, dtp)
                + np.outer(0.5 * ff, dtn)
            )
            dxb_ = (
                pb @ dxn
                + np.outer(-pb @ CTL.accel(p, xn, tn) + 0.5 * fb, dtn)
                + np.outer(0.5 * fb, dtp)
            )
            rows.append(CTL.W[:, None] * (dxf - dxb_))
    x, dx, t, dt = nd[0]
    info["gauge"].append(CTL.gauge(p, x, t))
    res.append(CTL.W_GAUGE * info["gauge"][0])
    if want_jac:
        gx, gt = CTL.gauge_grad(p, x, t)
        rows.append(CTL.W_GAUGE * (gx @ dx + gt * dt)[None, :])
    return np.asarray(res), (np.vstack(rows) if want_jac else np.zeros((0, N))), info


def stage_a(max_nfev: int) -> None:
    p, seed = build()
    path = OUT / "a_state.json"
    z = np.asarray(json.loads(path.read_text())["z"]) if path.exists() else seed.copy()
    if not path.exists():
        _, jac, _ = residual_and_jac(p, z)
        errs = {}
        for j in range(N):
            h = 1e-6 * max(1.0, abs(float(z[j])))
            zp, zm = z.copy(), z.copy()
            zp[j] += h
            zm[j] -= h
            fd = (residual_and_jac(p, zp, False)[0] - residual_and_jac(p, zm, False)[0]) / (2 * h)
            errs[f"col{j}"] = float(np.linalg.norm(fd - jac[:, j]) / max(np.linalg.norm(fd), 1e-30))
        _log(f"stage a: Jacobian FD check {errs}")
    cache: dict[str, Any] = {}

    def fun(zz: Arr) -> Arr:
        r, j, info = residual_and_jac(p, zz)
        cache["z"], cache["j"] = zz.copy(), j
        _log(
            f"stage a |r| {np.linalg.norm(r):.4e} dr {max(info['match_dr_km']):.3e} km "
            f"dv {max(info['match_dv_kms']):.3e} km/s gauge {abs(info['gauge'][0]):.2e}"
        )
        return r

    def jac(zz: Arr) -> Arr:
        if "z" in cache and np.array_equal(cache["z"], zz):
            return np.asarray(cache["j"])
        return residual_and_jac(p, zz)[1]

    sol = least_squares(
        fun,
        z,
        jac=jac,
        method="trf",
        x_scale="jac",
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=max_nfev,
    )
    r, _, info = residual_and_jac(p, np.asarray(sol.x), False)
    ok = (
        max(info["match_dr_km"]) < 1e-3
        and max(info["match_dv_kms"]) < 1e-6
        and abs(info["gauge"][0]) < 1e-9
    )
    path.write_text(
        json.dumps(
            {"z": np.asarray(sol.x).tolist(), "norm": float(np.linalg.norm(r)), "converged": ok}
        )
    )
    _log(f"stage a: end |r| {np.linalg.norm(r):.4e} converged={ok} nfev {sol.nfev}")


def stage_check() -> None:
    from scipy.integrate import solve_ivp
    from scipy.optimize import minimize_scalar

    p, _ = build()
    z = np.asarray(json.loads((OUT / "a_state.json").read_text())["z"])
    _, _, info = residual_and_jac(p, z, False)
    nd = nodes(p, z)
    out: dict[str, Any] = {
        "note": "docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 5.3",
        **info,
    }
    re, ve = p.circ.state(TARGET, nd[1][2])
    out["europa_hit_km"] = float(np.linalg.norm(nd[1][0][:3] - re))
    out["europa_speed_kms"] = float(np.linalg.norm(nd[1][0][3:] - ve))
    vb, db = CTL.vinf_gan(p, nd[0][0], nd[0][2])
    out["ganymede"] = {"vinf_kms": vb, "alt_rs_km": db - CTL.R_GAN_RS_KM, "t_days": nd[0][2] / DAY}
    t0 = nd[0][2]

    def rhs(t: float, y: Arr) -> Arr:
        return np.asarray(CTL.accel(p, y, t))

    sol = solve_ivp(
        rhs,
        (t0, t0 + p.period_s),
        nd[0][0],
        method="DOP853",
        rtol=CTL.RTOL,
        atol=CTL.ATOL,
        dense_output=True,
    )
    ts = np.arange(t0, t0 + p.period_s, 0.005 * DAY)
    ys = sol.sol(ts)
    rr = np.linalg.norm(ys[:3], axis=0)
    dg = np.array(
        [np.linalg.norm(ys[:3, i] - p.circ.state(GAN, float(t))[0]) for i, t in enumerate(ts)]
    )
    hill = 31_715.0
    mins = [
        {"t_days": float(ts[i]) / DAY, "dist_km": float(dg[i])}
        for i in range(1, len(ts) - 1)
        if dg[i] <= dg[i - 1] and dg[i] <= dg[i + 1]
    ]
    out["ganymede_local_minima"] = mins
    out["ganymede_inside_hill_other_than_node"] = [m for m in mins if m["dist_km"] < hill]
    i0, i1 = int(np.argmin(rr)), int(np.argmax(rr))
    out["r_min_km"] = float(
        minimize_scalar(
            lambda t: float(np.linalg.norm(sol.sol(t)[:3])),
            bounds=(ts[max(i0 - 1, 0)], ts[min(i0 + 1, len(ts) - 1)]),
            method="bounded",
        ).fun
    )
    out["r_max_km"] = -float(
        minimize_scalar(
            lambda t: -float(np.linalg.norm(sol.sol(t)[:3])),
            bounds=(ts[max(i1 - 1, 0)], ts[min(i1 + 1, len(ts) - 1)]),
            method="bounded",
        ).fun
    )
    xcheck = []
    for k in range(2):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for lab, x0, ts0 in (("fwd", xp, tp), ("bwd", xn, tn)):
            xd, _ = CTL.arc(p, x0, ts0, tm)
            tw = time.monotonic()
            xi, status = CTL.ias15(p, x0, ts0, tm)
            row: dict[str, Any] = {
                "leg": k,
                "dir": lab,
                "status": status,
                "wall_s": time.monotonic() - tw,
            }
            if xi is not None:
                row["dr_km"] = float(np.linalg.norm(xi[:3] - xd[:3]))
                row["dv_kms"] = float(np.linalg.norm(xi[3:] - xd[3:]))
            xcheck.append(row)
    out["ias15_halfarc"] = xcheck
    xi, status = CTL.ias15(p, nd[0][0], t0, t0 + p.period_s)
    tgt = p.q6 @ nd[0][0]
    out["ias15_full_period"] = (
        {
            "dr_km": float(np.linalg.norm(xi[:3] - tgt[:3])),
            "dv_kms": float(np.linalg.norm(xi[3:] - tgt[3:])),
        }
        if xi is not None
        else {"status": status}
    )
    (OUT / "a_check.json").write_text(json.dumps(out, indent=1))
    _log(
        f"stage check: {json.dumps({k: v for k, v in out.items() if k != 'ganymede_local_minima'})}"
    )


def stage_lane(max_nfev: int) -> None:
    """Criterion 5: the lane's residual (fixed wrap), seeded from corrector A's solution."""
    p, _ = build()
    z = np.asarray(json.loads((OUT / "a_state.json").read_text())["z"])
    nd = nodes(p, z)
    zero = np.zeros(3)
    seed = ShootingSeed(
        node_states=[np.asarray(q[0]) for q in nd],
        epochs=[float(q[2]) for q in nd],
        tofs=[(nd[i + 1][2] - nd[i][2]) / DAY for i in range(2)],
        sequence=(GAN, TARGET, GAN),
        slack_leg=0,
        period_days=p.period_s / DAY,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )
    moons = (GAN,)
    cache = JovianRailsCache(moons, p.circ, min(seed.epochs), max(seed.epochs))  # type: ignore[arg-type]

    def res(xv: Arr) -> Arr:
        r = jovian_defect_residual(
            _seed_with_states(seed, _x_to_states(xv, 3)),
            ephem=p.circ,
            cache=cache,
            moons=moons,
            max_wall_sec=3000.0,  # type: ignore[arg-type]
        )
        if np.any(r == 1e9):
            _log("lane: timeout/non-finite sentinel in a leg")
        _log(f"lane |r| {np.linalg.norm(r):.4e}")
        return r

    def jac(xv: Arr) -> Arr:
        return jovian_stm_jacobian(seed, xv, ephem=p.circ, moons=moons)  # type: ignore[arg-type]

    x0 = _states_to_x(seed.node_states)
    r0 = res(x0)
    sol = least_squares(res, x0, jac=jac, method="trf", x_scale="jac", max_nfev=max_nfev)
    rf = res(np.asarray(sol.x))
    legs = rf[:12].reshape(2, 6)
    states = _x_to_states(np.asarray(sol.x), 3)
    vb, db = CTL.vinf_gan(p, np.asarray(states[0]), seed.epochs[0])
    out = {
        "seed_norm": float(np.linalg.norm(r0)),
        "norm": float(np.linalg.norm(rf)),
        "leg_dr_km": np.linalg.norm(legs[:, :3], axis=1).tolist(),
        "leg_dv_kms": (np.linalg.norm(legs[:, 3:], axis=1) / _W_VEL).tolist(),
        "wrap_dr_km": float(np.linalg.norm(rf[-6:-3])),
        "wrap_dv_kms": float(np.linalg.norm(rf[-3:]) / _W_VEL),
        "ganymede_vinf_kms": vb,
        "ganymede_node_dist_km": db,
    }
    (OUT / "lane_state.json").write_text(json.dumps(out, indent=1))
    _log(f"lane end: {json.dumps(out)}")


def stage_gm(s_list: list[float], max_nfev: int) -> None:
    """Amendment 7: Ganymede-GM continuation of the pinned GanEur#43 solution (s_G from 1 down)."""
    p, _ = build()
    mu0, r0 = p.mus[GAN], p.surf[GAN]
    path = OUT / "gm.json"
    rec = json.loads(path.read_text()) if path.exists() else {"points": []}
    good = [q for q in rec["points"] if q["converged"]]
    hist = [(float(q["s"]), np.asarray(q["z"])) for q in good[-2:]]
    if not hist:
        hist = [(1.0, np.asarray(json.loads((OUT / "a_state.json").read_text())["z"]))]
    for s_g in s_list:
        if len(hist) == 2:
            (s0, z0), (s1, z1) = hist
            zs = z1 + (z1 - z0) * (math.log(s_g) - math.log(s1)) / (math.log(s1) - math.log(s0))
        else:
            s1, z1 = hist[-1]
            zs = z1.copy()
            rg, _ = p.circ.state(GAN, float(zs[6]))
            zs[0:3] = rg + (zs[0:3] - rg) * (s_g / s1)
        p.mus[GAN] = s_g * mu0
        p.surf[GAN] = s_g * r0
        sol = least_squares(
            lambda zz: residual_and_jac(p, zz, False)[0],
            zs,
            jac=lambda zz: residual_and_jac(p, zz)[1],
            method="trf",
            x_scale="jac",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev,
        )
        zn = np.asarray(sol.x)
        r, _, info = residual_and_jac(p, zn, False)
        ok = (
            max(info["match_dr_km"]) < 1e-3
            and max(info["match_dv_kms"]) < 1e-6
            and abs(info["gauge"][0]) < 1e-9
        )
        pt: dict[str, Any] = {
            "s": s_g,
            "converged": ok,
            "norm": float(np.linalg.norm(r)),
            "z": zn.tolist(),
        }
        if ok:
            nd = nodes(p, zn)
            v, d = CTL.vinf_gan(p, nd[0][0], nd[0][2])
            _, ve = p.circ.state(TARGET, nd[1][2])
            pt.update(
                {
                    "vinf_G": v,
                    "rp_over_s_km": d / s_g,
                    "europa_speed_kms": float(np.linalg.norm(nd[1][0][3:] - ve)),
                }
            )
            hist = [*hist[-1:], (s_g, zn)]
        rec["points"].append(pt)
        path.write_text(json.dumps(rec))
        _log(
            f"gm s={s_g:.4g}: converged={ok} |r| {pt['norm']:.2e} "
            f"{pt.get('vinf_G')} {pt.get('rp_over_s_km')}"
        )
        if not ok:
            break


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="ganeur43-continuous-ideal-positive-control",
        method=MethodCapability(
            genome="R-S 2009 GanEur#43 (one structure), circular coplanar, Europa massless",
            corrector="forward-backward multiple shooting, DOP853 + analytic STM",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["a", "check", "lane", "gm"], required=True)
    ap.add_argument("--max-nfev", type=int, default=60)
    ap.add_argument("--s-list", default="")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "a":
        stage_a(args.max_nfev)
    elif args.stage == "check":
        stage_check()
    elif args.stage == "gm":
        stage_gm([float(v) for v in args.s_list.split(",")], args.max_nfev)
    else:
        stage_lane(args.max_nfev)


if __name__ == "__main__":
    main()
