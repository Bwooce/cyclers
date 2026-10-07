"""#1004: gc-2 continued down in Callisto's GM in continuous gravity (pre-registration:
docs/notes/2026-10-07-1004-gc2-gancal5-callisto-gm-continuation.md).

R-S 2009 circular model, Ganymede and Callisto point masses (``jovian_stm.propagate_with_stm`` with
GM/radius overrides). Forward-backward multiple shooting with four massive-body nodes
B (G) -> C1 (C) -> C2 (C) -> A (G) -> B' = (Q x_B, t_B + T); periapsis gauge at each node.

    uv run python scripts/run_1004_gc2.py --stage sigma --list 0.05,0.07,...   # joint GM scale up
    uv run python scripts/run_1004_gc2.py --stage cont --variant S --points 20  # s_C down
    uv run python scripts/run_1004_gc2.py --stage cont --variant P --points 20

Checkpoints in data/1004_gc2/; live log in data/1004_gc2/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import periapsis_node
from cyclerfinder.nbody.jovian_stm import _accel_and_gradient, propagate_with_stm
from cyclerfinder.search.two_working_body import CircularSystem, cycle_flybys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1004_gc2"
DAY = 86400.0
KEY = "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s"
GAN, CAL = "Ganymede", "Callisto"
MOONS = (GAN, CAL)
NODE_MOONS = (GAN, CAL, CAL, GAN)
NZ = 28
W = np.array([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
W_GAUGE = 1e4
RTOL, ATOL = 1e-12, 1e-10
R_CAL_RS_KM = 2408.0
R_GAN_RS_KM = 2634.0
SOI_CAL_KM = 37_681.0
Arr = NDArray[np.float64]


def _log(msg: str) -> None:
    (OUT / "live").mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with (OUT / "live" / "runlog.txt").open("a") as fh:
        fh.write(line + "\n")


def _enum() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_942_enumerate", ROOT / "scripts" / "run_942_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_942_enumerate"] = mod
    spec.loader.exec_module(mod)
    return mod


def rot_z(a: float) -> Arr:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


@dataclass
class P:
    circ: CircularSystem
    period: float
    q6: Arr
    seed: Arr
    mus: dict[str, float] = field(default_factory=dict)
    surf: dict[str, float] = field(default_factory=dict)

    def scale(self, s_g: float, s_c: float, rad_g: float, rad_c: float) -> None:
        self.mus = {GAN: s_g * SATELLITES[GAN].mu_km3_s2, CAL: s_c * SATELLITES[CAL].mu_km3_s2}
        self.surf = {GAN: rad_g, CAL: rad_c}


def build() -> P:
    enum = _enum()
    circ, a, b = enum.cell_system("gc")
    cands = {
        c["key"]: c
        for c in json.loads((ROOT / "data/943_cell_gc_gauntlet.json").read_text())["candidates"]
    }
    _, cyc = enum.parse_cycle_key(KEY, circ, a, b)
    fl = cycle_flybys(circ, cyc, np.asarray(cands[KEY]["x_days"]) * DAY)
    assert fl is not None
    fl = sorted(fl, key=lambda f: f.t_s)
    assert [f.body for f in fl] == list(NODE_MOONS), [f.body for f in fl]
    period = cyc.period_s
    q = rot_z(2.0 * math.pi * period / circ.period_s(GAN))
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    seed = []
    for f in fl:
        r, v, _ = periapsis_node(f.body, f.t_s, f.vinf_in, f.vinf_out, circ)  # type: ignore[arg-type]
        seed.extend([*r, *v, f.t_s])
    p = P(circ, period, q6, np.asarray(seed, dtype=np.float64))
    p.scale(1.0, 1.0, SATELLITES[GAN].radius_eq_km, SATELLITES[CAL].radius_eq_km)
    return p


def nodes(p: P, z: Arr) -> list[tuple[Arr, Arr, float, Arr]]:
    out = []
    for k in range(4):
        dx = np.zeros((6, NZ))
        dx[:, 7 * k : 7 * k + 6] = np.eye(6)
        dt = np.zeros(NZ)
        dt[7 * k + 6] = 1.0
        out.append((z[7 * k : 7 * k + 6].copy(), dx, float(z[7 * k + 6]), dt))
    x0, dx0, t0, dt0 = out[0]
    out.append((p.q6 @ x0, p.q6 @ dx0, t0 + p.period, dt0.copy()))
    return out


def accel(p: P, x: Arr, t: float) -> Arr:
    a, _ = _accel_and_gradient(x[:3], t, ephem=p.circ, moons=MOONS, mus=p.mus, surf=p.surf)  # type: ignore[arg-type]
    return np.concatenate([x[3:], a])


def arc(p: P, x: Arr, t0: float, t1: float) -> tuple[Arr, Arr]:
    rf, vf, phi = propagate_with_stm(
        x[:3],
        x[3:],
        t0,
        t1,
        ephem=p.circ,  # type: ignore[arg-type]
        moons=MOONS,
        rtol=RTOL,
        atol=ATOL,
        mu_overrides=p.mus,
        radius_overrides=p.surf,
    )
    return np.concatenate([rf, vf]), phi


def gauge_and_grad(p: P, moon: str, x: Arr, t: float) -> tuple[float, Arr, float]:
    rm, vm = p.circ.state(moon, t)
    n_m = 2.0 * math.pi / p.circ.period_s(moon)
    am = -(n_m**2) * rm
    dr, dv = x[:3] - rm, x[3:] - vm
    a, b = float(np.linalg.norm(dr)), float(np.linalg.norm(dv))
    g = float(dr @ dv) / (a * b)
    g_dr = dv / (a * b) - g * dr / a**2
    g_dv = dr / (a * b) - g * dv / b**2
    return g, np.concatenate([g_dr, g_dv]), float(-(g_dr @ vm) - (g_dv @ am))


def residual_and_jac(p: P, z: Arr, want_jac: bool = True) -> tuple[Arr, Arr, dict[str, Any]]:
    nd = nodes(p, z)
    res: list[float] = []
    rows: list[Arr] = []
    info: dict[str, Any] = {"dr": [], "dv": [], "gauge": []}
    for k in range(4):
        xp, dxp, tp, dtp = nd[k]
        xn, dxn, tn, dtn = nd[k + 1]
        tm = 0.5 * (tp + tn)
        xf, pf = arc(p, xp, tp, tm)
        xb, pb = arc(p, xn, tn, tm)
        d = xf - xb
        info["dr"].append(float(np.linalg.norm(d[:3])))
        info["dv"].append(float(np.linalg.norm(d[3:])))
        res.extend(W * d)
        if want_jac:
            ff, fb = accel(p, xf, tm), accel(p, xb, tm)
            dxf = (
                pf @ dxp
                + np.outer(-pf @ accel(p, xp, tp) + 0.5 * ff, dtp)
                + np.outer(0.5 * ff, dtn)
            )
            dxb = (
                pb @ dxn
                + np.outer(-pb @ accel(p, xn, tn) + 0.5 * fb, dtn)
                + np.outer(0.5 * fb, dtp)
            )
            rows.append(W[:, None] * (dxf - dxb))
    for k in range(4):
        x, dx, t, dt = nd[k]
        g, gx, gt = gauge_and_grad(p, NODE_MOONS[k], x, t)
        info["gauge"].append(g)
        res.append(W_GAUGE * g)
        if want_jac:
            rows.append(W_GAUGE * (gx @ dx + gt * dt)[None, :])
    return np.asarray(res), (np.vstack(rows) if want_jac else np.zeros((0, NZ))), info


def converged(info: dict[str, Any]) -> bool:
    return (
        max(info["dr"]) < 1e-3
        and max(info["dv"]) < 1e-6
        and max(abs(g) for g in info["gauge"]) < 1e-9
    )


def describe(p: P, z: Arr) -> dict[str, Any]:
    nd = nodes(p, z)
    out: dict[str, Any] = {"nodes": []}
    for k in range(4):
        moon = NODE_MOONS[k]
        x, _, t, _ = nd[k]
        rm, vm = p.circ.state(moon, t)
        d = float(np.linalg.norm(x[:3] - rm))
        v2 = float(np.linalg.norm(x[3:] - vm)) ** 2
        mu = p.mus[moon]
        vinf2 = v2 - 2.0 * mu / d
        vinf = math.sqrt(max(vinf2, 0.0))
        e = 1.0 + d * vinf2 / mu if mu > 0 else float("inf")
        turn = math.degrees(2.0 * math.asin(1.0 / e)) if e > 1.0 else float("nan")
        out["nodes"].append(
            {"moon": moon, "t_days": t / DAY, "vinf_kms": vinf, "rp_km": d, "turn_deg": turn}
        )
    # C1 -> C2 arc osculating Jupiter-centred elements at mid-time
    x1, _, t1, _ = nd[1]
    tm = 0.5 * (t1 + nd[2][2])
    xm, _ = arc(p, x1, t1, tm)
    mu_j = p.circ.mu
    r, v = xm[:3], xm[3:]
    rn = float(np.linalg.norm(r))
    a = 1.0 / (2.0 / rn - float(v @ v) / mu_j)
    ev = ((float(v @ v) - mu_j / rn) * r - float(r @ v) * v) / mu_j
    out["c1c2_mid_a_km"] = a
    out["c1c2_mid_e"] = float(np.linalg.norm(ev))
    return out


class _DoneError(Exception):
    """Raised inside the residual once the floors are met (trf would otherwise keep polishing)."""


def solve(p: P, z0: Arr, max_nfev: int, tag: str) -> tuple[Arr, dict[str, Any]]:
    cache: dict[str, Any] = {}

    def fun(z: Arr) -> Arr:
        r, j, info = residual_and_jac(p, z)
        cache["z"], cache["j"] = z.copy(), j
        _log(f"{tag} |r| {np.linalg.norm(r):.3e} dr {max(info['dr']):.2e} dv {max(info['dv']):.2e}")
        if converged(info):
            cache["done"] = z.copy()
            raise _DoneError
        return r

    def jac(z: Arr) -> Arr:
        if "z" in cache and np.array_equal(cache["z"], z):
            return np.asarray(cache["j"])
        return residual_and_jac(p, z)[1]

    # Damped Newton first (the system is square); trf only if Newton stalls.
    z = np.asarray(z0, dtype=np.float64).copy()
    for _it in range(min(25, max_nfev)):
        r, j, info = residual_and_jac(p, z)
        _log(
            f"{tag} newton |r| {np.linalg.norm(r):.3e} dr {max(info['dr']):.2e} "
            f"dv {max(info['dv']):.2e}"
        )
        if converged(info):
            return z, info
        d = np.linalg.norm(j, axis=0)
        d[d == 0.0] = 1.0
        step = np.linalg.lstsq(j / d, -r, rcond=None)[0] / d
        f0 = float(np.linalg.norm(r))
        for alpha in (1.0, 0.5, 0.25, 0.125, 0.0625):
            rn = residual_and_jac(p, z + alpha * step, False)[0]
            if np.all(np.isfinite(rn)) and float(np.linalg.norm(rn)) < f0:
                z = z + alpha * step
                break
        else:
            break
    z0 = z
    try:
        sol = least_squares(
            fun,
            z0,
            jac=jac,
            method="trf",
            x_scale="jac",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev,
        )
        zf = np.asarray(sol.x)
    except _DoneError:
        zf = np.asarray(cache["done"])
    _, _, info = residual_and_jac(p, zf, False)
    return zf, info


def scale_offsets(p: P, z: Arr, f: float) -> Arr:
    """Scale every node's moon-relative position offset by f (patched-conic r_p is linear in GM)."""
    z = z.copy()
    for k in range(4):
        rm, _ = p.circ.state(NODE_MOONS[k], float(z[7 * k + 6]))
        z[7 * k : 7 * k + 3] = rm + (z[7 * k : 7 * k + 3] - rm) * f
    return z


def stage_sigma(sig_list: list[float], max_nfev: int) -> None:
    """Stage 1: joint GM scale sigma on both moons (radii scaled by sigma), upward to 1."""
    p = build()
    path = OUT / "sigma.json"
    rec = json.loads(path.read_text()) if path.exists() else {"points": []}
    good = [q for q in rec["points"] if q["converged"]]
    hist = [(float(q["sigma"]), np.asarray(q["z"])) for q in good[-2:]]
    for sg in sig_list:
        if len(hist) == 2:
            (s0, z0), (s1, z1) = hist
            zs = z1 + (z1 - z0) * (math.log(sg) - math.log(s1)) / (math.log(s1) - math.log(s0))
        elif hist:
            zs = scale_offsets(p, hist[-1][1], sg / hist[-1][0])
        else:
            zs = scale_offsets(p, p.seed, sg)
        p.scale(sg, sg, sg * SATELLITES[GAN].radius_eq_km, sg * SATELLITES[CAL].radius_eq_km)
        zn, info = solve(p, zs, max_nfev, f"sigma={sg:.4g}")
        ok = converged(info)
        pt: dict[str, Any] = {
            "sigma": sg,
            "converged": ok,
            "z": zn.tolist(),
            "dr": max(info["dr"]),
            "dv": max(info["dv"]),
        }
        if ok:
            pt.update(describe(p, zn))
            hist = [*hist[-1:], (sg, zn)]
        rec["points"].append(pt)
        path.write_text(json.dumps(rec))
        _log(
            f"sigma={sg:.4g} converged={ok} "
            + (
                json.dumps(
                    [
                        (
                            n["moon"][0],
                            round(n["vinf_kms"], 4),
                            round(n["rp_km"] / sg, 1),
                            round(n["turn_deg"], 3),
                        )
                        for n in pt["nodes"]
                    ]
                )
                if ok
                else ""
            )
        )
        if not ok:
            break


def j_s(p: P, z: Arr, s_c: float, variant: str) -> Arr:
    """d residual / d s_C by central differences (no change of state)."""
    out = []
    for sgn in (1.0, -1.0):
        sv = s_c * (1.0 + sgn * 1e-5)
        set_sc(p, sv, variant)
        out.append(residual_and_jac(p, z, False)[0])
    set_sc(p, s_c, variant)
    return (out[0] - out[1]) / (2e-5 * s_c)


def set_sc(p: P, s_c: float, variant: str) -> None:
    """Continuation parameter: "S"/"P" = Callisto GM scale (radius scaled / physical);
    "J" = joint scale sigma on both moons with radii scaled (stage 1 by pseudo-arclength)."""
    if variant == "J":
        p.scale(s_c, s_c, s_c * SATELLITES[GAN].radius_eq_km, s_c * SATELLITES[CAL].radius_eq_km)
        return
    rad_c = s_c * SATELLITES[CAL].radius_eq_km if variant == "S" else SATELLITES[CAL].radius_eq_km
    p.scale(1.0, s_c, SATELLITES[GAN].radius_eq_km, rad_c)


def stage_cont(variant: str, n_points: int, ds_target: float, max_nfev: int) -> None:
    """Stage 2: pseudo-arclength continuation in (z, s_C), s_C decreasing from 1."""
    p = build()
    path = OUT / f"cont_{variant}.json"
    if path.exists():
        rec = json.loads(path.read_text())
    else:
        sig = [q for q in json.loads((OUT / "sigma.json").read_text())["points"] if q["converged"]]
        if variant == "J":
            s_start = float(sig[-1]["sigma"])
        else:
            assert sig and abs(sig[-1]["sigma"] - 1.0) < 1e-12, "stage sigma must reach 1 first"
            s_start = 1.0
        z0 = np.asarray(sig[-1]["z"])
        set_sc(p, s_start, variant)
        _, jz, _ = residual_and_jac(p, z0)
        dz = np.linalg.norm(jz, axis=0)
        dz[dz == 0.0] = 1.0
        ks = float(np.linalg.norm(j_s(p, z0, s_start, variant)))
        rec = {
            "variant": variant,
            "dz": dz.tolist(),
            "ks": ks,
            "dstarget": ds_target,
            "points": [{"s": s_start, "z": z0.tolist(), **describe(p, z0)}],
        }
    dz = np.asarray(rec["dz"])
    ks = float(rec["ks"])
    t_prev = np.asarray(rec["tangent"]) if "tangent" in rec else None
    for _ in range(n_points):
        last = rec["points"][-1]
        z, s_c = np.asarray(last["z"]), float(last["s"])
        if s_c <= 1e-3 or (variant == "J" and s_c >= 1.0):
            _log(f"cont {variant}: parameter limit reached ({s_c:.4g}), stop")
            break
        set_sc(p, s_c, variant)
        _, jz, _ = residual_and_jac(p, z)
        js = j_s(p, z, s_c, variant)
        jaug = np.hstack([jz / dz, (js / ks)[:, None]])
        _, sv, vt = np.linalg.svd(jaug)
        svz = np.linalg.svd(jz / dz, compute_uv=False)
        tan = vt[-1]
        if t_prev is None:  # stage 2 starts downward in s_C; stage 1 ("J") upward in sigma
            want_up = variant == "J"
            tan = tan if (tan[-1] > 0) == want_up else -tan
        elif float(tan @ t_prev) < 0.0:
            tan = -tan
        gap = float(sv[-2] / max(sv[-1], 1e-300))
        u0 = np.concatenate([z * dz, [s_c * ks]])
        dst = float(rec["dstarget"])
        ok = False
        for _h in range(6):
            ds = dst * ks / max(abs(tan[-1]), 1e-6)
            ds = min(ds, 50.0 * dst * ks)

            def fun(u: Arr, u0: Arr = u0, tan: Arr = tan, ds: float = ds) -> Arr:
                set_sc(p, float(u[-1]) / ks, variant)
                r = residual_and_jac(p, u[:-1] / dz, False)[0]
                return np.concatenate([r, [1e3 * (float(tan @ (u - u0)) - ds)]])

            def jac(u: Arr, tan: Arr = tan) -> Arr:
                sc = float(u[-1]) / ks
                set_sc(p, sc, variant)
                jz_ = residual_and_jac(p, u[:-1] / dz)[1]
                js_ = j_s(p, u[:-1] / dz, sc, variant)
                return np.vstack([np.hstack([jz_ / dz, (js_ / ks)[:, None]]), 1e3 * tan[None, :]])

            # Damped Newton on the square augmented system (residual rows + arclength row).
            un = u0 + ds * tan
            for _it in range(max_nfev):
                fu = fun(un)
                zt, st = un[:-1] / dz, float(un[-1]) / ks
                set_sc(p, st, variant)
                if converged(residual_and_jac(p, zt, False)[2]) and abs(fu[-1]) < 1e-6:
                    break
                ju = jac(un)
                step = np.linalg.lstsq(ju, -fu, rcond=None)[0]
                f0 = float(np.linalg.norm(fu))
                for alpha in (1.0, 0.5, 0.25, 0.125):
                    if float(np.linalg.norm(fun(un + alpha * step))) < f0:
                        un = un + alpha * step
                        break
                else:
                    break
            zn, sn = un[:-1] / dz, float(un[-1]) / ks
            set_sc(p, sn, variant)
            _, _, info = residual_and_jac(p, zn, False)
            if converged(info) and sn > 0.0:
                ok = True
                break
            dst *= 0.5
            _log(
                f"cont {variant}: step failed at s~{sn:.4g} (dr {max(info['dr']):.2e}); "
                f"target ds_C -> {dst:.3g}"
            )
        if not ok:
            _log(f"cont {variant}: 6 halvings without convergence at s_C {s_c:.5g}: NUMERICAL stop")
            break
        d = describe(p, zn)
        pt = {
            "s": sn,
            "z": zn.tolist(),
            **d,
            "sv_aug_min2": sv[-2:].tolist(),
            "sv_z_min2": svz[-2:].tolist(),
            "gap_aug": gap,
            "ds_dl_sign": float(np.sign(tan[-1])),
        }
        rec["points"].append(pt)
        rec["tangent"] = tan.tolist()
        rec["dstarget"] = min(dst * 1.3, 0.05)
        t_prev = tan
        path.write_text(json.dumps(rec))
        cn = [
            (round(n["rp_km"], 1), round(n["turn_deg"], 3), round(n["vinf_kms"], 4))
            for n in d["nodes"]
        ]
        _log(
            f"cont {variant} s_C={sn:.5g} gap {gap:.2e} svz_min {svz[-1]:.2e} "
            f"sign {pt['ds_dl_sign']:+.0f} nodes {cn}"
        )
        if variant == "P" and min(d["nodes"][1]["rp_km"], d["nodes"][2]["rp_km"]) <= R_CAL_RS_KM:
            _log("cont P: a Callisto periapsis reached the surface (2408 km); stop")
            break


def main() -> None:
    preflight_search(
        task_no=1004,
        region_id="gc2-callisto-gm-continuation-continuous",
        method=MethodCapability(
            genome="gc-2 (one structure), R-S circular coplanar, Ganymede + Callisto massive",
            corrector="forward-backward shooting, DOP853 + STM; pseudo-arclength in GM",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["sigma", "cont"], required=True)
    ap.add_argument("--list", default="")
    ap.add_argument("--variant", choices=["S", "P", "J"], default="S")
    ap.add_argument("--points", type=int, default=10)
    ap.add_argument("--ds", type=float, default=0.02)
    ap.add_argument("--max-nfev", type=int, default=60)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "sigma":
        stage_sigma([float(v) for v in args.list.split(",")], args.max_nfev)
    else:
        stage_cont(args.variant, args.points, args.ds, args.max_nfev)


if __name__ == "__main__":
    main()
