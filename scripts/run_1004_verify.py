"""#1004 amendment 1 (note sec. 8.2): checks of the fold in the joint mass scale sigma.

    uv run python scripts/run_1004_verify.py --check two      # 1: two solutions at one sigma
    uv run python scripts/run_1004_verify.py --check bracket  # 2: fold bracket on the lower branch
    uv run python scripts/run_1004_verify.py --check ias15    # 3: IAS15 re-fly of the two solutions
    uv run python scripts/run_1004_verify.py --check sigma1   # 4: direct attempts at sigma = 1
    uv run python scripts/run_1004_verify.py --check noise    # 5: J_z small singular values vs rtol

Results append to data/1004_gc2/verify.json.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

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
VER = G.OUT / "verify.json"


def _save(key: str, val: Any) -> None:
    rec = json.loads(VER.read_text()) if VER.exists() else {}
    rec[key] = val
    VER.write_text(json.dumps(rec, indent=1))


def branches() -> tuple[list[tuple[float, Any]], list[tuple[float, Any]]]:
    lower = [
        (q["sigma"], np.asarray(q["z"]))
        for q in json.loads((G.OUT / "sigma.json").read_text())["points"]
        if q["converged"]
    ]
    cont = json.loads((G.OUT / "cont_J.json").read_text())["points"]
    upper = [(q["s"], np.asarray(q["z"])) for q in cont[2:]]  # after the turning point
    return lower, upper


def solve_at(p: Any, sigma: float, z0: Any, tag: str) -> tuple[Any, dict[str, Any], bool]:
    G.set_sc(p, sigma, "J")
    z, info = G.solve(p, z0, 40, tag)
    return z, info, bool(G.converged(info))


def summary(p: Any, z: Any, sigma: float) -> dict[str, Any]:
    d = G.describe(p, z)
    return {
        "vinf": [n["vinf_kms"] for n in d["nodes"]],
        "rp_over_sigma": [n["rp_km"] / sigma for n in d["nodes"]],
        "turn_deg": [n["turn_deg"] for n in d["nodes"]],
        "c1c2_a_km": d["c1c2_mid_a_km"],
        "c1c2_e": d["c1c2_mid_e"],
    }


def check_two() -> None:
    p = G.build()
    lower, upper = branches()
    out = {}
    for sg in (0.150, 0.160):
        zl = min(lower, key=lambda q: abs(q[0] - sg))[1]
        zu = min(upper, key=lambda q: abs(q[0] - sg))[1]
        a, _ia, oka = solve_at(p, sg, zl, f"two lower {sg}")
        b, _ib, okb = solve_at(p, sg, zu, f"two upper {sg}")
        row = {"lower_conv": oka, "upper_conv": okb}
        if oka:
            row["lower"] = summary(p, a, sg)
            row["lower_z"] = a.tolist()
        if okb:
            row["upper"] = summary(p, b, sg)
            row["upper_z"] = b.tolist()
        if oka and okb:
            row["d_vinf_C"] = abs(row["lower"]["vinf"][1] - row["upper"]["vinf"][1])
            row["d_rp_over_sigma_C"] = abs(
                row["lower"]["rp_over_sigma"][1] - row["upper"]["rp_over_sigma"][1]
            )
            row["fold_supported"] = row["d_vinf_C"] > 0.02 and row["d_rp_over_sigma_C"] > 100.0
        out[str(sg)] = row
        G._log(
            f"verify two sigma={sg}: "
            f"{json.dumps({k: v for k, v in row.items() if not k.endswith('_z')})}"
        )
    _save("two", out)


def check_bracket() -> None:
    p = G.build()
    lower, _ = branches()
    (s0, z0), (s1, z1) = lower[-2], lower[-1]
    last_ok = (s1, z1)
    prev = (s0, z0)
    step = 0.0005
    tried = []
    while step >= 1e-5:
        sg = last_ok[0] + step
        zp = last_ok[1] + (last_ok[1] - prev[1]) * (sg - last_ok[0]) / (last_ok[0] - prev[0])
        z, info, ok = solve_at(p, sg, zp, f"bracket {sg:.6f}")
        tried.append({"sigma": sg, "converged": ok, "dr": max(info["dr"]), "dv": max(info["dv"])})
        _save("bracket", {"tried": tried, "last_converged": last_ok[0]})
        if ok:
            prev, last_ok = last_ok, (sg, z)
        else:
            step = step / 5.0
    G._log(f"verify bracket: last converged sigma {last_ok[0]:.6f}")
    _save("bracket", {"tried": tried, "last_converged": last_ok[0], "last_z": last_ok[1].tolist()})


def check_ias15() -> None:
    from cyclerfinder.nbody.jovian import JovianRestrictedNBody

    p = G.build()
    rec = json.loads(VER.read_text())["two"]["0.15"]
    out = {}

    class Exact:
        def position(self, moon: str, t: float) -> Any:
            return p.circ.state(moon, t)[0]

    for lab in ("lower", "upper"):
        z = np.asarray(rec[f"{lab}_z"])
        G.set_sc(p, 0.150, "J")
        nd = G.nodes(p, z)
        rows = []
        for k in range(4):
            xp, _, tp, _ = nd[k]
            xn, _, tn, _ = nd[k + 1]
            tm = 0.5 * (tp + tn)
            for d_, x0, t0 in (("fwd", xp, tp), ("bwd", xn, tn)):
                xd, _ = G.arc(p, x0, t0, tm)
                # IAS15 with the same scaled GMs (the lane propagator reads registry GMs).
                arc = _ias15_scaled(p, x0, t0, tm, Exact())
                rows.append(
                    {
                        "leg": k,
                        "dir": d_,
                        "dr_km": float(np.linalg.norm(arc[:3] - xd[:3])),
                        "dv_kms": float(np.linalg.norm(arc[3:] - xd[3:])),
                    }
                )
        out[lab] = rows
        G._log(
            f"verify ias15 {lab}: max dr {max(r['dr_km'] for r in rows):.2e} "
            f"dv {max(r['dv_kms'] for r in rows):.2e}"
        )
    _ = JovianRestrictedNBody
    _save("ias15", out)


def _ias15_scaled(p: Any, x0: Any, t0: float, t1: float, cache: Any) -> Any:
    """REBOUND IAS15 with the scaled moon GMs and softening radii of ``p`` (second integrator)."""
    import rebound

    from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2

    sim = rebound.Simulation()
    sim.G = MU_JUPITER_KM3_S2
    sim.integrator = "ias15"
    sim.integrator.epsilon = 1e-11
    sim.add(m=1.0)
    sim.add(m=0.0, x=x0[0], y=x0[1], z=x0[2], vx=x0[3], vy=x0[4], vz=x0[5])
    sim.t = float(t0)
    mus, surf = dict(p.mus), dict(p.surf)

    def forces(ptr: Any) -> None:
        s = ptr.contents
        sc = s.particles[1]
        r = np.array([sc.x, sc.y, sc.z])
        acc = np.zeros(3)
        for m in G.MOONS:
            rm = cache.position(m, float(s.t))
            d = rm - r
            dn = max(float(np.linalg.norm(d)), surf[m])
            acc += mus[m] * (d / dn**3 - rm / float(np.linalg.norm(rm)) ** 3)
        sc.ax += acc[0]
        sc.ay += acc[1]
        sc.az += acc[2]

    sim.additional_forces = forces
    sim.force_is_velocity_dependent = 0
    sim.integrate(float(t1))
    sc = sim.particles[1]
    return np.array([sc.x, sc.y, sc.z, sc.vx, sc.vy, sc.vz])


def check_sigma1() -> None:
    p = G.build()
    lower, upper = branches()
    out = {}
    starts = {
        "patched_conic_seed": p.seed,
        "lower_last": G.scale_offsets(p, lower[-1][1], 1.0 / lower[-1][0]),
        "upper_last": G.scale_offsets(p, upper[-1][1], 1.0 / upper[-1][0]),
    }
    for lab, z0 in starts.items():
        z, info, ok = solve_at(p, 1.0, z0, f"sigma1 {lab}")
        row: dict[str, Any] = {"converged": ok, "dr": max(info["dr"]), "dv": max(info["dv"])}
        if ok:
            row.update(summary(p, z, 1.0))
            row["z"] = z.tolist()
        out[lab] = row
        G._log(f"verify sigma1 {lab}: {json.dumps({k: v for k, v in row.items() if k != 'z'})}")
    _save("sigma1", out)


def check_noise() -> None:
    p = G.build()
    lower, _ = branches()
    sg, z = lower[-1]
    G.set_sc(p, sg, "J")
    out = {}
    for rtol in (1e-12, 1e-13):
        G.RTOL = rtol
        G.ATOL = 1e-10 * rtol / 1e-12
        _, jz, _ = G.residual_and_jac(p, z)
        d = np.linalg.norm(jz, axis=0)
        out[str(rtol)] = np.linalg.svd(jz / d, compute_uv=False)[-6:].tolist()
        G._log(f"verify noise rtol {rtol}: smallest sv {out[str(rtol)]}")
    _save("noise", out)


def main() -> None:
    preflight_search(
        task_no=1004,
        region_id="gc2-sigma-fold-verification",
        method=MethodCapability(
            genome="gc-2 (one structure), R-S circular coplanar, Ganymede + Callisto massive",
            corrector="forward-backward shooting, DOP853 + STM; Newton at fixed sigma",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--check", choices=["two", "bracket", "ias15", "sigma1", "noise"], required=True
    )
    args = ap.parse_args()
    {
        "two": check_two,
        "bracket": check_bracket,
        "ias15": check_ias15,
        "sigma1": check_sigma1,
        "noise": check_noise,
    }[args.check]()


if __name__ == "__main__":
    main()
