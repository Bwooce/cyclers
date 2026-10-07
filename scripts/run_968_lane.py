"""#968 corrector B: the Jovian lane's own residual (fixed wrap) on GanCal#5.

Pre-registration: docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 3.3. The lane's
``jovian_defect_residual`` (REBOUND IAS15 on a rails spline) with ``jovian_stm_jacobian``, driven
by ``least_squares`` trf exactly as ``jovian_shoot`` does, but on the R-S circular ephemeris (the
lane's ``jovian_shoot`` hard-wires jup365). Nodes: the patched-conic periapsis nodes B, C, A, B' at
the patched-conic epochs (fixed); moons = ("Ganymede",) (Callisto massless, as in R-S).

    uv run python scripts/run_968_lane.py --max-nfev 40   # resumable, data/968_control/b_state.json
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

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
OUT = ROOT / "data" / "968_control"
DAY = 86400.0
MOONS = ("Ganymede",)


def _control() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_968_control", ROOT / "scripts" / "run_968_control.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_968_control"] = mod
    spec.loader.exec_module(mod)
    return mod


def build_seed(ctl: Any) -> tuple[Any, ShootingSeed]:
    p = ctl.build_problem()
    enum = ctl._load_enum()
    circ, a, b = enum.cell_system("gc1")
    cands = {
        c["key"]: c
        for c in json.loads((ROOT / "data/943_cell_gc_gauntlet.json").read_text())["candidates"]
    }
    _, cyc = enum.parse_cycle_key(ctl.KEY, circ, a, b)
    fl = cycle_flybys(circ, cyc, np.asarray(cands[ctl.KEY]["x_days"]) * DAY)
    fb, fc, fa = fl
    rb, vb, _ = periapsis_node("Ganymede", fb.t_s, fb.vinf_in, fb.vinf_out, p.circ)
    ra, va, _ = periapsis_node("Ganymede", fa.t_s, fa.vinf_in, fa.vinf_out, p.circ)
    rc, vc = p.circ.state("Callisto", fc.t_s)
    xb = np.concatenate([rb, vb])
    states = [xb, np.concatenate([rc, vc + fc.vinf_in]), np.concatenate([ra, va]), p.q6 @ xb]
    epochs = [fb.t_s, fc.t_s, fa.t_s, fb.t_s + p.period_s]
    rot = p.q6[:3, :3]
    seed = ShootingSeed(
        node_states=states,
        epochs=epochs,
        tofs=[(epochs[i + 1] - epochs[i]) / DAY for i in range(3)],
        sequence=("Ganymede", "Callisto", "Ganymede", "Ganymede"),
        slack_leg=0,
        period_days=p.period_s / DAY,
        vinf_in=[rot @ fb.vinf_in, fc.vinf_in, fa.vinf_in, rot @ fb.vinf_in],
        vinf_out=[fb.vinf_out, fc.vinf_out, fa.vinf_out, rot @ fb.vinf_out],
    )
    return p, seed


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="gancal5-continuous-ideal-positive-control-lane",
        method=MethodCapability(
            genome="R-S 2009 GanCal#5 (one structure), circular coplanar, Callisto massless",
            corrector="jovian_defect_residual (fixed wrap) + jovian_stm_jacobian, trf",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-nfev", type=int, default=40)
    ap.add_argument("--seed", choices=["pc", "a"], default="pc")
    args = ap.parse_args()
    ctl = _control()
    p, seed = build_seed(ctl)
    if args.seed == "a":  # amendment 4: corrector A's converged pinned orbit, its epochs
        za = np.asarray(json.loads((OUT / "a_state.json").read_text())["z"])
        nd = ctl.nodes(p, za)
        seed = ShootingSeed(
            node_states=[np.asarray(q[0]) for q in nd],
            epochs=[float(q[2]) for q in nd],
            tofs=[(nd[i + 1][2] - nd[i][2]) / DAY for i in range(3)],
            sequence=seed.sequence,
            slack_leg=0,
            period_days=seed.period_days,
            vinf_in=seed.vinf_in,
            vinf_out=seed.vinf_out,
        )
    n = len(seed.sequence)
    cache = JovianRailsCache(MOONS, p.circ, min(seed.epochs), max(seed.epochs))  # type: ignore[arg-type]
    state_path = OUT / f"b_state_{args.seed}seed.json"
    x = (
        np.asarray(json.loads(state_path.read_text())["x"])
        if state_path.exists()
        else _states_to_x(seed.node_states)
    )

    def res(xv: np.ndarray) -> np.ndarray:
        r = jovian_defect_residual(
            _seed_with_states(seed, _x_to_states(xv, n)),
            ephem=p.circ,  # type: ignore[arg-type]
            cache=cache,
            moons=MOONS,
            max_wall_sec=3000.0,
        )
        if np.any(r == 1e9):
            ctl._log("lane: a leg reported timeout/non-finite (1e9 sentinel)")
        legs = r[: 6 * (n - 1)].reshape(n - 1, 6)
        dr_max = float(np.max(np.linalg.norm(legs[:, :3], axis=1)))
        ctl._log(
            f"lane |r| {np.linalg.norm(r):.4e} leg dr max {dr_max:.3e} km "
            f"dv max {np.max(np.linalg.norm(legs[:, 3:], axis=1)) / _W_VEL:.3e} km/s "
            f"wrap dr {np.linalg.norm(r[-6:-3]):.3e} dv {np.linalg.norm(r[-3:]) / _W_VEL:.3e}"
        )
        return r

    def jac(xv: np.ndarray) -> np.ndarray:
        return jovian_stm_jacobian(seed, xv, ephem=p.circ, moons=MOONS)  # type: ignore[arg-type]

    sol = least_squares(res, x, jac=jac, method="trf", x_scale="jac", max_nfev=args.max_nfev)
    xf = np.asarray(sol.x)
    rf = res(xf)
    legs = rf[: 6 * (n - 1)].reshape(n - 1, 6)
    out: dict[str, Any] = {
        "x": xf.tolist(),
        "norm": float(np.linalg.norm(rf)),
        "leg_dr_km": np.linalg.norm(legs[:, :3], axis=1).tolist(),
        "leg_dv_kms": (np.linalg.norm(legs[:, 3:], axis=1) / _W_VEL).tolist(),
        "wrap_dr_km": float(np.linalg.norm(rf[-6:-3])),
        "wrap_dv_kms": float(np.linalg.norm(rf[-3:]) / _W_VEL),
    }
    states = _x_to_states(xf, n)
    gan = []
    for i in (0, 2):
        v, d = ctl.vinf_gan(p, np.asarray(states[i]), seed.epochs[i])
        gan.append({"vinf_kms": v, "alt_rs_km": d - ctl.R_GAN_RS_KM})
    out["ganymede_nodes"] = gan
    rc, _ = p.circ.state("Callisto", seed.epochs[1])
    out["callisto_node_dist_km"] = float(np.linalg.norm(np.asarray(states[1])[:3] - rc))
    state_path.write_text(json.dumps(out))
    ctl._log(f"lane end: {json.dumps({k: v for k, v in out.items() if k != 'x'})}")


if __name__ == "__main__":
    main()
