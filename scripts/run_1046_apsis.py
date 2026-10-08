"""#1046: full-revolution-aware open chain: apojove nodes (pre-registration:
docs/notes/2026-10-08-1046-apsis-node-open-chain.md).

The #1039 / #1044 open chain (``run_968_rungb.py``) with one extra multiple-shooting node at every
apojove of the seed trajectory between two flyby nodes (at least 1 day from both). An apsis node's
unknowns are its absolute Jupiter-centred state and epoch; its row is the Jupiter-centred apsis
gauge r.v / (|r||v|) = 0. It is carried as a pseudo-moon "APO" whose ephemeris state is zero, so
the rungb node, gauge and leg code is reused unchanged (moon-relative = absolute for APO).

    uv run python scripts/run_1046_apsis.py --target eggie|gancal1|gc1
        --stage build|sigma|down|ias15

Targets: the one-cycle EGGIE ideal control (#1039 form (b), damped Newton), then GanCal#1 at 2013
and gc-1 e2 on jup365 (#1044 form (b'), LM then Newton). Checkpoints in data/1046_<target>/; live
log in .../live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2
from cyclerfinder.search.two_working_body import kepler_step

ROOT = Path(__file__).resolve().parents[1]
APO = "APO"
DAY = 86400.0
MIN_GAP_S = 1.0 * DAY


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class ApsisEphem:
    """Forwards every call to the lane ephemeris; the pseudo-moon APO sits at Jupiter's centre."""

    def __init__(self, inner: Any) -> None:
        self._inner = inner

    def __getattr__(self, name: str) -> Any:
        return getattr(self._inner, name)

    def state(self, moon: str, t_sec: float) -> tuple[Any, Any]:
        if moon == APO:
            return np.zeros(3), np.zeros(3)
        return self._inner.state(moon, t_sec)  # type: ignore[no-any-return]

    def position(self, moon: str, t_sec: float) -> Any:
        if moon == APO:
            return np.zeros(3)
        return self._inner.position(moon, t_sec)


def apojoves(r: Any, v: Any, t0: float, t1: float) -> list[tuple[float, Any]]:
    """Apojove passages of the Jupiter-centred conic through (r, v) at t0, in
    (t0 + 1 d, t1 - 1 d)."""
    mu = MU_JUPITER_KM3_S2
    rn = float(np.linalg.norm(r))
    a = 1.0 / (2.0 / rn - float(v @ v) / mu)
    if a <= 0.0:
        return []
    ev = ((float(v @ v) - mu / rn) * r - float(r @ v) * v) / mu
    e = float(np.linalg.norm(ev))
    nu = math.atan2(
        float(np.cross(ev, r) @ np.cross(r, v)) / (e * float(np.linalg.norm(np.cross(r, v)))),
        float(ev @ r) / e,
    )
    ea = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(nu / 2.0), math.sqrt(1.0 + e) * math.cos(nu / 2.0)
    )
    m0 = ea - e * math.sin(ea)
    n = math.sqrt(mu / a**3)
    dt = ((math.pi - m0) % (2.0 * math.pi)) / n
    out = []
    while t0 + dt < t1 - MIN_GAP_S:
        if dt > MIN_GAP_S:
            rr, vv = kepler_step(r, v, dt, mu)
            out.append((t0 + dt, np.concatenate([rr, vv])))
        dt += 2.0 * math.pi / n
    return out


def install(rb: Any) -> None:
    """Patch the rungb module: build() gains apsis nodes; scale_offsets and describe skip them."""
    build0, scale0, describe0 = rb.build, rb.scale_offsets, rb.describe

    def build(n_cycles: int) -> Any:
        c = build0(n_cycles)
        c.eph = ApsisEphem(c.eph)
        moons, z = [], []
        inserted = []
        for k in range(c.m):
            moon = c.moons[k]
            zk = c.seed[7 * k : 7 * k + 7]
            moons.append(moon)
            z.extend(zk.tolist())
            if k == c.m - 1:
                break
            t = float(zk[6])
            t1 = float(c.seed[7 * (k + 1) + 6])
            rm, vm = c.eph.state(moon, t)
            vout = rb.asymptote(SATELLITES[moon].mu_km3_s2, zk[0:3], zk[3:6], "out")
            for ta, xa in apojoves(np.asarray(rm), np.asarray(vm) + vout, t, t1):
                moons.append(APO)
                z.extend([*xa, ta])
                inserted.append(
                    {
                        "after_node": k,
                        "t_days_after": (ta - t) / DAY,
                        "r_km": float(np.linalg.norm(xa[:3])),
                    }
                )
        c.moons = moons
        c.seed = np.asarray(z, dtype=np.float64)
        c.apsis_inserted = inserted
        return c

    def scale_offsets(c: Any, z: Any, f: float) -> Any:
        y = scale0(c, z, f)
        for k, moon in enumerate(c.moons):
            if moon == APO:
                y[7 * k : 7 * k + 3] = z[7 * k : 7 * k + 3]
        return y

    def describe(c: Any, z: Any, sigma: float) -> dict[str, Any]:
        keep = [k for k, moon in enumerate(c.moons) if moon != APO]
        apo = [k for k, moon in enumerate(c.moons) if moon == APO]
        moons0 = c.moons
        zz = np.concatenate([z[7 * k : 7 * k + 7] for k in keep])
        c.moons = [moons0[k] for k in keep]
        try:
            out: dict[str, Any] = describe0(c, zz, sigma)
        finally:
            c.moons = moons0
        out["apsis"] = [
            {"r_km": float(np.linalg.norm(z[7 * k : 7 * k + 3])), "t_s": float(z[7 * k + 6])}
            for k in apo
        ]
        return out

    rb.build, rb.scale_offsets, rb.describe = build, scale_offsets, describe


def main() -> None:
    preflight_search(
        task_no=1046,
        region_id="apsis-node-open-chain-controls",
        method=MethodCapability(
            genome="one-cycle open chains (EGGIE ideal; GanCal#1 2013, gc-1 e2 on jup365)",
            corrector="forward-backward shooting + apojove nodes; sigma continuation; (b), (b')",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["eggie", "gancal1", "gc1"], required=True)
    ap.add_argument("--stage", choices=["build", "sigma", "down", "ias15"], required=True)
    args = ap.parse_args()
    out = ROOT / "data" / f"1046_{args.target}"
    out.mkdir(parents=True, exist_ok=True)
    if args.target == "eggie":
        w = _load("run_1043_eggie", ROOT / "scripts" / "run_1043_eggie.py")
        rb = w.Q.R
        w.OUT = w.Q.OUT = rb.OUT = out
        rb.END_MODE = "direction"
        seed = out / "pc_chain_n1.json"
        if not seed.exists():
            shutil.copy(ROOT / "data" / "1043_eggie" / "pc_chain_n1.json", seed)
        ias15 = w.Q.stage_ias15
        install(rb)
        w.Q.build = rb.build  # Q.stage_ias15 calls its module-level build()
    else:
        w = _load("run_1044_run", ROOT / "scripts" / "run_1044_run.py")
        rb = _load("run_968_rungb", ROOT / "scripts" / "run_968_rungb.py")
        rb.OUT = out
        rb.END_MODE = "direction_seedmag"
        rb.PER_CYCLE = 4
        rb.FORCE = ("Ganymede", "Callisto")
        rb.SCALED = ("Ganymede", "Callisto")
        rb.SEED_CAP_KM = 1.0e12
        seed = out / "seed_chain.json"
        if not seed.exists():
            shutil.copy(ROOT / "data" / f"1044_{args.target}" / "seed_chain.json", seed)
        w._install_lm(rb)
        ias15 = rb.stage_ias15
        install(rb)
    if args.stage == "build":
        c = rb.build(1)
        info = {"moons": c.moons, "apsis_inserted": c.apsis_inserted}
        (out / "build.json").write_text(json.dumps(info, indent=1))
        rb._log(f"build: nodes {c.moons}; apsis {c.apsis_inserted}")
        for sg in (0.02,):
            c.set_sigma(sg)
            r, _, inf = rb.residual_and_jac(c, rb.scale_offsets(c, c.seed, sg), False)
            rb._log(
                f"build: seed at sigma {sg}: |r| {np.linalg.norm(r):.3e} dr {max(inf['dr']):.3e}"
            )
        return
    {"sigma": rb.stage_sigma, "down": rb.stage_down, "ias15": ias15}[args.stage](1)


if __name__ == "__main__":
    main()
