"""#1044: GanCal#1 at 2013 and gc-1 (e2) on jup365 under #1039 formulation (b') (pre-registration:
docs/notes/2026-10-08-1044-realeph-controls-gancal1-gc1.md). A wrapper of ``run_968_rungb.py``
with four encounters per cycle and Ganymede + Callisto massive.

    uv run python scripts/run_1044_run.py --target gancal1|gc1 --stage sigma|down|ias15
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from typing import Any

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


def _install_lm(rb: Any) -> None:
    """Amendment 2: Levenberg-Marquardt (scipy, analytic Jacobian) before the damped-Newton polish.
    The line-search Newton crawls on the full-revolution resonant legs (strongly nonlinear in the
    node epochs). LM state is checkpointed so that a call can stop and the next resume."""
    import numpy as np
    from scipy.optimize import least_squares

    newton0 = rb.newton
    ckpt = rb.OUT / "lm_checkpoint.npy"

    def newton_lm(c: Any, z: Any, iters: int, tag: str) -> tuple[Any, dict[str, Any]]:
        if ckpt.exists():
            zc = np.load(ckpt)
            if zc.shape == z.shape and float(np.max(np.abs(zc[6::7] - z[6::7]))) < 3600.0:
                z = zc
        cache: dict[str, Any] = {}

        def fun(zz: Any) -> Any:
            try:
                r, j, info = rb.residual_and_jac(c, zz)
            except RuntimeError:
                return np.full(7 * c.m, 1e6)
            cache["z"], cache["j"] = zz.copy(), j
            np.save(ckpt, zz)
            rb._log(
                f"{tag} lm |r| {np.linalg.norm(r):.3e} dr {max(info['dr']):.2e} "
                f"dv {max(info['dv']):.2e}"
            )
            return r

        def jac(zz: Any) -> Any:
            if "z" in cache and np.array_equal(cache["z"], zz):
                return cache["j"]
            return rb.residual_and_jac(c, zz)[1]

        sol = least_squares(
            fun,
            z,
            jac=jac,
            method="lm",
            x_scale="jac",
            max_nfev=25,
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
        )
        zf, info = newton0(c, np.asarray(sol.x), 6, tag + " polish")
        if rb.converged(info) or (max(info["dr"]) < 1e-2 and max(info["dv"]) < 1e-6):
            ckpt.unlink(missing_ok=True)
        return zf, info

    rb.newton = newton_lm


def main() -> None:
    preflight_search(
        task_no=1044,
        region_id="gancal1-gc1-jup365-continuous",
        method=MethodCapability(
            genome="one-cycle open chains on jup365 (GanCal#1 2013; gc-1 e2)",
            corrector="forward-backward shooting, DOP853 + STM; sigma continuation; (b') ends",
            capability_tags=frozenset({"ballistic", "n-body", "real-ephemeris"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["gancal1", "gc1"], required=True)
    ap.add_argument("--stage", choices=["sigma", "down", "ias15"], required=True)
    ap.add_argument("--solver", choices=["newton", "lm"], default="lm")
    args = ap.parse_args()
    rb = _load("run_968_rungb", ROOT / "scripts" / "run_968_rungb.py")
    rb.OUT = ROOT / "data" / f"1044_{args.target}"
    rb.END_MODE = "direction_seedmag"
    rb.PER_CYCLE = 4
    rb.FORCE = ("Ganymede", "Callisto")
    rb.SCALED = ("Ganymede", "Callisto")
    rb.SEED_CAP_KM = 1.0e12  # amendment 1: small turns need r_p beyond 0.6 SOI (as #1043)
    if args.solver == "lm":
        _install_lm(rb)
    {"sigma": rb.stage_sigma, "down": rb.stage_down, "ias15": rb.stage_ias15}[args.stage](1)


if __name__ == "__main__":
    main()
