"""#496 Step 2 — feasibility-first joint corrector for the #411 cross-system cycle.

The bounded_ls solver (#496 first fix) was blocked because the seed-convergence gate
requires BOTH legs to converge at the seed (c_em, c_se) before the optimizer runs.
The new solver="feasibility_ls" relaxes that gate: it folds the leg POSITION GAPS into
a 4-component residual [fwd_km/scale, ret_km/scale, R1, R2] so trf drives leg closure
and phase closure JOINTLY from an infeasible starting point.

Seed: use the best (c_em, c_se) from the Step-1 scan (scan_411_leg_closing.py).  If
the scan finds a both-legs-converging seed, that seed is passed directly to bounded_ls
(the existing phase solver) for the cleanest path.  feasibility_ls is the fallback when
no such seed exists.

Run (background, runlog):
  uv run python scripts/close_411_feasibility_ls.py 2>&1 | tee runlogs/close_411_feasibility_ls.log
"""

from __future__ import annotations

import datetime as _dt

from cyclerfinder.genome.cross_system_cycle import (
    FrameBridge,
    correct_cross_cycle,
    crosscheck_cross_cycle,
    em_moon_system,
    se_earth_system,
)

# Seed: start from the Canalias bifurcation point (c_se=3.00086 is the last stable
# node below the bifurcation).  c_em=3.150 is the working EM-L2 value from Task 3.
# These are adjusted by the solver.  Use the scan output to update these if a
# both-legs-closing seed is found.
C_EM0 = 3.150
C_SE0 = 3.00086

# Bounds for the optimization: c_se BELOW Canalias; c_em inside EM-L2 family.
CANALIAS_C_SE = 3.000863625
C_EM_BOUNDS = (3.110, 3.152)
C_SE_BOUNDS = (3.0000, CANALIAS_C_SE - 1e-9)  # strict below bifurcation

# Denser return-leg scan — the default (4,2) is too coarse to find the convergence
# basin shown by test_se_to_em_return_leg_converges (which uses (8,3)).
RETURN_SCAN_N = 8
RETURN_SCAN_N_TAU = 3


def _ts() -> str:
    return _dt.datetime.now(tz=_dt.UTC).isoformat(timespec="seconds")


def _on_iter(k: int, c_em: float, c_se: float, r1: float, r2: float) -> None:
    import math

    norm = math.hypot(r1, r2)
    print(
        f"[{_ts()}]   nfev={k:03d}  c_em={c_em:.6f}  c_se={c_se:.9f}  "
        f"R=({r1:+.4f},{r2:+.4f})  |R|={norm:.4e}",
        flush=True,
    )


def _run_solver(bridge: FrameBridge, solver: str) -> None:
    print(f"[{_ts()}] === solver={solver} seed=(c_em={C_EM0}, c_se={C_SE0}) ===", flush=True)
    kw: dict[str, object] = {}
    if solver in {"bounded_ls", "feasibility_ls"}:
        kw = {
            "c_em_bounds": C_EM_BOUNDS,
            "c_se_bounds": C_SE_BOUNDS,
        }
    cyc = correct_cross_cycle(
        bridge,
        em_lib="EM-L2",
        se_lib="SE-L2",
        c_em0=C_EM0,
        c_se0=C_SE0,
        n_em=1,
        n_se=1,
        max_iter=16,
        solver=solver,
        on_iter=_on_iter,
        return_scan_n=RETURN_SCAN_N,
        return_scan_n_tau=RETURN_SCAN_N_TAU,
        **kw,  # type: ignore[arg-type]
    )
    print(
        f"[{_ts()}] {solver}: closed={cyc.closed} "
        f"theta_norm={cyc.theta_residual_norm:.4e} rad "
        f"max_leg={cyc.max_leg_residual_km:.2e} km "
        f"patch_dv={cyc.total_patch_dv_kms:.4f} km/s "
        f"c_em={cyc.c_em:.6f} c_se={cyc.c_se:.9f} "
        f"nfev/iter={cyc.n_iter} | {cyc.notes}",
        flush=True,
    )
    if cyc.closed:
        print(f"[{_ts()}] #411 CLOSED via {solver}!  Running Radau cross-check ...", flush=True)
        # Construct a CrossCycle wrapper the crosscheck function expects.
        from cyclerfinder.genome.cross_system_cycle import CrossCycle

        cc_wrap = CrossCycle(
            connections=[cyc.forward, cyc.ret],
            c_em=cyc.c_em,
            c_se=cyc.c_se,
            libration_pair=cyc.libration_pair,
            theta_closure_residual=cyc.theta_residual_norm,
            closed=True,
            max_leg_residual=cyc.max_leg_residual_km,
            independent_residual=float("nan"),
            notes="",
        )
        checked = crosscheck_cross_cycle(bridge, cc_wrap)
        print(
            f"[{_ts()}] Radau independent_residual = {checked.independent_residual:.3e} km",
            flush=True,
        )


def main() -> None:
    print(f"[{_ts()}] #496 feasibility-first joint corrector for #411 cross-system cycle")
    bridge = FrameBridge(se=se_earth_system(), em=em_moon_system())
    # Run feasibility_ls (the new solver that relaxes the seed gate).
    _run_solver(bridge, "feasibility_ls")
    # Also run bounded_ls with denser return-leg scan (denser than previous attempt)
    # for comparison: if both legs converge with return_scan_n=8, bounded_ls should work.
    _run_solver(bridge, "bounded_ls")
    print(f"[{_ts()}] done", flush=True)


if __name__ == "__main__":
    main()
