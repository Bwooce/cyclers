"""#496 — break the #411 cross-system cycle stall with the bounded-least_squares solver.

The #411 time-consistent corrector (`correct_cross_cycle`, solver="newton") STALLS at
|R| ~ 0.59 rad: its damped FD-Newton step pushes c_se off the razor-thin SE-L2 Lyapunov
family (the family lives at C_se in ~[3.0000, 3.0008] with the Canalias bifurcation right
at 3.000863625), so no improving step keeps both legs converged and it breaks.

#496 fix (Braik-Ross 2025 / Ross-Scheeres blueprint: bound + feasibility-first): the new
solver="bounded_ls" path runs a bounded trust-region least_squares on (R1, R2) with c_se
CLAMPED below the bifurcation, so the solver cannot leave the stable family; the EM-L2
amplitude (the covering knob, per the #411 amplitude analysis) closes theta within bounds.

This driver runs BOTH solvers from the same seed and prints theta_residual_norm +
total_patch_dv so the stall-break is evidenced (or honestly not). Compute-heavy (the SE
return leg integrates many ~58-day SE periods); run detached + read the runlog.

Run: uv run python scripts/close_411_bounded_ls.py
"""

from __future__ import annotations

import datetime as _dt

from cyclerfinder.genome.cross_system_cycle import (
    FrameBridge,
    correct_cross_cycle,
    em_moon_system,
    se_earth_system,
)

CANALIAS_C_SE = 3.000863625  # SE-L2 saddle-center bifurcation; stay strictly BELOW it
C_EM0 = 3.150  # EM-L2 Lyapunov seed (family ~[3.108, 3.153]); the covering knob
C_SE0 = 3.00086  # SE-L2 seed at the Canalias node (the only known leg-converging region)
C_EM_BOUNDS = (3.112, 3.152)  # keep c_em inside the EM-L2 family
C_SE_BOUNDS = (3.00050, 3.00086)  # node-converging band up to the Canalias bifurcation


def _ts() -> str:
    return _dt.datetime.now(tz=_dt.UTC).isoformat(timespec="seconds")


def _on_iter(k: int, c_em: float, c_se: float, r1: float, r2: float) -> None:
    print(f"[{_ts()}]   it={k:02d}  c_em={c_em:.6f}  c_se={c_se:.9f}  R=({r1:+.4f},{r2:+.4f})")


def _run(bridge: FrameBridge, solver: str) -> None:
    print(f"[{_ts()}] === solver={solver} ===", flush=True)
    kw: dict[str, object] = {}
    if solver == "bounded_ls":
        kw = {"c_em_bounds": C_EM_BOUNDS, "c_se_bounds": C_SE_BOUNDS}
    cyc = correct_cross_cycle(
        bridge,
        em_lib="EM-L2",
        se_lib="SE-L2",
        c_em0=C_EM0,
        c_se0=C_SE0,
        n_em=1,
        n_se=1,
        max_iter=12,
        solver=solver,
        on_iter=_on_iter,
        **kw,  # type: ignore[arg-type]
    )
    print(
        f"[{_ts()}] {solver}: closed={cyc.closed} "
        f"theta_residual_norm={cyc.theta_residual_norm:.4e} rad "
        f"(seed stall ~0.59) | max_leg={cyc.max_leg_residual_km:.2e} km | "
        f"patch_dv={cyc.total_patch_dv_kms:.4f} km/s | c_em={cyc.c_em:.6f} "
        f"c_se={cyc.c_se:.9f} | nfev/iter={cyc.n_iter} | {cyc.notes}",
        flush=True,
    )


def main() -> None:
    print(f"[{_ts()}] #496 bounded-ls vs newton on the #411 cross-system cycle", flush=True)
    bridge = FrameBridge(se=se_earth_system(), em=em_moon_system())
    _run(bridge, "newton")  # baseline: reproduce the ~0.59 rad stall
    _run(bridge, "bounded_ls")  # the fix: c_se clamped below the Canalias bifurcation
    print(f"[{_ts()}] done", flush=True)


if __name__ == "__main__":
    main()
