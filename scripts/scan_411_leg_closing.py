"""#496 Step 1 — leg-closing seed scan over (c_em, c_se) grid.

Tests both the forward EM->SE leg AND the return SE->EM leg at each (c_em, c_se)
grid point using a DENSER scan than the cycle-corrector default (return_scan_n=4),
which the 2026-06-30 verdict identified as a likely contributor to the "legs don't
converge" wall.  The test `test_se_to_em_return_leg_converges` shows the return leg
DOES converge with scan_n=8, scan_n_tau=3 and branch_u=-1, branch_s=+1 at
c_em=3.15, c_se=3.00086; this scan checks whether that convergence extends across
the full (c_em, c_se) family range.

Run (background, runlog):
  uv run python scripts/scan_411_leg_closing.py 2>&1 | tee runlogs/scan_411_leg_closing.log
"""

from __future__ import annotations

import datetime as _dt
import sys

from cyclerfinder.genome.cross_system_cycle import (
    _RETURN_VARIANTS,
    FrameBridge,
    _build_em_node,
    _build_se_node,
    correct_cross_connection,
    em_moon_system,
    se_earth_system,
)

# Grid: c_em covers the EM-L2 Lyapunov family (~3.108-3.153).
# c_se covers only the KNOWN node-converging values of the SE-L2 family
# (patchy: nodes build at these but NOT at 3.0004/3.0005; the family ends at
# the Canalias bifurcation 3.000863625).
C_EM_VALUES = [3.110, 3.120, 3.130, 3.140, 3.148, 3.150, 3.152]
C_SE_VALUES = [3.0000, 3.0002, 3.0006, 3.00086]

# Denser scan than the corrector default (return_scan_n=4, return_scan_n_tau=2).
# test_se_to_em_return_leg_converges uses scan_n=8, scan_n_tau=3 and passes.
SCAN_N_FWD = 12  # forward leg: theta grid size
SCAN_N_TAU_FWD = 4  # forward leg: tau_u/tau_s grid size
SCAN_N_RET = 8  # return leg: theta grid size (denser than default 4)
SCAN_N_TAU_RET = 3  # return leg: tau grid size (denser than default 2)
MAX_TIME_FACTOR = 6.0  # manifold integration horizon (x orbit period)
TOL_KM = 1e2  # position-gap convergence threshold


def _ts() -> str:
    return _dt.datetime.now(tz=_dt.UTC).isoformat(timespec="seconds")


def main() -> None:
    n_em, n_se = len(C_EM_VALUES), len(C_SE_VALUES)
    print(f"[{_ts()}] #496 leg-closing seed scan: c_em x{n_em} c_se x{n_se}")
    n_fwd = f"({SCAN_N_FWD},{SCAN_N_TAU_FWD})"
    n_ret = f"({SCAN_N_RET},{SCAN_N_TAU_RET})"
    print(f"[{_ts()}] dense scan: fwd={n_fwd} ret={n_ret}")
    sys.stdout.flush()

    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)

    best_both_km: float = float("inf")
    best_key: tuple[float, float] | None = None
    any_both_converged = False

    header = (
        f"{'c_em':>8}  {'c_se':>10}  "
        f"{'fwd_km':>12}  {'fwd_ok':>6}  "
        f"{'ret_km':>12}  {'ret_ok':>6}  "
        f"{'both':>6}  note"
    )
    print(header)
    sys.stdout.flush()

    for c_em in C_EM_VALUES:
        for c_se in C_SE_VALUES:
            # Build EM-L2 and SE-L2 Lyapunov nodes.
            try:
                em_node = _build_em_node(em, "EM-L2", c_em)
                se_node = _build_se_node(se, "SE-L2", c_se)
            except Exception as exc:
                print(
                    f"[{_ts()}] {c_em:.4f}  {c_se:.7f}  "
                    f"{'---':>12}  {'---':>6}  {'---':>12}  {'---':>6}  "
                    f"{'SKIP':>6}  node build failed: {exc}"
                )
                sys.stdout.flush()
                continue

            if not em_node.converged:
                print(f"[{_ts()}] {c_em:.4f}  {c_se:.7f}  EM node failed")
                sys.stdout.flush()
                continue
            if not se_node.converged:
                print(f"[{_ts()}] {c_em:.4f}  {c_se:.7f}  SE node failed")
                sys.stdout.flush()
                continue

            # Forward leg: EM-L2 unstable -> SE-L2 stable.
            fwd = correct_cross_connection(
                bridge,
                em_node,
                se_node,
                label_from="EM-L2",
                label_to="SE-L2",
                tol_km=TOL_KM,
                scan_n=SCAN_N_FWD,
                scan_n_tau=SCAN_N_TAU_FWD,
                max_time_factor=MAX_TIME_FACTOR,
            )

            # Return leg: SE-L2 unstable -> EM-L2 stable.
            # Try all branch variants; keep the one with smallest position gap.
            ret_best = None
            for v in _RETURN_VARIANTS:
                ret_cand = correct_cross_connection(
                    bridge,
                    se_node,
                    em_node,
                    label_from="SE-L2",
                    label_to="EM-L2",
                    tol_km=TOL_KM,
                    branch_u=v["branch_u"],
                    branch_s=v["branch_s"],
                    scan_n=SCAN_N_RET,
                    scan_n_tau=SCAN_N_TAU_RET,
                    max_time_factor=MAX_TIME_FACTOR,
                )
                if ret_best is None or ret_cand.residual < ret_best.residual:
                    ret_best = ret_cand
                if ret_cand.converged:
                    break
            assert ret_best is not None

            both_ok = fwd.converged and ret_best.converged
            max_gap = max(fwd.residual, ret_best.residual)
            if both_ok:
                any_both_converged = True
                if max_gap < best_both_km:
                    best_both_km = max_gap
                    best_key = (c_em, c_se)

            note = (
                "BOTH CONVERGE"
                if both_ok
                else (
                    "fwd only"
                    if fwd.converged
                    else ("ret only" if ret_best.converged else "neither")
                )
            )
            print(
                f"[{_ts()}] {c_em:8.4f}  {c_se:10.7f}  "
                f"{fwd.residual:12.3e}  {fwd.converged!s:>6}  "
                f"{ret_best.residual:12.3e}  {ret_best.converged!s:>6}  "
                f"{both_ok!s:>6}  {note}"
            )
            sys.stdout.flush()

    print(f"\n[{_ts()}] === SUMMARY ===")
    if any_both_converged:
        print(
            f"[{_ts()}] BOTH-LEG CLOSURE FOUND: best at c_em={best_key[0]:.4f}, "
            f"c_se={best_key[1]:.7f}, max_leg_gap={best_both_km:.3e} km"
        )
        print(f"[{_ts()}] → feed this seed into feasibility_ls / bounded_ls driver")
    else:
        print(f"[{_ts()}] NO both-leg closure found across all (c_em, c_se) grid points.")
        print(f"[{_ts()}] Planar CR3BP wall confirmed → next step: 3D z-slicing (Gómez 2004)")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
