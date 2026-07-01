"""#512 — (n_em, n_se) integer resonance sweep for the #411 cross-system cycle.

Follow-up to #496 (docs/notes/2026-07-01-496-feasibility-first-verdict.md), which
CONFIRMED the phase-closure wall at n_em=1, n_se=1: R2 = -0.72 rad at the seed
(c_em=3.150, c_se=3.00086) and cannot be zeroed within the EM-L2 family (c_em max
~3.152-3.153) or below the SE-L2 manifold convergence threshold (c_se ~3.000854).

#496's next-steps note observes that each n_em/n_se increment shifts R2/R1 by a
FIXED phase offset (omega_rel * T_em / omega_rel * T_se respectively), independent
of the amplitude knobs (c_em, c_se) -- because T_em, T_se and omega_rel are properties
of the *orbits themselves* (built once per (c_em, c_se) in `_solve`), not of the
phase bookkeeping in `_resid`. This means the analytic wrap-shift prediction is
EXACT (not approximate) at any fixed seed:

    R1(n_se) = wrap[ R1(1) - omega_rel*(n_se-1)*T_se ]   (independent of n_em)
    R2(n_em) = wrap[ R2(1) - omega_rel*(n_em-1)*T_em ]   (independent of n_se)

because wrap[x - c] = wrap[wrap(x) - c] for any x, c (x and wrap(x) differ by an
integer multiple of 2*pi). So the (n_em, n_se) sweep decouples into two independent
1-D scans at the fixed seed -- no re-solving needed to RANK candidates, only to
verify the top ones actually close after re-optimizing (c_em, c_se).

Step 1 (cheap, this script's default action): print the analytic wrap-prediction
table for the full practical grid (n_em in 1..8, n_se in 1..4) plus a wider
diagnostic scan (n_em, n_se up to 120/60) to locate the true nearest-to-zero
crossings for completeness.

Step 2 (compute-heavy, gated by --solve): run `correct_cross_cycle(solver=
"bounded_ls")` -- the #496 fix, using its bounds -- at the top analytic candidates
NOT already characterized by #496 ((1,1) is already fully characterized: best
|R|=0.517-0.518 rad, closed=False; see the #496 verdict note; reused here rather
than re-run, per the "pull existing logs before rerunning known failures" policy).

Run (background, runlog):
  uv run python scripts/sweep_411_resonance.py --solve 2>&1 | tee runlogs/sweep_411_resonance.log
"""

from __future__ import annotations

import argparse
import datetime as _dt
import math
import sys

from cyclerfinder.genome.cross_system_cycle import (
    FrameBridge,
    correct_cross_cycle,
    em_moon_system,
    se_earth_system,
)

C_EM0 = 3.150
C_SE0 = 3.00086
# #496 bounded_ls family bounds (kept identical so results are comparable).
C_EM_BOUNDS = (3.112, 3.152)
C_SE_BOUNDS = (3.00050, 3.00086)

N_EM_MAX_PRACTICAL = 8
N_SE_MAX_PRACTICAL = 4
N_EM_MAX_WIDE = 120
N_SE_MAX_WIDE = 60

# Top candidates (beyond the already-characterized (1,1)) to actually re-solve with
# bounded_ls, ranked by the analytic |R| prediction over the practical grid.
CANDIDATES_TO_SOLVE = [(1, 2), (3, 1), (1, 3)]

MAX_ITER = 8  # max_nfev = max(8, max_iter*4) = 32 per candidate (bounded_ls, #496 pace)


def _ts() -> str:
    return _dt.datetime.now(tz=_dt.UTC).isoformat(timespec="seconds")


def _wrap(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", flush=True)


def analytic_baseline(bridge: FrameBridge) -> tuple[float, float, float, float, float]:
    """Exact seed R1(1), R2(1), omega_rel, T_em (s), T_se (s) via one cheap solve.

    Uses max_iter=0 so `correct_cross_cycle` (solver="newton") evaluates the seed
    residual ONCE and returns without iterating -- no optimization, just the exact
    (r1, r2) at (c_em0, c_se0, n_em=1, n_se=1), plus the orbit periods needed to
    extrapolate to other (n_em, n_se) analytically.
    """
    cyc = correct_cross_cycle(
        bridge,
        em_lib="EM-L2",
        se_lib="SE-L2",
        c_em0=C_EM0,
        c_se0=C_SE0,
        n_em=1,
        n_se=1,
        max_iter=0,
        solver="newton",
        return_scan_n=8,
        return_scan_n_tau=3,
    )
    em, se = bridge.em, bridge.se
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s
    # Get T_em, T_se independently from the node periods (exact, no re-solve), then
    # cross-check against cycle_time_s = t_fwd + n_se*T_se + t_ret + n_em*T_em (n=1 here).
    from cyclerfinder.genome.cross_system_cycle import _build_em_node, _build_se_node

    em_node = _build_em_node(em, "EM-L2", C_EM0)
    se_node = _build_se_node(se, "SE-L2", C_SE0)
    t_em = em_node.period * em.t_s
    t_se = se_node.period * se.t_s
    _cross_check = (
        cyc.forward.transit_time + 1 * t_se + cyc.ret.transit_time + 1 * t_em - cyc.cycle_time_s
    )
    assert abs(_cross_check) < 1.0, f"t_em/t_se reconstruction mismatch: {_cross_check} s"
    return cyc.r1_rad, cyc.r2_rad, omega_rel, t_em, t_se


def print_table(r1_1: float, r2_1: float, omega_rel: float, t_em: float, t_se: float) -> None:
    _log(f"seed R1(n_se=1)={r1_1:+.5f} rad  R2(n_em=1)={r2_1:+.5f} rad")
    _log(
        f"omega_rel={omega_rel:.6e} rad/s  T_em={t_em:.1f} s ({t_em / 86400:.3f} d)  "
        f"T_se={t_se:.1f} s ({t_se / 86400:.3f} d)"
    )
    _log(
        f"per-step shift: R1 step = {_wrap(-omega_rel * t_se):+.5f} rad/n_se  "
        f"R2 step = {_wrap(-omega_rel * t_em):+.5f} rad/n_em"
    )

    r1_of = {n: _wrap(r1_1 - omega_rel * (n - 1) * t_se) for n in range(1, N_SE_MAX_PRACTICAL + 1)}
    r2_of = {n: _wrap(r2_1 - omega_rel * (n - 1) * t_em) for n in range(1, N_EM_MAX_PRACTICAL + 1)}

    rows = []
    for n_em in range(1, N_EM_MAX_PRACTICAL + 1):
        for n_se in range(1, N_SE_MAX_PRACTICAL + 1):
            r1, r2 = r1_of[n_se], r2_of[n_em]
            rows.append((math.hypot(r1, r2), n_em, n_se, r1, r2))
    rows.sort(key=lambda t: t[0])

    _log(
        f"--- analytic wrap-prediction table, practical grid "
        f"n_em=1..{N_EM_MAX_PRACTICAL} x n_se=1..{N_SE_MAX_PRACTICAL} "
        f"({len(rows)} pairs), sorted by predicted |R| ---"
    )
    for mag, n_em, n_se, r1, r2 in rows:
        _log(f"  n_em={n_em}  n_se={n_se}  R1={r1:+.4f}  R2={r2:+.4f}  |R|={mag:.4f} rad")

    # Wide diagnostic scan: where are the true nearest-to-zero crossings for each
    # independent dimension (informational only -- these n are impractically large;
    # NOT solved by this driver).
    r1_wide = sorted(
        (abs(_wrap(r1_1 - omega_rel * (n - 1) * t_se)), n) for n in range(1, N_SE_MAX_WIDE + 1)
    )
    r2_wide = sorted(
        (abs(_wrap(r2_1 - omega_rel * (n - 1) * t_em)), n) for n in range(1, N_EM_MAX_WIDE + 1)
    )
    _log(f"--- wide scan (n_se up to {N_SE_MAX_WIDE}): best 5 |R1| minima ---")
    for a, n in r1_wide[:5]:
        _log(f"  n_se={n}  |R1|={a:.5f} rad  (cycle SE-dwell = {n * t_se / 86400:.1f} d)")
    _log(f"--- wide scan (n_em up to {N_EM_MAX_WIDE}): best 5 |R2| minima ---")
    for a, n in r2_wide[:5]:
        _log(f"  n_em={n}  |R2|={a:.5f} rad  (cycle EM-dwell = {n * t_em / 86400:.1f} d)")


def solve_candidate(bridge: FrameBridge, n_em: int, n_se: int) -> None:
    _log(f"=== solving n_em={n_em}, n_se={n_se} (bounded_ls, #496 bounds) ===")

    def _on_iter(k: int, c_em: float, c_se: float, r1: float, r2: float) -> None:
        _log(
            f"  n_em={n_em} n_se={n_se}  it={k:02d}  c_em={c_em:.6f}  c_se={c_se:.9f}  "
            f"R=({r1:+.4f},{r2:+.4f})  |R|={math.hypot(r1, r2):.4e}"
        )

    cyc = correct_cross_cycle(
        bridge,
        em_lib="EM-L2",
        se_lib="SE-L2",
        c_em0=C_EM0,
        c_se0=C_SE0,
        n_em=n_em,
        n_se=n_se,
        solver="bounded_ls",
        c_em_bounds=C_EM_BOUNDS,
        c_se_bounds=C_SE_BOUNDS,
        max_iter=MAX_ITER,
        return_scan_n=8,
        return_scan_n_tau=3,
        on_iter=_on_iter,
    )
    _log(
        f"RESULT n_em={n_em} n_se={n_se}: closed={cyc.closed} "
        f"|R|={cyc.theta_residual_norm:.4e} rad  R1={cyc.r1_rad:+.4f}  R2={cyc.r2_rad:+.4f}  "
        f"max_leg_residual_km={cyc.max_leg_residual_km:.2e}  "
        f"patch_dv_kms={cyc.total_patch_dv_kms:.4f}  "
        f"c_em={cyc.c_em:.6f}  c_se={cyc.c_se:.9f}  n_iter={cyc.n_iter}  notes={cyc.notes!r}"
    )
    if cyc.closed:
        _log(f"*** #411 CLOSE at n_em={n_em}, n_se={n_se} -- running crosscheck ***")
        try:
            from cyclerfinder.genome.cross_system_cycle import CrossCycle, crosscheck_cross_cycle

            as_cycle = CrossCycle(
                connections=[cyc.forward, cyc.ret],
                c_em=cyc.c_em,
                c_se=cyc.c_se,
                libration_pair=cyc.libration_pair,
                theta_closure_residual=cyc.theta_residual_norm,
                closed=cyc.closed,
                max_leg_residual=cyc.max_leg_residual_km,
                independent_residual=float("nan"),
            )
            checked = crosscheck_cross_cycle(bridge, as_cycle)
            _log(f"crosscheck independent_residual={checked.independent_residual:.4e} km")
        except Exception as exc:  # pragma: no cover - diagnostic path only
            _log(f"crosscheck FAILED: {exc!r}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solve", action="store_true", help="also run bounded_ls on top candidates")
    ap.add_argument(
        "--candidates",
        type=str,
        default=None,
        help="override candidate list, e.g. '1,2;3,1' (n_em,n_se pairs, ; separated)",
    )
    args = ap.parse_args()

    _log("#512 (n_em, n_se) resonance sweep for #411 cross-system cycle")
    bridge = FrameBridge(se=se_earth_system(), em=em_moon_system())

    _log("--- Step 1: analytic seed baseline + wrap-prediction table ---")
    r1_1, r2_1, omega_rel, t_em, t_se = analytic_baseline(bridge)
    print_table(r1_1, r2_1, omega_rel, t_em, t_se)

    if not args.solve:
        _log("(--solve not given: analytic table only. Exiting.)")
        return

    candidates = CANDIDATES_TO_SOLVE
    if args.candidates:
        candidates = []
        for pair in args.candidates.split(";"):
            a, b = pair.split(",")
            candidates.append((int(a), int(b)))

    _log(
        f"--- Step 2: re-solving top candidates {candidates} with bounded_ls "
        f"(NOTE: (1,1) already characterized by #496 -- not re-run) ---"
    )
    for n_em, n_se in candidates:
        solve_candidate(bridge, n_em, n_se)

    _log("done")


if __name__ == "__main__":
    sys.exit(main() or 0)
