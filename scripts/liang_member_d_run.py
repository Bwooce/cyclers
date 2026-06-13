"""Liang et al. 2024 Member D same-family n-body re-propagation lane (#223).

Drives the :mod:`cyclerfinder.nbody.jovian` toolkit end to end:

1. SEED   — phase the validated idealized CGCEC scaffold (Member A first-cycle
            leg ToFs, #222) to Member D's published Callisto departure epoch
            2033-09-25T18:04:43 (TDB assumed) on real JUP365 geometry.
2. CHAIN  — patched-conic CGCEC chain on real moon positions (multi-rev Lambert
            legs), per-cycle local optimisation of the flyby epochs to minimise
            the powered-flyby defect dv, chained cycle-by-cycle.
3. SHOOT  — re-propagate each converged conic cycle in the REBOUND restricted
            n-body model and attempt ballistic (continuity) closure.
4. COMPARE— qualitative signature vs the published traces: sequence achieved,
            per-cycle ToF inside the 99.4-100.5 d figure band, flybys above
            100 km altitude, residual-dv scale vs their 1.0383e-7 m/s.

HONEST CEILING (note docs/notes/2026-06-13-liang-member-d-nbody.md): this is a
same-family candidate of OUR construction, not their member. A close is at most
V1-class existence evidence, never V3. A clean failure to close is a finding.

Run (kernel not committed; ~1.14 GB):

    CYCLERFINDER_JUP365=~/dev/references/kernels/jup365.bsp \
        uv run python scripts/liang_member_d_run.py --n-cycles 3

If ``CYCLERFINDER_JUP365`` is unset this script falls back to the conventional
local path ``~/dev/references/kernels/jup365.bsp`` (still never committed).
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

# Conventional local kernel location (NOT committed; see the results note).
_DEFAULT_KERNEL = os.path.expanduser("~/dev/references/kernels/jup365.bsp")
if not os.environ.get("CYCLERFINDER_JUP365") and os.path.exists(_DEFAULT_KERNEL):
    os.environ["CYCLERFINDER_JUP365"] = _DEFAULT_KERNEL

import cyclerfinder.nbody.jovian as jovian  # noqa: E402
from cyclerfinder.core.constants import SECONDS_PER_DAY  # noqa: E402

# Member D published anchors (Liang et al. 2024, p. 19; figure traces 8d-e).
DEPART_ISO = "2033-09-25T18:04:43"
RETURN_ISO = "2036-06-22T01:44:39"
TOF_BAND_DAYS = (99.4, 100.5)
PUBLISHED_MAX_DV_MS = 1.0383e-7
MIN_ALT_KM = 100.0  # paper: all flybys above 100 km real altitude

# Idealized Member A first-cycle leg ToFs (Table 3) — the validated scaffold
# seed (#222), phased to the real epoch as the chain's per-cycle ToF seed.
TOF_SEED_DAYS = (31.8973, 18.1697, 29.9343, 19.9747)


def _fmt(xs: object, prec: int = 4) -> str:
    if isinstance(xs, (list, tuple)):
        return "[" + ", ".join(f"{float(v):.{prec}f}" for v in xs) + "]"
    return f"{float(xs):.{prec}f}"  # type: ignore[arg-type]


def _build_nodes(
    c: jovian.ConicCycle,
    ephem: jovian.JovianEphemeris,
) -> tuple[list[tuple[np.ndarray, np.ndarray]] | None, list[float]]:
    """Periapsis nodes + SOI-fraction encounter caps for one conic cycle."""
    legs = jovian._solve_cycle_legs(c.epochs_sec, ephem, jovian.BRANCH_PLAN)
    if legs is None:
        return None, []
    vinf_out, vinf_in = legs
    nodes: list[tuple[np.ndarray, np.ndarray]] = []
    # node 0: departure Callisto — only the outbound asymptote is known.
    r0, v0, _ = jovian.periapsis_node(
        jovian.CGCEC[0], c.epochs_sec[0], vinf_out[0], vinf_out[0], ephem
    )
    nodes.append((r0, v0))
    for k in range(1, 5):
        vo = vinf_out[k] if k < 4 else vinf_in[k - 1]
        r, v, _ = jovian.periapsis_node(jovian.CGCEC[k], c.epochs_sec[k], vinf_in[k - 1], vo, ephem)
        nodes.append((r, v))
    d_caps = [0.0] + [
        jovian.SATELLITES[jovian.CGCEC[k]].sma_km
        * (jovian.SATELLITES[jovian.CGCEC[k]].mu_km3_s2 / (3.0 * jovian.MU_JUPITER_KM3_S2))
        ** (1.0 / 3.0)
        for k in range(1, 5)
    ]
    return nodes, d_caps


def _report_seed_defects(
    c: jovian.ConicCycle,
    nodes: list[tuple[np.ndarray, np.ndarray]],
    d_caps: list[float],
    ephem: jovian.JovianEphemeris,
    cache: jovian.JovianRailsCache,
    vinf_prev_mag: float | None,
) -> None:
    """One n-body propagation per leg: the raw patched-conic -> n-body gap.

    This is the lane's headline measurement (no optimisation): how far the
    idealized patched-conic node states miss n-body continuity on the real
    moon-gravity model, per leg, in km (position) and m/s (velocity).
    """
    prop = jovian.JovianRestrictedNBody()
    leg_dr: list[float] = []
    leg_dv: list[float] = []
    for k in range(4):
        arc = prop.propagate(
            nodes[k][0], nodes[k][1], c.epochs_sec[k], c.epochs_sec[k + 1], cache=cache
        )
        if not arc.converged:
            leg_dr.append(float("nan"))
            leg_dv.append(float("nan"))
            continue
        leg_dr.append(float(np.linalg.norm(arc.r_km - nodes[k + 1][0])))
        leg_dv.append(float(np.linalg.norm(arc.v_km_s - nodes[k + 1][1])))
    print(f"cycle {c.index}: seed |dr|(km) ={_fmt(leg_dr, 1)}")
    print(f"          seed |dv|(m/s)={_fmt([d * 1e3 for d in leg_dv], 2)}")


def main() -> int:
    # Line-buffer stdout so per-cycle progress is visible during the long
    # n-body shooting (otherwise everything flushes only at process exit).
    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[union-attr]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cycles", type=int, default=3, help="conic cycles to chain")
    ap.add_argument("--bound-days", type=float, default=3.0, help="epoch search bound")
    ap.add_argument("--shoot", action="store_true", help="also run n-body shooting")
    ap.add_argument("--shoot-cycles", type=int, default=1, help="how many cycles to shoot")
    ap.add_argument("--shoot-max-nfev", type=int, default=60, help="least_squares nfev budget")
    ap.add_argument(
        "--diag",
        action="store_true",
        help="seed-only: report the patched-conic -> n-body leg defects without optimising",
    )
    args = ap.parse_args()

    kernel = jovian.jup365_kernel_path()
    if kernel is None:
        print(
            "ERROR: JUP365 kernel not found. Set CYCLERFINDER_JUP365 to the "
            f"jup365.bsp path (tried {_DEFAULT_KERNEL}).",
            file=sys.stderr,
        )
        return 2

    print(f"[setup] kernel = {kernel}")
    ephem = jovian.JovianEphemeris(kernel)
    t0 = jovian.tdb_sec_from_iso(DEPART_ISO)
    t_ret = jovian.tdb_sec_from_iso(RETURN_ISO)
    span_days = (t_ret - t0) / SECONDS_PER_DAY
    print(f"[setup] depart {DEPART_ISO} -> t0 = {t0:.3f} s (TDB since J2000)")
    print(f"[setup] return {RETURN_ISO} -> {t_ret:.3f} s  span = {span_days:.2f} d")

    # Departure-epoch moon geometry (sanity: real JUP365 placement).
    print("[setup] moon |r| (km) at departure epoch:")
    for m in jovian.GALILEAN:
        r, v = ephem.state(m, t0)
        print(f"          {m:9s} |r|={np.linalg.norm(r):11.2f}  |v|={np.linalg.norm(v):7.4f}")

    # --- 2. patched-conic CGCEC chain on real geometry --------------------------
    print(f"\n[chain] chaining {args.n_cycles} CGCEC cycles (seed ToF {_fmt(TOF_SEED_DAYS)} d) ...")
    t_chain = time.monotonic()
    cycles = jovian.chain_cycles(
        t0,
        ephem,
        n_cycles=args.n_cycles,
        tof_seed_days=TOF_SEED_DAYS,
        bound_days=args.bound_days,
        min_alt_km=MIN_ALT_KM,
        progress=True,
    )
    print(f"[chain] done in {time.monotonic() - t_chain:.1f} s ({len(cycles)} cycles)")

    # --- 4. signature comparison ------------------------------------------------
    print("\n=== PER-CYCLE SIGNATURE (conic chain on real JUP365) ===")
    hdr = (
        f"{'cyc':>3} {'tof_d':>9} {'in_band':>7} {'sum_dv_ms':>11} "
        f"{'max_dv_ms':>11} {'min_alt_km':>11} {'alts>100km':>10} {'conv':>5}"
    )
    print(hdr)
    print("-" * len(hdr))
    for c in cycles:
        in_band = TOF_BAND_DAYS[0] <= c.cycle_tof_days <= TOF_BAND_DAYS[1]
        max_dv = max(c.defects_ms) if c.defects_ms else 0.0
        finite_alts = [a for a in c.altitudes_km if np.isfinite(a)]
        min_alt = min(finite_alts) if finite_alts else float("nan")
        alts_ok = all(a > MIN_ALT_KM for a in finite_alts) if finite_alts else False
        print(
            f"{c.index:>3} {c.cycle_tof_days:>9.4f} {in_band!s:>7} "
            f"{c.sum_defect_ms:>11.3e} {max_dv:>11.3e} {min_alt:>11.1f} "
            f"{alts_ok!s:>10} {c.converged!s:>5}"
        )

    print("\n=== PER-CYCLE LEG DETAIL ===")
    for c in cycles:
        print(f"cycle {c.index}: leg ToFs(d)={_fmt(c.tofs_days)}  vinf(km/s)={_fmt(c.vinf_kms)}")
        print(f"          defects(m/s)={_fmt(c.defects_ms, 3)}  alts(km)={_fmt(c.altitudes_km, 1)}")

    # --- 3. optional n-body shooting / diagnostic ------------------------------
    if args.shoot or args.diag:
        mode = "DIAGNOSTIC (seed defects only)" if args.diag else "SHOOTING"
        print(f"\n=== N-BODY {mode} (REBOUND restricted model) ===")
        n_shoot = min(args.shoot_cycles, len(cycles))
        # One rails cache spanning the cycles we touch (plus pad inside cache).
        last = cycles[n_shoot - 1]
        cache = jovian.JovianRailsCache(jovian.GALILEAN, ephem, t0, last.epochs_sec[-1])
        vinf_prev_mag: float | None = None
        for c in cycles[:n_shoot]:
            nodes, d_caps = _build_nodes(c, ephem)
            if nodes is None:
                print(f"cycle {c.index}: conic legs unsolvable at converged epochs; skip")
                continue
            if args.diag:
                _report_seed_defects(c, nodes, d_caps, ephem, cache, vinf_prev_mag)
            else:
                t_shoot = time.monotonic()
                res = jovian.shoot_cycle(
                    nodes,
                    c.epochs_sec,
                    ephem,
                    cache,
                    cycle_index=c.index,
                    d_caps_km=d_caps,
                    vinf_in_mag_prev=vinf_prev_mag,
                    max_nfev=args.shoot_max_nfev,
                )
                print(
                    f"cycle {c.index}: converged={res.converged} nfev={res.nfev} "
                    f"({time.monotonic() - t_shoot:.1f} s)"
                )
                fin_dv_ms = [d * 1e3 for d in res.final_leg_defects_kms]
                print(f"          seed |dr|(km)  ={_fmt(res.seed_leg_defects_km, 3)}")
                print(f"          final |dr|(km) ={_fmt(res.final_leg_defects_km, 6)}")
                print(f"          final |dv|(m/s)={_fmt(fin_dv_ms, 6)}")
                print(f"          moon dist(km)  ={_fmt(res.moon_distances_km, 1)}")
                print(f"          boundary dv(m/s)={res.boundary_dv_ms:.6e}")
            vinf_prev_mag = c.vinf_kms[-1]

    # --- verdict scaffolding ----------------------------------------------------
    seq_ok = all(c.cycle_tof_days > 0.0 for c in cycles) and len(cycles) >= 1
    bands = [TOF_BAND_DAYS[0] <= c.cycle_tof_days <= TOF_BAND_DAYS[1] for c in cycles]
    print("\n=== SUMMARY ===")
    print(f"cycles run        : {len(cycles)} (requested {args.n_cycles})")
    print(f"sequence achieved : {seq_ok} (CGCEC topology constructed each cycle)")
    print(f"ToF in 99.4-100.5 : {sum(bands)}/{len(cycles)} cycles")
    print(f"published max dv   : {PUBLISHED_MAX_DV_MS:.4e} m/s")
    best_dv = min((c.sum_defect_ms for c in cycles), default=float("nan"))
    print(f"best conic sum_dv  : {best_dv:.3e} m/s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
