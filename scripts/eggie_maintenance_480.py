"""#480 — EGGIE maintenance-ΔV curve (real ephemeris, Approach A).

The remaining #480 quantitative direction: quantify the per-cycle MAINTENANCE ΔV of a
feasible real-ephemeris EGGIE (Europa-Ganymede-Ganymede-Io-Europa, 4-synodic/5-rev) member
by chaining it forward cycle-by-cycle in the real DE440/JUP365 ephemeris with the validated
per-cycle re-targeting lane (nbody.jovian.chain_cycles, validated on Liang Member D #223).

Honest framing: the paper prints NO EGGIE maintenance number (only "ballistic for ~2
cycles, then large impulses"); its one quantified maintenance figure (~30 m/s/10) is the
EIGE example, which the patched-conic lane cannot reproduce feasibly (the 1-rev B-plane
wall — see 2026-06-30-480-eige-realeph-maintenance-verdict.md). So this is a NOVEL
maintenance curve for a feasibly-discovered EGGIE member, method-validated by Liang #223,
NOT a paper reproduction. No catalogue impact.

Stage 1: scan branch-plan x departure epoch near the paper era for a feasible ballistic
EGGIE cycle-1 (V∞ an OUTPUT — no steering). Stage 2: chain n_cycles there; report the
maintenance ΔV per cycle + cumulative.
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from cyclerfinder.nbody import jovian
from cyclerfinder.search.eggie_ballistic import EGGIE_SEQUENCE
from cyclerfinder.search.resonant_conic import EGGIE_TOFS_TABLE4_DAYS

SECONDS_PER_DAY = 86400.0
MIN_ALT_KM = 50.0  # Liang optimizer floor (p.11); reported altitudes still checked vs window
ALT_MAX_KM = 70000.0  # paper feasibility ceiling (p.7)

# Candidate 4-leg (n_revs, branch) plans. Leg structure E->G (sub-rev), G->G, G->I, I->E
# (each ~1 spacecraft rev at the 4:5 resonant sma) — the ideal-model feasible structure
# (eggie_ballistic.FEASIBLE_BALLISTIC_PLAN) plus branch variations.
CANDIDATE_PLANS: tuple[tuple[tuple[int, str], ...], ...] = (
    ((0, "single"), (1, "high"), (1, "low"), (1, "high")),
    ((0, "single"), (1, "low"), (1, "high"), (1, "low")),
    ((0, "single"), (1, "high"), (1, "high"), (1, "high")),
    ((0, "single"), (1, "low"), (1, "low"), (1, "low")),
    ((0, "single"), (1, "high"), (1, "low"), (1, "low")),
    ((0, "single"), (1, "low"), (1, "high"), (1, "high")),
)
# Paper EGGIE departs 29-Sep / 02-Oct-2020; scan a few synodic periods around it.
SCAN_START_ISO = "2020-09-15T00:00:00"
SCAN_DAYS = 30.0


def _ts() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def eval_cycle1(
    t0: float, ephem: jovian.JovianEphemeris, plan: tuple[tuple[int, str], ...], bound_days: float
) -> tuple[float, bool, object]:
    """Optimised first-cycle (defect_ms, feasible, cycle) for a branch plan at epoch t0."""
    cyc, _ = jovian.optimize_cycle(
        t0,
        list(EGGIE_TOFS_TABLE4_DAYS),
        ephem,
        cycle_index=1,
        vinf_in_prev=None,
        bound_days=bound_days,
        min_alt_km=MIN_ALT_KM,
        sequence=EGGIE_SEQUENCE,
        branch_plan=plan,
    )
    finite_alts = [a for a in cyc.altitudes_km if np.isfinite(a)]
    feasible = (
        cyc.cycle_tof_days > 0
        and bool(finite_alts)
        and all(MIN_ALT_KM <= a <= ALT_MAX_KM for a in finite_alts)
    )
    return cyc.sum_defect_ms, feasible, cyc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-step-days", type=float, default=1.0)
    ap.add_argument("--n-cycles", type=int, default=10)
    ap.add_argument("--bound-days", type=float, default=4.0)
    args = ap.parse_args()

    kernel = jovian.jup365_kernel_path()
    if kernel is None:
        from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

        kernel = ensure_jup365_kernel()
    print(f"[{_ts()}] kernel = {kernel}", flush=True)
    ephem = jovian.JovianEphemeris(kernel)
    print(f"[{_ts()}] sequence = {EGGIE_SEQUENCE}", flush=True)
    print(f"[{_ts()}] ToF seed (Table-4, d) = {list(EGGIE_TOFS_TABLE4_DAYS)}", flush=True)

    # --- Stage 1: branch-plan x epoch scan for a feasible ballistic cycle-1 -----
    t_scan0 = jovian.tdb_sec_from_iso(SCAN_START_ISO)
    n_steps = int(SCAN_DAYS / args.scan_step_days) + 1
    print(
        f"\n[{_ts()}] STAGE 1: {len(CANDIDATE_PLANS)} plans x {n_steps} epochs "
        f"from {SCAN_START_ISO} (+{SCAN_DAYS} d)",
        flush=True,
    )
    best: tuple[float, float, tuple[tuple[int, str], ...]] | None = None
    best_feasible: tuple[float, float, tuple[tuple[int, str], ...]] | None = None
    t_start = time.monotonic()
    for p_idx, plan in enumerate(CANDIDATE_PLANS):
        for i in range(n_steps):
            t0 = t_scan0 + i * args.scan_step_days * SECONDS_PER_DAY
            defect, feasible, _ = eval_cycle1(t0, ephem, plan, args.bound_days)
            if best is None or defect < best[0]:
                best = (defect, t0, plan)
            if feasible and (best_feasible is None or defect < best_feasible[0]):
                best_feasible = (defect, t0, plan)
        bf = f"{best_feasible[0]:.3e}" if best_feasible else "none"
        print(f"  [{_ts()}] plan {p_idx} {plan} done  best_feasible_defect={bf}", flush=True)
    print(f"[{_ts()}] STAGE 1 done in {time.monotonic() - t_start:.1f}s", flush=True)

    chosen = best_feasible if best_feasible is not None else best
    assert chosen is not None
    label = "feasible-ballistic" if best_feasible is not None else "best-defect (NOT feasible)"
    defect0, t0_best, plan_best = chosen
    iso = "+%.2f d" % ((t0_best - t_scan0) / SECONDS_PER_DAY)
    print(
        f"\n[{_ts()}] chosen seed ({label}): defect={defect0:.3e} m/s  epoch {iso}  "
        f"plan={plan_best}",
        flush=True,
    )

    # --- Stage 2: chain n_cycles at the chosen seed -----------------------------
    print(f"\n[{_ts()}] STAGE 2: chain {args.n_cycles} EGGIE cycles", flush=True)
    t_chain = time.monotonic()
    cycles = jovian.chain_cycles(
        t0_best,
        ephem,
        n_cycles=args.n_cycles,
        tof_seed_days=EGGIE_TOFS_TABLE4_DAYS,
        bound_days=args.bound_days,
        min_alt_km=MIN_ALT_KM,
        sequence=EGGIE_SEQUENCE,
        branch_plan=plan_best,
        progress=True,
    )
    print(
        f"[{_ts()}] chained {len(cycles)} cycles in {time.monotonic() - t_chain:.1f}s", flush=True
    )

    print("\n=== EGGIE per-cycle maintenance ΔV (real JUP365, conic chain) ===", flush=True)
    print(
        f"{'cyc':>3} {'tof_d':>9} {'sum_dv_ms':>11} {'cumulative_ms':>13} "
        f"{'min_alt_km':>11} {'feasible':>9}"
    )
    cum = 0.0
    for c in cycles:
        cum += c.sum_defect_ms
        finite_alts = [a for a in c.altitudes_km if np.isfinite(a)]
        min_alt = min(finite_alts) if finite_alts else float("nan")
        feas = bool(finite_alts) and all(MIN_ALT_KM <= a <= ALT_MAX_KM for a in finite_alts)
        print(
            f"{c.index:>3} {c.cycle_tof_days:>9.4f} {c.sum_defect_ms:>11.3e} "
            f"{cum:>13.3e} {min_alt:>11.1f} {feas!s:>9}",
            flush=True,
        )
    print(
        f"\n[{_ts()}] cumulative maintenance ΔV over {len(cycles)} cycles = {cum:.3f} m/s "
        f"(NOVEL result; paper prints no EGGIE number — see module docstring)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
