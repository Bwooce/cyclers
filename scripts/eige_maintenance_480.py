"""#480 — EIGE maintenance-ΔV positive control (real ephemeris).

Drive the per-cycle chained maintenance lane (nbody.jovian.chain_cycles) on the EIGE
(Europa-Io-Ganymede-Europa, 1-synodic/1-rev) tour in the real DE440/JUP365 ephemeris,
seeded by the ideal-model ballistic construction (search.eige_ballistic). The paper's
Fig-5 EIGE is ballistic for the first cycle then grows to ~30 m/s over 10 repeat cycles
(AAS 17-608 pp.10-11). This reproduces that as the positive control for the maintenance
method BEFORE quoting an EGGIE number ([[feedback_verify_gauntlet_with_positive_control]]).

Stage 1: scan the departure epoch over a synodic period to find a near-ballistic first
cycle (the paper optimises cycle-1 to ballistic). Stage 2: chain 10 cycles there and
report the maintenance-ΔV curve. The paper prints no EIGE epoch, so the control is
class-level: PASS iff cycle-1 is ~ballistic and the 10-cycle growth is the ~30 m/s order.
"""

from __future__ import annotations

import argparse
import sys
import time

from cyclerfinder.nbody import jovian
from cyclerfinder.search.eige_ballistic import EIGE_SEQUENCE, eige_tof_seed_days
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

SECONDS_PER_DAY = 86400.0
# EIGE: 3 sub-revolution legs, each a single-rev (n_revs=0) Lambert arc.
EIGE_BRANCH_PLAN = ((0, "single"), (0, "single"), (0, "single"))
# Paper-era scan window (EGGIE/EGIEIE depart late Sep / early Oct 2020); EIGE epoch
# is unprinted, so we scan a full synodic period for a near-ballistic first cycle.
SCAN_START_ISO = "2020-09-20T00:00:00"
SCAN_DAYS = 8.0  # one synodic period (~7.05 d) + margin
MIN_ALT_KM = 25.0  # paper window floor (p.7)


def _ts() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def cycle1_defect(
    t0: float, ephem: jovian.JovianEphemeris, seed: tuple[float, ...], bound_days: float
) -> float:
    """Optimised first-cycle summed defect (m/s) at departure epoch ``t0``."""
    cyc, _ = jovian.optimize_cycle(
        t0,
        list(seed),
        ephem,
        cycle_index=1,
        vinf_in_prev=None,
        bound_days=bound_days,
        min_alt_km=MIN_ALT_KM,
        sequence=EIGE_SEQUENCE,
        branch_plan=EIGE_BRANCH_PLAN,
    )
    return cyc.sum_defect_ms if cyc.cycle_tof_days > 0 else float("inf")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-step-days", type=float, default=0.1)
    ap.add_argument("--n-cycles", type=int, default=10)
    ap.add_argument("--bound-days", type=float, default=1.0)
    args = ap.parse_args()

    kernel = ensure_jup365_kernel()
    print(f"[{_ts()}] kernel = {kernel}", flush=True)
    ephem = jovian.JovianEphemeris(kernel)
    seed = eige_tof_seed_days()
    seed_str = [f"{t:.4f}" for t in seed]
    print(f"[{_ts()}] EIGE ToF seed (d) = {seed_str}  sum={sum(seed):.4f}", flush=True)
    print(f"[{_ts()}] sequence = {EIGE_SEQUENCE}  branch_plan = {EIGE_BRANCH_PLAN}", flush=True)

    # --- Stage 1: epoch scan for a near-ballistic first cycle -------------------
    t_scan0 = jovian.tdb_sec_from_iso(SCAN_START_ISO)
    n_steps = int(SCAN_DAYS / args.scan_step_days) + 1
    print(
        f"\n[{_ts()}] STAGE 1: scan {n_steps} epochs from {SCAN_START_ISO} "
        f"(+{SCAN_DAYS} d, step {args.scan_step_days} d)",
        flush=True,
    )
    best_t0 = None
    best_def = float("inf")
    t_start = time.monotonic()
    for i in range(n_steps):
        t0 = t_scan0 + i * args.scan_step_days * SECONDS_PER_DAY
        d = cycle1_defect(t0, ephem, seed, args.bound_days)
        if d < best_def:
            best_def = d
            best_t0 = t0
        if i % 10 == 0:
            print(
                f"  [{_ts()}] step {i}/{n_steps}  +{i * args.scan_step_days:.2f}d  "
                f"defect={d:.3e} m/s  (best {best_def:.3e})",
                flush=True,
            )
    print(
        f"[{_ts()}] STAGE 1 done in {time.monotonic() - t_start:.1f}s. "
        f"best cycle-1 defect = {best_def:.3e} m/s at t0={best_t0:.1f}s "
        f"(+{(best_t0 - t_scan0) / SECONDS_PER_DAY:.3f} d)",
        flush=True,
    )

    # --- Stage 2: chain n_cycles at the best epoch ------------------------------
    print(f"\n[{_ts()}] STAGE 2: chain {args.n_cycles} EIGE cycles at best epoch", flush=True)
    t_chain = time.monotonic()
    cycles = jovian.chain_cycles(
        best_t0,
        ephem,
        n_cycles=args.n_cycles,
        tof_seed_days=seed,
        bound_days=args.bound_days,
        min_alt_km=MIN_ALT_KM,
        sequence=EIGE_SEQUENCE,
        branch_plan=EIGE_BRANCH_PLAN,
        progress=True,
    )
    print(
        f"[{_ts()}] chained {len(cycles)} cycles in {time.monotonic() - t_chain:.1f}s", flush=True
    )

    print("\n=== EIGE per-cycle maintenance ΔV (real JUP365, conic chain) ===", flush=True)
    hdr = f"{'cyc':>3} {'tof_d':>9} {'sum_dv_ms':>11} {'cumulative_ms':>13} {'vinf(km/s)':>40}"
    print(hdr)
    cum = 0.0
    for c in cycles:
        cum += c.sum_defect_ms
        vs = " ".join(f"{v:.2f}" for v in c.vinf_kms)
        print(
            f"{c.index:>3} {c.cycle_tof_days:>9.4f} {c.sum_defect_ms:>11.3e} {cum:>13.3e}   {vs}",
            flush=True,
        )
    print(
        f"\n[{_ts()}] cumulative maintenance ΔV over {len(cycles)} cycles = {cum:.3f} m/s "
        f"(paper EIGE: ~30 m/s over 10)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
