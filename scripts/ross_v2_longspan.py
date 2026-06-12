"""Long-span DISCRIMINATING inertial propagation for the 5 Ross rows (#229).

Executes the #224 recommendation (docs/notes/2026-06-13-ross-v2-assessment.md):
the campaign's 5-period inertial gate is weakly discriminating for
residual-zero ICs (seed error delta_1 ~ 1e-9..1e-8 cannot amplify past the 3A
band in 5 periods even for |lambda| ~ 10^2..10^3). This script extends the SAME
inertial REBOUND/IAS15 harness (`cr3bp-inertial-rebound-ias15-v1`, reused
verbatim from scripts/cr3bp_family_search.py -- frame conversion, Jacobi
bookkeeping, integrator settings all unchanged) to a span where the instrument
has teeth:

    T_span >= T * ln(3A / delta_1) / ln(|lambda_hypo|),  |lambda_hypo| = 2

evaluated with the ACTUAL measured delta_1 per row. With delta_1 ~ 1e-9 and
A ~ O(1) nd that is ~32 periods; we run N_PERIODS = 100 so even a hypothetical
|lambda| ~ 2 instability would visibly amplify past the 3A departure band.

Per row this records: measured delta_1, span used + per-row justification,
max per-period rotating-frame recurrence drift (first half vs second half of
the span -- bounded oscillation vs secular growth), Jacobi drift over the span,
and a bounded-band verdict.

MODEL SCOPE (like-for-like): the harness integrates the CR3BP idealisation in
inertial coordinates (primaries on the exact circular rail, massless
spacecraft) -- the rows' DEFINING model, NOT a real-ephemeris claim. V2-ballistic
is this lane's ceiling.

Output: human-readable progress + a markdown per-row table + a JSON blob on
stdout. NO catalogue writeback (held for user review, #229).

Usage:
    uv run python scripts/ross_v2_longspan.py
    uv run python scripts/ross_v2_longspan.py --n-periods 100
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp

# Reuse the campaign harness verbatim (seeds + inertial cross-check machinery).
# Import via the repo-root namespace package so the module identity matches the
# test suite's ``import scripts.cr3bp_family_search``.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scripts.cr3bp_family_search as fs

METHOD_TAG = "ross-v2-longspan-v1"
DEFAULT_N_PERIODS = 100
LAMBDA_HYPO = 2.0  # conservative hypothetical unstable multiplier (assessment note)
# Secular-growth gate: a stable member's recurrence drift oscillates in a band;
# slow linear phase drift gives second-half/first-half max ratio ~2, while a
# hypothetical |lambda|=2 instability over 50 periods gives ~2^50. 10 separates
# them by many orders of magnitude without loosening the band.
HALF_RATIO_MAX = 10.0

ROW_IDS = {
    "ross-(1,1)": "ross-rt-em-cycler-11-2025",
    "ross-(2,1)": "ross-rt-em-cycler-21-2025",
    "ross-(3,1)": "ross-rt-em-cycler-31-2025",
    "ross-(3,2)": "ross-rt-em-cycler-32-2025",
    "ross-(3,3)": "ross-rt-em-cycler-33-2025",
}


@dataclass(frozen=True)
class LongSpanResult:
    """Per-row long-span discriminating-run record."""

    row_id: str
    label: str
    jacobi: float
    period_nd: float
    nu: float
    amplitude_nd: float
    delta1_nd: float
    noise_floor_nd: float
    n_periods: int
    n_periods_required: float  # ln(3A/delta1)/ln(2) from the MEASURED delta1
    max_drift_first_half_nd: float
    max_drift_second_half_nd: float
    half_ratio: float
    max_drift_span_nd: float
    departure_band_nd: float  # 3A
    jacobi_drift_span: float
    diverged: bool
    verdict: str  # "BOUNDED" | "NOT-BOUNDED"
    elapsed_s: float


def run_row(seed: fs.FamilySeed, n_periods: int) -> LongSpanResult:
    """Correct the member, measure nu/A, run the long-span inertial propagation."""
    t0 = time.monotonic()
    system = cr3bp.CR3BPSystem(
        mu=seed.mu,
        primary=seed.primary,
        secondary=seed.secondary,
        l_km=seed.l_km,
        t_s=seed.t_s,
    )
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        seed.x0_seed,
        seed.jacobi,
        seed.period_guess,
        ydot0_sign=seed.ydot0_sign,
        half_crossings=seed.half_crossings,
        tol=1e-10,
    )
    if not orbit.converged:
        raise RuntimeError(f"{seed.label}: corrector did not converge")
    nu, lam = cp.barden_stability(system, orbit)
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    max_speed, amplitude = cc._gate_metrics(system, state0, orbit.period, rtol=1e-12, atol=1e-12)
    member = cc.BranchMember(
        state0=state0,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        nu=float(nu),
        abs_lambda=float(abs(lam)),
        crossing_residual=orbit.crossing_residual,
        radau_djacobi=0.0,  # not re-derived here; campaign artifact carries it
        max_speed_nd=max_speed,
        amplitude_nd=amplitude,
        n_iter=orbit.n_iter,
        stable=abs(float(nu)) < 1.0,
    )

    # The campaign harness, verbatim, at the long span. Its 5-period R1/R2/R3
    # verdict string is NOT used here -- only the measured metrics.
    iv = fs.inertial_crosscheck(member, seed, n_periods=n_periods)

    deltas = np.asarray(iv.delta_per_period_nd, dtype=float)
    finite = np.isfinite(deltas)
    half = n_periods // 2
    first = deltas[:half][finite[:half]]
    second = deltas[half:][finite[half:]]
    max_first = float(np.max(first)) if first.size else float("nan")
    max_second = float(np.max(second)) if second.size else float("nan")
    half_ratio = max_second / max_first if max_first > 0.0 else float("inf")
    max_span = float(np.max(deltas[finite])) if finite.any() else float("nan")

    band = 3.0 * amplitude
    delta1 = iv.delta1_nd
    n_req = math.log(band / delta1) / math.log(LAMBDA_HYPO) if 0.0 < delta1 < band else float("nan")

    bounded = (
        not iv.diverged and bool(finite.all()) and max_span <= band and half_ratio <= HALF_RATIO_MAX
    )
    return LongSpanResult(
        row_id=ROW_IDS[seed.label],
        label=seed.label,
        jacobi=seed.jacobi,
        period_nd=member.period,
        nu=float(nu),
        amplitude_nd=amplitude,
        delta1_nd=delta1,
        noise_floor_nd=iv.noise_floor_nd,
        n_periods=n_periods,
        n_periods_required=n_req,
        max_drift_first_half_nd=max_first,
        max_drift_second_half_nd=max_second,
        half_ratio=half_ratio,
        max_drift_span_nd=max_span,
        departure_band_nd=band,
        jacobi_drift_span=iv.jacobi_drift_bounded,
        diverged=iv.diverged,
        verdict="BOUNDED" if bounded else "NOT-BOUNDED",
        elapsed_s=time.monotonic() - t0,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-periods", type=int, default=DEFAULT_N_PERIODS)
    args = ap.parse_args()

    print(f"# {METHOD_TAG}: long-span discriminating inertial run (#229)")
    print(
        f"# harness: {fs.INERTIAL_TAG} (REBOUND/IAS15, epsilon={fs.IAS15_EPSILON},"
        f" {fs.N_SAMP_PER_PERIOD} samples/period); span = {args.n_periods} periods"
    )
    results: list[LongSpanResult] = []
    for seed in fs.ROSS_SEEDS:
        print(f"[{time.strftime('%H:%M:%S')}] {seed.label} ...", flush=True)
        res = run_row(seed, args.n_periods)
        results.append(res)
        print(
            f"  -> {res.verdict}: delta1={res.delta1_nd:.3e}, N_req={res.n_periods_required:.1f},"
            f" max drift {res.max_drift_first_half_nd:.3e}/{res.max_drift_second_half_nd:.3e}"
            f" (1st/2nd half, ratio {res.half_ratio:.2f}), dJ={res.jacobi_drift_span:.3e},"
            f" {res.elapsed_s:.1f}s",
            flush=True,
        )

    print("\n## Per-row table (markdown)\n")
    print(
        "| row | family | C | T (nd) | nu | A (nd) | delta_1 | N_req | span "
        "| max drift 1st half | max drift 2nd half | ratio | dJ span | verdict |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        print(
            f"| `{r.row_id}` | {r.label} | {r.jacobi:.15f} | {r.period_nd:.8f} "
            f"| {r.nu:+.5f} | {r.amplitude_nd:.4f} | {r.delta1_nd:.3e} "
            f"| {r.n_periods_required:.1f} | {r.n_periods} "
            f"| {r.max_drift_first_half_nd:.3e} | {r.max_drift_second_half_nd:.3e} "
            f"| {r.half_ratio:.2f} | {r.jacobi_drift_span:.3e} | {r.verdict} |"
        )

    print("\n## JSON\n")
    print(json.dumps([r.__dict__ for r in results], indent=2))


if __name__ == "__main__":
    main()
