"""Long-span DISCRIMINATING inertial propagation for the PC (3,2) row (#505).

Executes the V2-ballistic evidence run for catalogue row
``ross-rt-pc-cycler-32-2026`` (Pluto-Charon (3,2) CR3BP cycler, mu=0.10876,
V1 since #494).  Architecture mirrors ``scripts/ross_v2_longspan.py`` (the EM
five-row evidence package, #229) verbatim -- same harness, same discrimination
criterion, same bounded-band verdict gates:

    T_span >= T * ln(3A / delta_1) / ln(|lambda_hypo|),  |lambda_hypo| = 2

evaluated with the ACTUAL measured delta_1.  With delta_1 ~ 3.5e-10 and
A ~ 1.64 nd that is N_req ~ 33.7 periods; the default run uses N_PERIODS = 100
(~3x the requirement).

MODEL SCOPE (like-for-like): the harness integrates the CR3BP idealisation in
inertial coordinates -- the row's DEFINING model.  NOT a real-ephemeris claim.
V2-ballistic is this lane's ceiling.

NO catalogue writeback (held for user adjudication, #505).

Usage:
    uv run python scripts/pc_v2_longspan.py
    uv run python scripts/pc_v2_longspan.py --n-periods 100
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scripts.cr3bp_family_search as fs

ROW_ID = "ross-rt-pc-cycler-32-2026"
METHOD_TAG = "pc-v2-longspan-v1"
DEFAULT_N_PERIODS = 100

# Conservative hypothetical unstable multiplier (same as the EM evidence package)
LAMBDA_HYPO = 2.0
# Secular-growth gate: bounded oscillation gives half-ratio ~ 2; exponential
# |lambda|=2 over 50 periods gives ~2^50.  10 separates them by >14 orders.
HALF_RATIO_MAX = 10.0

# Catalogue parameters for ross-rt-pc-cycler-32-2026
# (state_nd DERIVED from #494 C-sweep nu=0 corrector at mu=0.10876)
PC_MU = 0.10876473603280369  # SOURCED: 106.1/975.5, satellites.py
PC_C = 3.57951501972907  # DERIVED: nu=0 midpoint C
PC_X0 = -0.693198287043369  # DERIVED: corrected x0
PC_T_GUESS = 11.8334625170346  # DERIVED: period in TU

# Corrector settings (from #494 Phase-3b and #504 confirmation)
_HALF_CROSSINGS = 6
_YDOT0_SIGN = -1.0


@dataclass(frozen=True)
class PCLongSpanResult:
    """Result record for the PC (3,2) long-span discriminating run."""

    row_id: str
    jacobi: float
    period_nd: float
    nu: float
    amplitude_nd: float
    delta1_nd: float
    noise_floor_nd: float
    n_periods: int
    n_periods_required: float
    max_drift_first_half_nd: float
    max_drift_second_half_nd: float
    half_ratio: float
    max_drift_span_nd: float
    departure_band_nd: float
    jacobi_drift_span: float
    diverged: bool
    verdict: str
    elapsed_s: float


def run_pc_longspan(n_periods: int = DEFAULT_N_PERIODS) -> PCLongSpanResult:
    """Correct the PC (3,2) IC, verify positive control, run long-span inertial."""
    t0 = time.monotonic()

    system = cr3bp.cr3bp_system("Pluto", "Charon")

    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        PC_X0,
        PC_C,
        PC_T_GUESS,
        ydot0_sign=_YDOT0_SIGN,
        half_crossings=_HALF_CROSSINGS,
        tol=1e-10,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"PC (3,2) corrector did not converge (res={orbit.crossing_residual:.2e})"
        )

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
        radau_djacobi=0.0,
        max_speed_nd=max_speed,
        amplitude_nd=amplitude,
        n_iter=orbit.n_iter,
        stable=abs(float(nu)) < 1.0,
    )

    # Build a FamilySeed for inertial_crosscheck (only mu/l_km/t_s are read,
    # but the type annotation requires the full FamilySeed dataclass).
    seed = fs.FamilySeed(
        label="pc-(3,2)",
        kind="pc",
        mu=system.mu,
        l_km=system.l_km,
        t_s=system.t_s,
        primary="Pluto",
        secondary="Charon",
        x0_seed=PC_X0,
        jacobi=PC_C,
        period_guess=PC_T_GUESS,
        ydot0_sign=_YDOT0_SIGN,
        half_crossings=_HALF_CROSSINGS,
        note="Catalogue row ross-rt-pc-cycler-32-2026, nu=0 midpoint",
    )

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

    return PCLongSpanResult(
        row_id=ROW_ID,
        jacobi=orbit.jacobi,
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

    print(f"# {METHOD_TAG}: PC (3,2) long-span discriminating inertial run (#505)")
    print(
        f"# harness: {fs.INERTIAL_TAG} (REBOUND/IAS15, epsilon={fs.IAS15_EPSILON},"
        f" {fs.N_SAMP_PER_PERIOD} samples/period); span = {args.n_periods} periods"
    )
    print(f"# row: {ROW_ID}")
    print(f"[{time.strftime('%H:%M:%S')}] running ...", flush=True)
    res = run_pc_longspan(args.n_periods)

    print(
        f"  -> {res.verdict}: delta1={res.delta1_nd:.3e}, N_req={res.n_periods_required:.1f},"
        f" max drift {res.max_drift_first_half_nd:.3e}/{res.max_drift_second_half_nd:.3e}"
        f" (1st/2nd half, ratio {res.half_ratio:.2f}), dJ={res.jacobi_drift_span:.3e},"
        f" {res.elapsed_s:.1f}s"
    )
    print(
        f"  -> amplitude={res.amplitude_nd:.6f}, band(3A)={res.departure_band_nd:.6f},"
        f" nu={res.nu:.3e}, noise_floor={res.noise_floor_nd:.3e}"
    )
    margin_orders = -math.log10(res.max_drift_span_nd / res.departure_band_nd)
    print(
        f"  -> drift/band = {res.max_drift_span_nd / res.departure_band_nd:.2e}"
        f" ({margin_orders:.1f} orders inside the band)"
    )


if __name__ == "__main__":
    main()
