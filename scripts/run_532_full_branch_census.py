"""#532 completion: census ALL distinct periodic-orbit branches in the
2:1/3:2 shared x-range (not just the two Kepler-seeded branches the prior
continuation script traced), then trace each one via continuation and check
every pairwise combination for a genuinely distinct shared-Jacobi unstable
overlap.

`run_532_resonance_family_continuation.py` traced exactly ONE branch per
resonance family -- the one continuously connected to each family's own
Kepler-derived circular-orbit seed -- and found no overlap. But the ORIGINAL
box-enumeration scan (`run_532_resonance_overlap_pilot.py`) surfaced several
OTHER distinct x0 candidates in the same x in [0.5, 0.9] region at a single
C (e.g. x0 around 0.50, 0.62, 0.67, 0.68, 0.76, 0.88 at various C) that were
never traced -- they were the source of the original box-overlap false
positive, but are themselves potentially GENUINE, DISTINCT periodic-orbit
families the #532 hypothesis cares about (any of them could turn out to be
the true "2:1" or "3:2" (or a THIRD, unrelated) resonance family member).

This script:
1. Runs ONE wide DA/HOTM enumeration at a single representative C=3.15
   across the FULL combined x in [0.45, 0.95] range (covering both original
   domain boxes), for n=1 and n=2, to find every distinct fixed point in
   this region -- not per-family, just "every periodic orbit here."
2. Deduplicates by corrected x0 (0.005 clustering radius, well above the
   1e-11 Newton tolerance).
3. Continues EACH distinct branch found in both C-directions (natural
   parameter continuation, warm-started, reusing the #530 monodromy-
   stability classification), covering C=[2.90, 3.25].
4. Checks every pairwise combination of branches for a shared-C point where
   BOTH are unstable AND their x0 values are genuinely distinct (>1e-3
   margin) -- the same distinctness check the earlier fix introduced, now
   applied exhaustively instead of just to the two seeded branches.
"""

from __future__ import annotations

import dataclasses
import datetime
import itertools
import pathlib
import sys
import time

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO))

import numpy as np  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap  # noqa: E402
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402
from scripts.run_532_resonance_overlap_pilot import (  # noqa: E402
    UNSTABLE_THRESHOLD,
    _classify_stability,
)

CENSUS_C = 3.15  # representative C: inside both families' original unstable-looking region
CENSUS_T_MAX = 20.0  # covers both Hilda (~12.6) and interior (~6.26) single-rev return times
CENSUS_DOMAIN = DomainBox(x_lo=0.45, x_hi=0.95, xdot_lo=-0.08, xdot_hi=0.08)
CENSUS_GRID = (41, 21)
CENSUS_N_RANGE = (1, 2)
RESIDUAL_TOL = 1e-2
DEDUP_X0_RADIUS = 0.005
DISTINCTNESS_MARGIN = 1e-3
C_LO, C_HI, C_STEP = 2.90, 3.25, 0.002

_METHOD = MethodCapability(
    genome="Single wide-box DA/HOTM census at C=3.15 across the combined 2:1/3:2 x-range, "
    "then natural-parameter continuation of every distinct branch found",
    corrector="correct_general_periodic (warm-started continuation) + monodromy stability "
    "classification (#530's method)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "natural-parameter-continuation"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


@dataclasses.dataclass(frozen=True)
class TracePoint:
    x0: float
    xdot0: float
    period: float
    lead_eig_mag: float | None
    unstable: bool


def _census(system: cr3bp.CR3BPSystem) -> list[dict]:
    """Find every distinct certified orbit in CENSUS_DOMAIN at CENSUS_C."""
    backend = SamplingSectionMap(system, c_target=CENSUS_C, ydot_sign=1.0, t_max=CENSUS_T_MAX)
    seeds = []
    for n in CENSUS_N_RANGE:
        cands = enumerate_fixed_points(
            backend,
            CENSUS_DOMAIN,
            n,
            residual_tol=RESIDUAL_TOL,
            grid=CENSUS_GRID,
            dedup_radius=0.01,
        )
        for cand in cands:
            orbit = correct_general_periodic(
                system,
                cand.x0,
                cand.xdot0,
                CENSUS_C,
                period_guess=10.0 * n,
                half_crossings=2 * n,
                ydot0_sign=1.0,
                tol=1e-11,
            )
            if orbit.converged and orbit.residual <= 1e-9:
                seeds.append(
                    {
                        "x0": orbit.x0,
                        "xdot0": orbit.xdot0,
                        "ydot0": orbit.ydot0,
                        "period": orbit.period,
                        "n": n,
                    }
                )

    # Deduplicate by x0 (independent of which n found it -- n=1 and n=2 finding the
    # same x0 means the n=2 "orbit" is just 2 reps of the n=1 one).
    distinct: list[dict] = []
    for s in sorted(seeds, key=lambda d: d["x0"]):
        if not any(abs(s["x0"] - d["x0"]) < DEDUP_X0_RADIUS for d in distinct):
            distinct.append(s)
    return distinct


def _trace_branch(
    system: cr3bp.CR3BPSystem,
    x0_seed: float,
    xdot0_seed: float,
    ydot0_seed: float,
    period_seed: float,
) -> dict[float, TracePoint]:
    points: dict[float, TracePoint] = {}
    c_grid = sorted(np.arange(C_LO, C_HI + C_STEP / 2, C_STEP))
    seed_idx = round((CENSUS_C - C_LO) / C_STEP)

    def _walk(indices, x0_start, xdot0_start, period_start) -> None:
        x0_cur, xdot0_cur, period_cur = x0_start, xdot0_start, period_start
        for idx in indices:
            c_target = float(c_grid[idx])
            orbit = correct_general_periodic(
                system,
                x0_cur,
                xdot0_cur,
                c_target,
                period_guess=period_cur,
                half_crossings=2,
                ydot0_sign=1.0,
                tol=1e-11,
            )
            if not (orbit.converged and orbit.residual <= 1e-9):
                break
            state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0], dtype=np.float64)
            lead_mag = _classify_stability(system, state0, orbit.period)
            is_unstable = lead_mag is not None and lead_mag > UNSTABLE_THRESHOLD
            points[round(c_target, 6)] = TracePoint(
                x0=orbit.x0,
                xdot0=orbit.xdot0,
                period=orbit.period,
                lead_eig_mag=lead_mag,
                unstable=is_unstable,
            )
            x0_cur, xdot0_cur, period_cur = orbit.x0, orbit.xdot0, orbit.period

    state0 = np.array([x0_seed, 0.0, 0.0, xdot0_seed, ydot0_seed, 0.0], dtype=np.float64)
    lead_mag0 = _classify_stability(system, state0, period_seed)
    points[round(CENSUS_C, 6)] = TracePoint(
        x0=x0_seed,
        xdot0=xdot0_seed,
        period=period_seed,
        lead_eig_mag=lead_mag0,
        unstable=lead_mag0 is not None and lead_mag0 > UNSTABLE_THRESHOLD,
    )
    _walk(range(seed_idx + 1, len(c_grid)), x0_seed, xdot0_seed, period_seed)
    _walk(range(seed_idx - 1, -1, -1), x0_seed, xdot0_seed, period_seed)
    return points


def main() -> None:
    print(f"[{_ts()}] #532 full-branch census starting.")
    system = cr3bp.cr3bp_system("Sun", "Jupiter")

    n_census_points = CENSUS_GRID[0] * CENSUS_GRID[1] * len(CENSUS_N_RANGE)
    n_continuation_points_est = 6 * int((C_HI - C_LO) / C_STEP)  # ~6 branches expected
    preflight_search(
        task_no=532,
        region_id="sun-jupiter-21-32-mmr-full-branch-census",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_census_points + n_continuation_points_est,
        timing_pilot_seconds_per_point=0.15,
    )

    t0 = time.time()
    branches = _census(system)
    print(f"[{_ts()}] Census at C={CENSUS_C}: {len(branches)} distinct branch(es) found:")
    for b in branches:
        print(f"[{_ts()}]     x0={b['x0']:.6f} (found via n={b['n']}), period={b['period']:.4f}")

    traces: dict[float, dict[float, TracePoint]] = {}
    for b in branches:
        trace = _trace_branch(system, b["x0"], b["xdot0"], b["ydot0"], b["period"])
        traces[b["x0"]] = trace
        n_unstable = sum(1 for p in trace.values() if p.unstable)
        c_lo_actual, c_hi_actual = min(trace), max(trace)
        print(
            f"[{_ts()}] branch x0_seed={b['x0']:.6f}: traced {len(trace)} points over "
            f"C=[{c_lo_actual:.4f},{c_hi_actual:.4f}], {n_unstable} unstable."
        )

    dt = time.time() - t0
    print(f"[{_ts()}] Census + continuation complete in {dt:.1f}s.")

    print()
    print(f"[{_ts()}] Pairwise cross-branch shared-Jacobi unstable check:")
    genuine_pairs = []
    for (x0_a, trace_a), (x0_b, trace_b) in itertools.combinations(traces.items(), 2):
        shared_c = sorted(set(trace_a) & set(trace_b))
        for c in shared_c:
            pa, pb = trace_a[c], trace_b[c]
            if pa.unstable and pb.unstable:
                distinct = abs(pa.x0 - pb.x0) > DISTINCTNESS_MARGIN
                print(
                    f"[{_ts()}]   C={c:.4f}: branch[{x0_a:.4f}] x0={pa.x0:.6f} "
                    f"(|lam|={pa.lead_eig_mag:.4f}) vs branch[{x0_b:.4f}] x0={pb.x0:.6f} "
                    f"(|lam|={pb.lead_eig_mag:.4f}) -- {'DISTINCT' if distinct else 'SAME ORBIT'}"
                )
                if distinct:
                    genuine_pairs.append((c, x0_a, x0_b))

    print()
    if genuine_pairs:
        print(
            f"[{_ts()}] VERDICT: GO -- {len(genuine_pairs)} genuinely distinct shared-Jacobi "
            f"unstable pair(s) found across the full branch census: {genuine_pairs}"
        )
    else:
        print(
            f"[{_ts()}] VERDICT: NO-GO -- no genuinely distinct shared-Jacobi unstable pair "
            f"found among any of the {len(branches)} censused branches over "
            f"C=[{C_LO},{C_HI}]. This is now a much more exhaustive negative than the "
            f"two-seed-only continuation."
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
