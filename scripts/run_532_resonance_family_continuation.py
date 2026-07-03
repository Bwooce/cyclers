"""#532 fix: replace box-enumeration "family" identification with genuine
natural-parameter continuation from each family's own seed.

`run_532_resonance_overlap_pilot.py` / `_finegrid.py` identified "family
membership" only by which overlapping DomainBox an independent per-C DA/HOTM
enumeration + Newton correction happened to converge into. That is NOT a
family-membership test: at C=3.14-3.16 BOTH the "3:2 Hilda" and "2:1
interior" scans converged to the IDENTICAL x0 (to 5 decimal places) with
IDENTICAL eigenvalues -- since (x0, xdot0=0, C) uniquely fixes the full
state via ydot0, an identical x0 at the same C is mathematically the SAME
orbit, not two distinct family members that happen to share energy. The
prior "GO" verdict was a false positive: the two overlapping search boxes
(Hilda x in [0.61,0.91], interior x in [0.50,0.76]) both converged onto one
single periodic-orbit branch passing through their shared overlap region.

The only rigorous definition of "family membership" is a connected
component of the solution manifold: a point belongs to family F iff it is
reachable from F's own defining seed by a continuous path of converged
periodic orbits. This script traces each family via NATURAL-PARAMETER
CONTINUATION IN C, starting from its own Kepler-derived seed (already
verified in the coarse pilot to reproduce #527's own committed C_SEED to 4
decimals), stepping C and re-using the previous step's converged (x0, xdot0,
period) as the next step's initial guess -- the same warm-start pattern
#534's L2 seed continuation used. This guarantees every point genuinely
traces back to its own family's seed by construction, unlike independent
per-C grid re-enumeration.

At any C where BOTH continuations are simultaneously unstable, this script
explicitly checks x0_hilda(C) != x0_interior(C) by a margin far above Newton
tolerance before calling it a genuine two-orbit candidate pair -- the exact
distinctness check the prior scripts skipped.
"""

from __future__ import annotations

import dataclasses
import datetime
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
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402
from scripts.run_532_resonance_overlap_pilot import (  # noqa: E402
    UNSTABLE_THRESHOLD,
    _classify_stability,
    _seed_from_period_ratio,
)

DISTINCTNESS_MARGIN = 1e-3  # far above the 1e-11 Newton tolerance
C_LO, C_HI, C_STEP = 3.00, 3.20, 0.002

_METHOD = MethodCapability(
    genome="Natural-parameter continuation in C from each family's own Kepler-derived seed "
    "(fixes the #532 box-enumeration family-identification bug)",
    corrector="correct_general_periodic, warm-started from the previous C step's converged "
    "solution + monodromy stability classification (#530's method)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "natural-parameter-continuation"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


@dataclasses.dataclass(frozen=True)
class TracePoint:
    c: float
    x0: float
    xdot0: float
    period: float
    lead_eig_mag: float | None
    unstable: bool


def _trace_family(
    system: cr3bp.CR3BPSystem,
    label: str,
    x0_seed: float,
    xdot0_seed: float,
    c_seed: float,
    period_guess_unit: float,
    n: int,
) -> dict[float, TracePoint]:
    """Continue from (x0_seed, xdot0_seed, c_seed) in both C-directions."""
    points: dict[float, TracePoint] = {}

    def _converge_one(x0: float, xdot0: float, c_target: float, period_guess: float):
        return correct_general_periodic(
            system,
            x0,
            xdot0,
            c_target,
            period_guess=period_guess,
            half_crossings=2 * n,
            ydot0_sign=1.0,
            tol=1e-11,
        )

    # Anchor point (the seed itself, or nearest grid value to it).
    seed_orbit = _converge_one(x0_seed, xdot0_seed, c_seed, period_guess_unit * n)
    if not (seed_orbit.converged and seed_orbit.residual <= 1e-9):
        raise RuntimeError(
            f"POSITIVE CONTROL FAILED: {label} seed itself did not converge "
            f"(residual={seed_orbit.residual:.3e}). Do not trust this family's trace."
        )
    print(
        f"[{_ts()}] {label}: seed converged at C={c_seed:.6f}, x0={seed_orbit.x0:.6f}, "
        f"residual={seed_orbit.residual:.3e}."
    )

    c_grid = sorted({round(c, 6) for c in np.arange(C_LO, C_HI + C_STEP / 2, C_STEP)})
    # Snap the seed's own C onto the grid so continuation starts exactly there.
    c_grid = sorted(set(c_grid) | {round(c_seed, 6)})
    seed_idx = c_grid.index(round(c_seed, 6))

    def _walk(indices) -> None:
        x0_cur, xdot0_cur, period_cur = seed_orbit.x0, seed_orbit.xdot0, seed_orbit.period
        for idx in indices:
            c_target = c_grid[idx]
            orbit = _converge_one(x0_cur, xdot0_cur, c_target, period_cur)
            if not (orbit.converged and orbit.residual <= 1e-9):
                break  # fold / boundary of this family's continuable range
            state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0], dtype=np.float64)
            lead_mag = _classify_stability(system, state0, orbit.period)
            is_unstable = lead_mag is not None and lead_mag > UNSTABLE_THRESHOLD
            points[c_target] = TracePoint(
                c=c_target,
                x0=orbit.x0,
                xdot0=orbit.xdot0,
                period=orbit.period,
                lead_eig_mag=lead_mag,
                unstable=is_unstable,
            )
            x0_cur, xdot0_cur, period_cur = orbit.x0, orbit.xdot0, orbit.period

    points[round(c_seed, 6)] = TracePoint(
        c=round(c_seed, 6),
        x0=seed_orbit.x0,
        xdot0=seed_orbit.xdot0,
        period=seed_orbit.period,
        lead_eig_mag=None,
        unstable=False,
    )
    _walk(range(seed_idx + 1, len(c_grid)))  # upward
    _walk(range(seed_idx - 1, -1, -1))  # downward

    n_unstable = sum(1 for p in points.values() if p.unstable)
    print(
        f"[{_ts()}] {label}: traced {len(points)} points over C in "
        f"[{min(points):.4f}, {max(points):.4f}], {n_unstable} unstable."
    )
    return points


def main() -> None:
    print(f"[{_ts()}] #532 rigorous family-continuation fix starting.")
    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    mu = float(system.mu)

    n_points = 2 * int((C_HI - C_LO) / C_STEP)
    preflight_search(
        task_no=532,
        region_id="sun-jupiter-21-32-mmr-family-continuation-verified-overlap",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Warm-started Newton correction from an adjacent converged point is far cheaper
        # than the box-enumeration pilot's cold-start grid search; conservative estimate.
        timing_pilot_seconds_per_point=0.15,
    )

    x0_h, _, c_seed_h = _seed_from_period_ratio(mu, 2.0 / 3.0)
    x0_i, _, c_seed_i = _seed_from_period_ratio(mu, 1.0 / 2.0)
    if abs(c_seed_h - 3.0613) > 1e-3:
        raise RuntimeError(
            f"POSITIVE CONTROL FAILED: re-derived Hilda seed C={c_seed_h:.6f} does not "
            f"match #527's committed C_SEED=3.0613."
        )
    print(f"[{_ts()}] Positive control PASSED: Hilda seed C={c_seed_h:.6f} matches #527.")

    t0 = time.time()
    hilda_points = _trace_family(system, "3:2 Hilda", x0_h, 0.0, c_seed_h, 12.6, n=1)
    interior_points = _trace_family(system, "2:1 interior", x0_i, 0.0, c_seed_i, 6.26, n=1)
    dt = time.time() - t0

    shared_c = sorted(set(hilda_points) & set(interior_points))
    print()
    print(f"[{_ts()}] Continuation complete in {dt:.1f}s. Shared-C grid points: {len(shared_c)}.")

    genuine_pairs = []
    for c in shared_c:
        hp, ip = hilda_points[c], interior_points[c]
        if hp.unstable and ip.unstable:
            distinct = abs(hp.x0 - ip.x0) > DISTINCTNESS_MARGIN
            print(
                f"[{_ts()}]   C={c:.4f}: Hilda x0={hp.x0:.6f} (|lam|={hp.lead_eig_mag:.4f}) "
                f"vs interior x0={ip.x0:.6f} (|lam|={ip.lead_eig_mag:.4f}) "
                f"-- {'DISTINCT' if distinct else 'SAME ORBIT (bug/overlap)'}"
            )
            if distinct:
                genuine_pairs.append(c)

    print()
    if genuine_pairs:
        print(
            f"[{_ts()}] VERDICT: GO -- {len(genuine_pairs)} Jacobi constant(s) host "
            f"GENUINELY DISTINCT unstable orbits in BOTH families (verified by seed-"
            f"traced continuation + explicit x0 distinctness check): {genuine_pairs}"
        )
    else:
        print(
            f"[{_ts()}] VERDICT: NO-GO -- no genuinely distinct shared-Jacobi unstable "
            f"pair found in C=[{C_LO},{C_HI}] at this continuation step size ({C_STEP}). "
            f"Any earlier apparent overlap was the box-enumeration same-orbit artifact."
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
