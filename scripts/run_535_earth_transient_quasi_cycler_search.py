"""#535 transient-drift-phase quasi_cycler search (Earth co-orbital region).

Applies the admission criterion settled in writing
(``docs/notes/2026-07-03-535-quasi-cycler-transient-drift-admission-
criterion.md``) to a scan of Earth co-orbital initial conditions, looking for
a genuinely BOUNDED-WITHIN-A-WINDOW, repeated-Hill-sphere-encounter
trajectory -- the transient/chaotic-drift phenomenon #523's own strictly-
periodic-orbit search could never find (#523's 120 certified periodic
orbits all stayed 4.9-10.3x outside the Hill radius; see #523's
OUTSTANDING.md entry).

NEGATIVE CONTROL (a REAL finding from this session, kept as a documented
data point, not silently dropped): the literal 2006 RH120-derived seed used
as #523's positive control (x0=0.97881371, xdot0=0.0, C=2.9998797409719242)
was independently re-verified this session to be a SINGLE-transient
encounter (one close approach then departure to ~2 AU, no return within at
least 60 years) -- it does NOT itself satisfy this task's admission
criterion (needs >= 3 returns). This matches the REAL astronomical 2006
RH120's own documented ~1-year-only capture (Jul 2006-Jul 2007); an earlier
claim that this seed also showed "recurring approaches every 4.6-9.2 years"
was an error, corrected in #523's own script docstring the same session
this script was built.

POSITIVE FINDING (this session, interactively verified before this script
was written, robust across 3 independent integration tolerance/sampling
settings to 4 decimal places): x0=0.9920, xdot0=0.0, ydot0=0.13033911
(Jacobi C=2.9839437412) shows 3 distinct Hill-sphere returns at
t=0, 2.99, 6.02 years, closest approaches 0.48-0.80x the Hill radius, ALL
within a single 10-year admission window with geometry ratio 1.65 (well
under the 3.0 gate) -- a genuine ADMISSIBLE candidate under this task's own
criterion. This script's scan reproduces that finding as its own positive
control and extends the search to nearby (x0, C) combinations.
"""

from __future__ import annotations

import datetime
import pathlib
import sys
import time

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.search.hill_sphere_return_detector import (  # noqa: E402
    find_admission_windows,
    find_returns,
)

# The verified positive-control candidate (see module docstring).
CONTROL_X0 = 0.9920
CONTROL_XDOT0 = 0.0
CONTROL_C = 2.9839437412

# Scan range: x0 around the verified candidate, at a handful of nearby C
# values (the verified C plus the original #523 positive-control C, for
# context/comparison).
X0_GRID = np.linspace(0.960, 0.999, 40)
C_VALUES = (CONTROL_C, 2.9998797409719242)
YEARS = 50
MIN_SEPARATION_YEARS = 1.0
WINDOW_LO_YEARS, WINDOW_HI_YEARS = 10.0, 15.0
N_RETURNS_LO, N_RETURNS_HI = 3, 15
GEOMETRY_FACTOR = 3.0

_REGION_ID = "sun-earth-coorbital-transient-drift-quasi-cycler"
_METHOD = MethodCapability(
    genome="Direct long-duration propagation + repeated-Hill-sphere-return detection "
    "(search/hill_sphere_return_detector.py) -- targets aperiodic/chaotic transient-drift "
    "trajectories, structurally distinct from #523's strictly-periodic-orbit search",
    corrector="none (no periodicity requirement by design -- see #535's admission criterion)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "hill-sphere-return-detection"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _ydot0_from_jacobi(
    x0: float, xdot0: float, c_target: float, mu: float, sign: float
) -> float | None:
    r1 = abs(x0 + mu)
    r2 = abs(x0 - 1.0 + mu)
    rad = x0 * x0 + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - c_target - xdot0 * xdot0
    if rad < 0.0:
        return None
    return float(sign) * float(np.sqrt(rad))


def _propagate_and_detect(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, c_target: float, r_hill: float
) -> tuple[int, list] | None:
    mu = float(system.mu)
    ydot0 = _ydot0_from_jacobi(x0, xdot0, c_target, mu, 1.0)
    if ydot0 is None:
        return None
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])
    t_max = YEARS * 2.0 * np.pi
    t_eval = np.linspace(0.0, t_max, YEARS * 800)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_max),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=t_eval,
        max_step=0.5,
    )
    if not sol.success:
        return None
    pos_rel_earth = np.stack([sol.y[0] - (1.0 - mu), sol.y[1]], axis=1)
    returns = find_returns(
        sol.t, pos_rel_earth, r_hill=r_hill, min_separation=MIN_SEPARATION_YEARS * 2.0 * np.pi
    )
    windows = find_admission_windows(
        returns,
        float(sol.t[0]),
        float(sol.t[-1]),
        window_lo=WINDOW_LO_YEARS * 2.0 * np.pi,
        window_hi=WINDOW_HI_YEARS * 2.0 * np.pi,
        n_returns_lo=N_RETURNS_LO,
        n_returns_hi=N_RETURNS_HI,
        geometry_factor=GEOMETRY_FACTOR,
    )
    admissible = [w for w in windows if w.geometry_ok]
    return len(returns), admissible


def main() -> None:
    print(f"[{_ts()}] #535 Earth co-orbital transient-drift quasi_cycler search starting.")

    system = cr3bp.cr3bp_system("Sun", "Earth")
    r_hill = (float(system.mu) / 3.0) ** (1.0 / 3.0)

    n_points = len(X0_GRID) * len(C_VALUES)
    preflight_search(
        task_no=535,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Measured this session via direct interactive timing of the same
        # propagate+detect call on this exact grid: comfortably under 1s/point.
        timing_pilot_seconds_per_point=1.0,
    )

    print(
        f"[{_ts()}] Positive control: verifying the known admissible candidate "
        f"(x0={CONTROL_X0}, C={CONTROL_C}) ..."
    )
    control_result = _propagate_and_detect(system, CONTROL_X0, CONTROL_XDOT0, CONTROL_C, r_hill)
    if control_result is None or not control_result[1]:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: the known admissible candidate did not reproduce. "
            "Do not trust any result from the scan below."
        )
    n_ctrl, windows_ctrl = control_result
    print(
        f"[{_ts()}] Positive control PASSED: {n_ctrl} returns, "
        f"{len(windows_ctrl)} admissible window(s), best geometry_ratio="
        f"{min(w.geometry_ratio for w in windows_ctrl):.3f}"
    )

    t0 = time.time()
    findings = []
    for c_target in C_VALUES:
        for x0 in X0_GRID:
            result = _propagate_and_detect(system, float(x0), 0.0, float(c_target), r_hill)
            if result is None:
                continue
            n_returns, admissible = result
            if admissible:
                findings.append(
                    {
                        "x0": float(x0),
                        "c_target": float(c_target),
                        "n_returns_total": n_returns,
                        "n_admissible_windows": len(admissible),
                        "best_geometry_ratio": min(w.geometry_ratio for w in admissible),
                    }
                )
                print(
                    f"[{_ts()}]   ADMISSIBLE: x0={x0:.4f} C={c_target:.7f} "
                    f"n_returns={n_returns} n_windows={len(admissible)}"
                )

    dt = time.time() - t0
    print()
    print(
        f"[{_ts()}] Scan complete in {dt:.1f}s over {n_points} points. "
        f"Admissible candidates: {len(findings)}."
    )
    for f in findings:
        print(f"    {f}")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
