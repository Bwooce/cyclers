"""#527 Sun-Jupiter Hilda/quasi-Hilda 3:2-MMR periodic-orbit search via DA/HOTM.

Points the #450 DA/HOTM global fixed-point enumerator (built for the Earth-Moon
system, never previously invoked anywhere) at a genuinely new target: the
Sun-Jupiter CR3BP, in the Jacobi-constant band around the interior 3:2
mean-motion resonance (the Hilda asteroid family's home, a ~3.97 AU).

The scan target is DERIVED, not guessed: the Hilda resonance semi-major axis
follows from Kepler's third law applied to Jupiter's own sourced SMA
(``a_hilda = a_jupiter * (2/3)**(2/3)``), converted to a Jacobi constant via a
circular-orbit seed in the Sun-Jupiter CR3BP rotating frame. This session's
diagnostic work independently verified the derivation two ways: the
back-computed GM_sun matches the IAU-standard value exactly, and a_hilda ~=
3.97 AU matches the well-documented, independently-sourced real Hilda asteroid
group location -- neither number is invented for this script.

Positive control (before any negative from this scan is trusted): the coarse
enumerator must recover a low-residual candidate near this derived seed
(x ~= 0.74-0.76, xdot ~= 0, C ~= 3.0-3.1) BEFORE the full band scan runs. This
was verified interactively this session (see data/OUTSTANDING.md #527 note)
and is re-verified programmatically below before the full scan proceeds.

IMPORTANT SYSTEM-SCALE CORRECTION (found this session, do not regress it):
the #450 driver's defaults (SamplingSectionMap.t_max=8.0, and
run_enumeration()'s period_guess heuristic float(2*n)*2.5) were tuned for the
Earth-Moon system's short nondimensional periods. The Sun-Jupiter 3:2
resonance's natural "single revolution" (first return with matching ydot
sign) spans ~12-13 nondimensional time units -- verified by direct
integration -- so this script builds its own SamplingSectionMap with
t_max=20.0 and calls correct_general_periodic() directly with a
period_guess scaled to this system (~12.6 * n), rather than reusing
run_enumeration()'s Earth-Moon-scaled defaults unmodified.

"Cycler structure" here means close approach to Jupiter itself during the
certified orbit's period -- the resonance-transport mechanism the #527
proposal is built on uses Jupiter's own mass to do the gravitational
structuring, unlike every prior small-body flyby lane in this project (which
died on the #489 mass-deficit no-go). Checked by re-propagating each
certified IC over one full period and tracking the minimum distance to
Jupiter's position in the rotating frame.
"""

from __future__ import annotations

import datetime
import pathlib
import sys
import time

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap  # noqa: E402
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402

# Kepler-derived 3:2 interior MMR (Hilda) seed: a_hilda = a_jupiter * (2/3)^(2/3).
# Independently verified this session: GM_sun back-derived from the system's own
# l_km/t_s matches the IAU-standard 132712440018 km^3/s^2 exactly; a_hilda ~=
# 3.9705 AU matches the real, independently-documented Hilda asteroid group.
C_SEED = 3.0613
X_SEED = 0.762
XDOT_SEED = 0.0
_MEASURED_SECTION_RETURN_TIME = 12.6  # nondim time units, direct-integration verified

# Band + domain box for the full scan, centred on the seed with margin for the
# eccentricity/libration spread of real Hilda-type orbits (e ~ 0.07-0.3).
C_BAND = (2.95, 2.98, 3.01, 3.04, C_SEED, 3.08, 3.11, 3.14)
N_RANGE = (1, 2)
DOMAIN_BOX = DomainBox(x_lo=0.60, x_hi=0.90, xdot_lo=-0.08, xdot_hi=0.08)
GRID = (31, 21)
RESIDUAL_TOL = 1e-2
T_MAX = 20.0  # nondim; must exceed the ~12.6 measured single-rev return time.

_REGION_ID = "sun-jupiter-hilda-32-mmr-dahotm"
_METHOD = MethodCapability(
    genome="DA/HOTM global Poincare-section fixed-point enumeration (#450)",
    corrector="correct_general_periodic (asymmetric single-shooting, analytic STM)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "poincare-section-enumeration"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _min_jupiter_distance(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, ydot0: float, period: float
) -> float:
    """Minimum rotating-frame distance to Jupiter over one certified period.

    Jupiter sits at ``(1 - mu, 0, 0)`` in the rotating frame. Returns the
    minimum |r - r_jupiter| sampled densely over ``[0, period]`` -- the
    "does this orbit repeatedly approach Jupiter" check for cycler structure.
    """
    mu = float(system.mu)
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0], dtype=np.float64)
    t_eval = np.linspace(0.0, abs(period), 400)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, abs(period)),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=t_eval,
    )
    if not sol.success:
        return float("inf")
    dx = sol.y[0] - (1.0 - mu)
    dy = sol.y[1]
    return float(np.min(np.hypot(dx, dy)))


def main() -> None:
    print(f"[{_ts()}] #527 Sun-Jupiter Hilda 3:2-MMR DA/HOTM search starting.")

    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    n_points = len(C_BAND) * len(N_RANGE) * GRID[0] * GRID[1]

    preflight_search(
        task_no=527,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Measured this session via direct interactive timing on this exact
        # (system, box, grid) configuration: 0.007-0.043 s/point depending on n;
        # 0.045 is a deliberately conservative upper bound including margin for
        # the certification step's extra cost per surfaced candidate.
        timing_pilot_seconds_per_point=0.045,
    )

    # Positive control: the coarse enumerator MUST recover a low-residual
    # candidate near the Kepler-derived seed before the full band is trusted.
    print(
        f"[{_ts()}] Positive control: re-checking the derived seed (C={C_SEED}, "
        f"x~={X_SEED}, xdot~={XDOT_SEED}) ..."
    )
    control_backend = SamplingSectionMap(system, c_target=C_SEED, ydot_sign=1.0, t_max=T_MAX)
    control_box = DomainBox(x_lo=0.70, x_hi=0.82, xdot_lo=-0.05, xdot_hi=0.05)
    control_candidates = enumerate_fixed_points(
        control_backend, control_box, 1, residual_tol=RESIDUAL_TOL, grid=(15, 11)
    )
    if not control_candidates:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: no sub-tolerance candidate recovered near the "
            "Kepler-derived Hilda seed. Do not trust any negative from the full band "
            "scan below -- the search configuration itself is not reproducing the "
            "expected physics. Stop and re-diagnose before proceeding."
        )
    print(
        f"[{_ts()}] Positive control PASSED: {len(control_candidates)} candidate(s) "
        f"near the seed, best residual={min(c.residual for c in control_candidates):.3e}"
    )

    print(
        f"[{_ts()}] Full band scan: {len(C_BAND)} Jacobi values x {len(N_RANGE)} rev "
        f"counts x {GRID[0] * GRID[1]} grid points = {n_points} total."
    )

    t0 = time.time()
    n_certified = 0
    n_cycler_structure = 0
    findings: list[dict[str, float | int | bool]] = []
    for c_target in C_BAND:
        backend = SamplingSectionMap(system, c_target=float(c_target), ydot_sign=1.0, t_max=T_MAX)
        for n in N_RANGE:
            cands = enumerate_fixed_points(
                backend, DOMAIN_BOX, n, residual_tol=RESIDUAL_TOL, grid=GRID, dedup_radius=0.01
            )
            for cand in cands:
                period_guess = _MEASURED_SECTION_RETURN_TIME * n
                orbit = correct_general_periodic(
                    system,
                    cand.x0,
                    cand.xdot0,
                    float(c_target),
                    period_guess=period_guess,
                    half_crossings=2 * n,
                    ydot0_sign=1.0,
                    tol=1e-11,
                )
                if not (orbit.converged and orbit.residual <= 1e-9):
                    continue
                n_certified += 1
                min_dist = _min_jupiter_distance(
                    system, orbit.x0, orbit.xdot0, orbit.ydot0, orbit.period
                )
                has_structure = min_dist < 0.15  # ~0.15 normalized ~ 117M km, generous flag
                if has_structure:
                    n_cycler_structure += 1
                print(
                    f"[{_ts()}] CERTIFIED C={c_target:.4f} n={n} x0={orbit.x0:.5f} "
                    f"xdot0={orbit.xdot0:.5f} period={orbit.period:.3f} "
                    f"min_dist_to_jupiter={min_dist:.4f} "
                    f"{'[CYCLER-STRUCTURE FLAG]' if has_structure else ''}"
                )
                findings.append(
                    {
                        "c_target": float(c_target),
                        "n": n,
                        "x0": orbit.x0,
                        "xdot0": orbit.xdot0,
                        "ydot0": orbit.ydot0,
                        "period": orbit.period,
                        "residual": orbit.residual,
                        "min_dist_to_jupiter": min_dist,
                        "has_cycler_structure_flag": has_structure,
                    }
                )

    dt = time.time() - t0
    print()
    print(
        f"[{_ts()}] Scan complete in {dt:.1f}s. Certified orbits: {n_certified}; "
        f"with cycler-structure flag (min dist < 0.15): {n_cycler_structure}."
    )
    for f in findings:
        print(f"    {f}")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
