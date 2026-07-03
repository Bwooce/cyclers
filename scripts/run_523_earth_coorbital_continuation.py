"""#523 rework: continuation-based Earth co-orbital horseshoe/QS orbit sweep.

Replaces ``run_523_earth_coorbital_search.py``'s brute-force per-candidate
certification (60-100s each via chained Newton passes from a cold, possibly
far coarse-enumerator guess) with pseudo-arclength continuation along the
family curve in ``(x0, xdot0, C)``
(``search/cr3bp_general_periodic_free_c.py``): certify ONE orbit, then walk
outward in both directions, each step a cheap Newton correction from an
already-good tangent-based predictor. Measured this session:
~1.4s/continuation-step vs. 60-100s/candidate for the original approach --
see ``search/cr3bp_general_periodic_free_c.py``'s module docstring for the
full derivation and validation of the extra analytic ``dR/dC`` Jacobian
column this relies on.

Same physical target, positive control, and encounter criterion as the
original script (read its docstring for the full justification): Sun-Earth
co-orbital horseshoe/quasi-satellite periodic orbits, Hill-sphere encounter
criterion, seeded from 2006 RH120's real published elements.
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
from cyclerfinder.search.cr3bp_general_periodic import (  # noqa: E402
    _ydot0_general,
    correct_general_periodic,
)
from cyclerfinder.search.cr3bp_general_periodic_free_c import (  # noqa: E402
    make_residual_and_jacobian_fns,
)
from cyclerfinder.search.pseudo_arclength import (  # noqa: E402
    continue_curve,
)

# Same #523 positive-control seed (2006 RH120-derived elements) and the same
# coarse guess already confirmed (this session) to certify at it.
SEED_C = 2.9998797409719242
SEED_X0_GUESS = 0.93000
SEED_XDOT0_GUESS = 0.03600

# Covers the original script's C_BAND (2.9990-3.0010) with margin.
C_LO, C_HI = 2.9980, 3.0020
STEP_SIZE = 0.0005
MAX_STEPS_PER_DIRECTION = 200  # generous headroom; the walk stops at C_LO/C_HI well before this
STEP_CAPS = [0.02, 0.02, 0.001]
T_HI_FRAC = 1.15
HALF_CROSSINGS = 2

_REGION_ID = "sun-earth-coorbital-horseshoe-qsat-dahotm"
_METHOD = MethodCapability(
    genome="Pseudo-arclength continuation in (x0,xdot0,C) via an analytic STM "
    "Jacobian (search/cr3bp_general_periodic_free_c.py) -- #523 rework, strictly "
    "stronger coverage than the original DA/HOTM grid enumeration at the same region",
    corrector="correct_general_periodic seed + continue_curve (analytic Jacobian)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "pseudo-arclength-continuation"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _min_dist_to_earth(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, ydot0: float, period: float
) -> float:
    mu = float(system.mu)
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0], dtype=np.float64)
    t_eval = np.linspace(0.0, abs(period), 600)
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
    print(f"[{_ts()}] #523 rework: continuation-based Sun-Earth co-orbital sweep starting.")

    system = cr3bp.cr3bp_system("Sun", "Earth")
    r_hill = (float(system.mu) / 3.0) ** (1.0 / 3.0)

    n_walk_points = 2 * MAX_STEPS_PER_DIRECTION
    preflight_search(
        task_no=523,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_walk_points,
        # Measured this session via direct interactive timing of continue_curve
        # on this exact residual/Jacobian pair: ~1.4 s/continuation-step
        # (dominated by two STM propagations over the ~80-nondim-TU horizon).
        timing_pilot_seconds_per_point=1.4,
    )

    print(
        f"[{_ts()}] Positive control: certifying the 2006-RH120-derived seed "
        f"(C={SEED_C:.6f}) via the existing corrector ..."
    )
    seed_orbit = correct_general_periodic(
        system,
        SEED_X0_GUESS,
        SEED_XDOT0_GUESS,
        SEED_C,
        period_guess=100.0,
        half_crossings=HALF_CROSSINGS,
        ydot0_sign=1.0,
        tol=1e-11,
        t_hi_frac=T_HI_FRAC,
        max_iter=60,
    )
    if not (seed_orbit.converged and seed_orbit.residual <= 1e-9):
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: could not re-certify the 2006-RH120-derived seed. "
            "Do not trust any result from the continuation walk below."
        )
    print(
        f"[{_ts()}] Positive control PASSED: x0={seed_orbit.x0:.6f} xdot0={seed_orbit.xdot0:.6f} "
        f"period={seed_orbit.period:.4f} residual={seed_orbit.residual:.3e}"
    )

    t0 = time.time()
    n_certified = 0
    n_encounter = 0
    findings: list[dict[str, float | int | bool]] = []
    mu = float(system.mu)

    for direction, step_size in ((1, STEP_SIZE), (-1, -STEP_SIZE)):
        residual_fn, jacobian_fn = make_residual_and_jacobian_fns(
            mu, 1.0, HALF_CROSSINGS, T_HI_FRAC, seed_orbit.period, rtol=1e-12, atol=1e-12
        )
        z0 = np.array([seed_orbit.x0, seed_orbit.xdot0, SEED_C])
        curve = continue_curve(
            residual_fn,
            z0,
            jacobian_fn=jacobian_fn,
            step_size=step_size,
            max_steps=MAX_STEPS_PER_DIRECTION,
            tol=1e-9,
            max_iter=30,
            step_caps=STEP_CAPS,
            target_index=2,
            target_value=C_HI if direction > 0 else C_LO,
        )
        print(
            f"[{_ts()}] direction={direction:+d}: {len(curve.points)} points, "
            f"stop_reason={curve.stop_reason}"
        )
        for p in curve.points:
            x0_p, xdot0_p, c_p = float(p.z[0]), float(p.z[1]), float(p.z[2])
            try:
                _ydot0_general(x0_p, xdot0_p, c_p, mu, 1.0)
            except ValueError:
                continue  # infeasible point (shouldn't happen post-continuation, but be safe)
            # continue_curve's own points carry (x0,xdot0,C) but not period/ydot0/
            # residual in the GeneralPeriodicOrbit shape the encounter check
            # needs -- a cheap few-iteration polish at this already-converged
            # point recovers those fields (this is NOT a fresh cold-start
            # certification; it should converge in 1-2 iterations).
            polished = correct_general_periodic(
                system,
                x0_p,
                xdot0_p,
                c_p,
                period_guess=seed_orbit.period,
                half_crossings=HALF_CROSSINGS,
                ydot0_sign=1.0,
                tol=1e-9,
                t_hi_frac=T_HI_FRAC,
                max_iter=5,
            )
            if not polished.converged:
                continue
            n_certified += 1
            min_dist = _min_dist_to_earth(
                system, polished.x0, polished.xdot0, polished.ydot0, polished.period
            )
            is_encounter = min_dist < r_hill
            if is_encounter:
                n_encounter += 1
            findings.append(
                {
                    "c_target": c_p,
                    "x0": polished.x0,
                    "xdot0": polished.xdot0,
                    "ydot0": polished.ydot0,
                    "period": polished.period,
                    "residual": polished.residual,
                    "min_dist_to_earth": min_dist,
                    "is_encounter": is_encounter,
                    "direction": direction,
                }
            )

    dt = time.time() - t0
    print()
    print(
        f"[{_ts()}] Continuation sweep complete in {dt:.1f}s. Certified points: {n_certified}; "
        f"meeting the Hill-sphere encounter criterion: {n_encounter}."
    )
    c_values = sorted({round(f["c_target"], 4) for f in findings})
    print(
        f"    C coverage: {len(c_values)} distinct values from "
        f"{min(c_values):.4f} to {max(c_values):.4f}"
    )
    for f in findings:
        print(f"    {f}")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
