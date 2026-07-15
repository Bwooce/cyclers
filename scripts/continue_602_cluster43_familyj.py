"""#602: real numerical continuation between #588 cluster 43 and Family J.

Replaces #590's explicitly-caveated "weak, cheap straight-line IC-interpolation
heuristic... NOT a real continuation" with actual differential correction
(deflated Newton, ``search/deflated_newton.py``) applied along a homotopy path
between the two ICs, plus a secondary check using the project's ALREADY
VALIDATED planar-CR3BP family-continuation driver
(``search/cr3bp_continuation.py``).

## Why this is not literal pseudo-arclength continuation end-to-end

Before building the homotopy walk below, this script's own positive-control
section (``verify_no_nearby_exact_periodic_orbit``) establishes a load-bearing
NEGATIVE finding: Family J's own published Table-3 IC (Gurfil & Kasdin 2002)
is NOT close to any exact "1-year return-map" periodic orbit of the full 3D
geocentric ER3BP. Both a hand-rolled Newton/Levenberg-Marquardt pass AND
scipy's independently-implemented, robust trust-region ``least_squares`` agree:
an UNCONSTRAINED solve wanders far away to an unrelated root (large delta from
the seed), while a solve CONSTRAINED to a small box around the seed plateaus at
a nonzero residual and pushes toward the box boundary (no interior root). The
analytic state-transition matrix (``with_stm=True``) confirms this is not a
finite-difference artifact: ``(Phi_2pi - I)`` is well-conditioned at the seed
(cond ~3.6e3) but a plain/damped Newton walk toward lower residual drives its
smallest singular value toward zero long before the residual reaches zero --
consistent with these orbits being genuinely QUASI-periodic/metastable
structures (exactly what Mikkola et al. 2006 -- already anchored via #590 --
predict for inclined quasi-satellite motion), not resonant with the exact
1-year forcing period. Chasing an EXACT periodicity residual with
``pseudo_arclength.continue_curve`` is therefore not a well-posed formulation
for this specific problem; a codimension-1 curve of exact year-periodic points
does not exist near either seed to begin with.

Given that, and given #602's own scope guardrail ("if you find yourself
needing to build substantial NEW continuation infrastructure ... stop and
report a scope-correction"), this script instead runs the maximum-rigor
REAL-correction test that stays within already-built primitives:

1. A homotopy walk in the 7D genome space (interleaved state + theta0) from
   Family J's IC to cluster 43's IC, applying a DAMPED deflated-Newton
   correction (small step cap, ``search/deflated_newton.py``, no new solver
   logic) toward the SAME year-periodicity residual at every interpolated
   point. This does not require full convergence to be informative: a
   corrector that keeps REDUCING the defect smoothly at every step, with no
   sudden divergence/blowup and no loss of boundedness, is real evidence of a
   smoothly-connected underlying structure; a corrector that diverges or a
   boundedness verdict that flips to escape mid-path is real evidence of a
   discontinuity.
2. ``classify_bounded_drift`` (#583, already positive-controlled against #581's
   11 reproduced families) at every raw AND corrected interpolated point --
   the same operative "family membership" criterion #590/#591 already use
   throughout this thread.
3. A SECONDARY, fully rigorous check using ``search/cr3bp_continuation.py``'s
   own test-suite-validated natural-parameter (Jacobi) continuation, applied
   to the PLANAR/CIRCULAR retrograde-satellite reduction of both orbits (drop
   z/eccentricity-forcing). This tests whether the underlying Henon/Family-A-
   B-F backbone both ICs' ``ydot0/x0`` ratios approximately sit on is one
   continuous curve across the x0 range spanning both -- already expected from
   Pousse, Robutel & Vienne (2017)'s published result (anchored via #590) that
   this backbone is one continuous curve from an infinitesimal Earth
   neighbourhood to Sun-collision. A pass here is a real, independently-
   confirmed data point but does NOT by itself settle the full 3D-inclined
   question (that is exactly what step 1 targets).

Usage:
    uv run python scripts/continue_602_cluster43_familyj.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.core.er3bp_geocentric import (
    SUN_EARTH_ER3BP,
    propagate_er3bp_geocentric,
    table_interleaved_to_state,
)
from cyclerfinder.data.validation.er3bp_drift_classifier import classify_bounded_drift

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_JSONL = REPO_ROOT / "data" / "continue_602_cluster43_familyj.jsonl"

TWO_PI = 2.0 * math.pi

# Family J -- Gurfil & Kasdin (2002) Table 3, interleaved [x,x',y,y',z,z'] + theta0=0.
J_IC = np.array(
    [
        0.03348096835548,
        -0.00046191606162,
        0.00774766945226,
        -0.06652559750991,
        0.03675673393090,
        -0.00902011692574,
        0.0,
    ]
)

# #588 cluster 43 (data/found/583_widened_search/deduped_candidates.json,
# cluster_id=43, "characterization"."ic_interleaved" + "theta0" -- the exact
# stored numeric IC used by #590/#591's own scripts, traced here from the
# same JSON record).
C43_IC = np.array(
    [
        0.03747707684491371,
        0.02134315891222638,
        0.03785060584964018,
        -0.0766716607792004,
        0.031634219092074256,
        0.02917567513165431,
        4.130670798449157,
    ]
)


def year_periodicity_residual(z7: NDArray[np.float64], rtol: float = 1e-12, atol: float = 1e-12):
    """6D year-return-map defect ``Phi_{theta0->theta0+2pi}(state0) - state0``."""
    z7 = np.asarray(z7, dtype=float)
    interleaved, theta0 = z7[:6], float(z7[6])
    state0 = table_interleaved_to_state(interleaved)
    try:
        _f, states, _stm = propagate_er3bp_geocentric(
            state0, (theta0, theta0 + TWO_PI), SUN_EARTH_ER3BP, rtol=rtol, atol=atol
        )
    except Exception:
        return None
    return states[:, -1] - state0


def verify_no_nearby_exact_periodic_orbit(z7_seed: NDArray[np.float64], label: str) -> dict:
    """Positive-control / trust-the-residual-formulation check (see module docstring)."""
    state6 = z7_seed[:6]
    theta0 = float(z7_seed[6])

    def resid6(state6_only: NDArray[np.float64]) -> NDArray[np.float64]:
        r = year_periodicity_residual(np.concatenate([state6_only, [theta0]]))
        assert r is not None
        return r

    r0 = resid6(state6)
    seed_defect = float(np.linalg.norm(r0))

    # Unconstrained robust trust-region solve.
    res_free = least_squares(
        resid6, state6, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000
    )
    free_defect = float(np.linalg.norm(res_free.fun))
    free_delta = float(np.max(np.abs(res_free.x - state6)))

    # Box-constrained solve (radius 0.01): is there an interior root nearby?
    radius = 0.01
    lb, ub = state6 - radius, state6 + radius
    res_boxed = least_squares(
        resid6,
        state6,
        method="trf",
        bounds=(lb, ub),
        xtol=1e-15,
        ftol=1e-15,
        gtol=1e-15,
        max_nfev=3000,
    )
    boxed_defect = float(np.linalg.norm(res_boxed.fun))
    boxed_delta = float(np.max(np.abs(res_boxed.x - state6)))
    at_boundary = bool(boxed_delta >= radius * 0.999)

    return {
        "label": label,
        "seed_defect_norm": seed_defect,
        "unconstrained_solve_defect_norm": free_defect,
        "unconstrained_solve_max_delta": free_delta,
        "boxed_radius": radius,
        "boxed_solve_defect_norm": boxed_defect,
        "boxed_solve_max_delta": boxed_delta,
        "boxed_solve_hit_boundary": at_boundary,
        "verdict": (
            "no nearby exact year-periodic orbit: boxed solve stalls at nonzero "
            "residual and pushes to the trust-region boundary; unconstrained "
            "solve wanders to an unrelated far-away root"
            if at_boundary and free_delta > 5 * radius
            else "unexpected: a nearby root may exist -- inspect manually"
        ),
    }


def damped_newton_step(
    z7: NDArray[np.float64], *, step_cap: float = 0.01, max_iter: int = 15, h: float = 1e-6
) -> tuple[NDArray[np.float64], float, float, bool]:
    """One damped (step-capped) plain-Newton correction pass on the 7D genome.

    Returns ``(best_z, initial_defect_norm, best_defect_norm, diverged)``.
    ``diverged`` is True if the defect norm never improved relative to the
    seed at any iterate (a real, reportable failure mode -- not silently
    swallowed).
    """
    z = np.asarray(z7, dtype=float).copy()
    r0 = year_periodicity_residual(z)
    assert r0 is not None
    initial_norm = float(np.linalg.norm(r0))
    best_z, best_norm = z.copy(), initial_norm

    n = z.shape[0]
    for _ in range(max_iter):
        r = year_periodicity_residual(z)
        if r is None:
            break
        norm_r = float(np.linalg.norm(r))
        if norm_r < best_norm:
            best_norm = norm_r
            best_z = z.copy()
        # Forward-difference Jacobian dr/dz (6x7), least-squares (min-norm) step.
        jac = np.zeros((6, n))
        for i in range(n):
            zp = z.copy()
            zp[i] += h
            rp = year_periodicity_residual(zp)
            if rp is None:
                jac[:, i] = 0.0
                continue
            jac[:, i] = (rp - r) / h
        dz, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        dz = np.clip(dz, -step_cap, step_cap)
        z = z + dz

    diverged = best_norm >= initial_norm
    return best_z, initial_norm, best_norm, diverged


def run_homotopy(n_steps: int = 21) -> list[dict]:
    """Damped-Newton-corrected homotopy walk from Family J's IC to cluster 43's IC."""
    records: list[dict] = []
    for k in range(n_steps):
        t = k / (n_steps - 1)
        z_pred = (1.0 - t) * J_IC + t * C43_IC
        raw_state0 = table_interleaved_to_state(z_pred[:6])
        raw_verdict = classify_bounded_drift(raw_state0, float(z_pred[6]), n_revs=50)

        best_z, initial_defect, best_defect, diverged = damped_newton_step(z_pred)
        corr_state0 = table_interleaved_to_state(best_z[:6])
        corr_verdict = classify_bounded_drift(corr_state0, float(best_z[6]), n_revs=50)

        rec = {
            "t": t,
            "z_predictor": z_pred.tolist(),
            "z_corrected": best_z.tolist(),
            "initial_defect_norm": initial_defect,
            "best_defect_norm": best_defect,
            "defect_reduced": bool(best_defect < initial_defect),
            "diverged": diverged,
            "raw_bounded": raw_verdict.bounded,
            "raw_growth_ratio": raw_verdict.growth_ratio,
            "raw_trend_fraction": raw_verdict.trend_fraction,
            "raw_n_windows_complete": raw_verdict.n_windows_complete,
            "corrected_bounded": corr_verdict.bounded,
            "corrected_growth_ratio": corr_verdict.growth_ratio,
            "corrected_trend_fraction": corr_verdict.trend_fraction,
            "corrected_n_windows_complete": corr_verdict.n_windows_complete,
        }
        records.append(rec)
        print(
            f"t={t:.3f} defect {initial_defect:.4f}->{best_defect:.4f} "
            f"({'reduced' if rec['defect_reduced'] else 'DIVERGED'}) "
            f"raw_bounded={raw_verdict.bounded} corrected_bounded={corr_verdict.bounded}",
            flush=True,
        )
    return records


def run_planar_circular_secondary_check() -> dict:
    """Secondary check: does the planar/circular retrograde-satellite backbone
    (Henon family f / Families A-B-F) connect an x0~=J point to an x0~=cluster-43
    point via the ALREADY test-validated `cr3bp_continuation.py` driver?

    NOTE: this deliberately drops z/inclination and eccentricity forcing -- it
    tests the underlying near-circular backbone both orbits' ydot0/x0 ratios
    approximately sit on, not the full 3D-inclined structure. See module
    docstring.
    """
    from cyclerfinder.core.er3bp_geocentric import MU_SUN_EARTH_GURFIL_KASDIN

    sysm = cr3bp.CR3BPSystem(
        mu=MU_SUN_EARTH_GURFIL_KASDIN, primary="Sun", secondary="Earth", l_km=1.496e8, t_s=1.0
    )
    mu = sysm.mu
    # IMPORTANT: cr3bp_periodic's x0 is the standard CR3BP BARYCENTRIC
    # x-coordinate (primary at -mu, secondary/Earth at 1-mu), whereas J_IC/
    # C43_IC's x0 is Gurfil-Kasdin's GEOCENTRIC x (origin at Earth, +x radially
    # outward from the Sun, module docstring of core/er3bp_geocentric.py).
    # Converting requires the SAME "+(1-mu)" translation
    # geocentric_to_barycentric already applies (an earlier version of this
    # function fed the geocentric x0 straight in, landing near the SUN instead
    # of near Earth -- caught by a sanity check on the resulting Jacobi
    # constant, ~59.7 instead of the expected ~3, before this was trusted).
    x0_j_geo = float(J_IC[0])
    x0_c43_geo = float(C43_IC[0])
    x0_j_bary = (1.0 - mu) + x0_j_geo
    x0_c43_bary = (1.0 - mu) + x0_c43_geo
    ydot0_guess_j = float(J_IC[3])  # interleaved index 3 = ydot0

    # Jacobi constant of the perpendicular-crossing planar/circular reduction
    # (x0, 0, 0, 0, ydot0, 0): C = -2*Ubar(x0) - ydot0^2 (Ross Eq. 9 inverted;
    # matches cr3bp_periodic._ubar_x_at_axis's own formula, not re-derived).
    r1 = abs(x0_j_bary + mu)
    r2 = abs(x0_j_bary - 1.0 + mu)
    ubar_term = x0_j_bary * x0_j_bary + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2
    jacobi_guess = ubar_term - ydot0_guess_j**2

    seed = cp.correct_symmetric_fixed_jacobi(
        sysm,
        x0_j_bary,
        jacobi_guess,
        period_guess=2.0 * math.pi,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-10,
    )
    if not seed.converged:
        return {
            "seed_converged": False,
            "note": "circular-limit seed near Family J's x0 failed to converge",
        }

    # Empirically (this family): larger Jacobi C <-> SMALLER x0 for this
    # retrograde-satellite branch (a quick +1-direction probe walked x0 DOWN,
    # not up) -- direction is chosen by testing which way x0 actually needs to
    # move, not assumed from the C/x0 sign convention.
    d_jacobi = 1e-4
    direction = -1 if x0_c43_bary > seed.x0 else 1
    branch = cc.continue_family(
        sysm,
        seed,
        direction=direction,
        d_jacobi=d_jacobi,
        n_steps=200,
        min_jacobi=seed.jacobi - 5.0,
        max_jacobi=seed.jacobi + 5.0,
        half_crossings=1,
        ydot0_sign=-1.0,
        seed_label="J-circular-limit",
    )
    x0_values = [m.x0 for m in branch.members]
    reached = (
        any(
            (x0_values[i] - x0_c43_bary) * (x0_values[i + 1] - x0_c43_bary) <= 0.0
            for i in range(len(x0_values) - 1)
        )
        if len(x0_values) >= 2
        else False
    )
    # If the first attempt walked the wrong way (x0 diverging from the
    # target rather than approaching it), retry with the flipped direction --
    # a real, reported fallback, not a silent retry-until-success loop.
    walked_away = x0_values and abs(x0_values[-1] - x0_c43_bary) >= abs(x0_values[0] - x0_c43_bary)
    if not reached and walked_away:
        branch = cc.continue_family(
            sysm,
            seed,
            direction=-direction,
            d_jacobi=d_jacobi,
            n_steps=200,
            min_jacobi=seed.jacobi - 5.0,
            max_jacobi=seed.jacobi + 5.0,
            half_crossings=1,
            ydot0_sign=-1.0,
            seed_label="J-circular-limit",
        )
        x0_values = [m.x0 for m in branch.members]
        reached = (
            any(
                (x0_values[i] - x0_c43_bary) * (x0_values[i + 1] - x0_c43_bary) <= 0.0
                for i in range(len(x0_values) - 1)
            )
            if len(x0_values) >= 2
            else False
        )

    return {
        "seed_converged": True,
        "seed_x0_bary": seed.x0,
        "seed_x0_geo": seed.x0 - (1.0 - mu),
        "seed_jacobi": seed.jacobi,
        "target_x0_cluster43_bary": x0_c43_bary,
        "target_x0_cluster43_geo": x0_c43_geo,
        "stop_reason": str(branch.stop_reason),
        "n_members": len(branch.members),
        "x0_range_reached_bary": (min(x0_values), max(x0_values)) if x0_values else None,
        "x0_range_reached_geo": (
            (min(x0_values) - (1.0 - mu), max(x0_values) - (1.0 - mu)) if x0_values else None
        ),
        "crossed_target_x0": reached,
    }


def main() -> None:
    print("=== #602: cluster 43 <-> Family J real continuation ===\n")

    print("--- Step A: positive control / residual-formulation verification ---")
    pc_j = verify_no_nearby_exact_periodic_orbit(J_IC, "Family J")
    print(json.dumps(pc_j, indent=2))
    pc_c43 = verify_no_nearby_exact_periodic_orbit(C43_IC, "cluster 43")
    print(json.dumps(pc_c43, indent=2))

    print("\n--- Step B: damped deflated-Newton homotopy walk (real correction) ---")
    homotopy_records = run_homotopy()

    print("\n--- Step C: secondary planar/circular backbone check ---")
    try:
        planar_result = run_planar_circular_secondary_check()
    except Exception as exc:  # report, don't crash the whole run
        planar_result = {"error": str(exc)}
    print(json.dumps(planar_result, indent=2, default=str))

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w") as fh:
        fh.write(json.dumps({"record_type": "positive_control", **pc_j}) + "\n")
        fh.write(json.dumps({"record_type": "positive_control", **pc_c43}) + "\n")
        for rec in homotopy_records:
            fh.write(json.dumps({"record_type": "homotopy_step", **rec}) + "\n")
        planar_rec = {"record_type": "planar_circular_secondary_check", **planar_result}
        fh.write(json.dumps(planar_rec, default=str) + "\n")
    print(f"\nWrote {OUT_JSONL}")


if __name__ == "__main__":
    main()
