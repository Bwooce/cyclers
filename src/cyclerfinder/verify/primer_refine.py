"""Add-an-impulse refinement of a non-optimal impulsive coast (DIAGNOSTIC).

This module turns the *qualitative* primer-vector verdict of
:mod:`cyclerfinder.verify.primer` (``IMPROVABLE_ADD_IMPULSE`` when ``max|p| >
1`` on a coast) into a *quantitative* recoverable-ΔV estimate, via the classic
Lion & Handelsman (1968) **add-an-impulse** step.

Theory (citations)
------------------
* Lion, P. M., & Handelsman, M. (1968). *Primer Vector on Fixed-Time Impulsive
  Trajectories.* AIAA Journal 6(1), 127-132. DOI 10.2514/3.4452. — when ``|p|``
  exceeds unity on a coast, total ΔV is reduced by adding a midcourse impulse;
  the first-order best place/direction to seed it is the primer peak time
  ``t*`` and the unit primer direction ``p̂(t*)`` there.
* Prussing, J. E. (2010). *Primer Vector Theory and Applications*, Ch. 2 in
  Conway (ed.), *Spacecraft Trajectory Optimization*, Cambridge UP. — the
  standard textbook treatment: seed the added impulse from the primer peak,
  then re-optimise the (now N+1)-impulse fixed-endpoint schedule.

The refinement (fixed boundary states — like-for-like)
------------------------------------------------------
A single coast is a two-impulse transfer that carries the spacecraft from a
fixed boundary state at A — position ``r_A`` at time ``t_A``, with an
externally-imposed pre-departure velocity ``v_A_before`` — to a fixed boundary
state at B (``r_B`` at ``t_B``, post-arrival velocity ``v_B_after``). The
original schedule pays

    ΔV_orig = |v_dep - v_A_before| + |v_B_after - v_arr|

where ``(v_dep, v_arr)`` is the single Lambert arc A→B over ``[t_A, t_B]``.

The add-an-impulse refinement inserts a midcourse impulse at an interior time
``t_m`` and position ``r_m``, splitting the coast into two Lambert arcs A→m and
m→B. The **boundary states (endpoints, encounter times, and the imposed
boundary velocities) are held fixed**, so the comparison is like-for-like:

    ΔV_refined = |v1_Am - v_A_before|            (departure burn at A)
               + |v1_mB - v2_Am|                 (midcourse burn at m)
               + |v_B_after - v2_mB|             (arrival burn at B)

The free parameters are ``(τ, δr)`` where ``τ ∈ (0, 1)`` is the fractional
midcourse time (``t_m = t_A + τ (t_B - t_A)``) and ``δr ∈ ℝ³`` is the offset of
``r_m`` from the *original* (single-arc) position at that time. ``δr = 0,
τ → anything`` reproduces the original arc to numerical precision (the midcourse
burn vanishes), so the optimum can never be worse than the original. Seeded from
the primer peak (``τ* = t*/T`` and ``δr ∝ p̂(t*)``), scipy then minimises the
small 4-vector.

Honesty (binding)
-----------------
The recoverable-ΔV number this module returns is **OUR computation**
(diagnostic, provisional): first-order primer theory predicts that *some*
improvement exists, but the realised magnitude is produced by the optimiser
here, not by primer theory. Report it as "realised improvement (our optimiser)".
The Guzman 2002 multi-rev caveats of :mod:`cyclerfinder.verify.primer` apply
unchanged. This must NOT be published to the site or written to the catalogue.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.verify.primer import _coast_stm, _solve_primer_rate

Vec3 = NDArray[np.float64]

_REFINE_MAXITER: int = 1000
_REFINE_FTOL: float = 1.0e-10
_REFINE_XTOL: float = 1.0e-6
"""Nelder-Mead is used (the objective has |·| kinks at the burns, so a
gradient method is unreliable); these cap its iterations / convergence. The
maxiter is generous because the 4-D simplex needs several hundred evaluations
to converge on the kinked objective."""

_SEED_OFFSET_FRAC: float = 0.02
"""Initial midcourse position offset, as a fraction of |r_A|, along the primer
peak direction. Small but non-zero so the optimiser starts off the (zero-burn)
original arc; the optimiser then sizes it. Pure seeding only."""

_PENALTY: float = 1.0e9
"""Finite penalty returned when a trial midcourse split is Lambert-pathological,
so Nelder-Mead does not choke on inf/nan."""


@dataclass(frozen=True)
class AddImpulseRefinement:
    """Outcome of one add-an-impulse step on a single coast.

    Attributes
    ----------
    original_dv_kms:
        Two-impulse ΔV of the original single-arc coast (departure + arrival),
        km/s.
    refined_dv_kms:
        Three-impulse ΔV after inserting and optimising one midcourse impulse,
        km/s. ``<= original_dv_kms`` by construction (zero offset reproduces the
        original).
    recoverable_dv_kms:
        ``original_dv_kms - refined_dv_kms`` (km/s). OUR computed value — the
        realised improvement from this optimiser, NOT a primer-predicted bound.
    recoverable_fraction:
        ``recoverable_dv_kms / original_dv_kms`` (dimensionless). ``0`` when the
        original is degenerate (≈0 km/s).
    midcourse_time_frac:
        Optimised fractional midcourse time ``τ = (t_m - t_A)/(t_B - t_A)``.
    seed_peak_time_frac:
        Primer-peak fractional time used to seed the step (for the note).
    converged:
        ``True`` when the optimiser ran to a feasible (Lambert-buildable)
        optimum.
    """

    original_dv_kms: float
    refined_dv_kms: float
    recoverable_dv_kms: float
    recoverable_fraction: float
    midcourse_time_frac: float
    seed_peak_time_frac: float
    converged: bool


def _primer_peak(
    r0: Vec3,
    v0: Vec3,
    p0_hat: Vec3,
    p1_hat: Vec3,
    duration_s: float,
    mu: float,
    *,
    n_samples: int,
) -> tuple[float, Vec3]:
    """Time and unit direction of the primer peak on a coast.

    Reconstructs the full primer *vector* ``p(t)`` (not just its magnitude)
    along the coast via the same STM/BVP machinery as
    :func:`cyclerfinder.verify.primer.primer_on_coast`, and returns the
    fractional time of ``max|p|`` together with the unit primer direction there
    — the Lion & Handelsman seed for the added impulse.
    """
    p0 = np.asarray(p0_hat, dtype=np.float64)
    p1 = np.asarray(p1_hat, dtype=np.float64)
    p0 = p0 / np.linalg.norm(p0)
    p1 = p1 / np.linalg.norm(p1)

    times, _ref, stms = _coast_stm(r0, v0, duration_s, mu, n_samples=n_samples)
    phi_end = stms[-1]
    phi_rr = phi_end[:3, :3]
    phi_rv = phi_end[:3, 3:]
    # Robust to near-integer-rev singularity of Φ_rv (Saloglu 2023 Sec. III.F):
    # continuity (truncated-SVD min-norm) primer rate instead of a direct solve.
    pdot0, _ill, _rcond = _solve_primer_rate(phi_rr, phi_rv, p0, p1)
    state0 = np.concatenate([p0, pdot0])

    p_vectors = np.empty((times.shape[0], 3), dtype=np.float64)
    for k in range(times.shape[0]):
        p_vectors[k] = stms[k][:3, :] @ state0
    mags = np.linalg.norm(p_vectors, axis=1)
    k_max = int(np.argmax(mags))
    peak_dir = p_vectors[k_max]
    peak_dir = peak_dir / np.linalg.norm(peak_dir)
    peak_frac = float(times[k_max] / duration_s)
    return peak_frac, peak_dir


def _two_impulse_dv(
    r_a: Vec3,
    r_b: Vec3,
    tof_s: float,
    v_a_before: Vec3,
    v_b_after: Vec3,
    mu: float,
) -> float:
    """Departure + arrival ΔV of the single Lambert arc A→B (km/s)."""
    sol = lambert(r_a, r_b, tof_s, mu=mu)[0]
    v_dep = np.asarray(sol.v1, dtype=np.float64)
    v_arr = np.asarray(sol.v2, dtype=np.float64)
    return float(np.linalg.norm(v_dep - v_a_before) + np.linalg.norm(v_b_after - v_arr))


def add_impulse_refine_coast(
    r_a: Vec3,
    r_b: Vec3,
    tof_s: float,
    v_a_before: Vec3,
    v_b_after: Vec3,
    *,
    primer_peak_frac: float,
    primer_peak_dir: Vec3,
    mu: float = MU_SUN_KM3_S2,
) -> AddImpulseRefinement:
    """One Lion & Handelsman add-an-impulse step on a single fixed-endpoint coast.

    Holds the boundary states ``(r_a, v_a_before)`` and ``(r_b, v_b_after)`` and
    the encounter times fixed, inserts a midcourse impulse seeded at the primer
    peak, and numerically minimises the three-impulse total ΔV over the
    midcourse ``(τ, δr)``. Returns the original vs refined ΔV and the realised
    (our-optimiser) recoverable ΔV.

    Parameters
    ----------
    r_a, r_b:
        Fixed coast endpoints (km), in the same inertial frame.
    tof_s:
        Coast duration (s), ``> 0``.
    v_a_before:
        Imposed pre-departure boundary velocity at A (km/s) — for a cycler leg
        this is the spacecraft arrival velocity from the previous leg.
    v_b_after:
        Imposed post-arrival boundary velocity at B (km/s) — the spacecraft
        departure velocity for the next leg.
    primer_peak_frac:
        Fractional time ``t*/T`` of the primer peak (the seed midcourse time).
    primer_peak_dir:
        Unit primer direction ``p̂(t*)`` (the seed offset direction).
    mu:
        Central-body gravitational parameter, km³/s².

    Returns
    -------
    AddImpulseRefinement
    """
    if tof_s <= 0.0:
        raise ValueError(f"coast duration must be positive, got {tof_s}")

    r_a = np.asarray(r_a, dtype=np.float64)
    r_b = np.asarray(r_b, dtype=np.float64)
    v_a_before = np.asarray(v_a_before, dtype=np.float64)
    v_b_after = np.asarray(v_b_after, dtype=np.float64)
    peak_dir = np.asarray(primer_peak_dir, dtype=np.float64)
    peak_dir = peak_dir / np.linalg.norm(peak_dir)

    original_dv = _two_impulse_dv(r_a, r_b, tof_s, v_a_before, v_b_after, mu)

    # Reference (single-arc) position at the midcourse fractional time, so δr is
    # measured as an offset from the original trajectory (δr = 0 reproduces it).
    sol0 = lambert(r_a, r_b, tof_s, mu=mu)[0]
    v_dep0 = np.asarray(sol0.v1, dtype=np.float64)

    def _ref_position(tau: float) -> Vec3:
        r_mid, _v = propagate(r_a, v_dep0, tau * tof_s, mu)
        return np.asarray(r_mid, dtype=np.float64)

    def _objective(x: NDArray[np.float64]) -> float:
        tau = float(x[0])
        if not (1.0e-3 < tau < 1.0 - 1.0e-3):
            return _PENALTY
        delta_r = x[1:4]
        r_m = _ref_position(tau) + delta_r
        t1 = tau * tof_s
        t2 = (1.0 - tau) * tof_s
        try:
            sol_am = lambert(r_a, r_m, t1, mu=mu)[0]
            sol_mb = lambert(r_m, r_b, t2, mu=mu)[0]
        except (LambertConvergenceError, LambertGeometryError, ValueError):
            return _PENALTY
        v1_am = np.asarray(sol_am.v1, dtype=np.float64)
        v2_am = np.asarray(sol_am.v2, dtype=np.float64)
        v1_mb = np.asarray(sol_mb.v1, dtype=np.float64)
        v2_mb = np.asarray(sol_mb.v2, dtype=np.float64)
        dv_dep = float(np.linalg.norm(v1_am - v_a_before))
        dv_mid = float(np.linalg.norm(v1_mb - v2_am))
        dv_arr = float(np.linalg.norm(v_b_after - v2_mb))
        return dv_dep + dv_mid + dv_arr

    # Seed: primer-peak time, small offset along p̂(t*).
    tau_seed = float(np.clip(primer_peak_frac, 0.05, 0.95))
    offset_seed = _SEED_OFFSET_FRAC * float(np.linalg.norm(r_a)) * peak_dir
    x0 = np.concatenate([[tau_seed], offset_seed])

    res = minimize(
        _objective,
        x0,
        method="Nelder-Mead",
        options={"maxiter": _REFINE_MAXITER, "fatol": _REFINE_FTOL, "xatol": _REFINE_XTOL},
    )

    refined_dv = float(res.fun)
    converged = bool(res.success) and refined_dv < _PENALTY
    # Never report worse than the original: the original arc is always a
    # feasible point (the zero-offset limit), so clamp to it defensively.
    if refined_dv > original_dv:
        refined_dv = original_dv
        tau_opt = primer_peak_frac
    else:
        tau_opt = float(res.x[0])

    recoverable = original_dv - refined_dv
    frac = recoverable / original_dv if original_dv > 1.0e-12 else 0.0
    return AddImpulseRefinement(
        original_dv_kms=original_dv,
        refined_dv_kms=refined_dv,
        recoverable_dv_kms=recoverable,
        recoverable_fraction=frac,
        midcourse_time_frac=tau_opt,
        seed_peak_time_frac=float(primer_peak_frac),
        converged=converged,
    )


def refine_coast_from_states(
    r0: Vec3,
    v0_depart: Vec3,
    r_b: Vec3,
    tof_s: float,
    v_a_before: Vec3,
    v_b_after: Vec3,
    p0_hat: Vec3,
    p1_hat: Vec3,
    *,
    mu: float = MU_SUN_KM3_S2,
    n_samples: int = 300,
) -> AddImpulseRefinement:
    """Compute the primer peak on a coast, then add-an-impulse refine it.

    Convenience wrapper: runs the primer BVP (to find the peak time/direction
    seed) and then :func:`add_impulse_refine_coast`. The coast's heliocentric
    departure state ``(r0, v0_depart)`` defines the original arc; ``(p0_hat,
    p1_hat)`` are the bounding unit ΔV directions used to pin the primer BVP.

    Parameters
    ----------
    r0, v0_depart:
        Heliocentric state at the start of the coast (just after departure
        burn). ``r0`` is also the fixed endpoint A.
    r_b:
        Fixed endpoint B (km).
    tof_s:
        Coast duration (s).
    v_a_before, v_b_after:
        Imposed boundary velocities at A and B (km/s); see
        :func:`add_impulse_refine_coast`.
    p0_hat, p1_hat:
        Unit ΔV directions at the bounding impulses (the primer BVP pins).
    mu:
        Central-body gravitational parameter, km³/s².
    n_samples:
        Grid resolution for locating the primer peak.

    Returns
    -------
    AddImpulseRefinement
    """
    peak_frac, peak_dir = _primer_peak(
        r0, v0_depart, p0_hat, p1_hat, tof_s, mu, n_samples=n_samples
    )
    return add_impulse_refine_coast(
        r0,
        r_b,
        tof_s,
        v_a_before,
        v_b_after,
        primer_peak_frac=peak_frac,
        primer_peak_dir=peak_dir,
        mu=mu,
    )


__all__ = [
    "AddImpulseRefinement",
    "add_impulse_refine_coast",
    "refine_coast_from_states",
]
