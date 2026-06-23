"""#307 Task 2: automated DSM placement (Vasile-Conway 2006 MGA-DSM transcription).

Position-closing deep-space-manoeuvre (DSM) placement for repairing a deficient
gravity-assist leg. The standard MGA-DSM transcription (Vasile & Conway 2006 §3.2)
splits a leg in two:

* **Arc 1** — depart the departure body with ``v_body + vinf_dep`` and propagate
  ballistically (Kepler initial-value) for ``eta * tof`` to the DSM point.
* **Arc 2** — a **Lambert** boundary-value solve from the DSM point to the arrival
  body over the remaining ``(1 - eta) * tof``. Because arc 2 *targets* the arrival
  body, the arrival POSITION closes by construction (unlike the forward-propagate
  DSM executor in :mod:`cyclerfinder.genome.epoch_aware_genome`, which leaves a
  position miss that only shows up in the independent cross-check).

The DSM Δv is the velocity discontinuity at the split point
(``|v_arc2_depart - v_arc1_arrive|``) — the optimised actuator cost. The free
variables are the DSM timing ``eta`` and the departure ``vinf`` vector (magnitude
+ two direction angles); the optimiser drives the *arrival* V∞ to a target while
minimising Δv.

Standalone by design: this does NOT modify the epoch-locked closure driver's
existing (forward-propagate) DSM semantics; it is a position-closing repair lane
the precursor matcher can call. (The DSMSpec docstring in epoch_aware_genome
describes this Lambert-arc-2 model but the code there forward-propagates — a
docstring/code mismatch logged for separate cleanup.)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import Bounds, differential_evolution

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler import KeplerConvergenceError, propagate
from cyclerfinder.core.lambert import lambert


@dataclass(frozen=True)
class DsmLegResult:
    """Outcome of a single Vasile-Conway DSM leg evaluation.

    All velocities km/s, positions km, ``eta`` dimensionless in (0, 1).
    ``arrival_pos_residual_km`` is the Lambert/propagate consistency check (the
    arrival position closes by construction, so this is ~0 to integrator
    precision — a non-zero value flags a degenerate/near-singular Lambert).
    """

    eta: float
    vinf_dep_kms: float
    dsm_dv_kms: float
    arrival_vinf_kms: float
    arrival_pos_residual_km: float


def _vinf_vec(mag: float, theta: float, phi: float) -> NDArray[np.float64]:
    """Spherical (magnitude, longitude, latitude) -> Cartesian km/s vector."""
    cphi = np.cos(phi)
    return mag * np.array(
        [cphi * np.cos(theta), cphi * np.sin(theta), np.sin(phi)],
        dtype=np.float64,
    )


def evaluate_dsm_leg(
    r_dep_km: NDArray[np.float64],
    v_body_dep_kms: NDArray[np.float64],
    r_arr_km: NDArray[np.float64],
    v_body_arr_kms: NDArray[np.float64],
    tof_s: float,
    *,
    eta: float,
    vinf_dep_vec_kms: NDArray[np.float64],
    mu_km3_s2: float = MU_SUN_KM3_S2,
) -> DsmLegResult:
    """Evaluate one Vasile-Conway DSM leg (arc 1 propagate + arc 2 Lambert)."""
    if not (0.0 < eta < 1.0):
        raise ValueError(f"eta must be in (0, 1), got {eta}")
    if tof_s <= 0.0:
        raise ValueError(f"tof_s must be positive, got {tof_s}")
    r_dep = np.asarray(r_dep_km, dtype=np.float64)
    r_arr = np.asarray(r_arr_km, dtype=np.float64)
    v_body_arr = np.asarray(v_body_arr_kms, dtype=np.float64)
    v_dep = np.asarray(v_body_dep_kms, dtype=np.float64) + np.asarray(
        vinf_dep_vec_kms, dtype=np.float64
    )

    t1 = eta * tof_s
    t2 = tof_s - t1
    r_dsm, v_dsm_minus = propagate(r_dep, v_dep, t1, mu=mu_km3_s2)
    r_dsm = np.asarray(r_dsm, dtype=np.float64)

    sols = lambert(r_dsm, r_arr, t2, mu=mu_km3_s2, prograde=True, max_revs=0)
    arc2 = sols[0]
    dsm_dv = float(np.linalg.norm(arc2.v1 - np.asarray(v_dsm_minus, dtype=np.float64)))
    arrival_vinf = float(np.linalg.norm(arc2.v2 - v_body_arr))

    # Consistency: re-propagate arc 2 forward and confirm it lands on the arrival
    # body (Lambert boundary-value vs Kepler initial-value agreement). A strongly
    # hyperbolic arc-2 solution can defeat the universal-variable Newton; treat
    # that as a non-closing (degenerate) leg so callers/optimisers reject it.
    try:
        r_check, _ = propagate(r_dsm, arc2.v1, t2, mu=mu_km3_s2)
        pos_res = float(np.linalg.norm(np.asarray(r_check, dtype=np.float64) - r_arr))
    except KeplerConvergenceError:
        pos_res = float("inf")

    return DsmLegResult(
        eta=float(eta),
        vinf_dep_kms=float(np.linalg.norm(vinf_dep_vec_kms)),
        dsm_dv_kms=dsm_dv,
        arrival_vinf_kms=arrival_vinf,
        arrival_pos_residual_km=pos_res,
    )


def optimize_dsm_leg(
    r_dep_km: NDArray[np.float64],
    v_body_dep_kms: NDArray[np.float64],
    r_arr_km: NDArray[np.float64],
    v_body_arr_kms: NDArray[np.float64],
    tof_s: float,
    target_arrival_vinf_kms: float,
    *,
    mu_km3_s2: float = MU_SUN_KM3_S2,
    vinf_dep_max_kms: float = 8.0,
    vinf_match_weight: float = 5.0,
    seed: int = 0,
    maxiter: int = 40,
    popsize: int = 15,
) -> DsmLegResult:
    """Place a DSM minimising Δv while driving arrival V∞ to a target.

    Decision variables: ``eta`` (DSM timing), departure ``|vinf|`` + two direction
    angles. Objective ``dsm_dv + vinf_match_weight * |arrival_vinf - target|``
    (position closure is automatic via the Lambert arc 2). Uses
    ``differential_evolution`` (global, gradient-free); ``seed`` makes it
    deterministic.
    """
    bounds = Bounds(
        [0.05, 0.0, 0.0, -0.5 * np.pi],  # eta, |vinf_dep|, longitude, latitude (lo)
        [0.95, vinf_dep_max_kms, 2.0 * np.pi, 0.5 * np.pi],  # (hi)
    )

    def objective(x: NDArray[np.float64]) -> float:
        eta = float(x[0])
        vinf = _vinf_vec(float(x[1]), float(x[2]), float(x[3]))
        try:
            res = evaluate_dsm_leg(
                r_dep_km,
                v_body_dep_kms,
                r_arr_km,
                v_body_arr_kms,
                tof_s,
                eta=eta,
                vinf_dep_vec_kms=vinf,
                mu_km3_s2=mu_km3_s2,
            )
        except (ValueError, ZeroDivisionError, FloatingPointError):
            return 1.0e9
        if not np.isfinite(res.dsm_dv_kms) or not np.isfinite(res.arrival_vinf_kms):
            return 1.0e9
        return res.dsm_dv_kms + vinf_match_weight * abs(
            res.arrival_vinf_kms - target_arrival_vinf_kms
        )

    result = differential_evolution(
        objective,
        bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        tol=1e-4,
        polish=True,
    )
    eta = float(result.x[0])
    vinf = _vinf_vec(float(result.x[1]), float(result.x[2]), float(result.x[3]))
    return evaluate_dsm_leg(
        r_dep_km,
        v_body_dep_kms,
        r_arr_km,
        v_body_arr_kms,
        tof_s,
        eta=eta,
        vinf_dep_vec_kms=vinf,
        mu_km3_s2=mu_km3_s2,
    )
