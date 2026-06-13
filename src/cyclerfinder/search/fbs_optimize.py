r"""ΔV-minimising multi-leg DSM-chain optimiser driven by the FBS analytic Jacobian.

This is the **decision-gate** lane for #243: it tests the Ellison-2018
forward-backward-shooting (FBS) transcription in its ACTUAL intended role — as the
gradient source for gradient-based ΔV OPTIMISATION of multi-DSM / multi-gravity-assist
trajectories (the EMTG/SNOPT use case) — NOT as a single-leg feasibility solver (that
was #242, the mis-scoped test).

The transcription (Ellison Sec. II, our massless Path-B)
--------------------------------------------------------
A chain of ``M`` two-body legs over a FIXED body sequence with FIXED epochs (departure
epoch and per-leg ToFs are inputs; the gradient question is isolated to the
velocity/impulse NLP). The decision vector is the EMTG-style impulses + shared
match-point boundary velocities:

    x = [ Δv_0..Δv_{M-1} (3M) | v_0..v_M (3(M+1)) ]

where ``v_j`` is the SHARED heliocentric boundary velocity at body ``j`` (Ellison's
match-point variables): leg ``i`` uses ``v_i`` as its left-boundary velocity and
``v_{i+1}`` as its arrival velocity. ``v_0`` is the heliocentric departure velocity
(its v∞ rides the launch energy) and ``v_M`` is the terminal arrival velocity; the
interior ``v_j`` (0 < j < M) are the flyby match points the optimiser is free to move.

Objective:  f(x) = Σ_i ‖Δv_i‖   (+ ‖v_M - v_planet_arr‖ when ``rendezvous``)
Constraints (equality, = 0):  the stacked per-leg match-point defects ``chain_defect``,
    with the analytic block-sparse ``chain_defect_jacobian`` (#226) supplying ∂c/∂x
    exactly.

This is genuinely UNDER-determined: ``3M + 3(M+1)`` variables versus ``6M`` equality
constraints leaves ``3M + 3`` optimisation degrees of freedom (the v∞ vectors plus the
impulse null space), so SLSQP minimises ΔV over a real null space rather than just
root-finding — exactly the optimisation problem FBS analytic gradients exist to speed up.

Head-to-head subject for ``scripts/fbs_optimizer_fair_trial.py``: the SAME NLP is handed
to SLSQP twice — once with the analytic constraint Jacobian (``use_analytic_jac=True``,
the FBS lane) and once with an SLSQP finite-difference Jacobian
(``use_analytic_jac=False``, the Lambert+FD lane's gradient style) — so robustness /
cost / optimum-quality differences are attributable to the GRADIENT SOURCE alone, which
is precisely the FBS claimed advantage.

NO catalogue writeback: this is a method evaluation. Pure w.r.t. ``data/``.

Source: D. H. Ellison et al., "Analytic Gradient Computation for Bounded-Impulse
Trajectory Models Using Two-Sided Shooting," JGCD 41(7), 2018, doi:10.2514/1.G003077.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import NonlinearConstraint, minimize

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.fbs_match_point import (
    FbsLeg,
    chain_defect,
    chain_defect_jacobian,
)
from cyclerfinder.core.kepler import KeplerError

Vec3 = NDArray[np.float64]

# Large finite scaled-defect penalty returned when a boundary velocity drives the
# two-body propagation past Newton convergence (hyperbolic probe). Big enough to be
# unambiguously infeasible, finite so SLSQP does not choke on inf/nan.
_DEFECT_PENALTY: float = 1.0e3


@dataclass(frozen=True)
class ChainLegSpec:
    """Fixed (non-decision) geometry of one chain leg.

    Attributes
    ----------
    r0, rf:
        Left / right boundary heliocentric positions (km). FIXED (body positions at
        the leg's departure / arrival epochs).
    tof_s, alpha:
        Total ToF (s) and burn fraction. FIXED.

    The boundary VELOCITIES are decision variables of the chain (the shared ``v_j``),
    not stored here.
    """

    r0: Vec3
    rf: Vec3
    tof_s: float
    alpha: float


@dataclass(frozen=True)
class FbsOptimizeResult:
    """Outcome of one FBS-gradient (or FD) chain ΔV optimisation.

    Attributes
    ----------
    total_dv_kms:
        The minimised objective Σ‖Δv_i‖ (+ arrival v∞ when rendezvous), km/s.
    dv_per_leg_kms:
        Per-leg converged impulse magnitudes, km/s.
    dvs:
        Converged per-leg impulse vectors (km/s).
    boundary_vs:
        Converged shared boundary velocities ``v_0..v_M`` (M+1 vectors, km/s).
    max_defect:
        ``max|chain match-point defect|`` at the solution on the NON-DIMENSIONALISED
        scale (position rows / AU, velocity rows / v_circ). Feasibility quantity.
    feasible:
        ``max_defect < feas_tol``.
    success:
        SLSQP's own convergence flag (audit).
    nfev, njev, nit:
        SLSQP objective-eval, objective-jac-eval and iteration counts.
    constr_nfev, constr_njev:
        Constraint function- and (analytic) jacobian-eval counts (tallied here; SLSQP
        does not expose these). In the FD lane ``constr_njev`` is 0 — the constraint
        Jacobian is finite-differenced inside SLSQP, costing extra ``constr_nfev``.
    wall_s:
        Wall-clock seconds for the solve (the COST metric).
    used_analytic_jac:
        Whether the analytic FBS Jacobian was supplied (True) or finite difference
        was used (False).
    """

    total_dv_kms: float
    dv_per_leg_kms: tuple[float, ...]
    dvs: tuple[Vec3, ...]
    boundary_vs: tuple[Vec3, ...]
    max_defect: float
    feasible: bool
    success: bool
    nfev: int
    njev: int
    nit: int
    constr_nfev: int
    constr_njev: int
    wall_s: float
    used_analytic_jac: bool


def _defect_scale(legs: tuple[ChainLegSpec, ...], mu: float) -> NDArray[np.float64]:
    """Per-row non-dimensionalisation of the stacked defect (pos/AU, vel/v_circ).

    Mirrors ``dsm_leg_correct_fbs``: position rows by AU, velocity rows by the local
    circular speed at each leg's left boundary, so SLSQP sees a well-scaled equality
    constraint vector and ``feas_tol`` is a meaningful threshold across legs.
    """
    rows: list[float] = []
    for leg in legs:
        v_scale = float(np.sqrt(mu / float(np.linalg.norm(leg.r0))))
        rows.extend([AU_KM, AU_KM, AU_KM, v_scale, v_scale, v_scale])
    return np.asarray(rows, dtype=np.float64)


def _fbs_legs(
    legs: tuple[ChainLegSpec, ...],
    boundary_vs: tuple[Vec3, ...],
    mu: float,
) -> tuple[FbsLeg, ...]:
    """Build the FBS legs from the fixed geometry + the shared boundary velocities."""
    return tuple(
        FbsLeg(
            r0=leg.r0,
            v0=np.asarray(boundary_vs[i], dtype=np.float64),
            rf=leg.rf,
            vf=np.asarray(boundary_vs[i + 1], dtype=np.float64),
            tof_s=leg.tof_s,
            alpha=leg.alpha,
            mu=mu,
        )
        for i, leg in enumerate(legs)
    )


def optimize_chain_fbs(
    legs: tuple[ChainLegSpec, ...],
    dv0_per_leg: tuple[Vec3, ...],
    boundary_v0: tuple[Vec3, ...],
    *,
    mu: float = MU_SUN_KM3_S2,
    rendezvous_vplanet: Vec3 | None = None,
    use_analytic_jac: bool = True,
    feas_tol: float = 1.0e-8,
    maxiter: int = 200,
    ftol: float = 1.0e-9,
) -> FbsOptimizeResult:
    r"""Minimise chain ΔV with the per-leg match-point defects as equality constraints.

    Decision vector ``x = [Δv_0..Δv_{M-1} (3M) | v_0..v_M (3(M+1))]`` — the impulses
    and the shared match-point boundary velocities (Ellison's match-point variables).
    Objective ``Σ‖Δv_i‖`` (plus the terminal arrival v∞ when ``rendezvous_vplanet``).
    The equality constraint is the stacked non-dimensionalised chain match-point defect
    (:func:`...chain_defect`); its Jacobian is the ANALYTIC
    :func:`...chain_defect_jacobian` when ``use_analytic_jac`` (the FBS lane) or SLSQP
    finite-difference when not (the FD baseline). The ONLY difference between the two
    lanes is the constraint-gradient source — so robustness / cost / optimum differences
    are attributable to it.

    Parameters
    ----------
    legs:
        Fixed per-leg geometry (boundary positions, ToF, alpha). ``M`` legs.
    dv0_per_leg:
        Seeds for the ``M`` impulse vectors (km/s).
    boundary_v0:
        Seeds for the ``M+1`` shared boundary velocities ``v_0..v_M`` (km/s).
    mu:
        Central-body gravitational parameter.
    rendezvous_vplanet:
        If given, the terminal body's heliocentric velocity (km/s); the arrival v∞
        magnitude ``‖v_M - v_planet‖`` is added to the objective (a rendezvous cost).
        ``None`` (default) optimises only the interior impulses (flyby finish).
    use_analytic_jac:
        Analytic FBS constraint Jacobian (True) vs SLSQP finite difference (False).
    feas_tol, maxiter, ftol:
        Feasibility threshold on the scaled defect, SLSQP iteration cap, objective
        tolerance.
    """
    m = len(legs)
    if m < 1:
        raise ValueError("chain must contain at least one leg")
    if len(dv0_per_leg) != m:
        raise ValueError(f"dv0_per_leg must have {m} entries, got {len(dv0_per_leg)}")
    if len(boundary_v0) != m + 1:
        raise ValueError(f"boundary_v0 must have {m + 1} entries, got {len(boundary_v0)}")

    scale = _defect_scale(legs, mu)
    n_dv = 3 * m
    n_v = 3 * (m + 1)
    counts = {"obj": 0, "grad": 0, "con": 0, "cjac": 0}

    def _split(x: NDArray[np.float64]) -> tuple[tuple[Vec3, ...], tuple[Vec3, ...]]:
        dvs = tuple(np.asarray(x[3 * i : 3 * i + 3], dtype=np.float64) for i in range(m))
        bvs = tuple(
            np.asarray(x[n_dv + 3 * j : n_dv + 3 * j + 3], dtype=np.float64) for j in range(m + 1)
        )
        return dvs, bvs

    def _objective(x: NDArray[np.float64]) -> float:
        counts["obj"] += 1
        dvs, bvs = _split(x)
        total = float(sum(float(np.linalg.norm(dv)) for dv in dvs))
        if rendezvous_vplanet is not None:
            vinf_arr = np.asarray(bvs[-1], dtype=np.float64) - np.asarray(
                rendezvous_vplanet, dtype=np.float64
            )
            total += float(np.linalg.norm(vinf_arr))
        return total

    def _objective_grad(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["grad"] += 1
        dvs, bvs = _split(x)
        g = np.zeros(n_dv + n_v, dtype=np.float64)
        for i, dv in enumerate(dvs):
            n = float(np.linalg.norm(dv))
            if n > 0.0:
                g[3 * i : 3 * i + 3] = dv / n
        if rendezvous_vplanet is not None:
            vinf_arr = np.asarray(bvs[-1], dtype=np.float64) - np.asarray(
                rendezvous_vplanet, dtype=np.float64
            )
            n = float(np.linalg.norm(vinf_arr))
            if n > 0.0:
                g[n_dv + 3 * m : n_dv + 3 * (m + 1)] = vinf_arr / n
        return g

    def _constraint(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["con"] += 1
        dvs, bvs = _split(x)
        fbs_legs = _fbs_legs(legs, bvs, mu)
        try:
            c = chain_defect(fbs_legs, dvs)  # (6M,)
        except KeplerError:
            # A bad boundary velocity (e.g. SLSQP probing a hyperbolic v0/vf) can
            # drive the universal-variable Newton past convergence. Return a large
            # finite defect of the right shape so the solver rejects the direction
            # rather than crashing — the SAME handling both lanes see (the Lambert
            # lane's dsm_chain_correct does the equivalent penalty substitution).
            return np.full(6 * m, _DEFECT_PENALTY, dtype=np.float64)
        return c / scale

    def _constraint_jac(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["cjac"] += 1
        dvs, bvs = _split(x)
        fbs_legs = _fbs_legs(legs, bvs, mu)
        try:
            # chain_defect_jacobian already returns exactly the (6M) x (3M + 3(M+1))
            # block-sparse Jacobian in the decision layout [Δv | v_0..v_M].
            jac = chain_defect_jacobian(fbs_legs, dvs)
        except KeplerError:
            return np.zeros((6 * m, n_dv + n_v), dtype=np.float64)
        return jac / scale[:, None]

    x0 = np.concatenate(
        [np.asarray(d, dtype=np.float64) for d in dv0_per_leg]
        + [np.asarray(v, dtype=np.float64) for v in boundary_v0]
    )

    jac_arg: Any = _constraint_jac if use_analytic_jac else "2-point"
    constraint = NonlinearConstraint(_constraint, 0.0, 0.0, jac=jac_arg)

    t_start = time.perf_counter()
    res: Any = minimize(
        _objective,
        x0,
        jac=_objective_grad,
        method="SLSQP",
        constraints=[constraint],
        options={"maxiter": maxiter, "ftol": ftol},
    )
    wall = time.perf_counter() - t_start

    x_sol = np.asarray(res.x, dtype=np.float64)
    dvs_sol, bvs_sol = _split(x_sol)
    fbs_legs = _fbs_legs(legs, bvs_sol, mu)
    try:
        defect = chain_defect(fbs_legs, dvs_sol) / scale
        max_defect = float(np.max(np.abs(defect)))
    except KeplerError:
        # SLSQP landed on a hyperbolic-propagation point (infeasible): report the
        # penalty defect so the run is correctly classed infeasible, not crashed.
        max_defect = _DEFECT_PENALTY
    total_dv = _objective(x_sol)
    dv_per_leg = tuple(float(np.linalg.norm(dv)) for dv in dvs_sol)

    return FbsOptimizeResult(
        total_dv_kms=total_dv,
        dv_per_leg_kms=dv_per_leg,
        dvs=dvs_sol,
        boundary_vs=bvs_sol,
        max_defect=max_defect,
        feasible=bool(max_defect < feas_tol),
        success=bool(res.success),
        nfev=int(res.get("nfev", 0)),
        njev=int(res.get("njev", 0)),
        nit=int(res.get("nit", 0)),
        constr_nfev=int(counts["con"]),
        constr_njev=int(counts["cjac"]),
        wall_s=float(wall),
        used_analytic_jac=use_analytic_jac,
    )
