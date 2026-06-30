"""V3 independent-integrator QP-torus invariance gauntlet (#319 Phase 2 / #320 follow-on).

Spec reinterpretation for QP-tori
---------------------------------
The strict-periodic V3 (:mod:`v3_3d_periodic`) re-propagates a CR3BP periodic-orbit
IC under an INDEPENDENT integrator (REBOUND IAS15 Gauss-Radau) and checks the orbit
signature survives the integrator swap — "same model, different integrator". V1_qp /
V2_qp establish a QP-torus's Fourier-mode invariance under the project's scipy
``cr3bp.propagate`` (DOP853) flow. **V3_qp is the same invariance equation under IAS15**:

    Does ``phi^{IAS15}_{(k+1)·t_strob}(u(theta))`` still lie on
    ``u(theta + (k+1)·rho)`` — and does IAS15 AGREE with the DOP853 result V2 used?

If the IAS15 invariance residual stays bounded AND IAS15 reproduces DOP853's
propagated endpoint to a tight floor, then the V1/V2 torus signature is a REAL
dynamical property of the CR3BP, not an artifact of the DOP853-based corrector +
invariance stack. If IAS15 and DOP853 DISAGREE, V2's bounded invariance was
integrator noise and the candidate must be retired.

This is the QP-torus analogue of :func:`v3_3d_periodic.run_v3_periodic_3d`'s
IAS15-vs-DOP853 cross-check. V4 (real-ephemeris HFEM) is a separate, later gate.

Independent-integrator architecture (the whole point)
-----------------------------------------------------
* V1_qp / V2_qp propagate with :func:`cyclerfinder.core.cr3bp.propagate` — scipy
  ``solve_ivp`` DOP853 of the analytic CR3BP rotating-frame EOM.
* V3_qp propagates with :func:`v3_3d_periodic._ias15_propagate_cr3bp_rotating` —
  REBOUND IAS15 (rotating-frame callback; scipy LSODA fallback if REBOUND absent).

A fundamentally different integrator family. The torus invariant circle
(:func:`evaluate_invariant_circle`) and the rotation-number comparison are identical
to V2_qp — only the propagator changes, isolating the integrator contribution.

Floor rationale
---------------
* ``drift_floor`` (the IAS15 invariance residual bar) reuses V2_qp's empirically-
  calibrated ``V2_QP_DRIFT_FLOOR`` (5e-2 nondim; Olikara 2016 §4 + #319 Phase-1
  N=2 calibration) — V3 should not be MORE permissive than V2.
* ``agreement_floor`` (the IAS15-vs-DOP853 terminal disagreement bar) is the
  integrator cross-check: two 1e-12-tol integrators of the same CR3BP EOM over
  ~3·t_strob (~30 TU) of a near-Neimark-Sacker torus agree to well under 1e-3
  nondim absent a real model discrepancy. Default 1e-3 nondim (judgment call,
  one order under the V2 drift floor); calibrated against the #320 SILVER tori.

A V3_qp PASS does NOT admit to the catalogue (V4 real-eph + human review follow).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v1_qp import _off_grid_thetas
from cyclerfinder.data.validation.v2_qp import (
    V2_QP_DRIFT_FLOOR,
    V2_QP_N_CYCLES_MIN,
    V2_QP_OFF_GRID_PER_CYCLE,
)
from cyclerfinder.data.validation.v3_3d_periodic import _ias15_propagate_cr3bp_rotating
from cyclerfinder.genome.qp_tori import QPTorus, evaluate_invariant_circle

#: V3_qp integrator-agreement floor (nondim): max |IAS15 - DOP853| terminal
#: disagreement over the off-grid samples. Calibrated against the #320 SILVERs.
V3_QP_AGREEMENT_FLOOR: float = 1.0e-3


@dataclass(frozen=True)
class V3VerdictQP:
    """Result of the V3 independent-integrator (IAS15) QP-torus invariance gauntlet."""

    candidate_id: str
    integrator: str  # the IAS15 integrator label (or LSODA fallback)
    n_cycles_requested: int
    n_cycles_longitudinal_propagated: int
    #: Per-cycle IAS15 invariance residual (nondim L_inf over off-grid samples).
    per_cycle_invariance_residual_ias15: tuple[float, ...]
    max_invariance_drift_ias15: float
    max_invariance_drift_ias15_km: float
    #: Per-cycle max |IAS15 - DOP853| terminal disagreement (the integrator cross-check).
    per_cycle_integrator_disagreement: tuple[float, ...]
    max_integrator_disagreement: float
    max_integrator_disagreement_km: float
    drift_floor: float
    agreement_floor: float
    n_cycles_min: int
    converged_each_cycle: bool
    passes_v3_qp: bool
    n_off_grid_samples_per_cycle: int
    n_modes: int
    notes: str


def _propagate_one_cycle_v3(
    torus: QPTorus,
    off_thetas: NDArray[np.float64],
    cycles_done: int,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, float, str, bool]:
    """IAS15 invariance residual + IAS15-vs-DOP853 disagreement for one cycle.

    Returns ``(max_invariance_err_ias15, max_integrator_disagreement, integrator_label,
    ok)``. ``ok`` is False if any propagation raised (the candidate fails V3).
    """
    t_total = (cycles_done + 1) * torus.t_strob
    shift = (cycles_done + 1) * torus.rho
    mu = float(torus.system.mu)
    max_err = 0.0
    max_disagree = 0.0
    label = ""
    for theta in off_thetas:
        u0 = evaluate_invariant_circle(torus.fourier_coeffs, float(theta))
        try:
            state_ias15, label = _ias15_propagate_cr3bp_rotating(u0, t_total, mu)
            arc_dop = cr3bp.propagate(
                torus.system, u0, t_total, with_stm=False, rtol=rtol, atol=atol
            )
        except Exception:
            return float("inf"), float("inf"), label, False
        u_target = evaluate_invariant_circle(torus.fourier_coeffs, float(theta) + shift)
        err = float(np.linalg.norm(state_ias15 - u_target))
        disagree = float(np.linalg.norm(state_ias15 - arc_dop.state_f))
        if err > max_err:
            max_err = err
        if disagree > max_disagree:
            max_disagree = disagree
    return max_err, max_disagree, label, True


def run_v3_qp(
    candidate_id: str,
    torus: QPTorus,
    *,
    n_cycles: int = V2_QP_N_CYCLES_MIN,
    drift_floor: float = V2_QP_DRIFT_FLOOR,
    agreement_floor: float = V3_QP_AGREEMENT_FLOOR,
    n_off_grid_samples_per_cycle: int = V2_QP_OFF_GRID_PER_CYCLE,
    rng_seed_base: int = 0x3DCAFE,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V3VerdictQP:
    """Run the V3 independent-integrator (REBOUND IAS15) QP-torus invariance gauntlet.

    For ``k`` in ``range(n_cycles)``: draw a fresh off-grid angle batch, propagate each
    by ``(k+1)·t_strob`` from the invariant circle under IAS15, and (a) check the result
    lies on the rotated circle ``u(theta + (k+1)·rho)`` (the invariance residual), (b)
    compare the IAS15 endpoint to the DOP853 endpoint (the integrator cross-check).

    ``passes_v3_qp`` iff: all ``n_cycles`` complete, ``n_cycles >= V2_QP_N_CYCLES_MIN``,
    the max IAS15 invariance residual ``<= drift_floor``, AND the max IAS15-vs-DOP853
    disagreement ``<= agreement_floor`` (the V1/V2 signature is integrator-independent).

    A V3_qp PASS does NOT admit to the catalogue.
    """
    if not isinstance(torus, QPTorus):
        raise TypeError(f"torus must be a QPTorus instance; got {type(torus).__name__}")
    if n_cycles < V2_QP_N_CYCLES_MIN:
        raise ValueError(f"V3_qp requires n_cycles >= {V2_QP_N_CYCLES_MIN}; got {n_cycles}")
    if drift_floor <= 0.0:
        raise ValueError(f"drift_floor must be > 0; got {drift_floor}")
    if agreement_floor <= 0.0:
        raise ValueError(f"agreement_floor must be > 0; got {agreement_floor}")
    if n_off_grid_samples_per_cycle < 1:
        raise ValueError(
            f"n_off_grid_samples_per_cycle must be >= 1; got {n_off_grid_samples_per_cycle}"
        )
    if torus.system.l_km <= 0.0 or torus.system.t_s <= 0.0:
        raise ValueError(
            f"invalid CR3BP system for V3-qp km conversion: "
            f"l_km={torus.system.l_km} t_s={torus.system.t_s}"
        )

    n_grid = torus.n_samples
    grid_thetas = 2 * np.pi * np.arange(n_grid) / n_grid

    per_cycle_inv: list[float] = []
    per_cycle_dis: list[float] = []
    integrator = ""
    converged = True
    n_done = 0
    for k in range(n_cycles):
        cycle_seed = rng_seed_base ^ k
        off_thetas = _off_grid_thetas(
            n_off_grid_samples_per_cycle, grid_thetas, rng_seed=cycle_seed
        )
        inv_err, dis_err, label, ok = _propagate_one_cycle_v3(
            torus, off_thetas, cycles_done=k, rtol=rtol, atol=atol
        )
        if label:
            integrator = label
        per_cycle_inv.append(inv_err)
        per_cycle_dis.append(dis_err)
        if not ok:
            converged = False
            break
        n_done += 1

    max_inv = max(per_cycle_inv) if per_cycle_inv else float("inf")
    max_dis = max(per_cycle_dis) if per_cycle_dis else float("inf")
    l_km = float(torus.system.l_km)

    passes = bool(
        converged
        and n_done >= V2_QP_N_CYCLES_MIN
        and max_inv <= drift_floor
        and max_dis <= agreement_floor
    )

    return V3VerdictQP(
        candidate_id=candidate_id,
        integrator=integrator,
        n_cycles_requested=int(n_cycles),
        n_cycles_longitudinal_propagated=int(n_done),
        per_cycle_invariance_residual_ias15=tuple(per_cycle_inv),
        max_invariance_drift_ias15=float(max_inv),
        max_invariance_drift_ias15_km=float(max_inv * l_km),
        per_cycle_integrator_disagreement=tuple(per_cycle_dis),
        max_integrator_disagreement=float(max_dis),
        max_integrator_disagreement_km=float(max_dis * l_km),
        drift_floor=float(drift_floor),
        agreement_floor=float(agreement_floor),
        n_cycles_min=V2_QP_N_CYCLES_MIN,
        converged_each_cycle=bool(converged and n_done == n_cycles),
        passes_v3_qp=passes,
        n_off_grid_samples_per_cycle=int(n_off_grid_samples_per_cycle),
        n_modes=int(torus.n_modes),
        notes=notes,
    )


__all__ = [
    "V3_QP_AGREEMENT_FLOOR",
    "V3VerdictQP",
    "run_v3_qp",
]
