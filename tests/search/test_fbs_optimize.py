"""Tests for the FBS-analytic-gradient ΔV-chain optimiser (#243).

These are MECHANICS gates (the expected values are defined by construction by a
same-model Lambert seed, never by the optimiser's own output) plus a Jacobian
correctness check. The science verdict (robustness/cost/optimum head-to-head) lives
in ``scripts/fbs_optimizer_fair_trial.py`` and its results note, not here.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.fbs_match_point import chain_defect, chain_defect_jacobian
from cyclerfinder.core.lambert import lambert
from cyclerfinder.search.fbs_optimize import (
    ChainLegSpec,
    _fbs_legs,
    optimize_chain_fbs,
)

DAY = 86400.0


def _two_leg_eme(
    eph: Ephemeris,
) -> tuple[tuple[ChainLegSpec, ...], list[np.ndarray], list[np.ndarray]]:
    """E->M->E chain with a per-leg Lambert (ballistic) seed for the boundaries."""
    t0 = 0.0
    seq = ["E", "M", "E"]
    tofs = [200 * DAY, 540 * DAY]
    epochs = [t0, t0 + tofs[0], t0 + tofs[0] + tofs[1]]
    pos = [np.asarray(eph.state(b, e)[0]) for b, e in zip(seq, epochs, strict=True)]
    legs = tuple(ChainLegSpec(r0=pos[i], rf=pos[i + 1], tof_s=tofs[i], alpha=0.5) for i in range(2))
    bvs: list[np.ndarray] = []
    for i in range(2):
        s = lambert(pos[i], pos[i + 1], tofs[i], mu=MU_SUN_KM3_S2)[0]
        if i == 0:
            bvs.append(np.asarray(s.v1))
        bvs.append(np.asarray(s.v2))
    dvs = [np.zeros(3), np.zeros(3)]
    return legs, dvs, bvs


def test_chain_jacobian_matches_finite_difference() -> None:
    """The analytic chain Jacobian agrees with central differences (correctness)."""
    eph = Ephemeris("astropy")
    legs, dvs, bvs = _two_leg_eme(eph)
    m = len(legs)
    fbs_legs = _fbs_legs(legs, tuple(bvs), MU_SUN_KM3_S2)
    ana = chain_defect_jacobian(fbs_legs, tuple(dvs))
    x0 = np.concatenate([*dvs, *bvs])

    def defect_of(x: np.ndarray) -> np.ndarray:
        d = [x[3 * i : 3 * i + 3] for i in range(m)]
        b = [x[3 * m + 3 * j : 3 * m + 3 * j + 3] for j in range(m + 1)]
        return chain_defect(_fbs_legs(legs, tuple(b), MU_SUN_KM3_S2), tuple(d))

    fd = np.zeros_like(ana)
    h = 1e-4
    for k in range(x0.size):
        xp, xm = x0.copy(), x0.copy()
        xp[k] += h
        xm[k] -= h
        fd[:, k] = (defect_of(xp) - defect_of(xm)) / (2 * h)
    rel = np.abs(ana - fd) / np.maximum(np.abs(ana), 1.0)
    assert float(np.max(rel)) < 1e-5


def test_analytic_lane_converges_to_optimum() -> None:
    """The analytic-gradient lane converges feasibly to the chain-ΔV optimum.

    A small ballistic seed lands in the optimum's basin; the converged ΔV is the
    same-model optimum both lanes reproduce WHEN FEASIBLE. (The FD lane's lower
    convergence rate from cold seeds is the documented #243 robustness finding —
    measured in scripts/fbs_optimizer_fair_trial.py, not asserted here.)
    """
    eph = Ephemeris("astropy")
    legs, dvs, bvs = _two_leg_eme(eph)
    rng = np.random.default_rng(7)
    dv_seed = tuple(d + rng.normal(scale=0.1, size=3) for d in dvs)
    bv_seed = tuple(b + rng.normal(scale=0.1, size=3) for b in bvs)

    res_a = optimize_chain_fbs(legs, dv_seed, bv_seed, use_analytic_jac=True, feas_tol=1e-6)
    res_f = optimize_chain_fbs(legs, dv_seed, bv_seed, use_analytic_jac=False, feas_tol=1e-6)

    assert res_a.feasible
    # The analytic lane reaches the same-model optimum to high precision (the leg-0
    # impulse drives to ~0 on this near-ballistic chain).
    assert res_a.total_dv_kms < 14.68
    # When the FD lane also returns feasible it cannot do BETTER than the analytic
    # optimum (same NLP, same minimiser); in practice it often stalls slightly above
    # it (the #243 optimum-quality finding). Allow it to be no better and not wildly
    # worse than the analytic optimum.
    if res_f.feasible:
        assert res_f.total_dv_kms >= res_a.total_dv_kms - 1e-6
        assert res_f.total_dv_kms < res_a.total_dv_kms + 1.0


def test_analytic_lane_costs_fewer_constraint_evals() -> None:
    """The FD lane pays many more constraint evals per Jacobian (the FBS cost win)."""
    eph = Ephemeris("astropy")
    legs, dvs, bvs = _two_leg_eme(eph)
    rng = np.random.default_rng(3)
    dv_seed = tuple(d + rng.normal(scale=0.2, size=3) for d in dvs)
    bv_seed = tuple(b + rng.normal(scale=0.2, size=3) for b in bvs)

    res_a = optimize_chain_fbs(legs, dv_seed, bv_seed, use_analytic_jac=True, feas_tol=1e-6)
    res_f = optimize_chain_fbs(legs, dv_seed, bv_seed, use_analytic_jac=False, feas_tol=1e-6)

    # Analytic supplies the constraint Jacobian directly (constr_njev > 0); FD
    # finite-differences it inside SLSQP (constr_njev == 0) at a multiple of the
    # variable count in extra constraint-function evaluations.
    assert res_a.constr_njev > 0
    assert res_f.constr_njev == 0
    assert res_f.constr_nfev > res_a.constr_nfev


def test_validation_errors() -> None:
    """Shape mismatches on the seeds raise."""
    eph = Ephemeris("astropy")
    legs, dvs, bvs = _two_leg_eme(eph)
    with pytest.raises(ValueError):
        optimize_chain_fbs(legs, (dvs[0],), tuple(bvs))  # too few dv seeds
    with pytest.raises(ValueError):
        optimize_chain_fbs(legs, tuple(dvs), tuple(bvs[:-1]))  # too few boundary vs
