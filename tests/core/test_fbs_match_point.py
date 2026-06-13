"""Tests for cyclerfinder.core.fbs_match_point (Ellison 2018 FBS, massless DSM leg).

CONSISTENCY-TEST PATTERN (same discipline as ``tests/core/test_kepler_stm.py``
and ``tests/nbody`` for ``flyby_gradients``): Ellison 2018 prints no unit-level
numeric gradient (mining note ``docs/notes/2026-06-10-ellison-2018-analytic-
gradients-mining.md`` §6), so the analytic match-point-defect Jacobian is
validated against central differences of the defect itself, never a sourced
golden. The zero-defect test below uses a forward-CONSTRUCTED self-consistent leg
(no external numbers), not a published value.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.fbs_match_point import (
    FbsLeg,
    match_point_defect,
    match_point_defect_jacobian,
)
from cyclerfinder.core.kepler import coe_to_rv, propagate

_V_CIRC_1AU = sqrt(MU_SUN_KM3_S2 / AU_KM)

# Central-difference step + tolerance, mirroring tests/core/test_kepler_stm.py
# (cbrt(eps)-scaled step; truncation ~ h^2 and roundoff ~ eps/h both ~ few e-11).
_FD_STEP = 1.0e-6
_FD_RTOL = 1.0e-6


def _self_consistent_leg(
    r0: np.ndarray, v0: np.ndarray, tof_s: float, alpha: float, dv: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Build ``(rf, vf)`` so the leg with this ``(Δv, alpha)`` has a ~zero defect."""
    t_burn = alpha * tof_s
    r_b, v_b = propagate(r0, v0, t_burn)
    v_b_post = v_b + dv
    rf, vf = propagate(r_b, v_b_post, tof_s - t_burn)
    return rf, vf


def _fd_defect_jacobian(leg: FbsLeg, dv: np.ndarray) -> np.ndarray:
    """Central-difference 6x9 Jacobian of the defect w.r.t. ``[Δv; v0; vf]``."""
    base_dv = np.asarray(dv, dtype=np.float64)
    cols: list[np.ndarray] = []

    # Δv columns
    for j in range(3):
        h = _FD_STEP * max(abs(base_dv[j]), 1.0)
        e = np.zeros(3)
        e[j] = h
        cp = match_point_defect(leg, base_dv + e)
        cm = match_point_defect(leg, base_dv - e)
        cols.append((cp - cm) / (2.0 * h))

    # v0 columns
    v0 = np.asarray(leg.v0, dtype=np.float64)
    n_v0 = float(np.linalg.norm(v0))
    for j in range(3):
        h = _FD_STEP * max(abs(v0[j]), n_v0)
        e = np.zeros(3)
        e[j] = h
        lp = FbsLeg(
            r0=leg.r0, v0=v0 + e, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu
        )
        lm = FbsLeg(
            r0=leg.r0, v0=v0 - e, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu
        )
        cols.append((match_point_defect(lp, base_dv) - match_point_defect(lm, base_dv)) / (2.0 * h))

    # vf columns
    vf = np.asarray(leg.vf, dtype=np.float64)
    n_vf = float(np.linalg.norm(vf))
    for j in range(3):
        h = _FD_STEP * max(abs(vf[j]), n_vf)
        e = np.zeros(3)
        e[j] = h
        lp = FbsLeg(
            r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=vf + e, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu
        )
        lm = FbsLeg(
            r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=vf - e, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu
        )
        cols.append((match_point_defect(lp, base_dv) - match_point_defect(lm, base_dv)) / (2.0 * h))

    return np.column_stack(cols)


def _block_rel_err(a: np.ndarray, b: np.ndarray) -> float:
    """Block-wise (r-rows / v-rows) Frobenius-relative error.

    Block-wise (not element-wise) because the position and velocity rows carry
    different units and individual entries pass through zero.
    """
    errs: list[float] = []
    for r0 in (0, 3):
        aa, bb = a[r0 : r0 + 3, :], b[r0 : r0 + 3, :]
        denom = float(np.linalg.norm(bb))
        if denom > 0.0:
            errs.append(float(np.linalg.norm(aa - bb)) / denom)
    return max(errs)


def test_defect_zero_on_self_consistent_leg() -> None:
    """A forward-constructed self-consistent leg has a ~zero match-point defect."""
    r0 = np.array([AU_KM, 0.0, 0.0])
    v0 = np.array([0.0, _V_CIRC_1AU, 0.03 * _V_CIRC_1AU])
    tof_s = 250.0 * SECONDS_PER_DAY
    alpha = 0.4
    dv = np.array([0.2, -0.1, 0.05])
    rf, vf = _self_consistent_leg(r0, v0, tof_s, alpha, dv)
    leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)
    c = match_point_defect(leg, dv)
    assert c.shape == (6,)
    assert float(np.linalg.norm(c[0:3])) / AU_KM < 1.0e-10
    assert float(np.linalg.norm(c[3:6])) / _V_CIRC_1AU < 1.0e-10


def test_jacobian_fd_consistency_random() -> None:
    """Analytic defect Jacobian matches central differences (consistency test)."""
    rng = np.random.default_rng(20260613)
    for _ in range(6):
        a_km = float(rng.uniform(0.6, 3.0)) * AU_KM
        e = float(rng.uniform(0.0, 0.6))
        nu = float(rng.uniform(0.0, 2.0 * pi))
        argp = float(rng.uniform(0.0, 2.0 * pi))
        r0, v0 = coe_to_rv(a_km, e, nu, arg_peri_rad=argp)
        period = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
        tof_s = float(rng.uniform(0.15, 0.6)) * period
        alpha = float(rng.uniform(0.25, 0.75))
        dv = rng.uniform(-0.3, 0.3, size=3)
        # Arbitrary (NOT self-consistent) right boundary -> generic non-zero defect.
        rf, vf = coe_to_rv(
            float(rng.uniform(0.6, 3.0)) * AU_KM,
            float(rng.uniform(0.0, 0.6)),
            float(rng.uniform(0.0, 2.0 * pi)),
            arg_peri_rad=float(rng.uniform(0.0, 2.0 * pi)),
        )
        leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)
        j_an = match_point_defect_jacobian(leg, dv)
        j_fd = _fd_defect_jacobian(leg, dv)
        assert j_an.shape == (6, 9)
        assert _block_rel_err(j_an, j_fd) < _FD_RTOL, (a_km, e, nu, alpha, tof_s)
