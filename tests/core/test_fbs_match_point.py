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
    BodyKinematics,
    FbsLeg,
    match_point_defect,
    match_point_defect_epoch_column,
    match_point_defect_jacobian,
    match_point_defect_vinf_jacobian,
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


def _fd_phase_columns(leg: FbsLeg, dv: np.ndarray) -> np.ndarray:
    """Central-difference 6x2 phase columns of the defect w.r.t. ``[tof_s; alpha]``."""
    base_dv = np.asarray(dv, dtype=np.float64)
    cols: list[np.ndarray] = []

    # tof_s column
    h = _FD_STEP * leg.tof_s
    lp = FbsLeg(
        r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s + h, alpha=leg.alpha, mu=leg.mu
    )
    lm = FbsLeg(
        r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s - h, alpha=leg.alpha, mu=leg.mu
    )
    cols.append((match_point_defect(lp, base_dv) - match_point_defect(lm, base_dv)) / (2.0 * h))

    # alpha column
    h = _FD_STEP * max(leg.alpha, 1.0 - leg.alpha)
    lp = FbsLeg(
        r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha + h, mu=leg.mu
    )
    lm = FbsLeg(
        r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha - h, mu=leg.mu
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


def test_phase_jacobian_fd_consistency_random() -> None:
    """Analytic phase-TOF columns match central differences (consistency test).

    Extends the Phase 1 FD-vs-analytic discipline (Pitkin Eqs. 43-44 / Ellison
    Eq. 58) to the ``[∂c/∂tof_s | ∂c/∂alpha]`` columns of the 6x11 Jacobian.
    """
    rng = np.random.default_rng(20260614)
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
        rf, vf = coe_to_rv(
            float(rng.uniform(0.6, 3.0)) * AU_KM,
            float(rng.uniform(0.0, 0.6)),
            float(rng.uniform(0.0, 2.0 * pi)),
            arg_peri_rad=float(rng.uniform(0.0, 2.0 * pi)),
        )
        leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)
        j_an = match_point_defect_jacobian(leg, dv, include_phase=True)
        assert j_an.shape == (6, 11)
        # The first 9 columns must be byte-identical to the default Jacobian.
        np.testing.assert_array_equal(j_an[:, 0:9], match_point_defect_jacobian(leg, dv))
        j_phase_fd = _fd_phase_columns(leg, dv)
        assert _block_rel_err(j_an[:, 9:11], j_phase_fd) < _FD_RTOL, (a_km, e, alpha, tof_s)


# --- Phase 3: boundary v-infinity + moving-endpoint epoch partials -----------


def _accel(r: np.ndarray, mu: float = MU_SUN_KM3_S2) -> np.ndarray:
    """Two-body acceleration ``-mu r / |r|^3`` (km/s^2)."""
    r = np.asarray(r, dtype=np.float64)
    rn = float(np.linalg.norm(r))
    return -mu * r / rn**3


def _moving_body_leg(
    r_b0_ref: np.ndarray,
    v_b0_ref: np.ndarray,
    r_bf_ref: np.ndarray,
    v_bf_ref: np.ndarray,
    vinf_out: np.ndarray,
    vinf_in: np.ndarray,
    t0: float,
    tof_s: float,
    alpha: float,
) -> FbsLeg:
    """Build a leg whose endpoints ride Keplerian bodies at epoch ``t0`` / ``t0+tof``.

    Body states are propagated from a reference epoch (=0) so the leg is a true
    function of ``t0``: the FD epoch test perturbs ``t0`` and rebuilds the leg.
    """
    r_b0, v_b0 = propagate(r_b0_ref, v_b0_ref, t0)
    r_bf, v_bf = propagate(r_bf_ref, v_bf_ref, t0 + tof_s)
    r0 = r_b0
    v0 = v_b0 + vinf_out
    rf = r_bf
    vf = v_bf + vinf_in
    return FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)


def test_vinf_columns_equal_boundary_velocity_columns() -> None:
    """``∂c/∂v∞`` equals the v0/vf columns (Eq. 57: v∞ map is additive, consistency)."""
    r0 = np.array([AU_KM, 0.0, 0.0])
    v0 = np.array([0.0, _V_CIRC_1AU, 0.02 * _V_CIRC_1AU])
    rf = np.array([0.1 * AU_KM, 1.4 * AU_KM, -0.05 * AU_KM])
    vf = np.array([-0.7 * _V_CIRC_1AU, 0.3 * _V_CIRC_1AU, 0.0])
    leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=200.0 * SECONDS_PER_DAY, alpha=0.45)
    dv = np.array([0.1, 0.2, -0.05])
    j_full = match_point_defect_jacobian(leg, dv)
    j_vinf = match_point_defect_vinf_jacobian(leg, dv)
    assert j_vinf.shape == (6, 6)
    np.testing.assert_array_equal(j_vinf, j_full[:, 3:9])


def test_epoch_column_fd_consistency_random() -> None:
    """Analytic ``∂c/∂t0`` matches central differences over moving bodies (consistency).

    Ellison Eqs. 59-61: with TOF fixed, both boundary bodies translate in epoch.
    The leg endpoints are Keplerian bodies propagated from a reference epoch, so
    perturbing ``t0`` rebuilds a genuinely epoch-dependent leg; the analytic
    ``Φ_bwd·[v_bf;a_bf] - Φ_fwd·[v_b0;a_b0]`` is checked against the FD of the
    defect.
    """
    rng = np.random.default_rng(20260615)
    for _ in range(6):
        a0 = float(rng.uniform(0.8, 1.2)) * AU_KM
        r_b0_ref, v_b0_ref = coe_to_rv(
            a0,
            float(rng.uniform(0.0, 0.2)),
            float(rng.uniform(0.0, 2.0 * pi)),
            arg_peri_rad=float(rng.uniform(0.0, 2.0 * pi)),
        )
        af = float(rng.uniform(1.3, 1.8)) * AU_KM
        r_bf_ref, v_bf_ref = coe_to_rv(
            af,
            float(rng.uniform(0.0, 0.2)),
            float(rng.uniform(0.0, 2.0 * pi)),
            arg_peri_rad=float(rng.uniform(0.0, 2.0 * pi)),
        )
        vinf_out = rng.uniform(-3.0, 3.0, size=3)
        vinf_in = rng.uniform(-3.0, 3.0, size=3)
        t0 = float(rng.uniform(0.0, 300.0)) * SECONDS_PER_DAY
        tof_s = float(rng.uniform(150.0, 350.0)) * SECONDS_PER_DAY
        alpha = float(rng.uniform(0.3, 0.7))
        dv = rng.uniform(-0.3, 0.3, size=3)

        leg = _moving_body_leg(
            r_b0_ref, v_b0_ref, r_bf_ref, v_bf_ref, vinf_out, vinf_in, t0, tof_s, alpha
        )
        # Body kinematics at the two epochs (v = dr/dt, a = dv/dt of the body).
        r_b0, v_b0 = propagate(r_b0_ref, v_b0_ref, t0)
        r_bf, v_bf = propagate(r_bf_ref, v_bf_ref, t0 + tof_s)
        body0 = BodyKinematics(v=v_b0, a=_accel(r_b0))
        bodyf = BodyKinematics(v=v_bf, a=_accel(r_bf))
        col_an = match_point_defect_epoch_column(leg, dv, body0=body0, bodyf=bodyf)

        h = _FD_STEP * max(t0, tof_s)
        lp = _moving_body_leg(
            r_b0_ref, v_b0_ref, r_bf_ref, v_bf_ref, vinf_out, vinf_in, t0 + h, tof_s, alpha
        )
        lm = _moving_body_leg(
            r_b0_ref, v_b0_ref, r_bf_ref, v_bf_ref, vinf_out, vinf_in, t0 - h, tof_s, alpha
        )
        col_fd = (match_point_defect(lp, dv) - match_point_defect(lm, dv)) / (2.0 * h)

        assert col_an.shape == (6,)
        err = _block_rel_err(col_an.reshape(6, 1), col_fd.reshape(6, 1))
        assert err < _FD_RTOL, (a0, af, t0, tof_s, alpha, err)
