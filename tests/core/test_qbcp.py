"""QBCP core EOM / STM / propagator / coordinate-mapping tests (#533)."""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.qbcp as qbcp

# ---------------------------------------------------------------------------
# Sample states and time anchors
# ---------------------------------------------------------------------------
_SAMPLE_STATE_PV = np.array([0.5, 0.1, 0.05, 0.02, 0.3, -0.01], dtype=np.float64)
_SAMPLE_STATE_PLANAR_PV = np.array([0.5, 0.1, 0.0, 0.02, 0.3, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# 1. Coordinate Mapping Round-Trip
# ---------------------------------------------------------------------------
def test_coordinate_mapping_roundtrip() -> None:
    """Verify that mapping between PV and PM coordinates is a perfect round-trip."""
    system = qbcp.qbcp_default()
    times = [0.0, 1.25, 4.5, 10.0]

    # Generate some random states around the sample
    np.random.seed(42)
    for t in times:
        for _ in range(5):
            state_pv = _SAMPLE_STATE_PV + np.random.normal(0.0, 0.05, 6)
            state_pm = qbcp.state_pv_to_pm(state_pv, t, system)
            state_pv_back = qbcp.state_pm_to_pv(state_pm, t, system)

            assert np.allclose(state_pv, state_pv_back, rtol=0.0, atol=1e-15), (
                f"PV-to-PM round-trip failed at t={t}. Max delta: "
                f"{np.max(np.abs(state_pv - state_pv_back)):.3e}"
            )


def test_jacobian_inverse_relation() -> None:
    """Verify that transformation_jacobian and transformation_jacobian_inverse are inverses."""
    system = qbcp.qbcp_default()
    times = [0.0, 1.25, 4.5, 10.0]
    for t in times:
        jac_m = qbcp.transformation_jacobian(t, system)
        jac_minv = qbcp.transformation_jacobian_inverse(t, system)

        # Check that jac_m * jac_minv is identity
        prod = jac_m @ jac_minv
        eye = np.eye(6)
        assert np.allclose(prod, eye, rtol=0.0, atol=1e-15), (
            f"M * Minv is not identity at t={t}. Max delta: {np.max(np.abs(prod - eye)):.3e}"
        )


# ---------------------------------------------------------------------------
# 2. Symplectic Conservation
# ---------------------------------------------------------------------------
def test_symplectic_conservation() -> None:
    """Verify that the propagated canonical STM remains symplectic.

    For a symplectic matrix Phi, det(Phi) = 1.0.
    """
    system = qbcp.qbcp_default()
    state_pv0 = _SAMPLE_STATE_PV.copy()
    t_horizon = 2.0

    # Propagate with STM
    times, states_pv = qbcp.propagate_qbcp_pv(state_pv0, (0.0, t_horizon), system, with_stm=True)

    # Check at the final state
    stm_pv = states_pv[-1, 6:].reshape((6, 6))
    # The coordinate transformation jacobian is not necessarily symplectic,
    # but the canonical STM (in PM variables) must be.
    # Let's map STM back to canonical variables to verify det(Phi_PM) = 1.0.
    t0, tf = 0.0, times[-1]
    # stm_pv = M_tf @ stm_pm @ Minv_t0 -> stm_pm = M_tf_inv @ stm_pv @ Minv_t0_inv
    jac_tf_inv = qbcp.transformation_jacobian_inverse(tf, system)
    jac_t0 = qbcp.transformation_jacobian(t0, system)
    stm_pm = jac_tf_inv @ stm_pv @ jac_t0

    det_pm = np.linalg.det(stm_pm)
    assert abs(det_pm - 1.0) < 1e-10, f"det(Phi_PM) = {det_pm:.6f} deviates from 1.0"


# ---------------------------------------------------------------------------
# 3. STM Finite Difference Validation
# ---------------------------------------------------------------------------
def test_stm_finite_difference_consistency() -> None:
    """Verify the analytic QBCP STM EOM against finite differences in PV coordinates."""
    system = qbcp.qbcp_default()
    state_pv0 = _SAMPLE_STATE_PV.copy()
    t_horizon = 0.2

    # Propagate nominal state and get STM
    _times, states_pv = qbcp.propagate_qbcp_pv(state_pv0, (0.0, t_horizon), system, with_stm=True)
    stm_pv = states_pv[-1, 6:].reshape((6, 6))

    # Finite differences
    eps = 1e-6
    fd_jac = np.zeros((6, 6), dtype=np.float64)

    for i in range(6):
        perturb = np.zeros(6, dtype=np.float64)
        perturb[i] = eps

        # Upper perturb
        _, spv_up = qbcp.propagate_qbcp_pv(
            state_pv0 + perturb, (0.0, t_horizon), system, with_stm=False
        )
        # Lower perturb
        _, spv_down = qbcp.propagate_qbcp_pv(
            state_pv0 - perturb, (0.0, t_horizon), system, with_stm=False
        )

        fd_jac[:, i] = (spv_up[-1] - spv_down[-1]) / (2.0 * eps)

    # Check max relative/absolute difference
    max_diff = np.max(np.abs(stm_pv - fd_jac))
    assert max_diff < 1e-5, f"STM does not match finite differences. Max diff: {max_diff:.3e}"


# ---------------------------------------------------------------------------
# 4. Circular limit (structural match to BCR4BP)
# ---------------------------------------------------------------------------
def test_qbcp_circular_limit_eom() -> None:
    """Verify that when QBCP Fourier terms are set to circular limits, EOMs match BCR4BP."""
    # We patch evaluate_alphas to return circular limits:
    # alpha_1 = 1, alpha_2 = 0, alpha_3 = 1, alpha_4 = 0, alpha_5 = 0, alpha_6 = 1
    # alpha_7 = -a_S * cos(theta), alpha_8 = -a_S * sin(theta)
    system = qbcp.qbcp_default()

    # We will compute the EOM using QBCP's EOM at t=0.5 with patched alphas,
    # and compare it to BCR4BP's EOM.
    t = 0.5
    theta = system.theta_sun0 + system.omega_sun_nondim * t
    a_s = system.a_sun_nondim

    # Patched alphas
    alphas_patched = np.zeros(9, dtype=np.float64)
    alphas_patched[1] = 1.0
    alphas_patched[2] = 0.0
    alphas_patched[3] = 1.0
    alphas_patched[4] = system.mu_sun / (a_s**2) * math.cos(theta)
    alphas_patched[5] = system.mu_sun / (a_s**2) * math.sin(theta)
    alphas_patched[6] = 1.0
    # Match BCR4BP's circular Sun position: (a_S * cos(theta), a_S * sin(theta))
    alphas_patched[7] = a_s * math.cos(theta)
    alphas_patched[8] = a_s * math.sin(theta)

    # We temporarily patch qbcp.evaluate_alphas to return alphas_patched
    original_evaluate_alphas = qbcp.evaluate_alphas
    try:
        qbcp.evaluate_alphas = lambda _t, _sys: alphas_patched  # type: ignore[assignment]

        # Test EOM equivalence at sample state
        state_pv = _SAMPLE_STATE_PV.copy()
        state_pm = qbcp.state_pv_to_pm(state_pv, t, system)

        # Calculate d(state_pm)/dt from QBCP
        deriv_pm = qbcp.qbcp_eom(t, state_pm, system)

        # Convert deriv_pm to deriv_pv using M:
        # deriv_pv = M * deriv_pm + dM/dt * state_pm
        # Since alphas are constant, dM/dt = 0, so deriv_pv = M * deriv_pm.
        jac_m = qbcp.transformation_jacobian(t, system)
        deriv_pv_qbcp = jac_m @ deriv_pm

        # Compare to BCR4BP EOM (which uses the circular approximation)
        sys_bcr = bcr4bp.BCR4BPSystem(
            mu=system.mu,
            mu_sun=system.mu_sun,
            a_sun_nondim=system.a_sun_nondim,
            omega_sun_nondim=system.omega_sun_nondim,
            theta_sun0=system.theta_sun0,
        )
        deriv_pv_bcr = bcr4bp.bcr4bp_eom(t, state_pv, sys_bcr)

        # Compare
        assert np.allclose(deriv_pv_qbcp, deriv_pv_bcr, rtol=0.0, atol=1e-12), (
            f"QBCP circular EOM limit deviates from BCR4BP. Max diff: "
            f"{np.max(np.abs(deriv_pv_qbcp - deriv_pv_bcr)):.3e}"
        )

    finally:
        qbcp.evaluate_alphas = original_evaluate_alphas
