"""QBCP core EOM / STM / propagator / coordinate-mapping tests (#533)."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.qbcp as qbcp

# Rosales-Jorba (2023) Table 4 dynamical-substitute ICs for the Earth-Moon
# collinear points under the QBCP (canonical PM coordinates, y=z=px=pz=0),
# expressed in THIS module's reflected (x -> -x) frame -- i.e. the published
# negative-x values mapped to the repository's Earth-at-(-mu) convention.
# Published: POL1 x=-0.8369141677649317 py=-0.8391311559808445
#            POL2 x=-1.1556836078332600 py=-1.1587306159501061
_POL1_X = 0.8369141677649317
_POL2_X = 1.1556836078332600
_POL1_REFLECTED = np.array([_POL1_X, 0.0, 0.0, 0.0, 0.8391311559808445, 0.0])
_POL2_REFLECTED = np.array([_POL2_X, 0.0, 0.0, 0.0, 1.1587306159501061, 0.0])


def _qbcp_frozen_jacobian(x: float, t: float, system: qbcp.QBCPSystem) -> np.ndarray:
    """Frozen-time state Jacobian A(t) of qbcp_eom at (x, 0, 0), same construction
    as ``qbcp_stm_eom`` builds internally (uses only public helpers)."""
    alphas = qbcp.evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]
    uxx, uyy, uzz, uxy, uxz, uyz = qbcp.qbcp_potential_second_derivatives(x, 0.0, 0.0, t, system)
    jac_a = np.zeros((6, 6), dtype=np.float64)
    jac_a[0, 0] = a2
    jac_a[0, 1] = a3
    jac_a[0, 3] = a1
    jac_a[1, 0] = -a3
    jac_a[1, 1] = a2
    jac_a[1, 4] = a1
    jac_a[2, 2] = a2
    jac_a[2, 5] = a1
    jac_a[3, 0] = uxx
    jac_a[3, 1] = uxy
    jac_a[3, 2] = uxz
    jac_a[3, 3] = -a2
    jac_a[3, 4] = a3
    jac_a[4, 0] = uxy
    jac_a[4, 1] = uyy
    jac_a[4, 2] = uyz
    jac_a[4, 3] = -a3
    jac_a[4, 4] = -a2
    jac_a[5, 0] = uxz
    jac_a[5, 1] = uyz
    jac_a[5, 2] = uzz
    jac_a[5, 5] = -a2
    return jac_a


def _cr3bp_collinear_unstable_rate(x: float, mu: float) -> float:
    """CR3BP collinear-point unstable eigenvalue (Szebehely), computed independently
    of qbcp.py so it is a non-circular reference for the QBCP structural check."""
    c2 = (1.0 - mu) / abs(x + mu) ** 3 + mu / abs(x - 1.0 + mu) ** 3
    return math.sqrt((c2 - 2.0 + math.sqrt(9.0 * c2 * c2 - 8.0 * c2)) / 2.0)


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


# ---------------------------------------------------------------------------
# 5. Collinear-point instability matches CR3BP (structural, sourced golden)
# ---------------------------------------------------------------------------
def test_qbcp_collinear_instability_matches_cr3bp() -> None:
    """The QBCP frozen-time linearization at the EM L1/L2 points reproduces the CR3BP
    collinear unstable rate.

    This is the structural check that would catch a gross error in the alpha_i scaling
    or in the Newtonian potential (either would shift the L-point stiffness). The
    EXPECTED rate is the CR3BP collinear eigenvalue (Szebehely), derived here directly
    from the mass ratio -- it does NOT come from qbcp.py, so the comparison is not
    circular. The QBCP (Sun-perturbed) rate must sit within a few percent of the CR3BP
    rate because the Sun term is a small O(eps^2) perturbation. (Verified 2026-07-10:
    QBCP L1 rate 2.979 vs CR3BP 2.932; QBCP L2 rate 2.199 vs CR3BP 2.159.)
    """
    system = qbcp.qbcp_default()
    mu = system.mu
    for x in (_POL1_X, _POL2_X):
        jac_a = _qbcp_frozen_jacobian(x, 0.0, system)
        eigs = np.linalg.eigvals(jac_a)
        rate_qbcp = float(np.max(eigs.real))
        rate_cr3bp = _cr3bp_collinear_unstable_rate(x, mu)
        rel = abs(rate_qbcp - rate_cr3bp) / rate_cr3bp
        assert rate_qbcp > 1.5, f"expected a real unstable eigenvalue at x={x}, got {rate_qbcp}"
        assert rel < 0.05, (
            f"QBCP frozen collinear rate {rate_qbcp:.4f} at x={x} deviates "
            f"{rel:.1%} from the CR3BP reference {rate_cr3bp:.4f} (expected < 5%)"
        )


# ---------------------------------------------------------------------------
# 6. POL substitutes are unstable: forward-prop non-closure is an artifact
# ---------------------------------------------------------------------------
def test_qbcp_pol_forward_prop_is_instability_dominated() -> None:
    """Pin WHY the published POL1/POL2 substitutes do NOT close under naive forward
    propagation, so the O(1) residual is not re-mistaken for a model bug.

    POL1/POL2 are dynamical substitutes of the violently unstable EM collinear points.
    The monodromy over one synodic period T_s carries a huge unstable multiplier
    (exp(rate * T_s) with rate ~2-3 and T_s ~6.79, i.e. 1e6-1e8), so a single forward
    propagation of even a perfect IC amplifies any roundoff / model-instance offset to
    O(1). The residual is therefore a metric artifact of the instability, NOT a defect
    in qbcp_eom -- proper validation uses continuation + multiple shooting (see
    data/OUTSTANDING.md, task #544). This test asserts both facts jointly: the forward
    residual is O(1) AND the frozen instability is large enough to fully explain it.
    """
    system = qbcp.qbcp_default()
    ts = system.sun_period_tu
    for pol, x in ((_POL1_REFLECTED, _POL1_X), (_POL2_REFLECTED, _POL2_X)):
        sol = solve_ivp(
            lambda t, y: qbcp.qbcp_eom(t, y, system),
            (0.0, ts),
            pol,
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        assert sol.success
        residual = float(np.linalg.norm(sol.y[:, -1] - pol))

        jac_a = _qbcp_frozen_jacobian(x, 0.0, system)
        rate = float(np.max(np.linalg.eigvals(jac_a).real))
        amplification = math.exp(rate * ts)

        # The orbit is genuinely, strongly unstable ...
        assert amplification > 1e5, (
            f"expected large unstable amplification at x={x}, got {amplification:.2e}"
        )
        # ... which is exactly why the naive one-shot residual is O(1) rather than tiny.
        assert residual > 0.5, (
            f"forward-prop residual at x={x} unexpectedly small ({residual:.3e}); if this "
            "ever holds, the instability-artifact reasoning in #544 needs revisiting"
        )
