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
    a_s = system.a_sun_nondim
    sys_bcr = bcr4bp.BCR4BPSystem(
        mu=system.mu,
        mu_sun=system.mu_sun,
        a_sun_nondim=system.a_sun_nondim,
        omega_sun_nondim=system.omega_sun_nondim,
        theta_sun0=system.theta_sun0,
    )
    # The circular Sun, wherever the bicircular model puts it (its sense is tested separately
    # in tests/core/test_bcr4bp_sun_sense.py; this test is about the structure of the equations).
    sun_x, sun_y, _ = bcr4bp._sun_position(t, sys_bcr)

    # Patched alphas
    alphas_patched = np.zeros(9, dtype=np.float64)
    alphas_patched[1] = 1.0
    alphas_patched[2] = 0.0
    alphas_patched[3] = 1.0
    alphas_patched[4] = system.mu_sun * sun_x / a_s**3
    alphas_patched[5] = system.mu_sun * sun_y / a_s**3
    alphas_patched[6] = 1.0
    alphas_patched[7] = sun_x
    alphas_patched[8] = sun_y

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
def test_qbcp_pol_forward_prop_closes_as_far_as_the_instability_allows() -> None:
    """One forward period from the published POL1/POL2 states.

    History (#544, #892): this test used to assert that the one-period residual is of order 1
    and to explain that by the instability alone. The residual was of order 1 because the model
    was wrong (two series evaluated with exchanged parity, alpha_6 on the Sun term only). With
    the model corrected, the published states agree with the module's own periodic orbits to
    about 2e-8, and the one-period residual is that mismatch times the unstable multiplier:
    measured 2.1e-2 at POL2 (multiplier about 1e6). At POL1 the multiplier is about 1e8, so
    2e-8 is amplified to order 1 and the one-shot residual says nothing; the quarter-period
    test above is the meaningful check there.
    """
    system = qbcp.qbcp_default()
    ts = system.sun_period_tu
    residuals = {}
    for name, pol, x in (("POL1", _POL1_REFLECTED, _POL1_X), ("POL2", _POL2_REFLECTED, _POL2_X)):
        sol = solve_ivp(
            lambda t, y: qbcp.qbcp_eom(t, y, system),
            (0.0, ts),
            pol,
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        assert sol.success
        residuals[name] = float(np.linalg.norm(sol.y[:, -1] - pol))

        jac_a = _qbcp_frozen_jacobian(x, 0.0, system)
        rate = float(np.max(np.linalg.eigvals(jac_a).real))
        assert math.exp(rate * ts) > 1e5, f"expected large unstable amplification at x={x}"

    assert residuals["POL2"] < 0.1


# Andreu (1998), "The Quasi-bicircular Problem", PhD thesis, Table 1.5 (printed page 41), column
# alpha_1 (cosine series), as printed. Filed in the private paper corpus as
# andreu-1998-quasi-bicircular-problem-phd-thesis.pdf.
_ANDREU_1998_TABLE_1_5_ALPHA1 = [
    1.00184160892484e00,
    5.76751772619840e-04,
    1.43877702550763e-02,
    -2.63036297497202e-06,
    1.17627835611893e-04,
    -8.06858139100555e-08,
    9.84324976650129e-07,
    -1.17205439441820e-09,
    8.31190597087959e-09,
    -1.40858423869539e-11,
    7.05071378646684e-11,
    -1.49425963491046e-13,
    5.98241897945123e-13,
]


def test_alpha1_coefficients_match_andreu_1998_table_1_5() -> None:
    """The alpha_1 Fourier table agrees with the printed source to its 15 digits.

    Guards a transcription slip found 2026-10-04 (#884): the j = 5 entry was typed as
    -38.068581391005552e-08, 4.7 times the printed -8.06858139100555e-08.
    """
    code = qbcp._COEFFS_ALPHA1
    assert len(code) == len(_ANDREU_1998_TABLE_1_5_ALPHA1)
    for j, (got, printed) in enumerate(zip(code, _ANDREU_1998_TABLE_1_5_ALPHA1, strict=True)):
        assert math.isclose(got, printed, rel_tol=1e-14), (j, got, printed)


# ---------------------------------------------------------------------------
# #892: the model against published orbits and its own reversing symmetry
# ---------------------------------------------------------------------------

# Rosales, Jorba & Jorba-Cusco (2023), CMDA 135, Table 4: the dynamical substitutes of L1 and L2 at
# t = 0 (y = z = px = pz = 0), here in this module's reflected frame. Jorba-Cusco, Farres & Jorba
# (2018), section 5.1: "the orbits replacing L1 and L2 are small, their maximal distance to the
# corresponding equilibrium point is of order O(10^-6)".
_POL1_PUBLISHED = np.array([0.8369141677649317, 0.0, 0.0, 0.0, 0.8391311559808445, 0.0])
_POL2_PUBLISHED = np.array([1.1556836078332600, 0.0, 0.0, 0.0, 1.1587306159501061, 0.0])


def _max_position_excursion(state0: np.ndarray, system: qbcp.QBCPSystem) -> float:
    quarter = 0.25 * system.sun_period_tu
    worst = 0.0
    for sign in (1.0, -1.0):
        sol = solve_ivp(
            qbcp.qbcp_eom,
            (0.0, sign * quarter),
            state0,
            args=(system,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
            dense_output=True,
        )
        samples = sol.sol(np.linspace(0.0, sign * quarter, 400))
        worst = max(worst, float(np.max(np.hypot(samples[0] - state0[0], samples[1] - state0[1]))))
    return worst


def test_published_substitutes_stay_within_1e5_of_their_point() -> None:
    """From the published POL1 and POL2 states the position stays put for a quarter period
    each way (2.6e-6 and 3.6e-6 measured), as the published description says it must.

    The orbits are unstable by 1e8 and 1e6 per period, so the test stops at a quarter period.
    With alpha_2 and alpha_3 evaluated with exchanged parity the excursion is larger than 1;
    with alpha_6 on the Sun term alone it is 0.17 and 0.10.
    """
    system = qbcp.qbcp_default()
    assert _max_position_excursion(_POL1_PUBLISHED, system) < 1e-5
    assert _max_position_excursion(_POL2_PUBLISHED, system) < 1e-5


def test_reversing_symmetry() -> None:
    """(theta, x, y, z, px, py, pz) -> (-theta, x, -y, z, -px, py, -pz) is a symmetry of the
    Hamiltonian (Jorba-Cusco, Farres & Jorba 2018, section 5). It fails by order 1 if a sine
    series is evaluated as a cosine series or the reverse."""
    system = qbcp.qbcp_default()
    mirror = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
    state0 = np.array([0.83, 0.02, 0.01, 0.01, 0.82, 0.02])
    fwd = solve_ivp(
        qbcp.qbcp_eom, (0.0, 2.0), state0, args=(system,), method="DOP853", rtol=1e-13, atol=1e-13
    ).y[:, -1]
    bwd = solve_ivp(
        qbcp.qbcp_eom,
        (0.0, -2.0),
        mirror * state0,
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    ).y[:, -1]
    assert np.linalg.norm(mirror * fwd - bwd) < 1e-9


def test_alpha_parities_match_andreu_table_1_5() -> None:
    """alpha_2, alpha_5, alpha_8 vanish at theta = 0 (sine series); the others do not."""
    alphas = qbcp.evaluate_alphas(0.0, qbcp.qbcp_default())
    for k in (2, 5, 8):
        assert alphas[k] == 0.0
    assert math.isclose(alphas[3], sum(qbcp._COEFFS_ALPHA3), rel_tol=1e-14)
    assert alphas[3] > 1.019
