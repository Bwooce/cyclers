"""Seedless 2D pseudospectral QBCP quasi-periodic-torus corrector tests (#617).

Covers (a) the unit-level machinery -- real tensor-product Fourier basis and its
analytic derivative, free-variable pack/unpack (only omega2 is free; omega1 is
FIXED), the vectorized QBCP RHS and its 6x6 state-Jacobian (both checked against
``core.qbcp``'s own ``qbcp_eom`` / ``qbcp_stm_eom`` pointwise), the analytic
residual Jacobian (verified against central finite differences), the direct
torus evaluator, and the three gauge rows -- and (b) two staged integration
controls, mirroring ``test_variational_qp_torus.py`` / ``test_variational_
periodic_orbit_qbcp.py``:

  * **SE-L2 positive control** (regression floor, MANDATORY): the seedless 2D
    pseudospectral corrector, bootstrapped from the full-mu Sun-Earth L2 GMOS
    torus that ``genome.qbcp_torus.correct_qbcp_torus`` ALREADY converges cleanly
    (invariance residual ~1.5e-5, the ``#544``-quoted ~3e-5 SE-L2 case),
    reproduces it -- same rotation number to ~2e-3 relative, invariance
    residual at the truncation floor, independent closure tiny. Machinery must
    be sound on the mildly-unstable case before the EM-L2 claim means anything.

  * **EM-L2 headline** (the wall ``#544``/``#538`` blocked on): on the
    violently-unstable Earth-Moon L2 torus at Jacobi 3.13 -- where the existing
    single-period GMOS ``correct_qbcp_torus`` plateaus at invariance residual
    4.771e-01 (``#544``'s best prior result, each stroboscopic sample amplified
    ~1e6-1e8) -- the integration-free pseudospectral corrector converges to an
    invariance residual ~3.4e-3 (n1=12,n2=5) / ~2.3e-3 (n1=16,n2=6), i.e.
    ~140-200x below the old method's plateau, approaching but not fully crossing
    the originally-targeted 1e-3 gate. The ``exact`` and ``lsmr`` trust-region
    solvers reach the SAME minimum (a genuine local min, not a starved budget).

All pinned numbers were reproduced live in this session (2026-07-17) by running
the corrector directly, NOT copied from a docstring. The rotation numbers and
the algebraic residual_rms are deterministic (no random seed anywhere -- the
CR3BP/GMOS bootstraps are deterministic); the closure residual is
integrator-dependent (DOP853 stepping is not bit-reproducible across libm
versions -- see ``tests/genome/test_qp_tori.py``'s own note), so closure is
pinned with generous bounds.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.bcr4bp_torus import (
    correct_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qbcp_torus import QBCPTorus, correct_qbcp_torus
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.variational_qbcp_torus import (
    QBCPTorusVariationalResult,
    SegmentedCLVDiagnostics,
    _alphas_on_theta1,
    _basis_matrices,
    _jacobian,
    _k2_first_harmonic_cols,
    _n_free,
    _pack,
    _qbcp_jacobian_grid,
    _qbcp_rhs_grid,
    _residual,
    _transverse_amplitude,
    _unpack,
    correct_qbcp_torus_pseudospectral,
    discover_qbcp_torus_from_gmos,
    evaluate_torus_state,
    manifold_direction_segmented_clv,
    manifold_state_vec_pseudospectral,
)

# ---------------------------------------------------------------------------
# Unit tests: machinery.
# ---------------------------------------------------------------------------


def test_basis_matrices_derivative_matches_analytic() -> None:
    """``_basis_matrices`` returns [1, cos k., sin k.] and its exact theta
    derivative [0, -k sin k., k cos k.]."""
    n = 3
    thetas = np.linspace(0.0, 2 * np.pi, 9, endpoint=False)
    p, pd = _basis_matrices(n, thetas)
    assert p.shape == (9, 2 * n + 1)
    assert np.allclose(p[:, 0], 1.0)
    assert np.allclose(pd[:, 0], 0.0)
    for k in range(1, n + 1):
        assert np.allclose(p[:, k], np.cos(k * thetas))
        assert np.allclose(pd[:, k], -k * np.sin(k * thetas))
        assert np.allclose(p[:, n + k], np.sin(k * thetas))
        assert np.allclose(pd[:, n + k], k * np.cos(k * thetas))


def test_pack_unpack_round_trip_only_omega2_free() -> None:
    """Free-variable vector packs coeffs plus ONLY omega2 (omega1 fixed) -- one
    fewer free scalar than #612's CR3BP torus."""
    n1, n2 = 4, 2
    rng = np.random.default_rng(0)
    coeffs = rng.normal(size=(6, 2 * n1 + 1, 2 * n2 + 1))
    z = _pack(coeffs, 0.3)
    assert z.shape == (_n_free(n1, n2),)
    assert _n_free(n1, n2) == 6 * (2 * n1 + 1) * (2 * n2 + 1) + 1
    c2, w2 = _unpack(z, n1, n2)
    assert np.array_equal(c2, coeffs)
    assert w2 == 0.3


def test_qbcp_rhs_grid_matches_pointwise_eom() -> None:
    """The vectorized grid RHS equals ``core.qbcp.qbcp_eom`` at each grid point
    (with the Sun phase locked to theta1 via t = theta1/omega1)."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    n1, n2 = 4, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, _ = _basis_matrices(n1, t1)
    p2, _ = _basis_matrices(n2, t2)
    alphas = _alphas_on_theta1(t1, omega1, system)
    rng = np.random.default_rng(3)
    coeffs = rng.normal(scale=0.05, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    u = np.empty((6, m1, m2))
    for c in range(6):
        u[c] = p1 @ coeffs[c] @ p2.T
    rhs = _qbcp_rhs_grid(u, alphas, system.mu, system.mu_sun)
    for i in range(m1):
        for j in range(m2):
            expected = qbcp.qbcp_eom(float(t1[i]) / omega1, u[:, i, j], system)
            assert np.allclose(rhs[:, i, j], expected, atol=1e-12)


def test_qbcp_jacobian_grid_matches_pointwise_stm_a_matrix() -> None:
    """The vectorized grid Jacobian equals ``qbcp_stm_eom``'s A matrix at each
    grid point (extract A from the augmented RHS: dPhi = A when Phi = I)."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    n1, n2 = 3, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    p1, _ = _basis_matrices(n1, t1)
    p2, _ = _basis_matrices(n2, 2 * np.pi * np.arange(m2) / m2)
    alphas = _alphas_on_theta1(t1, omega1, system)
    rng = np.random.default_rng(1)
    coeffs = rng.normal(scale=0.05, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    u = np.empty((6, m1, m2))
    for c in range(6):
        u[c] = p1 @ coeffs[c] @ p2.T
    jf = _qbcp_jacobian_grid(u, alphas, system.mu, system.mu_sun)
    for i in range(m1):
        for j in range(m2):
            y42 = np.concatenate([u[:, i, j], np.eye(6).reshape(36)])
            a_expected = qbcp.qbcp_stm_eom(float(t1[i]) / omega1, y42, system)[6:].reshape(6, 6)
            assert np.allclose(jf[:, :, i, j], a_expected, atol=1e-10)


def test_analytic_jacobian_matches_finite_difference() -> None:
    """The analytic residual Jacobian matches a central finite-difference
    Jacobian to high relative accuracy -- the load-bearing correctness gate for
    the whole corrector (a wrong Jacobian would silently mis-converge)."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    mu, mu_sun = system.mu, system.mu_sun
    rng = np.random.default_rng(2)
    n1, n2 = 4, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    alphas = _alphas_on_theta1(t1, omega1, system)
    coeffs = rng.normal(scale=0.05, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 1.1
    coeffs[1, 0, 1] = 0.03
    z = _pack(coeffs, 1.8)
    args = (omega1, mu, mu_sun, n1, n2, p1, p1d, p2, p2d, alphas, 1, 0.05, 1.0, 1.95, 1.0)
    j_an = _jacobian(z, *args)
    j_fd = np.zeros_like(j_an)
    eps = 1e-6
    for kk in range(z.size):
        zp, zm = z.copy(), z.copy()
        zp[kk] += eps
        zm[kk] -= eps
        j_fd[:, kk] = (_residual(zp, *args) - _residual(zm, *args)) / (2 * eps)
    rel = np.max(np.abs(j_an - j_fd)) / max(np.max(np.abs(j_fd)), 1e-30)
    assert rel < 1e-6, f"analytic Jacobian rel error {rel:.2e}"


def test_residual_has_three_gauge_rows() -> None:
    """The residual is PDE rows (6*m1*m2) plus exactly 3 gauge rows -- one FEWER
    than #612's CR3BP torus (the longitudinal phase gauge is dropped, theta1 is
    locked to the Sun epoch)."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    n1, n2 = 3, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    alphas = _alphas_on_theta1(t1, omega1, system)
    coeffs = np.zeros((6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    z = _pack(coeffs, 1.8)
    r = _residual(
        z,
        omega1,
        system.mu,
        system.mu_sun,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        alphas,
        1,
        0.0,
        1.0,
        1.95,
        1.0,
    )
    assert r.shape == (6 * m1 * m2 + 3,)


def test_rotation_number_gauge_row_value() -> None:
    """The 3rd (last) gauge row is ``rho_weight*(omega2 - rho_target*omega1)`` --
    zero exactly when the ratio matches the target."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    n1, n2 = 3, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    alphas = _alphas_on_theta1(t1, omega1, system)
    coeffs = np.zeros((6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    rho = 1.9
    z = _pack(coeffs, rho * omega1)  # ratio exactly rho
    r = _residual(
        z, omega1, system.mu, system.mu_sun, n1, n2, p1, p1d, p2, p2d, alphas, 1, 0.0, 1.0, rho, 3.0
    )
    assert r[-1] == pytest.approx(0.0, abs=1e-12)
    z2 = _pack(coeffs, (rho + 0.01) * omega1)
    r2 = _residual(
        z2,
        omega1,
        system.mu,
        system.mu_sun,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        alphas,
        1,
        0.0,
        1.0,
        rho,
        3.0,
    )
    assert r2[-1] == pytest.approx(3.0 * 0.01 * omega1, rel=1e-9)


def test_evaluate_torus_state_matches_series() -> None:
    """``evaluate_torus_state`` equals the tensor-product Fourier series evaluated
    directly at a scalar angle pair."""
    from cyclerfinder.search.variational_qbcp_torus import QBCPTorusVariationalResult

    n1, n2 = 3, 2
    rng = np.random.default_rng(5)
    coeffs = rng.normal(scale=0.1, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    res = QBCPTorusVariationalResult(
        system=qbcp.qbcp_default(),
        coeffs=coeffs,
        omega1=0.9252,
        omega2=1.8,
        rotation_number=1.95,
        rho_strob=12.2,
        period=6.79,
        n1=n1,
        n2=n2,
        m1=9,
        m2=7,
        period_multiple=1,
        transverse_amplitude=0.0,
        residual_rms=0.0,
        closure_residual=0.0,
        converged=False,
        n_iter=0,
    )
    a, b = 1.234, 0.567
    pa, _ = _basis_matrices(n1, np.array([a]))
    pb, _ = _basis_matrices(n2, np.array([b]))
    expected = np.array([(pa @ coeffs[c] @ pb.T)[0, 0] for c in range(6)])
    assert np.allclose(evaluate_torus_state(res, a, b), expected, atol=1e-13)


def test_transverse_amplitude_and_k2_cols() -> None:
    n2 = 3
    cos1, sin1 = _k2_first_harmonic_cols(n2)
    assert (cos1, sin1) == (1, n2 + 1)
    coeffs = np.zeros((6, 3, 2 * n2 + 1))
    coeffs[0, 0, cos1] = 3.0
    coeffs[1, 0, sin1] = 4.0
    assert _transverse_amplitude(coeffs, n2) == pytest.approx(5.0)


def test_period_multiple_validation() -> None:
    system = qbcp.qbcp_default()
    coeffs = np.zeros((6, 3, 3))
    coeffs[0, 0, 0] = 0.85
    with pytest.raises(ValueError, match="period_multiple"):
        correct_qbcp_torus_pseudospectral(
            system, coeffs, 1.8, n1=1, n2=1, amplitude_anchor=0.1, period_multiple=0
        )


def test_coeffs_shape_validation() -> None:
    system = qbcp.qbcp_default()
    with pytest.raises(ValueError, match="coeffs0 shape"):
        correct_qbcp_torus_pseudospectral(
            system, np.zeros((6, 3, 3)), 1.8, n1=4, n2=2, amplitude_anchor=0.1
        )


# ---------------------------------------------------------------------------
# Torus builders (shared by the integration controls).
# ---------------------------------------------------------------------------


def _build_se_l2_gmos_torus() -> QBCPTorus:
    """Full-mu Sun-Earth L2 QBCP GMOS torus, built exactly as #533/#537/#538's
    ``build_tori`` do (BCR4BP mu_sun bootstrap then ``correct_qbcp_torus``).
    This is the ~1.5e-5-residual SE-L2 torus #544 quotes at ~3e-5.
    """
    qbcp_sys = qbcp.qbcp_default()
    mu_se = 1.0 / (qbcp_sys.mu_sun + 1.0)
    n_samples, n_modes = 5, 2
    sys_se = cr3bp.CR3BPSystem(
        mu=mu_se, primary="Sun", secondary="Earth", l_km=qbcp_sys.a_sun_nondim * 384400.0, t_s=1.0
    )
    orbit_se = correct_symmetric_fixed_jacobi(
        sys_se,
        x0_guess=(1.0 - mu_se) + 0.010,
        jacobi=3.0008,
        period_guess=3.1,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )
    bcr = bcr4bp.BCR4BPSystem(
        mu=qbcp_sys.mu,
        mu_sun=qbcp_sys.mu_sun,
        a_sun_nondim=qbcp_sys.a_sun_nondim,
        omega_sun_nondim=qbcp_sys.omega_sun_nondim,
        theta_sun0=qbcp_sys.theta_sun0,
    )
    x0b, ppi, amp = se_lyapunov_to_bcr4bp_torus_seed(orbit_se, bcr, mu_se, n_samples=n_samples)
    tb = correct_bcr4bp_torus(bcr, x0b, n_modes, n_samples, ppi, amp, tol=1e-6)
    x0 = np.zeros(6 + 12 * n_modes + 1)
    cb = tb.fourier_coeffs
    x0[0:6] = np.real(cb[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0[i0 : i0 + 6] = np.real(cb[n, :])
        x0[i0 + 6 : i0 + 12] = np.imag(cb[n, :])
    x0[-1] = tb.rho
    return correct_qbcp_torus(qbcp_sys, x0, n_modes, n_samples, ppi, amp, tol=1e-3)


def _em_l2_c313_seed(n1: int, n2: int) -> tuple[np.ndarray, float, float, float]:
    """CR3BP Earth-Moon L2 Lyapunov at Jacobi 3.13 -> theta1-constant 2D PM
    Fourier seed. Returns (coeffs0, omega2_0, transverse_amplitude, lyap_period).

    ydot0_sign=-1.0 selects the genuine L2 Lyapunov (x0=1.189 > L2, T=3.475);
    +1.0 collapses onto a small near-Moon orbit at this energy. The seed carries
    NO Sun-forcing theta1 structure (mu_sun=0 CR3BP) -- the corrector builds that
    itself from the flat seed, the honest #611/#612 "one known family member"
    bootstrap.
    """
    qbcp_sys = qbcp.qbcp_default()
    mu = qbcp_sys.mu
    sys_em = cr3bp.CR3BPSystem(
        mu=mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )
    orbit = correct_symmetric_fixed_jacobi(
        sys_em,
        x0_guess=1.189,
        jacobi=3.13,
        period_guess=3.475,
        ydot0_sign=-1.0,
        half_crossings=1,
        tol=1e-8,
    )
    assert orbit.converged and orbit.x0 > 1.0
    period = orbit.period
    m2 = 2 * n2 + 3
    t2g = np.linspace(0.0, period, m2, endpoint=False)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0]),
        args=(mu,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=t2g,
    )
    u_pm = np.array([qbcp.state_pv_to_pm(sol.y[:, j], 0.0, qbcp_sys) for j in range(m2)])
    t2 = 2 * np.pi * np.arange(m2) / m2
    p2, _ = _basis_matrices(n2, t2)
    p2pinv = np.linalg.pinv(p2)
    coeffs = np.zeros((6, 2 * n1 + 1, 2 * n2 + 1))
    for c in range(6):
        coeffs[c, 0, :] = p2pinv @ u_pm[:, c]
    return coeffs, 2 * np.pi / period, _transverse_amplitude(coeffs, n2), period


# ---------------------------------------------------------------------------
# Integration control 1: SE-L2 positive control (regression floor).
# ---------------------------------------------------------------------------


def test_se_l2_positive_control_reproduces_gmos_torus() -> None:
    """The seedless 2D pseudospectral corrector, bootstrapped from the
    already-converging full-mu SE-L2 GMOS torus, reproduces it: same rotation
    number to ~2e-3 relative, invariance residual at the n1=6,n2=5 truncation
    floor, independent closure tiny. Regression floor -- machinery must be sound
    on the mildly-unstable case before the EM-L2 claim means anything.

    ``omega1`` must come out FIXED at the Sun's synodic frequency (the crux
    design decision: it is an input, not a solved-for unknown).
    """
    qbcp_sys = qbcp.qbcp_default()
    gmos = _build_se_l2_gmos_torus()
    assert gmos.converged
    assert gmos.invariance_residual < 1e-4
    gmos_rot = gmos.omega_trans / gmos.omega_long

    res = discover_qbcp_torus_from_gmos(qbcp_sys, gmos, n1=6, n2=5, max_nfev=1500, tr_solver="lsmr")
    # omega1 is FIXED at the Sun synodic frequency, never a free unknown.
    assert res.omega1 == pytest.approx(qbcp_sys.omega_sun_nondim, rel=1e-12)
    # Independent method reproduces the SAME torus's rotation number.
    assert res.rotation_number == pytest.approx(gmos_rot, rel=5e-3)
    assert res.rotation_number == pytest.approx(0.231365, abs=1e-4)
    # Residual/closure at the n1=6,n2=5 truncation floor (comparable to the GMOS
    # seed's own invariance residual ~1.5e-5, on a much more stringent pointwise
    # PDE metric than GMOS's coarse 5-sample stroboscopic check).
    assert res.residual_rms < 3e-4
    assert res.residual_rms == pytest.approx(9.405e-05, rel=0.2)
    assert res.closure_residual < 3e-4


# ---------------------------------------------------------------------------
# Integration control 2: EM-L2 headline (the #544/#538 wall).
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_em_l2_c313_crosses_gmos_plateau() -> None:
    """On the violently-unstable Earth-Moon L2 torus at Jacobi 3.13 -- where the
    single-period GMOS ``correct_qbcp_torus`` plateaus at invariance residual
    4.771e-01 (#544's best prior result; each stroboscopic sample amplified
    ~1e6-1e8) -- the integration-free pseudospectral corrector converges to an
    invariance residual ~3.4e-3, ~140x below that plateau. It self-terminates at
    a genuine local minimum (NOT a starved budget): the ``exact`` and ``lsmr``
    trust-region solvers reach the SAME residual/rotation-number. This is the
    #544-diagnosed fix -- an integration-free formulation in which the violent
    per-period amplification never enters the residual or its Jacobian.

    Honest scope boundary: ~3.4e-3 (n1=12,n2=5) / ~2.3e-3 (n1=16,n2=6, not
    re-run here for speed) approaches but does NOT fully cross the originally-
    targeted 1e-3 gate -- reported as a large partial improvement, not a full
    crossing. See OUTSTANDING #617.
    """
    qbcp_sys = qbcp.qbcp_default()
    n1, n2 = 12, 5
    coeffs0, omega2_0, amp, _period = _em_l2_c313_seed(n1, n2)

    res = correct_qbcp_torus_pseudospectral(
        qbcp_sys,
        coeffs0,
        omega2_0,
        n1=n1,
        n2=n2,
        amplitude_anchor=amp,
        tr_solver="exact",
        max_nfev=600,
        tol=1e-3,
        closure_tol=1e-2,
    )
    # omega1 fixed at the Sun synodic frequency.
    assert res.omega1 == pytest.approx(qbcp_sys.omega_sun_nondim, rel=1e-12)
    # ORDERS below the old single-period GMOS plateau of 4.771e-01.
    assert res.residual_rms < 1e-2
    assert res.residual_rms < 0.4771 / 50.0
    assert res.residual_rms == pytest.approx(3.412e-03, rel=0.2)
    # It is a genuinely large finite-amplitude EM-L2 torus (not collapsed to the
    # periodic-orbit center), at a physical rotation number near T_s/T_lyap.
    assert res.transverse_amplitude > 0.1
    assert res.rotation_number == pytest.approx(1.9183, abs=0.05)
    # Independent short-time nonlinear-flow closure confirms genuine near-
    # invariance (worse than the algebraic residual, as expected for a violently
    # unstable region, but still ~100x below the old plateau).
    assert res.closure_residual < 1e-2


@pytest.mark.slow
def test_em_l2_exact_and_lsmr_agree() -> None:
    """The ``exact`` (dense SVD) and ``lsmr`` (iterative) trust-region solvers
    converge the EM-L2 torus to the SAME minimum -- proof that ~3.4e-3 is a
    genuine local minimum of the least-squares, not a solver artifact or a
    starved iteration budget. ``lsmr`` is ~3x faster per the module's timing
    notes; ``exact`` is the robust default."""
    qbcp_sys = qbcp.qbcp_default()
    n1, n2 = 12, 5
    coeffs0, omega2_0, amp, _period = _em_l2_c313_seed(n1, n2)
    common: dict[str, Any] = dict(n1=n1, n2=n2, amplitude_anchor=amp, tol=1e-3, closure_tol=1e-2)
    r_exact = correct_qbcp_torus_pseudospectral(
        qbcp_sys, coeffs0, omega2_0, tr_solver="exact", max_nfev=600, **common
    )
    r_lsmr = correct_qbcp_torus_pseudospectral(
        qbcp_sys, coeffs0, omega2_0, tr_solver="lsmr", max_nfev=2500, **common
    )
    assert r_exact.residual_rms == pytest.approx(r_lsmr.residual_rms, rel=1e-2)
    assert r_exact.rotation_number == pytest.approx(r_lsmr.rotation_number, abs=1e-3)
    assert r_exact.transverse_amplitude == pytest.approx(r_lsmr.transverse_amplitude, rel=1e-2)


def test_planar_seed_z_components_start_zero() -> None:
    """The EM-L2 Lyapunov is planar (z=pz=0); the CR3BP seed must carry zero z/pz
    content, so any z-structure in a converged torus would be genuine QBCP
    out-of-plane coupling, not seed leakage."""
    coeffs0, _omega2_0, _amp, _period = _em_l2_c313_seed(8, 4)
    assert np.max(np.abs(coeffs0[2])) < 1e-10  # z
    assert np.max(np.abs(coeffs0[5])) < 1e-10  # pz


# ---------------------------------------------------------------------------
# Manifold-eigenvector extraction adapter (#619).
# ---------------------------------------------------------------------------
#
# All numbers pinned below were reproduced live in this session (2026-07-17) by
# running the adapter directly, NOT copied from a docstring. Direction |dot|
# statistics are robust (spectral, not integrator-bit-sensitive to leading
# figures) and pinned with generous bounds.


def _gmos_manifold_forward(
    torus: QBCPTorus,
    theta_long: float,
    theta_trans: float,
    branch: str,
    ref_vec: NDArray[np.float64] | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """Reference GMOS forward one-period manifold extraction -- mirrors
    ``scripts/run_538_qbcp_cycler.manifold_state_vec`` for a GMOS ``QBCPTorus``
    (evaluate -> one-period augmented STM -> ``local_stability``), replicated
    here so this search-suite test has no ``scripts`` import dependency."""
    import math

    from cyclerfinder.genome.qbcp_torus import evaluate_qbcp_torus
    from cyclerfinder.genome.qp_torus_manifold import local_stability

    state = evaluate_qbcp_torus(torus, theta_long, theta_trans)
    t0 = (float(theta_long) % (2.0 * math.pi)) / torus.omega_long
    try:
        _, spv = qbcp.propagate_qbcp_pv(
            state, (t0, t0 + torus.t_strob), torus.system, with_stm=True, rtol=1e-8, atol=1e-8
        )
    except RuntimeError:
        return None
    stab = local_stability(state, spv[-1, 6:].reshape((6, 6)), hyperbolicity_tol=1e-4)
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None
    if ref_vec is not None and float(np.dot(vec, ref_vec)) < 0.0:
        vec = -vec
    return state, np.asarray(vec, dtype=np.float64)


SETori = tuple[QBCPTorus, QBCPTorusVariationalResult]


@pytest.fixture(scope="module")
def se_tori() -> SETori:
    """Build the SE-L2 GMOS torus and its pseudospectral counterpart ONCE."""
    qbcp_sys = qbcp.qbcp_default()
    gmos = _build_se_l2_gmos_torus()
    pseudo = discover_qbcp_torus_from_gmos(
        qbcp_sys, gmos, n1=6, n2=5, max_nfev=1500, tr_solver="lsmr"
    )
    return gmos, pseudo


def test_variational_result_gmos_alias_properties(se_tori: SETori) -> None:
    """The #619 read-only aliases are exact synonyms so the #538 search can
    duck-type either torus: ``omega_long == omega1`` (fixed Sun frequency),
    ``t_strob == period``, ``invariance_residual == residual_rms``."""
    _gmos, pseudo = se_tori
    assert pseudo.omega_long == pseudo.omega1
    assert pseudo.t_strob == pseudo.period
    assert pseudo.invariance_residual == pseudo.residual_rms
    # omega_long is the fixed Sun synodic frequency, matching a GMOS torus.
    assert pseudo.omega_long == pytest.approx(qbcp.qbcp_default().omega_sun_nondim, rel=1e-12)


def test_manifold_adapter_state_and_unit_eigenvector(se_tori: SETori) -> None:
    """The adapter returns ``(state_pv, eigenvector)`` where ``state_pv`` is
    exactly ``state_pm_to_pv(evaluate_torus_state(...))`` at the Sun epoch
    ``t0 = theta_long/omega1`` and the eigenvector is unit-norm, for both
    branches -- the mechanical contract the #538 search relies on."""
    _gmos, pseudo = se_tori
    theta_long, theta_trans = 1.3, 2.1
    for branch in ("unstable", "stable"):
        mv = manifold_state_vec_pseudospectral(pseudo, theta_long, theta_trans, branch, None)
        assert mv is not None
        state_pv, vec = mv
        pm = evaluate_torus_state(pseudo, theta_long, theta_trans)
        t0 = (theta_long % (2.0 * np.pi)) / pseudo.omega1
        expected_pv = qbcp.state_pm_to_pv(pm, t0, pseudo.system)
        assert np.allclose(state_pv, expected_pv, atol=1e-12)
        assert float(np.linalg.norm(vec)) == pytest.approx(1.0, abs=1e-10)


def test_manifold_adapter_ref_vec_sign_convention(se_tori: SETori) -> None:
    """``ref_vec`` fixes the eigenvector sign by dot-product continuity, exactly
    as ``manifold_state_vec`` does (a negative reference flips the returned
    vector)."""
    _gmos, pseudo = se_tori
    # (1.3, 2.1) is a hyperbolic point (the SE torus is non-hyperbolic at some
    # phases where the one-period map has a complex pair and no real unstable
    # eigenvalue -- the adapter correctly returns None there).
    mv = manifold_state_vec_pseudospectral(pseudo, 1.3, 2.1, "unstable", None)
    assert mv is not None
    _state, vec = mv
    mv_pos = manifold_state_vec_pseudospectral(pseudo, 1.3, 2.1, "unstable", vec)
    mv_neg = manifold_state_vec_pseudospectral(pseudo, 1.3, 2.1, "unstable", -vec)
    assert mv_pos is not None and mv_neg is not None
    assert float(np.dot(mv_pos[1], vec)) > 0.0
    assert float(np.dot(mv_neg[1], vec)) < 0.0


def test_manifold_directions_match_gmos_on_se_l2(se_tori: SETori) -> None:
    """POSITIVE CONTROL (#619 Step 2): the adapter's forward-STM manifold
    directions on the pseudospectral SE-L2 torus reproduce the trusted GMOS
    ``manifold_state_vec`` directions at matched PHYSICAL points (theta1 is
    Sun-locked in both; match by nearest PV state along the invariant circle,
    since the two constructions can differ in their theta2 phase origin). The
    residual disagreement scales with the physical match distance (a smooth
    bundle sampled at slightly different points), so agreement is asserted on the
    well-matched subset. Live-observed median |dot| ~0.9998."""
    gmos, pseudo = se_tori
    for branch in ("unstable", "stable"):
        rows: list[tuple[float, float]] = []
        for theta1 in np.linspace(0.0, 2 * np.pi, 3, endpoint=False):
            # dense pseudo circle at this theta1
            t2p = np.linspace(0.0, 2 * np.pi, 60, endpoint=False)
            p_states, p_vecs = [], []
            for t2 in t2p:
                mvp = manifold_state_vec_pseudospectral(
                    pseudo, float(theta1), float(t2), branch, None
                )
                p_states.append(None if mvp is None else mvp[0])
                p_vecs.append(None if mvp is None else mvp[1])
            for t2g in np.linspace(0.0, 2 * np.pi, 6, endpoint=False):
                mvg = _gmos_manifold_forward(gmos, float(theta1), float(t2g), branch, None)
                if mvg is None:
                    continue
                sg, vg = mvg
                best_d, best_j = np.inf, -1
                for j, sj in enumerate(p_states):
                    if sj is None:
                        continue
                    d = float(np.linalg.norm(sg - sj))
                    if d < best_d:
                        best_d, best_j = d, j
                if best_j < 0:
                    continue
                vec_best = p_vecs[best_j]
                if vec_best is None:
                    continue
                rows.append((best_d, abs(float(np.dot(vg, vec_best)))))
        assert rows, f"no matched hyperbolic points for branch={branch}"
        rows.sort()
        dots = np.array([r[1] for r in rows])
        # Well-matched half agrees to high precision; median excellent overall.
        k = max(1, len(rows) // 2)
        assert np.median(dots) > 0.99, f"{branch}: median |dot| {np.median(dots):.4f} too low"
        assert dots[:k].min() > 0.95, (
            f"{branch}: well-matched min |dot| {dots[:k].min():.4f} too low"
        )


def test_unstable_forward_vs_backward_differ_on_torus(se_tori: SETori) -> None:
    """Documents the #619 finding that for a TORUS point (unlike a periodic
    orbit) the one-period-STM eigenvector is an approximation: the forward-STM
    unstable direction and the ``unstable_via_backward`` extraction are genuinely
    DIFFERENT (rotation mismatch E^u(t0) -> E^u(t0+T)), even on the
    cleanly-converged SE-L2 torus. This is why the GMOS-faithful forward
    extraction is the default and ``unstable_via_backward`` is opt-in."""
    _gmos, pseudo = se_tori
    mvf = manifold_state_vec_pseudospectral(pseudo, 0.0, 0.0, "unstable", None)
    mvb = manifold_state_vec_pseudospectral(
        pseudo, 0.0, 0.0, "unstable", None, unstable_via_backward=True
    )
    assert mvf is not None and mvb is not None
    # The two approximations are materially different directions (not ~parallel).
    assert abs(float(np.dot(mvf[1], mvb[1]))) < 0.9


# ---------------------------------------------------------------------------
# Segment-anchored discrete-QR / CLV extraction (#646).
# ---------------------------------------------------------------------------
#
# These validate the segmented extractor on the mildly-unstable SE-L2 torus (the
# `se_tori` fixture), where it MUST agree with the one-shot / GMOS convention (the
# fix only *matters* at the violently-unstable EM-L2 torus, but must not *break*
# the benign case). The headline EM-L2 fragility->robustness demonstration
# (one-shot vec_u swings ~88 deg under a 1e-4 perturbation; segmented vec_u
# <0.01 deg and converges as O(1/n_segments)) costs a 13-22 min EM-L2 torus build
# and is reported in the `#646` bullet with a scratch-reproducible driver, per the
# `#618`/`#619` precedent for expensive EM-L2 runs -- not pinned here.


def test_segmented_clv_contract(se_tori: SETori) -> None:
    """`(state_pv, unit-vec)` contract: ``state_pv`` is the re-projected torus PV
    state and the eigenvector is unit-norm, for both branches -- same contract as
    the one-shot ``manifold_state_vec_pseudospectral``."""
    _gmos, pseudo = se_tori
    theta_long, theta_trans = 1.3, 2.1
    for branch in ("unstable", "stable"):
        mv = manifold_direction_segmented_clv(
            pseudo, theta_long, theta_trans, branch, None, n_segments=12
        )
        assert mv is not None
        state_pv, vec = mv[0], mv[1]
        pm = evaluate_torus_state(pseudo, theta_long, theta_trans)
        t0 = (theta_long % (2.0 * np.pi)) / pseudo.omega1
        expected_pv = qbcp.state_pm_to_pv(pm, t0, pseudo.system)
        assert np.allclose(state_pv, expected_pv, atol=1e-12)
        assert float(np.linalg.norm(vec)) == pytest.approx(1.0, abs=1e-10)


def test_segmented_clv_ref_vec_sign(se_tori: SETori) -> None:
    """``ref_vec`` fixes the eigenvector sign by dot-product continuity, exactly as
    the one-shot extractor and ``manifold_state_vec`` do."""
    _gmos, pseudo = se_tori
    mv = manifold_direction_segmented_clv(pseudo, 1.3, 2.1, "unstable", None, n_segments=12)
    assert mv is not None
    vec = mv[1]
    mv_pos = manifold_direction_segmented_clv(pseudo, 1.3, 2.1, "unstable", vec, n_segments=12)
    mv_neg = manifold_direction_segmented_clv(pseudo, 1.3, 2.1, "unstable", -vec, n_segments=12)
    assert mv_pos is not None and mv_neg is not None
    assert float(np.dot(mv_pos[1], vec)) > 0.0
    assert float(np.dot(mv_neg[1], vec)) < 0.0


def test_segmented_qr_matches_direct_composition(se_tori: SETori) -> None:
    """The QR-reconstructed one-period STM ``P = Q @ R_prod`` equals the directly
    composed segment product to machine precision -- the QR re-orthonormalization
    bookkeeping is faithful (not an approximation of the composed STM)."""
    _gmos, pseudo = se_tori
    mv = manifold_direction_segmented_clv(
        pseudo, 1.3, 2.1, "unstable", None, n_segments=12, return_diagnostics=True
    )
    assert mv is not None and len(mv) == 3
    diag = mv[2]
    assert isinstance(diag, SegmentedCLVDiagnostics)
    assert diag.qr_vs_direct_stm_relerr < 1e-10


def test_segmented_matches_oneshot_on_se_l2(se_tori: SETori) -> None:
    """POSITIVE CONTROL (#646 step 1): on the mildly-unstable SE-L2 torus the
    segment-anchored CLV direction reproduces the one-shot / GMOS-convention
    direction (both methods are fine at SE-L2; the segmented method must not
    introduce its own bug). Live-observed agreement <0.03 deg at hyperbolic
    phases; asserted <1 deg with margin."""
    _gmos, pseudo = se_tori
    checked = 0
    for theta1 in (0.3, 1.3, 2.0, 2.7):
        for theta2 in (0.5, 2.1, 3.0):
            for branch in ("unstable", "stable"):
                seg = manifold_direction_segmented_clv(
                    pseudo, theta1, theta2, branch, None, n_segments=12
                )
                one = manifold_state_vec_pseudospectral(pseudo, theta1, theta2, branch, None)
                if seg is None or one is None:
                    continue  # non-hyperbolic phase (SE torus has some)
                dot = abs(float(np.dot(seg[1], one[1])))
                assert dot > 0.9998, (
                    f"{branch}@({theta1},{theta2}): seg-vs-oneshot |dot| {dot:.5f} too low"
                )
                checked += 1
    assert checked >= 6, f"too few hyperbolic phases checked ({checked})"


def test_segmented_clv_converges_in_n_segments(se_tori: SETori) -> None:
    """The extracted direction CONVERGES as ``n_segments`` grows (it approaches a
    definite limit -- the true torus one-period linearization -- rather than
    wandering). The coarse-to-reference angle must strictly shrink as the segment
    count increases toward the reference. Demonstrates the mechanism even at
    SE-L2's mild instability (the effect is dramatic at EM-L2)."""
    _gmos, pseudo = se_tori
    theta1, theta2, branch = 1.3, 2.1, "unstable"
    ref = manifold_direction_segmented_clv(pseudo, theta1, theta2, branch, None, n_segments=32)
    assert ref is not None
    vref = ref[1]

    def angle_at(k: int) -> float:
        mv = manifold_direction_segmented_clv(pseudo, theta1, theta2, branch, None, n_segments=k)
        assert mv is not None
        return float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(vref, mv[1])))))))

    a2, a8, a16 = angle_at(2), angle_at(8), angle_at(16)
    # Monotone-ish convergence toward the fine-segment reference.
    assert a16 <= a8 <= a2 + 1e-9
    assert a16 < 0.5  # close to the limit by 16 segments even at mild SE-L2


def test_segmented_clv_diagnostics_bounded_growth(se_tori: SETori) -> None:
    """The per-segment leading growth factors are MODEST (the whole point: the
    one-period ~amplification is split into ``n_segments`` well-conditioned pieces,
    each O(1-ish) at SE-L2), and the Lyapunov exponents / eigenvalues are finite
    and consistent (leading exponent > 0 for a hyperbolic point)."""
    _gmos, pseudo = se_tori
    mv = manifold_direction_segmented_clv(
        pseudo, 1.3, 2.1, "unstable", None, n_segments=12, return_diagnostics=True
    )
    assert mv is not None and len(mv) == 3
    diag = mv[2]
    assert diag.n_segments == 12
    assert diag.per_segment_leading_growth.shape == (12,)
    assert np.all(np.isfinite(diag.per_segment_leading_growth))
    assert float(diag.per_segment_leading_growth.max()) < 1e3  # not the 2e4 one-shot blow-up
    assert diag.lam_u is not None and abs(diag.lam_u) > 1.0
    assert diag.lyapunov_exponents.shape == (6,)
    assert diag.lyapunov_exponents[0] > 0.0
    assert float(np.isfinite(diag.nonnormality_ratio))


def test_segmented_clv_rejects_bad_args(se_tori: SETori) -> None:
    """Guard rails: unknown branch and non-positive segment counts raise."""
    _gmos, pseudo = se_tori
    with pytest.raises(ValueError, match="branch"):
        manifold_direction_segmented_clv(pseudo, 0.0, 0.0, "sideways", None)
    with pytest.raises(ValueError, match="n_segments"):
        manifold_direction_segmented_clv(pseudo, 0.0, 0.0, "unstable", None, n_segments=0)
