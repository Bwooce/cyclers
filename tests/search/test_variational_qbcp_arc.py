"""Integration-free Chebyshev-collocation transfer-arc corrector tests (#620).

Covers (a) the unit machinery -- the Chebyshev-Gauss-Lobatto differentiation
matrix vs an analytic derivative, the QBCP state-Jacobian vs ``core.qbcp``'s own
``qbcp_stm_eom``, spectral convergence of the single-arc ODE defect, pack/unpack,
the residual-shape/velocity-match guard, and the analytic connection-Jacobian
node block vs central finite differences -- and (b) the physics positive controls:

  * **Full-corrector positive control (periodic-orbit split):** a genuine `#611`
    QBCP periodic orbit is split into two arcs whose endpoints are fixed physical
    states (:class:`ConstantTorusPoint`); the full two-arc corrector drives the
    combined residual to MACHINE precision and the INDEPENDENT Radau closure check
    confirms it (loop defect limited only by the input orbit's own ~3e-6 closure,
    not the corrector). This validates the entire residual + analytic-Jacobian +
    independent-check path end-to-end against a known closed trajectory.

  * **Ghost rejection (the mandatory-check demonstration):** from a poor (linear-
    interpolation) seed the SAME split converges to a ghost -- a degree-``K``
    polynomial with MACHINE-ZERO nodal residual that is NOT a real trajectory
    (the independent Radau loop defect is O(1), ~20,000 km arrival error). This
    is the concrete form of this project's "it closed is the danger signal"
    discipline: the algebraic residual alone is never sufficient here.

  * **Real-torus endpoint control (SE-L2 GMOS):** ``torus_point_pv`` reproduces
    the trusted ``evaluate_qbcp_torus`` on a genuinely-built Sun-Earth L2 GMOS
    torus, and a collocation arc between two of its points reproduces the real
    on-torus segment (independent Radau tiny).

All pinned numbers were reproduced live in this session (2026-07-17) by running
the corrector directly, NOT copied from a docstring. Algebraic residuals are
deterministic; the Radau closure defects are integrator-dependent and pinned with
generous bounds (see ``test_variational_qbcp_torus.py``'s identical note).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.bcr4bp_torus import (
    correct_bcr4bp_torus,
    se_lyapunov_to_bcr4bp_torus_seed,
)
from cyclerfinder.genome.qbcp_torus import QBCPTorus, correct_qbcp_torus, evaluate_qbcp_torus
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.variational_periodic_orbit_qbcp import discover_qbcp_periodic_orbit
from cyclerfinder.search.variational_qbcp_arc import (
    ConstantTorusPoint,
    QBCPTorusVariationalResult,
    _qbcp_state_jacobian,
    arc_node_times,
    arc_ode_defects,
    build_connection_jacobian,
    build_connection_residual,
    cheb_diff_matrix,
    connection_shape,
    correct_qbcp_arc_connection,
    n_residuals,
    n_unknowns,
    pack_unknowns,
    torus_point_pv,
    unpack_unknowns,
)

# ---------------------------------------------------------------------------
# Unit tests: machinery.
# ---------------------------------------------------------------------------


def test_cheb_diff_matrix_derivative_matches_analytic() -> None:
    """The CGL differentiation matrix differentiates a smooth function to
    spectral accuracy, and its nodes are the Gauss-Lobatto points x_j =
    cos(pi j / K) with endpoints +1 and -1."""
    x, d = cheb_diff_matrix(24)
    assert x[0] == pytest.approx(1.0)
    assert x[-1] == pytest.approx(-1.0)
    f = np.exp(x) * np.sin(3.0 * x)
    fp_exact = np.exp(x) * np.sin(3.0 * x) + 3.0 * np.exp(x) * np.cos(3.0 * x)
    assert np.max(np.abs(d @ f - fp_exact)) < 1e-11


def test_cheb_diff_matrix_rejects_order_zero() -> None:
    with pytest.raises(ValueError, match="order must be >= 1"):
        cheb_diff_matrix(0)


def test_state_jacobian_matches_stm_eom() -> None:
    """``_qbcp_state_jacobian`` equals the A matrix ``qbcp_stm_eom`` builds
    (extract A from the augmented RHS with Phi = I)."""
    system = qbcp.qbcp_default()
    rng = np.random.default_rng(0)
    for _ in range(4):
        t = float(rng.uniform(0.0, 7.0))
        u = np.array([0.6, 0.1, 0.02, 0.01, 0.5, -0.01]) + 0.05 * rng.standard_normal(6)
        a = _qbcp_state_jacobian(t, u, system)
        y42 = np.concatenate([u, np.eye(6).reshape(36)])
        a_ref = qbcp.qbcp_stm_eom(t, y42, system)[6:].reshape(6, 6)
        assert np.max(np.abs(a - a_ref)) < 1e-12


def test_arc_ode_defects_spectral_convergence() -> None:
    """A genuine integrated QBCP trajectory segment has an ODE collocation defect
    that decreases GEOMETRICALLY as the node count rises (the spectral-convergence
    signature) -- confirming the differentiation-matrix time-scaling and RHS
    wiring are correct. (A real solve, with the interior nodes free, drives this
    to machine zero; see the positive-control tests.)"""
    system = qbcp.qbcp_default()
    t_start, tau = 0.3, 0.6
    state0_pm = qbcp.state_pv_to_pm(np.array([0.5, 0.1, 0.0, 0.02, 0.55, 0.0]), t_start, system)
    sol = solve_ivp(
        lambda t, y: qbcp.qbcp_eom(t, y, system),
        (t_start, t_start + tau),
        state0_pm,
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert sol.success
    # dense_output=True above (+ the success check just above) guarantees a
    # dense interpolant -- see the same pattern's comment in
    # src/cyclerfinder/search/perimoon_passage.py.
    assert sol.sol is not None
    prev = np.inf
    for order in (12, 20, 28):
        nodes, diff = cheb_diff_matrix(order)
        t_j = arc_node_times(t_start, tau, nodes)
        node_states = np.array([sol.sol(t) for t in t_j])
        defects = arc_ode_defects(node_states, tau, t_start, diff, nodes, system)
        rms = float(np.sqrt(np.mean(defects[1:order] ** 2)))
        assert rms < 0.5 * prev  # geometric decrease
        prev = rms
    assert prev < 1e-2


def test_pack_unpack_round_trip() -> None:
    order = 6
    rng = np.random.default_rng(1)
    nf = rng.normal(size=(order + 1, 6))
    nr = rng.normal(size=(order + 1, 6))
    th = [rng.normal(size=2) for _ in range(4)]
    z = pack_unknowns(nf, nr, 1.1, 2.2, *th)
    assert z.shape == (n_unknowns(order),)
    nf2, nr2, tf, tr, t0, t1, t2, t3 = unpack_unknowns(z, order)
    assert np.array_equal(nf2, nf)
    assert np.array_equal(nr2, nr)
    assert (tf, tr) == (1.1, 2.2)
    for a, b in zip((t0, t1, t2, t3), th, strict=True):
        assert np.array_equal(a, b)


def test_connection_shape_velocity_match_present() -> None:
    """The raw two-arc count is under-determined by 4 (genuine torus solution-
    family freedom), but the load-bearing anti-#537 guarantee -- the full velocity
    half of every endpoint match is present -- holds. Interior-only ODE gives
    n_residuals = 2*(6*(K-1) + 12 + 1) + 4."""
    for order in (10, 24, 40):
        nu, nr, vel = connection_shape(order)
        assert nu == 12 * (order + 1) + 10
        assert nr == 2 * (6 * (order - 1) + 13) + 4
        assert nu - nr == 4  # under-determined by exactly 4 (family freedom)
        assert vel is True
        assert nr == n_residuals(order) and nu == n_unknowns(order)


def test_connection_jacobian_node_block_matches_finite_difference() -> None:
    """The analytic node-value columns of the connection Jacobian (the large,
    well-conditioned bulk) match a central finite-difference Jacobian to high
    relative accuracy -- the load-bearing correctness gate. The ten special
    columns (tau_f, tau_r, 8 phases) are FD-by-construction and excluded here."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    # Cheap constant-coeff pseudospectral tori (no build) as endpoint providers.
    c0 = np.zeros((6, 3, 3))
    c0[:, 0, 0] = np.array([0.6, 0.05, 0.0, 0.01, 0.5, 0.0])
    c1 = np.zeros((6, 3, 3))
    c1[:, 0, 0] = np.array([0.9, -0.03, 0.0, -0.01, 0.4, 0.0])

    def _tor(coeffs: np.ndarray) -> QBCPTorusVariationalResult:
        return QBCPTorusVariationalResult(
            system=system,
            coeffs=coeffs,
            omega1=omega1,
            omega2=1.8,
            rotation_number=1.9,
            rho_strob=12.0,
            period=2 * np.pi / omega1,
            n1=1,
            n2=1,
            m1=5,
            m2=5,
            period_multiple=1,
            transverse_amplitude=0.0,
            residual_rms=0.0,
            closure_residual=0.0,
            converged=True,
            n_iter=0,
        )

    order = 8
    res = build_connection_residual(_tor(c0), _tor(c1), order, omega1)
    jacf = build_connection_jacobian(_tor(c0), _tor(c1), order, omega1, res)
    rng = np.random.default_rng(2)
    z = rng.normal(scale=0.1, size=n_unknowns(order))
    nb = 6 * (order + 1)
    z[2 * nb], z[2 * nb + 1] = 1.2, 1.5
    z[2 * nb + 2 :] = rng.uniform(-0.5, 0.5, 8)
    j_an = jacf(z)
    eps = 1e-7
    j_fd = np.zeros_like(j_an)
    for k in range(z.size):
        zp, zm = z.copy(), z.copy()
        zp[k] += eps
        zm[k] -= eps
        j_fd[:, k] = (res(zp) - res(zm)) / (2 * eps)
    node_cols = slice(0, 2 * nb)
    rel = np.max(np.abs(j_an[:, node_cols] - j_fd[:, node_cols])) / max(
        np.max(np.abs(j_fd[:, node_cols])), 1e-30
    )
    assert rel < 1e-7, f"analytic node-block Jacobian rel error {rel:.2e}"


def test_torus_point_pv_constant_and_pseudospectral() -> None:
    """``torus_point_pv`` returns the fixed PV for a ConstantTorusPoint and the
    ``state_pm_to_pv(evaluate_torus_state(...))`` value for a pseudospectral
    torus."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    pv = np.array([0.7, 0.2, 0.0, 0.03, 0.4, 0.0])
    const = ConstantTorusPoint(pv=pv, system=system, omega_long=omega1)
    assert np.array_equal(torus_point_pv(const, 1.3, 0.4, omega1), pv)
    c = np.zeros((6, 3, 3))
    c[:, 0, 0] = np.array([0.6, 0.05, 0.0, 0.01, 0.5, 0.0])
    tor = QBCPTorusVariationalResult(
        system=system,
        coeffs=c,
        omega1=omega1,
        omega2=1.8,
        rotation_number=1.9,
        rho_strob=12.0,
        period=2 * np.pi / omega1,
        n1=1,
        n2=1,
        m1=5,
        m2=5,
        period_multiple=1,
        transverse_amplitude=0.0,
        residual_rms=0.0,
        closure_residual=0.0,
        converged=True,
        n_iter=0,
    )
    got = torus_point_pv(tor, 1.3, 0.4, omega1)
    expect = qbcp.state_pm_to_pv(np.array([0.6, 0.05, 0.0, 0.01, 0.5, 0.0]), 1.3 / omega1, system)
    assert np.allclose(got, expect, atol=1e-13)


# ---------------------------------------------------------------------------
# Physics positive controls (fixtures build a #611 orbit and an SE-L2 torus).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def periodic_split() -> tuple[
    qbcp.QBCPSystem, ConstantTorusPoint, ConstantTorusPoint, dict[str, Any]
]:
    """A genuine #611 QBCP periodic orbit split into two arcs, with the true
    split-trajectory node states available for a truth seed."""
    system = qbcp.qbcp_default()
    omega1 = system.omega_sun_nondim
    orb = discover_qbcp_periodic_orbit(
        system, n_harmonics=32, n_restarts=4, rng=np.random.default_rng(0), tol=1e-6, max_nfev=1500
    )
    assert orb.residual_rms < 1e-6
    sol = solve_ivp(
        lambda t, y: qbcp.qbcp_eom(t, y, system),
        (0.0, orb.period),
        orb.state0_pm,
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert sol.success
    # dense_output=True above (+ the success check just above) guarantees a
    # dense interpolant -- see the same pattern's comment in
    # src/cyclerfinder/search/perimoon_passage.py.
    assert sol.sol is not None
    tau_f = 0.42 * orb.period
    tau_r = orb.period - tau_f
    p0 = qbcp.state_pm_to_pv(sol.sol(0.0), 0.0, system)
    p1 = qbcp.state_pm_to_pv(sol.sol(tau_f), tau_f, system)
    tor0 = ConstantTorusPoint(pv=p0, system=system, omega_long=omega1)
    tor1 = ConstantTorusPoint(pv=p1, system=system, omega_long=omega1)
    info = {"sol": sol, "tau_f": tau_f, "tau_r": tau_r, "omega1": omega1}
    return system, tor0, tor1, info


def test_full_corrector_periodic_split_positive_control(
    periodic_split: tuple[qbcp.QBCPSystem, ConstantTorusPoint, ConstantTorusPoint, dict[str, Any]],
) -> None:
    """POSITIVE CONTROL: the full two-arc corrector, truth-seeded from a real
    #611 periodic orbit split, drives the combined residual to MACHINE precision
    and the INDEPENDENT Radau loop defect below 1e-3 (limited by the input orbit's
    own ~3e-6 closure, not the corrector). Proves the whole residual + analytic-
    Jacobian + independent-check path against a known closed trajectory."""
    _system, tor0, tor1, info = periodic_split
    sol, tau_f, tau_r, omega1 = info["sol"], info["tau_f"], info["tau_r"], info["omega1"]
    order = 28
    nodes, _ = cheb_diff_matrix(order)
    nf = np.array([sol.sol(t) for t in arc_node_times(0.0, tau_f, nodes)])
    nr = np.array([sol.sol(t) for t in arc_node_times(tau_f, tau_r, nodes)])
    th0 = np.array([0.0, 0.0])
    th1 = np.array([omega1 * tau_f, 0.0])
    z0 = pack_unknowns(nf, nr, tau_f, tau_r, th0, th1, th1.copy(), th0.copy())
    res = correct_qbcp_arc_connection(
        tor0, tor1, z0, order, tol=1e-6, closure_tol=1e-3, max_nfev=200, notes="posctl"
    )
    assert res.residual_rms < 1e-10
    assert res.match_norm_f < 1e-9 and res.match_norm_r < 1e-9
    assert res.phase_close_norm < 1e-9
    assert res.closure_loop_defect < 1e-3  # independent Radau, not the optimizer


def test_ghost_solution_rejected_by_independent_closure(
    periodic_split: tuple[qbcp.QBCPSystem, ConstantTorusPoint, ConstantTorusPoint, dict[str, Any]],
) -> None:
    """MANDATORY-CHECK DEMONSTRATION: from a poor (linear-interpolation) seed the
    SAME split converges to a GHOST -- a machine-zero nodal residual that is NOT a
    real trajectory. The algebraic residual is tiny yet the INDEPENDENT Radau loop
    defect is O(1) (~20,000 km arrival error), and ``converged`` is False. A small
    residual alone is never sufficient here."""
    _system, tor0, tor1, info = periodic_split
    sol, tau_f, tau_r, omega1 = info["sol"], info["tau_f"], info["tau_r"], info["omega1"]
    order = 28
    nodes, _ = cheb_diff_matrix(order)
    pm0, pm1 = sol.sol(0.0), sol.sol(tau_f)
    # linear interpolation in PM (node idx K = s=-1 = departure, idx 0 = arrival)
    nf = np.array([pm0 + (pm1 - pm0) * (1 - (s + 1) / 2) for s in nodes])
    nr = np.array([pm1 + (pm0 - pm1) * (1 - (s + 1) / 2) for s in nodes])
    th0 = np.array([0.0, 0.0])
    th1 = np.array([omega1 * tau_f, 0.0])
    z0 = pack_unknowns(nf, nr, tau_f, tau_r, th0, th1, th1.copy(), th0.copy())
    res = correct_qbcp_arc_connection(
        tor0, tor1, z0, order, tol=1e-6, closure_tol=1e-3, max_nfev=300, notes="ghost"
    )
    assert res.residual_rms < 1e-8  # machine-zero nodal residual (a ghost)
    assert res.closure_loop_defect > 1e-1  # but NOT a real trajectory
    assert res.converged is False  # the independent check vetoes it


# ---------------------------------------------------------------------------
# Real-torus endpoint control (SE-L2 GMOS torus).
# ---------------------------------------------------------------------------


def _build_se_l2_gmos_torus() -> QBCPTorus:
    """Full-mu Sun-Earth L2 QBCP GMOS torus (built exactly as #533/#538 do)."""
    system = qbcp.qbcp_default()
    mu_se = 1.0 / (system.mu_sun + 1.0)
    n_samples, n_modes = 5, 2
    sys_se = cr3bp.CR3BPSystem(
        mu=mu_se, primary="Sun", secondary="Earth", l_km=system.a_sun_nondim * 384400.0, t_s=1.0
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
        mu=system.mu,
        mu_sun=system.mu_sun,
        a_sun_nondim=system.a_sun_nondim,
        omega_sun_nondim=system.omega_sun_nondim,
        theta_sun0=system.theta_sun0,
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
    return correct_qbcp_torus(system, x0, n_modes, n_samples, ppi, amp, tol=1e-3)


@pytest.fixture(scope="module")
def se_gmos() -> QBCPTorus:
    return _build_se_l2_gmos_torus()


def test_torus_point_pv_gmos_matches_evaluate(se_gmos: QBCPTorus) -> None:
    """``torus_point_pv`` on a real GMOS QBCPTorus dispatches to the trusted
    ``evaluate_qbcp_torus`` (byte-identical PV)."""
    omega1 = se_gmos.omega_long
    for tl, tt in ((0.0, 0.0), (1.3, 2.1), (2.7, 0.9)):
        got = torus_point_pv(se_gmos, tl, tt, omega1)
        expect = evaluate_qbcp_torus(se_gmos, tl, tt)
        assert np.allclose(got, expect, atol=1e-14)


def test_on_torus_arc_reproduced_between_gmos_points(se_gmos: QBCPTorus) -> None:
    """A single collocation arc between two real SE-L2 GMOS torus points
    reproduces the genuine on-torus QBCP segment (residual to machine precision;
    the reconstructed nodes match a direct integration; an independent Radau
    re-propagation of the departure node lands on the collocation arrival to < 1
    km). The task's literal 'short arc between two nearby points on the same
    already-converged torus, confirmed by an independent integrator' positive
    control, exercising the real GMOS torus endpoint via ``torus_point_pv``."""
    from scipy.optimize import least_squares

    system = se_gmos.system
    omega1 = se_gmos.omega_long
    tl_a, tt_a = 0.6, 1.1
    p_a = torus_point_pv(se_gmos, tl_a, tt_a, omega1)
    t_a = tl_a / omega1
    tau = 0.5  # short on-torus segment
    pm_a = qbcp.state_pv_to_pm(p_a, t_a, system)
    seg = solve_ivp(
        lambda t, y: qbcp.qbcp_eom(t, y, system),
        (t_a, t_a + tau),
        pm_a,
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert seg.success
    # dense_output=True above (+ the success check just above) guarantees a
    # dense interpolant -- see the same pattern's comment in
    # src/cyclerfinder/search/perimoon_passage.py.
    assert seg.sol is not None
    pm_b = seg.sol(t_a + tau)
    p_b = qbcp.state_pm_to_pv(pm_b, t_a + tau, system)
    order = 20
    nodes, diff = cheb_diff_matrix(order)
    t_j = arc_node_times(t_a, tau, nodes)
    truth = np.array([seg.sol(t) for t in t_j])

    def single_arc_res(z: np.ndarray) -> np.ndarray:
        x = z.reshape(order + 1, 6)
        defects = arc_ode_defects(x, tau, t_a, diff, nodes, system)[1:order].reshape(-1)
        dep = qbcp.state_pm_to_pv(x[order], t_a, system) - p_a
        arr = qbcp.state_pm_to_pv(x[0], t_a + tau, system) - p_b
        return np.concatenate([defects, dep, arr])

    # Seed with a straight PM interpolation -- benign short arc, unique BVP.
    seed = np.array([pm_a + (pm_b - pm_a) * (1.0 - (s + 1.0) / 2.0) for s in nodes])
    sol = least_squares(
        single_arc_res, seed.reshape(-1), method="trf", xtol=1e-14, ftol=1e-14, gtol=1e-14
    )
    x = sol.x.reshape(order + 1, 6)
    assert float(np.sqrt(np.mean(single_arc_res(sol.x) ** 2))) < 1e-9
    assert np.max(np.abs(x - truth)) < 1e-6  # reconstructed nodes match integration
    # Independent Radau re-propagation of the departure node -> arrival node.
    prop = solve_ivp(
        lambda t, y: qbcp.qbcp_eom(t, y, system),
        (t_a, t_a + tau),
        x[order],
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    ).y[:, -1]
    assert float(np.linalg.norm(prop[:3] - x[0][:3])) * 384400.0 < 1.0  # < 1 km
