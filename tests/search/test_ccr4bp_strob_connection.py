"""#882 gates for the stroboscopic-map connection module.

Every expected value is an identity of the model (inverse maps, symplecticity,
exact invariance of a periodic orbit's circle in the unperturbed problem,
Floquet theory, first-order manifold scaling) or an independent computation
(core ``ccr4bp_eom`` / ``cr3bp_stm_eom`` integrations, finite differences).
None is a number produced earlier by the code under test.

Test system: Jupiter-Europa(-Ganymede), the UNSTABLE x-axis-symmetric 3:4
exterior resonant orbit with periapsis on the +x axis (Keplerian seed e = 0.1),
apojove about 1.33 Europa radii, well inside Ganymede's orbit (1.595). Its
Floquet multiplier is real, positive and about 4.0.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_strob_connection as sc

N_NODES = 201


@pytest.fixture(scope="module")
def physical() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


@pytest.fixture(scope="module")
def unperturbed(physical: ccr4bp.CCR4BPSystem) -> ccr4bp.CCR4BPSystem:
    return dataclasses.replace(physical, mu_gan=0.0)


@pytest.fixture(scope="module")
def orbit(physical: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    mu = physical.mu
    a = (4.0 / 3.0) ** (2.0 / 3.0)
    r = a * 0.9
    vin = math.sqrt((1.0 - mu) * (2.0 / r - 1.0 / a))
    x0 = r - mu
    s4, period, res = sc.symmetric_periodic_orbit(mu, x0, vin - x0, 4.0 * math.pi)
    assert res < 1e-11
    return s4, period


@pytest.fixture(scope="module")
def circle0(
    unperturbed: ccr4bp.CCR4BPSystem, orbit: tuple[np.ndarray, float]
) -> sc.InvariantCircle:
    s4, period = orbit
    return sc.seed_circle_from_periodic_orbit(unperturbed, s4, period, N_NODES)


@pytest.fixture(scope="module")
def bundles0(circle0: sc.InvariantCircle) -> sc.HyperbolicBundles:
    return sc.hyperbolic_bundles(circle0)


@pytest.fixture(scope="module")
def circle_phys(physical: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle) -> sc.InvariantCircle:
    return sc.correct_invariant_circle(physical, circle0.nodes, circle0.rho, t0=0.0, tol=1e-10)


@pytest.fixture(scope="module")
def bundles_phys(circle_phys: sc.InvariantCircle) -> sc.HyperbolicBundles:
    return sc.hyperbolic_bundles(circle_phys)


# ---------------------------------------------------------------------------
# Own-EOM gate: the vectorised planar RHS against the core module.
# ---------------------------------------------------------------------------


def test_batch_rhs_matches_core_eom(physical: ccr4bp.CCR4BPSystem) -> None:
    rng = np.random.default_rng(882)
    m = 7
    states = np.column_stack(
        [
            rng.uniform(0.6, 1.4, m),
            rng.uniform(-0.5, 0.5, m),
            rng.normal(0, 0.3, m),
            rng.normal(0, 0.3, m),
        ]
    )
    phi = rng.normal(size=(m, 4, 4))
    for t in (0.0, 3.7, -11.2):
        y = np.concatenate([states.T, phi.transpose(1, 2, 0).reshape(16, m)]).reshape(-1)
        out = sc.planar_rhs_batch(t, y, physical, m, True).reshape(20, m)
        for i in range(m):
            s6 = sc.to_state6(states[i])
            phi6 = np.zeros((6, 6))
            phi6[np.ix_(sc._PLANAR, sc._PLANAR)] = phi[i]
            ref = ccr4bp.ccr4bp_stm_eom(t, np.concatenate([s6, phi6.reshape(-1)]), physical)
            ref_state = sc.to_state4(ref[:6])
            ref_phi = ref[6:].reshape(6, 6)[np.ix_(sc._PLANAR, sc._PLANAR)]
            np.testing.assert_allclose(out[:4, i], ref_state, rtol=0, atol=1e-13)
            np.testing.assert_allclose(out[4:, i].reshape(4, 4), ref_phi, rtol=0, atol=1e-12)
            a = sc.planar_jacobian(t, states[i], physical)
            np.testing.assert_allclose(a @ phi[i], ref_phi, rtol=0, atol=1e-12)


def test_batch_map_matches_core_propagator(
    physical: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle
) -> None:
    period = physical.ganymede_synodic_period
    orb = sc.strob_iterates(physical, circle0.nodes, n=1, with_stm=True)
    assert orb.stms is not None
    for j in range(0, circle0.n_nodes, 25):
        arc = ccr4bp.propagate_ccr4bp(
            physical, sc.to_state6(circle0.nodes[j]), period, with_stm=True, rtol=1e-13, atol=1e-13
        )
        assert arc.stm is not None
        np.testing.assert_allclose(orb.states[1, j], sc.to_state4(arc.state_f), rtol=0, atol=1e-11)
        np.testing.assert_allclose(
            orb.stms[1, j], arc.stm[np.ix_(sc._PLANAR, sc._PLANAR)], rtol=0, atol=1e-9
        )


# ---------------------------------------------------------------------------
# G1: the stroboscopic map.
# ---------------------------------------------------------------------------


def test_g1_inverse_and_phase_periodicity(
    physical: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle
) -> None:
    x = circle0.nodes[17] + np.array([1e-3, -2e-3, 5e-4, 1e-3])
    for t0 in (0.0, 2.3):
        fx, _ = sc.strob_map(physical, x, n=1, t0=t0)
        back, _ = sc.strob_map(physical, fx, n=-1, t0=t0 + physical.ganymede_synodic_period)
        np.testing.assert_allclose(back, x, rtol=0, atol=1e-9)
    # the map is the same at t0 and t0 + P (time-periodic equations)
    a, _ = sc.strob_map(physical, x, n=1, t0=2.3)
    b, _ = sc.strob_map(physical, x, n=1, t0=2.3 + physical.ganymede_synodic_period)
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)
    # backward map from t0: F^-1 started at t0 inverts F started at t0 - P
    c, _ = sc.strob_map(physical, x, n=-1, t0=2.3)
    d, _ = sc.strob_map(physical, c, n=1, t0=2.3 - physical.ganymede_synodic_period)
    np.testing.assert_allclose(d, x, rtol=0, atol=1e-9)


def test_g1_stm_matches_central_differences(
    physical: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle
) -> None:
    x = circle0.nodes[40].copy()
    _, phi = sc.strob_map(physical, x, n=1, t0=1.0, with_stm=True)
    assert phi is not None
    h = 1e-6
    fd = np.empty((4, 4))
    for k in range(4):
        dx = np.zeros(4)
        dx[k] = h
        fp, _ = sc.strob_map(physical, x + dx, n=1, t0=1.0)
        fm, _ = sc.strob_map(physical, x - dx, n=1, t0=1.0)
        fd[:, k] = (fp - fm) / (2 * h)
    assert np.max(np.abs(fd - phi)) <= 1e-6 * np.max(np.abs(phi))


def test_g1_symplectic(physical: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle) -> None:
    """Canonical momenta px = vx - y, py = vy + x: the map is symplectic in the
    standard J after the constant linear change of variables T."""
    _, phi = sc.strob_map(physical, circle0.nodes[5], n=2, t0=0.7, with_stm=True)
    assert phi is not None
    t = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 1, 0], [1, 0, 0, 1]], dtype=float)
    m = t @ phi @ np.linalg.inv(t)
    j = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    assert abs(np.linalg.det(phi) - 1.0) <= 1e-8
    assert np.max(np.abs(m.T @ j @ m - j)) <= 1e-8 * max(1.0, np.max(np.abs(m)) ** 2)


# ---------------------------------------------------------------------------
# G2 / G3: unperturbed identities.
# ---------------------------------------------------------------------------


def test_g2_unperturbed_circle_is_invariant(
    unperturbed: ccr4bp.CCR4BPSystem, circle0: sc.InvariantCircle
) -> None:
    res, _ = sc.circle_residual(circle0)
    assert res <= 1e-9
    corrected = sc.correct_invariant_circle(unperturbed, circle0.nodes, circle0.rho, tol=1e-9)
    assert corrected.residual <= 1e-9
    assert np.max(np.abs(corrected.nodes - circle0.nodes)) <= 1e-9
    jac = [cr3bp.jacobi_constant(sc.to_state6(u), unperturbed.mu) for u in circle0.nodes]
    assert max(jac) - min(jac) <= 1e-10


def test_g3_unperturbed_bundles_match_floquet(
    unperturbed: ccr4bp.CCR4BPSystem,
    orbit: tuple[np.ndarray, float],
    bundles0: sc.HyperbolicBundles,
) -> None:
    s4, period = orbit
    y0 = np.concatenate([sc.to_state6(s4), np.eye(6).reshape(-1)])
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom,
        (0.0, period),
        y0,
        args=(unperturbed.mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    mono = sol.y[6:, -1].reshape(6, 6)[np.ix_(sc._PLANAR, sc._PLANAR)]
    ev = np.linalg.eigvals(mono)
    lam_floquet = float(np.max(np.abs(ev)))
    expected = lam_floquet ** (unperturbed.ganymede_synodic_period / period)
    assert abs(bundles0.lam_u * bundles0.lam_s - 1.0) <= 1e-6
    assert abs(bundles0.lam_u / expected - 1.0) <= 1e-5
    assert bundles0.residual_u <= 1e-7
    assert bundles0.residual_s <= 1e-7


# ---------------------------------------------------------------------------
# G4: first-order manifold consistency.
# ---------------------------------------------------------------------------


def test_g4_unstable_manifold_fundamental_domain_consistency(
    circle0: sc.InvariantCircle, bundles0: sc.HyperbolicBundles
) -> None:
    """Wu(theta, s, n+1) = Wu(theta + rho, lam_u*s, n) + O((eps*s)^2 * growth)."""
    theta, s, n = 0.9, 1.3, 3
    mism = []
    for eps in (1e-4, 5e-5, 2.5e-5):
        a, _ = sc.manifold_point(circle0, bundles0, "unstable", theta, s, n + 1, eps=eps)
        b, _ = sc.manifold_point(
            circle0, bundles0, "unstable", theta + circle0.rho, bundles0.lam_u * s, n, eps=eps
        )
        mism.append(float(np.linalg.norm(a - b)))
    r1 = mism[0] / mism[1]
    r2 = mism[1] / mism[2]
    assert 3.2 <= r1 <= 4.8 and 3.2 <= r2 <= 4.8, (mism, r1, r2)
    # bound: (eps*s)^2 times the growth over n periods, times an O(10) curvature constant
    assert mism[0] <= 10.0 * (1e-4 * s) ** 2 * bundles0.lam_u ** (n + 1)


# ---------------------------------------------------------------------------
# G5 (fast part): the hyperbolic circle at physical Ganymede mass.
# ---------------------------------------------------------------------------


def test_g5_perturbed_circle_converges_with_independent_closure(
    physical: ccr4bp.CCR4BPSystem, circle_phys: sc.InvariantCircle
) -> None:
    assert circle_phys.converged, circle_phys.residual_history
    assert circle_phys.residual <= 1e-9
    orb = sc.strob_iterates(physical, circle_phys.nodes, n=3, t0=circle_phys.t0)
    target = circle_phys.state(circle_phys.thetas + 3.0 * circle_phys.rho)
    closure = float(np.max(np.abs(orb.states[3] - target)))
    assert closure <= 1e-8
    assert circle_phys.fourier_tail() <= 1e-8


# ---------------------------------------------------------------------------
# G6: verification logic.
# ---------------------------------------------------------------------------


def test_g6_trivial_connection_rejected(
    circle0: sc.InvariantCircle, bundles0: sc.HyperbolicBundles
) -> None:
    eps = 1e-6
    x0, _ = sc.manifold_departure(circle0, bundles0, "unstable", 0.4, 1.0, eps=eps)
    target, _ = sc.strob_map(circle0.system, x0, n=1, t0=circle0.t0)
    ver = sc.verify_trajectory(
        circle0.system, x0, circle0.t0, 1, 1, circle0, junction_target=target, offset_size=eps
    )
    assert ver.junction_residual <= 1e-9
    assert not ver.genuine
    assert any("excursion" in r for r in ver.reasons)


def test_g6_phase_mismatched_connection_rejected(
    physical: ccr4bp.CCR4BPSystem,
    circle_phys: sc.InvariantCircle,
    bundles_phys: sc.HyperbolicBundles,
) -> None:
    """Regression test for defect 1: a junction that matches (x, y, vx, vy)
    exactly but at a different forcing phase is not one trajectory."""
    eps, n_u, n_s = 1e-5, 3, 10
    period = physical.ganymede_synodic_period
    z, _ = sc.manifold_point(circle_phys, bundles_phys, "stable", 1.0, 1.3, n_s, eps=eps)
    outcomes = {}
    for delta in (0.0, period / 3.0):
        # pull the honest stable point back n_u periods as if it were reached at phase t0 + delta
        d, _ = sc.strob_map(physical, z, n=-n_u, t0=circle_phys.t0 + delta + n_u * period)
        ver = sc.verify_trajectory(
            physical,
            d,
            circle_phys.t0 + delta,
            n_u,
            n_s,
            circle_phys,
            junction_target=z,
            offset_size=eps * 1.3,
        )
        outcomes[delta] = ver
    honest = outcomes[0.0]
    fake = outcomes[period / 3.0]
    # control: at the right phase the continued trajectory does fall onto the circle
    assert honest.phase_consistent
    assert honest.final_distance <= 1e-3 * honest.max_excursion
    # fabricated: 4D junction is exact, but phase flag and the physics both fail
    assert fake.junction_residual <= 1e-9
    assert not fake.phase_consistent
    assert fake.final_distance > 1e-3 * fake.max_excursion
    assert not fake.genuine
