"""Gates for the #884 Sun-forced continuation module.

Every EXPECTED value is an identity (CR3BP limit, map inverse, symplecticity,
Jacobi conservation), an independent computation (finite differences, a second
integrator, an existing independent EOM implementation), or a prediction of the
published theory (Brown, Peterson, Henry & Scheeres 2025, Proposition 1 and
Proposition 3: zeros of the Melnikov function at the symmetric phases). No value
produced earlier by this module is on the expected side.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search import sun_forced_periodic_884 as sf

MODEL = sf.default_model()
TG = MODEL.tg
# Braik-Ross 2026 Table 2 LL1 representative (data/golden/braik_ross_2026_em_family_ics.yaml).
LL1 = (0.8115256290557147, 0.2561843220006502, 2.946253150022597)
STATES = [
    np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0]),
    np.array([0.8, -0.2, 0.05, 0.1, 0.2, -0.03]),
    np.array([-0.3, 0.7, -0.1, -0.4, 0.05, 0.02]),
]


@pytest.fixture(scope="module")
def l1_member() -> sf.SymmetricOrbit:
    """L1 planar Lyapunov member at T* = Tg/2 (a = 2), selected by period."""
    member, info = sf.continue_to_period(MODEL, np.array([LL1[0], 0.0, LL1[1]]), LL1[2], TG / 2)
    assert member is not None, info
    assert member.residual < 1e-11
    return member


def test_tg_is_the_models_synodic_period() -> None:
    sysd = bcr4bp.andreu_default()
    assert pytest.approx(2 * math.pi / sysd.omega_sun_nondim, rel=1e-15) == TG
    assert pytest.approx(bcr4bp.sun_commensurate_period(sysd.omega_sun_nondim, 1)) == TG


@pytest.mark.parametrize("state", STATES)
def test_rhs_matches_independent_implementations(state: np.ndarray) -> None:
    """eps = 0 equals cr3bp_eom; eps = 1 equals bcr4bp_eom (separately written code)."""
    t, th0 = 1.234, 0.7
    f0 = sf._rhs(t, state, MODEL.mu, 0.0, MODEL.a_sun, MODEL.omega_sun, th0)
    np.testing.assert_allclose(f0, cr3bp.cr3bp_eom(t, state, MODEL.mu), rtol=0, atol=1e-14)
    sysd = bcr4bp.BCR4BPSystem(
        mu=MODEL.mu,
        mu_sun=MODEL.mu_sun_phys,
        a_sun_nondim=MODEL.a_sun,
        omega_sun_nondim=MODEL.omega_sun,
        theta_sun0=th0,
    )
    f1 = sf._rhs(t, state, MODEL.mu, MODEL.mu_sun_phys, MODEL.a_sun, MODEL.omega_sun, th0)
    np.testing.assert_allclose(f1, bcr4bp.bcr4bp_eom(t, state, sysd), rtol=1e-13, atol=1e-13)
    y = np.concatenate([state, np.eye(6).reshape(36)])
    a_ref = bcr4bp.bcr4bp_stm_eom(t, y, sysd)[6:]
    a_new = sf._rhs_var(
        t,
        np.concatenate([y, np.zeros(6)]),
        MODEL.mu,
        MODEL.mu_sun_phys,
        1.0,
        MODEL.a_sun,
        MODEL.omega_sun,
        th0,
    )[6:42]
    np.testing.assert_allclose(a_new, a_ref, rtol=1e-12, atol=1e-12)


def test_jacobi_conserved_with_sun_off(l1_member: sf.SymmetricOrbit) -> None:
    sol = sf.propagate(MODEL, 0.0, 0.3, l1_member.state, 0.0, 3 * TG, dense=True)
    c0 = sf.jacobi(l1_member.state, MODEL.mu)
    cs = [sf.jacobi(sol.sol(t), MODEL.mu) for t in np.linspace(0, 3 * TG, 50)]
    assert max(abs(c - c0) for c in cs) < 1e-10


def test_map_inverse() -> None:
    x0 = np.array([-0.5, 0.0, 0.0, 0.0, -1.2, 0.0])
    fwd = sf.propagate(MODEL, 1.0, 0.9, x0, 0.0, TG).state_f
    back = sf.propagate(MODEL, 1.0, 0.9, fwd, TG, 0.0).state_f
    np.testing.assert_allclose(back, x0, rtol=0, atol=1e-10)


def test_stm_and_eps_sensitivity_against_finite_differences() -> None:
    x0 = STATES[1]
    t1, th0, eps = 2.0, 1.1, 0.6
    arc = sf.propagate(MODEL, eps, th0, x0, 0.0, t1, variational=True)
    assert arc.stm is not None and arc.dxdeps is not None
    h = 1e-6
    fd = np.empty((6, 6))
    for j in range(6):
        dx = np.zeros(6)
        dx[j] = h
        fp = sf.propagate(MODEL, eps, th0, x0 + dx, 0.0, t1).state_f
        fm = sf.propagate(MODEL, eps, th0, x0 - dx, 0.0, t1).state_f
        fd[:, j] = (fp - fm) / (2 * h)
    np.testing.assert_allclose(arc.stm, fd, rtol=1e-6, atol=1e-6)
    he = 1e-3
    fp = sf.propagate(MODEL, eps + he, th0, x0, 0.0, t1).state_f
    fm = sf.propagate(MODEL, eps - he, th0, x0, 0.0, t1).state_f
    np.testing.assert_allclose(arc.dxdeps, (fp - fm) / (2 * he), rtol=1e-5, atol=1e-8)


def test_monodromy_symplectic_in_canonical_coordinates() -> None:
    arc = sf.propagate(MODEL, 1.0, 0.4, STATES[1], 0.0, TG, variational=True)
    assert arc.stm is not None
    scale = float(np.linalg.norm(arc.stm)) ** 2
    assert sf.symplectic_defect(arc.stm) / scale < 1e-10
    # Cartesian rotating-frame STM is NOT symplectic with the standard J: the
    # canonical transformation is what makes the identity hold.
    j = np.zeros((6, 6))
    j[:3, 3:] = np.eye(3)
    j[3:, :3] = -np.eye(3)
    assert float(np.max(np.abs(arc.stm.T @ j @ arc.stm - j))) / scale > 1e-6


def test_melnikov_quadrature_matches_variational_and_finite_eps(
    l1_member: sf.SymmetricOrbit,
) -> None:
    """Three independent evaluations of the bifurcation function agree."""
    th0 = 0.4
    quad = float(sf.melnikov_scan(MODEL, l1_member.state, TG, np.array([th0]))[0])
    var = sf.melnikov_variational(MODEL, l1_member.state, TG, th0)
    assert abs(quad - var) < 1e-8 * max(1.0, abs(var))
    # Finite-eps derivative of the Jacobi change. The identity
    # dC/deps = grad C . dX/deps = -2 int v . a_sun holds on any arc; it is checked
    # on a mildly unstable arc because on the L1 member (max |lambda| ~ 6e5) the
    # O(eps^2) terms are amplified far beyond the first-order term at any usable eps.
    x0 = np.array([-0.5, 0.0, 0.0, 0.0, -1.2, 0.0])
    var2 = sf.melnikov_variational(MODEL, x0, TG, th0)
    quad2 = float(sf.melnikov_scan(MODEL, x0, TG, np.array([th0]))[0])
    eps = 1e-4
    xp = sf.propagate(MODEL, eps, th0, x0, 0.0, TG).state_f
    xm = sf.propagate(MODEL, -eps, th0, x0, 0.0, TG).state_f
    fin = (sf.jacobi(xp, MODEL.mu) - sf.jacobi(xm, MODEL.mu)) / (2 * eps)
    assert abs(quad2 - var2) < 1e-8 * abs(var2)
    assert abs(fin - var2) < 1e-6 * abs(var2)


def test_melnikov_period_two_pi_over_a(l1_member: sf.SymmetricOrbit) -> None:
    """With T* = (n/a) Tg the Melnikov function is 2 pi / a periodic (here a = 2)."""
    th = np.array([0.2, 0.9, 1.3])
    m0 = sf.melnikov_scan(MODEL, l1_member.state, TG, th)
    m1 = sf.melnikov_scan(MODEL, l1_member.state, TG, th + math.pi)
    assert float(np.max(np.abs(m1 - m0))) < 1e-9 * float(np.max(np.abs(m0)))


def test_melnikov_zeros_at_symmetric_phases(l1_member: sf.SymmetricOrbit) -> None:
    """Brown et al. Prop. 3 (via the reversing symmetry): x-axis start, zeros at 0, pi/2."""
    _, _, zeros = sf.find_zeros(MODEL, l1_member.state, TG, 2, n_grid=61)
    assert len(zeros) == 2
    assert min(abs(zeros[0]), abs(zeros[0] - math.pi)) < 1e-7
    assert abs(zeros[1] - math.pi / 2) < 1e-7


def test_forced_orbit_closure_independent_integrator_and_reversibility(
    l1_member: sf.SymmetricOrbit,
) -> None:
    prob = sf.ShootingProblem(MODEL, TG, 0.0, 2)
    xs0 = prob.nodes_from_orbit(l1_member.state)
    br = sf.continue_in_eps(prob, xs0, eps_target=0.2, ds0=0.05)
    assert br.stop_reason == "reached_target"
    path = br.eps[:-1]  # the last entry is the exact-target correction after an overshoot
    assert not br.folds
    assert all(b > a for a, b in itertools.pairwise(path))  # monotone in eps
    xs1 = br.nodes[-1]
    res, _, _, _ = prob.evaluate(xs1, 0.2)
    assert float(np.max(np.abs(res))) < 1e-10
    for k in range(2):
        t = prob.times()
        arc = sf.propagate(MODEL, 0.2, 0.0, xs1[k], t[k], t[k + 1], method="Radau")
        assert float(np.max(np.abs(arc.state_f - xs1[(k + 1) % 2]))) < 1e-7
    back = sf.continue_in_eps(prob, xs0, reverse_from=(xs1, 0.2), ds0=0.05)
    assert back.stop_reason == "reached_target"
    assert back.eps[-1] == 0.0
    assert float(np.max(np.abs(back.nodes[-1] - xs0))) < 1e-6


def test_symmetric_forced_orbit_matches_half_period_shooting(
    l1_member: sf.SymmetricOrbit,
) -> None:
    """Independent code path: perpendicular-crossing shooting in the forced model.

    With theta0 = 0 the reversing symmetry (y, vx, t, theta) -> (-y, -vx, -t, -theta)
    makes a periodic orbit through a perpendicular x-axis crossing at t = 0 cross
    perpendicularly again at t = P/2 (Sun phase pi there).
    """
    eps = 0.2
    prob = sf.ShootingProblem(MODEL, TG, 0.0, 2)
    br = sf.continue_in_eps(prob, prob.nodes_from_orbit(l1_member.state), eps_target=eps)
    xa = br.nodes[-1][0]
    assert abs(xa[1]) < 1e-9 and abs(xa[3]) < 1e-9
    u = np.array([l1_member.x0, l1_member.vy0])
    for _ in range(20):
        x0 = np.array([u[0], 0.0, 0.0, 0.0, u[1], 0.0])
        arc = sf.propagate(MODEL, eps, 0.0, x0, 0.0, TG / 2, variational=True)
        assert arc.stm is not None
        r = arc.state_f[[1, 3]]
        if float(np.max(np.abs(r))) < 1e-12:
            break
        u = u - np.linalg.solve(arc.stm[np.ix_([1, 3], [0, 4])], r)
    np.testing.assert_allclose([xa[0], xa[4]], u, rtol=0, atol=1e-9)
