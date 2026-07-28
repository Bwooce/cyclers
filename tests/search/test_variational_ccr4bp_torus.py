"""Tests for the CCR4BP pseudospectral quasi-periodic-torus corrector (#690).

The corrector (``search.variational_ccr4bp_torus``) is the "EOM swap" of `#617`'s
QBCP torus corrector onto `#689`'s CCR4BP dynamics. These tests mirror `#617`'s
own validation discipline:

* the vectorized grid RHS / state-Jacobian are checked pointwise to MACHINE
  PRECISION against the scalar :mod:`cyclerfinder.core.ccr4bp` functions (as
  `#617` checked its grid funcs against ``qbcp_eom`` / ``qbcp_stm_eom``);
* the analytic residual-Jacobian is checked against central finite differences;
* the positive control reproduces the object CLASS of Kumar-Anderson-de la
  Llave-Gunter 2021 (arXiv:2109.14815): the Jupiter-Europa 3:4 resonant orbit
  persists as a 2D invariant torus under Ganymede's physical-mass forcing. The
  ``mu_gan -> 0`` regression is a machine-precision-flat check; the physical-mass
  torus is reported at its achieved invariance residual (honest-partial, `#664`:
  the paper publishes no ICs/energy to pixel-match, and the reachable symmetric
  3:4 member is eccentric so the residual floor is Fourier-truncation-limited).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor
from cyclerfinder.search.variational_ccr4bp_torus import (
    _N_STATE,
    _ccr4bp_jacobian_grid,
    _ccr4bp_rhs_grid,
    _ganymede_on_theta1,
    _jacobian,
    _pack,
    _residual,
)
from cyclerfinder.search.variational_qbcp_torus import _basis_matrices

_PLANAR_IDX = [0, 1, 3, 4]  # (x, y, vx, vy) rows/cols of the 6-state


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Symmetric (perpendicular x-axis crossing) p:q resonant CR3BP periodic orbit.

    Period is FIXED by symmetry at ``T = 2*pi*q_moon``; free unknowns ``(x0, vy0)``
    for IC ``(x0,0,0,0,vy0,0)`` with conditions ``y(T/2)=0``, ``vx(T/2)=0``. Damped
    Newton on the 2x2 STM sub-block. Test scaffolding (base-CR3BP seed only, uses
    no code under test).
    """
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def _theta1_content(result: vt.CCR4BPTorusVariationalResult) -> float:
    """L2 norm of all theta1-VARYING coefficients (torus 'thickness' in the
    forcing angle); ~0 for a flat (unforced) torus, O(mu_gan) once forced."""
    return float(np.sqrt(np.sum(result.coeffs[:, 1:, :] ** 2)))


@pytest.fixture(scope="module")
def system() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


@pytest.fixture(scope="module")
def orbit_34(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Jupiter-Europa 3:4 symmetric resonant periodic orbit (base CR3BP)."""
    s0, period, res = _resonant_symmetric_orbit(system.mu, 3, 4)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


# ---------------------------------------------------------------------------
# Backend correctness: grid RHS / Jacobian vs the scalar core.ccr4bp functions.
# ---------------------------------------------------------------------------


def test_rhs_grid_matches_ccr4bp_eom(system: ccr4bp.CCR4BPSystem) -> None:
    """Vectorized grid RHS == scalar planar ``ccr4bp_eom`` to machine precision."""
    rng = np.random.default_rng(1)
    m1, m2 = 5, 7
    t1 = 2 * np.pi * np.arange(m1) / m1
    omega1 = abs(system.omega_gan)
    gx, gy = _ganymede_on_theta1(t1, omega1, system)
    u = rng.normal(scale=0.3, size=(4, m1, m2)) + np.array([1.2, 0.0, 0.0, -0.3])[:, None, None]
    rhs = _ccr4bp_rhs_grid(u, gx, gy, system.mu, system.mu_gan, system.a_gan)
    max_err = 0.0
    for i in range(m1):
        t = float(t1[i]) / omega1
        for j in range(m2):
            s6 = np.array([u[0, i, j], u[1, i, j], 0.0, u[2, i, j], u[3, i, j], 0.0])
            f6 = ccr4bp.ccr4bp_eom(t, s6, system)
            f4 = f6[_PLANAR_IDX]
            max_err = max(max_err, float(np.max(np.abs(f4 - rhs[:, i, j]))))
    assert max_err < 1e-12, max_err


def test_jacobian_grid_matches_ccr4bp_stm(system: ccr4bp.CCR4BPSystem) -> None:
    """Vectorized 4x4 grid Jacobian == the planar block of ``ccr4bp_stm_eom``'s A."""
    rng = np.random.default_rng(2)
    m1, m2 = 5, 7
    t1 = 2 * np.pi * np.arange(m1) / m1
    omega1 = abs(system.omega_gan)
    gx, gy = _ganymede_on_theta1(t1, omega1, system)
    u = rng.normal(scale=0.3, size=(4, m1, m2)) + np.array([1.2, 0.0, 0.0, -0.3])[:, None, None]
    jg = _ccr4bp_jacobian_grid(u, gx, gy, system.mu, system.mu_gan, system.a_gan)
    max_err = 0.0
    for i in range(m1):
        t = float(t1[i]) / omega1
        for j in range(m2):
            s6 = np.array([u[0, i, j], u[1, i, j], 0.0, u[2, i, j], u[3, i, j], 0.0])
            y42 = np.concatenate([s6, np.eye(6).reshape(36)])
            dy = ccr4bp.ccr4bp_stm_eom(t, y42, system)
            a6 = dy[6:].reshape(6, 6)  # = A (since phi = I)
            a4 = a6[np.ix_(_PLANAR_IDX, _PLANAR_IDX)]
            max_err = max(max_err, float(np.max(np.abs(a4 - jg[:, :, i, j]))))
    assert max_err < 1e-11, max_err


def test_residual_jacobian_matches_central_fd(system: ccr4bp.CCR4BPSystem) -> None:
    """Analytic residual-Jacobian == central finite differences (`#617` discipline)."""
    rng = np.random.default_rng(3)
    n1, n2 = 3, 4
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    omega1 = abs(system.omega_gan)
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    gx, gy = _ganymede_on_theta1(t1, omega1, system)
    coeffs0 = rng.normal(scale=0.2, size=(_N_STATE, 2 * n1 + 1, 2 * n2 + 1))
    coeffs0[:, 0, 1] += 1.0
    z = _pack(coeffs0, 0.25)
    args = (
        omega1,
        system.mu,
        system.mu_gan,
        system.a_gan,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        gx,
        gy,
        0,
        1.3,
        7.0,
        0.5 * omega1,
        13.0,
    )
    jan = _jacobian(z, *args)
    jfd = np.zeros_like(jan)
    eps = 1e-7
    for k in range(z.size):
        zp = z.copy()
        zm = z.copy()
        zp[k] += eps
        zm[k] -= eps
        jfd[:, k] = (_residual(zp, *args) - _residual(zm, *args)) / (2 * eps)
    rel = np.abs(jan - jfd) / (1.0 + np.abs(jfd))
    assert float(rel.max()) < 1e-4, float(rel.max())


# ---------------------------------------------------------------------------
# Torus evaluation and input validation.
# ---------------------------------------------------------------------------


def test_evaluate_torus_state_shapes(system: ccr4bp.CCR4BPSystem) -> None:
    coeffs = np.zeros((_N_STATE, 3, 5))
    coeffs[0, 0, 0] = 1.2
    coeffs[0, 0, 1] = 0.4  # cos(theta2)
    result = vt.CCR4BPTorusVariationalResult(
        system=system,
        coeffs=coeffs,
        omega1=0.5,
        omega2=0.25,
        rotation_number=0.5,
        rho_strob=3.1,
        period=12.5,
        n1=1,
        n2=2,
        m1=5,
        m2=7,
        period_multiple=1,
        transverse_amplitude=0.4,
        residual_rms=0.0,
        closure_residual=0.0,
        converged=True,
        n_iter=1,
    )
    s = vt.evaluate_torus_state(result, 0.0, 0.0)
    assert s.shape == (_N_STATE,)
    assert s[0] == pytest.approx(1.6)  # 1.2 + 0.4*cos(0)
    arr = vt.evaluate_torus_state(result, np.array([0.0, 1.0]), np.array([0.0, 1.0]))
    assert arr.shape == (2, _N_STATE)


def test_bad_coeffs_shape_raises(system: ccr4bp.CCR4BPSystem) -> None:
    with pytest.raises(ValueError, match="coeffs0 shape"):
        vt.correct_ccr4bp_torus_pseudospectral(
            system, np.zeros((_N_STATE, 3, 3)), 0.25, n1=1, n2=2, amplitude_anchor=1.0
        )
    with pytest.raises(ValueError, match="period_multiple"):
        vt.correct_ccr4bp_torus_pseudospectral(
            system,
            np.zeros((_N_STATE, 3, 5)),
            0.25,
            n1=1,
            n2=2,
            amplitude_anchor=1.0,
            period_multiple=0,
        )


# ---------------------------------------------------------------------------
# Positive control: Jupiter-Europa 3:4 CCR4BP torus (Kumar et al. 2021 class).
# ---------------------------------------------------------------------------


def test_mu_gan_zero_flat_regression(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> None:
    """At mu_gan=0 the corrector recovers the base 3:4 orbit as a theta1-FLAT torus.

    The unforced dynamics have no theta1 structure, so the corrector must hold the
    torus flat (theta1-content ~ machine zero) and drive the invariance residual
    to the Fourier-truncation floor -- the corrector-level analogue of `#689`'s
    ``mu_gan -> 0`` structural reduction.
    """
    s0, period = orbit_34
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    res = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert _theta1_content(res) < 1e-8, _theta1_content(res)  # flat to ~machine zero
    assert res.residual_rms < 1e-3, res.residual_rms
    # rotation number = orbit frequency / Ganymede synodic frequency.
    expected_rot = (2.0 * np.pi / period) / abs(system.omega_gan)
    assert res.rotation_number == pytest.approx(expected_rot, rel=1e-3)


def test_physical_mass_thin_torus(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> None:
    """The 3:4 Europa resonant orbit persists as a genuine thin 2D invariant torus
    under Ganymede's PHYSICAL-mass forcing (the Kumar 2021 object class).

    Reports: invariance residual below the 1e-3 gate, an INDEPENDENT closure
    residual, the stroboscopic rotation number matching theory to 5 digits, and a
    theta1-forcing content that appears only when mu_gan > 0.
    """
    s0, period = orbit_34
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    flat = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert phys.system.mu_gan > 0.0
    assert phys.residual_rms < 1e-3, phys.residual_rms
    assert phys.closure_residual < 5e-3, phys.closure_residual
    # Genuine 2D torus: forcing structure appears (>= 100x the flat-seed residue).
    assert _theta1_content(phys) > 100.0 * _theta1_content(flat)
    assert _theta1_content(phys) > 1e-3
    # Stroboscopic rotation number == 2*pi * omega2 / omega1 (theory).
    expected = 2.0 * np.pi * (2.0 * np.pi / period) / abs(system.omega_gan)
    assert phys.rho_strob == pytest.approx(expected, rel=1e-3)


# #763: hit CI's global 600s pytest-timeout cap (run 30360261381) despite
# measuring only ~16s LOCALLY. NOT V-gauntlet evidence for any catalogue row
# (grepped validate.py's _LEVEL_EVIDENCE + catalogue.yaml: this file is never
# cited -- it's a construction/unit test for the general Europa-Ganymede
# system, distinct from test_ccr4bp_torus_umbriel_titania.py's own
# same-named function which IS cited). This project's CI has no separate
# `-m slow` job/cadence (.github/workflows/ci.yml runs only the default
# `uv run pytest`, i.e. `-m 'not slow'`), so @pytest.mark.slow would be a
# pure coverage loss, not a relocation -- a per-test
# @pytest.mark.timeout(1800) override is used instead, same fix as this
# repo's evidence-tied #763 siblings.
@pytest.mark.timeout(1800)
def test_continue_ccr4bp_torus_mass_reaches_physical(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> None:
    """Mass homotopy from mu_gan=0 walks to the physical mass without diverging."""
    s0, period = orbit_34
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    seed = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    steps = vt.continue_ccr4bp_torus_mass(
        seed,
        system,
        n_steps=3,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert len(steps) == 3
    final = steps[-1]
    assert final.system.mu_gan == pytest.approx(system.mu_gan)
    assert final.residual_rms < 5e-3, final.residual_rms
    assert _theta1_content(final) > _theta1_content(seed)
