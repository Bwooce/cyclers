"""Tests for the N=5 CRNBP pseudospectral quasi-periodic-torus corrector (`#720`).

`#720` continues `#690`'s own validated Jupiter-Europa 3:4 CCR4BP torus in
``mu_Io: 0 -> physical`` (~4.7e-5, `#717`'s Laplace-projected registry value),
using ``core.crnbp``'s N=5 EOM (`#717`). Mirrors `#690`'s own validation
discipline:

* the vectorized grid RHS / state-Jacobian (which deliberately OMITS the
  Io-Ganymede mutual coupling term, per `#717`'s proof that it contributes
  exactly zero) are checked pointwise to MACHINE PRECISION against the FULL,
  coupling-INCLUDED scalar :func:`cyclerfinder.core.crnbp.crnbp_eom` /
  :func:`cyclerfinder.core.crnbp.crnbp_stm_eom` -- the sharpest possible
  confirmation that the omission is exact, not approximate;
* the analytic residual-Jacobian is checked against central finite differences;
* the ``mu_io -> 0`` Tier-0 gate: the N=5 corrector at ``mu_io = 0``, warm-
  started from `#690`'s own converged Ganymede-only torus, lands at
  essentially the SAME invariance residual (the residual FUNCTION is bit-
  identical at the shared warm-start point -- checked directly -- so any
  divergence downstream is solver-path noise near the Fourier-truncation
  floor, the exact same noise `#690`'s own corrector exhibits under a
  warm-restart, checked here too);
* the physical finding: at ``n1 = 1`` (`#690`'s own default, tuned for
  Ganymede's forcing, which is PURELY first-harmonic in ``theta1``) the
  continuation lands on a bad branch once ``mu_Io > 0`` -- because the exact
  Laplace lock (``omega_Io = -2 * omega_Gan``) forces Io's own forcing to be
  PURELY SECOND-HARMONIC in ``theta1`` (``theta_Io(theta1) = theta_Io0 +
  2*theta1``), which an ``n1 = 1`` Fourier ansatz cannot represent at all.
  ``n1 = 2`` (the minimal basis that can represent Io's forcing) reproduces
  the ``mu_io = 0`` residual almost exactly, and the mu_Io continuation is
  then smooth and monotonic all the way to the physical value -- the FIRST
  known N=5 CRNBP invariant torus (subject to the mandatory literature-
  novelty gate before that word is used in a catalogue writeback).
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.variational_ccr4bp_torus as vt
import cyclerfinder.search.variational_crnbp_torus as vc
from cyclerfinder.genome.composed_moon_map import resonance_semimajor
from cyclerfinder.search.variational_crnbp_torus import (
    _N_STATE,
    _crnbp_jacobian_grid,
    _crnbp_rhs_grid,
    _jacobian,
    _pack,
    _perturber_grids,
    _residual,
)
from cyclerfinder.search.variational_qbcp_torus import _basis_matrices, _k2_first_harmonic_cols

_PLANAR_IDX = [0, 1, 3, 4]  # (x, y, vx, vy) rows/cols of the 6-state


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Symmetric p:q resonant CR3BP periodic orbit -- copied verbatim (test
    scaffolding, no code under test) from
    ``tests/search/test_variational_ccr4bp_torus.py``'s own helper of the
    same name, so this module has no cross-test-file import dependency.
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


def _k1_harmonic_content(coeffs: np.ndarray, n1: int, k: int) -> float:
    """L2 norm of the ``k``-th theta1 harmonic's coefficients (``k=0`` -> the
    theta1-constant row)."""
    if k == 0:
        return float(np.sqrt(np.sum(coeffs[:, 0, :] ** 2)))
    return float(np.sqrt(np.sum(coeffs[:, [k, n1 + k], :] ** 2)))


@pytest.fixture(scope="module")
def ccr4bp_system() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


@pytest.fixture(scope="module")
def crnbp_target() -> crnbp.CRNBPSystem:
    return crnbp.jupiter_europa_io_ganymede_default()


@pytest.fixture(scope="module")
def orbit_34(ccr4bp_system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Jupiter-Europa 3:4 symmetric resonant periodic orbit (base CR3BP)."""
    s0, period, res = _resonant_symmetric_orbit(ccr4bp_system.mu, 3, 4)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


@pytest.fixture(scope="module")
def system5(
    ccr4bp_system: ccr4bp.CCR4BPSystem, crnbp_target: crnbp.CRNBPSystem
) -> crnbp.CRNBPSystem:
    """A physical-mass N=5 system for grid-level RHS/Jacobian checks (not a
    converged torus -- just the parameter set)."""
    io = crnbp_target.perturbers[0]
    gan = crnbp.CRNBPPerturber(
        mu=ccr4bp_system.mu_gan,
        a=ccr4bp_system.a_gan,
        omega=ccr4bp_system.omega_gan,
        theta0=ccr4bp_system.theta_gan0,
    )
    return crnbp.CRNBPSystem(mu=ccr4bp_system.mu, perturbers=(io, gan))


# ---------------------------------------------------------------------------
# Backend correctness: grid RHS / Jacobian vs the FULL (coupling-included)
# scalar core.crnbp functions -- proves the coupling-term omission is exact.
# ---------------------------------------------------------------------------


def test_rhs_grid_matches_crnbp_eom(system5: crnbp.CRNBPSystem) -> None:
    """Vectorized grid RHS (coupling OMITTED) == scalar planar ``crnbp_eom``
    (coupling INCLUDED), to machine precision -- confirms `#717`'s
    zero-net-coupling-contribution proof transfers exactly into this grid
    context, at PHYSICAL Io+Ganymede masses (both nonzero simultaneously)."""
    rng = np.random.default_rng(1)
    m1, m2 = 5, 7
    t1 = 2 * np.pi * np.arange(m1) / m1
    clock = system5.perturbers[-1]
    omega1 = abs(clock.omega)
    perturbers_grid = _perturber_grids(system5, t1, omega1)
    u = rng.normal(scale=0.3, size=(4, m1, m2)) + np.array([1.2, 0.0, 0.0, -0.3])[:, None, None]
    rhs = _crnbp_rhs_grid(u, system5.mu, perturbers_grid)
    max_err = 0.0
    for i in range(m1):
        t = float(t1[i]) / omega1
        for j in range(m2):
            s6 = np.array([u[0, i, j], u[1, i, j], 0.0, u[2, i, j], u[3, i, j], 0.0])
            f6 = crnbp.crnbp_eom(t, s6, system5)
            f4 = f6[_PLANAR_IDX]
            max_err = max(max_err, float(np.max(np.abs(f4 - rhs[:, i, j]))))
    assert max_err < 1e-12, max_err


def test_jacobian_grid_matches_crnbp_stm(system5: crnbp.CRNBPSystem) -> None:
    """Vectorized 4x4 grid Jacobian == the planar block of ``crnbp_stm_eom``'s
    A matrix (coupling included there, provably zero-derivative, omitted here)."""
    rng = np.random.default_rng(2)
    m1, m2 = 5, 7
    t1 = 2 * np.pi * np.arange(m1) / m1
    clock = system5.perturbers[-1]
    omega1 = abs(clock.omega)
    perturbers_grid = _perturber_grids(system5, t1, omega1)
    u = rng.normal(scale=0.3, size=(4, m1, m2)) + np.array([1.2, 0.0, 0.0, -0.3])[:, None, None]
    jg = _crnbp_jacobian_grid(u, system5.mu, perturbers_grid)
    max_err = 0.0
    for i in range(m1):
        t = float(t1[i]) / omega1
        for j in range(m2):
            s6 = np.array([u[0, i, j], u[1, i, j], 0.0, u[2, i, j], u[3, i, j], 0.0])
            y42 = np.concatenate([s6, np.eye(6).reshape(36)])
            dy = crnbp.crnbp_stm_eom(t, y42, system5)
            a6 = dy[6:].reshape(6, 6)  # = A (since phi = I)
            a4 = a6[np.ix_(_PLANAR_IDX, _PLANAR_IDX)]
            max_err = max(max_err, float(np.max(np.abs(a4 - jg[:, :, i, j]))))
    assert max_err < 1e-11, max_err


def test_residual_jacobian_matches_central_fd(system5: crnbp.CRNBPSystem) -> None:
    """Analytic residual-Jacobian == central finite differences (`#617`/`#690` discipline)."""
    rng = np.random.default_rng(3)
    n1, n2 = 3, 4
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    clock = system5.perturbers[-1]
    omega1 = abs(clock.omega)
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    perturbers_grid = _perturber_grids(system5, t1, omega1)
    coeffs0 = rng.normal(scale=0.2, size=(_N_STATE, 2 * n1 + 1, 2 * n2 + 1))
    coeffs0[:, 0, 1] += 1.0
    z = _pack(coeffs0, 0.25)
    args = (
        omega1,
        system5.mu,
        perturbers_grid,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
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
# Tier-0 gate: mu_Io = 0 reduces to #690's own converged CCR4BP torus.
# ---------------------------------------------------------------------------


def _ccr4bp_seed_n1_2(
    ccr4bp_system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    s0, period = orbit_34
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        ccr4bp_system,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )


# #763: all four tests below (the file cited whole in the europa-3-4-crnbp-
# torus-jupiter-2026 row's own V1 evidence, validate.py's _LEVEL_EVIDENCE)
# hit the global 600s pytest-timeout cap in CI (run 30360261381) despite
# measuring only 12-124s LOCALLY -- a CI runner that is both slower per-core
# (documented BLAS/DOP853 platform gap) and contended (2-core runner, -n
# auto xdist) for these SVD-heavy least_squares(tr_solver="exact") solves.
# Per-test @pytest.mark.timeout(1800) overrides raise each test's own budget
# (established pattern: tests/search/test_png_lane_recovery.py's
# @pytest.mark.timeout(1200)) with real margin over the worst local
# measurement (124s, ~14.5x) rather than touching the pinned solver settings
# these tests exist to regression-guard.
@pytest.mark.timeout(1800)
def test_mu_io_zero_residual_function_bit_identical_to_ccr4bp(
    ccr4bp_system: ccr4bp.CCR4BPSystem,
    crnbp_target: crnbp.CRNBPSystem,
    orbit_34: tuple[np.ndarray, float],
) -> None:
    """At the SAME warm-start point, the N=5 residual (``mu_io=0``) and #690's
    own N=4 residual are bit-identical (max abs diff ~1e-16) -- the sharpest
    possible Tier-0 structural check, independent of any solver-path noise."""
    phys = _ccr4bp_seed_n1_2(ccr4bp_system, orbit_34)
    io = crnbp_target.perturbers[0]
    system5 = vc.crnbp_system_at_mu_io(
        ccr4bp_system.mu,
        ccr4bp_system.mu_gan,
        ccr4bp_system.a_gan,
        ccr4bp_system.omega_gan,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        theta_gan0=ccr4bp_system.theta_gan0,
    )
    n1, n2, m1, m2 = phys.n1, phys.n2, phys.m1, phys.m2
    omega1 = phys.omega1
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    z0 = _pack(phys.coeffs, phys.omega2)
    amp0 = vt._transverse_amplitude(phys.coeffs, n2)
    cos1, sin1 = _k2_first_harmonic_cols(n2)
    trans_by_coord = np.sqrt(
        np.sum(phys.coeffs[:, :, cos1] ** 2 + phys.coeffs[:, :, sin1] ** 2, axis=1)
    )
    phase2_coord = int(np.argmax(trans_by_coord))

    gx, gy = vt._ganymede_on_theta1(t1, omega1, ccr4bp_system)
    res_ccr4bp = vt._residual(
        z0,
        omega1,
        ccr4bp_system.mu,
        ccr4bp_system.mu_gan,
        ccr4bp_system.a_gan,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        gx,
        gy,
        phase2_coord,
        amp0,
        30.0,
        phys.rotation_number,
        100.0,
    )
    perturbers_grid = _perturber_grids(system5, t1, omega1)
    res_crnbp = _residual(
        z0,
        omega1,
        system5.mu,
        perturbers_grid,
        n1,
        n2,
        p1,
        p1d,
        p2,
        p2d,
        phase2_coord,
        amp0,
        30.0,
        phys.rotation_number,
        100.0,
    )
    assert float(np.max(np.abs(res_ccr4bp - res_crnbp))) < 1e-12


def test_mu_io_zero_seed_reproduces_ccr4bp_residual(
    ccr4bp_system: ccr4bp.CCR4BPSystem,
    crnbp_target: crnbp.CRNBPSystem,
    orbit_34: tuple[np.ndarray, float],
) -> None:
    """:func:`discover_crnbp_torus_from_ccr4bp_seed` at ``mu_io=0`` re-converges
    to essentially the SAME residual as #690's own torus (both below the
    project's established 1e-3 gate; the two need not match to machine
    precision after a re-solve, per the bit-identical-residual-function check
    above already covering that ground)."""
    phys = _ccr4bp_seed_n1_2(ccr4bp_system, orbit_34)
    io = crnbp_target.perturbers[0]
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert seed.system.perturbers[0].mu == 0.0
    assert seed.residual_rms < 1e-3, seed.residual_rms
    assert seed.residual_rms == pytest.approx(phys.residual_rms, rel=0.1)
    assert seed.closure_residual < 5e-3, seed.closure_residual


# ---------------------------------------------------------------------------
# The n1=1 harmonic-capacity finding (documented, not silently worked around).
# ---------------------------------------------------------------------------


@pytest.mark.timeout(1800)  # #763: see the file-level note above
def test_n1_1_cannot_represent_ios_second_harmonic_forcing(
    ccr4bp_system: ccr4bp.CCR4BPSystem,
    crnbp_target: crnbp.CRNBPSystem,
    orbit_34: tuple[np.ndarray, float],
) -> None:
    """Locks in the load-bearing methodological finding: the exact Laplace lock
    (``omega_Io = -2 * omega_Gan``) makes Io's synodic-frame position a PURE
    second-harmonic function of ``theta1`` (``theta_Io(theta1) = theta_Io0 +
    2*theta1``), which an ``n1=1`` Fourier ansatz (#690's own default, correct
    for Ganymede's purely-first-harmonic forcing) cannot represent at all.
    Turning on ANY nonzero mu_Io at ``n1=1`` throws the corrector onto a badly
    different branch (residual/closure orders of magnitude above the ``n1=1``
    ``mu_io=0`` floor) REGARDLESS of how small mu_Io is -- a representation-
    capacity wall, not a magnitude-proportional physical effect (confirmed by
    the residual staying roughly flat across the whole mu_Io sweep at n1=1,
    rather than shrinking as mu_Io -> 0). ``n1=2`` fixes this cleanly (see
    ``test_continuation_reaches_physical_mu_io_at_n1_2`` below).
    """
    s0, period = orbit_34
    phys1 = vt.discover_ccr4bp_torus_from_resonant_orbit(
        ccr4bp_system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    io = crnbp_target.perturbers[0]
    seed1 = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys1,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    target1 = crnbp.CRNBPSystem(mu=crnbp_target.mu, perturbers=(io, crnbp_target.perturbers[1]))
    steps1 = vc.continue_crnbp_torus_mu_io(
        seed1,
        target1,
        n_steps=3,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert len(steps1) == 3
    # Every step lands far above both the mu_io=0 seed's own residual AND the
    # project's 1e-3 physical-torus gate -- a genuine branch jump, not a
    # small forced correction.
    for step in steps1:
        assert step.residual_rms > 20.0 * seed1.residual_rms, (
            step.residual_rms,
            seed1.residual_rms,
        )
        assert step.residual_rms > 1e-3, (
            step.residual_rms
        )  # above the project's own physical-torus gate
    # Flat across the sweep (not shrinking toward the small-mu_io end) --
    # the signature of a representation-capacity wall, not a magnitude effect.
    assert steps1[0].residual_rms == pytest.approx(steps1[-1].residual_rms, rel=0.1)


# ---------------------------------------------------------------------------
# The n1=2 mu_Io continuation: the actual #720 discovery result.
# ---------------------------------------------------------------------------


@pytest.mark.timeout(1800)  # #763: see the file-level note above
def test_continuation_reaches_physical_mu_io_at_n1_2(
    ccr4bp_system: ccr4bp.CCR4BPSystem,
    crnbp_target: crnbp.CRNBPSystem,
    orbit_34: tuple[np.ndarray, float],
) -> None:
    """The #720 headline result: at ``n1=2`` (the minimal basis able to
    represent Io's second-harmonic forcing), the mu_Io continuation from
    #690's own validated Europa 3:4 CCR4BP torus (``mu_io=0``) reaches the
    FULL physical Io mass parameter (~4.7e-5) with a SMOOTH, monotonic,
    small (~1.5%) residual/closure increase and a stable rotation number --
    the first known N=5 CRNBP invariant torus (self-continuation-validated;
    literature-novelty status is NOT adjudicated by this test)."""
    phys = _ccr4bp_seed_n1_2(ccr4bp_system, orbit_34)
    io = crnbp_target.perturbers[0]
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        crnbp_target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert len(steps) == 8, "continuation diverged before reaching the physical mu_io"
    final = steps[-1]
    assert final.system.perturbers[0].mu == pytest.approx(crnbp_target.perturbers[0].mu)
    assert final.system.perturbers[1].mu == pytest.approx(crnbp_target.perturbers[1].mu)

    # Residuals stay well inside #690's own established gates at every step.
    for step in steps:
        assert step.residual_rms < 1e-3, step.residual_rms
        assert step.closure_residual < 5e-3, step.closure_residual

    # Smooth, small, monotonic (to within re-solve noise) drift from the
    # mu_io=0 seed to the physical endpoint (NOT a branch jump like the n1=1
    # case above). `#723`: at the corrected physical theta_io0 (`pi -
    # 2*theta_gan0`, was the buggy all-aligned 0.0 -- see crnbp.py's module
    # docstring), the very first continuation step dips by ~9.5e-10 relative
    # (~9.5e-9 absolute) below the seed before resuming its monotonic climb --
    # solver re-solve noise, not a branch jump (two orders of magnitude below
    # the ~1.33e-6 absolute drift over the WHOLE continuation, and four+
    # orders below the >20x-seed branch-jump signature the n1=1 test above
    # locks in). The per-step tolerance is loosened from the prior
    # (theta_io0=0-only-valid) 1e-9 to 2e-8 to absorb this while still
    # tightly bounding genuine divergence.
    assert final.residual_rms < 5.0 * seed.residual_rms
    assert final.closure_residual < 5.0 * seed.closure_residual
    residuals = [seed.residual_rms] + [s.residual_rms for s in steps]
    assert all(b >= a - 2e-8 for a, b in itertools.pairwise(residuals)), residuals

    # Rotation number essentially unperturbed by Io's (tiny) forcing.
    assert final.rotation_number == pytest.approx(seed.rotation_number, rel=1e-3)

    # Io's own second-harmonic content is present at the physical endpoint and
    # is LARGER than at the mu_io=0 seed (a genuine, if small, forcing signature).
    k2_seed = _k1_harmonic_content(seed.coeffs, seed.n1, 2)
    k2_final = _k1_harmonic_content(final.coeffs, final.n1, 2)
    assert k2_final > k2_seed


@pytest.mark.timeout(1800)  # #763: see the file-level note above
def test_one_shot_direct_to_physical_matches_continuation(
    ccr4bp_system: ccr4bp.CCR4BPSystem,
    crnbp_target: crnbp.CRNBPSystem,
    orbit_34: tuple[np.ndarray, float],
) -> None:
    """Robustness check: correcting DIRECTLY at the physical mu_Io (no
    intermediate continuation steps, cold-started from the mu_io=0 seed's own
    coefficients) converges to essentially the same torus as the stepped
    continuation -- confirming the result is not an artifact of a fragile,
    finely-stepped homotopy path."""
    phys = _ccr4bp_seed_n1_2(ccr4bp_system, orbit_34)
    io = crnbp_target.perturbers[0]
    oneshot = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=io.mu,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert oneshot.residual_rms < 1e-3, oneshot.residual_rms
    assert oneshot.closure_residual < 5e-3, oneshot.closure_residual
