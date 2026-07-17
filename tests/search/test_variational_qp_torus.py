"""Seedless 2D pseudospectral quasi-periodic-torus corrector tests (#612).

Covers (a) the unit-level machinery -- real tensor-product Fourier basis and
its analytic derivative, free-variable pack/unpack, the analytic Jacobian
(verified against a finite-difference Jacobian), and the direct torus
evaluator -- and (b) two staged integration controls, mirroring
``test_variational_periodic_orbit_qbcp.py``'s pattern:

  * **L2 positive control** (regression floor, MANDATORY): the seedless 2D
    pseudospectral corrector, bootstrapped from a small-amplitude GMOS torus on
    the ALREADY-converging EM L2 near-bifurcation quasi-halo at C=3.15
    (`#555`'s working L2 case, latitudinal frequency in Owen & Baresi's ~0.0216
    regime), reproduces that torus -- same rotation number to ~1e-6, invariance
    residual at the truncation floor, independent closure tiny.

  * **L1 wall crossing** (the headline): on the EM L1 quasi-halo at C=3.15 --
    where ``genome.qp_tori.correct_qp_torus`` (GMOS, longitudinal direction
    resolved by stroboscopic FLOW INTEGRATION over a parent orbit with monodromy
    spectral radius ~1540) cannot converge above amp~0.01 and times out at
    amp~0.015 -- the integration-free pseudospectral corrector continues the
    transverse amplitude cleanly PAST the wall (to ``initial_torus_amplitude``-
    equivalent > 0.015), holding an independent closure residual ~1e-8 at every
    step. (`#612`'s live diagnosis, re-confirmed here, measured GMOS timing out
    at amp=0.015 in 60s and returning invariance residual 2.2e-3 at amp=0.02;
    that slow-failure comparison is documented in the module/OUTSTANDING, not
    re-run here to keep the suite fast.)

All pinned numbers were reproduced live in this session (2026-07-16) by running
the corrector directly, NOT copied from a docstring. The rotation numbers are
energy-pinned and stable; the residual/closure floors are truncation- and
integrator-dependent, so they are pinned with generous bounds (DOP853 stepping
is not bit-reproducible across libm versions -- see ``tests/genome/
test_qp_tori.py``'s own note).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import correct_qp_torus
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.nrho_continuation import SymmetricNRHO, correct_symmetric_nrho
from cyclerfinder.search.variational_qp_torus import (
    _basis_matrices,
    _cr3bp_jacobian_grid,
    _jacobian,
    _k2_first_harmonic_cols,
    _n_free,
    _pack,
    _residual,
    _transverse_amplitude,
    _unpack,
    continue_qp_torus_amplitude,
    discover_qp_torus,
    evaluate_torus_state,
)

MU = 0.012153643
SYS = cr3bp.CR3BPSystem(mu=MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0)


# ---------------------------------------------------------------------------
# Unit tests: machinery.
# ---------------------------------------------------------------------------


def test_basis_matrices_derivative_matches_analytic() -> None:
    """``_basis_matrices`` returns [1, cos k., sin k.] and its exact theta
    derivative [0, -k sin k., k cos k.] -- checked at several angles."""
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


def test_pack_unpack_round_trip() -> None:
    n1, n2 = 4, 2
    rng = np.random.default_rng(0)
    coeffs = rng.normal(size=(6, 2 * n1 + 1, 2 * n2 + 1))
    z = _pack(coeffs, 1.7, 0.3)
    assert z.shape == (_n_free(n1, n2),)
    c2, w1, w2 = _unpack(z, n1, n2)
    assert np.array_equal(c2, coeffs)
    assert (w1, w2) == (1.7, 0.3)


def test_cr3bp_jacobian_grid_matches_pointwise_stm_a_matrix() -> None:
    """The vectorized grid Jacobian equals ``cr3bp_stm_eom``'s A matrix at each
    grid point (extract A from the augmented RHS: dPhi = A when Phi = I)."""
    rng = np.random.default_rng(1)
    u = rng.normal(scale=0.3, size=(6, 3, 2))
    u[0] += 0.85  # keep away from primaries
    jf = _cr3bp_jacobian_grid(u, MU)
    for i in range(3):
        for j in range(2):
            state = u[:, i, j]
            y42 = np.concatenate([state, np.eye(6).reshape(36)])
            a_expected = cr3bp.cr3bp_stm_eom(0.0, y42, MU)[6:].reshape(6, 6)
            assert np.allclose(jf[:, :, i, j], a_expected, atol=1e-10)


def test_analytic_jacobian_matches_finite_difference() -> None:
    """The analytic residual Jacobian matches a central finite-difference
    Jacobian to high relative accuracy -- the load-bearing correctness gate for
    the whole corrector (a wrong Jacobian would silently mis-converge)."""
    rng = np.random.default_rng(2)
    n1, n2 = 4, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    coeffs = rng.normal(scale=0.05, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    coeffs[1, 0, 0] = 0.05
    z = _pack(coeffs, 1.7, 0.3)
    args = (MU, n1, n2, p1, p1d, p2, p2d, 0, 2, 0.05, 1.0, 0.176, 1.0)
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


def test_residual_has_four_gauge_rows() -> None:
    """The residual is PDE rows (6*m1*m2) plus exactly 4 gauge rows."""
    n1, n2 = 3, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    coeffs = np.zeros((6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    z = _pack(coeffs, 1.7, 0.3)
    r = _residual(z, MU, n1, n2, p1, p1d, p2, p2d, 0, 2, 0.0, 1.0, 0.176, 1.0)
    assert r.shape == (6 * m1 * m2 + 4,)


def test_rotation_number_gauge_row_value() -> None:
    """The 4th gauge row is ``rho_weight*(omega2 - rho_target*omega1)`` -- zero
    exactly when the ratio matches the target."""
    n1, n2 = 3, 2
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    coeffs = np.zeros((6, 2 * n1 + 1, 2 * n2 + 1))
    coeffs[0, 0, 0] = 0.85
    omega1, rho = 2.0, 0.074
    z = _pack(coeffs, omega1, rho * omega1)  # ratio exactly rho
    r = _residual(z, MU, n1, n2, p1, p1d, p2, p2d, 0, 2, 0.0, 1.0, rho, 3.0)
    assert r[-1] == pytest.approx(0.0, abs=1e-12)
    z2 = _pack(coeffs, omega1, (rho + 0.01) * omega1)
    r2 = _residual(z2, MU, n1, n2, p1, p1d, p2, p2d, 0, 2, 0.0, 1.0, rho, 3.0)
    assert r2[-1] == pytest.approx(3.0 * 0.01 * omega1, rel=1e-9)


def test_evaluate_torus_state_matches_series() -> None:
    """``evaluate_torus_state`` equals the tensor-product Fourier series
    evaluated directly at a scalar angle pair."""
    from cyclerfinder.search.variational_qp_torus import QPTorusVariationalResult

    n1, n2 = 3, 2
    rng = np.random.default_rng(5)
    coeffs = rng.normal(scale=0.1, size=(6, 2 * n1 + 1, 2 * n2 + 1))
    res = QPTorusVariationalResult(
        system=SYS,
        coeffs=coeffs,
        omega1=2.0,
        omega2=0.15,
        rotation_number=0.075,
        n1=n1,
        n2=n2,
        m1=9,
        m2=7,
        transverse_amplitude=0.0,
        residual_rms=0.0,
        closure_residual=0.0,
        converged=False,
        n_iter=0,
        jacobi=0.0,
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


# ---------------------------------------------------------------------------
# Halo builders (shared by the integration controls).
# ---------------------------------------------------------------------------


def _best_k(phi: float) -> int:
    bk, bd = 5, math.inf
    for kk in range(3, 81):
        for j in range(1, kk):
            if math.gcd(j, kk) == 1 and abs(phi - 2 * math.pi * j / kk) < bd:
                bd, bk = abs(phi - 2 * math.pi * j / kk), kk
    return bk


def _center_pair(
    r: SymmetricNRHO,
) -> tuple[NDArray[np.float64], list[np.complex128]]:
    s0 = np.array([r.x0, 0.0, r.z0, 0.0, r.ydot0, 0.0])
    eigs = floquet_multipliers(monodromy(SYS, s0, r.T_TU))
    cands = [e for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    return s0, cands


def _l1_halo_at_315() -> SymmetricNRHO | None:
    x0, z0, ydot0, t = 0.82431, -0.06058, 0.17154, 2.7647
    prev = None
    for _ in range(500):
        r = correct_symmetric_nrho(SYS, float(x0), z0, ydot0, t, with_monodromy=False)
        if not r.converged or abs(r.z0) < 5e-3:
            return None
        z0, ydot0, t = r.z0, r.ydot0, r.T_TU
        if prev is not None and (prev - 3.15) * (r.jacobi - 3.15) <= 0:
            return r
        prev = r.jacobi
        x0 -= 3e-5
    return None


def _l2_halo_at_315() -> SymmetricNRHO | None:
    lyap = None
    for xg, tg in ((1.1225, 3.43), (1.1204, 3.42), (1.1180, 3.40), (1.1250, 3.45)):
        cand = correct_symmetric_fixed_jacobi(
            SYS,
            x0_guess=xg,
            jacobi=3.15,
            period_guess=tg,
            ydot0_sign=+1.0,
            half_crossings=1,
            tol=1e-10,
            x0_bounds=(1.10, 1.20),
        )
        if cand.converged and abs(cand.period - 3.42) < 0.15:
            lyap = cand
            break
    if lyap is None:
        return None
    for z0s in (0.006, 0.008, 0.012, 0.016, 0.02):
        r = correct_symmetric_nrho(
            SYS,
            float(lyap.x0),
            z0s,
            float(lyap.ydot0),
            float(lyap.period),
            with_monodromy=False,
        )
        if r.converged and abs(r.z0) > 3e-3 and abs(r.jacobi - 3.15) < 3e-3:
            return r
    return None


# ---------------------------------------------------------------------------
# Integration control 1: L2 positive control (regression floor).
# ---------------------------------------------------------------------------


def test_l2_positive_control_reproduces_gmos_torus() -> None:
    """The seedless 2D pseudospectral corrector, bootstrapped from the
    already-converging small-amplitude GMOS L2 quasi-halo torus, reproduces it:
    same rotation number to ~1e-5, invariance residual at the truncation floor,
    independent closure tiny. This is the regression floor -- machinery must be
    sound on the EASY case before the L1 wall claim means anything."""
    l2 = _l2_halo_at_315()
    assert l2 is not None, "L2 halo build failed"
    s0, cands = _center_pair(l2)
    assert len(cands) >= 2
    k = _best_k(abs(math.atan2(cands[0].imag, cands[0].real)))
    gmos = correct_qp_torus(
        SYS,
        s0,
        l2.T_TU,
        (cands[0], cands[1]),
        k=k,
        n_trans=4,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
    )
    gmos_rot = gmos.omega_trans / gmos.omega_long
    # L2 near-bifurcation rotation number sits in Owen & Baresi's ~0.0216
    # latitudinal-frequency regime (#555 reports the family bottoms at ~0.0214).
    assert gmos_rot == pytest.approx(0.023272, rel=5e-3)

    res = discover_qp_torus(SYS, gmos, n1=10, n2=4, max_nfev=300)
    # Independent method reproduces the SAME torus.
    assert res.rotation_number == pytest.approx(gmos_rot, rel=1e-3)
    assert abs(res.omega2 - gmos.omega_trans) < 1e-4
    assert abs(res.omega1 - gmos.omega_long) < 1e-2
    # Residual/closure at the n1=10,n2=4 truncation floor (better than the GMOS
    # seed's own invariance residual of ~5.8e-6).
    assert res.residual_rms < 1e-5
    assert res.closure_residual < 1e-5


# ---------------------------------------------------------------------------
# Integration control 2: L1 wall crossing (headline).
# ---------------------------------------------------------------------------


def test_l1_crosses_gmos_amplitude_wall() -> None:
    """On the EM L1 quasi-halo at C=3.15 -- where GMOS cannot converge above
    amp~0.01 and times out at amp~0.015 (parent monodromy spectral radius
    ~1540; longitudinal direction resolved by stroboscopic flow integration) --
    the integration-free pseudospectral corrector continues the transverse
    amplitude cleanly PAST the wall, holding an independent closure residual
    ~1e-8 at every step. The rotation number stays energy-pinned at ~-0.074
    (confirming #555: at C=3.15 the L1 latitudinal frequency is NOT Owen &
    Baresi's 0.2739; that is an energy fact, not a corrector limitation)."""
    l1 = _l1_halo_at_315()
    assert l1 is not None, "L1 halo build failed"
    assert l1.jacobi == pytest.approx(3.15, abs=1e-3)
    s0, cands = _center_pair(l1)
    assert len(cands) >= 2
    # Parent halo is violently unstable -- the whole reason GMOS's stroboscopic
    # residual degrades. Confirm the monodromy amplification directly.
    sr = max(abs(e) for e in floquet_multipliers(monodromy(SYS, s0, l1.T_TU)))
    assert sr > 500.0, f"expected violently-unstable parent, got spectral radius {sr:.1f}"

    k = _best_k(abs(math.atan2(cands[0].imag, cands[0].real)))
    gmos = correct_qp_torus(
        SYS,
        s0,
        l1.T_TU,
        (cands[0], cands[1]),
        k=k,
        n_trans=4,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
    )
    seed = discover_qp_torus(SYS, gmos, n1=12, n2=6, max_nfev=300)
    assert seed.converged
    assert seed.rotation_number == pytest.approx(-0.074024, abs=1e-4)
    # Bootstrap amplitude scale: transverse_amplitude per GMOS initial amplitude.
    amp_scale = seed.transverse_amplitude / 5e-4
    assert amp_scale == pytest.approx(3.04, rel=0.1)

    # Continue PAST the GMOS wall. Target transverse amp 0.05 == GMOS
    # initial_torus_amplitude-equivalent ~0.0164, well beyond GMOS's amp=0.015
    # timeout / amp=0.02 invres=2.2e-3 failure (measured in #612's diagnosis).
    steps = continue_qp_torus_amplitude(seed, 0.05, n_steps=8, max_nfev=300)
    assert all(s.converged for s in steps), (
        "some continuation step failed to converge below the wall-crossing target"
    )
    last = steps[-1]
    gmos_equiv = last.transverse_amplitude / amp_scale
    assert gmos_equiv > 0.015, (
        f"did not cross the GMOS wall: GMOS-equivalent amplitude only {gmos_equiv:.4f}"
    )
    # Independent closure (short-time nonlinear flow, NOT the algebraic residual)
    # stays tiny at the large amplitude where GMOS diverges.
    assert last.closure_residual < 1e-6
    assert last.residual_rms < 1e-5
    # Rotation number held energy-pinned across the whole continuation.
    assert last.rotation_number == pytest.approx(-0.074024, abs=1e-4)
    # The large-amplitude torus is a genuinely big quasi-halo: it librates in z
    # across a range wider than the parent halo's own |z0|.
    t1 = np.linspace(0, 2 * np.pi, 25, endpoint=False)
    t2 = np.linspace(0, 2 * np.pi, 17, endpoint=False)
    zvals = np.array([evaluate_torus_state(last, float(a), float(b))[2] for a in t1 for b in t2])
    assert (zvals.max() - zvals.min()) > abs(l1.z0)
