"""Tests for :mod:`cyclerfinder.core.kepler_stm` (analytic two-body STM).

CONSISTENCY-TEST PATTERN (same discipline as ``tests/nbody`` for
``flyby_gradients``): Shepperd 1985 was read 2026-06-13 (#233) and prints NO
worked numeric example — it is a pure theory/algorithm paper (symbolic
universal-variable formulation, M-matrix STM, U_n functions, and the
continued-fraction evaluation; Appendix A is the algorithm summary, not a
numeric case). There is therefore no wireable printed golden, so the analytic
STM is validated against central-difference finite differences of the
INDEPENDENT existing propagator :func:`cyclerfinder.core.kepler.propagate`
across all conic regimes, plus the STM group properties (identity at
``dt == 0``, symplecticity, unit determinant, composition). These are
internal-consistency checks, never sourced goldens.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.kepler import coe_to_rv, propagate
from cyclerfinder.core.kepler_stm import Mat6, Vec3, shepperd_stm

_V_CIRC_1AU = sqrt(MU_SUN_KM3_S2 / AU_KM)
_V_ESC_1AU = sqrt(2.0) * _V_CIRC_1AU
_PERIOD_1AU_S = 2.0 * pi * sqrt(AU_KM**3 / MU_SUN_KM3_S2)

# Regime sweep: (id, r0 [km], v0 [km/s], dt [s]). All states are 3-dimensional
# (out-of-plane components included) except where the regime is the point.
_REGIME_CASES: list[tuple[str, Vec3, Vec3, float]] = [
    (
        "elliptic-low-e",
        np.array([AU_KM, 0.0, 0.0]),
        np.array([0.0, 1.05 * _V_CIRC_1AU, 0.02 * _V_CIRC_1AU]),
        200.0 * SECONDS_PER_DAY,
    ),
    (
        # rp = 0.2 AU, a = 1 AU -> e = 0.8, started at periapsis.
        "elliptic-high-e",
        np.array([0.2 * AU_KM, 0.0, 0.0]),
        np.array([0.0, sqrt(MU_SUN_KM3_S2 * (2.0 / (0.2 * AU_KM) - 1.0 / AU_KM)), 0.0]),
        120.0 * SECONDS_PER_DAY,
    ),
    (
        "hyperbolic",
        np.array([AU_KM, 0.0, 0.0]),
        np.array([0.0, 1.5 * _V_ESC_1AU, 0.1 * _V_CIRC_1AU]),
        120.0 * SECONDS_PER_DAY,
    ),
    (
        # 1 + 1e-6 of escape speed: |alpha| ~ 2.7e-14 1/km, the parabolic
        # bootstrap branch of the Newton solve.
        "near-parabolic",
        np.array([AU_KM, 0.0, 0.0]),
        np.array([0.0, (1.0 + 1.0e-6) * _V_ESC_1AU, 0.0]),
        100.0 * SECONDS_PER_DAY,
    ),
    (
        "multi-rev-elliptic",
        np.array([AU_KM, 0.0, 0.0]),
        np.array([0.0, _V_CIRC_1AU, 0.0]),
        3.7 * _PERIOD_1AU_S,
    ),
    (
        "backward",
        np.array([AU_KM, 0.0, 0.0]),
        np.array([0.0, 0.98 * _V_CIRC_1AU, 0.05 * _V_CIRC_1AU]),
        -150.0 * SECONDS_PER_DAY,
    ),
]

_REGIME_IDS = [case[0] for case in _REGIME_CASES]

# Central-difference step scale ~ cbrt(machine eps): truncation ~ h^2 ~ 4e-11
# relative, roundoff ~ eps/h ~ 4e-11 relative -- comfortably inside _FD_RTOL.
_FD_STEP_REL = 6.0e-6
_FD_RTOL = 1.0e-6

_J_CANONICAL = np.block(
    [
        [np.zeros((3, 3)), np.eye(3)],
        [-np.eye(3), np.zeros((3, 3))],
    ]
)


def _fd_stm(r0: Vec3, v0: Vec3, dt: float) -> Mat6:
    """Central-difference STM from the independent ``kepler.propagate``."""
    phi = np.zeros((6, 6), dtype=np.float64)
    h_r = _FD_STEP_REL * float(np.linalg.norm(r0))
    h_v = _FD_STEP_REL * float(np.linalg.norm(v0))
    for j in range(3):
        e = np.zeros(3)
        e[j] = h_r
        r_p, v_p = propagate(r0 + e, v0, dt)
        r_m, v_m = propagate(r0 - e, v0, dt)
        phi[0:3, j] = (r_p - r_m) / (2.0 * h_r)
        phi[3:6, j] = (v_p - v_m) / (2.0 * h_r)
    for j in range(3):
        e = np.zeros(3)
        e[j] = h_v
        r_p, v_p = propagate(r0, v0 + e, dt)
        r_m, v_m = propagate(r0, v0 - e, dt)
        phi[0:3, 3 + j] = (r_p - r_m) / (2.0 * h_v)
        phi[3:6, 3 + j] = (v_p - v_m) / (2.0 * h_v)
    return phi


def _max_block_rel_err(phi_a: Mat6, phi_b: Mat6) -> float:
    """Max Frobenius-relative error over the four 3x3 quadrants.

    Block-wise (not element-wise) because the quadrants carry different units
    (s, 1/s, dimensionless) and individual elements pass through zero.
    """
    errs = []
    for i0 in (0, 3):
        for j0 in (0, 3):
            a = phi_a[i0 : i0 + 3, j0 : j0 + 3]
            b = phi_b[i0 : i0 + 3, j0 : j0 + 3]
            errs.append(float(np.linalg.norm(a - b)) / float(np.linalg.norm(b)))
    return max(errs)


def _nondimensionalise(phi: Mat6, length_km: float) -> Mat6:
    """Symplectic similarity scaling to balanced (canonical) units.

    ``S = diag(1/L, 1/L, 1/L, T/L, T/L, T/L)`` with ``T = sqrt(L^3/mu)``
    satisfies ``S^T J S = (T/L^2) J``, so ``S Phi S^-1`` is symplectic iff
    ``Phi`` is — but with O(1) entries, making the residual meaningful.
    """
    t_scale = sqrt(length_km**3 / MU_SUN_KM3_S2)
    s = np.diag([1.0 / length_km] * 3 + [t_scale / length_km] * 3)
    s_inv = np.diag([length_km] * 3 + [length_km / t_scale] * 3)
    result: Mat6 = s @ phi @ s_inv
    return result


def test_zero_dt_identity() -> None:
    """``dt == 0`` returns the original state (as copies) and Phi = I."""
    r0 = np.array([AU_KM, 0.0, 0.0])
    v0 = np.array([0.0, _V_CIRC_1AU, 0.0])
    r, v, phi = shepperd_stm(r0, v0, 0.0)
    assert np.array_equal(r, r0)
    assert np.array_equal(v, v0)
    assert np.array_equal(phi, np.eye(6))
    r[0] = 12345.0
    assert r0[0] != 12345.0


@pytest.mark.parametrize(("r0", "v0", "dt"), [c[1:] for c in _REGIME_CASES], ids=_REGIME_IDS)
def test_state_matches_existing_propagator(r0: Vec3, v0: Vec3, dt: float) -> None:
    """The propagated state agrees with ``kepler.propagate`` to numerical noise."""
    r, v, _ = shepperd_stm(r0, v0, dt)
    r_ref, v_ref = propagate(r0, v0, dt)
    assert float(np.linalg.norm(r - r_ref)) / float(np.linalg.norm(r_ref)) < 1.0e-12
    assert float(np.linalg.norm(v - v_ref)) / float(np.linalg.norm(v_ref)) < 1.0e-12


@pytest.mark.parametrize(("r0", "v0", "dt"), [c[1:] for c in _REGIME_CASES], ids=_REGIME_IDS)
def test_fd_consistency(r0: Vec3, v0: Vec3, dt: float) -> None:
    """Analytic STM matches central differences of the independent propagator."""
    _, _, phi = shepperd_stm(r0, v0, dt)
    phi_fd = _fd_stm(r0, v0, dt)
    assert _max_block_rel_err(phi, phi_fd) < _FD_RTOL


def test_fd_consistency_random_elliptic() -> None:
    """Seeded random planar ellipses (a, e, nu, argp, dt) all pass the FD check."""
    rng = np.random.default_rng(20260613)
    for _ in range(6):
        a_km = float(rng.uniform(0.4, 4.0)) * AU_KM
        e = float(rng.uniform(0.0, 0.9))
        nu = float(rng.uniform(0.0, 2.0 * pi))
        argp = float(rng.uniform(0.0, 2.0 * pi))
        period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
        dt = float(rng.uniform(-1.0, 1.0)) * period_s
        r0, v0 = coe_to_rv(a_km, e, nu, arg_peri_rad=argp)
        _, _, phi = shepperd_stm(r0, v0, dt)
        phi_fd = _fd_stm(r0, v0, dt)
        assert _max_block_rel_err(phi, phi_fd) < _FD_RTOL, (a_km, e, nu, argp, dt)


@pytest.mark.parametrize(("r0", "v0", "dt"), [c[1:] for c in _REGIME_CASES], ids=_REGIME_IDS)
def test_symplectic(r0: Vec3, v0: Vec3, dt: float) -> None:
    """Phi^T J Phi = J with the canonical J (two-body flow is symplectic)."""
    _, _, phi = shepperd_stm(r0, v0, dt)
    phi_nd = _nondimensionalise(phi, float(np.linalg.norm(r0)))
    residual = phi_nd.T @ _J_CANONICAL @ phi_nd - _J_CANONICAL
    rel = float(np.linalg.norm(residual)) / float(np.linalg.norm(_J_CANONICAL))
    assert rel < 1.0e-10


@pytest.mark.parametrize(("r0", "v0", "dt"), [c[1:] for c in _REGIME_CASES], ids=_REGIME_IDS)
def test_det_one(r0: Vec3, v0: Vec3, dt: float) -> None:
    """det(Phi) = 1 (Liouville; implied by symplecticity, checked directly)."""
    _, _, phi = shepperd_stm(r0, v0, dt)
    phi_nd = _nondimensionalise(phi, float(np.linalg.norm(r0)))
    assert abs(float(np.linalg.det(phi_nd)) - 1.0) < 1.0e-10


@pytest.mark.parametrize(("r0", "v0", "dt"), [c[1:] for c in _REGIME_CASES], ids=_REGIME_IDS)
def test_composition(r0: Vec3, v0: Vec3, dt: float) -> None:
    """Phi(t2, t0) = Phi(t2, t1) Phi(t1, t0), with t1 an interior point."""
    t1 = 0.4 * dt
    r1, v1, phi_10 = shepperd_stm(r0, v0, t1)
    _, _, phi_21 = shepperd_stm(r1, v1, dt - t1)
    _, _, phi_20 = shepperd_stm(r0, v0, dt)
    assert _max_block_rel_err(phi_21 @ phi_10, phi_20) < 1.0e-12
