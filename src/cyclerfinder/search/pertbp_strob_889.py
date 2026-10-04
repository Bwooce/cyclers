"""Published positive control for the torus-connection lane: the planar elliptic RTBP.

#889. Reproduces (or fails to reproduce, and says where) the refined heteroclinic
connection between the Jupiter-Europa 3:4 and 5:6 resonant tori printed in Kumar,
Anderson & de la Llave, SIAM J. Applied Dynamical Systems 24(1):219-258 (2025),
filed in the private paper corpus as
``kumar-anderson-delallave-2025-gpu-connections-tori-perturbed-crtbp-siam-ads-24-219-arxiv-2109.14814.pdf``
(here "SIAM"), with the model and torus machinery of the companion Kumar, Anderson &
de la Llave, Celest. Mech. Dyn. Astron. 134:3 (2022), filed as
``kumar-anderson-delallave-2021-whiskered-tori-manifolds-cmda-arxiv-2105.11100.pdf``
("CMDA"), and the periodic orbits of Kumar, Anderson & de la Llave (2021),
arXiv:2109.14800, filed as
``kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf``
("CNSNS").

Model (SIAM Eq. 3.5 = CMDA Eq. 5, verbatim)
--------------------------------------------
``H(x, y, px, py, t) = (px^2 + py^2)/2 + n(t) (px y - py x) - (1-mu)/r1 - mu/r2`` with
``r1 = |(x + mu rho(t), y)|``, ``r2 = |(x - (1-mu) rho(t), y)|``, ``rho = 1 - e cos E(t)``,
``E - e sin E = t`` (periapse at ``t = 0``) and ``n(t) = sqrt(1-e^2)/rho^2`` (the rate of
the true anomaly). Coordinates: barycentric, rotating with the primaries' line (NOT
pulsating), ``(px, py) = (xdot - n y, ydot + n x)`` canonical momenta, time ``t`` the
independent variable, forcing period ``2 pi``. Jupiter-Europa: ``mu ~ 2.527e-5``,
``e = 0.0094``.

Stroboscopic map ``F``: the time-``2 pi`` flow map started at ``t = 0 mod 2 pi`` (CMDA
Sec. 3.1 fixes ``theta_p = 0``). An invariant circle satisfies ``F(K(theta)) =
K(theta + omega)`` with ``omega = 2 pi Omega_1 / Omega_p = 4 pi^2 / T1`` (CMDA Sec. 4.9,
``T1`` the parent periodic orbit's period). Bundles ``DF(K) v = lam v(theta + omega)``
with CONSTANT ``lam`` (CMDA Sec. 4.8); the manifold parameterisation ``W(theta, s)``
satisfies ``F(W(theta, s)) = W(theta + omega, lam s)`` (SIAM Eq. 4.6) and is normalised
as in SIAM Sec. 8.2-8.3: ``|W_1(0)| = 1`` and ``K(0)`` on the negative x-axis.

What is reused from :mod:`cyclerfinder.search.ccr4bp_strob_connection` (#882): the
odd-node trigonometric interpolation helpers (``node_angles``, ``fourier_eval``,
``shift_matrix``, ``fourier_tail``) and the PCRTBP symmetric periodic-orbit
corrector ``symmetric_periodic_orbit`` (whose model is exactly the PCRTBP when the
perturber mass is zero). Everything that integrates is rewritten here because the
#882 functions take a ``CCR4BPSystem`` (four-body, velocity coordinates, constant
rotation) and cannot carry the elliptic forcing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numba
import numpy as np
from numba import njit, prange
from numpy.typing import NDArray
from scipy.integrate._ivp import dop853_coefficients as _dc

from cyclerfinder.search.ccr4bp_strob_connection import (
    fourier_eval,
    fourier_tail,
    node_angles,
    shift_matrix,
    symmetric_periodic_orbit,
)

FloatArray = NDArray[np.float64]

TWO_PI = 2.0 * math.pi

# ---------------------------------------------------------------------------
# Printed data (golden values, with source).
# ---------------------------------------------------------------------------

#: Jupiter-Europa mass ratio printed in CNSNS Sec. 3 ("we use mu_E = ...").
MU_JE_CNSNS = 2.5266448850435028e-5
#: Europa mass ratio printed in Kumar et al. AAS 21-651 (arXiv:2109.14815) Sec. 5.
MU_JE_AAS21651 = 2.5265115494603433e-5
#: Europa's eccentricity used in SIAM Sec. 3.3 and CMDA Sec. 4.12.
ECC_EUROPA = 0.0094

#: CNSNS Table 1 (C = 3.0024): (x, y, xdot, ydot), period, lam_s, lam_u.
TABLE1_5_6 = {
    "x": -1.231240907544348,
    "ydot": 0.371411618064504,
    "T": 38.328135171743014,
    "lam_s": 0.001256465177783,
    "lam_u": 795.8835769446018,
}
TABLE1_3_4 = {
    "x": -1.391929713356257,
    "ydot": 0.609863420586548,
    "T": 25.338526603095760,
    "lam_s": 0.011341070996024,
    "lam_u": 88.175093899915780,
}
TABLE1_JACOBI = 3.0024

#: SIAM Sec. 8.3, verbatim.
PRINTED_POINT = np.array([-0.96064, 0.88783, -0.51377, -0.64714])
PRINTED_THETA_U = 2.39703
PRINTED_S_U = -77.73428
PRINTED_THETA_S = 1.83093
PRINTED_S_S = 202.62277
PRINTED_OMEGA_U = 1.558039
PRINTED_OMEGA_S = 1.030011

#: Reversing involution of the PERTBP at t = 0 (CMDA Sec. 5.3): M = diag(1, -1, -1, 1).
REVERSOR = np.array([1.0, -1.0, -1.0, 1.0])

_J4 = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])

# DOP853 tableau (Hairer's coefficients, as shipped by scipy).
_A = np.ascontiguousarray(_dc.A[: _dc.N_STAGES, : _dc.N_STAGES], dtype=np.float64)
_B = np.ascontiguousarray(_dc.B, dtype=np.float64)
_C = np.ascontiguousarray(_dc.C[: _dc.N_STAGES], dtype=np.float64)
_E3 = np.ascontiguousarray(_dc.E3, dtype=np.float64)
_E5 = np.ascontiguousarray(_dc.E5, dtype=np.float64)


def omega_from_period(period: float) -> float:
    """Rotation number of the stroboscopic circle born from a period-``T1`` orbit.

    CMDA Sec. 4.9: ``omega = 4 pi^2 / (T1 Omega_p)`` with ``Omega_p = 1``.
    """
    return (4.0 * math.pi**2 / period) % TWO_PI


def period_from_omega(omega: float) -> float:
    """Inverse of :func:`omega_from_period` on the branch ``omega = 4 pi^2 / T1``."""
    return 4.0 * math.pi**2 / omega


# ---------------------------------------------------------------------------
# Equations of motion (numba).
# ---------------------------------------------------------------------------


@njit(cache=True)  # type: ignore[untyped-decorator]
def _n_rho(t: float, e: float) -> tuple[float, float]:
    """``(n(t), rho(t))``: true-anomaly rate and primaries' separation."""
    if e == 0.0:
        return 1.0, 1.0
    m = t - TWO_PI * math.floor(t / TWO_PI)
    ecc = m + e * math.sin(m)
    for _ in range(50):
        d = (ecc - e * math.sin(ecc) - m) / (1.0 - e * math.cos(ecc))
        ecc -= d
        if abs(d) < 1e-16:
            break
    rho = 1.0 - e * math.cos(ecc)
    return math.sqrt(1.0 - e * e) / (rho * rho), rho


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs(t: float, z: FloatArray, mu: float, e: float, out: FloatArray, stm: bool) -> None:
    nn, rho = _n_rho(t, e)
    x = z[0]
    y = z[1]
    px = z[2]
    py = z[3]
    om = 1.0 - mu
    dx1 = x + mu * rho
    dx2 = x - om * rho
    r1sq = dx1 * dx1 + y * y
    r2sq = dx2 * dx2 + y * y
    r1 = math.sqrt(r1sq)
    r2 = math.sqrt(r2sq)
    i13 = 1.0 / (r1sq * r1)
    i23 = 1.0 / (r2sq * r2)
    gx = om * dx1 * i13 + mu * dx2 * i23
    gy = om * y * i13 + mu * y * i23
    out[0] = px + nn * y
    out[1] = py - nn * x
    out[2] = nn * py - gx
    out[3] = -nn * px - gy
    if stm:
        i15 = i13 / r1sq
        i25 = i23 / r2sq
        gxx = om * (i13 - 3.0 * dx1 * dx1 * i15) + mu * (i23 - 3.0 * dx2 * dx2 * i25)
        gyy = om * (i13 - 3.0 * y * y * i15) + mu * (i23 - 3.0 * y * y * i25)
        gxy = -3.0 * om * dx1 * y * i15 - 3.0 * mu * dx2 * y * i25
        # A = [[0, n, 1, 0], [-n, 0, 0, 1], [-gxx, -gxy, 0, n], [-gxy, -gyy, -n, 0]]
        for j in range(4):
            p0 = z[4 + j]
            p1 = z[8 + j]
            p2 = z[12 + j]
            p3 = z[16 + j]
            out[4 + j] = nn * p1 + p2
            out[8 + j] = -nn * p0 + p3
            out[12 + j] = -gxx * p0 - gxy * p1 + nn * p3
            out[16 + j] = -gxy * p0 - gyy * p1 - nn * p2


@njit(cache=True)  # type: ignore[untyped-decorator]
def _dop853(
    z0: FloatArray,
    t0: float,
    t1: float,
    mu: float,
    e: float,
    stm: bool,
    rtol: float,
    atol: float,
    max_steps: int,
    a: FloatArray,
    b: FloatArray,
    c: FloatArray,
    e3: FloatArray,
    e5: FloatArray,
) -> tuple[FloatArray, int]:
    """Single-trajectory DOP853 (a port of scipy's step-size control, no dense output)."""
    n = z0.size
    y = z0.copy()
    t = t0
    if t1 == t0:
        return y, 0
    direction = 1.0 if t1 > t0 else -1.0
    k = np.empty((13, n))
    f = np.empty(n)
    _rhs(t, y, mu, e, f, stm)
    ytmp = np.empty(n)
    ynew = np.empty(n)
    h_abs = min(0.01, abs(t1 - t0))
    nsteps = 0
    expo = -1.0 / 8.0
    while direction * (t1 - t) > 0.0:
        if nsteps >= max_steps:
            return y, -1
        rejected = False
        while True:
            if h_abs < 1e-13:
                return y, -2
            h = h_abs * direction
            t_new = t + h
            if direction * (t_new - t1) > 0.0:
                t_new = t1
            h = t_new - t
            h_abs = abs(h)
            for i in range(n):
                k[0, i] = f[i]
            for s in range(1, 12):
                for i in range(n):
                    acc = 0.0
                    for j in range(s):
                        acc += a[s, j] * k[j, i]
                    ytmp[i] = y[i] + h * acc
                _rhs(t + c[s] * h, ytmp, mu, e, k[s], stm)
            for i in range(n):
                acc = 0.0
                for j in range(12):
                    acc += b[j] * k[j, i]
                ynew[i] = y[i] + h * acc
            _rhs(t + h, ynew, mu, e, k[12], stm)
            err5 = 0.0
            err3 = 0.0
            for i in range(n):
                sc = atol + max(abs(y[i]), abs(ynew[i])) * rtol
                s5 = 0.0
                s3 = 0.0
                for j in range(13):
                    s5 += k[j, i] * e5[j]
                    s3 += k[j, i] * e3[j]
                err5 += (s5 / sc) ** 2
                err3 += (s3 / sc) ** 2
            if err5 == 0.0 and err3 == 0.0:
                en = 0.0
            else:
                en = h_abs * err5 / math.sqrt((err5 + 0.01 * err3) * n)
            if en < 1.0:
                factor = 10.0 if en == 0.0 else min(10.0, 0.9 * en**expo)
                if rejected:
                    factor = min(1.0, factor)
                h_abs *= factor
                break
            h_abs *= max(0.2, 0.9 * en**expo)
            rejected = True
        t = t_new
        for i in range(n):
            y[i] = ynew[i]
            f[i] = k[12, i]
        nsteps += 1
    return y, nsteps


@njit(cache=True, parallel=True)  # type: ignore[untyped-decorator]
def _batch(
    z0s: FloatArray,
    t0: float,
    t1: float,
    mu: float,
    e: float,
    stm: bool,
    rtol: float,
    atol: float,
    max_steps: int,
    a: FloatArray,
    b: FloatArray,
    c: FloatArray,
    e3: FloatArray,
    e5: FloatArray,
) -> tuple[FloatArray, NDArray[np.int64]]:
    m = z0s.shape[0]
    out = np.empty_like(z0s)
    status = np.empty(m, dtype=np.int64)
    for i in prange(m):
        yi, st = _dop853(z0s[i], t0, t1, mu, e, stm, rtol, atol, max_steps, a, b, c, e3, e5)
        out[i] = yi
        status[i] = st
    return out, status


def set_threads(n: int = 4) -> None:
    """Cap numba's worker threads (the machine is shared; the task allows 4)."""
    numba.set_num_threads(max(1, min(n, numba.config.NUMBA_NUM_THREADS)))


def rhs(t: float, state: FloatArray, mu: float, e: float) -> FloatArray:
    """Hamiltonian vector field ``(dx, dy, dpx, dpy)/dt`` (SIAM Eq. 3.3 with 3.5)."""
    out = np.empty(4)
    _rhs(float(t), np.asarray(state, dtype=np.float64), float(mu), float(e), out, False)
    return out


def flow(
    states: FloatArray,
    t0: float,
    t1: float,
    mu: float,
    e: float,
    *,
    with_stm: bool = False,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    max_steps: int = 2_000_000,
) -> tuple[FloatArray, FloatArray | None]:
    """Flow ``(m, 4)`` (or ``(4,)``) states from ``t0`` to ``t1``; optional 4x4 STMs."""
    s = np.asarray(states, dtype=np.float64)
    single = s.ndim == 1
    s2 = np.atleast_2d(s)
    m = s2.shape[0]
    eye = np.tile(np.eye(4).reshape(-1), (m, 1))
    z0 = np.hstack([s2, eye]) if with_stm else s2.copy()
    out, status = _batch(
        np.ascontiguousarray(z0),
        float(t0),
        float(t1),
        float(mu),
        float(e),
        bool(with_stm),
        float(rtol),
        float(atol),
        int(max_steps),
        _A,
        _B,
        _C,
        _E3,
        _E5,
    )
    if np.any(status < 0):
        raise RuntimeError(f"integration failed for {int(np.sum(status < 0))} of {m} states")
    st = out[:, :4].copy()
    phi = out[:, 4:].reshape(m, 4, 4).copy() if with_stm else None
    if single:
        return st[0], (phi[0] if phi is not None else None)
    return st, phi


def strob(
    states: FloatArray,
    n: int,
    mu: float,
    e: float,
    *,
    with_stm: bool = False,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> tuple[FloatArray, FloatArray | None]:
    """``F^n`` of the stroboscopic map at phase ``t = 0`` (``n < 0`` maps backwards)."""
    return flow(states, 0.0, TWO_PI * n, mu, e, with_stm=with_stm, rtol=rtol, atol=atol)


def n_rho(t: float, e: float) -> tuple[float, float]:
    """Python access to ``(n(t), rho(t))``."""
    nn, rho = _n_rho(float(t), float(e))
    return float(nn), float(rho)


def momentum_to_velocity(t: float, state: FloatArray, e: float) -> FloatArray:
    """``(x, y, px, py) -> (x, y, xdot, ydot)`` with ``xdot = px + n y``, ``ydot = py - n x``."""
    s = np.asarray(state, dtype=np.float64)
    nn, _ = n_rho(t, e)
    out = s.copy()
    out[..., 2] = s[..., 2] + nn * s[..., 1]
    out[..., 3] = s[..., 3] - nn * s[..., 0]
    return out


def velocity_to_momentum(t: float, state: FloatArray, e: float) -> FloatArray:
    """Inverse of :func:`momentum_to_velocity`."""
    s = np.asarray(state, dtype=np.float64)
    nn, _ = n_rho(t, e)
    out = s.copy()
    out[..., 2] = s[..., 2] - nn * s[..., 1]
    out[..., 3] = s[..., 3] + nn * s[..., 0]
    return out


def hamiltonian(t: float, state: FloatArray, mu: float, e: float) -> FloatArray:
    """SIAM Eq. 3.5 evaluated on ``(..., 4)`` momentum states."""
    s = np.asarray(state, dtype=np.float64)
    nn, rho = n_rho(t, e)
    x, y, px, py = s[..., 0], s[..., 1], s[..., 2], s[..., 3]
    r1 = np.hypot(x + mu * rho, y)
    r2 = np.hypot(x - (1.0 - mu) * rho, y)
    return np.asarray(
        0.5 * (px * px + py * py) + nn * (px * y - py * x) - (1.0 - mu) / r1 - mu / r2
    )


def jacobi(state: FloatArray, mu: float) -> FloatArray:
    """PCRTBP Jacobi constant ``C = -2 H0`` (CNSNS Eq. 3; no ``mu(1-mu)`` term)."""
    return np.asarray(-2.0 * hamiltonian(0.0, state, mu, 0.0))


# ---------------------------------------------------------------------------
# PCRTBP resonant periodic orbits.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PeriodicOrbit:
    """x-axis-symmetric PCRTBP orbit; ``state`` in MOMENTUM coordinates at ``t = 0``."""

    mu: float
    state: FloatArray
    period: float
    lam_s: float
    lam_u: float
    jacobi: float
    residual: float

    @property
    def omega(self) -> float:
        return omega_from_period(self.period)


def _monodromy(mu: float, state_mom: FloatArray, period: float) -> FloatArray:
    _, phi = flow(state_mom, 0.0, period, mu, 0.0, with_stm=True)
    assert phi is not None
    return phi


def _po_from_velocity(mu: float, x0: float, vy: float, period: float, res: float) -> PeriodicOrbit:
    sv = np.array([x0, 0.0, 0.0, vy])
    sm = velocity_to_momentum(0.0, sv, 0.0)
    mono = _monodromy(mu, sm, period)
    ev = np.linalg.eigvals(mono)
    real = np.sort(np.abs(ev.real[np.abs(ev.imag) < 1e-6]))
    return PeriodicOrbit(
        mu=mu,
        state=sm,
        period=period,
        lam_s=float(real[0]),
        lam_u=float(real[-1]),
        jacobi=float(jacobi(sm, mu)),
        residual=res,
    )


def po_at_fixed_x(mu: float, x0: float, vy_guess: float, period_guess: float) -> PeriodicOrbit:
    """Symmetric orbit through ``(x0, 0)`` (reuses the #882 PCRTBP corrector)."""
    s4, period, res = symmetric_periodic_orbit(mu, x0, vy_guess, 0.5 * period_guess, tol=1e-13)
    return _po_from_velocity(mu, float(s4[0]), float(s4[3]), float(period), float(res))


def po_at_fixed_period(
    mu: float, x_guess: float, vy_guess: float, period: float, *, tol: float = 1e-13
) -> PeriodicOrbit:
    """Symmetric orbit of prescribed period: unknowns ``(x0, vy0)``, conditions
    ``y(T/2) = vx(T/2) = 0`` (velocity frame; momentum ``px = vx - y`` at ``y = 0``)."""
    x0 = float(x_guess)
    vy = float(vy_guess)
    res = float("inf")
    for _ in range(40):
        sm = velocity_to_momentum(0.0, np.array([x0, 0.0, 0.0, vy]), 0.0)
        end, phi = flow(sm, 0.0, 0.5 * period, mu, 0.0, with_stm=True)
        assert phi is not None
        # d(state_mom)/d(x0, vy): px = vx - y = 0, py = vy + x
        dmom = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
        dend = phi @ dmom
        # conditions: y = 0 and vx = px + y = 0
        g = np.array([end[1], end[2] + end[1]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([dend[1], dend[2] + dend[1]])
        dz = np.linalg.solve(jac, -g)
        x0 += float(dz[0])
        vy += float(dz[1])
    return _po_from_velocity(mu, x0, vy, period, res)


# ---------------------------------------------------------------------------
# Invariant circles of the stroboscopic map.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Circle:
    """``F(K(theta)) = K(theta + omega)`` for the map at ``t = 0``; odd node count."""

    mu: float
    e: float
    omega: float
    nodes: FloatArray
    residual: float = float("nan")
    converged: bool = False
    history: tuple[float, ...] = ()

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    def state(self, theta: float | FloatArray) -> FloatArray:
        return fourier_eval(self.nodes, theta)

    def dstate(self, theta: float | FloatArray) -> FloatArray:
        return fourier_eval(self.nodes, theta, deriv=1)

    def tail(self) -> float:
        return fourier_tail(self.nodes)


def seed_circle(po: PeriodicOrbit, n_nodes: int) -> Circle:
    """``K(theta) = phi_PCRTBP(x0, T theta / 2 pi)`` (CMDA Sec. 4.9); exact at ``e = 0``.

    ``K(0) = x0`` is the perpendicular negative-x-axis crossing, the SIAM Sec. 8.3
    phase convention ("theta = 0 points ... on the negative x-axis").
    """
    if n_nodes % 2 != 1:
        raise ValueError("odd node count required")
    times = node_angles(n_nodes) * po.period / TWO_PI
    nodes = np.empty((n_nodes, 4))
    nodes[0] = po.state
    cur = po.state.copy()
    for j in range(1, n_nodes):
        cur, _ = flow(cur, float(times[j - 1]), float(times[j]), po.mu, 0.0)
        nodes[j] = cur
    return Circle(mu=po.mu, e=0.0, omega=po.omega, nodes=nodes)


def circle_residual(circle: Circle) -> tuple[float, float]:
    """(on-node, off-node) invariance residuals ``max |F(K(th)) - K(th + omega)|``."""
    n = circle.n_nodes
    fu, _ = strob(circle.nodes, 1, circle.mu, circle.e)
    on = float(np.max(np.abs(shift_matrix(n, -circle.omega) @ fu - circle.nodes)))
    th = node_angles(n) + math.pi / n
    fm, _ = strob(circle.state(th), 1, circle.mu, circle.e)
    off = float(np.max(np.abs(fm - circle.state(th + circle.omega))))
    return on, off


def correct_circle(
    circle: Circle,
    e: float,
    *,
    tol: float = 1e-11,
    max_iter: int = 12,
    verbose: bool = False,
) -> Circle:
    """Gauss-Newton on the node values at FIXED ``omega`` (the paper's continuation
    is also at fixed ``omega``, CMDA Sec. 4.9). Phase condition: ``y(K(0)) = 0``."""
    nodes = circle.nodes.copy()
    n = circle.n_nodes
    smat = shift_matrix(n, -circle.omega)
    prow = np.zeros(4 * n)
    prow[1] = 1.0
    hist: list[float] = []
    # Bordering column: the symplectic conjugate of the tangent, J^-1 DK. The fixed-
    # omega invariance operator has DK in its kernel (removed by the phase row) and a
    # one-dimensional cokernel (the averaged action equation, cf. CMDA Sec. 4.5), which
    # J^-1 DK spans and DK does not. Its multiplier is a translated-torus counterterm
    # that must vanish (quadratically) at a true solution; it is reported in ``history``.
    dk = fourier_eval(nodes, node_angles(n), deriv=1)
    tcol = np.column_stack([-dk[:, 2], -dk[:, 3], dk[:, 0], dk[:, 1]]).reshape(-1)
    tcol /= np.linalg.norm(tcol)

    def evaluate(u: FloatArray) -> tuple[FloatArray, FloatArray]:
        fu, dfs = strob(u, 1, circle.mu, e, with_stm=True)
        assert dfs is not None
        return smat @ fu - u, dfs

    g, dfs = evaluate(nodes)
    res = float(np.max(np.abs(g)))
    hist.append(res)
    it = 0
    while res > tol and it < max_iter:
        it += 1
        a = np.zeros((4 * n + 1, 4 * n + 1))
        a[: 4 * n, : 4 * n] = np.einsum("jl,lab->jalb", smat, dfs).reshape(4 * n, 4 * n)
        a[: 4 * n, : 4 * n] -= np.eye(4 * n)
        a[: 4 * n, -1] = tcol
        a[-1, : 4 * n] = prow
        rhs_v = -np.concatenate([g.reshape(-1), [nodes[0, 1]]])
        sol = np.linalg.solve(a, rhs_v)
        step = sol[: 4 * n].reshape(n, 4)
        counterterm = float(sol[-1])
        alpha = 1.0
        for _ in range(6):
            trial = nodes + alpha * step
            g_t, dfs_t = evaluate(trial)
            res_t = float(np.max(np.abs(g_t)))
            if res_t <= 2.0 * res:
                break
            alpha *= 0.5
        nodes, g, dfs, res = trial, g_t, dfs_t, res_t
        hist.append(res)
        if verbose:
            print(
                f"    circle GN {it}: residual {res:.3e} alpha {alpha:g} "
                f"counterterm {counterterm:.2e}",
                flush=True,
            )
        if len(hist) >= 4 and res > 0.5 * hist[-4]:
            break
    return Circle(
        mu=circle.mu,
        e=e,
        omega=circle.omega,
        nodes=nodes,
        residual=res,
        converged=res <= tol,
        history=tuple(hist),
    )


def continue_circle_in_e(
    circle: Circle,
    e_target: float,
    *,
    first_step: float = 2e-5,
    max_step: float = 1e-3,
    min_step: float = 1e-8,
    tol: float = 1e-11,
    verbose: bool = False,
) -> list[Circle]:
    """Adaptive continuation in eccentricity at fixed ``omega`` (secant predictor;
    step halved on failure, grown by 1.5 after an easy step). Returns the accepted
    chain; the last element is at ``e_target`` iff the continuation succeeded."""
    out: list[Circle] = [circle]
    h = first_step
    while out[-1].e < e_target - 1e-15:
        cur = out[-1]
        e_new = min(e_target, cur.e + h)
        nodes = cur.nodes
        if len(out) >= 2:
            prev = out[-2]
            r = (e_new - cur.e) / (cur.e - prev.e)
            nodes = cur.nodes + r * (cur.nodes - prev.nodes)
        guess = Circle(mu=cur.mu, e=cur.e, omega=cur.omega, nodes=nodes)
        nxt = correct_circle(guess, e_new, tol=tol, max_iter=8)
        if nxt.converged:
            out.append(nxt)
            if verbose:
                print(
                    f"  e = {e_new:.7f} (h {h:.2e}): {len(nxt.history) - 1} it, "
                    f"residual {nxt.residual:.2e}",
                    flush=True,
                )
            if len(nxt.history) - 1 <= 4:
                h = min(max_step, 1.5 * h)
        else:
            h *= 0.5
            if h < min_step:
                if verbose:
                    print(f"  continuation stalled at e = {cur.e:.7f}", flush=True)
                break
    return out


def reversibility_defect(circle: Circle, n_test: int = 257) -> float:
    """``max_theta |M K(-theta) - K(theta)|``: the circle is symmetric under the
    PERTBP reversor at ``t = 0`` (CMDA Sec. 5.3) iff this vanishes."""
    th = TWO_PI * np.arange(n_test) / n_test
    return float(np.max(np.abs(REVERSOR * circle.state(-th) - circle.state(th))))


def distance_to_circle(
    circle: Circle, state: FloatArray, n_fine: int = 4096
) -> tuple[float, float]:
    """``min_theta |state - K(theta)|`` (grid then golden-section) and its argmin."""
    s = np.asarray(state, dtype=np.float64)
    th = TWO_PI * np.arange(n_fine) / n_fine
    d = np.linalg.norm(circle.state(th) - s[None, :], axis=1)
    i = int(np.argmin(d))
    lo, hi = th[i] - TWO_PI / n_fine, th[i] + TWO_PI / n_fine
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(60):
        c1 = hi - gr * (hi - lo)
        c2 = lo + gr * (hi - lo)
        if np.linalg.norm(circle.state(c1) - s) < np.linalg.norm(circle.state(c2) - s):
            hi = c2
        else:
            lo = c1
    tm = 0.5 * (lo + hi)
    dm = float(np.linalg.norm(circle.state(tm) - s))
    return min(dm, float(d[i])), float(tm % TWO_PI)


# ---------------------------------------------------------------------------
# Bundles (CMDA Sec. 4.8 and 4.11) and second-order manifold terms.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Bundles:
    """Constant-multiplier bundles, ``|v(0)| = 1`` and ``v_x(0) > 0`` (pre-registered
    from SIAM Fig. 11, whose x-component of the unit-normalised 3:4 unstable bundle
    at ``theta = 0`` is positive, about 0.116). ``w_*`` are the ``s^2`` coefficients."""

    lam_u: float
    lam_s: float
    v_u: FloatArray
    v_s: FloatArray
    w_u: FloatArray
    w_s: FloatArray
    offgrid_u: float
    offgrid_s: float
    tail_u: float
    tail_s: float

    def lam(self, branch: str) -> float:
        return self.lam_u if branch == "unstable" else self.lam_s

    def v(self, branch: str) -> FloatArray:
        return self.v_u if branch == "unstable" else self.v_s

    def w(self, branch: str) -> FloatArray:
        return self.w_u if branch == "unstable" else self.w_s


def _wavenumbers(n: int) -> FloatArray:
    return np.asarray(np.rint(np.fft.fftfreq(n, 1.0 / n)), dtype=np.float64)


def _const_multiplier(
    v: FloatArray, dfs: FloatArray, omega: float
) -> tuple[FloatArray, float, float]:
    """Rescale a line field so that ``DF v = lam v(theta + omega)`` with constant
    ``lam`` (CMDA Eq. 51-52): solve ``log a(th + w) - log a(th) = log lam(th) - mean``."""
    n = v.shape[0]
    vp = shift_matrix(n, omega) @ v
    g = np.einsum("jab,jb->ja", dfs, v)
    lam_th = np.sum(g * vp, axis=1) / np.sum(vp * vp, axis=1)
    pointwise = float(np.max(np.abs(g - lam_th[:, None] * vp)))
    sign = float(np.sign(lam_th[0]))
    loglam = np.log(np.abs(lam_th))
    mean = float(loglam.mean())
    k = _wavenumbers(n)
    bh = np.fft.fft(loglam - mean)
    den = np.exp(1j * k * omega) - 1.0
    ah = np.zeros_like(bh)
    nz = k != 0
    ah[nz] = bh[nz] / den[nz]
    loga = np.real(np.fft.ifft(ah))
    return v * np.exp(loga)[:, None], sign * math.exp(mean), pointwise


def _power_line_field(
    dfs: FloatArray, omega: float, unstable: bool, n_iter: int = 200, tol: float = 1e-15
) -> FloatArray:
    """Projective power iteration for the line field (CMDA Eq. 76-77)."""
    n = dfs.shape[0]
    sm = shift_matrix(n, -omega)
    sp = shift_matrix(n, omega)
    inv = np.linalg.inv(dfs)

    def continuous(u: FloatArray) -> FloatArray:
        # orient each node like its predecessor so the field is a smooth function
        # of theta (the bundles are cylinders: positive multipliers, CMDA Sec. 4.9)
        u = u / np.linalg.norm(u, axis=1)[:, None]
        for j in range(1, u.shape[0]):
            if float(u[j] @ u[j - 1]) < 0.0:
                u[j] = -u[j]
        return u

    v = continuous(np.einsum("jab,b->ja", dfs if unstable else inv, np.array([1.0, 0.7, 0.4, 0.2])))
    for _ in range(n_iter):
        if unstable:
            # DF(K(th_j)) v(th_j) is the direction at th_j + omega: normalise it there
            # (a smooth unit field) BEFORE interpolating back onto the nodes.
            nv = sm @ continuous(np.einsum("jab,jb->ja", dfs, v))
        else:
            nv = np.einsum("jab,jb->ja", inv, sp @ v)
        nv = continuous(nv)
        if float(np.sum(nv * v)) < 0.0:
            nv = -nv
        d = float(np.max(np.abs(nv - v)))
        v = nv
        if d < tol:
            break
    return v


def hyperbolic_bundles(circle: Circle, *, fd_h: float = 1e-4) -> Bundles:
    """Bundles by projective power iteration, constant-multiplier rescaling, the SIAM
    Sec. 8.2 unit normalisation at ``theta = 0``, and FD second-order terms."""
    n = circle.n_nodes
    f0, dfs = strob(circle.nodes, 1, circle.mu, circle.e, with_stm=True)
    assert dfs is not None
    smat = shift_matrix(n, -circle.omega)
    op = np.einsum("jl,lab->jalb", smat, dfs).reshape(4 * n, 4 * n)
    out: dict[str, tuple[float, FloatArray, FloatArray, float, float]] = {}
    th_mid = node_angles(n) + math.pi / n
    _, df_mid = strob(circle.state(th_mid), 1, circle.mu, circle.e, with_stm=True)
    assert df_mid is not None
    for branch in ("unstable", "stable"):
        v = _power_line_field(dfs, circle.omega, branch == "unstable")
        v, lam, _ = _const_multiplier(v, dfs, circle.omega)
        v = v / np.linalg.norm(v[0])
        if v[0, 0] < 0.0:
            v = -v
        lhs = np.einsum("jab,jb->ja", df_mid, fourier_eval(v, th_mid))
        off = float(np.max(np.abs(lhs - lam * fourier_eval(v, th_mid + circle.omega))))
        # second order: DF w(th) + q/2 = lam^2 w(th + omega), q = D^2F[v, v] by FD
        fp, _ = strob(circle.nodes + fd_h * v, 1, circle.mu, circle.e)
        fm, _ = strob(circle.nodes - fd_h * v, 1, circle.mu, circle.e)
        q = (fp + fm - 2.0 * f0) / (fd_h * fd_h)
        w = np.linalg.solve(op - lam * lam * np.eye(4 * n), -0.5 * (smat @ q).reshape(-1))
        out[branch] = (lam, v, w.reshape(n, 4), off, fourier_tail(v))
    lu, vu, wu, ou, tu = out["unstable"]
    ls, vs, ws, os_, ts = out["stable"]
    return Bundles(
        lam_u=float(lu),
        lam_s=float(ls),
        v_u=vu,
        v_s=vs,
        w_u=wu,
        w_s=ws,
        offgrid_u=float(ou),
        offgrid_s=float(os_),
        tail_u=float(tu),
        tail_s=float(ts),
    )


# ---------------------------------------------------------------------------
# Manifolds: local parameterisation and globalisation (SIAM Eq. 4.9-4.10, 7.2-7.7).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Whisker:
    """One manifold of one circle: ``W(theta, s)``, ``F(W(th, s)) = W(th + omega, lam s)``."""

    circle: Circle
    bundles: Bundles
    branch: str  # "unstable" | "stable"
    order: int = 2

    @property
    def lam(self) -> float:
        return self.bundles.lam(self.branch)

    def local(self, theta: float, sig: float) -> tuple[FloatArray, FloatArray, FloatArray]:
        """``W_loc``, ``d/dtheta``, ``d/dsig`` at one point."""
        c = self.circle
        v = self.bundles.v(self.branch)
        w = self.bundles.w(self.branch)
        p = c.state(theta) + sig * fourier_eval(v, theta)
        dth = c.dstate(theta) + sig * fourier_eval(v, theta, deriv=1)
        ds = fourier_eval(v, theta).copy()
        if self.order >= 2:
            p = p + sig * sig * fourier_eval(w, theta)
            dth = dth + sig * sig * fourier_eval(w, theta, deriv=1)
            ds = ds + 2.0 * sig * fourier_eval(w, theta)
        return p, dth, ds

    def n_maps(self, s: float, delta: float) -> int:
        """Smallest ``m >= 0`` with ``|s| * |lam|^-m <= delta`` (unstable) or
        ``|s| * |lam|^m <= delta`` (stable)."""
        a = abs(self.lam)
        if abs(s) <= delta:
            return 0
        k = math.log(abs(s) / delta) / abs(math.log(a))
        return max(0, math.ceil(k - 1e-12))

    def evaluate(
        self, theta: float, s: float, m: int, *, with_jac: bool = False
    ) -> tuple[FloatArray, FloatArray | None, float, float]:
        """``W(theta, s)`` by ``m`` maps from the local chart.

        Returns ``(point, 4x2 jacobian [d/dtheta, d/ds] or None, local theta, local s)``.
        """
        lam = self.lam
        om = self.circle.omega
        if self.branch == "unstable":
            th_l = theta - m * om
            sig = s * lam ** (-m)
            dsig = lam ** (-m)
            nmap = m
        else:
            th_l = theta + m * om
            sig = s * lam**m
            dsig = lam**m
            nmap = -m
        p, dth, ds = self.local(th_l, sig)
        if nmap == 0:
            jac = np.column_stack([dth, ds * dsig]) if with_jac else None
            return p, jac, th_l, sig
        x, phi = strob(p, nmap, self.circle.mu, self.circle.e, with_stm=with_jac)
        jac = None
        if with_jac:
            assert phi is not None
            jac = np.column_stack([phi @ dth, phi @ ds * dsig])
        return x, jac, th_l, sig


@dataclass(frozen=True)
class Connection:
    """A zero of ``W_u(theta_u, s_u) - W_s(theta_s, s_s)`` (SIAM Eq. 7.1)."""

    z: FloatArray  # (theta_u, s_u, theta_s, s_s)
    point: FloatArray
    residual: float
    m_u: int
    m_s: int
    delta: float
    singular_values: FloatArray
    null_vector: FloatArray
    n_iter: int
    history: tuple[float, ...]


def refine_connection(
    wu: Whisker,
    ws: Whisker,
    z0: FloatArray,
    *,
    delta: float,
    tol: float = 1e-12,
    max_iter: int = 60,
    verbose: bool = False,
) -> Connection:
    """Damped Gauss-Newton (as in SIAM Eq. 7.11, with step halving) at fixed map
    counts chosen from the seed and the local-chart radius ``delta``."""
    z = np.asarray(z0, dtype=np.float64).copy()
    m_u = wu.n_maps(z[1], delta)
    m_s = ws.n_maps(z[3], delta)

    def ev(zz: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
        pu, ju, _, _ = wu.evaluate(zz[0], zz[1], m_u, with_jac=True)
        ps, js, _, _ = ws.evaluate(zz[2], zz[3], m_s, with_jac=True)
        assert ju is not None and js is not None
        return pu - ps, np.hstack([ju, -js]), pu

    f, jac, pu = ev(z)
    res = float(np.linalg.norm(f))
    hist = [res]
    it = 0
    while res > tol and it < max_iter:
        it += 1
        step = np.linalg.lstsq(jac, -f, rcond=None)[0]
        alpha = 1.0
        for _ in range(30):
            zt = z + alpha * step
            ft, jt, put = ev(zt)
            rt = float(np.linalg.norm(ft))
            if rt < res:
                break
            alpha *= 0.5
        if rt >= res:
            break
        z, f, jac, pu, res = zt, ft, jt, put, rt
        hist.append(res)
        if verbose:
            print(f"    refine {it}: |f| = {res:.3e} alpha {alpha:g}", flush=True)
    cols = jac / np.linalg.norm(jac, axis=0)[None, :]
    _, sv, vt = np.linalg.svd(cols)
    return Connection(
        z=z,
        point=pu,
        residual=res,
        m_u=m_u,
        m_s=m_s,
        delta=delta,
        singular_values=sv,
        null_vector=vt[-1] / np.linalg.norm(jac, axis=0),
        n_iter=it,
        history=tuple(hist),
    )


@dataclass(frozen=True)
class ValleyPoint:
    """Least-squares minimum of ``|W_u - W_s|`` at FIXED ``theta_u`` (3 unknowns)."""

    z: FloatArray
    point: FloatArray
    mismatch: FloatArray  # W_u - W_s at the minimum (4-vector)
    norm: float


def valley_scan(
    wu: Whisker,
    ws: Whisker,
    z_start: FloatArray,
    thetas_u: FloatArray,
    *,
    m_u: int,
    m_s: int,
    max_iter: int = 12,
) -> list[ValleyPoint]:
    """Follow the near-degenerate intersection valley by fixing ``theta_u`` and
    minimising ``|W_u(theta_u, s_u) - W_s(theta_s, s_s)|`` over the other three
    parameters (warm-started along ``thetas_u``, which must start next to
    ``z_start[0]``). In the circular problem the intersection is a whole curve; the
    elliptic forcing leaves a small signed mismatch along it whose zeros are the
    isolated connections (a Melnikov-type function, measured here directly)."""
    out: list[ValleyPoint] = []
    z = np.asarray(z_start, dtype=np.float64).copy()
    for th in thetas_u:
        z[0] = float(th)
        best = None
        for _ in range(max_iter):
            pu, ju, _, _ = wu.evaluate(z[0], z[1], m_u, with_jac=True)
            ps, js, _, _ = ws.evaluate(z[2], z[3], m_s, with_jac=True)
            assert ju is not None and js is not None
            f = pu - ps
            nrm = float(np.linalg.norm(f))
            if best is None or nrm < best.norm:
                best = ValleyPoint(z=z.copy(), point=pu, mismatch=f, norm=nrm)
            jac = np.column_stack([ju[:, 1], -js[:, 0], -js[:, 1]])
            step = np.linalg.lstsq(jac, -f, rcond=None)[0]
            z[1:] = z[1:] + step
            if float(np.max(np.abs(step) / np.maximum(1.0, np.abs(z[1:])))) < 1e-13:
                break
        assert best is not None
        out.append(best)
        z = best.z.copy()
    return out


def closest_on_whisker(
    wh: Whisker, target: FloatArray, theta: float, s: float, m: int, *, max_iter: int = 60
) -> tuple[float, float, FloatArray, float]:
    """Damped Gauss-Newton for ``min_(theta, s) |W(theta, s) - target|`` at a fixed map
    count ``m``. Returns ``(theta, s, W, distance)``."""
    p = np.array([theta, s], dtype=np.float64)
    w, jac, _, _ = wh.evaluate(p[0], p[1], m, with_jac=True)
    assert jac is not None
    d = float(np.linalg.norm(w - target))
    for _ in range(max_iter):
        step = np.linalg.lstsq(jac, target - w, rcond=None)[0]
        alpha = 1.0
        improved = False
        for _ in range(30):
            pt = p + alpha * step
            wt, jt, _, _ = wh.evaluate(pt[0], pt[1], m, with_jac=True)
            dt = float(np.linalg.norm(wt - target))
            if dt < d:
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break
        assert jt is not None
        rel = float(np.max(np.abs(pt - p) / np.maximum(1.0, np.abs(p))))
        p, w, jac, d = pt, wt, jt, dt
        if rel < 1e-14:
            break
    return float(p[0]), float(p[1]), w, d


def chart_coordinates(
    wh: Whisker, state: FloatArray, *, max_iter: int = 30
) -> tuple[float, float, float]:
    """Least-squares ``(theta, sig)`` with ``W_loc(theta, sig) ~ state`` (2 unknowns,
    4 equations); returns ``(theta, sig, residual norm)``. Seeded by the nearest
    circle point and the symplectic hyperbolic coordinate."""
    _, th = distance_to_circle(wh.circle, state)
    c_s, c_u = hyperbolic_coordinates(wh.circle, wh.bundles, th, state)
    sig = c_u if wh.branch == "unstable" else c_s
    p = np.array([th, sig])
    for _ in range(max_iter):
        w, dth, ds = wh.local(p[0], p[1])
        r = w - state
        step = np.linalg.lstsq(np.column_stack([dth, ds]), -r, rcond=None)[0]
        p = p + step
        if float(np.max(np.abs(step))) < 1e-14:
            break
    w, _, _ = wh.local(p[0], p[1])
    return float(p[0]), float(p[1]), float(np.linalg.norm(w - state))


def seed_from_point(
    wu: Whisker, ws: Whisker, point: FloatArray, k_u: int, k_s: int
) -> tuple[FloatArray, dict[str, float]]:
    """Convention-free seed for :func:`refine_connection` from a phase-space point:
    map it ``k_u`` periods back and ``k_s`` forward, read off the local chart
    coordinates there, and convert to ``(theta_u, s_u, theta_s, s_s)`` at the point."""
    yb, _ = strob(point, -k_u, wu.circle.mu, wu.circle.e)
    yf, _ = strob(point, k_s, ws.circle.mu, ws.circle.e)
    th_b, sg_b, r_b = chart_coordinates(wu, yb)
    th_f, sg_f, r_f = chart_coordinates(ws, yf)
    z = np.array(
        [
            (th_b + k_u * wu.circle.omega) % TWO_PI,
            sg_b * wu.lam**k_u,
            (th_f - k_s * ws.circle.omega) % TWO_PI,
            sg_f * ws.lam ** (-k_s),
        ]
    )
    return z, {"chart_residual_u": r_b, "chart_residual_s": r_f, "sig_u": sg_b, "sig_s": sg_f}


# ---------------------------------------------------------------------------
# Verification criteria (pre-registered by the #882 adversarial review, Sec. 4).
# ---------------------------------------------------------------------------


def hyperbolic_coordinates(
    circle: Circle, bundles: Bundles, theta: float, state: FloatArray
) -> tuple[float, float]:
    """``(c_s, c_u)`` of ``state - K(theta)`` via the symplectic form.

    The centre bundle is symplectically orthogonal to both hyperbolic bundles (same
    argument as CMDA Eq. 68), so ``c_u = Omega(v_s, d)/Omega(v_s, v_u)`` and
    ``c_s = Omega(v_u, d)/Omega(v_u, v_s)``.
    """
    d = np.asarray(state) - circle.state(theta)
    vs = fourier_eval(bundles.v_s, theta)
    vu = fourier_eval(bundles.v_u, theta)
    om_su = float(vs @ _J4 @ vu)
    c_u = float(vs @ _J4 @ d) / om_su
    c_s = float(vu @ _J4 @ d) / (-om_su)
    return c_s, c_u


def symplectic_defect(phi: FloatArray) -> float:
    """``max |phi^T J phi - J|``."""
    return float(np.max(np.abs(phi.T @ _J4 @ phi - _J4)))
