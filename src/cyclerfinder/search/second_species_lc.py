"""Planar moon-centred Levi-Civita propagator with transition matrices, compiled (#899 step 2).

A fast planar companion of :mod:`cyclerfinder.core.cr3bp_ks` (#928) for the mass continuation
of :mod:`cyclerfinder.search.second_species_continuation`, where thousands of one-period
transition matrices are needed. Same model and frame as the project's CR3BP (primary of mass
1 - mu at (-mu, 0), secondary of mass mu at (1 - mu, 0), unit rotation), regularised at the
SECONDARY by z = x - (1 - mu) + i y = w^2 and dt = |w|^2 ds. With w' = dw/ds and the Jacobi
constant C carried as a state (it is constant), the equations are

    w'' = -2 i r w' + w (2 Omega' - C) / 4 + (r conj(w) / 2) G,       t' = r = |w|^2,

where Omega' = |z + 1 - mu|^2 / 2 + (1 - mu)/|z + 1| is the potential without the secondary's
term and G = dOmega'/dx + i dOmega'/dy. The secondary's singular term cancels exactly against
the energy relation |w'|^2 = r (2 Omega' - C)/4 + mu/2 (derivation in the #899 step 2 record);
every term is regular at w = 0. Physical velocity: zdot = 2 w' / conj(w).

Transition matrix: the 6 x 6 variational equations of (a, b, p, q, C, t) (w = a + i b,
w' = p + i q) with the analytic Jacobian, lifted by the Jacobian of (x, y, vx, vy) -> (a, b, p,
q, C, 0), projected back, and corrected to fixed physical time by subtracting zdot (dt/dz0),
as in :mod:`cyclerfinder.core.cr3bp_ks`. The vertical (z, vz) variational block of the planar
orbit is integrated alongside in the regular form zeta' = r vz, vz' = -[(1 - mu) r / r1^3 +
mu / r^2] zeta (no fixed-time correction is needed: zdot = 0 on the planar orbit).

Integrator: Dormand & Prince DOP853 with scipy's coefficients and scipy's error norm and
step-size rule (``scipy.integrate._ivp.rk``), compiled with numba; the end point at the exact
physical time is found by Newton steps in s from the last accepted step. The periselene and
perigee are located inside each accepted step by an Illinois root of the distance rate on
state-only steps from the step's start.

Validated against :func:`cyclerfinder.core.cr3bp_ks.propagate_ks` (state, 6 x 6 transition
matrix, periselene) in ``tests/search/test_second_species_continuation.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import njit
from numpy.typing import NDArray
from scipy.integrate._ivp import dop853_coefficients as _dop

FloatArray = NDArray[np.float64]

_NS = 12
_A = np.ascontiguousarray(_dop.A[:_NS, :_NS], dtype=np.float64)
_B = np.ascontiguousarray(_dop.B, dtype=np.float64)
_E3 = np.ascontiguousarray(_dop.E3, dtype=np.float64)
_E5 = np.ascontiguousarray(_dop.E5, dtype=np.float64)

NSTATE = 6  # a, b, p, q, C, t
NPLANAR = NSTATE + 36
NFULL = NPLANAR + 4  # plus the vertical 2 x 2 block


@njit(cache=True)  # type: ignore[untyped-decorator]
def _potential_terms(a: float, b: float, mu: float) -> tuple[float, ...]:
    xr = a * a - b * b
    y = 2.0 * a * b
    d = xr + 1.0
    om = 1.0 - mu
    r1sq = d * d + y * y
    r1 = math.sqrt(r1sq)
    r13 = r1sq * r1
    bx = xr + om  # barycentric x
    omega = 0.5 * (bx * bx + y * y) + om / r1
    gx = bx - om * d / r13
    gy = y - om * y / r13
    r15 = r13 * r1sq
    hxx = 1.0 - om * (1.0 / r13 - 3.0 * d * d / r15)
    hyy = 1.0 - om * (1.0 / r13 - 3.0 * y * y / r15)
    hxy = 3.0 * om * d * y / r15
    return xr, y, r1, omega, gx, gy, hxx, hxy, hyy


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs(y: FloatArray, mu: float, mode: int, out: FloatArray) -> None:
    """mode 0: state only (6); 1: state + planar STM (42); 2: plus vertical block (46)."""
    a = y[0]
    b = y[1]
    p = y[2]
    q = y[3]
    cj = y[4]
    r = a * a + b * b
    _xr, _yy, r1, omega, gx, gy, hxx, hxy, hyy = _potential_terms(a, b, mu)
    k = (2.0 * omega - cj) * 0.25
    s1 = a * gx + b * gy
    s2 = a * gy - b * gx
    out[0] = p
    out[1] = q
    out[2] = 2.0 * r * q + a * k + 0.5 * r * s1
    out[3] = -2.0 * r * p + b * k + 0.5 * r * s2
    out[4] = 0.0
    out[5] = r
    if mode == 0:
        return
    # Jacobian (6 x 6) of the state equations
    xa = 2.0 * a
    xb = -2.0 * b
    ya = 2.0 * b
    yb = 2.0 * a
    ra = 2.0 * a
    rb = 2.0 * b
    gxa = hxx * xa + hxy * ya
    gxb = hxx * xb + hxy * yb
    gya = hxy * xa + hyy * ya
    gyb = hxy * xb + hyy * yb
    ka = 0.5 * (gx * xa + gy * ya)
    kb = 0.5 * (gx * xb + gy * yb)
    s1a = gx + a * gxa + b * gya
    s1b = a * gxb + gy + b * gyb
    s2a = gy + a * gya - b * gxa
    s2b = a * gyb - gx - b * gxb
    jac = np.zeros((6, 6))
    jac[0, 2] = 1.0
    jac[1, 3] = 1.0
    jac[2, 0] = 2.0 * ra * q + k + a * ka + 0.5 * ra * s1 + 0.5 * r * s1a
    jac[2, 1] = 2.0 * rb * q + a * kb + 0.5 * rb * s1 + 0.5 * r * s1b
    jac[2, 3] = 2.0 * r
    jac[2, 4] = -0.25 * a
    jac[3, 0] = -2.0 * ra * p + b * ka + 0.5 * ra * s2 + 0.5 * r * s2a
    jac[3, 1] = -2.0 * rb * p + k + b * kb + 0.5 * rb * s2 + 0.5 * r * s2b
    jac[3, 2] = -2.0 * r
    jac[3, 4] = -0.25 * b
    jac[5, 0] = ra
    jac[5, 1] = rb
    for i in range(6):
        for j in range(6):
            acc = 0.0
            for m in range(6):
                acc += jac[i, m] * y[NSTATE + 6 * m + j]
            out[NSTATE + 6 * i + j] = acc
    if mode == 1:
        return
    # vertical block: zeta' = r vz, vz' = -[(1 - mu) r / r1^3 + mu / r^2] zeta
    coef = (1.0 - mu) * r / (r1 * r1 * r1) + mu / (r * r)
    base = NPLANAR
    # row-major 2x2 [[zz, zv], [vz, vv]]
    out[base + 0] = r * y[base + 2]
    out[base + 1] = r * y[base + 3]
    out[base + 2] = -coef * y[base + 0]
    out[base + 3] = -coef * y[base + 1]


@njit(cache=True)  # type: ignore[untyped-decorator]
def _step(
    y: FloatArray, f0: FloatArray, h: float, mu: float, mode: int, kst: FloatArray
) -> tuple[FloatArray, FloatArray]:
    n = y.shape[0]
    for i in range(n):
        kst[0, i] = f0[i]
    tmp = np.empty(n)
    for s in range(1, _NS):
        for i in range(n):
            acc = 0.0
            for j in range(s):
                acc += _A[s, j] * kst[j, i]
            tmp[i] = y[i] + h * acc
        _rhs(tmp, mu, mode, kst[s])
    y_new = np.empty(n)
    for i in range(n):
        acc = 0.0
        for j in range(_NS):
            acc += _B[j] * kst[j, i]
        y_new[i] = y[i] + h * acc
    f_new = np.empty(n)
    _rhs(y_new, mu, mode, f_new)
    for i in range(n):
        kst[_NS, i] = f_new[i]
    return y_new, f_new


@njit(cache=True)  # type: ignore[untyped-decorator]
def _err_norm(
    kst: FloatArray, h: float, y: FloatArray, y_new: FloatArray, rtol: float, atol: float
) -> float:
    n = y.shape[0]
    e5 = 0.0
    e3 = 0.0
    for i in range(n):
        sc = atol + max(abs(y[i]), abs(y_new[i])) * rtol
        a5 = 0.0
        a3 = 0.0
        for j in range(_NS + 1):
            a5 += _E5[j] * kst[j, i]
            a3 += _E3[j] * kst[j, i]
        a5 /= sc
        a3 /= sc
        e5 += a5 * a5
        e3 += a3 * a3
    if e5 == 0.0 and e3 == 0.0:
        return 0.0
    den = e5 + 0.01 * e3
    return abs(h) * e5 / math.sqrt(den * n)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _earth_d2_and_rate(y: FloatArray, mu: float) -> tuple[float, float]:
    """Squared distance from the primary and its s-derivative."""
    a = y[0]
    b = y[1]
    p = y[2]
    q = y[3]
    xr = a * a - b * b
    yy = 2.0 * a * b
    dx = xr + 1.0
    # dz/ds = 2 w w'
    dxs = 2.0 * (a * p - b * q)
    dys = 2.0 * (a * q + b * p)
    return dx * dx + yy * yy, 2.0 * (dx * dxs + yy * dys)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _dist_and_rate(y: FloatArray, mu: float, kind: int) -> tuple[float, float]:
    if kind == 0:
        return y[0] * y[0] + y[1] * y[1], 2.0 * (y[0] * y[2] + y[1] * y[3])
    d2, rate = _earth_d2_and_rate(y, mu)
    return float(d2), float(rate)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _refine_min(
    y: FloatArray, h: float, mu: float, kind: int, d0: float, d1: float, fallback: float
) -> float:
    """Minimum of the distance (kind 0: |z|, secondary; kind 1: squared, primary) inside a step
    of size h from y whose rate goes from d0 < 0 to d1 >= 0: Illinois root of the rate on
    state-only steps."""
    ys = y[:NSTATE].copy()
    fs = np.empty(NSTATE)
    _rhs(ys, mu, 0, fs)
    kst = np.empty((_NS + 1, NSTATE))
    lo, hi, flo, fhi = 0.0, h, d0, d1
    best = fallback
    side = 0
    for _ in range(60):
        if fhi == flo:
            break
        m = hi - fhi * (hi - lo) / (fhi - flo)
        if not (min(lo, hi) < m < max(lo, hi)):
            m = 0.5 * (lo + hi)
        ym, _ = _step(ys, fs, m, mu, 0, kst)
        val, rate = _dist_and_rate(ym, mu, kind)
        if val < best:
            best = val
        if rate == 0.0 or abs(hi - lo) < 1e-15 * max(1.0, abs(h)):
            break
        if (rate < 0.0) == (flo < 0.0):
            lo, flo = m, rate
            if side == -1:
                fhi *= 0.5
            side = -1
        else:
            hi, fhi = m, rate
            if side == 1:
                flo *= 0.5
            side = 1
    return best


@njit(cache=True)  # type: ignore[untyped-decorator]
def _integrate(
    y0: FloatArray, mu: float, t_final: float, mode: int, rtol: float, atol: float, max_steps: int
) -> tuple[FloatArray, float, float, int, int]:
    """Integrate to physical time t_final (> y0[5]). Returns (y, r_min, d2_earth_min, nfev,
    status); status 0 ok, 1 too many steps, 2 step underflow."""
    n = y0.shape[0]
    kst = np.empty((_NS + 1, n))
    y = y0.copy()
    f = np.empty(n)
    _rhs(y, mu, mode, f)
    nfev = 1
    r0 = y[0] * y[0] + y[1] * y[1]
    r_min = r0
    d2_min, _ = _earth_d2_and_rate(y, mu)
    # initial step: a small fraction of the remaining time in s units
    h = min(1e-3, 0.01 * (t_final - y[5]) / max(r0, 1e-300))
    status = 1
    for _ in range(max_steps):
        while True:
            y_new, f_new = _step(y, f, h, mu, mode, kst)
            nfev += _NS
            err = _err_norm(kst, h, y, y_new, rtol, atol)
            if err < 1.0:
                fac = 10.0 if err == 0.0 else min(10.0, 0.9 * err ** (-1.0 / 8.0))
                break
            h *= max(0.2, 0.9 * err ** (-1.0 / 8.0))
            if abs(h) < 1e-300:
                return y, r_min, d2_min, nfev, 2
        # extrema inside the step
        ra = y[0] * y[0] + y[1] * y[1]
        rb = y_new[0] * y_new[0] + y_new[1] * y_new[1]
        da = 2.0 * (y[0] * y[2] + y[1] * y[3])
        db = 2.0 * (y_new[0] * y_new[2] + y_new[1] * y_new[3])
        crossed = y_new[5] >= t_final
        if not crossed:
            r_min = min(r_min, ra, rb)
            if da < 0.0 <= db:
                r_min = min(r_min, _refine_min(y, h, mu, 0, da, db, r_min))
            ea, eda = _earth_d2_and_rate(y, mu)
            eb, edb = _earth_d2_and_rate(y_new, mu)
            d2_min = min(d2_min, ea, eb)
            if eda < 0.0 <= edb:
                d2_min = min(d2_min, _refine_min(y, h, mu, 1, eda, edb, d2_min))
            y = y_new
            f = f_new
            h *= fac
            continue
        # end point: Newton in s from the last accepted point
        rr = ra
        hs = (t_final - y[5]) / rr
        ye = y
        for _k in range(10):
            ye, _fe = _step(y, f, hs, mu, mode, kst)
            nfev += _NS
            dt = t_final - ye[5]
            if abs(dt) <= 4e-16 * max(1.0, abs(t_final)):
                break
            hs += dt / (ye[0] * ye[0] + ye[1] * ye[1])
        rb = ye[0] * ye[0] + ye[1] * ye[1]
        db = 2.0 * (ye[0] * ye[2] + ye[1] * ye[3])
        r_min = min(r_min, ra, rb)
        if da < 0.0 <= db:
            r_min = min(r_min, _refine_min(y, hs, mu, 0, da, db, r_min))
        ea, eda = _earth_d2_and_rate(y, mu)
        eb, edb = _earth_d2_and_rate(ye, mu)
        d2_min = min(d2_min, ea, eb)
        if eda < 0.0 <= edb:
            d2_min = min(d2_min, _refine_min(y, hs, mu, 1, eda, edb, d2_min))
        return ye, r_min, d2_min, nfev, 0
    return y, r_min, d2_min, nfev, status


# ---------------------------------------------------------------------------------------------
# Python side: lift, projection, fixed-time correction


def _cmul(c: complex) -> FloatArray:
    return np.array([[c.real, -c.imag], [c.imag, c.real]])


def _cconj_mul(c: complex) -> FloatArray:
    """Real 2 x 2 matrix of dz -> c conj(dz)."""
    return np.array([[c.real, c.imag], [c.imag, -c.real]])


def jacobi_planar(state4: FloatArray, mu: float) -> float:
    x, y, vx, vy = (float(v) for v in state4)
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    return x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - vx * vx - vy * vy


def _jacobi_grad(state4: FloatArray, mu: float) -> FloatArray:
    x, y, vx, vy = (float(v) for v in state4)
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    ox = x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3
    oy = y - (1.0 - mu) * y / r1**3 - mu * y / r2**3
    return np.array([2.0 * ox, 2.0 * oy, -2.0 * vx, -2.0 * vy])


def lift(state4: FloatArray, mu: float) -> tuple[FloatArray, FloatArray]:
    """(a, b, p, q, C, t=0) of a planar state and its 6 x 4 Jacobian."""
    z = complex(state4[0] - (1.0 - mu), state4[1])
    zd = complex(state4[2], state4[3])
    w = complex(np.sqrt(z))
    if w == 0:
        raise ZeroDivisionError("lift: state at the secondary")
    wp = zd * w.conjugate() / 2.0
    cj = jacobi_planar(state4, mu)
    y = np.array([w.real, w.imag, wp.real, wp.imag, cj, 0.0])
    jac = np.zeros((6, 4))
    dw_dz = _cmul(1.0 / (2.0 * w))
    jac[0:2, 0:2] = dw_dz
    # dw' = dzdot conj(w)/2 + zdot conj(dw)/2, conj(dw) = conj(dz)/(2 conj(w))
    jac[2:4, 2:4] = _cmul(w.conjugate() / 2.0)
    jac[2:4, 0:2] = _cconj_mul(zd / (4.0 * w.conjugate()))
    jac[4, :] = _jacobi_grad(state4, mu)
    return y, jac


def project(y: FloatArray, mu: float) -> tuple[FloatArray, FloatArray]:
    """Planar physical state of (a, b, p, q, ...) and its 4 x 6 Jacobian."""
    w = complex(y[0], y[1])
    wp = complex(y[2], y[3])
    z = w * w
    zd = 2.0 * wp / w.conjugate()
    st = np.array([z.real + (1.0 - mu), z.imag, zd.real, zd.imag])
    jac = np.zeros((4, 6))
    jac[0:2, 0:2] = _cmul(2.0 * w)
    jac[2:4, 2:4] = _cmul(2.0 / w.conjugate())
    # dzdot = 2 dw'/conj(w) - 2 w' conj(dw)/conj(w)^2
    jac[2:4, 0:2] = _cconj_mul(-2.0 * wp / w.conjugate() ** 2)
    return st, jac


def planar_accel(state4: FloatArray, mu: float) -> FloatArray:
    x, y, vx, vy = (float(v) for v in state4)
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    ax = x + 2.0 * vy - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3
    ay = y - 2.0 * vx - (1.0 - mu) * y / r1**3 - mu * y / r2**3
    return np.array([vx, vy, ax, ay])


@dataclass(frozen=True)
class LCArc:
    end: FloatArray  # planar (x, y, vx, vy) at t
    stm4: FloatArray | None  # planar 4 x 4 fixed-time transition matrix
    stm_z: FloatArray | None  # vertical 2 x 2 (z, vz) block
    r_min: float  # closest approach to the secondary
    earth_min: float  # closest approach to the primary
    nfev: int
    jacobi_drift: float  # Jacobi constant at the end minus the carried constant


class LCPropagationError(RuntimeError):
    pass


def propagate_lc(
    mu: float,
    state4: FloatArray,
    t: float,
    *,
    with_stm: bool = True,
    rtol: float = 1e-12,
    atol: float = 1e-14,
    max_steps: int = 200_000,
) -> LCArc:
    """Propagate a planar state for physical time ``t > 0``."""
    if not t > 0.0:
        raise ValueError("propagate_lc: t must be positive")
    y6, j0 = lift(np.asarray(state4, dtype=np.float64), mu)
    mode = 2 if with_stm else 0
    y0 = np.concatenate([y6, np.eye(6).ravel(), np.eye(2).ravel()]) if with_stm else y6
    ye, r_min, d2_min, nfev, status = _integrate(y0, mu, t, mode, rtol, atol, max_steps)
    if status != 0:
        raise LCPropagationError(f"propagate_lc: integration failed (status {status})")
    end, pj = project(ye, mu)
    drift = jacobi_planar(end, mu) - float(ye[4])
    stm4 = None
    stm_z = None
    if with_stm:
        phi = ye[NSTATE:NPLANAR].reshape(6, 6)
        chain = phi @ j0  # d y_end / d z0 at fixed s
        stm4 = pj @ chain - np.outer(planar_accel(end, mu), chain[5])
        stm_z = ye[NPLANAR:NFULL].reshape(2, 2).copy()
    return LCArc(
        end=end,
        stm4=stm4,
        stm_z=stm_z,
        r_min=float(r_min),
        earth_min=math.sqrt(d2_min),
        nfev=int(nfev),
        jacobi_drift=float(drift),
    )
