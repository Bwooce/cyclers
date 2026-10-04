"""Rerun of the #884 Sun-forced Earth-Moon periodic-orbit search in the corrected Sun models (#905).

What this module does
---------------------
A three-body (CR3BP) periodic orbit of period ``T*`` whose period is commensurate with the
Sun's synodic period ``Tg = 2 pi / omega_S`` (``M T* = N Tg``, ``gcd(M, N) = 1``) is continued
into a Sun-forced model as a periodic orbit of period ``P = N Tg``: ``N`` Sun periods, so that the
Sun is back at its starting phase when the orbit closes (not the three-body period ``T*``; the
closure is checked explicitly by :func:`orbit_diagnostics`). Two forced models:

* ``bcr4bp``: the corrected bicircular model of :mod:`cyclerfinder.core.bcr4bp` (Sun regressing,
  ``theta_S = theta_S0 - omega_S t``). Homotopy ``eps`` multiplies the Sun's mass (direct and
  indirect term), as in Oshima (2022) eq. 2.
* ``qbcp``: the corrected coherent model of :mod:`cyclerfinder.core.qbcp` (Andreu's Fourier
  coefficients, canonical variables, alpha_6 on the whole potential). Homotopy
  ``H = H_RTBP + eps (H_QBCP - H_RTBP)``, which is Leiva & Briozzo (2008) eq. 5. ``H`` is affine in
  the eight alpha functions, so ``eps`` blends each alpha from its three-body value and scales the
  Sun's potential.

Both fields are linear in ``eps``: ``f_eps = f_0 + eps g``. They are re-implemented here in numba
(with the state transition matrix and the ``eps`` sensitivity) for speed; the tests gate them
against ``core.bcr4bp.bcr4bp_eom`` and ``core.qbcp.qbcp_eom``. Neither core module is modified.

Method (sources)
----------------
* Multiple shooting with the node states as unknowns over the whole forced period, nodes placed
  by the growth of the state transition matrix (Kallrath, Schloder & Bock 1993, CMDA 56:353:
  single shooting fails on chaotic systems; Leiva & Briozzo 2008 p239 report exactly that failure
  over five or more Sun periods).
* Screens before continuation (Rhouma & Chicone 2000, Theorem 3.1 and section 4 of the digest):
  commensurate period; ``dT*/dC != 0``; no nontrivial one-lap multiplier with ``lambda^M = 1``;
  Melnikov function not identically zero, with simple zeros.
* Melnikov function in the Sun phase: the first-order change of the Jacobi constant over ``P``,
  ``Mel(tau) = int_0^P grad C(y(s)) . g(tau + s, y(s)) ds``, where ``tau`` is the clock time at
  which the forced orbit passes the parent's reference point. Periodic in ``tau``, period
  ``Tg / M``.
  Evaluated by the periodic trapezoid rule on one lap of the parent (spectrally accurate for a
  periodic integrand).
* Continuation in ``eps`` from 0 to 1 by pseudo-arclength in (nodes, eps), with a predictor guard
  and a tangent-cosine guard; folds (sign change of the ``eps`` component of the tangent), branch
  points (sign change of the bordered determinant) and crossings of ``eps = 0`` are recorded, and
  a bounded excursion to negative ``eps`` is allowed (Oshima 2022 Fig. 4; Rosales et al. 2021).
* Family walk in the CR3BP with bifurcation detection: crossings of the stability indices
  through ``2 cos(2 pi k / m)`` (roots of unity, Sanaga & Howell 2025), and arrival at the plane
  of a spatial family (pitchfork) are recorded (#884 review section 14 item 5).
* Diagnostics: refined periselene and perigee with a surface exclusion, minimal period and
  planarity of the parent, Floquet magnitudes with the reciprocal-pair reliability check, a
  separate-integrator (Radau) closure, the Sun-phase closure check, and one orbit per symmetry
  class (time shift by whole Sun periods and the reversing symmetry).

Clocks. ``bcr4bp``: the Sun is at angle ``theta_S0 - omega_S t`` (on +x at t = 0 for the default
``theta_S0 = 0``). ``qbcp``: the alpha functions are evaluated at ``theta_S0 + omega_S t`` and the
Sun is on -x at t = 0 (Andreu's syzygy, the paper's +x after the rotation by pi).

Discipline: pure numerics; no catalogue writeback.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import field as dc_field
from fractions import Fraction
from typing import Any

import numba
import numpy as np
from numba import njit, prange
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.integrate._ivp import dop853_coefficients as _dc
from scipy.optimize import brentq, minimize_scalar

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.qbcp as qbcp

FloatArr = NDArray[np.float64]

#: Earth-Moon distance used to convert lengths to km (Leiva & Briozzo 2008 p228; catalogue unit).
EM_LENGTH_KM: float = 384_400.0
#: Mean lunar radius, km (surface exclusion; #884 review section 14 item 4).
MOON_RADIUS_KM: float = 1737.4
#: Earth equatorial radius, km.
EARTH_RADIUS_KM: float = 6378.1
#: Lunar Laplace sphere of influence quoted in the catalogue comments, km (cycler-class label).
LUNAR_SOI_KM: float = 66_182.9

KIND_BCR4BP = 0
KIND_QBCP = 1

# DOP853 tableau (Hairer's coefficients, as shipped by scipy).
_A = np.ascontiguousarray(_dc.A[: _dc.N_STAGES, : _dc.N_STAGES], dtype=np.float64)
_B = np.ascontiguousarray(_dc.B, dtype=np.float64)
_C = np.ascontiguousarray(_dc.C[: _dc.N_STAGES], dtype=np.float64)
_E3 = np.ascontiguousarray(_dc.E3, dtype=np.float64)
_E5 = np.ascontiguousarray(_dc.E5, dtype=np.float64)

# Andreu's alpha tables as held by core.qbcp (read only). Parity 1 = cosine, 0 = sine
# (Andreu 1998 Table 1.5; #892).
_QBCP_TABLES = (
    qbcp._COEFFS_ALPHA1,
    qbcp._COEFFS_ALPHA2,
    qbcp._COEFFS_ALPHA3,
    qbcp._COEFFS_ALPHA4,
    qbcp._COEFFS_ALPHA5,
    qbcp._COEFFS_ALPHA6,
    qbcp._COEFFS_ALPHA7,
    qbcp._COEFFS_ALPHA8,
)
_QBCP_PARITY = np.array([0, 1, 0, 1, 1, 0, 1, 1, 0], dtype=np.int64)  # index 1..8


def _qbcp_coeff_matrix() -> FloatArr:
    width = max(len(c) for c in _QBCP_TABLES)
    m = np.zeros((9, width), dtype=np.float64)
    for i, c in enumerate(_QBCP_TABLES, start=1):
        m[i, : len(c)] = c
    return m


# ---------------------------------------------------------------------------
# Model description.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SunModel:
    """A Sun-forced Earth-Moon model with a homotopy parameter ``eps`` (0 = CR3BP, 1 = full)."""

    kind: str  # "bcr4bp" or "qbcp"
    mu: float
    mu_sun: float
    a_sun: float
    omega_sun: float
    theta_sun0: float = 0.0

    @property
    def kind_id(self) -> int:
        return KIND_BCR4BP if self.kind == "bcr4bp" else KIND_QBCP

    @property
    def tg(self) -> float:
        """Forcing period: the Sun's synodic period ``2 pi / omega_S`` (TU)."""
        return 2.0 * math.pi / self.omega_sun

    @property
    def prm(self) -> FloatArr:
        return np.array(
            [self.mu, self.mu_sun, self.a_sun, self.omega_sun, self.theta_sun0], dtype=np.float64
        )

    @property
    def coeffs(self) -> FloatArr:
        return _COEFFS if self.kind == "qbcp" else _NO_COEFFS


_COEFFS = _qbcp_coeff_matrix()
_NO_COEFFS = np.zeros((9, 1), dtype=np.float64)


def bcr4bp_model(system: bcr4bp.BCR4BPSystem | None = None) -> SunModel:
    """The corrected bicircular model (default: ``core.bcr4bp.andreu_default()``)."""
    s = system if system is not None else bcr4bp.andreu_default()
    return SunModel("bcr4bp", s.mu, s.mu_sun, s.a_sun_nondim, s.omega_sun_nondim, s.theta_sun0)


def qbcp_model(system: qbcp.QBCPSystem | None = None) -> SunModel:
    """The corrected coherent model (default: ``core.qbcp.qbcp_default()``)."""
    s = system if system is not None else qbcp.qbcp_default()
    return SunModel("qbcp", s.mu, s.mu_sun, s.a_sun_nondim, s.omega_sun_nondim, s.theta_sun0)


# ---------------------------------------------------------------------------
# Fields (numba). State: 6 model coordinates; with mode 1 also the STM (36, row-major) and the
# eps sensitivity dX/deps (6), 48 in all.
# ---------------------------------------------------------------------------


@njit(cache=True)  # type: ignore[untyped-decorator]
def _alphas(theta: float, coeffs: FloatArr, parity: NDArray[np.int64], out: FloatArr) -> None:
    width = coeffs.shape[1]
    c1 = math.cos(theta)
    s1 = math.sin(theta)
    for i in range(1, 9):
        out[i] = coeffs[i, 0]
    ck = 1.0
    sk = 0.0
    for k in range(1, width):
        ck, sk = ck * c1 - sk * s1, sk * c1 + ck * s1
        for i in range(1, 9):
            if parity[i] == 1:
                out[i] += coeffs[i, k] * ck
            else:
                out[i] += coeffs[i, k] * sk


@njit(cache=True)  # type: ignore[untyped-decorator]
def _point_mass(
    m: float, dx: float, dy: float, dz: float, grad: FloatArr, hess: FloatArr, add: bool
) -> None:
    """Accumulate ``m d / r^3`` (gradient of -m/r) and its Jacobian ``m (I/r^3 - 3 d d^T/r^5)``."""
    r2 = dx * dx + dy * dy + dz * dz
    r = math.sqrt(r2)
    i3 = 1.0 / (r2 * r)
    i5 = i3 / r2
    grad[0] += m * dx * i3
    grad[1] += m * dy * i3
    grad[2] += m * dz * i3
    if add:
        d0 = dx
        d1 = dy
        d2 = dz
        hess[0, 0] += m * (i3 - 3.0 * d0 * d0 * i5)
        hess[1, 1] += m * (i3 - 3.0 * d1 * d1 * i5)
        hess[2, 2] += m * (i3 - 3.0 * d2 * d2 * i5)
        hess[0, 1] += -3.0 * m * d0 * d1 * i5
        hess[0, 2] += -3.0 * m * d0 * d2 * i5
        hess[1, 2] += -3.0 * m * d1 * d2 * i5
        hess[1, 0] = hess[0, 1]
        hess[2, 0] = hess[0, 2]
        hess[2, 1] = hess[1, 2]


@njit(cache=True)  # type: ignore[untyped-decorator]
def _field(
    t: float,
    y: FloatArr,
    eps: float,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
    mode: int,
    out: FloatArr,
    wk: FloatArr,
) -> None:
    """``out = f_eps(t, y)`` (and the STM / sensitivity derivatives when ``mode == 1``).

    ``wk`` is caller-owned scratch of length at least 80 (no allocation per call)."""
    mu = prm[0]
    mus = prm[1]
    a_s = prm[2]
    w_s = prm[3]
    th0 = prm[4]
    x = y[0]
    yy = y[1]
    z = y[2]
    om = 1.0 - mu
    wk[:80] = 0.0
    gp = wk[0:3]  # gradient of the primaries' potential term
    hp = wk[3:12].reshape((3, 3))
    gs = wk[12:15]  # gradient of the Sun's direct term (unit Sun mass)
    hs = wk[15:24].reshape((3, 3))
    want = mode == 1
    _point_mass(om, x + mu, yy, z, gp, hp, want)
    _point_mass(mu, x - 1.0 + mu, yy, z, gp, hp, want)
    a = wk[24:60].reshape((6, 6))  # Jacobian of f_eps
    g = wk[60:66]  # d f / d eps
    if kind == 0:
        vx = y[3]
        vy = y[4]
        vz = y[5]
        th = th0 - w_s * t
        sx = a_s * math.cos(th)
        sy = a_s * math.sin(th)
        _point_mass(1.0, x - sx, yy - sy, z, gs, hs, want)
        a3 = a_s * a_s * a_s
        g[3] = -mus * gs[0] - mus * sx / a3
        g[4] = -mus * gs[1] - mus * sy / a3
        g[5] = -mus * gs[2]
        out[0] = vx
        out[1] = vy
        out[2] = vz
        out[3] = 2.0 * vy + x - gp[0] + eps * g[3]
        out[4] = -2.0 * vx + yy - gp[1] + eps * g[4]
        out[5] = -gp[2] + eps * g[5]
        if want:
            a[0, 3] = 1.0
            a[1, 4] = 1.0
            a[2, 5] = 1.0
            for i in range(3):
                for j in range(3):
                    a[3 + i, j] = -hp[i, j] - eps * mus * hs[i, j]
            a[3, 0] += 1.0
            a[4, 1] += 1.0
            a[3, 4] = 2.0
            a[4, 3] = -2.0
    else:
        px = y[3]
        py = y[4]
        pz = y[5]
        al = wk[66:75]
        _alphas(th0 + w_s * t, coeffs, parity, al)
        b1 = 1.0 + eps * (al[1] - 1.0)
        b2 = eps * al[2]
        b3 = 1.0 + eps * (al[3] - 1.0)
        b4 = eps * al[4]
        b5 = eps * al[5]
        b6p = 1.0 + eps * (al[6] - 1.0)
        b6s = eps * al[6] * mus
        _point_mass(1.0, x - al[7], yy - al[8], z, gs, hs, want)
        out[0] = b1 * px + b2 * x + b3 * yy
        out[1] = b1 * py + b2 * yy - b3 * x
        out[2] = b1 * pz + b2 * z
        out[3] = -b2 * px + b3 * py - b4 - b6p * gp[0] - b6s * gs[0]
        out[4] = -b2 * py - b3 * px - b5 - b6p * gp[1] - b6s * gs[1]
        out[5] = -b2 * pz - b6p * gp[2] - b6s * gs[2]
        d1 = al[1] - 1.0
        d3 = al[3] - 1.0
        d6 = al[6] - 1.0
        g[0] = d1 * px + al[2] * x + d3 * yy
        g[1] = d1 * py + al[2] * yy - d3 * x
        g[2] = d1 * pz + al[2] * z
        g[3] = -al[2] * px + d3 * py - al[4] - d6 * gp[0] - al[6] * mus * gs[0]
        g[4] = -al[2] * py - d3 * px - al[5] - d6 * gp[1] - al[6] * mus * gs[1]
        g[5] = -al[2] * pz - d6 * gp[2] - al[6] * mus * gs[2]
        if want:
            for i in range(3):
                a[i, i] = b2
                a[i, 3 + i] = b1
                a[3 + i, 3 + i] = -b2
                for j in range(3):
                    a[3 + i, j] = -b6p * hp[i, j] - b6s * hs[i, j]
            a[0, 1] = b3
            a[1, 0] = -b3
            a[3, 4] = b3
            a[4, 3] = -b3
    if want:
        for i in range(6):
            for j in range(6):
                acc = 0.0
                for k in range(6):
                    acc += a[i, k] * y[6 + 6 * k + j]
                out[6 + 6 * i + j] = acc
            acc = 0.0
            for k in range(6):
                acc += a[i, k] * y[42 + k]
            out[42 + i] = acc + g[i]


@njit(cache=True)  # type: ignore[untyped-decorator]
def _gdot(
    t: float,
    y: FloatArr,
    w: FloatArr,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
    al: FloatArr,
) -> float:
    """``w . g(t, y)`` with ``g = d f / d eps`` (fields affine in eps). ``al`` is scratch."""
    mu = prm[0]
    mus = prm[1]
    x = y[0]
    yy = y[1]
    z = y[2]
    if kind == 0:
        a_s = prm[2]
        th = prm[4] - prm[3] * t
        sx = a_s * math.cos(th)
        sy = a_s * math.sin(th)
        dx = x - sx
        dy = yy - sy
        d2 = dx * dx + dy * dy + z * z
        i3 = 1.0 / (d2 * math.sqrt(d2))
        a3 = a_s * a_s * a_s
        return float(
            mus * (w[3] * (-dx * i3 - sx / a3) + w[4] * (-dy * i3 - sy / a3) + w[5] * (-z * i3))
        )
    _alphas(prm[4] + prm[3] * t, coeffs, parity, al)
    px = y[3]
    py = y[4]
    pz = y[5]
    om = 1.0 - mu
    r1sq = (x + mu) ** 2 + yy * yy + z * z
    r2sq = (x - 1.0 + mu) ** 2 + yy * yy + z * z
    i13 = om / (r1sq * math.sqrt(r1sq))
    i23 = mu / (r2sq * math.sqrt(r2sq))
    gpx = i13 * (x + mu) + i23 * (x - 1.0 + mu)
    gpy = (i13 + i23) * yy
    gpz = (i13 + i23) * z
    dx = x - al[7]
    dy = yy - al[8]
    d2 = dx * dx + dy * dy + z * z
    i3 = mus / (d2 * math.sqrt(d2))
    d1 = al[1] - 1.0
    d3 = al[3] - 1.0
    d6 = al[6] - 1.0
    g0 = d1 * px + al[2] * x + d3 * yy
    g1 = d1 * py + al[2] * yy - d3 * x
    g2 = d1 * pz + al[2] * z
    g3 = -al[2] * px + d3 * py - al[4] - d6 * gpx - al[6] * i3 * dx
    g4 = -al[2] * py - d3 * px - al[5] - d6 * gpy - al[6] * i3 * dy
    g5 = -al[2] * pz - d6 * gpz - al[6] * i3 * z
    return float(w[0] * g0 + w[1] * g1 + w[2] * g2 + w[3] * g3 + w[4] * g4 + w[5] * g5)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _forcing(
    t: float,
    y: FloatArr,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
) -> FloatArr:
    """``g = d f / d eps`` at ``(t, y)``."""
    out = np.empty(48)
    yy = np.zeros(48)
    yy[:6] = y[:6]
    _field(t, yy, 0.0, kind, prm, coeffs, parity, 1, out, np.zeros(80))
    return out[42:48].copy()


@njit(cache=True)  # type: ignore[untyped-decorator]
def _dop853(
    z0: FloatArr,
    t0: float,
    t1: float,
    eps: float,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
    mode: int,
    rtol: float,
    atol: float,
    max_steps: int,
    a: FloatArr,
    b: FloatArr,
    c: FloatArr,
    e3: FloatArr,
    e5: FloatArr,
) -> tuple[FloatArr, int]:
    """Single-trajectory DOP853 (scipy's step-size control, no dense output)."""
    n = z0.size
    y = z0.copy()
    t = t0
    if t1 == t0:
        return y, 0
    direction = 1.0 if t1 > t0 else -1.0
    k = np.empty((13, n))
    f = np.empty(n)
    wk = np.zeros(80)
    _field(t, y, eps, kind, prm, coeffs, parity, mode, f, wk)
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
            if h_abs < 1e-14:
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
                _field(t + c[s] * h, ytmp, eps, kind, prm, coeffs, parity, mode, k[s], wk)
            for i in range(n):
                acc = 0.0
                for j in range(12):
                    acc += b[j] * k[j, i]
                ynew[i] = y[i] + h * acc
            _field(t + h, ynew, eps, kind, prm, coeffs, parity, mode, k[12], wk)
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
    z0s: FloatArr,
    t0s: FloatArr,
    t1s: FloatArr,
    eps: float,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
    mode: int,
    rtol: float,
    atol: float,
    max_steps: int,
    a: FloatArr,
    b: FloatArr,
    c: FloatArr,
    e3: FloatArr,
    e5: FloatArr,
) -> tuple[FloatArr, NDArray[np.int64]]:
    m = z0s.shape[0]
    out = np.empty_like(z0s)
    status = np.empty(m, dtype=np.int64)
    for i in prange(m):
        yi, st = _dop853(
            z0s[i],
            t0s[i],
            t1s[i],
            eps,
            kind,
            prm,
            coeffs,
            parity,
            mode,
            rtol,
            atol,
            max_steps,
            a,
            b,
            c,
            e3,
            e5,
        )
        out[i] = yi
        status[i] = st
    return out, status


def _limit_blas_threads(n: int) -> None:
    """Best effort: cap the OpenBLAS threads bundled with numpy and scipy.

    The shooting systems are a few hundred unknowns; with numba's workers busy, OpenBLAS's own
    thread pool oversubscribes a shared machine and a 121 x 121 solve was measured at 0.34 s
    instead of under 1 ms."""
    import ctypes
    import glob
    import os

    import scipy

    dirs = [os.path.dirname(np.__file__) + ".libs", os.path.dirname(scipy.__file__) + ".libs"]
    names = (
        "scipy_openblas_set_num_threads64_",
        "scipy_openblas_set_num_threads",
        "openblas_set_num_threads64_",
        "openblas_set_num_threads",
    )
    for d in dirs:
        for lib in glob.glob(os.path.join(d, "*openblas*")):
            try:
                handle = ctypes.CDLL(lib)
            except OSError:
                continue
            for name in names:
                fn = getattr(handle, name, None)
                if fn is not None:
                    fn(ctypes.c_int(n))
                    break


def set_threads(n: int = 4, blas: int = 1) -> None:
    """Cap numba's worker threads and the BLAS threads (the machine is shared)."""
    numba.set_num_threads(max(1, min(n, numba.config.NUMBA_NUM_THREADS)))
    _limit_blas_threads(blas)


# ---------------------------------------------------------------------------
# Python access.
# ---------------------------------------------------------------------------

RTOL = 1e-12
ATOL = 1e-13


def vector_field(model: SunModel, eps: float, t: float, y: FloatArr) -> FloatArr:
    """The model's vector field ``f_eps(t, y)`` in model coordinates (6)."""
    out = np.empty(6)
    _field(
        float(t),
        np.asarray(y, dtype=np.float64)[:6].copy(),
        float(eps),
        model.kind_id,
        model.prm,
        model.coeffs,
        _QBCP_PARITY,
        0,
        out,
        np.zeros(80),
    )
    return out


def forcing(model: SunModel, t: float, y: FloatArr) -> FloatArr:
    """``g = d f / d eps`` (model coordinates)."""
    return _forcing(  # type: ignore[no-any-return]
        float(t),
        np.asarray(y, dtype=np.float64),
        model.kind_id,
        model.prm,
        model.coeffs,
        _QBCP_PARITY,
    )


def flow(
    model: SunModel,
    eps: float,
    states: FloatArr,
    t0s: FloatArr | float,
    t1s: FloatArr | float,
    *,
    variational: bool = False,
    rtol: float = RTOL,
    atol: float = ATOL,
    max_steps: int = 5_000_000,
) -> tuple[FloatArr, FloatArr | None, FloatArr | None]:
    """Flow ``(m, 6)`` states from clock times ``t0s`` to ``t1s`` (in parallel).

    Returns final states ``(m, 6)``, STMs ``(m, 6, 6)`` and ``dX/deps`` ``(m, 6)`` (the last two
    only when ``variational``).
    """
    s = np.atleast_2d(np.asarray(states, dtype=np.float64))
    m = s.shape[0]
    t0a = np.broadcast_to(np.asarray(t0s, dtype=np.float64), (m,)).copy()
    t1a = np.broadcast_to(np.asarray(t1s, dtype=np.float64), (m,)).copy()
    if variational:
        z0 = np.hstack([s, np.tile(np.eye(6).reshape(-1), (m, 1)), np.zeros((m, 6))])
    else:
        z0 = s.copy()
    out, status = _batch(
        np.ascontiguousarray(z0),
        t0a,
        t1a,
        float(eps),
        model.kind_id,
        model.prm,
        model.coeffs,
        _QBCP_PARITY,
        1 if variational else 0,
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
        raise RuntimeError(f"integration failed for {int(np.sum(status < 0))} of {m} arcs")
    if variational:
        return out[:, :6].copy(), out[:, 6:42].reshape(m, 6, 6).copy(), out[:, 42:48].copy()
    return out[:, :6].copy(), None, None


def dense(
    model: SunModel,
    eps: float,
    y0: FloatArr,
    t0: float,
    t1: float,
    *,
    method: str = "DOP853",
    rtol: float = RTOL,
    atol: float = ATOL,
) -> Any:
    """scipy ``solve_ivp`` solution with dense output (diagnostics and the Radau cross-check)."""
    prm = model.prm
    coeffs = model.coeffs
    kid = model.kind_id
    wk = np.zeros(80)

    def fun(t: float, y: FloatArr) -> FloatArr:
        out = np.empty(6)
        _field(t, y, eps, kid, prm, coeffs, _QBCP_PARITY, 0, out, wk)
        return out

    sol = solve_ivp(  # type: ignore[call-overload]
        fun,
        (t0, t1),
        np.asarray(y0, dtype=np.float64),
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"dense propagation failed: {sol.message}")
    return sol


# ---------------------------------------------------------------------------
# Coordinates, Jacobi constant.
# ---------------------------------------------------------------------------


def blended_alphas(model: SunModel, t: float, eps: float) -> FloatArr:
    """The QBCP alphas at clock ``t`` blended with their three-body values (index 1..8)."""
    al = np.zeros(9)
    _alphas(model.theta_sun0 + model.omega_sun * t, _COEFFS, _QBCP_PARITY, al)
    base = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    out = base + eps * (al - base)
    out[7], out[8] = al[7], al[8]  # the Sun's position is not blended (its mass is)
    return np.asarray(out, dtype=np.float64)


def pv_to_model(model: SunModel, pv: FloatArr, t: float, eps: float) -> FloatArr:
    """Position-velocity state to model coordinates (identity for the bicircular model)."""
    s = np.asarray(pv, dtype=np.float64).copy()
    if model.kind == "bcr4bp":
        return s
    a = blended_alphas(model, t, eps)
    x, y, z, vx, vy, vz = s
    return np.array(
        [
            x,
            y,
            z,
            (vx - a[2] * x - a[3] * y) / a[1],
            (vy - a[2] * y + a[3] * x) / a[1],
            (vz - a[2] * z) / a[1],
        ]
    )


def model_to_pv(model: SunModel, st: FloatArr, t: float, eps: float) -> FloatArr:
    """Model coordinates to position-velocity (inverse of :func:`pv_to_model`)."""
    s = np.asarray(st, dtype=np.float64).copy()
    if model.kind == "bcr4bp":
        return s
    a = blended_alphas(model, t, eps)
    x, y, z, px, py, pz = s
    return np.array(
        [
            x,
            y,
            z,
            a[1] * px + a[2] * x + a[3] * y,
            a[1] * py + a[2] * y - a[3] * x,
            a[1] * pz + a[2] * z,
        ]
    )


def jacobi_pv(pv: FloatArr, mu: float) -> float:
    """Jacobi constant ``C = x^2 + y^2 + 2(1-mu)/r1 + 2 mu/r2 - v^2`` (Earth at -mu)."""
    x, y, z, vx, vy, vz = (float(v) for v in pv)
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return x * x + y * y + 2 * (1 - mu) / r1 + 2 * mu / r2 - (vx * vx + vy * vy + vz * vz)


def grad_jacobi_pv(pv: FloatArr, mu: float) -> FloatArr:
    """Gradient of :func:`jacobi_pv` with respect to the position-velocity state."""
    x, y, z, vx, vy, vz = (float(v) for v in pv)
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    c1 = 2 * (1 - mu) / r1**3
    c2 = 2 * mu / r2**3
    return np.array(
        [
            2 * x - c1 * (x + mu) - c2 * (x - 1 + mu),
            2 * y - (c1 + c2) * y,
            -(c1 + c2) * z,
            -2 * vx,
            -2 * vy,
            -2 * vz,
        ]
    )


def grad_jacobi(model: SunModel, st: FloatArr) -> FloatArr:
    """Gradient of the three-body Jacobi constant in the model's eps = 0 coordinates."""
    mu = model.mu
    x, y, z, a, b, c = (float(v) for v in st)
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    c1 = 2 * (1 - mu) / r1**3
    c2 = 2 * mu / r2**3
    gx = -c1 * (x + mu) - c2 * (x - 1 + mu)
    gy = -(c1 + c2) * y
    gz = -(c1 + c2) * z
    if model.kind == "bcr4bp":  # (x, v): C = x^2 + y^2 + 2U - v^2
        return np.array([2 * x + gx, 2 * y + gy, gz, -2 * a, -2 * b, -2 * c])
    # canonical (x, p): C = -2 H, H = |p|^2/2 + y px - x py - U
    return np.array([2 * b + gx, -2 * a + gy, gz, -2 * (a + y), -2 * (b - x), -2 * c])


def sun_angle(model: SunModel, t: float) -> float:
    """Direction of the Sun in the rotating frame at clock ``t`` (rad, in [0, 2 pi))."""
    if model.kind == "bcr4bp":
        return (model.theta_sun0 - model.omega_sun * t) % (2 * math.pi)
    a = blended_alphas(model, t, 1.0)
    return math.atan2(a[8], a[7]) % (2 * math.pi)


# ---------------------------------------------------------------------------
# Multiple shooting.
# ---------------------------------------------------------------------------


@dataclass
class Shooting:
    """Cyclic multiple shooting over ``[t0, t0 + period]`` with node times ``t0 + s_k``.

    ``offsets`` are the node offsets ``s_k`` (``s_0 = 0``, increasing, below ``period``). The
    unknowns are the node states in model coordinates; the residual is
    ``phi(X_k; t_k -> t_{k+1}) - X_{k+1}`` cyclically. For a forced model ``period`` is a whole
    number of Sun periods, so the last segment ends with the Sun at its starting phase.
    """

    model: SunModel
    period: float
    t0: float
    offsets: FloatArr

    @property
    def n_seg(self) -> int:
        return len(self.offsets)

    def times(self) -> FloatArr:
        return np.concatenate([self.t0 + self.offsets, [self.t0 + self.period]])

    def evaluate(self, xs: FloatArr, eps: float) -> tuple[FloatArr, FloatArr, FloatArr, FloatArr]:
        """Residual (6N), Jacobian (6N x 6N), d residual / d eps (6N), segment STMs (N, 6, 6)."""
        n = self.n_seg
        tt = self.times()
        xf, stms, sens = flow(self.model, eps, xs, tt[:-1], tt[1:], variational=True)
        assert stms is not None and sens is not None
        res = np.empty(6 * n)
        jac = np.zeros((6 * n, 6 * n))
        for k in range(n):
            kn = (k + 1) % n
            res[6 * k : 6 * k + 6] = xf[k] - xs[kn]
            jac[6 * k : 6 * k + 6, 6 * k : 6 * k + 6] += stms[k]
            jac[6 * k : 6 * k + 6, 6 * kn : 6 * kn + 6] -= np.eye(6)
        return res, jac, sens.reshape(-1), stms


def monodromy(stms: FloatArr) -> FloatArr:
    m = np.eye(6)
    for s in stms:
        m = s @ m
    return m


def node_offsets(
    model: SunModel,
    x0: FloatArr,
    lap: float,
    laps: int,
    *,
    growth: float = 200.0,
    max_len: float = 1.5,
    min_per_lap: int = 4,
    n_probe: int = 400,
) -> FloatArr:
    """Node offsets over ``laps`` laps of a three-body orbit, by STM growth.

    One lap is cut into ``n_probe`` pieces; consecutive pieces are merged while the largest
    singular value of the merged STM stays below ``growth`` and the merged length below
    ``max_len``. The lap pattern is repeated ``laps`` times.
    """
    cuts = np.linspace(0.0, lap, n_probe + 1)
    starts = np.empty((n_probe, 6))
    starts[0] = x0
    for i in range(1, n_probe):
        starts[i] = flow(model, 0.0, starts[i - 1], cuts[i - 1], cuts[i])[0][0]
    _, stms, _ = flow(model, 0.0, starts, cuts[:-1], cuts[1:], variational=True)
    assert stms is not None
    marks = [0.0]
    acc = np.eye(6)
    t_start = 0.0
    for i in range(n_probe):
        trial = stms[i] @ acc
        too_big = float(np.linalg.norm(trial, 2)) > growth or cuts[i + 1] - t_start > max_len
        if too_big and cuts[i] > t_start:
            marks.append(float(cuts[i]))
            t_start = float(cuts[i])
            acc = stms[i].copy()
        else:
            acc = trial
    per_lap = np.array(marks)
    if len(per_lap) < min_per_lap:
        per_lap = np.linspace(0.0, lap, min_per_lap + 1)[:-1]
    return np.concatenate([per_lap + j * lap for j in range(laps)])


def newton(
    prob: Shooting,
    xs: FloatArr,
    eps: float,
    *,
    tol: float = 1e-11,
    max_iter: int = 15,
    max_step: float = 0.1,
) -> tuple[FloatArr, float, bool]:
    """Damped Newton (least squares) on the shooting residual at fixed ``eps``."""
    nrm = math.inf
    for _ in range(max_iter):
        res, jac, _, _ = prob.evaluate(xs, eps)
        nrm = float(np.max(np.abs(res)))
        if not math.isfinite(nrm):
            return xs, nrm, False
        if nrm < tol:
            return xs, nrm, True
        dx = np.linalg.lstsq(jac, -res, rcond=None)[0]
        step = float(np.max(np.abs(dx)))
        if step > max_step:
            dx *= max_step / step
        xs = xs + dx.reshape(xs.shape)
    res, _, _, _ = prob.evaluate(xs, eps)
    nrm = float(np.max(np.abs(res)))
    return xs, nrm, nrm < tol


# ---------------------------------------------------------------------------
# Floquet reporting.
# ---------------------------------------------------------------------------


def floquet_report(mono: FloatArr) -> dict[str, Any]:
    """Multipliers, largest magnitude, and the reciprocal-pair reliability check (review s7)."""
    eig = np.linalg.eigvals(mono)
    mags = np.abs(eig)
    order = np.argsort(-mags)
    eig = eig[order]
    mags = mags[order]
    recip = float(mags[0] * mags[-1])
    reliable = abs(recip - 1.0) < 1e-3
    return {
        "multipliers": [[float(e.real), float(e.imag)] for e in eig],
        "max_abs": float(mags[0]),
        "reciprocal_product": recip,
        "subdominant_reliable": bool(reliable),
        "stable": bool(mags[0] < 1.0 + 1e-6),
    }


# ---------------------------------------------------------------------------
# Three-body parent: correction, screens.
# ---------------------------------------------------------------------------


@dataclass
class Parent:
    """A three-body periodic orbit at a commensurate period.

    ``state`` is in position-velocity at the reference point (parent time 0); ``period`` is the
    one-lap period ``T*``; ``laps = M`` and ``n_sun = N`` with ``M T* = N Tg``.
    """

    label: str
    state: FloatArr
    period: float
    laps: int
    n_sun: int
    residual: float = math.nan
    info: dict[str, Any] = dc_field(default_factory=dict)

    @property
    def forced_period(self) -> float:
        return self.laps * self.period


def commensurability(period: float, tg: float, max_den: int = 6) -> tuple[int, int, float]:
    """``(N, M, defect)`` with ``period / tg ~ N / M``; defect is ``|M period - N tg| / tg``."""
    fr = Fraction(period / tg).limit_denominator(max_den)
    n, m = fr.numerator, fr.denominator
    return n, m, abs(m * period - n * tg) / tg


def correct_parent(
    model: SunModel,
    pv0: FloatArr,
    period: float,
    *,
    n_nodes: int | None = None,
    tol: float = 1e-12,
) -> tuple[FloatArr, float]:
    """Correct a three-body periodic orbit at fixed period by multiple shooting (eps = 0).

    The orbit need not be symmetric. Returns the corrected position-velocity state at the
    first node and the residual.
    """
    x0 = pv_to_model(model, pv0, 0.0, 0.0)
    if n_nodes is None:
        offs = node_offsets(model, x0, period, 1)
    else:
        offs = np.linspace(0.0, period, n_nodes + 1)[:-1]
    prob = Shooting(model, period, 0.0, offs)
    xs = np.empty((len(offs), 6))
    xs[0] = x0
    for k in range(1, len(offs)):
        xs[k] = flow(model, 0.0, xs[k - 1], offs[k - 1], offs[k])[0][0]
    xs, nrm, _ = newton(prob, xs, 0.0, tol=tol, max_iter=25, max_step=0.02)
    return model_to_pv(model, xs[0], 0.0, 0.0), nrm


def lap_monodromy(model: SunModel, parent: Parent) -> FloatArr:
    """One-lap monodromy of the parent (multiple segments, eps = 0, model coordinates)."""
    x0 = pv_to_model(model, parent.state, 0.0, 0.0)
    offs = node_offsets(model, x0, parent.period, 1)
    tt = np.concatenate([offs, [parent.period]])
    xs = np.empty((len(offs), 6))
    xs[0] = x0
    for k in range(1, len(offs)):
        xs[k] = flow(model, 0.0, xs[k - 1], tt[k - 1], tt[k])[0][0]
    _, stms, _ = flow(model, 0.0, xs, tt[:-1], tt[1:], variational=True)
    assert stms is not None
    return monodromy(stms)


def stability_indices(mono: FloatArr, planar: bool) -> dict[str, float | None]:
    """Stability indices ``s = lambda + 1/lambda`` of the nontrivial pairs of a one-lap monodromy.

    Planar orbit (state ordered x, y, z, vx, vy, vz): in-plane ``tr(M_4) - 2`` and vertical
    ``tr(M_2)`` from the decoupled blocks (traces, robust at large multipliers). Spatial orbit:
    the two roots of ``s^2 - A s + (A^2 - B) / 2`` with ``A = tr M - 2`` and ``B = tr M^2 + 2``;
    ``None`` for a complex pair (Krein quadruplet).
    """
    if planar:
        ip = [0, 1, 3, 4]
        vp = [2, 5]
        return {
            "in_plane": float(np.trace(mono[np.ix_(ip, ip)]) - 2.0),
            "vertical": float(np.trace(mono[np.ix_(vp, vp)])),
        }
    a = float(np.trace(mono)) - 2.0
    b = float(np.trace(mono @ mono)) + 2.0
    disc = a * a - 2.0 * (a * a - b)
    if disc < 0:
        return {"s1": None, "s2": None}
    r = math.sqrt(disc)
    s1, s2 = 0.5 * (a + r), 0.5 * (a - r)
    return {"s1": s1, "s2": s2}


def root_of_unity_screen(mono: FloatArr, laps: int) -> dict[str, Any]:
    """Rhouma & Chicone 2(b): no nontrivial one-lap multiplier with ``lambda^M = 1``.

    The trivial pair is taken as the two multipliers closest to 1. Returns the smallest
    ``|lambda^M - 1|`` over the other four.
    """
    eig = np.linalg.eigvals(mono)
    order = np.argsort(np.abs(eig - 1.0))
    nontrivial = eig[order[2:]]
    d = [float(abs(lam**laps - 1.0)) for lam in nontrivial]
    return {
        "min_abs_lambda_M_minus_1": min(d),
        "nontrivial": [[float(e.real), float(e.imag)] for e in nontrivial],
    }


def period_slope(model: SunModel, parent: Parent, *, rel: float = 2e-5) -> tuple[float | None, str]:
    """``dT*/dC`` by fixed-period corrections at ``T*(1 +- rel)``; ``None`` if either fails.

    A failure near a period extremum (a fold in T) is itself the degeneracy being screened for.
    Prefer the family-walk tangent value when one is available.
    """
    cs = []
    for sgn in (1.0, -1.0):
        tp = parent.period * (1.0 + sgn * rel)
        pv, res = correct_parent(model, parent.state, tp)
        if not res < 1e-10:
            return None, f"fixed-period correction failed at T*(1{sgn:+.0f} rel) (res {res:.1e})"
        cs.append(jacobi_pv(pv, model.mu))
    dc = cs[0] - cs[1]
    if dc == 0.0:
        return math.inf, "C does not change with T"
    return 2.0 * rel * parent.period / dc, "finite difference"


def minimal_period_and_planarity(
    model: SunModel, parent: Parent, *, k_max: int = 6, tol: float = 1e-6
) -> dict[str, Any]:
    """Minimal period (does the orbit return at ``T*/k``?) and planarity of the parent."""
    x0 = pv_to_model(model, parent.state, 0.0, 0.0)
    sol = dense(model, 0.0, x0, 0.0, parent.period)
    returns = {}
    minimal_k = 1
    for k in range(2, k_max + 1):
        d = float(np.max(np.abs(sol.sol(parent.period / k) - x0)))
        returns[str(k)] = d
        if d < tol:
            minimal_k = k
    ts = np.linspace(0.0, parent.period, 4001)
    st = sol.sol(ts)
    zmax = float(np.max(np.abs(st[2])))
    vzmax = float(np.max(np.abs(st[5])))
    tt = np.array([0.0, parent.period])
    peri = min(_close_approaches([sol], tt, 1.0 - model.mu)) * EM_LENGTH_KM
    peri_e = min(_close_approaches([sol], tt, -model.mu)) * EM_LENGTH_KM
    return {
        "minimal_period_divisor": minimal_k,
        "return_at_T_over_k": returns,
        "planar": bool(zmax < 1e-10 and vzmax < 1e-10),
        "max_abs_z": zmax,
        "parent_periselene_km": peri,
        "parent_perigee_km": peri_e,
        "parent_below_moon_surface": bool(peri < MOON_RADIUS_KM),
        "parent_below_earth_surface": bool(peri_e < EARTH_RADIUS_KM),
    }


def screens(
    model: SunModel, parent: Parent, *, dtdc: float | None = None, dtdc_source: str = ""
) -> dict[str, Any]:
    """All Rhouma-Chicone screens and the section-14 parent checks, as one record."""
    n, m, defect = commensurability(parent.period, model.tg)
    out: dict[str, Any] = {
        "N": n,
        "M": m,
        "commensurability_defect": defect,
        "commensurate": bool(defect < 1e-9 and n == parent.n_sun and m == parent.laps),
    }
    mono = lap_monodromy(model, parent)
    out["parent_floquet"] = floquet_report(mono)
    shape = minimal_period_and_planarity(model, parent)
    out.update(shape)
    out["stability_indices"] = stability_indices(mono, shape["planar"])
    out["root_of_unity"] = root_of_unity_screen(mono, parent.laps)
    if dtdc is None:
        dtdc, dtdc_source = period_slope(model, parent)
    out["dT_dC"] = dtdc
    out["dT_dC_source"] = dtdc_source
    out["twist_nondegenerate"] = bool(dtdc is not None and abs(dtdc) > 1e-3)
    out["no_root_of_unity"] = bool(out["root_of_unity"]["min_abs_lambda_M_minus_1"] > 1e-3)
    out["pass"] = bool(
        out["commensurate"]
        and out["twist_nondegenerate"]
        and out["no_root_of_unity"]
        and shape["minimal_period_divisor"] == 1
    )
    return out


# ---------------------------------------------------------------------------
# Melnikov function in the Sun phase.
# ---------------------------------------------------------------------------


@njit(cache=True)  # type: ignore[untyped-decorator]
def _mel_values(
    taus: FloatArr,
    ys: FloatArr,
    ss: FloatArr,
    gradc: FloatArr,
    lap: float,
    laps: int,
    weight: float,
    kind: int,
    prm: FloatArr,
    coeffs: FloatArr,
    parity: NDArray[np.int64],
) -> FloatArr:
    out = np.zeros(taus.size)
    al = np.zeros(9)
    for it in range(taus.size):
        acc = 0.0
        for lp in range(laps):
            for j in range(ss.size):
                acc += _gdot(
                    taus[it] + ss[j] + lp * lap, ys[j], gradc[j], kind, prm, coeffs, parity, al
                )
        out[it] = acc * weight
    return out


@dataclass
class Melnikov:
    """The parent sampled on one lap for the Melnikov quadrature."""

    model: SunModel
    parent: Parent
    ss: FloatArr
    ys: FloatArr
    gradc: FloatArr

    def __call__(self, taus: FloatArr | float) -> FloatArr:
        t = np.atleast_1d(np.asarray(taus, dtype=np.float64))
        p = self.parent
        weight = p.period / self.ss.size  # periodic trapezoid
        return _mel_values(  # type: ignore[no-any-return]
            t,
            self.ys,
            self.ss,
            self.gradc,
            p.period,
            p.laps,
            weight,
            self.model.kind_id,
            self.model.prm,
            self.model.coeffs,
            _QBCP_PARITY,
        )

    @property
    def span(self) -> float:
        """Period of ``Mel`` in ``tau``: ``Tg / M``."""
        return self.model.tg / self.parent.laps


def melnikov(model: SunModel, parent: Parent, *, samples_per_tu: int = 2000) -> Melnikov:
    x0 = pv_to_model(model, parent.state, 0.0, 0.0)
    sol = dense(model, 0.0, x0, 0.0, parent.period)
    k = max(64, int(samples_per_tu * parent.period))
    ss = np.asarray(np.arange(k) * (parent.period / k), dtype=np.float64)
    ys = np.ascontiguousarray(sol.sol(ss).T)
    gradc = np.ascontiguousarray(np.array([grad_jacobi(model, y) for y in ys]))
    return Melnikov(model, parent, ss, ys, gradc)


def melnikov_variational(model: SunModel, parent: Parent, tau: float) -> float:
    """Cross-check: ``grad C(x(P)) . dX/deps(P)`` from the sensitivity equation at eps = 0."""
    x0 = pv_to_model(model, parent.state, 0.0, 0.0)
    xf, _, sens = flow(model, 0.0, x0, tau, tau + parent.forced_period, variational=True)
    assert sens is not None
    return float(grad_jacobi(model, xf[0]) @ sens[0])


def melnikov_zeros(mel: Melnikov, *, n_grid: int = 240) -> dict[str, Any]:
    """Zeros of ``Mel`` over one period ``[0, Tg/M)``, with slopes and simplicity."""
    span = mel.span
    cell = span / n_grid
    th = np.asarray(-0.5 * cell + cell * np.arange(n_grid + 1), dtype=np.float64)
    vals = mel(th)
    amp = float(np.max(np.abs(vals)))
    zeros = []

    def f(t: float) -> float:
        return float(mel(t)[0])

    for i in range(n_grid):
        if vals[i] * vals[i + 1] < 0:
            z = float(brentq(f, th[i], th[i + 1], xtol=1e-13, rtol=1e-14))
            h = 1e-5 * span
            slope = (f(z + h) - f(z - h)) / (2 * h)
            simple = bool(abs(slope) * span / (2 * math.pi) > 1e-4 * amp)
            zeros.append(
                {
                    "tau": z % span,
                    "slope": slope,
                    "simple": simple,
                    "sun_angle": sun_angle(mel.model, z % span),
                }
            )
    return {
        "span": span,
        "amplitude": amp,
        "identically_zero": bool(amp < 1e-12),
        "zeros": zeros,
        "grid": [float(v) for v in vals],
    }


# ---------------------------------------------------------------------------
# Continuation in eps.
# ---------------------------------------------------------------------------


@dataclass
class Branch:
    """One continuation branch in ``eps`` at a fixed Sun phase (``tau``)."""

    tau: float
    eps: list[float] = dc_field(default_factory=list)
    max_abs_floquet: list[float] = dc_field(default_factory=list)
    nodes_at: dict[float, FloatArr] = dc_field(default_factory=dict)
    folds: list[float] = dc_field(default_factory=list)
    branch_points: list[float] = dc_field(default_factory=list)
    eps0_crossings: list[float] = dc_field(default_factory=list)
    stop_reason: str = ""
    final_nodes: FloatArr | None = None
    final_eps: float = math.nan
    steps: int = 0
    wall_s: float = 0.0

    def summary(self) -> dict[str, Any]:
        return {
            "tau": self.tau,
            "stop_reason": self.stop_reason,
            "final_eps": self.final_eps,
            "steps": self.steps,
            "eps_min": min(self.eps) if self.eps else math.nan,
            "eps_max": max(self.eps) if self.eps else math.nan,
            "folds": self.folds,
            "branch_points": self.branch_points,
            "eps0_crossings": self.eps0_crossings,
            "max_abs_floquet_path": [
                [e, m]
                for e, m in zip(
                    self.eps[:: max(1, len(self.eps) // 40)],
                    self.max_abs_floquet[:: max(1, len(self.eps) // 40)],
                    strict=False,
                )
            ],
            "wall_s": self.wall_s,
        }


def _tangent(jac: FloatArr, jeps: FloatArr, prev: FloatArr | None, sign_eps: float) -> FloatArr:
    full = np.column_stack([jac, jeps])
    if prev is not None:
        # Null vector of the full-rank 6N x (6N+1) system from the bordered solve.
        rhs = np.zeros(full.shape[0] + 1)
        rhs[-1] = 1.0
        try:
            t = np.linalg.solve(np.vstack([full, prev]), rhs)
            t = t / np.linalg.norm(t)
            return np.asarray(t, dtype=np.float64)
        except np.linalg.LinAlgError:
            pass
    t = np.linalg.svd(full)[2][-1]
    if prev is not None:
        if t @ prev < 0:
            t = -t
    elif t[-1] * sign_eps < 0:
        t = -t
    return np.asarray(t, dtype=np.float64)


def _bordered_sign(jac: FloatArr, jeps: FloatArr, tan: FloatArr) -> float:
    big = np.vstack([np.column_stack([jac, jeps]), tan])
    sign, _ = np.linalg.slogdet(big)
    return float(sign)


def continue_in_eps(
    prob: Shooting,
    xs0: FloatArr,
    *,
    eps_start: float = 1e-4,
    ds0: float = 0.01,
    ds_max: float = 0.05,
    ds_min: float = 1e-6,
    eps_min: float = -0.5,
    eps_target: float = 1.0,
    max_steps: int = 600,
    tol: float = 1e-10,
    cos_min: float = 0.9,
    wall_s: float = math.inf,
    log: Callable[[str], None] | None = None,
) -> Branch:
    """Pseudo-arclength in ``(nodes, eps)`` from the three-body orbit at a Melnikov zero.

    ``xs0`` are the parent's node states at the zero's phase (eps = 0). The first point is a
    fixed-eps Newton solve at ``eps_start`` from the first-order predictor ``x + eps xi`` with
    ``(Phi - I) xi = -y_eps`` (Rhouma & Chicone eq. 3.4, minimum-norm solution). Steps are
    accepted only if Newton converges within 8 iterations, the corrector lands within half a
    step of the predictor, and the new unit tangent has cosine above ``cos_min`` with the old.
    The branch is followed through folds and through ``eps = 0`` (recorded) down to ``eps_min``;
    it stops at ``eps_target`` (corrected exactly there), below ``eps_min``, on step collapse,
    or at ``max_steps`` / ``wall_s``.
    """
    t_start = time.monotonic()
    n6 = 6 * prob.n_seg
    br = Branch(tau=prob.t0)
    _, jac0, jeps0, _ = prob.evaluate(xs0, 0.0)
    y = np.linalg.lstsq(jac0, -jeps0, rcond=None)[0]
    ok = False
    nrm = math.inf
    xs = xs0
    e_try = eps_start
    for e_try in (eps_start, 0.1 * eps_start, 0.01 * eps_start):
        xs, nrm, ok = newton(prob, xs0 + e_try * y.reshape(xs0.shape), e_try, tol=0.1 * tol)
        if ok:
            break
    if not ok:
        br.stop_reason = f"start Newton failed down to eps={e_try:.1e} (res {nrm:.2e})"
        br.wall_s = time.monotonic() - t_start
        return br
    eps = e_try
    res, jac, jeps, stms = prob.evaluate(xs, eps)
    br.eps.append(eps)
    br.max_abs_floquet.append(float(np.max(np.abs(np.linalg.eigvals(monodromy(stms))))))
    z = np.concatenate([xs.reshape(n6), [eps]])
    tan = _tangent(jac, jeps, None, 1.0)
    bsign = _bordered_sign(jac, jeps, tan)
    ds = ds0
    step = 0
    for step in range(max_steps):
        if time.monotonic() - t_start > wall_s:
            br.stop_reason = f"wall-clock limit at eps={z[-1]:.6g}"
            break
        pred = z + ds * tan
        zz = pred.copy()
        conv = False
        for _ in range(8):
            res, jac, jeps, stms = prob.evaluate(zz[:n6].reshape(-1, 6), float(zz[-1]))
            nrm = float(np.max(np.abs(res)))
            arc_res = float(tan @ (zz - pred))
            if not math.isfinite(nrm):
                break
            if nrm < tol and abs(arc_res) < 1e-10:
                conv = True
                break
            big = np.vstack([np.column_stack([jac, jeps]), tan])
            try:
                dz = np.linalg.solve(big, -np.concatenate([res, [arc_res]]))
            except np.linalg.LinAlgError:
                break
            if not np.all(np.isfinite(dz)) or float(np.max(np.abs(dz))) > 0.5:
                break
            zz = zz + dz
        newtan = _tangent(jac, jeps, tan, 1.0) if conv else tan
        if conv and float(np.max(np.abs(zz - pred))) > 0.5 * ds:
            conv = False  # corrector left the predictor's neighbourhood
        if conv and float(newtan @ tan) < cos_min:
            conv = False  # tangent turned too fast: possible branch jump
        if not conv:
            ds *= 0.5
            if ds < ds_min:
                br.stop_reason = f"step collapse at eps={z[-1]:.6g}"
                break
            continue
        e_old, e_new = float(z[-1]), float(zz[-1])
        if newtan[-1] * tan[-1] < 0:
            br.folds.append(e_new)
        nsign = _bordered_sign(jac, jeps, newtan)
        if nsign != bsign:
            br.branch_points.append(0.5 * (e_old + e_new))
        bsign = nsign
        if e_old > 0.0 >= e_new or e_old < 0.0 <= e_new:
            br.eps0_crossings.append(e_new)
        z_prev = z
        z, tan = zz, newtan
        br.eps.append(e_new)
        br.max_abs_floquet.append(float(np.max(np.abs(np.linalg.eigvals(monodromy(stms))))))
        if log is not None and step % 10 == 0:
            log(
                f"    step {step}: eps={e_new:.6f} ds={ds:.3g}"
                f" max|lam|={br.max_abs_floquet[-1]:.3e} folds={len(br.folds)}"
                f" bp={len(br.branch_points)}"
            )
        if (e_old - eps_target) * (e_new - eps_target) <= 0.0 and e_old != eps_target:
            w = (eps_target - e_old) / (e_new - e_old) if e_new != e_old else 1.0
            guess = (z_prev[:n6] + w * (z[:n6] - z_prev[:n6])).reshape(-1, 6)
            xs_t, nrm_t, ok_t = newton(prob, guess, eps_target, tol=tol)
            if ok_t:
                br.final_nodes = xs_t
                br.final_eps = eps_target
                br.stop_reason = "reached_target"
            else:
                br.stop_reason = f"target correction failed (res {nrm_t:.2e})"
            break
        if e_new < eps_min:
            br.stop_reason = f"below eps_min at eps={e_new:.6g}"
            break
        if float(np.max(np.abs(z[:n6]))) > 20:
            br.stop_reason = f"left domain at eps={e_new:.6g}"
            break
        ds = min(ds * 1.5, ds_max)
    else:
        br.stop_reason = f"max_steps at eps={z[-1]:.6g}"
    if br.final_nodes is None:
        br.final_nodes = z[:n6].reshape(-1, 6).copy()
        br.final_eps = float(z[-1])
    br.steps = step + 1
    br.wall_s = time.monotonic() - t_start
    return br


# ---------------------------------------------------------------------------
# The forced problem for one parent at one Melnikov zero.
# ---------------------------------------------------------------------------


def forced_problem(
    model: SunModel, parent: Parent, tau: float, **node_kw: Any
) -> tuple[Shooting, FloatArr]:
    """Shooting problem over ``P = N Tg`` from clock ``tau`` and the parent's nodes there."""
    x0 = pv_to_model(model, parent.state, 0.0, 0.0)
    offs = node_offsets(model, x0, parent.period, parent.laps, **node_kw)
    lap_offs = offs[offs < parent.period]
    xs_lap = np.empty((len(lap_offs), 6))
    xs_lap[0] = x0
    for k in range(1, len(lap_offs)):
        xs_lap[k] = flow(model, 0.0, xs_lap[k - 1], lap_offs[k - 1], lap_offs[k])[0][0]
    xs = np.vstack([xs_lap] * parent.laps)
    prob = Shooting(model, parent.forced_period, tau, offs)
    return prob, xs


# ---------------------------------------------------------------------------
# Diagnostics of a converged forced orbit.
# ---------------------------------------------------------------------------


def _close_approaches(
    sols: list[Any], tt: FloatArr, centre: float, per_tu: int = 4000
) -> list[float]:
    """Every local minimum of the distance to ``(centre, 0, 0)`` over the closed orbit, refined
    by bounded minimisation on the dense output (length units). Review section 8: a sampled
    minimum missed a pass by 395 km."""
    ts_all: list[FloatArr] = []
    seg_all: list[NDArray[np.int64]] = []
    d_all: list[FloatArr] = []
    for k, sol in enumerate(sols):
        a, b = tt[k], tt[k + 1]
        n = max(200, int(per_tu * (b - a)))
        ts = np.linspace(a, b, n)[:-1]
        st = sol.sol(ts)
        ts_all.append(ts)
        seg_all.append(np.full(ts.size, k))
        d_all.append(np.sqrt((st[0] - centre) ** 2 + st[1] ** 2 + st[2] ** 2))
    ts_c = np.concatenate(ts_all)
    seg = np.concatenate(seg_all)
    d = np.concatenate(d_all)
    n = d.size
    mins = []
    for i in range(n):
        if d[i] <= d[i - 1] and d[i] <= d[(i + 1) % n]:
            k = int(seg[i])
            sol = sols[k]
            lo = ts_c[i - 1] if i > 0 and seg[i - 1] == k else tt[k]
            hi = ts_c[i + 1] if i + 1 < n and seg[i + 1] == k else tt[k + 1]

            def dist(t: float, s: Any = sol) -> float:
                v = s.sol(t)
                return float(math.sqrt((v[0] - centre) ** 2 + v[1] ** 2 + v[2] ** 2))

            r = minimize_scalar(dist, bounds=(lo, hi), method="bounded", options={"xatol": 1e-12})
            mins.append(min(float(r.fun), float(d[i])))
    return mins


def orbit_diagnostics(
    prob: Shooting, xs: FloatArr, eps: float, *, radau: bool = True, samples_per_tg: int = 128
) -> dict[str, Any]:
    """Closure (DOP853 and Radau), Sun-phase closure, Floquet, refined passes, samples."""
    model = prob.model
    mu = model.mu
    res, _, _, stms = prob.evaluate(xs, eps)
    tt = prob.times()
    out: dict[str, Any] = {"closure_dop853": float(np.max(np.abs(res)))}
    n_sun = prob.period / model.tg
    out["period_in_sun_periods"] = n_sun
    out["sun_angle_start"] = sun_angle(model, tt[0])
    out["sun_angle_end"] = sun_angle(model, tt[-1])
    d_ang = (out["sun_angle_end"] - out["sun_angle_start"] + math.pi) % (2 * math.pi) - math.pi
    out["sun_phase_returns"] = bool(abs(n_sun - round(n_sun)) < 1e-12 and abs(d_ang) < 1e-9)
    sols = []
    rad = 0.0
    for k in range(prob.n_seg):
        sols.append(dense(model, eps, xs[k], tt[k], tt[k + 1]))
        if radau:
            s = dense(model, eps, xs[k], tt[k], tt[k + 1], method="Radau", rtol=1e-12, atol=1e-13)
            rad = max(rad, float(np.max(np.abs(s.y[:, -1] - xs[(k + 1) % prob.n_seg]))))
    out["closure_radau"] = rad if radau else None
    out["floquet"] = floquet_report(monodromy(stms))
    moon = _close_approaches(sols, tt, 1.0 - mu)
    earth = _close_approaches(sols, tt, -mu)
    out["periselene_km"] = min(moon) * EM_LENGTH_KM
    out["perigee_km"] = min(earth) * EM_LENGTH_KM
    out["lunar_minima_km"] = sorted(m * EM_LENGTH_KM for m in moon)
    out["n_lunar_minima_in_soi"] = int(sum(m * EM_LENGTH_KM <= LUNAR_SOI_KM for m in moon))
    out["below_moon_surface"] = bool(out["periselene_km"] < MOON_RADIUS_KM)
    out["below_earth_surface"] = bool(out["perigee_km"] < EARTH_RADIUS_KM)
    out["cycler_class"] = bool(out["periselene_km"] <= LUNAR_SOI_KM)
    # Positions on the absolute clock grid j Tg / samples_per_tg (for the symmetry classes).
    n_samp = round(n_sun) * samples_per_tg
    grid = np.arange(n_samp) * (model.tg / samples_per_tg)
    pos = np.empty((n_samp, 3))
    for j, tj in enumerate(grid):
        tl = prob.t0 + ((tj - prob.t0) % prob.period)
        k = int(np.searchsorted(tt, tl, side="right") - 1)
        k = min(max(k, 0), prob.n_seg - 1)
        pos[j] = sols[k].sol(tl)[:3]
    out["_clock_samples"] = pos
    out["nodes_pv"] = [
        model_to_pv(model, x, t, eps).tolist() for x, t in zip(xs, tt[:-1], strict=False)
    ]
    return out


def same_orbit(pos_a: FloatArr, pos_b: FloatArr, per_tg: int, *, tol: float = 1e-6) -> str | None:
    """How orbit B relates to orbit A on the absolute clock grid, or None.

    ``"same"``: ``B(t) = A(t + k Tg)`` for a whole number ``k`` (the same orbit started some Sun
    periods later). The symmetries of both Sun models (Oshima 2022 eqs. 6-8; both clocks are
    symmetric about t = 0 and the Sun is in the plane):
    ``"s2"``: ``B(t) = (x, -y, z) A(-t + k Tg)``;
    ``"s1"``: ``B(t) = (x, -y, -z) A(-t + k Tg)``;
    ``"s3"``: ``B(t) = (x, y, -z) A(t + k Tg)``.
    """
    if pos_a.shape != pos_b.shape:
        return None
    n = pos_a.shape[0]
    n_sun = n // per_tg
    rev_idx = (-np.arange(n)) % n
    images = {
        "same": pos_a,
        "s2": (pos_a * np.array([1.0, -1.0, 1.0]))[rev_idx],
        "s1": (pos_a * np.array([1.0, -1.0, -1.0]))[rev_idx],
        "s3": pos_a * np.array([1.0, 1.0, -1.0]),
    }
    for k in range(n_sun):
        sh = k * per_tg
        for name, img in images.items():
            if float(np.max(np.abs(np.roll(img, -sh, axis=0) - pos_b))) < tol:
                return name
    return None


def symmetry_classes(samples: list[FloatArr], per_tg: int, *, tol: float = 1e-6) -> list[int]:
    """Class index per orbit: orbits related by :func:`same_orbit` share a class."""
    cls: list[int] = []
    reps: list[int] = []
    for i, s in enumerate(samples):
        for c, r in enumerate(reps):
            if same_orbit(samples[r], s, per_tg, tol=tol) is not None:
                cls.append(c)
                break
        else:
            reps.append(i)
            cls.append(len(reps) - 1)
    return cls


def value_at_clock(prob: Shooting, xs: FloatArr, eps: float, t: float) -> FloatArr:
    """Position-velocity state of the forced orbit at absolute clock ``t``."""
    tt = prob.times()
    tl = prob.t0 + ((t - prob.t0) % prob.period)
    k = int(np.searchsorted(tt, tl, side="right") - 1)
    k = min(max(k, 0), prob.n_seg - 1)
    st = flow(prob.model, eps, xs[k], tt[k], tl)[0][0]
    return model_to_pv(prob.model, st, tl, eps)


# ---------------------------------------------------------------------------
# Family walk in the three-body problem with bifurcation detection.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymmetricMember:
    """A symmetric three-body orbit ``(x0, 0, z0, 0, vy0, 0)`` at period ``period``."""

    x0: float
    z0: float
    vy0: float
    period: float
    residual: float

    @property
    def state(self) -> FloatArr:
        return np.array([self.x0, 0.0, self.z0, 0.0, self.vy0, 0.0])


def _half(model: SunModel, u: FloatArr, period: float, spatial: bool) -> tuple[Any, ...]:
    pv = np.array([u[0], 0.0, u[1] if spatial else 0.0, 0.0, u[-1], 0.0])
    x0 = pv_to_model(model, pv, 0.0, 0.0)
    xf, stm, _ = flow(model, 0.0, x0, 0.0, 0.5 * period, variational=True)
    assert stm is not None
    xf0 = xf[0]
    f = vector_field(model, 0.0, 0.5 * period, xf0)
    # Symmetric crossing conditions in position-velocity: y = vx = vz = 0. In canonical
    # coordinates at eps = 0, vx = px + y and vz = pz, so (y, px, pz) = 0 is equivalent.
    rows = [1, 3, 5] if spatial else [1, 3]
    cols = [0, 2, 4] if spatial else [0, 4]
    # d(model state)/d(pv inputs): in canonical coords px = vx - y, py = vy + x.
    conv = np.eye(6)
    if model.kind == "qbcp":
        conv[3, 1] = -1.0
        conv[4, 0] = 1.0
    jfull = stm[0] @ conv
    return xf0[rows], jfull[np.ix_(rows, cols)], 0.5 * f[rows]


def correct_symmetric(
    model: SunModel, guess: FloatArr, period: float, *, tol: float = 1e-12, max_iter: int = 30
) -> SymmetricMember:
    """Newton on ``(x0, [z0], vy0)`` at fixed period; ``guess = (x0, z0, vy0)``."""
    spatial = abs(float(guess[1])) > 0.0
    u = np.array([guess[0], guess[1], guess[2]]) if spatial else np.array([guess[0], guess[2]])
    res = np.array([np.inf])
    for _ in range(max_iter):
        res, jac, _ = _half(model, u, period, spatial)
        if float(np.max(np.abs(res))) < tol:
            break
        u = u - np.linalg.solve(jac, res)
    return SymmetricMember(
        x0=float(u[0]),
        z0=float(u[1]) if spatial else 0.0,
        vy0=float(u[-1]),
        period=float(period),
        residual=float(np.max(np.abs(res))),
    )


def _member_indices(model: SunModel, m: SymmetricMember) -> dict[str, float | None]:
    parent = Parent("walk", m.state, m.period, 1, 1)
    mono = lap_monodromy(model, parent)
    return stability_indices(mono, planar=m.z0 == 0.0)


def _thresholds(m_max: int) -> list[tuple[float, int, int]]:
    out = [(2.0, 1, 0)]
    for mm in range(2, m_max + 1):
        for k in range(1, mm // 2 + 1):
            if math.gcd(k, mm) == 1:
                out.append((2.0 * math.cos(2.0 * math.pi * k / mm), mm, k))
    return out


def walk_family(
    model: SunModel,
    seed: FloatArr,
    seed_period: float,
    targets: list[float],
    direction: float,
    *,
    ds: float = 0.005,
    ds_max: float = 0.03,
    max_steps: int = 400,
    tol: float = 1e-11,
    t_window: tuple[float, float] = (0.5, 60.0),
    m_max: int = 6,
    wall_s: float = math.inf,
    log: Callable[[str], None] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Pseudo-arclength along a symmetric family in ``(x0, [z0], vy0, T)``, collecting members.

    Every crossing of a target period is corrected there. Along the walk the stability indices
    are tracked and every crossing of ``2 cos(2 pi k / m)`` (``m <= m_max``) is recorded as a
    bifurcation event (``m = 1``: tangent / pitchfork; ``m = 2``: period doubling); a spatial
    family whose ``z0`` changes sign has met the plane (pitchfork, review section 4) and the walk
    stops there. Each member carries ``dT/dC`` from the walk tangent and the list of events
    passed before it. Steps are accepted only if the corrector lands within half a step of the
    predictor and the tangent cosine stays above 0.95.
    """
    spatial = abs(float(seed[1])) > 0.0
    first = correct_symmetric(model, seed, seed_period)
    found: list[dict[str, Any]] = []
    if first.residual > 1e-9:
        return found, {"reason": f"seed did not correct (res {first.residual:.2e})"}
    u = np.array([first.x0, first.z0, first.vy0]) if spatial else np.array([first.x0, first.vy0])
    z = np.concatenate([u, [seed_period]])
    _, jac, dt = _half(model, u, seed_period, spatial)
    tan = np.linalg.svd(np.column_stack([jac, dt]))[2][-1]
    if tan[-1] * direction < 0:
        tan = -tan
    thr = _thresholds(m_max)
    prev_idx = _member_indices(model, first)
    events: list[dict[str, Any]] = []
    t_min = t_max = seed_period
    info: dict[str, Any] = {"reason": "max_steps", "seed_corrected": first.state.tolist()}
    t_start = time.monotonic()
    step = 0

    def unpack(zv: FloatArr) -> SymmetricMember:
        return SymmetricMember(
            float(zv[0]), float(zv[1]) if spatial else 0.0, float(zv[-2]), float(zv[-1]), 0.0
        )

    for step in range(max_steps):
        if time.monotonic() - t_start > wall_s:
            info["reason"] = f"wall-clock limit at T={z[-1]:.6f}"
            break
        pred = z + ds * tan
        zz = pred.copy()
        ok = False
        try:
            for _ in range(8):
                res, jac, dt = _half(model, zz[:-1], zz[-1], spatial)
                if float(np.max(np.abs(res))) < tol:
                    ok = True
                    break
                big = np.vstack([np.column_stack([jac, dt]), tan])
                rhs = np.concatenate([res, [tan @ (zz - pred)]])
                zz = zz - np.linalg.solve(big, rhs)
        except (np.linalg.LinAlgError, RuntimeError, ValueError):
            ok = False
        newtan = tan
        if ok:
            newtan = np.linalg.svd(np.column_stack([jac, dt]))[2][-1]
            if newtan @ tan < 0:
                newtan = -newtan
            if float(np.max(np.abs(zz - pred))) > 0.5 * ds or float(newtan @ tan) < 0.95:
                ok = False
        if not ok:
            ds *= 0.5
            if ds < 1e-7:
                info["reason"] = f"step collapse at T={z[-1]:.6f}"
                break
            continue
        t_old, t_new = float(z[-1]), float(zz[-1])
        z_prev, tan_prev = z, tan
        z, tan = zz, newtan
        t_min, t_max = min(t_min, t_new), max(t_max, t_new)
        mem = unpack(z)
        idx = _member_indices(model, mem)
        for name, s_new in idx.items():
            s_old = prev_idx.get(name)
            if s_new is None or s_old is None:
                continue
            for c, mm, kk in thr:
                if (s_old - c) * (s_new - c) < 0:
                    events.append(
                        {
                            "index": name,
                            "m": mm,
                            "k": kk,
                            "T": [t_old, t_new],
                            "x0": float(z[0]),
                            "s": [s_old, s_new],
                        }
                    )
                    if log is not None:
                        log(f"  bifurcation event {name} m={mm} k={kk} near T={t_new:.6f}")
        prev_idx = idx
        if spatial and z_prev[1] * z[1] <= 0.0:
            events.append(
                {
                    "index": "z0",
                    "m": 1,
                    "k": 0,
                    "T": [t_old, t_new],
                    "x0": float(z[0]),
                    "s": None,
                    "plane": True,
                }
            )
            info["reason"] = f"spatial family met the plane at T={t_new:.6f} (pitchfork)"
            break
        for tp in targets:
            if (t_old - tp) * (t_new - tp) <= 0.0 and t_old != tp:
                w = (tp - t_old) / (t_new - t_old)
                zi = (1 - w) * z_prev + w * z
                guess = np.array([zi[0], zi[1] if spatial else 0.0, zi[-2]])
                try:
                    member = correct_symmetric(model, guess, tp)
                except (np.linalg.LinAlgError, RuntimeError):
                    continue
                near = abs(member.x0 - guess[0]) + abs(member.vy0 - guess[2]) < 1e-3
                if member.residual < 1e-9 and near:
                    tw = (1 - w) * tan_prev + w * tan
                    stw = np.array([zi[0], 0.0, zi[1] if spatial else 0.0, 0.0, zi[-2], 0.0])
                    g = grad_jacobi_pv(stw, model.mu)
                    dcds = g[0] * tw[0] + g[4] * tw[-2] + (g[2] * tw[1] if spatial else 0.0)
                    dtdc = float(tw[-1] / dcds) if dcds != 0.0 else math.inf
                    found.append(
                        {
                            "member": member,
                            "dT_dC": dtdc,
                            "events_before": [dict(e) for e in events],
                            # Only s = 2 crossings (m = 1) can swap the branch at the same
                            # period; period-multiplying ones branch off at m T.
                            "after_branch_point": any(e["m"] == 1 for e in events),
                        }
                    )
        if abs(z[0]) > 5 or abs(z[-2]) > 10 or not (t_window[0] < t_new < t_window[1]):
            info["reason"] = f"left domain at T={t_new:.6f}"
            break
        r1 = math.hypot(z[0] + model.mu, z[1] if spatial else 0.0)
        r2 = math.hypot(z[0] - 1 + model.mu, z[1] if spatial else 0.0)
        if min(r1, r2) < 0.005:
            info["reason"] = f"start point near a primary at T={t_new:.6f}"
            break
        if log is not None and step % 50 == 0:
            log(f"  family step {step}: T={t_new:.6f} x0={z[0]:.6f}")
        ds = min(ds * 1.3, ds_max)
    info.update({"T_min": t_min, "T_max": t_max, "steps": step, "events": events})
    return found, info
