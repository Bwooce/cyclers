"""Ballistic Titania-Oberon flyby arcs in a URA111 force model (#895).

Pre-registration and results: ``docs/notes/2026-10-04-895-titania-oberon-realeph.md``.

Force model
-----------
Uranus-centred (NAIF 799), non-rotating J2000 axes, km, km/s, TDB seconds measured from a model
reference epoch ``t_ref_et`` (SPICE ephemeris time). The spacecraft is massless and feels:

* Uranus as a point mass plus the zonal harmonics J2 and J4 about a fixed pole;
* each body ``k`` of a fixed list (Miranda, Ariel, Umbriel, Titania, Oberon, the Sun) as an
  UNSOFTENED point mass, direct term ``-GM_k (r - r_k)/|r - r_k|^3`` and indirect term
  ``-GM_k rho_k/|rho_k|^3`` (the acceleration of Uranus by that body). ``rho_k = r_k`` except in
  the `#890` provider, where Oberon circles the Uranus-Titania barycentre and its indirect anchor
  is its barycentric position (see :func:`model_890`).

A body's position is ``w_k`` times its URA111 position (cubic Hermite table sampled from SPICE)
plus ``1 - w_k`` times a prescribed circle; ``w_k`` is the homotopy blend of Bradley & Russell
(2014). The constants are the set URA111 was fitted with, read from the kernel's comment area.

Nothing here zeroes or softens a force near a body (the V4 Uranus lanes do; see #890 note 2.1).

Model checks: :func:`moon_positive_control` (each moon reproduced as a test body against
URA111) and :func:`model_890` (the planar `#890` four-body model rebuilt in this frame, against
its stored periodic orbit). Both are run by ``tests/search/test_titania_oberon_realeph_895.py``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numba import njit
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.integrate._ivp import dop853_coefficients as _dop
from scipy.optimize import brentq

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
DAY_S = 86400.0

# --------------------------------------------------------------------------------------------- #
# Constants (URA111 header, R. A. Jacobson, 2014-01-08; checked against the kernel by a test)
# --------------------------------------------------------------------------------------------- #

GM_URANUS = 5.793951322279009e06
GM_MIRANDA = 4.319516899232100e00
GM_ARIEL = 8.346344431770477e01
GM_UMBRIEL = 8.509338094489388e01
GM_TITANIA = 2.269437003741248e02
GM_OBERON = 2.053234302535623e02
GM_SUN = 1.327132332639000e11
J2_URA111 = 3.510685384697763e-03
J4_URA111 = -3.416639735448987e-05
R_REF_KM = 2.555900000000000e04

#: The repository's former Uranus J2 (paired with 25,559 km until #894) and French et al. 2024.
J2_FORMER_REPO = 3.34343e-3
J2_FRENCH_2024 = 3509.291e-6

#: IAU 2009 pole (GMAT PCK BODY799_POLE_RA/DEC, constant) and the URA111 header pole (ZACPL7,
#: ZDEPL7 at POLTIM = J2000, opposite sense).
POLE_IAU_RA_DEC_DEG = (257.311, -15.175)
POLE_URA111_RA_DEC_DEG = (77.30990252631723, 15.17245819840212)

BODY_NAMES: tuple[str, ...] = ("Miranda", "Ariel", "Umbriel", "Titania", "Oberon", "Sun")
NAIF_IDS: tuple[int, ...] = (705, 701, 702, 703, 704, 10)
BODY_GM: tuple[float, ...] = (GM_MIRANDA, GM_ARIEL, GM_UMBRIEL, GM_TITANIA, GM_OBERON, GM_SUN)
#: Mean radii, km (GMAT PCK BODY70x_RADII; Miranda and Ariel triaxial, mean of the axes).
BODY_RADIUS_KM: tuple[float, ...] = (235.8, 578.9, 584.7, 788.9, 761.4, 695700.0)
I_MIRANDA, I_ARIEL, I_UMBRIEL, I_TITANIA, I_OBERON, I_SUN = range(6)
INNER = (I_MIRANDA, I_ARIEL, I_UMBRIEL)
MOONS = (I_MIRANDA, I_ARIEL, I_UMBRIEL, I_TITANIA, I_OBERON)

TABLE_STEP_S = 600.0


def pole_vector(ra_deg: float, dec_deg: float) -> FloatArray:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])


POLE_IAU = pole_vector(*POLE_IAU_RA_DEC_DEG)

# --------------------------------------------------------------------------------------------- #
# Compiled core: body positions, acceleration, gravity gradient, right-hand side, DOP853
# --------------------------------------------------------------------------------------------- #

_A = np.ascontiguousarray(_dop.A[: _dop.N_STAGES, : _dop.N_STAGES])
_B = np.ascontiguousarray(_dop.B)
_C = np.ascontiguousarray(_dop.C[: _dop.N_STAGES])
_E3 = np.ascontiguousarray(_dop.E3)
_E5 = np.ascontiguousarray(_dop.E5)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _hermite(tab: Any, k: int, t0: float, h: float, t: float, out: Any, outv: Any) -> None:
    s = (t - t0) / h
    i = int(math.floor(s))  # noqa: RUF046 (numba typing)
    n = tab.shape[1]
    if i < 0 or i > n - 2:
        if i == n - 1 and s - i < 1e-9:
            i = n - 2
        else:
            raise ValueError("time outside the ephemeris table")
    u = s - i
    u2 = u * u
    u3 = u2 * u
    h00 = 2.0 * u3 - 3.0 * u2 + 1.0
    h10 = (u3 - 2.0 * u2 + u) * h
    h01 = -2.0 * u3 + 3.0 * u2
    h11 = (u3 - u2) * h
    d00 = (6.0 * u2 - 6.0 * u) / h
    d10 = 3.0 * u2 - 4.0 * u + 1.0
    d01 = (-6.0 * u2 + 6.0 * u) / h
    d11 = 3.0 * u2 - 2.0 * u
    for j in range(3):
        p0 = tab[k, i, j]
        v0 = tab[k, i, 3 + j]
        p1 = tab[k, i + 1, j]
        v1 = tab[k, i + 1, 3 + j]
        out[j] = h00 * p0 + h10 * v0 + h01 * p1 + h11 * v1
        outv[j] = d00 * p0 + d10 * v0 + d01 * p1 + d11 * v1


@njit(cache=True)  # type: ignore[untyped-decorator]
def _circle(circ: Any, k: int, t: float, out: Any, outv: Any) -> None:
    a = circ[k, 0]
    n = circ[k, 1]
    lon = circ[k, 2] + n * t
    c = math.cos(lon)
    s = math.sin(lon)
    for j in range(3):
        out[j] = a * (c * circ[k, 3 + j] + s * circ[k, 6 + j])
        outv[j] = a * n * (-s * circ[k, 3 + j] + c * circ[k, 6 + j])


@njit(cache=True)  # type: ignore[untyped-decorator]
def _bodies(
    t: float,
    w: Any,
    tabidx: Any,
    tab: Any,
    tab_t0: float,
    tab_h: float,
    circ: Any,
    cidx: Any,
    pos: Any,
    vel: Any,
    rho: Any,
) -> None:
    """Positions, velocities and indirect anchors of every body at time ``t``."""
    nb = w.shape[0]
    p = np.empty(3)
    v = np.empty(3)
    pc = np.empty(3)
    vc = np.empty(3)
    for k in range(nb):
        wk = w[k]
        for j in range(3):
            pos[k, j] = 0.0
            vel[k, j] = 0.0
            rho[k, j] = 0.0  # holds the barycentre shift until the end of the loop body
        if wk > 0.0:
            _hermite(tab, tabidx[k], tab_t0, tab_h, t, p, v)
            for j in range(3):
                pos[k, j] += wk * p[j]
                vel[k, j] += wk * v[j]
        if wk < 1.0:
            _circle(circ, k, t, p, v)
            for j in range(3):
                pos[k, j] += (1.0 - wk) * p[j]
                vel[k, j] += (1.0 - wk) * v[j]
            beta = circ[k, 9]
            if beta != 0.0:
                _circle(circ, cidx[k], t, pc, vc)
                for j in range(3):
                    pos[k, j] += (1.0 - wk) * beta * pc[j]
                    vel[k, j] += (1.0 - wk) * beta * vc[j]
                    rho[k, j] = (1.0 - wk) * beta * pc[j]
        for j in range(3):
            rho[k, j] = pos[k, j] - rho[k, j]


@njit(cache=True)  # type: ignore[untyped-decorator]
def _zonal(n: int, coef: float, r: Any, pole: Any, a: Any, g: Any, want_grad: bool) -> None:
    """Add the degree-``n`` zonal acceleration (and gradient) for ``coef = GM J_n R^n``."""
    rr = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    u = (r[0] * pole[0] + r[1] * pole[1] + r[2] * pole[2]) / rr
    if n == 2:
        pn = 0.5 * (3.0 * u * u - 1.0)
        dp = 3.0 * u
        ddp = 3.0
    else:
        pn = (35.0 * u**4 - 30.0 * u * u + 3.0) / 8.0
        dp = (140.0 * u**3 - 60.0 * u) / 8.0
        ddp = (420.0 * u * u - 60.0) / 8.0
    f = (n + 1) * pn + u * dp
    fp = (n + 2) * dp + u * ddp
    gg = dp
    gp = ddp
    rm3 = rr ** (-(n + 3))
    rm2 = rr ** (-(n + 2))
    for i in range(3):
        a[i] += coef * (f * rm3 * r[i] - gg * rm2 * pole[i])
    if want_grad:
        rhat = r / rr
        du = (pole - u * rhat) / rr
        for i in range(3):
            for j in range(3):
                d_ij = 1.0 if i == j else 0.0
                term = (
                    rm3 * (f * d_ij + r[i] * fp * du[j]) - (n + 3) * f * rm3 / rr * r[i] * rhat[j]
                )
                term += -gp * rm2 * pole[i] * du[j] + (n + 2) * gg * rm3 * pole[i] * rhat[j]
                g[i, j] += coef * term


@njit(cache=True)  # type: ignore[untyped-decorator]
def _accel(t: float, r: Any, prm: Any, a: Any, g: Any, want_grad: bool) -> None:
    gm, w, tabidx, tab, tab_t0, tab_h, circ, cidx, central = prm
    nb = gm.shape[0]
    pos = np.empty((nb, 3))
    vel = np.empty((nb, 3))
    rho = np.empty((nb, 3))
    _bodies(t, w, tabidx, tab, tab_t0, tab_h, circ, cidx, pos, vel, rho)
    gmc = central[0]
    rr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]
    rr = math.sqrt(rr2)
    r3 = rr2 * rr
    for i in range(3):
        a[i] = -gmc * r[i] / r3
    if want_grad:
        r5 = r3 * rr2
        for i in range(3):
            for j in range(3):
                d_ij = 1.0 if i == j else 0.0
                g[i, j] = -gmc * (d_ij / r3 - 3.0 * r[i] * r[j] / r5)
    pole = central[4:7]
    if central[1] != 0.0:
        _zonal(2, gmc * central[1] * central[3] ** 2, r, pole, a, g, want_grad)
    if central[2] != 0.0:
        _zonal(4, gmc * central[2] * central[3] ** 4, r, pole, a, g, want_grad)
    d = np.empty(3)
    for k in range(nb):
        gk = gm[k]
        if gk == 0.0:
            continue
        for j in range(3):
            d[j] = r[j] - pos[k, j]
        dd2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        dd = math.sqrt(dd2)
        d3 = dd2 * dd
        q2 = rho[k, 0] ** 2 + rho[k, 1] ** 2 + rho[k, 2] ** 2
        q3 = q2 * math.sqrt(q2)
        for i in range(3):
            a[i] += -gk * d[i] / d3 - gk * rho[k, i] / q3
        if want_grad:
            d5 = d3 * dd2
            for i in range(3):
                for j in range(3):
                    d_ij = 1.0 if i == j else 0.0
                    g[i, j] += -gk * (d_ij / d3 - 3.0 * d[i] * d[j] / d5)


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs(t: float, y: Any, prm: Any, out: Any) -> None:
    a = np.empty(3)
    g = np.empty((3, 3))
    stm = y.shape[0] > 6
    _accel(t, y[:3], prm, a, g, stm)
    for i in range(3):
        out[i] = y[3 + i]
        out[3 + i] = a[i]
    if stm:
        # Phi row-major 6x6 in y[6:42]; dPhi = [[0, I], [G, 0]] Phi
        for i in range(3):
            for j in range(6):
                out[6 + 6 * i + j] = y[6 + 6 * (3 + i) + j]
        for i in range(3):
            for j in range(6):
                s = 0.0
                for m in range(3):
                    s += g[i, m] * y[6 + 6 * m + j]
                out[6 + 6 * (3 + i) + j] = s


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs_ret(t: float, y: Any, prm: Any) -> Any:
    out = np.empty_like(y)
    _rhs(t, y, prm, out)
    return out


@njit(cache=True)  # type: ignore[untyped-decorator]
def _err_scale(y: Any, yn: Any, rtol: float, atol_r: float, atol_v: float, i: int) -> Any:
    m = max(abs(y[i]), abs(yn[i]))
    return (atol_r if i < 3 else atol_v) + rtol * m


@njit(cache=True)  # type: ignore[untyped-decorator]
def _dop853(
    t0: float,
    t1: float,
    y0: Any,
    prm: Any,
    rtol: float,
    atol_r: float,
    atol_v: float,
    max_steps: int,
    rec_t: Any,
    rec_y: Any,
) -> Any:
    """Adaptive DOP853 (scipy's coefficients and error estimator) on the six state components.

    Returns ``(y1, n_steps, n_rhs, status)``; ``status`` 0 = reached ``t1``, 1 = step limit.
    ``rec_t``/``rec_y`` (length >= max_steps + 1) receive the accepted steps' times and states
    (positions and velocities) when their length is positive.
    """
    n = y0.shape[0]
    record = rec_t.shape[0] > 0
    direction = 1.0 if t1 >= t0 else -1.0
    t = t0
    y = y0.copy()
    f = np.empty(n)
    _rhs(t, y, prm, f)
    nrhs = 1
    # initial step (Hairer's heuristic on the state components)
    d0 = 0.0
    d1 = 0.0
    for i in range(6):
        sc = _err_scale(y, y, rtol, atol_r, atol_v, i)
        d0 += (y[i] / sc) ** 2
        d1 += (f[i] / sc) ** 2
    d0 = math.sqrt(d0 / 6.0)
    d1 = math.sqrt(d1 / 6.0)
    h0 = 1e-6 if (d0 < 1e-5 or d1 < 1e-5) else 0.01 * d0 / d1
    h0 = min(h0, abs(t1 - t0))
    y1 = y + direction * h0 * f
    f1 = np.empty(n)
    _rhs(t + direction * h0, y1, prm, f1)
    nrhs += 1
    d2 = 0.0
    for i in range(6):
        sc = _err_scale(y, y, rtol, atol_r, atol_v, i)
        d2 += ((f1[i] - f[i]) / sc) ** 2
    d2 = math.sqrt(d2 / 6.0) / h0
    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1.0 / 8.0)
    h = min(100.0 * h0, h1, abs(t1 - t0))
    k = np.empty((13, n))
    ytmp = np.empty(n)
    ynew = np.empty(n)
    fnew = np.empty(n)
    nsteps = 0
    if record:
        rec_t[0] = t
        for i in range(6):
            rec_y[0, i] = y[i]
    rejected = False
    while direction * (t1 - t) > 0.0:
        if nsteps >= max_steps:
            return y, nsteps, nrhs, 1
        hs = min(h, abs(t1 - t))
        last = abs(t1 - t) - hs <= 1e-12 * max(1.0, abs(t1))
        hd = direction * (abs(t1 - t) if last else hs)
        for i in range(n):
            k[0, i] = f[i]
        for s in range(1, 12):
            for i in range(n):
                acc = 0.0
                for j in range(s):
                    acc += _A[s, j] * k[j, i]
                ytmp[i] = y[i] + hd * acc
            _rhs(t + _C[s] * hd, ytmp, prm, k[s])
        for i in range(n):
            acc = 0.0
            for j in range(12):
                acc += _B[j] * k[j, i]
            ynew[i] = y[i] + hd * acc
        _rhs(t + hd, ynew, prm, fnew)
        nrhs += 12
        for i in range(n):
            k[12, i] = fnew[i]
        e5 = 0.0
        e3 = 0.0
        for i in range(6):
            sc = _err_scale(y, ynew, rtol, atol_r, atol_v, i)
            s5 = 0.0
            s3 = 0.0
            for j in range(13):
                s5 += _E5[j] * k[j, i]
                s3 += _E3[j] * k[j, i]
            e5 += (s5 / sc) ** 2
            e3 += (s3 / sc) ** 2
        err = 0.0 if e5 == 0.0 and e3 == 0.0 else abs(hd) * e5 / math.sqrt((e5 + 0.01 * e3) * 6.0)
        if err < 1.0:
            t = t1 if last else t + hd
            for i in range(n):
                y[i] = ynew[i]
                f[i] = fnew[i]
            nsteps += 1
            if record:
                rec_t[nsteps] = t
                for i in range(6):
                    rec_y[nsteps, i] = y[i]
            fac = 10.0 if err == 0.0 else min(10.0, 0.9 * err ** (-1.0 / 8.0))
            if rejected:
                fac = min(1.0, fac)
            h = abs(hd) * fac
            rejected = False
        else:
            h = abs(hd) * max(0.2, 0.9 * err ** (-1.0 / 8.0))
            rejected = True
    return y, nsteps, nrhs, 0


# --------------------------------------------------------------------------------------------- #
# Model container
# --------------------------------------------------------------------------------------------- #


@dataclass
class EphemerisTable:
    """URA111 positions and velocities of the six bodies relative to Uranus on a uniform grid."""

    t0: float  # seconds from the model reference epoch
    h: float
    data: FloatArray  # (6, N, 6)
    t_ref_et: float

    @property
    def t1(self) -> float:
        return float(self.t0 + self.h * (self.data.shape[1] - 1))


def _empty_table() -> EphemerisTable:
    return EphemerisTable(t0=0.0, h=1.0, data=np.zeros((6, 2, 6)), t_ref_et=0.0)


@dataclass
class ForceModel:
    """One force model: central body, six bodies with GM, blend weight, table and circle."""

    gm: FloatArray  # (6,)
    w: FloatArray  # (6,) weight of the URA111 table position
    circ: FloatArray  # (6, 10): a, n, lon0, ex(3), ey(3), beta
    cidx: IntArray  # (6,) circle whose position times beta is this body's circle centre
    central: FloatArray  # gm, J2, J4, R_ref, pole(3)
    table: EphemerisTable = field(default_factory=_empty_table)
    t_ref_et: float = 0.0
    label: str = ""

    def params(self) -> tuple[Any, ...]:
        return (
            np.ascontiguousarray(self.gm, dtype=np.float64),
            np.ascontiguousarray(self.w, dtype=np.float64),
            np.arange(6, dtype=np.int64),
            self.table.data,
            float(self.table.t0),
            float(self.table.h),
            np.ascontiguousarray(self.circ, dtype=np.float64),
            np.ascontiguousarray(self.cidx, dtype=np.int64),
            np.ascontiguousarray(self.central, dtype=np.float64),
        )

    def body_states(self, t: float) -> tuple[FloatArray, FloatArray]:
        """Positions and velocities (6, 3) of the six bodies relative to Uranus."""
        p = self.params()
        pos = np.empty((6, 3))
        vel = np.empty((6, 3))
        rho = np.empty((6, 3))
        _bodies(float(t), p[1], p[2], p[3], p[4], p[5], p[6], p[7], pos, vel, rho)
        return pos, vel

    def anchors(self, t: float) -> FloatArray:
        p = self.params()
        pos = np.empty((6, 3))
        vel = np.empty((6, 3))
        rho = np.empty((6, 3))
        _bodies(float(t), p[1], p[2], p[3], p[4], p[5], p[6], p[7], pos, vel, rho)
        return rho

    def accel(self, t: float, r: Sequence[float] | FloatArray) -> FloatArray:
        a = np.empty(3)
        g = np.empty((3, 3))
        _accel(float(t), np.asarray(r, dtype=np.float64), self.params(), a, g, False)
        return a

    def accel_grad(
        self, t: float, r: Sequence[float] | FloatArray
    ) -> tuple[FloatArray, FloatArray]:
        a = np.empty(3)
        g = np.empty((3, 3))
        _accel(float(t), np.asarray(r, dtype=np.float64), self.params(), a, g, True)
        return a, g

    def rhs(self) -> Callable[[float, FloatArray], FloatArray]:
        prm = self.params()

        def f(t: float, y: FloatArray) -> FloatArray:
            return _rhs_ret(float(t), np.asarray(y, dtype=np.float64), prm)  # type: ignore[no-any-return]

        return f


def central_vector(
    gm: float, j2: float, j4: float, pole: FloatArray, r_ref: float = R_REF_KM
) -> FloatArray:
    return np.array([gm, j2, j4, r_ref, pole[0], pole[1], pole[2]], dtype=np.float64)


# --------------------------------------------------------------------------------------------- #
# Propagation
# --------------------------------------------------------------------------------------------- #


@dataclass
class Propagation:
    state: FloatArray
    stm: FloatArray | None
    n_steps: int
    n_rhs: int
    rec_t: FloatArray | None = None
    rec_y: FloatArray | None = None
    complete: bool = True


def propagate(
    model: ForceModel,
    t0: float,
    t1: float,
    y0: Sequence[float] | FloatArray,
    *,
    stm: bool = False,
    rtol: float = 1e-12,
    atol_r: float = 1e-8,
    atol_v: float = 1e-14,
    max_steps: int = 200_000,
    record: bool = False,
    allow_incomplete: bool = False,
) -> Propagation:
    """Integrate with the module's DOP853 from ``t0`` to ``t1`` (either direction).

    With ``allow_incomplete`` a step-limit stop returns the state reached (``complete`` False)
    instead of raising; the step limit is how a collision course ends (the steps shrink without
    bound at an unsoftened point mass).
    """
    y = np.asarray(y0, dtype=np.float64)[:6]
    if stm:
        y = np.concatenate([y, np.eye(6).reshape(-1)])
    if record:
        rt = np.empty(max_steps + 1)
        ry = np.empty((max_steps + 1, 6))
    else:
        rt = np.empty(0)
        ry = np.empty((0, 6))
    y1, ns, nr, status = _dop853(
        float(t0), float(t1), y, model.params(), rtol, atol_r, atol_v, max_steps, rt, ry
    )
    if status != 0 and not allow_incomplete:
        raise RuntimeError(f"DOP853 step limit at t={t0}..{t1}")
    return Propagation(
        complete=status == 0,
        state=np.array(y1[:6]),
        stm=np.array(y1[6:]).reshape(6, 6) if stm else None,
        n_steps=int(ns),
        n_rhs=int(nr),
        rec_t=rt[: ns + 1].copy() if record else None,
        rec_y=ry[: ns + 1].copy() if record else None,
    )


def propagate_scipy(
    model: ForceModel,
    t0: float,
    t1: float,
    y0: Sequence[float] | FloatArray,
    *,
    method: str = "DOP853",
    rtol: float = 1e-13,
    atol_r: float = 1e-9,
    atol_v: float = 1e-15,
    dense: bool = False,
) -> Any:
    """Integrate the same right-hand side with a scipy integrator (verification path)."""
    atol = np.array([atol_r] * 3 + [atol_v] * 3)
    sol = solve_ivp(  # type: ignore[call-overload]
        model.rhs(),
        (float(t0), float(t1)),
        np.asarray(y0, dtype=np.float64)[:6],
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=dense,
    )
    if not sol.success:
        raise RuntimeError(f"{method} failed: {sol.message}")
    return sol


# --------------------------------------------------------------------------------------------- #
# The #890 planar four-body model in this frame (kernel-free)
# --------------------------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Model890:
    """Constants of the `#890` model (registry values, as `two_moon_periodic_890.TwoMoonModel`)."""

    gm_sys: float
    gm_t: float
    gm_o: float
    a_t: float
    a_o: float

    @property
    def gm_planet(self) -> float:
        return self.gm_sys - self.gm_t - self.gm_o

    @property
    def n_t(self) -> float:
        return math.sqrt((self.gm_planet + self.gm_t) / self.a_t**3)

    @property
    def mu(self) -> float:
        return self.gm_t / (self.gm_planet + self.gm_t)

    @property
    def mu_o(self) -> float:
        return self.gm_o / (self.gm_planet + self.gm_t)

    @property
    def omega_o(self) -> float:
        """Oberon's rate in the Titania frame, nondimensional (`ccr4bp.two_body_synodic_rate`)."""
        a = self.a_o / self.a_t
        return math.sqrt((1.0 - self.mu + self.mu_o) / a**3) - 1.0

    @property
    def n_o(self) -> float:
        return (1.0 + self.omega_o) * self.n_t

    @property
    def cycle_s(self) -> float:
        """Five forcing periods, seconds."""
        return 5.0 * 2.0 * math.pi / abs(self.omega_o) / self.n_t


def registry_890() -> Model890:
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    return Model890(
        gm_sys=float(PRIMARIES["Uranus"]),
        gm_t=float(SATELLITES["Titania"].mu_km3_s2),
        gm_o=float(SATELLITES["Oberon"].mu_km3_s2),
        a_t=float(SATELLITES["Titania"].sma_km),
        a_o=float(SATELLITES["Oberon"].sma_km),
    )


def model_890(c: Model890, *, oberon_about_barycentre: bool = True) -> ForceModel:
    """The `#890` planar model as a Uranus-centred force model.

    Titania circles Uranus at ``a_t`` with rate ``n_t``; Oberon circles the Uranus-Titania
    barycentre (which sits at ``mu r_T`` from Uranus) at ``a_o`` with rate ``n_o``. Uranus's
    acceleration (the indirect term) is ``GM_T r_T/a_t^3`` (its motion about the barycentre) plus
    ``GM_O R_O/a_o^3`` (the barycentre pulled by Oberon). With ``oberon_about_barycentre=False``
    Oberon circles Uranus instead (the negative control).
    """
    circ = np.zeros((6, 10))
    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    circ[I_TITANIA, :3] = (c.a_t, c.n_t, 0.0)
    circ[I_OBERON, :3] = (c.a_o, c.n_o, 0.0)
    for k in (I_TITANIA, I_OBERON):
        circ[k, 3:6] = ex
        circ[k, 6:9] = ey
    for k in (I_MIRANDA, I_ARIEL, I_UMBRIEL, I_SUN):
        circ[k, :3] = (1.0, 0.0, 0.0)
        circ[k, 3:6] = ex
        circ[k, 6:9] = ey
    cidx = np.zeros(6, dtype=np.int64)
    if oberon_about_barycentre:
        circ[I_OBERON, 9] = c.mu
        cidx[I_OBERON] = I_TITANIA
    gm = np.zeros(6)
    gm[I_TITANIA] = c.gm_t
    gm[I_OBERON] = c.gm_o
    return ForceModel(
        gm=gm,
        w=np.zeros(6),
        circ=circ,
        cidx=cidx,
        central=central_vector(c.gm_planet, 0.0, 0.0, np.array([0.0, 0.0, 1.0])),
        label="890",
    )


def state_890_to_inertial(c: Model890, s4: Sequence[float] | FloatArray) -> FloatArray:
    """`#890` rotating nondimensional state at tau = 0 -> Uranus-centred inertial 6-state."""
    x, y, vx, vy = (float(v) for v in s4)
    xu = x + c.mu  # Uranus sits at (-mu, 0) in the barycentric rotating frame
    length = c.a_t
    vel = c.a_t * c.n_t
    return np.array([xu * length, y * length, 0.0, (vx - y) * vel, (vy + xu) * vel, 0.0])


def inertial_to_rot_890(c: Model890, t: float, s6: FloatArray) -> FloatArray:
    """Uranus-centred inertial state at time ``t`` -> `#890` barycentric rotating nondim 4-state."""
    th = c.n_t * t
    cs, sn = math.cos(th), math.sin(th)
    r = np.array([cs * s6[0] + sn * s6[1], -sn * s6[0] + cs * s6[1]])
    v_in = np.array([cs * s6[3] + sn * s6[4], -sn * s6[3] + cs * s6[4]])
    v = v_in - c.n_t * np.array([-r[1], r[0]])
    return np.array(
        [r[0] / c.a_t - c.mu, r[1] / c.a_t, v[0] / (c.a_t * c.n_t), v[1] / (c.a_t * c.n_t)]
    )


# --------------------------------------------------------------------------------------------- #
# Encounters
# --------------------------------------------------------------------------------------------- #


def hill_radius_km(k: int, a_km: float) -> float:
    return float(a_km * (BODY_GM[k] / (3.0 * GM_URANUS)) ** (1.0 / 3.0))


@dataclass
class Encounter:
    body: int
    t: float
    dist_km: float
    rel_pos: FloatArray
    rel_vel: FloatArray

    @property
    def name(self) -> str:
        return BODY_NAMES[self.body]

    def osculating(self, gm: float | None = None) -> dict[str, float]:
        mu = BODY_GM[self.body] if gm is None else gm
        r = float(np.linalg.norm(self.rel_pos))
        v = float(np.linalg.norm(self.rel_vel))
        energy = 0.5 * v * v - mu / r
        vinf2 = v * v - 2.0 * mu / r
        rv = float(np.dot(self.rel_pos, self.rel_vel))
        e = float(np.linalg.norm(((v * v - mu / r) * self.rel_pos - rv * self.rel_vel) / mu))
        turn = 2.0 * math.degrees(math.asin(1.0 / e)) if e > 1.0 else float("nan")
        return {
            "altitude_km": r - BODY_RADIUS_KM[self.body],
            "speed_rel_kms": v,
            "energy_km2s2": energy,
            "vinf_kms": math.sqrt(vinf2) if vinf2 > 0 else float("nan"),
            "ecc": e,
            "turn_deg": turn,
        }


def find_minima(
    traj: Callable[[float], FloatArray],
    body_state: Callable[[float], tuple[FloatArray, FloatArray]],
    t0: float,
    t1: float,
    *,
    dt: float = 600.0,
    body: int = -1,
    max_dist: float = math.inf,
) -> list[Encounter]:
    """Local minima of the distance between a trajectory and a body, by root-finding on range rate.

    ``traj(t)`` returns the 6-state; ``body_state(t)`` the body's position and velocity.
    """

    def rdot(t: float) -> float:
        s = traj(t)
        p, v = body_state(t)
        return float(np.dot(s[:3] - p, s[3:6] - v))

    ts = np.arange(t0, t1, dt)
    ts = np.append(ts, t1)
    vals = np.array([rdot(t) for t in ts])
    out: list[Encounter] = []
    for i in range(len(ts) - 1):
        if vals[i] < 0.0 <= vals[i + 1]:
            tm = brentq(rdot, ts[i], ts[i + 1], xtol=1e-6, rtol=1e-15)
            s = traj(tm)
            p, v = body_state(tm)
            d = float(np.linalg.norm(s[:3] - p))
            if d <= max_dist:
                out.append(Encounter(body, tm, d, s[:3] - p, s[3:6] - v))
    return out


def hermite_traj(
    ts: FloatArray, ys: FloatArray, model: ForceModel
) -> Callable[[float], FloatArray]:
    """Quintic Hermite interpolant of recorded steps (positions, velocities, accelerations)."""
    acc = np.array([model.accel(t, y[:3]) for t, y in zip(ts, ys, strict=True)])

    def f(t: float) -> FloatArray:
        i = int(np.searchsorted(ts, t, side="right") - 1)
        i = min(max(i, 0), len(ts) - 2)
        h = ts[i + 1] - ts[i]
        u = (t - ts[i]) / h
        p0, p1 = ys[i, :3], ys[i + 1, :3]
        v0, v1 = ys[i, 3:] * h, ys[i + 1, 3:] * h
        a0, a1 = acc[i] * h * h, acc[i + 1] * h * h
        u2, u3, u4, u5 = u * u, u**3, u**4, u**5
        b = (
            (
                1 - 10 * u3 + 15 * u4 - 6 * u5,
                u - 6 * u3 + 8 * u4 - 3 * u5,
                0.5 * (u2 - 3 * u3 + 3 * u4 - u5),
            ),
            (10 * u3 - 15 * u4 + 6 * u5, -4 * u3 + 7 * u4 - 3 * u5, 0.5 * (u3 - 2 * u4 + u5)),
        )
        db = (
            (
                -30 * u2 + 60 * u3 - 30 * u4,
                1 - 18 * u2 + 32 * u3 - 15 * u4,
                0.5 * (2 * u - 9 * u2 + 12 * u3 - 5 * u4),
            ),
            (
                30 * u2 - 60 * u3 + 30 * u4,
                -12 * u2 + 28 * u3 - 15 * u4,
                0.5 * (3 * u2 - 8 * u3 + 5 * u4),
            ),
        )
        pos = (
            b[0][0] * p0 + b[0][1] * v0 + b[0][2] * a0 + b[1][0] * p1 + b[1][1] * v1 + b[1][2] * a1
        )
        vel = (
            db[0][0] * p0
            + db[0][1] * v0
            + db[0][2] * a0
            + db[1][0] * p1
            + db[1][1] * v1
            + db[1][2] * a1
        ) / h
        return np.concatenate([pos, vel])

    return f


# --------------------------------------------------------------------------------------------- #
# SPICE (kernel-dependent; spiceypy imported lazily)
# --------------------------------------------------------------------------------------------- #

LSK_PATH = Path(__file__).resolve().parents[1] / "verify" / "kernels" / "naif0012.tls"


def ura_path() -> Path:
    from cyclerfinder.data.validation.v4_uranus_strict import DEFAULT_URA_PATH

    return Path(DEFAULT_URA_PATH)


def kernels_present() -> bool:
    try:
        import spiceypy  # noqa: F401
    except ImportError:
        return False
    return ura_path().exists() and LSK_PATH.exists()


_FURNISHED = False


def load_kernels() -> Any:
    """Furnish the URA111 kernel and the leap seconds once; return the spiceypy module."""
    global _FURNISHED
    import spiceypy as sp

    if not _FURNISHED:
        sp.furnsh(str(LSK_PATH))
        sp.furnsh(str(ura_path()))
        _FURNISHED = True
    return sp


def et_of(utc_or_tdb: str) -> float:
    """SPICE ephemeris time of a calendar string (``"... TDB"`` is read as TDB)."""
    return float(load_kernels().str2et(utc_or_tdb))


def tdb_string(et: float) -> str:
    return str(load_kernels().timout(et, "YYYY-MM-DD HR:MN:SC.### ::TDB"))


def spice_state(naif: int, et: float) -> FloatArray:
    sp = load_kernels()
    st, _ = sp.spkgeo(naif, float(et), "J2000", 799)
    return np.asarray(st, dtype=np.float64)


def kernel_header_text() -> str:
    sp = load_kernels()
    handle = sp.dafopr(str(ura_path()))
    try:
        lines: list[str] = []
        while True:
            _n, buf, done = sp.dafec(handle, 100)
            lines.extend(buf)
            if done:
                break
    finally:
        sp.dafcls(handle)
    return "\n".join(lines)


def build_table(
    t_ref_et: float,
    t0: float,
    t1: float,
    *,
    step: float = TABLE_STEP_S,
    cache_dir: Path | None = None,
) -> EphemerisTable:
    """Sample the six bodies' URA111 states relative to Uranus on ``[t0, t1]`` (model seconds)."""
    n = math.ceil((t1 - t0) / step) + 1
    key = f"tab_{t_ref_et:.3f}_{t0:.1f}_{n}_{step:.1f}.npy"
    if cache_dir is not None and (cache_dir / key).exists():
        data = np.load(cache_dir / key)
    else:
        sp = load_kernels()
        ets = t_ref_et + t0 + step * np.arange(n)
        data = np.empty((6, n, 6))
        for k, naif in enumerate(NAIF_IDS):
            for i, et in enumerate(ets):
                data[k, i] = sp.spkgeo(naif, float(et), "J2000", 799)[0]
        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
            np.save(cache_dir / key, data)
    return EphemerisTable(t0=t0, h=step, data=np.ascontiguousarray(data), t_ref_et=t_ref_et)


def full_model(
    table: EphemerisTable,
    *,
    j2: float = J2_URA111,
    j4: float = J4_URA111,
    pole: FloatArray = POLE_IAU,
    sun: bool = True,
    exclude: int | None = None,
    self_mass: bool = True,
) -> ForceModel:
    """The force model of the pre-registration (section 1.1) at lam = 1.

    ``exclude`` turns a body into the test body of the positive control: its GM is removed from
    the perturbers and (``self_mass``) added to the central term.
    """
    gm = np.array(BODY_GM, dtype=np.float64)
    if not sun:
        gm[I_SUN] = 0.0
    gmc = GM_URANUS
    if exclude is not None:
        if self_mass:
            gmc += gm[exclude]
        gm[exclude] = 0.0
    circ = np.zeros((6, 10))
    circ[:, 0] = 1.0
    circ[:, 3] = 1.0
    circ[:, 7] = 1.0
    return ForceModel(
        gm=gm,
        w=np.ones(6),
        circ=circ,
        cidx=np.zeros(6, dtype=np.int64),
        central=central_vector(gmc, j2, j4, pole),
        table=table,
        t_ref_et=table.t_ref_et,
        label="full",
    )


def moon_positive_control(
    table: EphemerisTable,
    moon: int,
    t_start: float,
    days: Sequence[float],
    *,
    rtol: float = 1e-13,
    **model_kwargs: Any,
) -> list[float]:
    """Propagate a moon as a test body from its URA111 state; position errors (km) at ``days``."""
    model = full_model(table, exclude=moon, **model_kwargs)
    y = spice_state(NAIF_IDS[moon], table.t_ref_et + t_start)
    out: list[float] = []
    t = t_start
    for d in sorted(days):
        t1 = t_start + d * DAY_S
        y = propagate(model, t, t1, y, rtol=rtol, atol_r=1e-9, atol_v=1e-15).state
        t = t1
        ref = spice_state(NAIF_IDS[moon], table.t_ref_et + t1)
        out.append(float(np.linalg.norm(y[:3] - ref[:3])))
    return out


# --------------------------------------------------------------------------------------------- #
# Reference circles, conjunctions, homotopy models
# --------------------------------------------------------------------------------------------- #


@dataclass
class Circles:
    """Mean circular coplanar Titania and Oberon fitted to URA111 over a window."""

    ex: FloatArray
    ey: FloatArray
    ez: FloatArray
    a_t: float
    n_t: float
    lon0_t: float
    a_o: float
    n_o: float
    lon0_o: float
    fit_window: tuple[float, float]

    def lon_t(self, t: float) -> float:
        return self.lon0_t + self.n_t * t

    @property
    def synodic_s(self) -> float:
        return 2.0 * math.pi / (self.n_t - self.n_o)

    @property
    def cycle_s(self) -> float:
        return 5.0 * self.synodic_s

    def conjunction_near(self, t: float) -> float:
        m = round(((self.lon0_t - self.lon0_o) + (self.n_t - self.n_o) * t) / (2.0 * math.pi))
        return (2.0 * math.pi * m - (self.lon0_t - self.lon0_o)) / (self.n_t - self.n_o)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ez": self.ez.tolist(),
            "ex": self.ex.tolist(),
            "a_t_km": self.a_t,
            "n_t_rad_s": self.n_t,
            "lon0_t": self.lon0_t,
            "a_o_km": self.a_o,
            "n_o_rad_s": self.n_o,
            "lon0_o": self.lon0_o,
            "period_t_d": 2 * math.pi / self.n_t / DAY_S,
            "period_o_d": 2 * math.pi / self.n_o / DAY_S,
            "synodic_d": self.synodic_s / DAY_S,
            "cycle_d": self.cycle_s / DAY_S,
            "fit_window_s": list(self.fit_window),
        }


def fit_circles(table: EphemerisTable, t0: float, t1: float, step: float = 6 * 3600.0) -> Circles:
    model = full_model(table)
    ts = np.arange(t0, t1, step)
    pos_t = np.empty((len(ts), 3))
    pos_o = np.empty((len(ts), 3))
    hsum = np.zeros(3)
    for i, t in enumerate(ts):
        p, v = model.body_states(t)
        pos_t[i], pos_o[i] = p[I_TITANIA], p[I_OBERON]
        hh = np.cross(p[I_TITANIA], v[I_TITANIA])
        hsum += hh / np.linalg.norm(hh)
    ez = hsum / np.linalg.norm(hsum)
    ex = np.cross([0.0, 0.0, 1.0], ez)
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)

    def fit(pos: FloatArray) -> tuple[float, float, float]:
        lon = np.unwrap(np.arctan2(pos @ ey, pos @ ex))
        n, lon0 = np.polyfit(ts, lon, 1)
        return float(np.mean(np.linalg.norm(pos, axis=1))), float(n), float(lon0)

    a_t, n_t, l_t = fit(pos_t)
    a_o, n_o, l_o = fit(pos_o)
    return Circles(ex, ey, ez, a_t, n_t, l_t, a_o, n_o, l_o, (t0, t1))


def homotopy_model(
    table: EphemerisTable,
    circles: Circles,
    lam: float,
    *,
    j2: float = J2_URA111,
    j4: float = J4_URA111,
    pole: FloatArray = POLE_IAU,
    sun: bool = True,
) -> ForceModel:
    """Pre-registration section 1.3: lam = 0 circular coplanar, lam = 1 the full model."""
    gm = np.array(BODY_GM, dtype=np.float64)
    for k in INNER:
        gm[k] *= lam
    gm[I_SUN] = lam * GM_SUN if sun else 0.0
    w = np.ones(6)
    w[I_TITANIA] = lam
    w[I_OBERON] = lam
    circ = np.zeros((6, 10))
    circ[:, 0] = 1.0
    circ[:, 3:6] = circles.ex
    circ[:, 6:9] = circles.ey
    circ[I_TITANIA, :3] = (circles.a_t, circles.n_t, circles.lon0_t)
    circ[I_OBERON, :3] = (circles.a_o, circles.n_o, circles.lon0_o)
    gmc = GM_URANUS + (1.0 - lam) * sum(BODY_GM[k] for k in INNER)
    return ForceModel(
        gm=gm,
        w=w,
        circ=circ,
        cidx=np.zeros(6, dtype=np.int64),
        central=central_vector(gmc, lam * j2, lam * j4, pole),
        table=table,
        t_ref_et=table.t_ref_et,
        label=f"lam={lam:.6g}",
    )


# --------------------------------------------------------------------------------------------- #
# Seed from the #890 periodic orbit
# --------------------------------------------------------------------------------------------- #


@dataclass
class Orbit890:
    c: Model890
    traj: Callable[[float], FloatArray]

    def rot_state(self, t: float) -> tuple[FloatArray, FloatArray]:
        """Uranus-centred planar state in the frame rotating with Titania, ``t`` mod cycle."""
        tt = t % self.c.cycle_s
        s = self.traj(tt)
        th = self.c.n_t * tt
        cs, sn = math.cos(th), math.sin(th)
        r = np.array([cs * s[0] + sn * s[1], -sn * s[0] + cs * s[1]])
        v = np.array([cs * s[3] + sn * s[4], -sn * s[3] + cs * s[4]])
        return r, v - self.c.n_t * np.array([-r[1], r[0]])


def orbit_890(refined_state: Sequence[float] | FloatArray, c: Model890 | None = None) -> Orbit890:
    c = registry_890() if c is None else c
    model = model_890(c)
    s0 = state_890_to_inertial(c, refined_state)
    p = propagate(model, 0.0, c.cycle_s, s0, rtol=1e-13, atol_r=1e-10, atol_v=1e-16, record=True)
    assert p.rec_t is not None and p.rec_y is not None
    return Orbit890(c, hermite_traj(p.rec_t, p.rec_y, model))


def seed_from_890(orb: Orbit890, circles: Circles, t_c: float, times: FloatArray) -> FloatArray:
    """Map the `#890` orbit onto the fitted circles with its Titania flyby at ``t_c``."""
    s = circles.a_t / orb.c.a_t
    q = orb.c.cycle_s / circles.cycle_s
    out = np.empty((len(times), 6))
    for i, t in enumerate(times):
        tau = ((t - t_c) / circles.cycle_s) * orb.c.cycle_s
        r, v = orb.rot_state(tau)
        r = s * r
        v = s * q * v
        lon = circles.lon_t(t)
        cs, sn = math.cos(lon), math.sin(lon)
        vi = v + circles.n_t * np.array([-r[1], r[0]])
        r2 = np.array([cs * r[0] - sn * r[1], sn * r[0] + cs * r[1]])
        v2 = np.array([cs * vi[0] - sn * vi[1], sn * vi[0] + cs * vi[1]])
        out[i, :3] = r2[0] * circles.ex + r2[1] * circles.ey
        out[i, 3:] = v2[0] * circles.ex + v2[1] * circles.ey
    return out


def arc_times(circles: Circles, t_c: float, n_cycles: int, per_cycle: int = 24) -> FloatArray:
    h = circles.cycle_s / per_cycle
    return np.asarray(t_c - 2.0 * h + h * np.arange(per_cycle * n_cycles + 5), dtype=np.float64)


# --------------------------------------------------------------------------------------------- #
# Multiple shooting
# --------------------------------------------------------------------------------------------- #

L_SCALE_KM = 436300.0
V_SCALE_KMS = 3.64
_D = np.array([L_SCALE_KM] * 3 + [V_SCALE_KMS] * 3)


@dataclass
class ShootEval:
    jumps: FloatArray  # (K, 6): propagated end minus next node
    stms: FloatArray | None  # (K, 6, 6)

    @property
    def max_r(self) -> float:
        return float(np.max(np.linalg.norm(self.jumps[:, :3], axis=1)))

    @property
    def max_v(self) -> float:
        return float(np.max(np.linalg.norm(self.jumps[:, 3:], axis=1)))

    @property
    def scaled(self) -> float:
        return float(np.max(np.abs(self.jumps / _D)))


def shoot_eval(
    model: ForceModel, times: FloatArray, x: FloatArray, *, stm: bool = True, rtol: float = 1e-12
) -> ShootEval:
    k = len(times) - 1
    jumps = np.empty((k, 6))
    stms = np.empty((k, 6, 6)) if stm else None
    for i in range(k):
        p = propagate(model, times[i], times[i + 1], x[i], stm=stm, rtol=rtol)
        jumps[i] = p.state - x[i + 1]
        if stms is not None and p.stm is not None:
            stms[i] = p.stm
    return ShootEval(jumps, stms)


def min_norm_step(ev: ShootEval) -> FloatArray:
    """Minimum-norm Newton step (scaled variables) for the continuity equations."""
    assert ev.stms is not None
    k = ev.jumps.shape[0]
    jac = np.zeros((6 * k, 6 * (k + 1)))
    dinv = 1.0 / _D
    for i in range(k):
        jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = (dinv[:, None] * ev.stms[i]) * _D[None, :]
        jac[6 * i : 6 * i + 6, 6 * i + 6 : 6 * i + 12] = -np.eye(6)
    rhs = -(ev.jumps * dinv).reshape(-1)
    dz, *_ = np.linalg.lstsq(jac, rhs, rcond=None)
    return np.asarray(dz.reshape(k + 1, 6) * _D, dtype=np.float64)


@dataclass
class NewtonResult:
    x: FloatArray
    converged: bool
    history: list[dict[str, float]]
    reason: str = ""


def newton(
    model: ForceModel,
    times: FloatArray,
    x0: FloatArray,
    *,
    tol_r: float = 1e-5,
    tol_v: float = 1e-8,
    max_iter: int = 12,
    rtol: float = 1e-12,
) -> NewtonResult:
    x = x0.copy()
    hist: list[dict[str, float]] = []
    try:
        ev = shoot_eval(model, times, x, rtol=rtol)
    except RuntimeError as exc:
        return NewtonResult(x, False, hist, f"propagation: {exc}")
    for it in range(max_iter + 1):
        hist.append({"iter": it, "max_r_km": ev.max_r, "max_v_kms": ev.max_v, "scaled": ev.scaled})
        if ev.max_r < tol_r and ev.max_v < tol_v:
            return NewtonResult(x, True, hist)
        if it == max_iter:
            break
        dx = min_norm_step(ev)
        alpha = 1.0
        for _ in range(5):
            xt = x + alpha * dx
            try:
                evt = shoot_eval(model, times, xt, rtol=rtol)
            except RuntimeError:
                alpha *= 0.5
                continue
            if evt.scaled < ev.scaled:
                break
            alpha *= 0.5
        else:
            return NewtonResult(x, False, hist, "step halving exhausted")
        x, ev = xt, evt
        hist[-1]["alpha"] = alpha
    return NewtonResult(x, False, hist, "iteration limit")


# --------------------------------------------------------------------------------------------- #
# Encounters along an arc
# --------------------------------------------------------------------------------------------- #


def arc_trajectory(
    model: ForceModel, times: FloatArray, x: FloatArray, *, rtol: float = 1e-12
) -> Callable[[float], FloatArray]:
    """Piecewise trajectory: each node propagated to the next, quintic Hermite between steps."""
    pieces: list[Callable[[float], FloatArray]] = []
    for i in range(len(times) - 1):
        p = propagate(model, times[i], times[i + 1], x[i], rtol=rtol, record=True)
        assert p.rec_t is not None and p.rec_y is not None
        pieces.append(hermite_traj(p.rec_t, p.rec_y, model))

    def f(t: float) -> FloatArray:
        i = int(np.searchsorted(times, t, side="right") - 1)
        i = min(max(i, 0), len(times) - 2)
        return pieces[i](t)

    return f


@dataclass
class FlybyInfo:
    body: int
    t: float
    dist_km: float
    alt_km: float
    side: int  # +1 outside the moon's orbit (away from Uranus), -1 inside
    sense: int  # sign of the flyby angular momentum on the reference +z
    osc: dict[str, float]
    z_km: float  # spacecraft distance from the reference plane

    def to_dict(self, t_ref_et: float | None = None) -> dict[str, Any]:
        d: dict[str, Any] = {
            "body": BODY_NAMES[self.body],
            "t_s": self.t,
            "dist_km": self.dist_km,
            "alt_km": self.alt_km,
            "side": self.side,
            "sense": self.sense,
            "z_ref_km": self.z_km,
        }
        d.update(self.osc)
        return d


def encounters(
    model: ForceModel,
    traj: Callable[[float], FloatArray],
    t0: float,
    t1: float,
    ez: FloatArray,
    bodies: Sequence[int] = (I_TITANIA, I_OBERON),
    hill_factor: float = 2.0,
    dt: float = 1800.0,
) -> list[FlybyInfo]:
    """All local minima of the distance to the given bodies within ``hill_factor`` Hill radii."""
    out: list[FlybyInfo] = []
    a_mean = {I_MIRANDA: 129900.0, I_ARIEL: 190900.0, I_UMBRIEL: 266000.0}
    a_mean.update({I_TITANIA: 436300.0, I_OBERON: 583500.0})
    for k in bodies:
        rh = hill_radius_km(k, a_mean[k])

        def bstate(t: float, k: int = k) -> tuple[FloatArray, FloatArray]:
            p, v = model.body_states(t)
            return p[k], v[k]

        for e in find_minima(traj, bstate, t0, t1, dt=dt, body=k, max_dist=hill_factor * rh):
            pm, _vm = bstate(e.t)
            side = 1 if float(np.dot(e.rel_pos, pm)) > 0 else -1
            sense = 1 if float(np.dot(np.cross(e.rel_pos, e.rel_vel), ez)) > 0 else -1
            z = float(np.dot(traj(e.t)[:3], ez))
            out.append(
                FlybyInfo(
                    k, e.t, e.dist_km, e.dist_km - BODY_RADIUS_KM[k], side, sense, e.osculating(), z
                )
            )
    out.sort(key=lambda f: f.t)
    return out


def branch_signature(flybys: Sequence[FlybyInfo]) -> list[tuple[str, int, int]]:
    return [(BODY_NAMES[f.body], f.side, f.sense) for f in flybys]
