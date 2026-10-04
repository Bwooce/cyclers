"""Second-species (consecutive-collision) arcs of the planar circular restricted
problem at mu = 0 (#899 step 1).

At mu = 0 the small body P3 moves on a Kepler conic about P1 and the secondary P2
(on the unit circle, unit angular rate) has no mass.  A *collision* of P3 with P2
breaks the conic into *arcs*; a piece that begins and ends at a collision is an
orbit with consecutive collisions (Henon 1968).  These arcs are the generating
objects of Poincare's second-species periodic solutions at small mu > 0.

Sources (all digested under ``docs/notes/2026-10-04-digest-*``):

* Henon 1968, Bull. astron. 3(3):377 -- the timing equation (eq. 30), the families
  A, B, C and Tables 1-9.
* Brjuno 1978 (Celest. Mech. 18:9 and 18:51) -- sign conventions, collision
  velocity (4.C), the asymmetric arcs T_N, the e* membership test (Theorem 2.2),
  the curve f of Table II.
* Bruno 1981 (Celest. Mech. 24:255) -- the first-order periapsis parameter W.
* Hitzl & Henon 1977 (Celest. Mech. 15:421 and 16:?) -- the critical-orbit
  function S (eq. 40 of the stability paper); S = 0 at critical orbits.
* Henon 2001, Generating families II part B, chapter 18 -- R-arcs and R-orbits;
  Devaney 1981 -- their count by the baker map.
* Gomez & Olle 1986, Celest. Mech. 39:33 -- the elliptic extension (primaries'
  eccentricity ``e_p``, sign ``eps_p``); only the timing equation, the arc
  elements, the parabolic arc and the existence rule of the C families are
  implemented for ``e_p > 0``.  The paper has no Jacobi integral, so anything
  that needs the velocity at collision raises ``NotImplementedError`` there.

Conventions (Henon 1968): ``eps = +1`` if the apse S of P3 on the x axis has
positive abscissa, ``eps1`` (the paper's eps') is +1 for direct motion, ``eps2``
(eps'') is +1 if S is a pericentre.  ``sigma = eps * eps2`` is the only
combination in the timing equation; it is constant along a family: -1 for A,
+1 for B, ``(-1)**(i+j)`` for C_ij.  Brjuno's ``eps`` equals Henon's ``eps``
times ``eps2`` (so Brjuno's eps is our ``sigma``).

Units: P2 on the circle of radius 1 at angular rate 1; the Jacobi constant is
``C = 2 eps1 sqrt(a (1 - e^2)) + 1/a`` and the speed relative to P2 at collision
is ``V = sqrt(3 - C)``.

Never expected to be "validated" beyond the sourced tables in
``tests/search/test_second_species_arcs*.py``; this module is a mu = 0
generator, not a periodic-orbit solver at mu > 0.

Notes on the R-region solvers.  The recurrence ``y[i-1] - 2 y[i] + y[i+1] +
1/y[i] = 0`` is the stationarity condition of
``Phi(y) = 1/2 sum (y[i+1] - y[i])^2 - sum ln|y[i]|``.  Inside one open orthant
(a fixed sign code) Phi is strictly convex, and it is coercive unless all signs
agree (periodic case) so there is exactly one critical point per admissible sign
code: ``2**(n-1)`` R-arcs and ``2**n - 2`` R-orbits (the Devaney count, obtained
here without the baker map).  The solvers are damped Newton iterations on that
convex function: deterministic, one root per code, never a random multi-start.
Because the minimisation is well conditioned (unlike the shooting recurrence,
which doubles errors per step), double precision suffices to well beyond n = 12;
``refine_r_orbit_mp`` polishes a root at extended precision when wanted.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

PI = math.pi
TWO_PI = 2.0 * math.pi

_TOL_ZERO = 1e-9  # default threshold for "indeterminate" flags


# --------------------------------------------------------------------------------------
# Primary's mean anomaly (elliptic extension, Gomez & Olle 1986 eqs. 13c, 16)
# --------------------------------------------------------------------------------------
def primary_mean_anomaly(tau: float, e_p: float = 0.0, eps_p: int = 1) -> float:
    """Mean anomaly of P2 at true anomaly ``tau`` (eq. 13c left side).

    ``M = E - eps_p e_p sin E`` with ``tan(E/2) = K^(-eps_p) tan(tau/2)`` and
    ``K = sqrt((1 + e_p)/(1 - e_p))``; the branch of E follows tau continuously.
    For ``e_p = 0`` this is ``tau``.
    """
    if e_p == 0.0:
        return tau
    if not 0.0 <= e_p < 1.0:
        raise ValueError("e_p must lie in [0, 1)")
    k = math.sqrt((1.0 + e_p) / (1.0 - e_p)) ** (-eps_p)
    n = math.floor((tau + PI) / TWO_PI)
    tr = tau - TWO_PI * n
    ecc_anom = 2.0 * math.atan(k * math.tan(tr / 2.0)) + TWO_PI * n
    return ecc_anom - eps_p * e_p * math.sin(ecc_anom)


def primary_radius(tau: float, e_p: float = 0.0, eps_p: int = 1) -> float:
    """r_p of eq. 16: ``(1 - e_p^2) / (1 + eps_p e_p cos tau)``."""
    return (1.0 - e_p * e_p) / (1.0 + eps_p * e_p * math.cos(tau))


# --------------------------------------------------------------------------------------
# Timing equation (Henon 1968 eq. 30; Gomez & Olle 1986 eq. 17)
# --------------------------------------------------------------------------------------
def timing_residual(tau: float, eta: float, sigma: int, e_p: float = 0.0, eps_p: int = 1) -> float:
    """Residual of the implicit equation relating tau and eta.

    ``sigma = eps * eps2``.  At ``e_p = 0`` this is Henon's eq. 30::

        sqrt(D) [eta D - sin(eta)(cos(eta) - s cos(tau))] - tau |sin(eta)|^3

    with ``D = 1 - s cos(tau) cos(eta)`` and ``s = sigma`` (``s = eps_p sigma``
    in the elliptic form, eq. 17 with ``r_p^(3/2)`` and the mean anomaly of P2).
    """
    sg = eps_p * sigma
    cos_t = math.cos(tau)
    cos_h = math.cos(eta)
    d = 1.0 - sg * cos_t * cos_h
    r_p = primary_radius(tau, e_p, eps_p)
    bracket = eta * d - math.sin(eta) * (cos_h - sg * cos_t)
    mp = primary_mean_anomaly(tau, e_p, eps_p)
    return float(r_p**1.5 * math.sqrt(max(d, 0.0)) * bracket - mp * abs(math.sin(eta)) ** 3)


def _timing_residual_vec(
    tau: float, etas: NDArray[np.float64], sigma: int, e_p: float, eps_p: int
) -> NDArray[np.float64]:
    sg = eps_p * sigma
    cos_t = math.cos(tau)
    cos_h = np.cos(etas)
    d = 1.0 - sg * cos_t * cos_h
    r_p = primary_radius(tau, e_p, eps_p)
    bracket = etas * d - np.sin(etas) * (cos_h - sg * cos_t)
    mp = primary_mean_anomaly(tau, e_p, eps_p)
    out: NDArray[np.float64] = (
        r_p**1.5 * np.sqrt(np.maximum(d, 0.0)) * bracket - mp * np.abs(np.sin(etas)) ** 3
    )
    return out


# --------------------------------------------------------------------------------------
# Arc elements
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SArc:
    """One symmetric arc: elements at the central apse S and at the collisions.

    ``eps1 = 0`` means undefined (rectilinear, or tangent ellipse where both
    ``sin tau`` and ``sin eta`` vanish and either sense is a solution).
    """

    tau: float
    eta: float
    sigma: int
    eps: int
    eps1: int
    eps2: int
    a: float
    e: float
    x0: float
    x1: float
    e_p: float = 0.0
    eps_p: int = 1

    # signed eccentricity e * eps2 (Henon's (cos eta - sigma cos tau) / D)
    @property
    def signed_e(self) -> float:
        return self.eps2 * self.e

    def _need_circular(self, what: str) -> None:
        if self.e_p != 0.0:
            raise NotImplementedError(
                f"{what} is defined here only for circular primaries (e_p = 0): "
                "Gomez & Olle 1986 give no Jacobi integral or collision velocity."
            )

    def _need_sense(self, what: str) -> None:
        if self.eps1 == 0 and self.e < 1.0 - 1e-12:
            raise ValueError(
                f"{what} needs eps1 (direct or retrograde); this arc has eps1 = 0 "
                "(tangent ellipse: pass eps1 to arc_elements)."
            )

    @property
    def jacobi(self) -> float:
        """C = 2 eps1 sqrt(a (1 - e^2)) + 1/a (Henon eq. 33)."""
        self._need_circular("the Jacobi constant")
        self._need_sense("the Jacobi constant")
        c = self.eps1 * math.sqrt(max(self.a * (1.0 - self.e * self.e), 0.0))
        return 2.0 * c + 1.0 / self.a

    @property
    def speed(self) -> float:
        """V = sqrt(3 - C), speed relative to P2 at collision (Henon eq. 14)."""
        return math.sqrt(max(3.0 - self.jacobi, 0.0))

    @property
    def v1(self) -> float:
        """Radial component at collision, ``e eps2 sqrt(a) sin(eta)`` (Bruno eq. 3)."""
        self._need_circular("the collision velocity")
        return self.signed_e * math.sqrt(self.a) * math.sin(self.eta)

    @property
    def v2(self) -> float:
        """Tangential component, ``eps1 sqrt(a (1 - e^2)) - 1`` (Bruno eq. 3)."""
        self._need_circular("the collision velocity")
        self._need_sense("the collision velocity")
        return self.eps1 * math.sqrt(max(self.a * (1.0 - self.e * self.e), 0.0)) - 1.0

    @property
    def periapsis_distance_primary(self) -> float:
        """a (1 - e), the closest approach to P1."""
        return self.a * (1.0 - self.e)


def classify_signs(tau: float, eta: float, sigma: int) -> tuple[int, int]:
    """(eps, eps2) from Henon's eq. 29: eps2 is the sign of ``cos eta - sigma cos tau``.

    Returns eps2 = 0 where that vanishes (circular orbit, e = 0).
    """
    num = math.cos(eta) - sigma * math.cos(tau)
    eps2 = 1 if num > 0 else (-1 if num < 0 else 0)
    return sigma * eps2, eps2


def arc_elements(
    tau: float,
    eta: float,
    sigma: int,
    eps1: int | None = None,
    *,
    e_p: float = 0.0,
    eps_p: int = 1,
    tol: float = 1e-9,
) -> SArc:
    """Elements of the arc at (tau, eta) with product sign ``sigma`` (eqs. 28-33).

    ``a = r_p D / sin^2(eta)``, ``e = |cos eta - s cos tau| / D`` with
    ``s = eps_p sigma`` and ``D = 1 - s cos tau cos eta``; the abscissae of the two
    axis crossings are ``x0 = sigma a (1 - eps2 e)``, ``x1 = -sigma a (1 + eps2 e)``.
    ``eps1`` is inferred from the second of eqs. 27 when ``sin tau sin eta`` is not
    zero, else it must be given (or is left 0 for a rectilinear arc).

    Only the product of sigma with eps_p matters in D, so for ``e_p = 0`` and
    ``eps_p = -1`` the arc is the ``eps_p = +1`` arc turned by pi (Gomez & Olle).
    """
    sg = eps_p * sigma
    sin_h = math.sin(eta)
    if abs(sin_h) < 1e-14:
        raise ValueError("eta is an integer multiple of pi: use the tangent-ellipse routines")
    cos_t = math.cos(tau)
    cos_h = math.cos(eta)
    d = 1.0 - sg * cos_t * cos_h
    r_p = primary_radius(tau, e_p, eps_p)
    a = r_p * d / (sin_h * sin_h)
    num = cos_h - sg * cos_t
    s_e = num / d
    if abs(s_e) > 1.0 + 1e-9:
        raise ValueError("no elliptic arc: |e| > 1 from eq. 29")
    e = min(abs(s_e), 1.0)
    eps2 = 1 if s_e > tol else (-1 if s_e < -tol else (0 if e < tol else (1 if s_e > 0 else -1)))
    eps = sigma * eps2 if eps2 else 0
    if eps1 is None:
        prod = math.sin(tau) * sin_h
        if abs(prod) > 1e-9 and e < 1.0 - 1e-9:
            eps1 = sigma * (1 if prod > 0 else -1)
            if e_p != 0.0:
                eps1 = sigma * (1 if prod > 0 else -1)
        else:
            eps1 = 0
    x0 = sigma * a * (1.0 - s_e)
    x1 = -sigma * a * (1.0 + s_e)
    return SArc(
        tau=tau,
        eta=eta,
        sigma=sigma,
        eps=eps,
        eps1=eps1,
        eps2=eps2,
        a=a,
        e=e,
        x0=x0,
        x1=x1,
        e_p=e_p,
        eps_p=eps_p,
    )


def tangent_arc(i: int, j: int, eps1: int = 0) -> SArc:
    """The tangent-ellipse arc at (tau/pi, eta/pi) = (i, j) (Henon eqs. 40-41).

    Needs ``sigma = (-1)**(i + j)`` (so both branches of the double point carry
    it), ``a = (i/j)^(2/3)`` and ``e = |1 - (j/i)^(2/3)|``; ``eps1`` is free (the two
    senses are the two branches through the double point) so it must be given for
    the Jacobi constant.  The signs eps, eps2 follow from requiring the collision
    point (E = eta = j pi) to be the tangency point (cos tau, 0), cos tau =
    (-1)^i.  Demanded turn: V1 = 0 identically (type II, a resonance).
    """
    sigma = (-1) ** (i + j)
    a, e = tangent_ellipse(i, j)
    target = (-1.0) ** i
    for eps2 in (1, -1):
        eps = sigma * eps2
        x_t = eps * a * (eps2 * (-1.0) ** j - e)
        if abs(x_t - target) < 1e-9:
            return SArc(
                tau=i * PI,
                eta=j * PI,
                sigma=sigma,
                eps=eps,
                eps1=eps1,
                eps2=eps2,
                a=a,
                e=e,
                x0=sigma * a * (1.0 - eps2 * e),
                x1=-sigma * a * (1.0 + eps2 * e),
            )
    raise ValueError("no sign assignment puts the tangency at the collision point")


def closure_residuals(arc: SArc) -> tuple[float, float, float]:
    """Residuals of Henon's three collision equations (27) for an arc (e_p = 0).

    Returns (cos tau, sin tau, time) residuals; zero for a true arc.  sin tau is
    skipped (returned as 0) when eps1 is undefined.
    """
    arc._need_circular("the closure residuals")
    eps = arc.eps if arc.eps else arc.sigma  # circular arcs: eps eps2 = sigma, eps2 irrelevant
    eps2 = arc.eps2 if arc.eps2 else 1
    a, e, eta = arc.a, arc.e, arc.eta
    r1 = eps * a * (eps2 * math.cos(eta) - e) - math.cos(arc.tau)
    if arc.eps1:
        r2 = eps * arc.eps1 * a * math.sqrt(max(1.0 - e * e, 0.0)) * eps2 * math.sin(
            eta
        ) - math.sin(arc.tau)
    else:
        r2 = 0.0
    r3 = a**1.5 * (eta - eps2 * e * math.sin(eta)) - arc.tau
    return r1, r2, r3


def tau_branches_from_a_e(
    a: float, e: float, eps: int, eps1: int, eps2: int, *, eta_max: float = 7.0 * PI
) -> list[tuple[float, float]]:
    """(eta, tau) pairs consistent with a symmetric arc of elements (a, e) and Brjuno's
    signs (``eps`` = sgn a-tilde is Brjuno's, i.e. Henon's eps times eps2).

    From eq. 2.3, ``cos eta = (1 - 1/a)/(eps2 e)``; eta runs over all branches
    ``+-arccos(.) + 2 pi m`` in (0, eta_max]; tau follows from the third of eqs. 3.9,
    ``tau = a^(3/2)(eta - eps2 e sin eta)``, and a branch is kept only if the first two
    of eqs. 3.9 (``cos tau = eps a (cos eta - eps2 e)`` and ``sin tau = eps eps1
    a sqrt(1 - e^2) sin eta``) close on it to 3e-3.  Needs e > 0, a != 1.
    """
    c0 = (1.0 - 1.0 / a) / (eps2 * e)
    if abs(c0) > 1.0 + 1e-4:  # tangent rows: |cos eta| = 1 up to the 5-digit print
        return []
    base = math.acos(max(-1.0, min(1.0, c0)))
    out: list[tuple[float, float]] = []
    m_max = int(eta_max / TWO_PI) + 1
    for mm in range(m_max + 1):
        for sgn in (1, -1):
            eta = sgn * base + TWO_PI * mm
            if eta <= 1e-9 or eta > eta_max:
                continue
            tau = a**1.5 * (eta - eps2 * e * math.sin(eta))
            r1 = eps * a * (math.cos(eta) - eps2 * e) - math.cos(tau)
            r2 = eps * eps1 * a * math.sqrt(1.0 - e * e) * math.sin(eta) - math.sin(tau)
            near_tangent = (
                abs(math.sin(eta)) < 3e-2
            )  # both senses are solutions at the double point
            if abs(r1) < 3e-3 * max(1.0, a) and (near_tangent or abs(r2) < 3e-3 * max(1.0, a)):
                out.append((eta, tau))
    return sorted(set(out))


# --------------------------------------------------------------------------------------
# Collision velocity, turn, W (Brjuno 1978 section 4.C; Bruno 1981 eqs. 3, 10-12)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class CollisionData:
    """Velocity at collision, demanded turn and indeterminacy flags (#906)."""

    v1: float
    v2: float
    speed: float  # V = sqrt(3 - C)
    jacobi: float
    turn_deg: float  # delta, sin(delta/2) = |V1| / V
    w_eq12: float  # Bruno 1981 eq. 12 (synodic radial speed)
    w_sidereal: float  # Bruno 1981 Table IV form (sidereal speed in the denominator)
    indeterminate: tuple[str, ...]

    @property
    def turn_indeterminate(self) -> bool:
        return bool(self.indeterminate)


def collision_data(arc: SArc, tol: float = _TOL_ZERO) -> CollisionData:
    """Collision velocities, demanded turn and the two forms of Bruno's W.

    Demanded turn: the in and out velocities are (V1, V2) and (-V1, V2), so
    ``cos(delta) = (V2^2 - V1^2)/V^2`` i.e. ``sin(delta/2) = |V1|/V``.

    Indeterminate flags (strings; the turn is still returned):
    ``"tangent_resonance"`` (V1 = 0 because sin eta = 0: type II tangent ellipse,
    integer eta/pi: demanded turn exactly 0), ``"circular"`` (e = 0: V1 = 0),
    ``"half_turn"`` (V2 = 0, demanded turn 180 degrees), ``"radial"`` (e = 1: the
    arc passes through P1, angular momentum zero, sense undefined),
    ``"zero_relative_speed"`` (V = 0, type IV*).

    W forms (both labelled; they coincide only at e = 1):

    * ``w_eq12 = (1/V^2)(V/sqrt(2 - 1/a + a e^2 - a) - 1) eps2 sgn(sin eta)``
      (Bruno 1981 eq. 12, built on the synodic radial speed |V1|);
    * ``w_sidereal = (1/V^2)(V/sqrt(2 - 1/a) - 1)`` (the form that reproduces the
      Table IV W column; sidereal speed in the denominator).
    """
    v1 = arc.v1
    v2 = arc.v2
    cj = arc.jacobi
    v = math.sqrt(max(3.0 - cj, 0.0))
    flags: list[str] = []
    if abs(math.sin(arc.eta)) < tol:
        flags.append("tangent_resonance")
    if arc.e < tol:
        flags.append("circular")
    if abs(arc.e - 1.0) < tol:
        flags.append("radial")
    if abs(v2) < tol:
        flags.append("half_turn")
    if v < tol:
        flags.append("zero_relative_speed")
    ratio = min(abs(v1) / v, 1.0) if v > 0 else 0.0
    turn = math.degrees(2.0 * math.asin(ratio))
    a, e = arc.a, arc.e
    den12 = 2.0 - 1.0 / a + a * e * e - a
    sg = (arc.eps2 or 1) * (1 if math.sin(arc.eta) > 0 else -1)
    w12 = _w_form(v, den12) * sg
    w_sid = _w_form(v, 2.0 - 1.0 / a)
    return CollisionData(v1, v2, v, cj, turn, w12, w_sid, tuple(flags))


def _w_form(v: float, den: float) -> float:
    if v <= 0.0:
        return math.nan
    if den <= 0.0:
        return math.inf
    return (v / math.sqrt(den) - 1.0) / (v * v)


def w_e1(a: float) -> tuple[float, float]:
    """(|W|, V) at e = 1 (Bruno 1981 Table I): ``V = sqrt(3 - 1/a)`` and
    ``|W| = (sqrt(3 - 1/a)/sqrt(2 - 1/a) - 1)/(3 - 1/a)``; infinity at a = 1/2."""
    v2 = 3.0 - 1.0 / a
    v = math.sqrt(v2)
    if 2.0 - 1.0 / a <= 0.0:
        return math.inf, v
    return (v / math.sqrt(2.0 - 1.0 / a) - 1.0) / v2, v


# --------------------------------------------------------------------------------------
# Root finding on the timing equation
# --------------------------------------------------------------------------------------
def find_etas(
    tau: float,
    sigma: int,
    *,
    eta_max: float,
    eta_min: float = 0.0,
    e_p: float = 0.0,
    eps_p: int = 1,
    grid_per_pi: int = 600,
    drop_coincident: bool = True,
) -> list[float]:
    """All roots eta in (eta_min, eta_max) of the timing equation at fixed tau.

    A fine grid brackets sign changes; each is polished with Brent's method.  Roots
    with |e| > 1 (eq. 29) are discarded, as is the trivial coincident orbit
    (eta = tau, a = 1, e = 0; sigma = +1 only) unless ``drop_coincident`` is False.
    Two roots closer than the grid spacing (near a turning point of the curve in
    the (tau, eta) plane) may be missed; refine with ``refine_point``.
    """
    n = max(int((eta_max - eta_min) / PI * grid_per_pi), 20)
    etas = np.linspace(eta_min + 1e-9, eta_max, n)
    f = _timing_residual_vec(tau, etas, sigma, e_p, eps_p)
    roots: list[float] = []
    sg = eps_p * sigma
    for i in range(n - 1):
        f0, f1 = float(f[i]), float(f[i + 1])
        if f0 == 0.0:
            cand = float(etas[i])
        elif f0 * f1 < 0.0:
            cand = brentq(
                lambda x: timing_residual(tau, x, sigma, e_p, eps_p),
                float(etas[i]),
                float(etas[i + 1]),
                xtol=1e-15,
                rtol=4 * np.finfo(float).eps,
            )
        else:
            continue
        d = 1.0 - sg * math.cos(tau) * math.cos(cand)
        if d <= 0.0 or abs(math.cos(cand) - sg * math.cos(tau)) > d * (1.0 + 1e-9):
            continue
        if drop_coincident and sigma == 1 and abs(cand - tau) < 1e-7:
            continue
        if roots and abs(cand - roots[-1]) < 1e-12:
            continue
        roots.append(cand)
    return roots


def refine_point(
    tau0: float, eta0: float, sigma: int, *, e_p: float = 0.0, eps_p: int = 1
) -> tuple[float, float]:
    """Move (tau0, eta0) onto the curve F = 0 along the better-conditioned axis.

    The root of eq. 30 in eta is ill-conditioned where the curve is horizontal in
    the (tau, eta) plane (and the root in tau where it is vertical); the gradient
    decides which variable to solve for.
    """
    h = 1e-6
    ft = (
        timing_residual(tau0 + h, eta0, sigma, e_p, eps_p)
        - timing_residual(tau0 - h, eta0, sigma, e_p, eps_p)
    ) / (2 * h)
    fe = (
        timing_residual(tau0, eta0 + h, sigma, e_p, eps_p)
        - timing_residual(tau0, eta0 - h, sigma, e_p, eps_p)
    ) / (2 * h)
    if abs(fe) >= abs(ft):
        eta = _bracketed_root(lambda x: timing_residual(tau0, x, sigma, e_p, eps_p), eta0, 0.05)
        return tau0, eta
    tau = _bracketed_root(lambda x: timing_residual(x, eta0, sigma, e_p, eps_p), tau0, 0.05)
    return tau, eta0


def _bracketed_root(fun: Callable[[float], float], x0: float, width: float) -> float:
    """Brent root of ``fun`` bracketing x0 +- growing width."""
    w = width
    f0 = fun(x0)
    if f0 == 0.0:
        return x0
    for _ in range(12):
        lo, hi = x0 - w, x0 + w
        flo, fhi = fun(lo), fun(hi)
        if flo * fhi < 0.0:
            return float(brentq(fun, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps))
        # try sub-intervals around x0
        if flo * f0 < 0.0:
            return float(brentq(fun, lo, x0, xtol=1e-15, rtol=4 * np.finfo(float).eps))
        if f0 * fhi < 0.0:
            return float(brentq(fun, x0, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps))
        w *= 0.5
    raise ValueError("no bracketed root near the starting point")


def enumerate_arcs(
    tau: float,
    *,
    eta_max: float,
    sigmas: Sequence[int] = (-1, 1),
    e_p: float = 0.0,
    eps_p: int = 1,
    drop_coincident: bool = True,
) -> list[SArc]:
    """All symmetric arcs at duration parameter ``tau`` (both sigma by default).

    Arcs whose sense eps1 is not determined are returned with ``eps1 = 0``.
    """
    out: list[SArc] = []
    for sg in sigmas:
        for eta in find_etas(
            tau, sg, eta_max=eta_max, e_p=e_p, eps_p=eps_p, drop_coincident=drop_coincident
        ):
            try:
                out.append(arc_elements(tau, eta, sg, e_p=e_p, eps_p=eps_p))
            except ValueError:
                continue
    return out


def e1_arcs(k: int, eps2: int, *, eta_max: float = 6.6 * PI) -> list[SArc]:
    """The e = 1 (radial, through P1) arcs.

    Henon eq. 29 gives ``e = 1`` iff ``(eps2 cos eta - 1)(1 + eps cos tau) = 0``;
    the second factor means ``tau = k pi`` with ``eps = (-1)**(k+1)`` (the first,
    sin eta = 0, is the tangent-ellipse family).  Then ``a = 1/(1 - eps2 cos eta)``
    and the timing equation fixes eta.  Returns the arcs for tau = k pi with
    eps2 = +-1 (so sigma = eps eps2), eta in (0, eta_max).
    """
    if eps2 not in (-1, 1):
        raise ValueError("eps2 must be +1 or -1")
    eps = (-1) ** (k + 1)
    sigma = eps * eps2
    tau = k * PI
    arcs = []
    for eta in find_etas(tau, sigma, eta_max=eta_max):
        s = abs(math.sin(eta))
        if s < 1e-6:
            continue
        a = 1.0 / (1.0 - eps2 * math.cos(eta))
        arc = SArc(
            tau=tau,
            eta=eta,
            sigma=sigma,
            eps=eps,
            eps1=0,
            eps2=eps2,
            a=a,
            e=1.0,
            x0=sigma * a * (1.0 - eps2),
            x1=-sigma * a * (1.0 + eps2),
        )
        arcs.append(arc)
    return arcs


# --------------------------------------------------------------------------------------
# Hyperbolic and parabolic arcs (Henon eqs. 5-22; Gomez & Olle eqs. 19, 20)
# --------------------------------------------------------------------------------------
def hyperbolic_residual(tau: float, eta: float) -> float:
    """Henon eq. 10 for the (eps = -1, eps1 = -1) hyperbolic family A0."""
    ct = math.cos(tau)
    d = 1.0 + ct * math.cosh(eta)
    return (
        math.sqrt(d) * (math.sinh(eta) * (math.cosh(eta) + ct) - eta * d)
        - tau * math.sinh(eta) ** 3
    )


@dataclass(frozen=True)
class HyperbolicArc:
    tau: float
    eta: float
    a: float
    e: float
    x0: float
    jacobi: float
    speed: float


def hyperbolic_arc(tau: float, eta: float | None = None) -> HyperbolicArc:
    """Hyperbolic arc of A0 at ``tau`` (eqs. 8, 12, 13, 14 with eps = eps1 = -1).

    If ``eta`` is None it is found from eq. 10 (unique for 0 < tau < tau_parabolic).
    """
    if eta is None:
        hi = 1.0
        while hyperbolic_residual(tau, hi) * hyperbolic_residual(tau, 1e-3) > 0 and hi < 60:
            hi *= 1.5
        eta = float(brentq(lambda x: hyperbolic_residual(tau, x), 1e-3, hi, xtol=1e-15))
    ct = math.cos(tau)
    a = (1.0 + ct * math.cosh(eta)) / math.sinh(eta) ** 2
    e = (math.cosh(eta) + ct) / (1.0 + ct * math.cosh(eta))
    x0 = -a * (e - 1.0)
    cj = -2.0 * math.sqrt(a * (e * e - 1.0)) - 1.0 / a
    return HyperbolicArc(tau, eta, a, e, x0, cj, math.sqrt(3.0 - cj))


def parabolic_tau(e_p: float = 0.0, eps_p: int = 1, eps: int | None = None) -> float:
    """Duration parameter tau of the unique parabolic arc (Henon eq. 20; Gomez &
    Olle eq. 19).

    ``eps`` defaults to the value that has a solution: -1 for eps_p = +1, +1 for
    eps_p = -1.  Solves ``M_p(tau) = (1/3) r_p^(3/2) sqrt(1 - s cos tau)(2 + s cos tau)``
    with ``s = eps_p eps``; at e_p = 0, eps_p = 1, eps = -1 this is ``tau =
    (1/3)(2 - cos tau) sqrt(1 + cos tau)`` and ``tau/pi = 0.163926``.
    """
    if eps is None:
        eps = -eps_p
    s = eps_p * eps

    def f(t: float) -> float:
        ct = math.cos(t)
        r_p = primary_radius(t, e_p, eps_p)
        return float(
            primary_mean_anomaly(t, e_p, eps_p)
            - r_p**1.5 * math.sqrt(1.0 - s * ct) * (2.0 + s * ct) / 3.0
        )

    xs = np.linspace(1e-6, PI, 4000)
    vals = [f(float(x)) for x in xs]
    for i in range(len(xs) - 1):
        if vals[i] * vals[i + 1] < 0.0:
            return float(brentq(f, float(xs[i]), float(xs[i + 1]), xtol=1e-15))
    raise ValueError("no parabolic arc for these signs")


@dataclass(frozen=True)
class ParabolicArc:
    tau: float
    p: float
    x0: float
    jacobi: float
    speed: float
    sigma_param: float


def parabolic_arc() -> ParabolicArc:
    """The circular-problem parabolic arc (Henon eqs. 16-22)."""
    tau = parabolic_tau(0.0, 1, -1)
    p = 1.0 - math.cos(tau)
    sig = math.sqrt((1.0 + math.cos(tau)) / (1.0 - math.cos(tau)))
    # eps = -1: p = 1 + eps cos tau; sigma^2 = (1 - eps cos tau)/(1 + eps cos tau)
    return ParabolicArc(
        tau, p, -p / 2.0, -2.0 * math.sqrt(p), math.sqrt(3.0 + 2.0 * math.sqrt(p)), sig
    )


# --------------------------------------------------------------------------------------
# Families and tangent ellipses (Henon eqs. 36, 40-43; Gomez & Olle eq. 21)
# --------------------------------------------------------------------------------------
def family_sigma(kind: str, i: int = 0, j: int = 0) -> int:
    """sigma of a family: A -> -1, B -> +1, C_ij -> (-1)**(i + j)."""
    if kind == "A":
        return -1
    if kind == "B":
        return 1
    if kind == "C":
        return int((-1) ** (i + j))
    raise ValueError("kind must be 'A', 'B' or 'C'")


def tangent_ellipse(i: int, j: int, e_p: float = 0.0, eps_p: int = 1) -> tuple[float, float]:
    """(a, e) of the ellipse tangent to the unit circle at (tau/pi, eta/pi) = (i, j).

    ``a = (i/j)^(2/3)``, ``e = |1 - (j/i)^(2/3) (1 + (-1)^(i+1) eps_p e_p)|``
    (Gomez & Olle eq. 20-21 text; Henon eq. 41 at e_p = 0).  Raises ValueError if
    e > 1 (the point does not exist).
    """
    a = (i / j) ** (2.0 / 3.0)
    e = abs(1.0 - (j / i) ** (2.0 / 3.0) * (1.0 + (-1) ** (i + 1) * eps_p * e_p))
    if e > 1.0 + 1e-12:
        raise ValueError("no tangent ellipse: e > 1")
    return a, e


def c_family_max_j(i: int, e_p: float = 0.0, eps_p: int = 1) -> float:
    """Upper bound on j for the C_ij family (Gomez & Olle eq. 21):
    ``j <= (2/(1 + (-1)^(i+1) eps_p e_p))^(3/2) i``; at e_p = 0 this is Henon's
    ``j/i <= 2 sqrt 2``."""
    return float((2.0 / (1.0 + (-1) ** (i + 1) * eps_p * e_p)) ** 1.5 * i)


def c_family_exists(i: int, j: int, e_p: float = 0.0, eps_p: int = 1) -> bool:
    """Whether C_ij exists: j > i and j <= ``c_family_max_j`` (eq. 21)."""
    return j > i and j <= c_family_max_j(i, e_p, eps_p) * (1.0 + 1e-12)


def c_family_indices(
    max_i: int, e_p: float = 0.0, eps_p: int = 1, *, max_j: int = 200
) -> list[tuple[int, int]]:
    """All (i, j) with 1 <= i <= max_i for which C_ij exists (eq. 21)."""
    out = []
    for i in range(1, max_i + 1):
        for j in range(i + 1, max_j + 1):
            if c_family_exists(i, j, e_p, eps_p):
                out.append((i, j))
    return out


def first_c_family(
    e_p: float,
    eps_p: int,
    i_parity: int,
    *,
    max_i: int = 2000,
    factor: float | None = None,
) -> tuple[int, int]:
    """Smallest (i, i + 1) with C_ij existing and i of the given parity (0 even, 1 odd).

    ``factor`` replaces the exact ratio ``c_family_max_j(i)/i`` of eq. 21 (use it to
    reproduce the rounded 1.015 of the paper's text).  Gomez & Olle p. 44 quote C67,68
    (eps_p = +1, i odd) and C68,69 (eps_p = -1, i even) at e_p = 0.98; with the exact
    eq. 21 the second is C66,67 (66 * 1.01523 = 67.005 >= 67), the printed pair following
    from the rounded 1.015 (66 * 1.015 = 66.99).
    """
    for i in range(1, max_i + 1):
        if i % 2 != i_parity % 2:
            continue
        j = i + 1
        ok = c_family_exists(i, j, e_p, eps_p) if factor is None else j <= factor * i
        if ok:
            return i, j
    raise ValueError("none found")


# --------------------------------------------------------------------------------------
# Junction closed forms (Brjuno 1978b section 3; Henon eq. 41)
# --------------------------------------------------------------------------------------
def type_ii_junction(p: int, q: int, eps1: int) -> dict[str, float]:
    """Type II (tangent-ellipse) junction at a = (p/(p+q))^(2/3): e = |a-1|/a,
    C = 1/a + 2 eps1 sqrt(2 - 1/a), V1 = 0 so the demanded turn is exactly 0."""
    a = (p / (p + q)) ** (2.0 / 3.0)
    e = abs(a - 1.0) / a
    cj = 1.0 / a + 2.0 * eps1 * math.sqrt(2.0 - 1.0 / a)
    return {"a": a, "e": e, "C": cj, "V": math.sqrt(3.0 - cj), "turn_deg": 0.0}


JUNCTION_CONSTANTS: dict[str, float] = {"type_III_C": -1.0, "type_IV_star_C": 3.0}


# --------------------------------------------------------------------------------------
# Asymmetric arcs T_N (Brjuno 1978 section 3.C and 4.C)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class TArc:
    """An asymmetric arc: a resonant ellipse of mean motion N = (p+q)/p with one
    collision deleted.  Departure and arrival velocities are equal so the demanded
    turn at mu = 0 is exactly zero."""

    p: int
    q: int
    a: float
    e: float
    eps1: int
    eps2: int
    eta: float
    v1: float
    v2: float
    jacobi: float
    speed: float
    duration: float  # 2 pi p, time between successive collisions at the same point
    turn_deg: float = 0.0


def t_n_exists(p: int, q: int) -> bool:
    """T_N exists for rational N = (p+q)/p < 2 sqrt 2 (Brjuno II p.32; N = 3 empty)."""
    return (p + q) / p < 2.0 * math.sqrt(2.0)


def t_n_eta_intervals(p: int, q: int, k_max: int) -> list[tuple[int, float, float]]:
    """Intervals J_k = [k pi - arccos|1 - 1/a|, k pi + arccos|1 - 1/a|] (eq. 3.8) of
    admissible eta, for k = 0 .. k_max."""
    a = (p / (p + q)) ** (2.0 / 3.0)
    w = math.acos(min(1.0, abs(1.0 - 1.0 / a)))
    return [(k, k * PI - w, k * PI + w) for k in range(k_max + 1)]


def t_n_arc(p: int, q: int, eps1: int, eta: float) -> TArc:
    """T_N arc at eccentric anomaly eta (eq. 3.8): ``eps2 e = (1 - 1/a)/cos eta``.

    ``a = (p/(p+q))^(2/3)``.  Raises ValueError if |cos eta| < |1 - 1/a| (eta
    outside every J_k) or the bound N < 2 sqrt 2 fails.
    """
    if not t_n_exists(p, q):
        raise ValueError("T_N does not exist for N >= 2 sqrt 2")
    if eps1 not in (-1, 1):
        raise ValueError("eps1 must be +-1")
    a = (p / (p + q)) ** (2.0 / 3.0)
    ce = math.cos(eta)
    if abs(ce) < 1e-14:
        raise ValueError("cos eta = 0 requires a = 1")
    s_e = (1.0 - 1.0 / a) / ce
    if abs(s_e) > 1.0 + 1e-12:
        raise ValueError("eta outside the admissible intervals J_k")
    e = min(abs(s_e), 1.0)
    eps2 = 1 if s_e >= 0 else -1
    v1 = s_e * math.sqrt(a) * math.sin(eta)
    c = eps1 * math.sqrt(max(a * (1.0 - e * e), 0.0))
    v2 = c - 1.0
    cj = 2.0 * c + 1.0 / a
    return TArc(p, q, a, e, eps1, eps2, eta, v1, v2, cj, math.sqrt(max(3.0 - cj, 0.0)), TWO_PI * p)


# --------------------------------------------------------------------------------------
# e* membership test (Brjuno 1978b Theorem 2.2; sign reconstruction of eq. 2.12)
# --------------------------------------------------------------------------------------
def e_star_abs(a: float, e: float, eps1: int) -> float:
    """|e*| = |Re sqrt(Q1)| from (a, e, eps1), per the reconstruction of eq. 2.12
    in the Brjuno 1978b digest (section 1.1).  Needs ``|a - 1| < a e`` (a collision
    is possible) and ``0 < e < 1``; the overall sign of e* is not returned.

    For a > 1: ``arg Q1 = arg(X + iY) + N^-1 arg(U - i eps1 V) + eps1 w``;
    for a < 1: ``arg Q1 = arg(X - i eps1 Y) + N^-1 arg(U + iV) - w + pi|N^-1 - 1|``;
    with ``X = (a - 1 - a e^2)/e``, ``Y = sqrt(1-e^2) sqrt(a^2 e^2 - (a-1)^2)/e``,
    ``U = (a - 1)/(a e)``, ``V = sqrt(a^2 e^2 - (a-1)^2)/(a e)``, ``w = sqrt(a)
    sqrt(a^2 e^2 - (a-1)^2)``, ``N^-1 = a^(3/2)``.
    """
    if eps1 not in (-1, 1):
        raise ValueError("eps1 must be +-1")
    rad = a * a * e * e - (a - 1.0) ** 2
    if rad < -1e-12:
        raise ValueError("no collision possible: |a - 1| > a e")
    rad = max(rad, 0.0)
    ninv = a**1.5
    x = (a - 1.0 - a * e * e) / e
    y = math.sqrt(1.0 - e * e) * math.sqrt(rad) / e
    u = (a - 1.0) / (a * e)
    v = math.sqrt(rad) / (a * e)
    w = math.sqrt(a) * math.sqrt(rad)
    if a > 1.0:
        arg = math.atan2(y, x) + ninv * math.atan2(-eps1 * v, u) + eps1 * w
    else:
        arg = math.atan2(-eps1 * y, x) + ninv * math.atan2(v, u) - w + PI * abs(ninv - 1.0)
    return abs(math.cos(arg / 2.0))


def theorem_2_2_values(kind: str, ninv: float, j: int = 0, k: int = 0) -> list[float]:
    """Candidate |e*| from Theorem 2.2 for a family and N^-1 = a^(3/2).

    * A_j: omega_3 ``|sin((2[j/2] + 1)(pi/2)(N^-1 - 1))|`` and omega_4
      ``|sin([(j+1)/2] pi (N^-1 - 1))|`` (both returned);
    * B_k: ``|cos(k (pi/2)(N^-1 - 1))|``;
    * C_jk: ``|cos(k (pi/2)(N^-1 - 1))|`` (omega_1, omega_2) and ``|sin(...)|``
      (omega_3, omega_4) (both returned).

    The domain omega_i of an arc is not computed here, so every admissible form is
    listed; membership means matching any of them.
    """
    d = ninv - 1.0
    if kind == "A":
        return [
            abs(math.sin((2 * (j // 2) + 1) * (PI / 2.0) * d)),
            abs(math.sin(((j + 1) // 2) * PI * d)),
        ]
    if kind == "B":
        return [abs(math.cos(k * (PI / 2.0) * d))]
    if kind == "C":
        return [abs(math.cos(k * (PI / 2.0) * d)), abs(math.sin(k * (PI / 2.0) * d))]
    raise ValueError("kind must be 'A', 'B' or 'C'")


def e_star_matches(
    a: float, e: float, eps1: int, kind: str, j: int = 0, k: int = 0, tol: float = 3e-3
) -> bool:
    """Whether (a, e, eps1) is on a characteristic of the named family (the
    ``e*`` membership test, independent of the timing equation)."""
    es = e_star_abs(a, e, eps1)
    return any(abs(es - c) < tol for c in theorem_2_2_values(kind, a**1.5, j, k))


# --------------------------------------------------------------------------------------
# The curve f of Brjuno 1978b Table II (P = 0, eps1 = -1, a < 1)
# --------------------------------------------------------------------------------------
def curve_f_p(a: float, e: float) -> float:
    """P(a, e) = a e^2 + a - 1 - a^(3/2) sqrt(1 - e^2)(1 - a + a e^2) (eq. 2.16'), zero on f."""
    return float(a * e * e + a - 1.0 - a**1.5 * math.sqrt(1.0 - e * e) * (1.0 - a + a * e * e))


def curve_f_from_c(c: float) -> tuple[float, float]:
    """(a, e) on f from the area integral c in [0, 1): ``a = (c^2 + 1)/(c^3 - c + 2)``,
    ``1 - e^2 = c^2/a`` (parametric form, section 2.E)."""
    a = (c * c + 1.0) / (c**3 - c + 2.0)
    return a, math.sqrt(1.0 - c * c / a)


Z0 = 1.0 / (2.0 * math.sqrt(2.0) - 1.0)


def curve_f_phi(a: float, e_star: float) -> tuple[float, float, float]:
    """(x, y, phi) on f: ``x = 1/(N - 1)`` with ``N = a^(-3/2)``, ``y = 2 arccos(e*)/
    (pi (1 - N^-1))`` and ``phi = y - (x - z0)`` (eq. 2.30 and Table II)."""
    ninv = a**1.5
    nn = 1.0 / ninv
    x = 1.0 / (nn - 1.0)
    y = 2.0 * math.acos(e_star) / (PI * (1.0 - ninv))
    return x, y, y - (x - Z0)


# --------------------------------------------------------------------------------------
# Critical-orbit function (Hitzl & Henon 1977b eqs. 40-42)
# --------------------------------------------------------------------------------------
def critical_function_s(tau: float, eta: float, sigma: int) -> float:
    """S of Hitzl & Henon 1977b eq. 40; zero at critical orbits (extrema of the
    Jacobi constant along a family).  ``G* = -(1/2) sigma rho sin^2(eta) S`` (eq. 42)
    with ``rho = sgn(sin eta) sqrt(1 - sigma cos tau cos eta)``.
    """
    ct, ch = math.cos(tau), math.cos(eta)
    st, sh = math.sin(tau), math.sin(eta)
    s2t = math.sin(2.0 * tau)
    d = 1.0 - sigma * ct * ch
    sgn_h = 1.0 if sh > 0 else -1.0
    t1 = 2.0 * s2t
    t2 = -6.0 * (ch - sigma * ct) ** 2 * tau / (sh * sh)
    bracket = ch * (2.0 * ch**2 + 3.0 * ct**2 - ct**4) - sigma * ct * (
        2.0 + 2.0 * ch**2 - ct**2 + ch**4
    )
    t3 = 2.0 * sgn_h * bracket / (sh * sh * math.sqrt(d))
    _ = st
    return t1 + t2 + t3


# --------------------------------------------------------------------------------------
# R-arcs and R-orbits (Henon 2001 chapter 18; Devaney 1981)
# --------------------------------------------------------------------------------------
def _damped_newton(
    y0: NDArray[np.float64],
    grad_hess: Callable[[NDArray[np.float64]], tuple[NDArray[np.float64], NDArray[np.float64]]],
    phi: Callable[[NDArray[np.float64]], float],
    signs: NDArray[np.float64],
    tol: float,
    max_iter: int = 5000,
) -> NDArray[np.float64]:
    y = y0.copy()
    for _ in range(max_iter):
        g, h = grad_hess(y)
        scale = max(1.0, float(np.max(np.abs(1.0 / y))))
        gmax = float(np.max(np.abs(g)))
        if gmax < tol * scale:
            return y
        step = np.linalg.solve(h, -g)
        if gmax < 1e-8 * scale:
            # quadratic-convergence regime: the Armijo test is below roundoff here
            yn = y + step
            if np.all(np.sign(yn) == signs):
                y = yn
                continue
        t = 1.0
        f0 = phi(y)
        slack = 8.0 * np.finfo(float).eps * (abs(f0) + 1.0)
        while True:
            yn = y + t * step
            if np.all(np.sign(yn) == signs) and phi(yn) <= f0 + 1e-4 * t * float(g @ step) + slack:
                break
            t *= 0.5
            if t < 1e-16:
                return y
        y = yn
    raise RuntimeError("damped Newton did not converge")


def r_arc(signs: Sequence[int], tol: float = 1e-12) -> NDArray[np.float64]:
    """The R-arc with y_0 = y_n = 0 and sign(y_i) = signs[i-1] (i = 1 .. n-1).

    Solves ``y[i-1] - 2 y[i] + y[i+1] + 1/y[i] = 0`` (Henon 2001 eq. 18.4).  The
    unique critical point of the convex function ``Phi`` in the orthant; returns
    the vector (y_1, ..., y_(n-1)).
    """
    s = np.array(signs, dtype=float)
    m = len(s)
    if m == 0 or not np.all(np.abs(s) == 1.0):
        raise ValueError("signs must be a non-empty sequence of +-1")
    lap = 2.0 * np.eye(m) - np.eye(m, k=1) - np.eye(m, k=-1)

    def phi(y: NDArray[np.float64]) -> float:
        # R-arc: sum_{j=0}^{n-1} (y_{j+1}-y_j)^2 / 2 with y_0 = y_n = 0
        return float(0.5 * y @ lap @ y - np.sum(np.log(np.abs(y))))

    def gh(y: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        g = lap @ y - 1.0 / y
        h = lap + np.diag(1.0 / (y * y))
        return g, h

    return _damped_newton(s.copy(), gh, phi, s, tol)


def r_arcs(n: int) -> dict[tuple[int, ...], NDArray[np.float64]]:
    """All 2^(n-1) R-arcs of order n >= 2, keyed by the sign code of (y_1 .. y_(n-1)).

    Hard assertion of the Devaney count: exactly ``2**(n-1)`` arcs, all distinct.
    """
    out: dict[tuple[int, ...], NDArray[np.float64]] = {}
    for code in itertools.product((1, -1), repeat=n - 1):
        out[code] = r_arc(code)
    assert len(out) == 2 ** (n - 1)
    return out


def r_orbit(signs: Sequence[int], tol: float = 1e-12) -> NDArray[np.float64]:
    """The R-orbit (period n, cyclic) with sign(y_i) = signs[i]; unique for every
    sign code except all + and all - (which have no finite solution)."""
    s = np.array(signs, dtype=float)
    n = len(s)
    if n < 2 or not np.all(np.abs(s) == 1.0):
        raise ValueError("signs must be a sequence of +-1 with length >= 2")
    if np.all(s == 1.0) or np.all(s == -1.0):
        raise ValueError("all-equal sign codes have no R-orbit (y diverges)")
    lap = 2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    lap[0, n - 1] -= 1.0
    lap[n - 1, 0] -= 1.0

    def phi(y: NDArray[np.float64]) -> float:
        return float(0.5 * y @ lap @ y - np.sum(np.log(np.abs(y))))

    def gh(y: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        return lap @ y - 1.0 / y, lap + np.diag(1.0 / (y * y))

    return _damped_newton(s.copy(), gh, phi, s, tol)


def r_orbits(n: int) -> dict[tuple[int, ...], NDArray[np.float64]]:
    """All 2^n - 2 R-orbits of order n, keyed by sign code (sub-period solutions
    included, as in Henon's count).  Hard assertion of the count."""
    out: dict[tuple[int, ...], NDArray[np.float64]] = {}
    for code in itertools.product((1, -1), repeat=n):
        if len(set(code)) == 1:
            continue
        out[code] = r_orbit(code)
    assert len(out) == 2**n - 2
    return out


def r_orbit_residual(y: NDArray[np.float64]) -> float:
    """Max |y[i-1] - 2 y[i] + y[i+1] + 1/y[i]| for a cyclic R-orbit."""
    r = np.roll(y, 1) - 2.0 * y + np.roll(y, -1) + 1.0 / y
    return float(np.max(np.abs(r)))


def r_arc_residual(y: NDArray[np.float64]) -> float:
    """Max residual of the recurrence for an R-arc vector (y_1 .. y_(n-1))."""
    full = np.concatenate(([0.0], y, [0.0]))
    r = full[:-2] - 2.0 * full[1:-1] + full[2:] + 1.0 / full[1:-1]
    return float(np.max(np.abs(r)))


def refine_r_orbit_mp(y: Sequence[float], dps: int = 50) -> list[str]:
    """Polish an R-orbit at ``dps`` digits by Newton in extended precision
    (mpmath); returns the decimal strings.  Not needed for correctness at moderate
    n (see the module docstring) but gives an independent high-precision check."""
    import mpmath as mp

    mp.mp.dps = dps
    n = len(y)
    yy = mp.matrix([mp.mpf(v) for v in y])
    for _ in range(8):
        f = mp.matrix(n, 1)
        jac = mp.matrix(n, n)
        for i in range(n):
            im, ip = (i - 1) % n, (i + 1) % n
            f[i] = yy[im] - 2 * yy[i] + yy[ip] + 1 / yy[i]
            jac[i, im] += 1
            jac[i, ip] += 1
            jac[i, i] += -2 - 1 / yy[i] ** 2
        yy = yy - mp.lu_solve(jac, f)
    return [mp.nstr(v, dps) for v in yy]
