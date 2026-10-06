"""Two-body Kepler propagator (universal-variable formulation).

Propagates an inertial state ``(r0, v0)`` by an arbitrary signed time interval
``dt`` under a central body with gravitational parameter ``mu``. The universal
formulation handles elliptic, parabolic, and hyperbolic regimes with one set of
equations, and supports backward propagation (``dt < 0``) natively.

References
----------
Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th ed.,
Microcosm Press, 2013, Algorithm 3.4 (KEPLER) and §2.2 (universal variables);
Bate, R. R., Mueller, D. D., and White, J. E., *Fundamentals of Astrodynamics*,
Dover, 1971, §4.4.

Tolerances
----------
* Newton convergence: relative chi-step ``|delta| / max(|chi|, 1) < 1e-12``.
  (Absolute residual ``|f|`` scales with ``sqrt(mu) * dt`` and is not a useful
  convergence criterion across short and long propagations.)
* Iteration cap: 200 (#963: the bracketed fallback bisects; plain Newton needs far fewer).
* Newton is safeguarded by a bracket on the monotone f(chi) (#963).

JIT acceleration (#475)
-----------------------
The universal-anomaly Newton iteration is compiled to native code via
``numba.njit(cache=True)`` in ``_kepler_chi_newton``.  The public
``propagate()`` function calls the JIT core and applies the Lagrange
coefficients in numpy (array ops are not the bottleneck).  The pure-Python
reference ``_kepler_chi_newton_py`` is retained permanently for parity testing.

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §3.3.
"""

from __future__ import annotations

from math import cos, inf, isfinite, log, nan, sin, sqrt

import numba as nb
import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core._stumpff import stumpff_c, stumpff_c_py, stumpff_s, stumpff_s_py
from cyclerfinder.core.constants import MU_SUN_KM3_S2

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

_NEWTON_TOL_DELTA_REL: float = 1.0e-12
_NEWTON_MAX_ITER: int = 200


class KeplerError(Exception):
    """Base class for errors raised by :func:`propagate`."""


class KeplerConvergenceError(KeplerError):
    """Newton iteration on the universal anomaly failed to converge.

    Attributes
    ----------
    chi:
        Last universal-anomaly iterate.
    residual:
        Last value of ``f(chi)``.
    """

    def __init__(self, chi: float, residual: float) -> None:
        super().__init__(
            f"Kepler universal-variable Newton failed to converge: "
            f"chi={chi:.6e}, residual={residual:.6e}"
        )
        self.chi = chi
        self.residual = residual


def coe_to_rv(
    a_km: float,
    e: float,
    true_anom_rad: float,
    mu: float = MU_SUN_KM3_S2,
    *,
    arg_peri_rad: float = 0.0,
) -> tuple[Vec3, Vec3]:
    """Planar orbital elements -> inertial heliocentric state (inclination 0).

    Converts a classical-element triple ``(a, e, nu)`` into an inertial
    position/velocity pair using the standard perifocal formulation rotated
    by the argument of periapsis about the +z axis. The orbit is treated as
    planar (inclination and longitude of ascending node both zero), so the
    inertial frame coincides with the ecliptic and the state lies in the
    ``z = 0`` plane. This is the circular-coplanar idealisation used by the
    resonance-anchored construction.

    Parameters
    ----------
    a_km:
        Semi-major axis, km. Must be positive (elliptic).
    e:
        Eccentricity, ``0 <= e < 1``.
    true_anom_rad:
        True anomaly ``nu``, radians, measured from periapsis.
    mu:
        Central-body gravitational parameter, km^3/s^2. Defaults to the
        heliocentric :data:`cyclerfinder.core.constants.MU_SUN_KM3_S2`.
    arg_peri_rad:
        Argument of periapsis ``omega``, radians; rotates the perifocal
        frame about +z into the inertial frame. Defaults to ``0.0``.

    Returns
    -------
    (r, v):
        Inertial position and velocity as two ``(3,)`` float64 arrays in
        km and km/s, lying in the ``z = 0`` plane.

    Notes
    -----
    Perifocal frame (Vallado §2.6, eqs. 2-110/2-115; Bate-Mueller-White
    §2.5): with semi-latus rectum ``p = a(1 - e^2)`` and radius
    ``r = p / (1 + e cos nu)``::

        r_pf = r * [cos nu, sin nu, 0]
        v_pf = sqrt(mu / p) * [-sin nu, e + cos nu, 0]

    rotated by ``omega`` about +z to obtain the inertial state.
    """
    p = a_km * (1.0 - e * e)
    nu = true_anom_rad
    cos_nu = cos(nu)
    sin_nu = sin(nu)
    r_mag = p / (1.0 + e * cos_nu)

    r_pf = np.array([r_mag * cos_nu, r_mag * sin_nu, 0.0], dtype=np.float64)
    sqrt_mu_p = sqrt(mu / p)
    v_pf = np.array([-sqrt_mu_p * sin_nu, sqrt_mu_p * (e + cos_nu), 0.0], dtype=np.float64)

    cos_w = cos(arg_peri_rad)
    sin_w = sin(arg_peri_rad)
    rot = np.array(
        [
            [cos_w, -sin_w, 0.0],
            [sin_w, cos_w, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )

    r = rot @ r_pf
    v = rot @ v_pf
    return r, v


# ---------------------------------------------------------------------------
# Universal-anomaly Newton core — pure-Python reference
# ---------------------------------------------------------------------------


def _kepler_chi_newton_py(
    r0_n: float,
    v0_n: float,
    rv_dot: float,
    alpha: float,
    dt: float,
    mu: float,
    chi0: float,
) -> tuple[float, float, float, float, float]:
    """Pure-Python reference for the universal-anomaly Newton solve.

    Retained permanently as oracle for JIT parity tests.

    Returns ``(chi, f, g, f_dot, g_dot)`` — the converged universal anomaly
    and the four Lagrange coefficients that reconstruct ``(r, v)`` from
    ``(r0, v0)``.  Raises on non-convergence by returning sentinel
    ``chi = nan`` (caller must detect and raise :class:`KeplerConvergenceError`).
    """
    sqrt_mu = sqrt(mu)
    # #963: Newton safeguarded by a bracket. f(chi) is strictly increasing (f'(chi) = r > 0), and
    # f(0) = -sqrt(mu) dt, so the root lies in (0, inf) for dt > 0 and (-inf, 0) for dt < 0.
    # Every evaluated iterate narrows the bracket; a Newton step that leaves it is replaced by a
    # bisection (or, while the far side is still open, by a step that doubles the distance from
    # 0). A Newton step that converges is accepted before the bracket test, so iterates that
    # the plain Newton solve accepted are unchanged.
    if dt > 0.0:
        lo = 0.0
        hi = inf
    elif dt < 0.0:
        lo = -inf
        hi = 0.0
    else:
        lo = 0.0
        hi = 0.0
    chi = chi0
    residual = 0.0
    converged = False
    dx_last = inf
    chi_par = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0
    if chi_par == 0.0:
        chi_par = 1.0 if dt > 0.0 else -1.0
    for _iteration in range(_NEWTON_MAX_ITER):
        z = chi * chi * alpha
        c = stumpff_c_py(z)
        s = stumpff_s_py(z)
        chi2 = chi * chi
        chi3 = chi2 * chi

        f_val = (
            (rv_dot / sqrt_mu) * chi2 * c
            + (1.0 - alpha * r0_n) * chi3 * s
            + r0_n * chi
            - sqrt_mu * dt
        )
        f_prime = (rv_dot / sqrt_mu) * chi * (1.0 - z * s) + (1.0 - alpha * r0_n) * chi2 * c + r0_n
        residual = f_val

        if f_val == 0.0:
            converged = True
            break
        # f is monotone, so the sign of f places chi on one side of the root (a non-finite f
        # comes from a chi far beyond it, on chi's own side of 0).
        if f_val < 0.0 or (f_val != f_val and chi < 0.0):
            if chi > lo:
                lo = chi
        elif chi < hi:
            hi = chi

        if f_prime > 0.0 and f_prime < inf and abs(f_val) < inf:
            chi_new = chi - f_val / f_prime
        else:
            chi_new = nan
        if chi_new == chi_new:
            delta = chi - chi_new
            if abs(delta) < _NEWTON_TOL_DELTA_REL * max(abs(chi_new), 1.0):
                chi = chi_new
                converged = True
                break
        # Bisect (rtsafe rule, Press et al., Numerical Recipes, sec. 9.4) when the Newton step
        # leaves the bracket or, once both sides are finite, is not at least halving the
        # previous step: on the exponential (hyperbolic) side Newton otherwise creeps by
        # about 1/sqrt(-alpha) per iteration.
        both_finite = lo > -inf and hi < inf
        if not (lo < chi_new < hi) or (both_finite and abs(chi_new - chi) > 0.5 * dx_last):
            if hi == inf:
                # far side still open: step out from lo by doubling, starting at the
                # parabolic bootstrap (which has the sign of dt)
                chi_new = max(2.0 * lo, chi_par)
            elif lo == -inf:
                chi_new = min(2.0 * hi, chi_par)
            else:
                chi_new = 0.5 * (lo + hi)
                if hi - lo <= _NEWTON_TOL_DELTA_REL * max(abs(chi_new), 1.0):
                    chi = chi_new
                    converged = True
                    break
        dx_last = abs(chi_new - chi)
        chi = chi_new
    if not converged:
        return nan, 0.0, 0.0, 0.0, 0.0

    z = chi * chi * alpha
    c = stumpff_c_py(z)
    s = stumpff_s_py(z)
    chi2 = chi * chi
    chi3 = chi2 * chi

    f_coef = 1.0 - (chi2 / r0_n) * c
    g_coef = dt - (chi3 / sqrt_mu) * s
    return chi, f_coef, g_coef, residual, z


# ---------------------------------------------------------------------------
# Universal-anomaly Newton core — JIT-compiled (#475)
# ---------------------------------------------------------------------------


@nb.njit(cache=True)  # type: ignore[untyped-decorator]
def _kepler_chi_newton(
    r0_n: float,
    v0_n: float,
    rv_dot: float,
    alpha: float,
    dt: float,
    mu: float,
    chi0: float,
) -> tuple[float, float, float, float, float]:
    """JIT-compiled universal-anomaly Newton solve.

    Identical logic to ``_kepler_chi_newton_py``; compiled with
    ``numba.njit(cache=True)`` for a 10-50x speedup on the Kepler propagation
    hot path.  Returns ``(chi, f, g, f_dot_z, g_dot_z)`` where the last two
    entries are the intermediate ``z = chi^2 * alpha`` and the Lagrange g
    coefficient; caller assembles the Lagrange f/g_dot from chi and r_n.

    Returns ``(nan, 0, 0, 0, 0)`` on non-convergence; caller raises
    :class:`KeplerConvergenceError`.
    """
    sqrt_mu = sqrt(mu)
    # #963: Newton safeguarded by a bracket. f(chi) is strictly increasing (f'(chi) = r > 0), and
    # f(0) = -sqrt(mu) dt, so the root lies in (0, inf) for dt > 0 and (-inf, 0) for dt < 0.
    # Every evaluated iterate narrows the bracket; a Newton step that leaves it is replaced by a
    # bisection (or, while the far side is still open, by a step that doubles the distance from
    # 0). A Newton step that converges is accepted before the bracket test, so iterates that
    # the plain Newton solve accepted are unchanged.
    if dt > 0.0:
        lo = 0.0
        hi = inf
    elif dt < 0.0:
        lo = -inf
        hi = 0.0
    else:
        lo = 0.0
        hi = 0.0
    chi = chi0
    converged = False
    dx_last = inf
    chi_par = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0
    if chi_par == 0.0:
        chi_par = 1.0 if dt > 0.0 else -1.0
    for _iteration in range(_NEWTON_MAX_ITER):
        z = chi * chi * alpha
        c = stumpff_c(z)
        s = stumpff_s(z)
        chi2 = chi * chi
        chi3 = chi2 * chi

        f_val = (
            (rv_dot / sqrt_mu) * chi2 * c
            + (1.0 - alpha * r0_n) * chi3 * s
            + r0_n * chi
            - sqrt_mu * dt
        )
        f_prime = (rv_dot / sqrt_mu) * chi * (1.0 - z * s) + (1.0 - alpha * r0_n) * chi2 * c + r0_n

        if f_val == 0.0:
            converged = True
            break
        # f is monotone, so the sign of f places chi on one side of the root (a non-finite f
        # comes from a chi far beyond it, on chi's own side of 0).
        if f_val < 0.0 or (f_val != f_val and chi < 0.0):
            if chi > lo:
                lo = chi
        elif chi < hi:
            hi = chi

        if f_prime > 0.0 and f_prime < inf and abs(f_val) < inf:
            chi_new = chi - f_val / f_prime
        else:
            chi_new = nan
        if chi_new == chi_new:
            delta = chi - chi_new
            if abs(delta) < _NEWTON_TOL_DELTA_REL * max(abs(chi_new), 1.0):
                chi = chi_new
                converged = True
                break
        # Bisect (rtsafe rule, Press et al., Numerical Recipes, sec. 9.4) when the Newton step
        # leaves the bracket or, once both sides are finite, is not at least halving the
        # previous step: on the exponential (hyperbolic) side Newton otherwise creeps by
        # about 1/sqrt(-alpha) per iteration.
        both_finite = lo > -inf and hi < inf
        if not (lo < chi_new < hi) or (both_finite and abs(chi_new - chi) > 0.5 * dx_last):
            if hi == inf:
                # far side still open: step out from lo by doubling, starting at the
                # parabolic bootstrap (which has the sign of dt)
                chi_new = max(2.0 * lo, chi_par)
            elif lo == -inf:
                chi_new = min(2.0 * hi, chi_par)
            else:
                chi_new = 0.5 * (lo + hi)
                if hi - lo <= _NEWTON_TOL_DELTA_REL * max(abs(chi_new), 1.0):
                    chi = chi_new
                    converged = True
                    break
        dx_last = abs(chi_new - chi)
        chi = chi_new
    if not converged:
        return nan, 0.0, 0.0, 0.0, 0.0

    z = chi * chi * alpha
    c = stumpff_c(z)
    s = stumpff_s(z)
    chi2 = chi * chi
    chi3 = chi2 * chi

    f_coef = 1.0 - (chi2 / r0_n) * c
    g_coef = dt - (chi3 / sqrt_mu) * s
    # Return (chi, f, g, unused_residual_placeholder, z) — residual not tracked
    # in the JIT path for speed; non-convergence returned nan above.
    return chi, f_coef, g_coef, 0.0, z


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _initial_chi_guess(r0_n: float, rv_dot: float, alpha: float, dt: float, mu: float) -> float:
    """Initial universal-anomaly guess per Vallado Algorithm 3.4.

    The hyperbolic log-form guess (eq. 3-66) is used only when it is finite and
    its logarithm argument is positive; otherwise the parabolic bootstrap is
    used. For near-parabolic orbits (``|alpha|`` tiny, ``sqrt(-a)`` huge) the
    log guess can land outside the Newton basin (#934); the callers retry once
    from the parabolic bootstrap in that case.
    """
    sqrt_mu = sqrt(mu)
    chi_par = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0
    if alpha > 1.0e-9:
        # Elliptic / near-circular: chi ~ sqrt(mu) * alpha * dt is a good start.
        return sqrt_mu * alpha * dt
    if alpha < -1.0e-9:
        a = 1.0 / alpha  # negative for hyperbolic
        sign_dt = 1.0 if dt >= 0.0 else -1.0
        denom = rv_dot + sign_dt * sqrt(-mu * a) * (1.0 - r0_n * alpha)
        arg = (-2.0 * mu * alpha * dt) / denom if denom != 0.0 else -1.0
        if arg > 0.0:
            chi_h = sign_dt * sqrt(-a) * log(arg)
            if isfinite(chi_h):
                return chi_h
        return chi_par
    # Parabolic: Newton recovers from the simple bootstrap.
    return chi_par


def propagate(
    r0: Vec3,
    v0: Vec3,
    dt: float,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[Vec3, Vec3]:
    """Two-body universal-variable propagation of ``(r0, v0)`` by ``dt`` seconds.

    Parameters
    ----------
    r0, v0:
        Inertial position and velocity at the start, ``(3,)`` float64 arrays
        in km and km/s.
    dt:
        Time interval in seconds. May be positive (forward), negative
        (backward), or zero (identity).
    mu:
        Central-body gravitational parameter, km^3/s^2. Defaults to the
        heliocentric :data:`cyclerfinder.core.constants.MU_SUN_KM3_S2`.

    Returns
    -------
    (r, v):
        Inertial state at ``t + dt`` as two ``(3,)`` float64 arrays.

    Raises
    ------
    KeplerConvergenceError
        If Newton iteration on the universal anomaly fails to converge.
    """
    if dt == 0.0:
        return r0.copy(), v0.copy()

    r0_arr = np.asarray(r0, dtype=np.float64)
    v0_arr = np.asarray(v0, dtype=np.float64)

    r0_n = float(np.linalg.norm(r0_arr))
    v0_n = float(np.linalg.norm(v0_arr))
    sqrt_mu = sqrt(mu)

    # alpha = 1 / a (reciprocal semi-major axis). Sign distinguishes regime:
    # > 0 elliptic, == 0 parabolic, < 0 hyperbolic.
    alpha = 2.0 / r0_n - (v0_n * v0_n) / mu

    rv_dot = float(np.dot(r0_arr, v0_arr))

    chi = _initial_chi_guess(r0_n, rv_dot, alpha, dt, mu)

    # JIT-compiled Newton iteration on f(chi) = 0 (Vallado eq. 3-65).
    chi_conv, f_coef, g_coef, _res, z = _kepler_chi_newton(r0_n, v0_n, rv_dot, alpha, dt, mu, chi)

    if chi_conv != chi_conv:  # nan check (numba-safe, no math.isnan needed)
        # #934: retry once from the parabolic bootstrap, which is in the basin
        # for near-parabolic orbits where the hyperbolic log guess is not.
        chi_p = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0
        if chi_p != chi:
            chi_conv, f_coef, g_coef, _res, z = _kepler_chi_newton(
                r0_n, v0_n, rv_dot, alpha, dt, mu, chi_p
            )
        if chi_conv != chi_conv:
            raise KeplerConvergenceError(chi, 0.0)

    r = f_coef * r0_arr + g_coef * v0_arr
    r_n = float(np.linalg.norm(r))

    chi2 = chi_conv * chi_conv
    c = float(stumpff_c(z))
    g_dot = 1.0 - (chi2 / r_n) * c
    f_dot = (sqrt_mu / (r_n * r0_n)) * chi_conv * (z * float(stumpff_s(z)) - 1.0)

    v = f_dot * r0_arr + g_dot * v0_arr

    return r, v
