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
* Iteration cap: 50.

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §3.3.
"""

from __future__ import annotations

from math import log, sqrt

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.constants import MU_SUN_KM3_S2

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

_NEWTON_TOL_DELTA_REL: float = 1.0e-12
_NEWTON_MAX_ITER: int = 50


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

    # Initial chi guess per Vallado Algorithm 3.4.
    chi: float
    if alpha > 1.0e-9:
        # Elliptic / near-circular: chi ~ sqrt(mu) * alpha * dt is a good start.
        chi = sqrt_mu * alpha * dt
    elif alpha < -1.0e-9:
        # Hyperbolic: log expression keeps the initial guess in the right basin.
        a = 1.0 / alpha  # negative for hyperbolic
        sign_dt = 1.0 if dt >= 0.0 else -1.0
        # Vallado eq. (3-66).
        arg = (-2.0 * mu * alpha * dt) / (rv_dot + sign_dt * sqrt(-mu * a) * (1.0 - r0_n * alpha))
        chi = sign_dt * sqrt(-a) * log(arg)
    else:
        # Parabolic: use a simple bootstrap; Newton recovers from here.
        # Semi-latus rectum p = |r0 x v0|^2 / mu.
        h_vec = np.cross(r0_arr, v0_arr)
        p = float(np.dot(h_vec, h_vec)) / mu
        # From Vallado, but simplified: use sqrt(p) * something. The iteration
        # is forgiving for parabolic; start from zero and Newton in.
        chi = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0
        # Silence "unused" lint hint while keeping `p` available for debugging.
        _ = p

    # Newton iteration on f(chi) = 0 (Vallado eq. 3-65). The residual is
    # measured in km (length scale of sqrt_mu * dt); convergence is judged
    # by the relative size of the chi step, which is dimensionless under
    # scaling and works equally well for short and long propagations.
    residual: float = 0.0
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
        # f'(chi); equivalent to r(chi) per universal-variable identities.
        f_prime = (rv_dot / sqrt_mu) * chi * (1.0 - z * s) + (1.0 - alpha * r0_n) * chi2 * c + r0_n

        if f_prime == 0.0:
            raise KeplerConvergenceError(chi, f_val)

        delta = f_val / f_prime
        chi -= delta
        residual = f_val

        # Relative-step convergence; absolute step for the chi=0 case.
        if abs(delta) < _NEWTON_TOL_DELTA_REL * max(abs(chi), 1.0):
            break
    else:
        raise KeplerConvergenceError(chi, residual)

    # Lagrange coefficients (Vallado eqs. 2-119 to 2-123).
    z = chi * chi * alpha
    c = stumpff_c(z)
    s = stumpff_s(z)
    chi2 = chi * chi
    chi3 = chi2 * chi

    f = 1.0 - (chi2 / r0_n) * c
    g = dt - (chi3 / sqrt_mu) * s

    r = f * r0_arr + g * v0_arr
    r_n = float(np.linalg.norm(r))

    g_dot = 1.0 - (chi2 / r_n) * c
    f_dot = (sqrt_mu / (r_n * r0_n)) * chi * (z * s - 1.0)

    v = f_dot * r0_arr + g_dot * v0_arr

    return r, v
