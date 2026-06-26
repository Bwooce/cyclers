"""Analytic two-body state-transition matrix (universal variables, Shepperd lane).

Propagates an inertial state ``(r0, v0)`` by a signed interval ``dt`` under a
central body ``mu`` — exactly like :func:`cyclerfinder.core.kepler.propagate` —
and additionally returns the exact 6x6 state-transition matrix (STM)::

    Phi(t0 + dt, t0) = [ dr/dr0  dr/dv0 ]   =   [ R~  R ]
                       [ dv/dr0  dv/dv0 ]       [ V~  V ]

in the quadrant notation of Ellison et al. 2018, Eq. (17) (Battin's notation,
their ref [51]). One universal-variable formulation covers elliptic, parabolic,
and hyperbolic regimes, multi-revolution arcs, and backward propagation. This is
the gating ingredient for the Path-B MGAnDSMs re-transcription of the DSM lane
(#225/#226) and gives exact analytic sparsity to any Kepler-arc corrector.

Sources
-------
Shepperd, S. W., "Universal Keplerian State Transition Matrix", *Celestial
Mechanics* 35 (1985), pp. 129-144, doi:10.1007/BF01227666 — the canonical
reference this lane standardises on, cited via Ellison, D. H., Conway, B. A.,
Englander, J. A., and Ozimek, M. T., "Analytic Gradient Computation for
Bounded-Impulse Trajectory Models Using Two-Sided Shooting," *JGCD* 41(7),
2018, Eq. (17) and ref [52] (mining note
``docs/notes/2026-06-10-ellison-2018-analytic-gradients-mining.md`` §3, §7).

SOURCE HONESTY: Ellison 2018 prints only the quadrant *partition* (Eq. 17),
not the closed forms — its appendix carries the flyby gradients (A1-A6) only.
The quadrant expressions implemented below are therefore the standard
universal-variable closed form of the Keplerian STM in Battin's development
(equivalent to Shepperd's), with the secular term
``C = (3*U5 - chi*U4)/sqrt(mu) - dt*U2``. Validation is the CONSISTENCY-test
pattern (exactly like ``nbody/flyby_gradients.py``): central-difference FD
agreement against the independent :func:`cyclerfinder.core.kepler.propagate`
across all conic regimes, plus the STM group properties (identity at
``dt = 0``, symplecticity ``Phi^T J Phi = J``, ``det Phi = 1``, composition).

Shepperd 1985 read 2026-06-13 (#233, closes the #116 acquisition item for this
paper): the paper is pure theory/algorithm — Sections 2-7 give the symbolic
universal-variable formulation, the M-matrix STM (Eq. 17), the U_n functions
(Eqs. 20-30), and the Gaussian-continued-fraction evaluation, and Appendix A
summarises the algorithm (Eqs. A.1-A.46). It prints NO worked numeric example:
no input state, no propagated output, no numeric STM cells, no tabulated
U-function values anywhere in pp. 129-144. There is therefore no wireable
printed golden; the consistency-tests above remain the validation. No goldens
are fabricated here.

Universal functions
-------------------
With ``alpha = 1/a = 2/r0 - v0^2/mu`` (regime sign: > 0 elliptic, = 0
parabolic, < 0 hyperbolic), universal anomaly ``chi`` (units sqrt(km), same
convention as ``core/kepler.py``), and ``z = alpha * chi**2``::

    U_n(chi; alpha) = chi**n * c_n(z)          n = 0..5

where ``c_n`` are the Stumpff functions: ``c0 = cos`` / ``c1 = sinc`` analogues
with ``c0(z) = 1 - z*c2(z)``, ``c1(z) = 1 - z*c3(z)``; ``c2``, ``c3`` come from
:mod:`cyclerfinder.core._stumpff` (closed trig/hyperbolic forms with a series
guard near ``z = 0``); ``c4``, ``c5`` are evaluated here (see below). The
identities used: ``r = r0*U0 + (r0.v0/sqrt(mu))*U1 + U2`` and the Lagrange
coefficients ``f = 1 - U2/r0``, ``g = dt - U3/sqrt(mu)``,
``fdot = -sqrt(mu)*U1/(r*r0)``, ``gdot = 1 - U2/r``.

c4/c5 evaluation and convergence guard
--------------------------------------
Shepperd's own paper evaluates the transcendental kernel via a Gautschi
continued fraction; here the equivalent ``c4``/``c5`` values come from the
downward Stumpff recursion ``c_{k+2}(z) = (1/k! - c_k(z))/z``, which is exact
algebra but loses relative precision to cancellation as ``z -> 0`` (error ~
``12*eps/|z|``). Inside ``|z| < 1e-3`` (matching the ``_stumpff`` cutoff) the
Maclaurin series ``c_k(z) = sum_i (-z)^i / (k + 2i)!`` is used instead,
truncated at ``z**3`` (next term < 1e-19 at the cutoff). At the cutoff the
recursion's cancellation error is ~2e-13 relative — both branches agree well
inside the FD validation tolerance.

Tolerances
----------
Newton solve for ``chi``: identical to ``core/kepler.py`` (relative chi-step
``< 1e-12``, 50-iteration cap, Vallado Alg. 3.4 initial guesses). The solver is
replicated inline because the STM needs ``chi`` itself, which
:func:`~cyclerfinder.core.kepler.propagate` does not expose (and that module is
not modified here).
"""

from __future__ import annotations

from math import log, sqrt

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler import KeplerConvergenceError

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64
Mat6 = NDArray[np.float64]  # shape (6, 6), dtype float64

_NEWTON_TOL_DELTA_REL: float = 1.0e-12
_NEWTON_MAX_ITER: int = 50

# Series cutoff for c4/c5; below this the downward recursion from c2/c3 loses
# ~12*eps/|z| relative precision to cancellation (see module docstring).
_C45_SERIES_CUTOFF: float = 1.0e-3


def _stumpff_c4(z: float) -> float:
    """Stumpff ``c4(z)``: series inside the guard, downward recursion outside."""
    if abs(z) < _C45_SERIES_CUTOFF:
        # c4(z) = 1/24 - z/720 + z**2/40320 - z**3/3628800 + ...
        return 1.0 / 24.0 - z / 720.0 + z * z / 40320.0 - (z * z * z) / 3628800.0
    # c4 = (1/2! - c2(z)) / z
    return float(0.5 - stumpff_c(z)) / z


def _stumpff_c5(z: float) -> float:
    """Stumpff ``c5(z)``: series inside the guard, downward recursion outside."""
    if abs(z) < _C45_SERIES_CUTOFF:
        # c5(z) = 1/120 - z/5040 + z**2/362880 - z**3/39916800 + ...
        return 1.0 / 120.0 - z / 5040.0 + z * z / 362880.0 - (z * z * z) / 39916800.0
    # c5 = (1/3! - c3(z)) / z
    return float(1.0 / 6.0 - stumpff_s(z)) / z


def shepperd_stm(
    r0: Vec3,
    v0: Vec3,
    dt: float,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[Vec3, Vec3, Mat6]:
    """Propagate ``(r0, v0)`` by ``dt`` and return the analytic 6x6 STM.

    Parameters
    ----------
    r0, v0:
        Inertial position and velocity at the start, ``(3,)`` float64 arrays
        in km and km/s.
    dt:
        Time interval in seconds. May be positive (forward), negative
        (backward), or zero (identity state, identity STM).
    mu:
        Central-body gravitational parameter, km^3/s^2. Defaults to the
        heliocentric :data:`cyclerfinder.core.constants.MU_SUN_KM3_S2`.

    Returns
    -------
    (r, v, phi):
        Inertial state at ``t + dt`` (two ``(3,)`` float64 arrays, matching
        :func:`cyclerfinder.core.kepler.propagate` to numerical noise) and the
        ``(6, 6)`` state-transition matrix ``Phi(t0 + dt, t0)`` ordered
        ``[r; v]`` — position rows/columns first. Units: the ``dr/dv0`` block
        carries seconds, the ``dv/dr0`` block 1/seconds.

    Raises
    ------
    KeplerConvergenceError
        If Newton iteration on the universal anomaly fails to converge
        (same failure mode and exception type as the plain propagator).
    """
    r0_arr = np.asarray(r0, dtype=np.float64)
    v0_arr = np.asarray(v0, dtype=np.float64)

    if dt == 0.0:
        return r0_arr.copy(), v0_arr.copy(), np.eye(6, dtype=np.float64)

    r0_n = float(np.linalg.norm(r0_arr))
    v0_n = float(np.linalg.norm(v0_arr))
    sqrt_mu = sqrt(mu)

    # alpha = 1 / a (reciprocal semi-major axis); sign distinguishes regime.
    alpha = 2.0 / r0_n - (v0_n * v0_n) / mu
    rv_dot = float(np.dot(r0_arr, v0_arr))

    # --- Newton solve for the universal anomaly chi (mirrors core/kepler.py,
    # Vallado Algorithm 3.4 initial guesses) -------------------------------
    chi: float
    if alpha > 1.0e-9:
        chi = sqrt_mu * alpha * dt
    elif alpha < -1.0e-9:
        a = 1.0 / alpha  # negative for hyperbolic
        sign_dt = 1.0 if dt >= 0.0 else -1.0
        arg = (-2.0 * mu * alpha * dt) / (rv_dot + sign_dt * sqrt(-mu * a) * (1.0 - r0_n * alpha))
        chi = sign_dt * sqrt(-a) * log(arg)
    else:
        chi = sqrt_mu * dt / r0_n if r0_n > 0.0 else 0.0

    residual: float = 0.0
    for _iteration in range(_NEWTON_MAX_ITER):
        z = chi * chi * alpha
        c2 = stumpff_c(z)
        c3 = stumpff_s(z)
        chi2 = chi * chi
        chi3 = chi2 * chi

        f_val = (
            (rv_dot / sqrt_mu) * chi2 * c2
            + (1.0 - alpha * r0_n) * chi3 * c3
            + r0_n * chi
            - sqrt_mu * dt
        )
        f_prime = (
            (rv_dot / sqrt_mu) * chi * (1.0 - z * c3) + (1.0 - alpha * r0_n) * chi2 * c2 + r0_n
        )

        if f_prime == 0.0:
            raise KeplerConvergenceError(chi, f_val)

        delta = f_val / f_prime
        chi -= delta
        residual = f_val

        if abs(delta) < _NEWTON_TOL_DELTA_REL * max(abs(chi), 1.0):
            break
    else:
        raise KeplerConvergenceError(chi, residual)

    # --- Universal functions U0..U5 at the converged chi ------------------
    z = chi * chi * alpha
    c2 = stumpff_c(z)
    c3 = stumpff_s(z)
    u0 = 1.0 - z * c2  # c0(z)
    u1 = chi * (1.0 - z * c3)  # chi * c1(z)
    u2 = chi * chi * c2
    u3 = chi * chi * chi * c3
    u4 = chi**4 * _stumpff_c4(z)
    u5 = chi**5 * _stumpff_c5(z)

    # --- Lagrange coefficients and propagated state ------------------------
    r_n = r0_n * u0 + (rv_dot / sqrt_mu) * u1 + u2
    f = 1.0 - u2 / r0_n
    g = dt - u3 / sqrt_mu
    f_dot = -sqrt_mu * u1 / (r_n * r0_n)
    g_dot = 1.0 - u2 / r_n

    r_arr = f * r0_arr + g * v0_arr
    v_arr = f_dot * r0_arr + g_dot * v0_arr

    # --- STM quadrants (Battin-form universal-variable closed form) --------
    # Secular term (km*s); carries the multi-revolution sensitivity growth.
    c_sec = (3.0 * u5 - chi * u4) / sqrt_mu - dt * u2

    dr = r_arr - r0_arr
    dv = v_arr - v0_arr
    eye3 = np.eye(3, dtype=np.float64)

    # R~ = dr/dr0
    r_tilde = (
        (r_n / mu) * np.outer(dv, dv)
        + (1.0 / r0_n**3)
        * (r0_n * (1.0 - f) * np.outer(r_arr, r0_arr) + c_sec * np.outer(v_arr, r0_arr))
        + f * eye3
    )
    # R = dr/dv0
    r_mat = (
        (r0_n / mu) * (1.0 - f) * (np.outer(dr, v0_arr) - np.outer(dv, r0_arr))
        + (c_sec / mu) * np.outer(v_arr, v0_arr)
        + g * eye3
    )
    # V = dv/dv0
    v_mat = (
        (r0_n / mu) * np.outer(dv, dv)
        + (1.0 / r_n**3)
        * (r0_n * (1.0 - f) * np.outer(r_arr, r0_arr) - c_sec * np.outer(r_arr, v0_arr))
        + g_dot * eye3
    )
    # V~ = dv/dr0 (the involved quadrant: the f_dot bracket carries the
    # angular-momentum-direction correction (r v^T - v r^T) r (dv)^T / (mu r)).
    v_tilde = (
        -np.outer(dv, r0_arr) / r0_n**2
        - np.outer(r_arr, dv) / r_n**2
        + f_dot
        * (
            eye3
            - np.outer(r_arr, r_arr) / r_n**2
            + (1.0 / (mu * r_n))
            * np.outer((np.outer(r_arr, v_arr) - np.outer(v_arr, r_arr)) @ r_arr, dv)
        )
        - (mu * c_sec / (r_n**3 * r0_n**3)) * np.outer(r_arr, r0_arr)
    )

    phi = np.empty((6, 6), dtype=np.float64)
    phi[0:3, 0:3] = r_tilde
    phi[0:3, 3:6] = r_mat
    phi[3:6, 0:3] = v_tilde
    phi[3:6, 3:6] = v_mat

    return r_arr, v_arr, phi
