"""Analytic CR3BP planar periodic-orbit seed generator at arbitrary mu (#435).

Builds a small-amplitude *linear* seed for a planar periodic orbit and refines
it with the existing fixed-Jacobi symmetric corrector
(:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`).
Two families:

* :func:`lyapunov_seed` -- planar Lyapunov orbit about a collinear libration
  point (L1/L2), from the standard collinear linearization (Szebehely 1967;
  Koon, Lo, Marsden & Ross, *Dynamical Systems, the Three-Body Problem and
  Space Mission Design*, 2011, Ch. 2.5-2.7).
* :func:`dro_seed` -- distant retrograde orbit about the secondary, from a
  retrograde circular two-body guess clamped to stay bound.

Frame convention (matches :mod:`cyclerfinder.core.cr3bp`): nondimensional
rotating frame, primary at ``(-mu, 0, 0)``, secondary at ``(1-mu, 0, 0)``.

This feeds the #435 high-eccentricity Sun-planet ER3BP discovery, where a
converged CR3BP planar orbit at arbitrary mu is the continuation anchor.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.cr3bp import CR3BPSystem, jacobi_constant
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x

# Amplitude ladder (smallest-first) for the Lyapunov linear seed. The corrector's
# convergence basin is non-monotonic in amplitude and μ-dependent: at some μ the
# nominal 1e-3 lands in a non-convergent gap while smaller amplitudes (deeper in
# the linear regime, genuine small-amplitude Lyapunov, period ~ 2π/ω) converge
# cleanly. We try the caller's amplitude first, then walk this ladder.
_LYAPUNOV_AMPLITUDE_LADDER = (1.0e-4, 5.0e-5, 1.0e-5, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)


def _collinear_c2(x_l: float, mu: float) -> float:
    """In-plane linear coefficient c2 at a collinear point.

    ``c2 = (1-mu)/|x_L + mu|^3 + mu/|x_L - (1-mu)|^3`` -- the standard
    collinear-linearization coefficient (Koon-Lo-Marsden-Ross 2011, Eq. 2.5.x):
    the sum of the two primaries' gravity-gradient terms at ``x_L``.
    """
    d1 = abs(x_l + mu)
    d2 = abs(x_l - (1.0 - mu))
    return (1.0 - mu) / d1**3 + mu / d2**3


def lyapunov_seed(
    system: CR3BPSystem,
    *,
    point: str = "L1",
    amplitude: float = 1e-3,
) -> tuple[NDArray[np.float64], float]:
    """Converged planar Lyapunov orbit about a collinear point, at arbitrary mu.

    Construct the textbook small-amplitude linear seed at the collinear point
    ``point`` ("L1" or "L2") and refine it with the fixed-Jacobi symmetric
    corrector. ``amplitude`` is the nondimensional in-plane x-amplitude ``Ax``
    of the linear seed (small, e.g. ``1e-3``).

    The collinear linearization gives the in-plane frequency

        omega = sqrt((2 - c2 + sqrt(9 c2^2 - 8 c2)) / 2)

    and a planar Lyapunov initial condition (x-axis perpendicular crossing)

        x0  = x_L - Ax,   y0=z0=vx0=vz0 = 0,
        vy0 = -Ax * omega * tau,   tau = -(omega^2 + 1 + 2 c2) / (2 omega).

    Returns
    -------
    (state0, period)
        The corrected 6-state (perpendicular x-axis crossing) and the FULL
        nondimensional period ``T = 2 * t_half``.

    Raises
    ------
    ValueError
        If ``point`` is not a collinear point with an in-plane oscillatory mode
        (radicand ``9 c2^2 - 8 c2 < 0``), or if the corrector fails to converge.
    """
    mu = system.mu
    x_l = lagrange_collinear_x(mu, point)
    c2 = _collinear_c2(x_l, mu)
    radicand = 9.0 * c2 * c2 - 8.0 * c2
    if radicand < 0.0:
        raise ValueError(
            f"lyapunov_seed: no in-plane oscillatory mode at {point} "
            f"(9 c2^2 - 8 c2 = {radicand:.3e} < 0, c2={c2:.6f}, mu={mu:.3e})"
        )
    omega = math.sqrt((2.0 - c2 + math.sqrt(radicand)) / 2.0)
    tau = -(omega * omega + 1.0 + 2.0 * c2) / (2.0 * omega)
    period_guess = 2.0 * math.pi / omega

    # Try the requested amplitude first, then the smallest-first ladder. The
    # corrector's basin is non-monotonic in (amplitude, μ); walking the ladder
    # makes the generator reliable across Sun-Mercury/Mars/Pluto/Earth-Moon.
    tried: list[float] = []
    last_residual = math.inf
    last_iter = 0
    for ax in (float(amplitude), *_LYAPUNOV_AMPLITUDE_LADDER):
        if any(abs(ax - t) <= 1e-15 for t in tried):
            continue
        tried.append(ax)
        x0 = x_l - ax
        vy0 = -ax * omega * tau
        seed_state = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0], dtype=np.float64)
        jacobi = jacobi_constant(seed_state, mu)
        orbit = correct_symmetric_fixed_jacobi(
            system,
            x0_guess=x0,
            jacobi=jacobi,
            period_guess=period_guess,
            ydot0_sign=math.copysign(1.0, vy0) if vy0 != 0.0 else 1.0,
            half_crossings=1,
        )
        if orbit.converged:
            state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
            return state0, orbit.period
        last_residual = orbit.crossing_residual
        last_iter = orbit.n_iter
    raise ValueError(
        f"lyapunov_seed: corrector did not converge at {point} for any amplitude "
        f"in {tried} (last residual={last_residual:.3e}, n_iter={last_iter}, mu={mu:.3e})"
    )


def dro_seed(
    system: CR3BPSystem,
    *,
    amplitude: float = 5e-2,
) -> tuple[NDArray[np.float64], float]:
    """Converged distant retrograde orbit (DRO) about the secondary.

    Retrograde co-orbital seed near the secondary (at ``(1-mu, 0, 0)``): a
    perpendicular x-axis crossing exterior to the secondary at radius
    ``amplitude``, with a retrograde circular two-body guess
    ``v ~ sqrt(mu / amplitude)`` (clamped to stay bound), refined with the
    fixed-Jacobi symmetric corrector.

    Returns
    -------
    (state0, period)
        The corrected 6-state and FULL nondimensional period.

    Raises
    ------
    ValueError
        If the corrector fails to converge.
    """
    mu = system.mu
    a = float(amplitude)
    x0 = (1.0 - mu) + a
    # Retrograde circular two-body speed about the secondary, clamped so the
    # guess stays bound (escape ~ sqrt(2 mu / a)); retrograde => negative vy.
    v_circ = math.sqrt(mu / a)
    vy0 = -v_circ
    seed_state = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0], dtype=np.float64)
    jacobi = jacobi_constant(seed_state, mu)
    period_guess = 2.0 * math.pi * a / v_circ

    orbit = correct_symmetric_fixed_jacobi(
        system,
        x0_guess=x0,
        jacobi=jacobi,
        period_guess=period_guess,
        ydot0_sign=-1.0,
        half_crossings=1,
    )
    if not orbit.converged:
        raise ValueError(
            f"dro_seed: corrector did not converge "
            f"(residual={orbit.crossing_residual:.3e}, n_iter={orbit.n_iter})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    return state0, orbit.period
