"""Patched-conic gravity-assist mechanics.

Pure functions on heliocentric :math:`V_\\infty` vectors. The math primitives
(:func:`max_bend`, :func:`bend_angle`, :func:`is_ballistic_feasible`,
:func:`flyby_dv`) take explicit ``mu_planet`` and ``rp_min`` arguments and do
not reach into :mod:`cyclerfinder.core.constants`; the planet-aware
convenience wrapper :func:`flyby_dv_for` is the single touch-point with the
planet table.

Bend formula (spec §3, §6 ; Bate-Mueller-White §6.4):

.. math::

    \\sin\\!\\left(\\frac{\\delta_{\\max}}{2}\\right)
    = \\frac{1}{1 + r_p V_\\infty^2 / \\mu}

Single-impulse powered-flyby surrogate (Strange & Longuski, JSR 2002,
eq. 9; matches the plan §3.1 sketch): the flyby cost has two non-negative
components — a magnitude change ``|V_out| - |V_in|`` paid as a tangential
periapsis burn, and a bend deficit ``max(0, delta - delta_max)`` paid as a
second impulse rotating the asymptote. The cost is exactly zero on a
ballistic-feasible pair (equal-magnitude in-bounds bend) and strictly
positive otherwise. M5 may refine this if the optimiser demands it.

Plan: ``docs/phases/m2-flyby-maps/plan.md`` §3.1.
"""

from __future__ import annotations

from math import asin, sin

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64 — matches lambert.py


def max_bend(mu_planet: float, rp_min: float, vinf: float) -> float:
    """Maximum ballistic deflection angle (rad) for a hyperbolic flyby.

    ``sin(delta_max / 2) = 1 / (1 + rp_min * vinf**2 / mu_planet)``.

    Limits: returns ``pi`` as ``vinf -> 0`` (the asymptote bends all the way
    around in the limit of an infinitely slow flyby), and approaches ``0``
    as ``vinf -> infinity`` (a fast flyby cannot be turned).

    Parameters
    ----------
    mu_planet:
        Planet gravitational parameter ``GM``, km**3/s**2.
    rp_min:
        Minimum safe periapsis radius, km (planet radius plus safety altitude).
    vinf:
        Hyperbolic excess speed, km/s. Must be non-negative.

    Returns
    -------
    float
        Maximum deflection angle in radians, in ``[0, pi]``.
    """
    if vinf < 0.0:
        raise ValueError(f"vinf must be non-negative, got {vinf}")
    if vinf == 0.0:
        # Limiting case: 1/(1 + 0) = 1, asin(1) = pi/2, doubled to pi.
        return float(np.pi)
    arg = 1.0 / (1.0 + rp_min * vinf * vinf / mu_planet)
    # arg is in (0, 1] by construction; clip defensively for FP noise.
    return 2.0 * asin(min(1.0, max(0.0, arg)))


def dv_from_turn_deficit(vinf: float, delta_required: float, delta_max: float) -> float:
    """``Delta V`` (km/s) to make up a turn-angle deficit at a flyby.

    When the geometry demands more bending than the planet can ballistically
    provide (``delta_required > delta_max``), the shortfall is paid as a
    single periapsis impulse rotating the :math:`V_\\infty` asymptote
    (Strange & Longuski, JSR 2002, eq. 9):

    .. math::

        \\Delta V = 2 V_\\infty \\sin\\!\\left(\\tfrac{1}{2}
                    \\max(0,\\, \\delta_\\text{req} - \\delta_{\\max})\\right)

    Returns ``0.0`` when the required turn is within the achievable cone.
    This is the magnitude reported as *interesting output* for powered
    cyclers such as the classic Aldrin (1L1): McConaghy 2002 Table 4 gives the
    **Earth** (geocentric) ``delta_required = 84 deg`` against ``delta_max =
    72 deg`` (a 200 km Earth flyby), a 12 deg deficit. NOTE: this is a
    single-impulse surrogate that *over*-estimates
    the optimal maintenance cost; it is our computed value, not a
    source-attested ``Delta V`` (no published Aldrin ``Delta V`` exists —
    McConaghy defers it, p.8). Never use it as a golden-test target.

    Parameters
    ----------
    vinf:
        Hyperbolic excess speed at the flyby, km/s. Non-negative.
    delta_required:
        Turn angle the trajectory geometry demands, rad.
    delta_max:
        Maximum ballistically achievable turn angle, rad (see
        :func:`max_bend`).

    Returns
    -------
    float
        ``Delta V`` in km/s. Non-negative; ``0.0`` when within the cone.
    """
    if vinf < 0.0:
        raise ValueError(f"vinf must be non-negative, got {vinf}")
    excess = max(0.0, delta_required - delta_max)
    return 2.0 * vinf * sin(0.5 * excess)


def bend_angle(vin_vec: Vec3, vout_vec: Vec3) -> float:
    """Angle (rad) between two :math:`V_\\infty` vectors.

    Numerically robust via clipping the ``acos`` argument to ``[-1, 1]``.
    Both vectors must be non-zero; the result is the unsigned angle between
    them (the planar bend in the b-plane sense).

    Parameters
    ----------
    vin_vec, vout_vec:
        Length-3 heliocentric :math:`V_\\infty` vectors, km/s. Magnitudes
        should be (ideally) equal for a ballistic flyby; this routine just
        reports the angle and makes no feasibility judgement.

    Returns
    -------
    float
        Angle in radians, in ``[0, pi]``.
    """
    vin_norm = float(np.linalg.norm(vin_vec))
    vout_norm = float(np.linalg.norm(vout_vec))
    if vin_norm == 0.0 or vout_norm == 0.0:
        raise ValueError("bend_angle requires non-zero V_inf vectors")
    cos_arg = float(np.dot(vin_vec, vout_vec)) / (vin_norm * vout_norm)
    return float(np.arccos(np.clip(cos_arg, -1.0, 1.0)))


def is_ballistic_feasible(
    vin_vec: Vec3,
    vout_vec: Vec3,
    mu_planet: float,
    rp_min: float,
    speed_tol: float = 1.0e-6,
) -> bool:
    """True iff a ballistic flyby can map ``vin_vec`` to ``vout_vec``.

    Requires (a) the speeds match within ``speed_tol`` (km/s), and (b) the
    rotation between the two vectors is within the achievable bend cone
    ``max_bend(mu_planet, rp_min, |vin|)`` for the (shared) speed.

    Parameters
    ----------
    vin_vec, vout_vec:
        Length-3 :math:`V_\\infty` vectors, km/s.
    mu_planet:
        Planet ``GM``, km**3/s**2.
    rp_min:
        Minimum safe periapsis radius, km.
    speed_tol:
        Allowed magnitude mismatch in km/s. Default ``1e-6`` matches the
        Lambert residual scale from M1.

    Returns
    -------
    bool
    """
    vin_norm = float(np.linalg.norm(vin_vec))
    vout_norm = float(np.linalg.norm(vout_vec))
    if abs(vin_norm - vout_norm) > speed_tol:
        return False
    delta = bend_angle(vin_vec, vout_vec)
    delta_max = max_bend(mu_planet, rp_min, vin_norm)
    return delta <= delta_max


def flyby_dv(
    vin_vec: Vec3,
    vout_vec: Vec3,
    mu_planet: float,
    rp_min: float,
) -> float:
    """Powered-flyby ``Delta V`` required to convert ``vin_vec`` into ``vout_vec``.

    Returns exactly ``0.0`` when the requested transformation is
    ballistic-feasible (see :func:`is_ballistic_feasible`). Otherwise the
    single-impulse surrogate is:

    .. math::

        \\Delta V = \\underbrace{\\bigl||V_\\text{out}| - |V_\\text{in}|\\bigr|}_{\\text{magnitude}}
                  + \\underbrace{2 \\bar V_\\infty \\sin\\!\\left(\\tfrac{1}{2}
                    \\max(0,\\, \\delta - \\delta_{\\max})\\right)}_{\\text{excess bend}}

    where ``\\bar V_\\infty`` is the mean of the in/out magnitudes and
    ``delta_max`` is evaluated at that mean. The decomposition matches the
    spec §6 sketch and is exactly zero on the ballistic-feasible branch.

    Parameters
    ----------
    vin_vec, vout_vec:
        Length-3 :math:`V_\\infty` vectors, km/s.
    mu_planet:
        Planet ``GM``, km**3/s**2.
    rp_min:
        Minimum safe periapsis radius, km.

    Returns
    -------
    float
        ``Delta V`` in km/s. Non-negative; ``0.0`` exactly on the
        ballistic-feasible branch.
    """
    if is_ballistic_feasible(vin_vec, vout_vec, mu_planet, rp_min):
        return 0.0

    vin_norm = float(np.linalg.norm(vin_vec))
    vout_norm = float(np.linalg.norm(vout_vec))
    v_mean = 0.5 * (vin_norm + vout_norm)

    mag_cost = abs(vout_norm - vin_norm)
    delta = bend_angle(vin_vec, vout_vec)
    delta_achievable = max_bend(mu_planet, rp_min, v_mean)
    bend_cost = dv_from_turn_deficit(v_mean, delta, delta_achievable)

    return mag_cost + bend_cost


def flyby_dv_for(
    code: str,
    vin_vec: Vec3,
    vout_vec: Vec3,
) -> float:
    """Planet-aware convenience wrapper around :func:`flyby_dv`.

    Looks up ``mu_planet`` from :data:`cyclerfinder.core.constants.PLANETS`
    and ``rp_min`` from :data:`cyclerfinder.core.constants.SAFE_PERIHELION_KM`
    for the given one-letter planet code.

    Parameters
    ----------
    code:
        One-letter planet code (``"V"``, ``"E"``, ``"M"`` for the M2 set).
    vin_vec, vout_vec:
        Length-3 :math:`V_\\infty` vectors, km/s.

    Returns
    -------
    float
        ``Delta V`` in km/s.

    Raises
    ------
    KeyError
        If ``code`` is not in :data:`PLANETS`.
    """
    mu = PLANETS[code].mu_km3_s2
    rp = SAFE_PERIHELION_KM[code]
    return flyby_dv(vin_vec, vout_vec, mu, rp)
