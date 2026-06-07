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

from math import asin, sin, sqrt

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


def dv_powered_flyby_periapsis(
    vinf: float,
    delta_required: float,
    delta_max: float,
    mu_planet: float,
    rp_min: float,
) -> float:
    """Oberth-credited periapsis powered-flyby ``Delta V`` (km/s) for a turn deficit.

    An **alternative** to :func:`dv_from_turn_deficit` for the equal-:math:`|V_\\infty|`
    turn-closure case (the actuator the Aldrin maintenance schedule needs). Where
    :func:`dv_from_turn_deficit` rotates the :math:`V_\\infty` asymptote *at
    infinity* (speed :math:`V_\\infty`, no Oberth credit), this model supplies the
    residual turn by a **tangential impulse at periapsis**, deep in the planet's
    well, where the spacecraft moves at
    :math:`v_p = \\sqrt{V_\\infty^2 + 2\\mu/r_p}`.

    Mechanism (Takao 2025 Eq. 11 / Russell 2004 Eq. 5.5 family). The achievable
    ballistic bend cone at ``rp_min`` is ``delta_max`` for the incoming
    :math:`V_\\infty`. To deliver the *full* ``delta_required`` ballistically the
    flyby must run at a **lower** excess speed ``vinf_target`` for which the cone
    just opens to ``delta_required``:

    .. math::

        \\sin\\!\\left(\\tfrac{1}{2}\\delta_\\text{req}\\right)
        = \\frac{1}{1 + r_p\\, V_{\\infty,\\text{target}}^2 / \\mu}
        \\;\\Rightarrow\\;
        V_{\\infty,\\text{target}}
        = \\sqrt{\\frac{\\mu}{r_p}\\!\\left(\\frac{1}{\\sin(\\delta_\\text{req}/2)} - 1\\right)}.

    The maneuver is two tangential periapsis impulses charged by the Oberth speed
    difference (Takao Eq. 11): slow from the incoming periapsis speed
    :math:`v_p(V_\\infty)` to :math:`v_p(V_{\\infty,\\text{target}})` so the
    widened cone delivers the whole turn, then restore the original
    :math:`|V_\\infty|` on the way out so the closure magnitude is preserved:

    .. math::

        \\Delta V_\\text{Oberth}
        = 2\\,\\bigl|\\,\\sqrt{V_\\infty^2 + 2\\mu/r_p}
                       - \\sqrt{V_{\\infty,\\text{target}}^2 + 2\\mu/r_p}\\,\\bigr|.

    Properties (asserted in :mod:`tests.core.test_flyby_oberth`): zero deficit ->
    exactly ``0.0``; monotone non-decreasing in ``delta_required``; and, **in the
    deep-well regime** (``2 mu / rp`` dominant — equivalently ``vinf`` below a
    body-specific threshold, which covers the entire physically-relevant deficit
    range for an Earth flyby up to ``vinf ~ 6.9 km/s``), strictly below
    :func:`dv_from_turn_deficit` for the same deficit (the Oberth credit). The
    guarantee is **not** universal: at high ``vinf`` the ballistic cone is already
    so narrow that opening it to ``delta_required`` demands a large magnitude
    excursion and the periapsis maneuver can exceed the asymptote rotation. See
    the module-level note ``docs/notes/2026-06-07-oberth-flyby-recost.md``.

    Units/frames match the rest of the module: ``vinf`` and the result in km/s,
    ``mu_planet`` in km**3/s**2, ``rp_min`` in km, angles in radians.

    Parameters
    ----------
    vinf:
        Hyperbolic excess speed at the flyby (shared in/out magnitude), km/s.
        Non-negative.
    delta_required:
        Turn angle the geometry demands, rad.
    delta_max:
        Maximum ballistically achievable turn at ``rp_min`` for ``vinf`` (rad);
        see :func:`max_bend`. Passed in so the caller controls the cone
        convention (the same value fed to :func:`dv_from_turn_deficit`).
    mu_planet:
        Planet gravitational parameter ``GM``, km**3/s**2.
    rp_min:
        Periapsis radius at which the maneuver is performed (the tightest safe
        flyby), km.

    Returns
    -------
    float
        ``Delta V`` in km/s. Non-negative; exactly ``0.0`` within the cone.
    """
    if vinf < 0.0:
        raise ValueError(f"vinf must be non-negative, got {vinf}")
    if mu_planet <= 0.0:
        raise ValueError(f"mu_planet must be positive, got {mu_planet}")
    if rp_min <= 0.0:
        raise ValueError(f"rp_min must be positive, got {rp_min}")
    if delta_required <= delta_max:
        return 0.0
    # vinf_target whose ballistic cone equals the required turn. delta_required
    # is in (delta_max, pi) here, so sin(delta_required / 2) in (0, 1] and the
    # bracket is non-negative.
    s = sin(0.5 * delta_required)
    inner = (1.0 / s) - 1.0
    if inner <= 0.0:
        # delta_required >= pi: the cone can never open this far; no finite
        # periapsis target exists. Fall back to the asymptote-rotation cost so
        # the model never returns a spuriously small number.
        return dv_from_turn_deficit(vinf, delta_required, delta_max)
    vinf_target = sqrt(mu_planet / rp_min * inner)
    vp_in = sqrt(vinf * vinf + 2.0 * mu_planet / rp_min)
    vp_target = sqrt(vinf_target * vinf_target + 2.0 * mu_planet / rp_min)
    return 2.0 * abs(vp_in - vp_target)


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


def bend_decompose(
    vin_vec: Vec3,
    vout_vec: Vec3,
    orbit_normal: Vec3,
) -> tuple[float, float]:
    """Attribute a flyby :math:`V_\\infty` bend to in/out-of-plane components.

    **Diagnostic only** (M-3D Phase 3, Approval Q4): this splits the
    :math:`V_\\infty`-in → :math:`V_\\infty`-out rotation into the part *in* the
    orbit plane (which changes :math:`a, e`) and the part *out of* the orbit
    plane (which changes :math:`i` and the node), so a reviewer can see how much
    plane-change work a gravity assist is doing. No optimiser degree of freedom
    is introduced; the cost model (:func:`flyby_dv`) is unchanged.

    The in-plane component is the unsigned angle between the projections of the
    two vectors onto the plane orthogonal to ``orbit_normal``. The out-of-plane
    component is the change in elevation above that plane: the difference of each
    vector's signed ``arcsin(v . n_hat / |v|)``.

    Parameters
    ----------
    vin_vec, vout_vec:
        Length-3 heliocentric :math:`V_\\infty` vectors, km/s. Non-zero.
    orbit_normal:
        Length-3 normal of the reference orbit plane (need not be unit; it is
        normalised internally). The ecliptic normal ``(0, 0, 1)`` recovers the
        in/out-of-ecliptic split.

    Returns
    -------
    (delta_inplane_rad, delta_outofplane_rad):
        In-plane bend (unsigned, ``[0, pi]``) and out-of-plane bend (signed).

    Raises
    ------
    ValueError
        If either :math:`V_\\infty` vector or ``orbit_normal`` is zero-length.
    """
    vin_norm = float(np.linalg.norm(vin_vec))
    vout_norm = float(np.linalg.norm(vout_vec))
    n_norm = float(np.linalg.norm(orbit_normal))
    if vin_norm == 0.0 or vout_norm == 0.0:
        raise ValueError("bend_decompose requires non-zero V_inf vectors")
    if n_norm == 0.0:
        raise ValueError("bend_decompose requires a non-zero orbit_normal")

    n_hat = np.asarray(orbit_normal, dtype=np.float64) / n_norm

    # Out-of-plane: signed elevation of each vector above the plane.
    sin_in = float(np.dot(vin_vec, n_hat)) / vin_norm
    sin_out = float(np.dot(vout_vec, n_hat)) / vout_norm
    elev_in = float(asin(min(1.0, max(-1.0, sin_in))))
    elev_out = float(asin(min(1.0, max(-1.0, sin_out))))
    delta_outofplane = elev_out - elev_in

    # In-plane: angle between the projections onto the plane orthogonal to n_hat.
    proj_in = np.asarray(vin_vec, dtype=np.float64) - float(np.dot(vin_vec, n_hat)) * n_hat
    proj_out = np.asarray(vout_vec, dtype=np.float64) - float(np.dot(vout_vec, n_hat)) * n_hat
    proj_in_norm = float(np.linalg.norm(proj_in))
    proj_out_norm = float(np.linalg.norm(proj_out))
    if proj_in_norm == 0.0 or proj_out_norm == 0.0:
        # A vector lying along the normal has no in-plane direction; the bend is
        # purely out-of-plane.
        delta_inplane = 0.0
    else:
        cos_arg = float(np.dot(proj_in, proj_out)) / (proj_in_norm * proj_out_norm)
        delta_inplane = float(np.arccos(np.clip(cos_arg, -1.0, 1.0)))

    return delta_inplane, delta_outofplane


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
