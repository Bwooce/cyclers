r"""Jones AAS 17-577 B-plane flyby-targeting kernel (the GMAT V4 #171 piece).

Pure geometry from a flyby's incoming/outgoing hyperbolic-excess pair
:math:`(\mathbf{v}_\infty^-, \mathbf{v}_\infty^+)` to the body-centered-equatorial
B-plane goal a GMAT ``Target``/``Achieve`` block aims at:
``(BdotR, BdotT, r_p, theta_B)``. The math is Jones, Hernandez & Jesick,
"Low Excess Speed Triple Cyclers of Venus, Earth, and Mars," AAS 17-577 (2017),
Eqs.1-5, transcribed in ``docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md``
§2.2:

* :math:`\hat S = \mathbf{v}_\infty^-/|\mathbf{v}_\infty^-|` (incoming asymptote);
* :math:`\hat T = (\hat S \times \hat k)/\|\hat S \times \hat k\|`,
  :math:`\hat k = (0,0,1)` (the body pole) — Eq.4;
* :math:`\hat R = \hat S \times \hat T` — Eq.4;
* turn :math:`\delta = \angle(\mathbf{v}_\infty^-, \mathbf{v}_\infty^+)` — Eq.1;
* periapsis radius :math:`r_p` solved from Eq.2
  :math:`\arcsin\!\frac{\mu}{\mu + r_p v_\infty^{-2}}
  + \arcsin\!\frac{\mu}{\mu + r_p v_\infty^{+2}} = \delta`;
* B-plane angle :math:`\theta_B
  = \operatorname{atan2}(\hat v_\infty^+\!\cdot\hat R,\ \hat v_\infty^+\!\cdot\hat T) - \pi`
  — Eq.5;
* impact-parameter magnitude
  :math:`|B| = r_p\sqrt{1 + 2\mu/(r_p v_\infty^{-2})}`;
* B-vector goal :math:`B\!\cdot\!R = |B|\sin\theta_B`,
  :math:`B\!\cdot\!T = |B|\cos\theta_B`.

Honesty / golden discipline. Jones tabulates **no** worked
:math:`(\mathbf{v}_\infty^-, \mathbf{v}_\infty^+) \to (\theta_B, r_p, BdotR, BdotT)`
example (deep-dive §7), so the unit tests are **self-consistency round-trips**
(the goal B-vector reproduces the intended turn), NOT goldens. The only sourced
side is the row's published :math:`v_\infty` nodes (fed in by the caller). This
module REUSES :mod:`cyclerfinder.core.flyby` (``max_bend``, ``bend_angle``,
registry-resolved ``mu``/``rp_min`` via :data:`PLANETS`/:data:`SAFE_PERIHELION_KM`)
rather than re-literal-ing constants or re-deriving the bend cone.

The deliverable serves the manual, out-of-CI GMAT V4 generator
(``scripts/gmat_v4_*.py``); see ``docs/notes/2026-06-08-gmat-v4-design.md`` §3.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan2, cos, pi, sin, sqrt

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.flyby import bend_angle, max_bend

Vec3 = NDArray[np.float64]

_K_HAT: Vec3 = np.array([0.0, 0.0, 1.0], dtype=np.float64)
"""Body pole :math:`\\hat k = (0,0,1)` (Jones Eq.4, body-centered equatorial)."""


@dataclass(frozen=True)
class BPlaneTarget:
    """The GMAT B-plane goal computed from a flyby :math:`(v_\\infty^-, v_\\infty^+)`.

    Attributes
    ----------
    s_hat, t_hat, r_hat:
        The orthonormal B-plane frame (Jones Eq.4), each a unit length-3 vector.
    turn_rad:
        The required turn angle :math:`\\delta = \\angle(v_\\infty^-, v_\\infty^+)`,
        radians (Jones Eq.1).
    rp_km:
        Periapsis radius solving Jones Eq.2, km.
    theta_b_rad:
        The B-plane angle :math:`\\theta_B` (Jones Eq.5), radians in :math:`(-\\pi,\\pi]`.
    b_mag_km:
        Impact-parameter magnitude :math:`|B|`, km.
    bdot_r_km, bdot_t_km:
        The B-vector goal components :math:`B\\cdot R`, :math:`B\\cdot T`, km — what
        the GMAT ``Achieve`` lines target.
    feasible:
        ``True`` iff the required turn is within the ballistic bend cone at the safe
        periapsis (``turn_rad <= max_bend``). When ``False`` the geometry needs a
        powered maneuver / TCM; the generator either rejects the node or adds the
        impulse (Phase 0 Task 0.3).
    """

    s_hat: Vec3
    t_hat: Vec3
    r_hat: Vec3
    turn_rad: float
    rp_km: float
    theta_b_rad: float
    b_mag_km: float
    bdot_r_km: float
    bdot_t_km: float
    feasible: bool


def bplane_frame(vinf_minus: Vec3) -> tuple[Vec3, Vec3, Vec3]:
    r"""Body-centered-equatorial B-plane frame :math:`(\hat S, \hat T, \hat R)`.

    Jones Eq.4: :math:`\hat S = \hat v_\infty^-`,
    :math:`\hat T = (\hat S \times \hat k)/\|\cdot\|`, :math:`\hat R = \hat S \times \hat T`,
    with :math:`\hat k = (0,0,1)` the body pole. The three returned vectors are
    unit-length and mutually orthogonal (a right-handed triad).

    Parameters
    ----------
    vinf_minus:
        Incoming hyperbolic-excess vector :math:`v_\infty^-`, length-3, km/s.
        Must be non-zero and not (anti)parallel to the pole :math:`\hat k`.

    Returns
    -------
    (s_hat, t_hat, r_hat):
        The orthonormal B-plane frame.

    Raises
    ------
    ValueError
        If ``vinf_minus`` is zero-length, or parallel to the pole (``S x k`` is
        degenerate — the equatorial B-plane frame is undefined there).
    """
    s = np.asarray(vinf_minus, dtype=np.float64)
    s_norm = float(np.linalg.norm(s))
    if s_norm == 0.0:
        raise ValueError("bplane_frame requires a non-zero v_inf_minus")
    s_hat = s / s_norm
    s_cross_k = np.cross(s_hat, _K_HAT)
    cross_norm = float(np.linalg.norm(s_cross_k))
    if cross_norm < 1.0e-12:
        raise ValueError(
            "v_inf_minus is (anti)parallel to the pole k=(0,0,1); the equatorial "
            "B-plane frame T=(S x k)/|S x k| is undefined"
        )
    t_hat = s_cross_k / cross_norm
    r_hat = np.cross(s_hat, t_hat)
    return s_hat, t_hat, r_hat


def _solve_rp(turn_rad: float, vinf_minus_mag: float, vinf_plus_mag: float, mu: float) -> float:
    r"""Solve Jones Eq.2 for periapsis radius :math:`r_p` (km).

    :math:`\arcsin\frac{\mu}{\mu + r_p v_-^2} + \arcsin\frac{\mu}{\mu + r_p v_+^2} = \delta`.

    The left side is strictly decreasing in :math:`r_p` (each term decreases), from
    :math:`\pi` at :math:`r_p=0` to :math:`0` as :math:`r_p\to\infty`, so for any
    :math:`\delta \in (0, \pi)` a unique root exists — found by bisection. (Jones
    allows subsurface roots here; the min-altitude filter / feasibility flag handles
    them downstream.)
    """

    def lhs(rp: float) -> float:
        a1 = asin(min(1.0, mu / (mu + rp * vinf_minus_mag * vinf_minus_mag)))
        a2 = asin(min(1.0, mu / (mu + rp * vinf_plus_mag * vinf_plus_mag)))
        return a1 + a2

    lo, hi = 0.0, 1.0
    # Expand hi until lhs(hi) < turn (lhs is decreasing in rp).
    while lhs(hi) > turn_rad:
        hi *= 2.0
        if hi > 1.0e15:  # pragma: no cover — defensive ceiling
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if lhs(mid) > turn_rad:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def bplane_target(
    vinf_minus: Vec3,
    vinf_plus: Vec3,
    mu: float,
    rp_min: float,
) -> BPlaneTarget:
    r"""Compute the GMAT B-plane goal from a flyby :math:`(v_\infty^-, v_\infty^+)`.

    Implements Jones Eqs.1,2,4,5 plus :math:`|B| = r_p\sqrt{1 + 2\mu/(r_p v_-^2)}`
    and :math:`B\cdot R = |B|\sin\theta_B`, :math:`B\cdot T = |B|\cos\theta_B`.

    Parameters
    ----------
    vinf_minus, vinf_plus:
        Incoming / outgoing hyperbolic-excess vectors, length-3, km/s, body-centered
        equatorial. ``vinf_minus`` fixes the frame; the turn is the angle between them.
    mu:
        Flyby-body gravitational parameter ``GM``, km**3/s**2.
    rp_min:
        Minimum safe periapsis radius, km (used only for the feasibility flag; the
        targeted ``rp_km`` is the Eq.2 root for the requested turn).

    Returns
    -------
    BPlaneTarget
        The full goal (frame, turn, ``r_p``, ``theta_B``, ``|B|``, ``BdotR``/``BdotT``,
        feasibility).

    Raises
    ------
    ValueError
        If either vector is zero-length / the frame is degenerate (see
        :func:`bplane_frame`), or ``mu``/``rp_min`` non-positive.
    """
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    if rp_min <= 0.0:
        raise ValueError(f"rp_min must be positive, got {rp_min}")

    s_hat, t_hat, r_hat = bplane_frame(vinf_minus)
    vinf_minus_mag = float(np.linalg.norm(vinf_minus))
    vinf_plus_mag = float(np.linalg.norm(vinf_plus))
    if vinf_plus_mag == 0.0:
        raise ValueError("bplane_target requires a non-zero v_inf_plus")

    turn = bend_angle(vinf_minus, vinf_plus)  # Jones Eq.1 (reuses core/flyby)
    rp = _solve_rp(turn, vinf_minus_mag, vinf_plus_mag, mu)

    # Jones Eq.5: theta_B = atan2(v+ . R, v+ . T) - pi, wrapped to (-pi, pi].
    vplus_hat = np.asarray(vinf_plus, dtype=np.float64) / vinf_plus_mag
    theta_b = atan2(float(np.dot(vplus_hat, r_hat)), float(np.dot(vplus_hat, t_hat))) - pi
    if theta_b <= -pi:
        theta_b += 2.0 * pi

    # |B| = r_p sqrt(1 + 2 mu / (r_p v_-^2)) — impact parameter from periapsis radius.
    b_mag = rp * sqrt(1.0 + 2.0 * mu / (rp * vinf_minus_mag * vinf_minus_mag))
    bdot_r = b_mag * sin(theta_b)
    bdot_t = b_mag * cos(theta_b)

    # Feasibility: required turn within the ballistic cone at the SAFE periapsis,
    # for the (worst-case) smaller |v_inf|. Reuses core/flyby.max_bend.
    v_for_bend = min(vinf_minus_mag, vinf_plus_mag)
    feasible = turn <= max_bend(mu, rp_min, v_for_bend)

    return BPlaneTarget(
        s_hat=s_hat,
        t_hat=t_hat,
        r_hat=r_hat,
        turn_rad=turn,
        rp_km=rp,
        theta_b_rad=theta_b,
        b_mag_km=b_mag,
        bdot_r_km=bdot_r,
        bdot_t_km=bdot_t,
        feasible=feasible,
    )


def bplane_target_for(code: str, vinf_minus: Vec3, vinf_plus: Vec3) -> BPlaneTarget:
    """Planet-aware wrapper around :func:`bplane_target`.

    Resolves ``mu`` from :data:`cyclerfinder.core.constants.PLANETS` and ``rp_min``
    from :data:`SAFE_PERIHELION_KM` for the one-letter planet ``code`` — the same
    registry touch-point :func:`cyclerfinder.core.flyby.flyby_dv_for` uses.

    Raises
    ------
    KeyError
        If ``code`` is not in :data:`PLANETS`.
    """
    mu = PLANETS[code].mu_km3_s2
    rp = SAFE_PERIHELION_KM[code]
    return bplane_target(vinf_minus, vinf_plus, mu, rp)


__all__ = [
    "BPlaneTarget",
    "bplane_frame",
    "bplane_target",
    "bplane_target_for",
]
