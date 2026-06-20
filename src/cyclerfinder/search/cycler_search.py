"""Russell (2004) global cycler search: Eq 3.1 TOF and AR/TR feasibility gates.

This module layers Russell's ballistic Earth-Mars cycler search (#388) on top
of the generic-return machinery in :mod:`cyclerfinder.search.generic_return`.

A candidate cycler is parameterised by Russell's integer quadruple
``(p, h, s, i)`` (App A.1):

* ``p`` -- number of E-M synodic periods spanned,
* ``h`` -- a half-integer numerator (the ``h / 2`` term carries the half-rev
  bookkeeping of the transfer),
* ``s`` -- the period divisor,
* ``i`` -- a member index within the resulting cycler family.

Two scalar feasibility gates are reported per candidate:

* **AR** (aphelion ratio): the largest cycler aphelion (AU) divided by Mars'
  semi-major axis ``1.52`` AU. ``AR <= 1`` means the trajectory never reaches
  Mars and is rejected.
* **TR** (turn ratio): the turn the geometry demands divided by the maximum
  ballistic turn a 200 km-altitude Earth flyby can deliver. ``TR <= 1`` means
  the flyby can close the geometry without a maneuver.

Equations
---------
* **Eq 3.1** (App A.1): ``TOF = (tau * p - h / 2) / s`` years, where ``tau`` is
  the Earth-Mars synodic period :meth:`RussellModel.synodic_yr("E", "M")`.
* Aphelion is reconstructed from the departure state of a
  :class:`~cyclerfinder.search.generic_return.GenericReturn` -- never from its
  ``a_au`` field, which holds the transfer *period* (``a_sma ** 1.5``), not the
  semi-major axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, degrees, inf, pi, radians, sin, sqrt

import numpy as np

from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.cycler_assembly import group_half_years
from cyclerfinder.search.generic_return import (
    GenericReturn,
    RussellModel,
    _russell_branch,
    psi_of_vinf_vec,
)

# Mars semi-major axis (AU) used as the aphelion-ratio reference (App A.1).
MARS_SMA_AU = 1.52
# Earth equatorial radius (km) and the safe 200 km flyby altitude (App A.1).
EARTH_RADIUS_KM = 6378.0
EARTH_FLYBY_ALTITUDE_KM = 200.0


def cycler_tof(model: RussellModel, p: int, h: int, s: int) -> float:
    """Cycler time of flight in years (Russell Eq 3.1, App A.1).

    ``TOF = (tau * p - h / 2) / s`` where ``tau`` is the Earth-Mars synodic
    period in years.

    Parameters
    ----------
    model:
        The :class:`RussellModel` providing the synodic period.
    p, h, s:
        Russell's integer parameters (synodic-period count, half-integer
        numerator, and period divisor).

    Returns
    -------
    float
        Time of flight in years.
    """
    tau = model.synodic_yr("E", "M")
    return (tau * p - h / 2.0) / s


def aphelion_ratio(aphelion_au: float) -> float:
    """Aphelion ratio: cycler aphelion (AU) over Mars' semi-major axis.

    ``AR = aphelion_au / 1.52``. ``AR <= 1`` means the cycler never reaches
    Mars' orbit.
    """
    return aphelion_au / MARS_SMA_AU


def turn_ratio(max_allowable: float, omega_max: float) -> float:
    """Turn ratio: demanded turn over the maximum ballistic turn.

    ``TR = max_allowable / omega_max``. ``TR <= 1`` means the flyby can deliver
    the required bend ballistically.
    """
    return max_allowable / omega_max


def max_earth_flyby_bend(model: RussellModel, vinf: float) -> float:
    """Maximum ballistic bend (rad) of a 200 km-altitude Earth flyby.

    Reuses :func:`cyclerfinder.core.flyby.max_bend`, which is unit-agnostic
    (consistent units in, radians out): here ``mu = model.mu_earth`` (canonical
    AU**3/TU**2), ``rp = (6378 + 200) km / model.au_km`` (AU), and ``vinf`` in
    canonical AU/TU. The underlying formula is
    ``2 * asin(1 / (1 + rp * vinf**2 / mu))``.

    Parameters
    ----------
    model:
        The :class:`RussellModel` providing ``mu_earth`` and ``au_km``.
    vinf:
        Hyperbolic excess speed at the flyby, canonical AU/TU. Non-negative.

    Returns
    -------
    float
        Maximum bend angle in radians, in ``(0, pi]``.
    """
    rp_au = (EARTH_RADIUS_KM + EARTH_FLYBY_ALTITUDE_KM) / model.au_km
    return max_bend(model.mu_earth, rp_au, vinf)


def generic_return_aphelion(model: RussellModel, body: str, gret: GenericReturn) -> float:
    """Heliocentric aphelion (AU) of a generic return's transfer ellipse.

    Reconstructs the departure state from the generic return's ``psi`` and
    ``vinf`` -- ``gret.a_au`` is the transfer *period*, not the semi-major axis,
    so the aphelion must be rebuilt from the Keplerian elements of the
    departure state, never read off ``a_au`` directly.

    The departure state is built in the same in-plane basis Russell uses for
    ``psi`` (see :func:`cyclerfinder.search.generic_return.psi_of_vinf_vec`):
    ``e1 = v_hat_B``, ``e2`` the component of ``r_hat_B`` orthogonal to ``e1``.
    Then ``v0 = v_B0 + vinf * (cos psi * e1 + sin psi * e2)`` and the orbit
    elements follow from the vis-viva and angular-momentum of ``(r1, v0)``.

    Parameters
    ----------
    model:
        The :class:`RussellModel`.
    body:
        Body key understood by :meth:`RussellModel.body_state`.
    gret:
        The :class:`GenericReturn` whose transfer aphelion is wanted.

    Returns
    -------
    float
        Aphelion radius ``a_sma * (1 + e)`` in AU; ``+inf`` for a
        non-elliptical (parabolic/hyperbolic) reconstructed transfer.
    """
    r1, v_b0 = model.body_state(body, 0.0)

    e1 = v_b0 / np.linalg.norm(v_b0)
    rhat = r1 / np.linalg.norm(r1)
    e2 = rhat - np.dot(rhat, e1) * e1
    e2 = e2 / np.linalg.norm(e2)

    psi = radians(gret.psi_deg)
    vinf_vec = gret.vinf * (cos(psi) * e1 + sin(psi) * e2)
    v0 = v_b0 + vinf_vec

    mu = model.mu_sun
    r1_n = float(np.linalg.norm(r1))
    speed = float(np.linalg.norm(v0))
    denom = 2.0 / r1_n - speed * speed / mu
    if denom <= 0.0:
        # Parabolic (denom == 0) or hyperbolic (denom < 0): no bound aphelion.
        return inf
    a_sma = 1.0 / denom

    h_vec = np.cross(r1, v0)
    h_n = float(np.linalg.norm(h_vec))
    ecc = float(np.sqrt(max(0.0, 1.0 - h_n * h_n / (mu * a_sma))))
    return float(a_sma * (1.0 + ecc))


@dataclass(frozen=True)
class Cycler:
    """A ballistic Earth-Mars cycler candidate (Russell App A.1).

    Fields
    ------
    p, h, s, i:
        Russell's integer parameter quadruple identifying the family member.
    generic_return:
        The underlying :class:`GenericReturn` (departure geometry / TOF).
    turn_angles:
        Per-flyby turn angles (rad) demanded by the geometry.
    ar:
        Aphelion ratio (see :func:`aphelion_ratio`).
    tr:
        Turn ratio (see :func:`turn_ratio`).
    vinf_e, vinf_m:
        Hyperbolic excess speeds at Earth and Mars (canonical AU/TU).
    """

    p: int
    h: int
    s: int
    i: int
    generic_return: GenericReturn
    turn_angles: tuple[float, ...]
    ar: float
    tr: float
    vinf_e: float
    vinf_m: float


def generic_returns_at_tof(
    model: RussellModel,
    body: str,
    tof_body_periods: float,
    *,
    max_revs_cap: int = 15,
) -> list[GenericReturn]:
    """Generic returns from a SINGLE Lambert solve at a fixed ToF (spec App A.1).

    Unlike :func:`~cyclerfinder.search.generic_return.generate_generic_returns`,
    which sweeps a grid of times of flight, this performs ONE multi-revolution
    Lambert solve at the supplied ``tof_body_periods`` -- the ``i`` solutions of
    Russell's Fig 3.9 loop -- and wraps each Lambert solution into a
    :class:`GenericReturn` the same way the grid generator does (vinf vector from
    ``v1 - v_B0``, ``psi`` via :func:`psi_of_vinf_vec`, ``a_au = a_sma ** 1.5``,
    slow/fast label via the shared :func:`_russell_branch`).

    Parameters
    ----------
    model:
        The :class:`RussellModel`.
    body:
        Body key understood by :meth:`RussellModel.body_state`.
    tof_body_periods:
        Time of flight in body orbital periods.
    max_revs_cap:
        Maximum heliocentric revolution count passed to :func:`lambert`.

    Returns
    -------
    list[GenericReturn]
        One entry per Lambert solution; an empty list if Lambert fails to
        converge for the requested ToF.
    """
    p_tu = 2.0 * pi * model.sma_au(body) ** 1.5
    r1, v_b0 = model.body_state(body, 0.0)
    r2, _ = model.body_state(body, 2.0 * pi * tof_body_periods)
    tof_canonical = tof_body_periods * p_tu
    r1_n = float(np.linalg.norm(r1))

    try:
        sols = lambert(r1, r2, tof_canonical, mu=model.mu_sun, prograde=True, max_revs=max_revs_cap)
    except LambertError:
        return []

    returns: list[GenericReturn] = []
    for sol in sols:
        vinf_vec = sol.v1 - v_b0
        vinf = float(np.linalg.norm(vinf_vec))
        if vinf <= 0.0:
            continue
        psi_deg = degrees(psi_of_vinf_vec(vinf_vec, r1, v_b0))

        speed = float(np.linalg.norm(sol.v1))
        denom = 2.0 / r1_n - speed * speed / model.mu_sun
        if abs(denom) < 1.0e-12:
            # Parabolic transfer -- semi-major axis undefined; skip.
            continue
        a_sma = 1.0 / denom
        if a_sma <= 0.0:
            # Hyperbolic transfer -- not a Russell return.
            continue
        a_au = a_sma**1.5

        returns.append(
            GenericReturn(
                psi_deg=psi_deg,
                tof_body_periods=tof_body_periods,
                a_au=a_au,
                n_revs=int(sol.n_revs),
                branch=_russell_branch(tof_canonical, a_sma, model.mu_sun),
                vinf=vinf,
            )
        )

    return returns


def _arrival_vinf_vec(model: RussellModel, body: str, ret: GenericReturn) -> np.ndarray:
    """Incoming v_inf vector at the Earth return of a generic return (canonical).

    Reconstructs the departure state from ``ret``'s ``psi``/``vinf`` in the same
    in-plane basis :func:`psi_of_vinf_vec` uses, propagates it forward over the
    transfer time of flight, and differences the arrival velocity against the
    body's velocity at the arrival angle to recover the arrival v_inf vector.
    """
    r1, v_b0 = model.body_state(body, 0.0)

    e1 = v_b0 / np.linalg.norm(v_b0)
    rhat = r1 / np.linalg.norm(r1)
    e2 = rhat - np.dot(rhat, e1) * e1
    e2 = e2 / np.linalg.norm(e2)

    psi = radians(ret.psi_deg)
    v0 = v_b0 + ret.vinf * (cos(psi) * e1 + sin(psi) * e2)

    tof_canon = ret.tof_body_periods * 2.0 * pi * model.sma_au(body) ** 1.5
    _, v_arr = propagate(r1, v0, tof_canon, model.mu_sun)

    theta_arr = 2.0 * pi * ret.tof_body_periods
    _, v_body_arr = model.body_state(body, theta_arr)
    out: np.ndarray = v_arr - v_body_arr
    return out


def _mars_vinf(model: RussellModel, ret: GenericReturn) -> float:
    """Rough scalar Mars v_inf evidence for a generic return (canonical AU/TU).

    Evidence only: the transfer speed at ``r = 1.52`` AU (vis-viva with the
    transfer semi-major axis ``a_sma = ret.a_au ** (2/3)``) differenced against
    Mars' circular speed. Returns ``0.0`` for a degenerate (non-positive)
    semi-major axis.
    """
    a_sma = ret.a_au ** (2.0 / 3.0)
    if a_sma <= 0.0:
        return 0.0
    mu = model.mu_sun
    v_at_mars = sqrt(max(0.0, mu * (2.0 / MARS_SMA_AU - 1.0 / a_sma)))
    mars_circ = model.body_circular_speed("M")
    return abs(v_at_mars - mars_circ)


def _signed_rev(ret: GenericReturn) -> int:
    """Russell's signed revolution index ``i``: ``+N`` slow, ``-N`` fast."""
    return ret.n_revs if ret.branch == "slow" else -ret.n_revs


def search_cyclers(
    model: RussellModel,
    *,
    p_max: int,
    ar_min: float = 0.9,
    tr_min: float = 0.85,
    max_revs_cap: int = 15,
) -> list[Cycler]:
    """Russell's Fig 3.9 global ``(p, h, s, i)`` cycler search.

    Sweeps ``p in 1..p_max``, ``h in 1..5*p_max``, ``s in 1..3*p_max``; at each
    triple computes the Eq 3.1 time of flight, solves a single-ToF generic return
    (:func:`generic_returns_at_tof`), and applies the AR (aphelion-ratio) and TR
    (turn-ratio) feasibility gates in that order. Surviving candidates are
    recorded as :class:`Cycler` instances.

    Parameters
    ----------
    model:
        The :class:`RussellModel`.
    p_max:
        Upper bound on the synodic-period count ``p`` (and, scaled, on ``h`` and
        ``s``).
    ar_min:
        Aphelion-ratio threshold; candidates need ``ar > ar_min`` to pass.
    tr_min:
        Turn-ratio threshold; candidates need ``tr > tr_min`` to pass.
    max_revs_cap:
        Maximum heliocentric revolution count passed to :func:`lambert`.

    Returns
    -------
    list[Cycler]
        Every ``(p, h, s, i)`` candidate that clears both gates.
    """
    body = "E"
    period_yr = model.period_yr(body)
    cyclers: list[Cycler] = []

    for p in range(1, p_max + 1):
        for h in range(1, 5 * p_max + 1):
            for s in range(1, 3 * p_max + 1):
                tof_yr = cycler_tof(model, p, h, s)
                if tof_yr <= 0.0:
                    continue
                tof_bp = tof_yr / period_yr
                for ret in generic_returns_at_tof(model, body, tof_bp, max_revs_cap=max_revs_cap):
                    if ret.n_revs < 1:
                        # Need a real multi-rev generic return.
                        continue
                    aph = generic_return_aphelion(model, body, ret)
                    ar = aphelion_ratio(aph)
                    if ar <= ar_min:
                        continue
                    vinf_in_vec = _arrival_vinf_vec(model, body, ret)
                    _hs, omega_max = group_half_years(model, body, ret.vinf, vinf_in_vec, h, s)
                    if omega_max <= 0.0:
                        continue
                    tr = turn_ratio(max_earth_flyby_bend(model, ret.vinf), omega_max)
                    if tr <= tr_min:
                        continue
                    vinf_m = _mars_vinf(model, ret)
                    cyclers.append(
                        Cycler(
                            p=p,
                            h=h,
                            s=s,
                            i=_signed_rev(ret),
                            generic_return=ret,
                            turn_angles=(),
                            ar=ar,
                            tr=tr,
                            vinf_e=ret.vinf,
                            vinf_m=vinf_m,
                        )
                    )

    return cyclers
