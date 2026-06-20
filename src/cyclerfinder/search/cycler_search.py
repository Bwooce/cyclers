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
from math import cos, inf, radians, sin

import numpy as np

from cyclerfinder.core.flyby import max_bend
from cyclerfinder.search.generic_return import GenericReturn, RussellModel

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
