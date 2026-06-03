"""Resonance-anchored construction of a cycler from its sourced orbit.

This module implements the *construction* half of the McConaghy/Russell cycler
discovery method (McConaghy et al. 2002/2006; Russell & Ocampo 2004; Hollister
& Menning 1970): rather than free-optimising legs to hit V_inf targets, it
builds the spacecraft's heliocentric ellipse from its *sourced* ``(a, e)``,
places the planetary encounters at the orbit's radial crossings, and computes
each encounter's V_inf analytically in the circular-coplanar model. The result
is a family-correct seed (and, on its own, a golden cross-check).

Provenance / golden discipline
------------------------------
The orbit elements ``(a, e)`` are an *input* — they are the sourced quantity
(e.g. S1L1: a=1.30 AU, e=0.257). Everything this module returns (crossing true
anomalies, leg ToFs, per-body V_inf) is COMPUTED from those elements plus the
planets' circular-coplanar velocities. We never assert against the ``(a, e)``
we were given; the tests cross-check the *computed* V_inf against *independently
sourced* V_inf anchors.

For S1L1 the coplanar construction reproduces the independently published
COPLANAR V_inf anchors: Russell & Ocampo 2004 tabulate 4.99 / 5.10 km/s
(Earth / Mars) for the McConaghy "Notable Two-Synodic", and McConaghy 2006's
abstract gives 4.7 / 5.0 km/s. These two pairs are separate publications and
agree only in the circular-coplanar model, so matching them from ``(a, e)`` is
a non-circular validation.

Note on the spec §9 values (5.65 / 3.05 km/s): those are a *higher-fidelity*
(real-ephemeris) figure, not the coplanar construction's output. In particular
the Mars 3.05 km/s requires an *eccentric* Mars orbit (Mars e=0.093), which the
circular-coplanar model here deliberately does not carry. Reproducing 5.65/3.05
is the real-ephemeris domain (Task 3), not this module's.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, atan2, pi, sqrt, tan

import numpy as np

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.kepler import coe_to_rv


@dataclass(frozen=True)
class ResonantConstruction:
    """Result of :func:`construct_resonant_cycler`.

    All fields are COMPUTED from the input orbit elements and the planets'
    circular-coplanar velocities.

    Attributes
    ----------
    a_au:
        Spacecraft heliocentric semi-major axis used, AU (echo of the input).
    e:
        Spacecraft eccentricity used (echo of the input).
    vinf_kms:
        Hyperbolic-excess speed at each body's crossing, km/s, keyed by body
        code. ``|v_sc - v_planet|`` with the planet on a circular coplanar
        orbit at its ``sma_au``.
    crossing_true_anom_rad:
        True anomaly (radians, outbound branch in ``[0, pi]``) at which the
        spacecraft orbit crosses each body's heliocentric radius, keyed by
        body code.
    leg_tofs_days:
        Keplerian time-of-flight (days) along the spacecraft ellipse between
        consecutive crossings, in the order the ``bodies`` were supplied (the
        last entry wraps from the final body's crossing back to the first via
        apoapsis). Computed from the eccentric/mean anomaly.
    """

    a_au: float
    e: float
    vinf_kms: dict[str, float]
    crossing_true_anom_rad: dict[str, float]
    leg_tofs_days: dict[str, float]


def _true_to_mean_anomaly(nu: float, e: float) -> float:
    """True anomaly -> mean anomaly for an elliptic orbit (radians).

    Standard relations (Vallado §2.2): eccentric anomaly ``E`` from
    ``tan(E/2) = sqrt((1-e)/(1+e)) tan(nu/2)``, then Kepler's equation
    ``M = E - e sin E``. Returns ``M`` in ``[0, 2pi)``.
    """
    ecc_anom = 2.0 * atan2(
        sqrt(1.0 - e) * tan(nu / 2.0),
        sqrt(1.0 + e),
    )
    mean_anom = ecc_anom - e * np.sin(ecc_anom)
    return float(mean_anom % (2.0 * pi))


def construct_resonant_cycler(
    a_au: float,
    e: float,
    bodies: tuple[str, ...] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
) -> ResonantConstruction:
    """Construct a cycler's encounter schedule + V_inf from a sourced orbit.

    Builds the heliocentric ellipse from ``(a_au, e)``, finds the outbound
    true anomaly where the orbit crosses each body's circular radius, and
    computes the per-body V_inf and the inter-crossing Keplerian ToFs. This is
    the circular-coplanar resonance-anchored construction; planets are treated
    as being on circular coplanar orbits at their ``sma_au``.

    Parameters
    ----------
    a_au:
        Spacecraft heliocentric semi-major axis, AU (sourced input).
    e:
        Spacecraft eccentricity, ``0 <= e < 1`` (sourced input).
    bodies:
        Body codes (keys into :data:`cyclerfinder.core.constants.PLANETS`)
        the orbit encounters, ordered along the trajectory. Default
        ``("E", "M")``.
    mu:
        Heliocentric gravitational parameter, km^3/s^2.

    Returns
    -------
    ResonantConstruction
        Frozen dataclass of COMPUTED quantities (see its docstring).

    Raises
    ------
    ValueError
        If the orbit does not reach a requested body (``|cos nu| > 1``), i.e.
        the body's radius lies outside ``[a(1-e), a(1+e)]``.
    """
    a_km = a_au * AU_KM
    p = a_km * (1.0 - e * e)

    vinf_kms: dict[str, float] = {}
    crossing_true_anom_rad: dict[str, float] = {}
    period_s = 2.0 * pi * sqrt(a_km**3 / mu)

    for body in bodies:
        r_body = PLANETS[body].sma_au * AU_KM
        # Radius equation r = p / (1 + e cos nu) solved for cos nu.
        cos_nu = (p / r_body - 1.0) / e if e != 0.0 else 0.0
        if abs(cos_nu) > 1.0:
            raise ValueError(
                f"orbit (a={a_au} AU, e={e}) does not reach body {body!r} "
                f"(r={r_body / AU_KM:.4f} AU); cos(nu)={cos_nu:.4f} out of range "
                f"[a(1-e), a(1+e)] = [{a_au * (1 - e):.4f}, {a_au * (1 + e):.4f}] AU"
            )
        # Outbound branch: nu in [0, pi].
        nu = acos(cos_nu)
        crossing_true_anom_rad[body] = nu

        # Spacecraft inertial state at the crossing.
        r_sc, v_sc = coe_to_rv(a_km=a_km, e=e, true_anom_rad=nu, mu=mu)

        # Planet on a circular coplanar orbit at the same heliocentric radius:
        # speed sqrt(mu / r_body), direction purely transverse (perpendicular
        # to r_sc, prograde). Build the transverse unit vector by rotating the
        # radial unit vector +90 deg about +z.
        r_hat = r_sc / np.linalg.norm(r_sc)
        t_hat = np.array([-r_hat[1], r_hat[0], 0.0], dtype=np.float64)
        v_planet = sqrt(mu / r_body) * t_hat

        vinf_kms[body] = float(np.linalg.norm(v_sc - v_planet))

    # Leg ToFs along the ellipse between consecutive crossings, in body order,
    # the last wrapping back to the first. Time from periapsis to a crossing is
    # M / n; the leg time is the forward (positive, mod period) difference.
    n = 2.0 * pi / period_s
    mean_anoms = {b: _true_to_mean_anomaly(crossing_true_anom_rad[b], e) for b in bodies}
    leg_tofs_days: dict[str, float] = {}
    n_bodies = len(bodies)
    for i, body in enumerate(bodies):
        nxt = bodies[(i + 1) % n_bodies]
        dm = (mean_anoms[nxt] - mean_anoms[body]) % (2.0 * pi)
        leg_seconds = dm / n
        leg_tofs_days[f"{body}->{nxt}"] = leg_seconds / SECONDS_PER_DAY

    return ResonantConstruction(
        a_au=a_au,
        e=e,
        vinf_kms=vinf_kms,
        crossing_true_anom_rad=crossing_true_anom_rad,
        leg_tofs_days=leg_tofs_days,
    )
