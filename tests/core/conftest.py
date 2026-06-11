"""Shared fixtures for the M1 ``tests/core`` suite.

The three canonical Lambert test legs (Aldrin medium, short Earth-to-Earth,
long Earth-to-Mars) are built once here using
:class:`cyclerfinder.core.ephemeris.Ephemeris` in its circular backend so the
individual gate tests don't repeat the boilerplate.

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §4.1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris

Vec3 = NDArray[np.float64]


def coe3d_to_rv(
    a_km: float,
    e: float,
    raan_rad: float,
    inc_rad: float,
    argp_rad: float,
    nu_rad: float,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[Vec3, Vec3]:
    """Full classical elements -> inertial state (perifocal + 3-1-3 rotation).

    Test-side helper (closed-form, independent of the Lambert solver) used by
    the task #205 pinning tests in ``test_lambert.py`` /
    ``test_lambert_multirev.py`` to construct exact boundary states on an
    inclined generating orbit. The production :func:`coe_to_rv` is planar-only.
    """
    from math import cos, sin, sqrt

    p = a_km * (1.0 - e * e)
    r_mag = p / (1.0 + e * cos(nu_rad))
    r_pf = np.array([r_mag * cos(nu_rad), r_mag * sin(nu_rad), 0.0], dtype=np.float64)
    sqrt_mu_p = sqrt(mu / p)
    v_pf = np.array(
        [-sqrt_mu_p * sin(nu_rad), sqrt_mu_p * (e + cos(nu_rad)), 0.0], dtype=np.float64
    )

    def _rz(t: float) -> NDArray[np.float64]:
        return np.array(
            [[cos(t), -sin(t), 0.0], [sin(t), cos(t), 0.0], [0.0, 0.0, 1.0]], dtype=np.float64
        )

    def _rx(t: float) -> NDArray[np.float64]:
        return np.array(
            [[1.0, 0.0, 0.0], [0.0, cos(t), -sin(t)], [0.0, sin(t), cos(t)]], dtype=np.float64
        )

    rot = _rz(raan_rad) @ _rx(inc_rad) @ _rz(argp_rad)
    return rot @ r_pf, rot @ v_pf


@dataclass(frozen=True)
class Leg:
    """A canonical Lambert test leg.

    Attributes
    ----------
    name:
        Short identifier for the leg (logged in failures).
    r1, r2:
        Endpoint position vectors in km, heliocentric ecliptic J2000.
    tof:
        Time of flight in seconds.
    prograde:
        Transfer sense; ``True`` for the standard short-way prograde transfer.
    """

    name: str
    r1: Vec3
    r2: Vec3
    tof: float
    prograde: bool = True


def _leg_from_bodies(
    name: str,
    body_from: str,
    body_to: str,
    t_from_days: float,
    t_to_days: float,
    prograde: bool = True,
) -> Leg:
    eph = Ephemeris(model="circular")
    r1, _ = eph.state(body_from, t_from_days * SECONDS_PER_DAY)
    r2, _ = eph.state(body_to, t_to_days * SECONDS_PER_DAY)
    return Leg(
        name=name,
        r1=r1,
        r2=r2,
        tof=(t_to_days - t_from_days) * SECONDS_PER_DAY,
        prograde=prograde,
    )


@pytest.fixture(scope="session")
def leg_aldrin() -> Leg:
    """Aldrin-style Earth -> Mars leg, ~146 d (spec §9, plan §4.1 leg A)."""
    return _leg_from_bodies("aldrin-EM-146d", "E", "M", 0.0, 146.0)


@pytest.fixture(scope="session")
def leg_short() -> Leg:
    """Earth -> Earth short arc, 50 d (plan §4.1 leg B).

    Earth advances ~49.3 deg in 50 d — small transfer angle, stresses sign
    handling and small-A behaviour in the universal-variable solver.
    """
    return _leg_from_bodies("short-EE-50d", "E", "E", 0.0, 50.0)


@pytest.fixture(scope="session")
def leg_long() -> Leg:
    """Earth -> Mars long arc, 500 d (plan §4.1 leg C).

    Mars advances ~262 deg in 500 d — large transfer angle, slow Newton
    convergence regime.
    """
    return _leg_from_bodies("long-EM-500d", "E", "M", 0.0, 500.0)


@pytest.fixture(scope="session")
def all_legs(leg_aldrin: Leg, leg_short: Leg, leg_long: Leg) -> list[Leg]:
    """All three canonical legs in a list for parametrised consistency tests."""
    return [leg_aldrin, leg_short, leg_long]
