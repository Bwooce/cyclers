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

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris

Vec3 = NDArray[np.float64]


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
