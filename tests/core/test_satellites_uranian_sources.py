"""#894: the Uranian rows of the satellite registry against the papers they come from.

Sources (both filed in the private paper corpus, read 2026-10-04):
  * Jacobson (2014), "The orbits of the Uranian satellites and rings, the gravity field of the
    Uranian system, and the orientation of the pole of Uranus", AJ 148:76, Table 12, column
    "Current Results" (the URA111 solution): system GM and satellite GMs.
  * Jacobson & Park (2025), "The Orbits of Uranus, Its Satellites and Rings, the Gravity Field of
    the Uranian System, and the Orientation of the Poles of Uranus and Its Satellites", AJ 169:65,
    Table 6 ("Satellite Equatorial Geometric Orbital Elements", the URA182 solution): semi-major
    axes and periods; Table 2: that solution's own GMs.

Until this date the registry's values carried no paper citation and nobody had checked them.
The tests pin each value to the printed table, record that the masses and the distances come
from two different solutions, and hold open (as a strict expected failure) the fact that the
registry's mean motions are not the moons' real ones.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.satellites import (
    PRIMARIES,
    SATELLITES,
    URANIAN_PRINTED_PERIODS_DAYS,
)

# Jacobson (2014) Table 12, "Current Results", km^3/s^2.
_JACOBSON_2014_GM = {
    "Miranda": 4.3,
    "Ariel": 83.5,
    "Umbriel": 85.1,
    "Titania": 226.9,
    "Oberon": 205.3,
}
_JACOBSON_2014_SYSTEM_GM = 5794556.4

# Jacobson & Park (2025) Table 6: geometric semi-major axis (km) and period (days).
_JACOBSON_PARK_2025_A_KM = {
    "Miranda": 129846.0,
    "Ariel": 190929.0,
    "Umbriel": 265986.0,
    "Titania": 436298.0,
    "Oberon": 583511.0,
}
_JACOBSON_PARK_2025_PERIOD_DAYS = {
    "Miranda": 1.413479,
    "Ariel": 2.520379,
    "Umbriel": 4.144177,
    "Titania": 8.705869,
    "Oberon": 13.463237,
}
# Jacobson & Park (2025) Table 2, "Current" column, km^3/s^2: NOT what the registry carries.
_JACOBSON_PARK_2025_GM = {
    "Miranda": 4.11,
    "Ariel": 83.43,
    "Umbriel": 85.40,
    "Titania": 222.80,
    "Oberon": 214.21,
}

_MOONS = sorted(_JACOBSON_2014_GM)


def test_system_gm_is_jacobson_2014() -> None:
    assert PRIMARIES["Uranus"] == pytest.approx(_JACOBSON_2014_SYSTEM_GM, rel=1e-12)


@pytest.mark.parametrize("moon", _MOONS)
def test_moon_gm_is_jacobson_2014_table_12(moon: str) -> None:
    assert SATELLITES[moon].mu_km3_s2 == pytest.approx(_JACOBSON_2014_GM[moon], rel=1e-12)


@pytest.mark.parametrize("moon", _MOONS)
def test_moon_semi_major_axis_is_jacobson_park_2025_table_6(moon: str) -> None:
    assert SATELLITES[moon].sma_km == pytest.approx(_JACOBSON_PARK_2025_A_KM[moon], rel=1e-12)


def test_masses_and_distances_come_from_two_different_solutions() -> None:
    """Recorded, not hidden: the registry's GMs are the 2014 solution's and differ from the
    2025 solution's by up to 4.3 percent (Oberon), more than either paper's stated uncertainty
    for Oberon's 2025 value (2.14)."""
    assert SATELLITES["Titania"].mu_km3_s2 / _JACOBSON_PARK_2025_GM["Titania"] - 1.0 == (
        pytest.approx(0.0184, abs=5e-4)
    )
    assert SATELLITES["Oberon"].mu_km3_s2 / _JACOBSON_PARK_2025_GM["Oberon"] - 1.0 == (
        pytest.approx(-0.0416, abs=5e-4)
    )


@pytest.mark.parametrize("moon", _MOONS)
def test_printed_periods_are_transcribed(moon: str) -> None:
    assert URANIAN_PRINTED_PERIODS_DAYS[moon] == _JACOBSON_PARK_2025_PERIOD_DAYS[moon]


def test_printed_periods_agree_with_jacobson_2014_mean_longitude_rates() -> None:
    """Jacobson (2014) Table 2 prints mean longitude rates of 41.3514187 deg/day (Titania) and
    26.7394835 deg/day (Oberon); the two papers agree on the periods to the printed digits."""
    assert pytest.approx(URANIAN_PRINTED_PERIODS_DAYS["Titania"], abs=2e-6) == 360.0 / 41.3514187
    assert pytest.approx(URANIAN_PRINTED_PERIODS_DAYS["Oberon"], abs=2e-6) == 360.0 / 26.7394835


@pytest.mark.parametrize("moon", ["Titania", "Oberon"])
def test_keplers_law_on_the_geometric_axis_does_not_give_the_period(moon: str) -> None:
    """The trap, pinned: Kepler's third law with the system GM and the geometric semi-major
    axis is off by 4e-5 (Titania) and 2e-4 (Oberon) from the printed period."""
    sat = SATELLITES[moon]
    kepler_days = 2.0 * math.pi * math.sqrt(sat.sma_km**3 / PRIMARIES["Uranus"]) / 86400.0
    relative = kepler_days / URANIAN_PRINTED_PERIODS_DAYS[moon] - 1.0
    assert 3e-5 < relative < 3e-4


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#894: SatelliteData.mean_motion_deg_day is Kepler's third law on the geometric "
        "semi-major axis, not the moon's real rate. To be fixed after #895 has finished "
        "(the registry must not change under a running build)."
    ),
)
@pytest.mark.parametrize("moon", ["Titania", "Oberon"])
def test_registry_mean_motion_matches_the_printed_period(moon: str) -> None:
    printed_rate = 360.0 / URANIAN_PRINTED_PERIODS_DAYS[moon]
    assert SATELLITES[moon].mean_motion_deg_day == pytest.approx(printed_rate, rel=1e-6)


def test_registry_uranus_gm_is_the_system_value_not_the_planet() -> None:
    """Jacobson (2014) Table 12 prints two values: "System" 5794556.4 and, as quoted in Jacobson &
    Park (2025) Table 2, "Uranus" 5793951.3 +/- 4.4 km^3/s^2. The registry carries the SYSTEM
    value, and it equals the planet's GM plus the five major moons' to within the planet's stated
    uncertainty.

    This matters: a propagator that uses ``PRIMARIES["Uranus"]`` as the central mass AND adds
    the moons as separate bodies counts their mass twice (605 km^3/s^2, 1.04e-4 of the centre).
    Measured 2026-10-04 (#894): Titania and Oberon propagated that way drift 1,952 and 1,583 km
    from the URA111 kernel in 30 days, against 0.0 and 0.8 km with the planet's GM at the centre.
    The Uranian V4 lanes (``data/validation/v4_uranus.py``, ``v4_uranus_strict.py``) do exactly
    that. (For Jupiter and Saturn the registry value is the planet's own GM, so their lanes do
    not double count; the registry's comments call all of them "system GM".)
    """
    moons = sum(_JACOBSON_2014_GM.values())
    planet_gm_jacobson_2014 = 5793951.3
    assert PRIMARIES["Uranus"] - moons == pytest.approx(planet_gm_jacobson_2014, abs=4.4)
    assert moons / PRIMARIES["Uranus"] == pytest.approx(1.044e-4, rel=1e-2)
