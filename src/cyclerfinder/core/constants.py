"""Physical constants used across cyclerfinder.

This module is the single source of truth for the numerical values consumed
by the M1+ mechanics code (ephemeris, Lambert, Kepler, flyby). Every value
is cited inline so a future reader can audit it against the original source.

Sources
-------

Universal:
    * ``MU_SUN_KM3_S2`` — IAU 2015 nominal solar mass parameter
      ``GM_S^N = 1.3271244e20 m^3/s^2`` (IAU 2015 Resolution B3), converted to
      km^3/s^2 with the extra digits from JPL DE440 (1.32712440018e11).
    * ``AU_KM`` — IAU 2012 Resolution B2, exact definition
      ``1 au = 149_597_870_700 m``.
    * ``SECONDS_PER_DAY`` — SI definition.
    * ``DAYS_PER_JULIAN_YEAR`` — 365.25 d, matches ``astropy.units.yr``
      (Julian year, the convention used throughout the cycler literature).

Planets:
    * ``mu_km3_s2`` — gravitational parameter ``GM`` of the planet.
        - Venus, Mars: JPL DE440 (Park et al., AJ 2021, table 5).
        - Earth: IERS 2010 conventions / EGM2008 (``GM_E = 3.986004418e14
          m^3/s^2``); the trailing digits ``...507`` track JPL's DE440 value
          for the Earth+Moon vs Earth-only choice. We use the Earth-only μ
          (no Moon) since the patched-conic flyby treats only the central
          body. Source: JPL DE440 table 5.
    * ``radius_eq_km`` — mean equatorial radius.
        - Venus, Mars: IAU 2015 Working Group on Cartographic Coordinates
          and Rotational Elements (Archinal et al., CMDA 2018).
        - Earth: WGS84 semi-major axis, 6378.137 km exactly.
    * ``sma_au`` — heliocentric semi-major axis at J2000, from Standish &
      Williams, "Approximate Positions of the Planets" (JPL Solar System
      Dynamics, 1992/2006 update), table 1, ``a_0`` column.
    * ``mean_motion_deg_day`` — 360 / orbital period in days, derived from
      ``sma_au`` and ``MU_SUN_KM3_S2`` via Kepler's third law and rounded
      to six decimals. (Recomputed at import time so the table is internally
      consistent rather than copy-pasted.)
    * ``safe_alt_km`` — conservative minimum flyby altitude (atmosphere top
      plus margin). 300 km for all three; Aldrin's original work used 200 km
      and later phases may override per body via config.

Out of scope here (deliberately, per ``docs/phases/m0-scaffold/plan.md``
§4.4): Sun radius, planet eccentricity/inclination (the circular-coplanar
model treats both as zero; the JPL ephemeris backend in M6 will supply
them), moons, asteroids.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Final

# ---------------------------------------------------------------------------
# Universal constants
# ---------------------------------------------------------------------------

MU_SUN_KM3_S2: Final[float] = 1.32712440018e11
"""Heliocentric gravitational parameter, km^3/s^2. JPL DE440 / IAU 2015."""

AU_KM: Final[float] = 1.49597870700e8
"""Astronomical unit, km. IAU 2012 Resolution B2 (exact)."""

SECONDS_PER_DAY: Final[float] = 86400.0
"""Seconds in a day (SI, exact)."""

DAYS_PER_JULIAN_YEAR: Final[float] = 365.25
"""Days in a Julian year (exact convention, matches ``astropy.units.yr``)."""


# ---------------------------------------------------------------------------
# Planet record
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PlanetData:
    """Per-planet physical constants needed by the M1+ mechanics code.

    Attributes
    ----------
    name:
        Full English name, e.g. ``"Venus"``.
    code:
        One-letter short code used as the key in :data:`PLANETS` and in
        sequence strings (e.g. ``"VEM"``).
    mu_km3_s2:
        Gravitational parameter ``GM`` of the planet, km^3/s^2.
    radius_eq_km:
        Mean equatorial radius, km.
    sma_au:
        Mean heliocentric semi-major axis at J2000, AU.
    mean_motion_deg_day:
        Mean motion ``n = 360°/P``, deg/day, with ``P`` derived from
        ``sma_au`` and the heliocentric ``GM_Sun``.
    safe_alt_km:
        Minimum allowable flyby altitude above the equatorial radius, km.
        Conservative default; later phases may override per body via config.
    """

    name: str
    code: str
    mu_km3_s2: float
    radius_eq_km: float
    sma_au: float
    mean_motion_deg_day: float
    safe_alt_km: float


def _mean_motion_deg_day(sma_au: float) -> float:
    """Derive mean motion (deg/day) from the heliocentric semi-major axis.

    Uses Kepler's third law with ``MU_SUN_KM3_S2`` so all three planets'
    mean motions are internally consistent with the AU and μ_Sun adopted
    above (rather than hand-copied from disparate sources).
    """
    a_km = sma_au * AU_KM
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
    period_days = period_s / SECONDS_PER_DAY
    return 360.0 / period_days


# ---------------------------------------------------------------------------
# Planet table — Venus, Earth, Mars (others added as later phases need them)
# ---------------------------------------------------------------------------

_VENUS_SMA_AU: Final[float] = 0.72333566
_EARTH_SMA_AU: Final[float] = 1.00000261
_MARS_SMA_AU: Final[float] = 1.52371034

PLANETS: Final[dict[str, PlanetData]] = {
    "V": PlanetData(
        name="Venus",
        code="V",
        mu_km3_s2=3.24858592e5,
        radius_eq_km=6051.8,
        sma_au=_VENUS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_VENUS_SMA_AU),
        safe_alt_km=300.0,
    ),
    "E": PlanetData(
        name="Earth",
        code="E",
        mu_km3_s2=3.98600435507e5,
        radius_eq_km=6378.137,
        sma_au=_EARTH_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_EARTH_SMA_AU),
        safe_alt_km=300.0,
    ),
    "M": PlanetData(
        name="Mars",
        code="M",
        mu_km3_s2=4.282837521e4,
        radius_eq_km=3396.19,
        sma_au=_MARS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_MARS_SMA_AU),
        safe_alt_km=300.0,
    ),
}


# ---------------------------------------------------------------------------
# Derived helpers
# ---------------------------------------------------------------------------

SAFE_PERIHELION_KM: Final[dict[str, float]] = {
    code: data.radius_eq_km + data.safe_alt_km for code, data in PLANETS.items()
}
"""Minimum safe flyby periapsis radius (km) per planet code.

Lives here (rather than in the future ``flyby.py``) so the constants module
remains the single source of truth for physical numbers.
"""
