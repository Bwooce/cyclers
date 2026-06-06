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
        - Mercury: JPL DE440 (Park et al., AJ 2021); Mercury has no
          satellites so its GM is unambiguously planet-only.
        - Jupiter, Saturn, Uranus, Neptune: JPL DE440 *system* GM (planet +
          satellites; Park et al., AJ 2021). Unlike Earth (where the Moon's
          mass is ~1.2% of Earth's and is deliberately excluded), each giant
          planet's satellites contribute <2e-4 of the system GM and orbit far
          above any physical flyby periapsis, so the system value is the
          correct central-body μ a patched-conic flyby senses at closest
          approach. JPL DE440 publishes only the system GM for these bodies.
    * ``radius_eq_km`` — mean equatorial radius.
        - Venus, Mars: IAU 2015 Working Group on Cartographic Coordinates
          and Rotational Elements (Archinal et al., CMDA 2018).
        - Earth: WGS84 semi-major axis, 6378.137 km exactly.
        - Mercury, Jupiter, Saturn, Uranus, Neptune: IAU 2015 WGCCRE
          (Archinal et al., CMDA 2018, 130:22) — same source as Venus/Mars.
          Giant-planet radii are the 1-bar equatorial level.
    * ``sma_au`` — heliocentric semi-major axis at J2000, from Standish &
      Williams, "Approximate Positions of the Planets" (JPL Solar System
      Dynamics, 1992/2006 update), table 1, ``a_0`` column.
    * ``mean_motion_deg_day`` — 360 / orbital period in days, derived from
      ``sma_au`` and ``MU_SUN_KM3_S2`` via Kepler's third law and rounded
      to six decimals. (Recomputed at import time so the table is internally
      consistent rather than copy-pasted.)
    * ``safe_alt_km`` — conservative minimum flyby altitude (atmosphere top
      plus margin). 300 km for all three of V/E/M; Aldrin's original work used
      200 km and later phases may override per body via config. These are
      ENGINEERING DEFAULTS (an operational convention, not sourced physics).
      The outer-body defaults are scaled for their distinct hazards and are
      likewise convention, not measured constants:
        - Mercury (1000 km): no atmosphere, but extreme solar thermal flux and
          a sparse-tracking nav environment motivate a larger standoff than the
          inner-planet 300 km.
        - Jupiter (5000 km): the intense radiation belts above the cloud tops
          dominate; a high standoff keeps periapsis well clear of the worst
          dose region (engineering convention, not a belt-model output).
        - Saturn (5000 km): keeps periapsis well above the main ring system and
          upper atmosphere (a flyby must not graze the rings); convention.
        - Uranus, Neptune (1000 km): conservative atmosphere-plus-margin
          defaults pending mission-specific analysis; convention.

Out of scope here (deliberately, per ``docs/phases/m0-scaffold/plan.md``
§4.4): Sun radius, moons, asteroids. Planet eccentricity (``ecc``) and
inclination (``inc_deg``/``lan_deg``) are now carried as sourced J2000
elements for the 3-D Tisserand predicate; the circular-coplanar ephemeris
backend still treats both as zero (it ignores these fields).
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

STANDARD_GRAVITY_M_S2: Final[float] = 9.80665
"""Standard gravitational acceleration ``g0``, m/s^2 (CGPM 1901, exact).

The conventional reference acceleration in the rocket equation
``m_{i+1} = m_i * exp(-dv / (g0 * Isp))`` (Yam, Di Lorenzo & Izzo 2010, Eq. 5;
see ``docs/v2-future-references.md`` §1). Used by the Sims-Flanagan low-thrust
leg model. Note: this is the propulsion ``g0`` convention, not a local
gravitational field strength. Value fixed by the 3rd CGPM (1901).
"""

STANDARD_GRAVITY_KM_S2: Final[float] = STANDARD_GRAVITY_M_S2 / 1000.0
"""Standard gravitational acceleration ``g0`` in km/s^2 (derived).

Convenience form for the rocket equation when speeds are in km/s and ``Isp`` is
in seconds, so ``g0 * Isp`` carries units of km/s like the per-segment ``dv``.
"""


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
    inc_deg:
        Orbital inclination wrt the J2000 ecliptic, deg. Defaults to ``0.0``
        (the coplanar idealisation). Sourced J2000 values from Standish &
        Williams, "Approximate Positions of the Planets", JPL Solar System
        Dynamics, Table 1 (valid 1800-2050 AD).
    lan_deg:
        Longitude of the ascending node wrt the J2000 ecliptic, deg. Defaults
        to ``0.0``. Same Standish & Williams Table 1 source as ``inc_deg``.
    ecc:
        Mean orbital eccentricity at J2000. Defaults to ``0.0`` (the coplanar/
        circular idealisation). Sourced J2000 values from Standish & Williams,
        "Approximate Positions of the Planets", JPL Solar System Dynamics,
        Table 1, ``e_0`` column (same source as ``sma_au``): Venus 0.00677727,
        Earth 0.01671123, Mars 0.09340065, Mercury 0.20563593,
        Jupiter 0.04838624, Saturn 0.05386179, Uranus 0.04725744,
        Neptune 0.00859048. Consumed by the 3-D Tisserand predicate; the
        circular ephemeris backend ignores it.
    """

    name: str
    code: str
    mu_km3_s2: float
    radius_eq_km: float
    sma_au: float
    mean_motion_deg_day: float
    safe_alt_km: float
    # Optional 3D orbital elements; zero-default keeps the coplanar model
    # byte-identical for callers (and bodies) that don't set them.
    inc_deg: float = 0.0
    lan_deg: float = 0.0
    ecc: float = 0.0


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
# Planet table — heliocentric bodies (all eight major planets).
# sma_au values: Standish & Williams, "Approximate Positions of the Planets",
# JPL Solar System Dynamics, Table 1, a_0 column (1800-2050 AD).
# ---------------------------------------------------------------------------

_VENUS_SMA_AU: Final[float] = 0.72333566
_EARTH_SMA_AU: Final[float] = 1.00000261
_MARS_SMA_AU: Final[float] = 1.52371034
_MERCURY_SMA_AU: Final[float] = 0.38709927
_JUPITER_SMA_AU: Final[float] = 5.20288700
_SATURN_SMA_AU: Final[float] = 9.53667594
_URANUS_SMA_AU: Final[float] = 19.18916464
_NEPTUNE_SMA_AU: Final[float] = 30.06992276

PLANETS: Final[dict[str, PlanetData]] = {
    "V": PlanetData(
        name="Venus",
        code="V",
        mu_km3_s2=3.24858592e5,
        radius_eq_km=6051.8,
        sma_au=_VENUS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_VENUS_SMA_AU),
        safe_alt_km=300.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.00677727,
        # inc_deg/lan_deg deliberately left at the coplanar default 0.0 here so
        # the live ``circular`` backend stays byte-identical for every existing
        # caller and golden. The sourced J2000 3D elements (Standish & Williams
        # Table 1: Venus inc=3.39467605 deg, lan=76.67984255 deg) are exercised
        # via the inclined backend with injected PlanetData; STAGE 4 / 3D work
        # opts in explicitly rather than mutating the shared coplanar model.
    ),
    "E": PlanetData(
        name="Earth",
        code="E",
        mu_km3_s2=3.98600435507e5,
        radius_eq_km=6378.137,
        sma_au=_EARTH_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_EARTH_SMA_AU),
        safe_alt_km=300.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.01671123,
    ),
    "M": PlanetData(
        name="Mars",
        code="M",
        mu_km3_s2=4.282837521e4,
        radius_eq_km=3396.19,
        sma_au=_MARS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_MARS_SMA_AU),
        safe_alt_km=300.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.09340065,
        # Left at coplanar default 0.0 (see the Venus note above). Sourced
        # J2000 3D elements (Standish & Williams Table 1: Mars inc=1.84969142
        # deg, lan=49.55953891 deg) are exercised via the inclined backend.
    ),
    "Me": PlanetData(
        name="Mercury",
        code="Me",
        # JPL DE440 (Park et al., AJ 2021); planet-only (Mercury has no moons).
        mu_km3_s2=2.2031868551e4,
        # IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22), mean eq. radius.
        radius_eq_km=2440.53,
        sma_au=_MERCURY_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_MERCURY_SMA_AU),
        # Engineering default (convention, not sourced physics): no atmosphere,
        # but solar-thermal + sparse-tracking nav margins → larger standoff.
        safe_alt_km=1000.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.20563593,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Mercury
        # inc=7.00497902 deg, lan=48.33076593 deg) live in the inclined backend.
    ),
    "J": PlanetData(
        name="Jupiter",
        code="J",
        # JPL DE440 system GM (Park et al., AJ 2021); satellites contribute
        # <2e-4 and orbit far above any flyby periapsis, so the system value is
        # the central-body mu sensed at closest approach (see module docstring).
        mu_km3_s2=1.267127641e8,
        # IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22), 1-bar eq. radius.
        radius_eq_km=71492.0,
        sma_au=_JUPITER_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_JUPITER_SMA_AU),
        # Engineering default (convention, not sourced physics): radiation belts
        # above the cloud tops dominate → high standoff keeps periapsis clear.
        safe_alt_km=5000.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.04838624,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Jupiter
        # inc=1.30439695 deg, lan=100.47390909 deg) live in the inclined backend.
    ),
    "S": PlanetData(
        name="Saturn",
        code="S",
        # JPL DE440 system GM (Park et al., AJ 2021); see Jupiter note.
        mu_km3_s2=3.79405848418e7,
        # IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22), 1-bar eq. radius.
        radius_eq_km=60268.0,
        sma_au=_SATURN_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_SATURN_SMA_AU),
        # Engineering default (convention, not sourced physics): keeps periapsis
        # well above the main ring system + upper atmosphere (must not graze).
        safe_alt_km=5000.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.05386179,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Saturn
        # inc=2.48599187 deg, lan=113.66242448 deg) live in the inclined backend.
    ),
    "U": PlanetData(
        name="Uranus",
        code="U",
        # JPL DE440 system GM (Park et al., AJ 2021); see Jupiter note.
        mu_km3_s2=5.7945564e6,
        # IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22), 1-bar eq. radius.
        radius_eq_km=25559.0,
        sma_au=_URANUS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_URANUS_SMA_AU),
        # Engineering default (convention, not sourced physics): conservative
        # atmosphere-plus-margin standoff pending mission-specific analysis.
        safe_alt_km=1000.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.04725744,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Uranus
        # inc=0.77263783 deg, lan=74.01692503 deg) live in the inclined backend.
    ),
    "N": PlanetData(
        name="Neptune",
        code="N",
        # JPL DE440 system GM (Park et al., AJ 2021); see Jupiter note.
        mu_km3_s2=6.83652710058e6,
        # IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22), 1-bar eq. radius.
        radius_eq_km=24764.0,
        sma_au=_NEPTUNE_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_NEPTUNE_SMA_AU),
        # Engineering default (convention, not sourced physics): conservative
        # atmosphere-plus-margin standoff pending mission-specific analysis.
        safe_alt_km=1000.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.00859048,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Neptune
        # inc=1.77004347 deg, lan=131.78422574 deg) live in the inclined backend.
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

SUPPORTED_BODIES: Final[tuple[str, ...]] = tuple(PLANETS.keys())
"""Body codes the compute machinery supports, derived from :data:`PLANETS`.

The single source of truth for "which bodies exist": adding a (sourced)
:class:`PlanetData` entry to ``PLANETS`` is all that is needed to make a new
body resolvable everywhere — ephemeris, frames, flyby, and maintenance all key
off this table rather than hardcoding their own V/E/M sets.
"""
