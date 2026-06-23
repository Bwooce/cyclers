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

Dwarf planets / planetoids (added 2026-06-14, #260):
    The eight major planets are joined by the IAU dwarf planets Pluto, Ceres,
    Eris, Makemake, Haumea and the two largest main-belt asteroids Vesta and
    Pallas, so every meaningful heliocentric body is available for the
    Tisserand/V_inf screen to prune (most are LOW-MASS -> negligible gravity
    assist; they are added for COMPLETENESS, not because each is a useful flyby
    body). Their sourcing:
      * ``mu_km3_s2`` (GM) —
          - Ceres, Pallas, Vesta, Pluto(system): JPL DE440 small-body /
            planetary constants kernel ``gm_de440.tpc`` (Park et al., AJ 2021).
            Ceres 62.6289, Pallas 13.6659, Vesta 17.2882, Pluto system 975.5.
          - Eris, Makemake, Haumea: JPL/SBDB publishes no GM, only a mass from
            satellite dynamics, so GM = G*M with G = 6.67430e-20 km^3 kg^-1 s^-2
            (CODATA 2018). Eris M=(1.638±0.014)e22 kg (Brown & Schaller 2007,
            Science 316:1585); Makemake M=(2.69±0.20)e21 kg (Bamberger et al.
            2025, from moon S/2015 (136472) 1); Haumea M=(3.952±0.011)e21 kg
            (Proudfoot et al. 2024). These are MASS-derived, hence less precise
            than the DE440 GMs — flagged inline. Pluto here is the planet-only
            heliocentric body (its system GM lives in satellites.PRIMARIES).
      * ``radius_eq_km`` — mean radii: Pluto 1188.3, Ceres 469.7 (IAU 2015
        WGCCRE / Dawn, Archinal et al. 2018), Eris 1163, Makemake 715,
        Haumea ~780 (volume-equivalent), Vesta 262.7, Pallas 256.5 (mean radii,
        JPL SSD / Russell et al. 2012). Used only for the flyby periapsis floor.
      * ``sma_au`` / ``ecc`` / time-phasing angles —
          - Pluto: Standish & Williams "Keplerian Elements for Approximate
            Positions of the Major Planets" (JPL SSD), Table 2a (3000 BC-3000
            AD) — the only table that retains Pluto (Table 1 dropped it post-IAU
            2006). a0=39.48686035, e0=0.24885238, I0=17.14104260,
            L0=238.96535011, varpi0=224.09702598, Omega0=110.30167986.
          - Ceres/Eris/Makemake/Haumea + Vesta/Pallas: JPL Small-Body Database
            osculating heliocentric elements (NASA/JPL SBDB ``full-prec=1``), all
            at the common epoch JD 2461200.5 (TDB), accessed 2026-06-14. These
            are osculating (not a mean-element table), so the time-phasing
            varpi/L0 are left at 0.0 — only sma/ecc/inc feed the Tisserand
            screen. Ceres a=2.76555 e=0.07969 i=10.588; Eris a=67.9339 e=0.43824
            i=43.926; Makemake a=45.5709 e=0.15889 i=29.028; Haumea a=43.0603
            e=0.19444 i=28.208; Vesta a=2.36137 e=0.09020 i=7.144; Pallas
            a=2.76956 e=0.23070 i=34.933.
      * ``safe_alt_km`` — engineering default 100 km for all; convention, not
        sourced physics — EXCEPT Pluto, whose 100 km is now SOURCED as a design
        floor (Stern-Tapley-Finley-Scherrer 2020 JSR A34658 tour spec 3: Pluto
        periapse 100-500 km for in-situ atmospheric measurement; #429). Same
        value -> no change, provenance upgrade only. See
        data/flyby_altitude_references.yaml + 2026-06-23-digest-stern-2020-pluto-orbiter.md.
    The non-coplanar inclinations (Pluto 17°, Eris 44°, Makemake 29°,
    Haumea 28°) matter for the 3-D Tisserand screen but, per the existing
    convention, ``inc_deg``/``lan_deg`` stay at the coplanar 0.0 default in
    PLANETS and are exposed via the inclined backend's element map instead.

Out of scope here (deliberately, per ``docs/phases/m0-scaffold/plan.md``
§4.4): Sun radius, moons, smaller asteroids. Planet eccentricity (``ecc``) and
inclination (``inc_deg``/``lan_deg``) are now carried as sourced J2000
elements for the 3-D Tisserand predicate; the circular-coplanar ephemeris
backend still treats both as zero (it ignores these fields). The
time-phasing angles ``varpi_deg`` (longitude of perihelion ϖ) and ``L0_deg``
(mean longitude at J2000) are likewise sourced from Standish & Williams
Table 1 for the time-true visualization (planet-elements.json emitter); no
ephemeris backend reads them either.
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
    varpi_deg:
        Longitude of perihelion ``ϖ = Ω + ω`` wrt the J2000 ecliptic, deg.
        Defaults to ``0.0``. Sourced J2000 values from Standish & Williams,
        "Approximate Positions of the Planets", JPL Solar System Dynamics,
        Table 1, ``varpi_0`` (longitude of perihelion) column — the SAME table
        as ``sma_au``/``ecc``/``inc_deg``. Needed to place a planet at a real
        date (``ω = ϖ - Ω``). The circular/inclined ephemeris backends ignore
        it; it is consumed by the time-true viz (planet-elements.json emitter).
    L0_deg:
        Mean longitude at J2000 ``L0 = ϖ + M0`` wrt the J2000 ecliptic, deg.
        Defaults to ``0.0``. Sourced J2000 values from Standish & Williams,
        Table 1, ``L_0`` (mean longitude) column. Needed to place a planet on
        its ellipse at a given date (``M0 = L0 - ϖ``). Same provenance and same
        consumer (time-true viz) as ``varpi_deg``; the ephemeris backends
        ignore it.
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
    # Time-phasing angles (Standish & Williams Table 1 varpi_0 / L_0 columns).
    # No ephemeris backend reads these; they feed the time-true viz emitter.
    varpi_deg: float = 0.0
    L0_deg: float = 0.0


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

# Dwarf planets / planetoids (added 2026-06-14, #260).
# Pluto: Standish & Williams Table 2a a0 (the only Standish table with Pluto).
# The rest: JPL SBDB osculating a at epoch JD 2461200.5 (see module docstring).
_PLUTO_SMA_AU: Final[float] = 39.48686035
_CERES_SMA_AU: Final[float] = 2.765552595034094
_ERIS_SMA_AU: Final[float] = 67.93394687853566
_MAKEMAKE_SMA_AU: Final[float] = 45.57093317300052
_HAUMEA_SMA_AU: Final[float] = 43.06029023650952
_VESTA_SMA_AU: Final[float] = 2.361365965127599
_PALLAS_SMA_AU: Final[float] = 2.769559010737709

# G (CODATA 2018), km^3 kg^-1 s^-2 — for the mass-derived GMs (Eris/Makemake/
# Haumea), which JPL/SBDB publish only as a mass, not a GM.
_G_KM3_KG_S2: Final[float] = 6.67430e-20

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
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=131.60246718,
        L0_deg=181.97909950,
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
        # SOURCED design floor: Russell 2004 p.165 r_p,min,Earth = 6578.0 km (= 200 km),
        # Table 3.4 footnote b ("all flybys min altitude > 200 km" = strictly ballistic).
        # The model the Earth-Mars cyclers were designed under; the prior unsourced 300 km
        # spuriously charged ~40 m/s to S1L1's one marginal flyby. Mission cross-check:
        # Galileo EGA2 flew 303 km (D'Amario-Bright-Wolf 1992) — above this physical floor.
        safe_alt_km=200.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.01671123,
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (Earth-Moon barycentre row).
        varpi_deg=102.93768193,
        L0_deg=100.46457166,
    ),
    "M": PlanetData(
        name="Mars",
        code="M",
        mu_km3_s2=4.282837521e4,
        radius_eq_km=3396.19,
        sma_au=_MARS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_MARS_SMA_AU),
        # SOURCED design floor: Russell 2004 p.165 r_p,min,Mars = 3598.5 km (= 200 km),
        # Table 3.4 footnote b. Parallels Earth; the model the cyclers were designed under.
        safe_alt_km=200.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.09340065,
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=-23.94362959,
        L0_deg=-4.55343205,
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
        # 200 km: DEMONSTRATED practical floor (#428) — BepiColombo flew Mercury at
        # 199/200/236 km (Mercury-1/2/3); Mariner-10 ~703 km. Mercury is airless so the
        # floor is surface + nav margin, not atmosphere. Not a published *design minimum*
        # (none in corpus); the prior 1000 km was an unsourced conservative convention.
        # docs/notes/2026-06-23-flyby-altitude-corpus-mining.md (Mercury section).
        safe_alt_km=200.0,
        # J2000 mean eccentricity, Standish & Williams Table 1, e_0 column.
        ecc=0.20563593,
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=77.45779628,
        L0_deg=252.25032350,
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
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=14.72847983,
        L0_deg=34.39644051,
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
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=92.59887831,
        L0_deg=49.95424423,
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
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=170.95427630,
        L0_deg=313.23810451,
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
        # J2000 longitude of perihelion / mean longitude, Standish & Williams
        # Table 1, varpi_0 / L_0 columns (time-true viz phasing angles).
        varpi_deg=44.96476227,
        L0_deg=-55.12002969,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced J2000 3D elements (Standish & Williams Table 1: Neptune
        # inc=1.77004347 deg, lan=131.78422574 deg) live in the inclined backend.
    ),
    # -----------------------------------------------------------------------
    # Dwarf planets / planetoids (added 2026-06-14, #260). LOW-MASS bodies
    # added for COMPLETENESS so the Tisserand/V_inf screen can self-prune;
    # most are poor flyby/cycler nodes. Two-letter codes (single letters are
    # taken). Full sourcing in the module docstring.
    # -----------------------------------------------------------------------
    "Pl": PlanetData(
        name="Pluto",
        code="Pl",
        # Planet-only heliocentric body. GM = Pluto+Charon system GM minus
        # Charon would be ideal, but DE440 publishes only the system value and
        # Charon is ~11% of it; for the heliocentric screen the system GM is the
        # mass that perturbs a flyby, so we use the DE440 system GM (975.5).
        mu_km3_s2=9.755e2,
        # Mean radius 1188.3 km (New Horizons / IAU 2015 WGCCRE, Archinal 2018).
        radius_eq_km=1188.3,
        sma_au=_PLUTO_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_PLUTO_SMA_AU),
        # 100 km: SOURCED design floor (Stern-Tapley-Finley-Scherrer 2020 JSR
        # A34658 tour spec 3 — Pluto periapse 100-500 km for in-situ ATMOSPHERIC
        # measurement; #429). Was an unsourced "no atmosphere" convention at the
        # same value; Pluto does have a thin atmosphere, so the 100 km is a real
        # atmospheric-pass design minimum. Provenance upgrade, no value change.
        safe_alt_km=100.0,
        # Standish & Williams Table 2a (3000 BC-3000 AD) — the only Standish
        # table retaining Pluto. e0/varpi0/L0 from that same row.
        ecc=0.24885238,
        varpi_deg=224.09702598,
        L0_deg=238.96535011,
        # inc_deg/lan_deg left at coplanar default 0.0 (see the Venus note).
        # Sourced Table 2a 3D elements (inc=17.14104260 deg, lan=110.30167986
        # deg) live in the inclined backend element map (ephemeris.py).
    ),
    "Ce": PlanetData(
        name="Ceres",
        code="Ce",
        # JPL DE440 small-body constants (gm_de440.tpc, body 2000001): 62.6289.
        mu_km3_s2=6.26288886444e1,
        # Mean radius 469.7 km (Dawn / IAU 2015 WGCCRE, Archinal et al. 2018).
        radius_eq_km=469.7,
        sma_au=_CERES_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_CERES_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5 (see module docstring).
        ecc=0.07969229514816586,
        # inc=10.58802780183462 deg lives in the inclined backend element map.
    ),
    "Er": PlanetData(
        name="Eris",
        code="Er",
        # No published GM; GM = G*M, M=(1.638e22) kg (Brown & Schaller 2007,
        # Science 316:1585, from Dysnomia's orbit). MASS-derived (less precise).
        mu_km3_s2=_G_KM3_KG_S2 * 1.638e22,
        # Mean radius 1163 km (Sicardy et al. 2011 stellar occultation).
        radius_eq_km=1163.0,
        sma_au=_ERIS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_ERIS_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5.
        ecc=0.4382385347971672,
        # inc=43.9258279471791 deg lives in the inclined backend element map.
    ),
    "Mk": PlanetData(
        name="Makemake",
        code="Mk",
        # No published GM; GM = G*M, M=(2.69e21) kg (Bamberger et al. 2025, from
        # the moon S/2015 (136472) 1). MASS-derived (less precise).
        mu_km3_s2=_G_KM3_KG_S2 * 2.69e21,
        # Mean radius 715 km (Ortiz et al. 2012 stellar occultation).
        radius_eq_km=715.0,
        sma_au=_MAKEMAKE_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_MAKEMAKE_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5.
        ecc=0.1588889953992523,
        # inc=29.02785603743067 deg lives in the inclined backend element map.
    ),
    "Ha": PlanetData(
        name="Haumea",
        code="Ha",
        # No published GM; GM = G*M, M=(3.952e21) kg (Proudfoot et al. 2024,
        # from the satellites Hi'iaka/Namaka). MASS-derived (less precise).
        mu_km3_s2=_G_KM3_KG_S2 * 3.952e21,
        # Volume-equivalent mean radius ~780 km (Ortiz et al. 2017 occultation;
        # Haumea is markedly triaxial, so this is the equivalent-sphere radius).
        radius_eq_km=780.0,
        sma_au=_HAUMEA_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_HAUMEA_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5.
        ecc=0.1944430148898797,
        # inc=28.20847393040364 deg lives in the inclined backend element map.
    ),
    "Ve": PlanetData(
        name="Vesta",
        code="Ve",
        # JPL DE440 small-body constants (gm_de440.tpc, body 2000004): 17.2882.
        mu_km3_s2=1.72882328792e1,
        # Mean radius 262.7 km (Dawn, Russell et al. 2012).
        radius_eq_km=262.7,
        sma_au=_VESTA_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_VESTA_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5.
        ecc=0.09020374382834395,
        # inc=7.143925545058711 deg lives in the inclined backend element map.
    ),
    "Pa": PlanetData(
        name="Pallas",
        code="Pa",
        # JPL DE440 small-body constants (gm_de440.tpc, body 2000002): 13.6659.
        mu_km3_s2=1.36658781460e1,
        # Mean radius 256.5 km (Marsset et al. 2020 / SPHERE imaging).
        radius_eq_km=256.5,
        sma_au=_PALLAS_SMA_AU,
        mean_motion_deg_day=_mean_motion_deg_day(_PALLAS_SMA_AU),
        safe_alt_km=100.0,
        # JPL SBDB osculating e at epoch JD 2461200.5.
        ecc=0.2307000995648547,
        # inc=34.93279321851542 deg lives in the inclined backend element map.
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


def vinf_ceiling_kms(body: str) -> float:
    """Elliptic-periodicity hyperbolic-excess ceiling at ``body`` (km/s).

    PHYSICS, not a convention. A heliocentric cycler orbit must be *bound*
    (closed/periodic), i.e. the spacecraft's heliocentric speed at the
    encounter radius ``r_B`` cannot exceed the local solar-escape speed
    ``v_esc_sun(r_B) = sqrt(2 * mu_Sun / r_B)`` — a faster vehicle is on a
    hyperbolic heliocentric trajectory and never returns, so it cannot be a
    periodic cycler. The hyperbolic excess at the body is the vector difference
    ``V_inf = v_sc - v_body``; its magnitude is maximised (worst case) when the
    spacecraft moves *anti-parallel* to the body, giving the bound

        |V_inf| <= v_esc_sun(r_B) + |v_body|

    where ``|v_body| = sqrt(mu_Sun / r_B)`` is the body's circular orbital
    speed (the planets are very nearly circular). This is a hard upper limit
    for ANY periodic heliocentric orbit at that body — exceeding it means the
    encounter cannot belong to a closed cycler. At Earth it evaluates to
    ``42.122 + 29.785 = 71.91 km/s``.

    Use as a BUG ceiling, never a family filter: legitimate high-energy
    catalogue rows (Russell-Ocampo reaches 20.3 km/s at Earth) sit far below
    it. A value above the ceiling signals a degenerate / unit-error / off-family
    solve, not a high-energy cycler.
    """
    p = PLANETS[body]
    r_b_km = p.sma_au * AU_KM
    v_body = sqrt(MU_SUN_KM3_S2 / r_b_km)
    v_esc_sun = sqrt(2.0 * MU_SUN_KM3_S2 / r_b_km)
    return v_esc_sun + v_body


VINF_CEILING_KMS: Final[dict[str, float]] = {code: vinf_ceiling_kms(code) for code in PLANETS}
"""Per-body elliptic-periodicity V_inf ceiling (km/s); see :func:`vinf_ceiling_kms`.

Precomputed for every supported body. Values (km/s): Venus 84.55, Earth 71.91,
Mars 58.25, Mercury 115.57, Jupiter 31.52, Saturn 23.29, Uranus 16.42,
Neptune 13.11. A periodic heliocentric encounter V_inf above the body's entry
here is physically impossible and marks a BUG, not a high-energy solution.
"""
