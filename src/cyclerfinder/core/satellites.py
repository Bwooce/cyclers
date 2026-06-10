"""Planet-centric satellite (moon) registry — Tier-1 patched-conic moon tours.

A moon cannot live in :data:`cyclerfinder.core.constants.PLANETS`:
``PlanetData.sma_au`` is intrinsically heliocentric (constants.py:139) and its
mean motion derives from ``MU_SUN``. This is the about-the-primary sibling.

Body-code scheme: FULL MOON NAMES (data/README.md:65-79) —
``"Io"``/``"Europa"``/``"Ganymede"``/``"Callisto"``/``"Titan"``/… — to avoid
collision with the heliocentric V/E/M planet codes. ``mean_motion_deg_day`` is
DERIVED at import from ``sma_km`` + the primary's mu (Kepler III), never
hand-copied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SatelliteData:
    name: str
    code: str  # full moon name (scheme: data/README.md:65-79)
    primary: str  # "Jupiter" | "Saturn" | "Earth" | "Mars"
    mu_km3_s2: float  # moon GM
    radius_eq_km: float
    sma_km: float  # SMA ABOUT THE PRIMARY (km, not AU, not Sun-relative)
    mean_motion_deg_day: float  # about the primary (use mean_motion_deg_day_about)
    safe_alt_km: float


def mean_motion_deg_day_about(sma_km: float, *, mu_primary: float) -> float:
    """Mean motion (deg/day) about a primary, Kepler III (cf. constants.py:149-159)."""
    period_s = 2.0 * math.pi * math.sqrt(sma_km**3 / mu_primary)
    return 360.0 / (period_s / 86400.0)


# Primary (central-body) gravitational parameters GM [km^3/s^2].
#
# SOURCING (golden discipline): these are the JPL Solar System Dynamics planetary
# *system* GM values (http://ssd.jpl.nasa.gov/, gm_de440 planetary constants),
# the SAME upstream the Endgame Part-1 Table 3 footnote cites. They are sourced
# INDEPENDENTLY of the paper's transcribed V_M anchors so the Task-1.3 golden
# (registry vs published Table 3) stays non-circular. Accessed 2026-06-07.
PRIMARIES: dict[str, float] = {
    # Earth system GM (JPL SSD gm_de440 planetary constants): 4.0350323562548019e5 km^3/s^2.
    # Used for the Earth-Moon CR3BP (Task 4, 2026-06-10).
    "Earth": 4.0350323562548019e5,
    # Jupiter system GM (JPL SSD gm_de440 planetary constants): 1.26686534e8.
    "Jupiter": 1.26686534e8,
    # Saturn system GM (JPL SSD gm_de440 planetary constants): 3.7931207e7.
    "Saturn": 3.7931207e7,
}


def _sat(
    name: str,
    primary: str,
    mu_km3_s2: float,
    radius_eq_km: float,
    sma_km: float,
    safe_alt_km: float,
) -> SatelliteData:
    """Build a SatelliteData with mean motion DERIVED from sma + primary mu."""
    return SatelliteData(
        name=name,
        code=name,  # full-name scheme
        primary=primary,
        mu_km3_s2=mu_km3_s2,
        radius_eq_km=radius_eq_km,
        sma_km=sma_km,
        mean_motion_deg_day=mean_motion_deg_day_about(sma_km, mu_primary=PRIMARIES[primary]),
        safe_alt_km=safe_alt_km,
    )


# Satellite registry keyed by full moon name.
#
# SOURCING (golden, non-negotiable): mu_km3_s2 / radius_eq_km / sma_km come from
# the JPL Solar System Dynamics satellite physical-parameter tables
# (http://ssd.jpl.nasa.gov/, satellite GMs from the jup365/sat441 ephemeris fits;
# mean radii + semi-major axes from the SSD satellite physical/orbital tables),
# accessed 2026-06-07. They are sourced INDEPENDENTLY of the Endgame Part-1
# Table 3 transcription — the paper's a_M / V_M are the cross-check anchor
# (Task 1.3), not the registry source.
#
# safe_alt_km = the altitude the paper used for its Delta-V tables so the Phase-5
# anchors are comparable: 100 km for all moons except Titan at 1500 km (mining
# note line 352).
SATELLITES: dict[str, SatelliteData] = {
    # --- Earth-Moon system ---
    # Moon: GM 4902.800118 km^3/s^2 (JPL SSD satellite GM table, gm_de440);
    #       mean radius 1737.4 km; SMA 384400 km (JPL SSD satellite physical/orbital
    #       tables); safe_alt 100 km (standard low lunar orbit altitude).
    #       Sources accessed 2026-06-10. Added for Earth-Moon CR3BP mu (Task 4).
    "Moon": _sat("Moon", "Earth", 4902.800118, 1737.4, 384400.0, 100.0),
    # --- Galilean (Jupiter) ---
    # Io: GM 5959.916, mean R 1821.49 km, a 421800 km (JPL SSD).
    "Io": _sat("Io", "Jupiter", 5959.916, 1821.49, 421800.0, 100.0),
    # Europa: GM 3202.739, mean R 1560.8 km, a 671100 km (JPL SSD).
    "Europa": _sat("Europa", "Jupiter", 3202.739, 1560.8, 671100.0, 100.0),
    # Ganymede: GM 9887.834, mean R 2631.2 km, a 1070400 km (JPL SSD).
    "Ganymede": _sat("Ganymede", "Jupiter", 9887.834, 2631.2, 1070400.0, 100.0),
    # Callisto: GM 7179.289, mean R 2410.3 km, a 1882700 km (JPL SSD).
    "Callisto": _sat("Callisto", "Jupiter", 7179.289, 2410.3, 1882700.0, 100.0),
    # --- Saturnian midsize + Titan (Saturn) ---
    # Mimas: GM 2.503, mean R 198.2 km, a 185540 km (JPL SSD).
    "Mimas": _sat("Mimas", "Saturn", 2.503, 198.2, 185540.0, 100.0),
    # Enceladus: GM 7.211, mean R 252.1 km, a 238040 km (JPL SSD).
    "Enceladus": _sat("Enceladus", "Saturn", 7.211, 252.1, 238040.0, 100.0),
    # Tethys: GM 41.21, mean R 533.0 km, a 294670 km (JPL SSD).
    "Tethys": _sat("Tethys", "Saturn", 41.21, 533.0, 294670.0, 100.0),
    # Dione: GM 73.116, mean R 561.4 km, a 377420 km (JPL SSD).
    "Dione": _sat("Dione", "Saturn", 73.116, 561.4, 377420.0, 100.0),
    # Rhea: GM 153.94, mean R 763.8 km, a 527070 km (JPL SSD).
    "Rhea": _sat("Rhea", "Saturn", 153.94, 763.8, 527070.0, 100.0),
    # Titan: GM 8978.14, mean R 2574.7 km, a 1221870 km (JPL SSD); 1500 km alt.
    "Titan": _sat("Titan", "Saturn", 8978.14, 2574.7, 1221870.0, 1500.0),
}
