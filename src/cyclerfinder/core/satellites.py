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
    primary: str  # "Earth" | "Mars" | "Jupiter" | "Saturn" | "Uranus" | "Neptune" | "Pluto"
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
    # Mars system GM (JPL DE440 planetary constants, ssd.jpl.nasa.gov astro_par /
    # gm_de440.tpc): 4.282837521e4 km^3/s^2. Phobos/Deimos add <2e-8 of this, so
    # the system value is the body GM the moon-frame Kepler sees. Accessed 2026-06-14.
    "Mars": 4.282837521e4,
    # Jupiter system GM (JPL SSD gm_de440 planetary constants): 1.26686534e8.
    "Jupiter": 1.26686534e8,
    # Saturn system GM (JPL SSD gm_de440 planetary constants): 3.7931207e7.
    "Saturn": 3.7931207e7,
    # Uranus system GM (JPL DE440 planetary constants, ssd.jpl.nasa.gov astro_par):
    # 5.794556400e6 km^3/s^2. Accessed 2026-06-14.
    "Uranus": 5.7945564e6,
    # Neptune system GM (JPL DE440 planetary constants, ssd.jpl.nasa.gov astro_par):
    # 6.836527100580e6 km^3/s^2. Accessed 2026-06-14.
    "Neptune": 6.836527100580e6,
    # Pluto system GM (JPL DE440 / gm_de440.tpc BODY9_GM, Brozovic et al. 2015):
    # 9.755e2 km^3/s^2. The Pluto-Charon pair is a near-binary; this is the
    # Pluto+Charon SYSTEM GM. Accessed 2026-06-14.
    "Pluto": 9.755e2,
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
# note line 352). For the bodies added 2026-06-14 (Mars/Uranus/Neptune/Pluto
# moons + Jupiter Amalthea + Saturn Iapetus/Hyperion) there is no paper-anchored
# altitude, so safe_alt_km is an ENGINEERING DEFAULT scaled to the moon's size:
# 10 km for the tiny irregulars (Phobos/Deimos/Amalthea/Hyperion/Nix/Hydra,
# mean radius < ~140 km), 100 km otherwise. This is a convention, not a sourced
# constant; the Tisserand/V_inf screen self-prunes the low-mass bodies regardless.
#
# CAVEAT (completeness, not advocacy): most of the 2026-06-14 additions are
# LOW-MASS and make poor gravity-assist / cycler bodies. Phobos/Deimos and the
# small irregulars (Amalthea, Hyperion, Nix, Hydra) have GM < 1 km^3/s^2 ->
# negligible bending. Triton is large but RETROGRADE and inclined (a hostile
# capture/flyby geometry). The genuinely interesting additions are the five
# regular Uranian moons (Miranda/Ariel/Umbriel/Titania/Oberon) and the
# Pluto-Charon pair (Charon's GM is ~11% of the system -> a real binary). They
# are added so the screen has the full body set to prune from, not because every
# one is a useful cycler node.
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
    # --- Mars system (added 2026-06-14) ---
    # All values JPL SSD satellite physical (phys_par, ref MAR097) + mean elements
    # (sats/elem), accessed 2026-06-14. Tiny bodies: GM ~ 7e-4 / 1e-4 -> screen
    # self-prunes. safe_alt 10 km (irregular, mean R ~ 6-11 km).
    # Phobos: GM 0.0007087 km^3/s^2, mean R 11.08 km, a 9375 km.
    "Phobos": _sat("Phobos", "Mars", 0.0007087, 11.08, 9375.0, 10.0),
    # Deimos: GM 0.0000962 km^3/s^2, mean R 6.2 km, a 23457 km.
    "Deimos": _sat("Deimos", "Mars", 0.0000962, 6.2, 23457.0, 10.0),
    # --- Uranus system (added 2026-06-14) ---
    # All values JPL SSD satellite physical (phys_par, ref URA111) + mean elements
    # (sats/elem), accessed 2026-06-14. The five classical regular moons; the most
    # interesting outer-planet additions (moderate GM, prograde, low inclination).
    # Miranda: GM 4.3 km^3/s^2, mean R 235.8 km, a 129846 km.
    "Miranda": _sat("Miranda", "Uranus", 4.3, 235.8, 129846.0, 100.0),
    # Ariel: GM 83.5 km^3/s^2, mean R 578.9 km, a 190929 km.
    "Ariel": _sat("Ariel", "Uranus", 83.5, 578.9, 190929.0, 100.0),
    # Umbriel: GM 85.1 km^3/s^2, mean R 584.7 km, a 265986 km.
    "Umbriel": _sat("Umbriel", "Uranus", 85.1, 584.7, 265986.0, 100.0),
    # Titania: GM 226.9 km^3/s^2, mean R 788.9 km, a 436298 km.
    "Titania": _sat("Titania", "Uranus", 226.9, 788.9, 436298.0, 100.0),
    # Oberon: GM 205.3 km^3/s^2, mean R 761.4 km, a 583511 km.
    "Oberon": _sat("Oberon", "Uranus", 205.3, 761.4, 583511.0, 100.0),
    # --- Neptune system (added 2026-06-14) ---
    # JPL SSD satellite physical (phys_par, ref NEP097/NEP101) + mean elements,
    # accessed 2026-06-14. Triton is large (GM ~ 1428) but RETROGRADE + inclined
    # (a ~ 354800 km on a near-circular retrograde orbit) -> hostile flyby
    # geometry; carried for completeness. Proteus added (regular, prograde).
    # Nereid OMITTED: JPL SSD lists its GM as 0.0 (mass not determined), so per
    # the sourcing discipline we omit it rather than guess from a size estimate.
    # Triton: GM 1428.495 km^3/s^2, mean R 1352.6 km, a 354800 km.
    "Triton": _sat("Triton", "Neptune", 1428.49546, 1352.6, 354800.0, 100.0),
    # Proteus: GM 2.58342 km^3/s^2, mean R 208.0 km, a 117600 km.
    "Proteus": _sat("Proteus", "Neptune", 2.58342, 208.0, 117600.0, 100.0),
    # --- Pluto system (added 2026-06-14) ---
    # JPL SSD satellite physical (phys_par, ref PLU060) + mean elements, accessed
    # 2026-06-14. Charon's GM (106.1) is ~11% of the Pluto+Charon system GM
    # (975.5) -> a genuine binary, the interesting addition here. Nix/Hydra are
    # tiny (GM ~ 1.5e-3 / 2e-3); included for completeness, screen self-prunes.
    # Charon: GM 106.1 km^3/s^2, mean R 606.0 km, a 19600 km.
    "Charon": _sat("Charon", "Pluto", 106.1, 606.0, 19600.0, 100.0),
    # Nix: GM 0.0015 km^3/s^2, mean R 18.0 km, a 49300 km.
    "Nix": _sat("Nix", "Pluto", 0.0015, 18.0, 49300.0, 10.0),
    # Hydra: GM 0.0020 km^3/s^2, mean R 18.5 km, a 65200 km.
    "Hydra": _sat("Hydra", "Pluto", 0.0020, 18.5, 65200.0, 10.0),
    # --- Completeness additions to existing systems (added 2026-06-14) ---
    # Jupiter Amalthea: JPL SSD phys_par ref JUP365 + mean elements, accessed
    # 2026-06-14. GM 0.16456 km^3/s^2, mean R 83.5 km, a 181400 km (tiny inner
    # irregular; screen self-prunes).
    "Amalthea": _sat("Amalthea", "Jupiter", 0.16456, 83.5, 181400.0, 10.0),
    # Saturn Iapetus: JPL SSD phys_par ref SAT441 + mean elements, accessed
    # 2026-06-14. GM 120.51511 km^3/s^2, mean R 734.3 km, a 3561700 km (regular
    # outer moon, moderate GM).
    "Iapetus": _sat("Iapetus", "Saturn", 120.51511, 734.3, 3561700.0, 100.0),
    # Saturn Hyperion: JPL SSD phys_par ref SAT441 + mean elements, accessed
    # 2026-06-14. GM 0.37049 km^3/s^2, mean R 135.0 km, a 1481500 km (small
    # chaotic-rotation moon; screen self-prunes).
    "Hyperion": _sat("Hyperion", "Saturn", 0.37049, 135.0, 1481500.0, 10.0),
}
