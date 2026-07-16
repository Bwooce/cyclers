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
    # RETROGRADE (#599): True iff this moon orbits its primary in the sense
    # OPPOSITE the primary's own rotation / the primary's other regular
    # moons -- e.g. Neptune's Triton. Default False (every other satellite
    # in this registry is prograde), so every pre-existing call site and
    # every pre-existing consumer of SatelliteData is unaffected. This flag
    # is load-bearing for two-body synodic-period + relative-phase-evolution
    # math (``synodic_period_days`` / the #558-#563 symmetric-closure
    # construction): a counter-orbiting pair conjuncts via ``1/(1/Ta+1/Tb)``,
    # not the same-sense ``1/|1/Ta-1/Tb|``, and its phase must be advanced
    # with a NEGATIVE mean motion, not a positive one.
    retrograde: bool = False


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
    # --- Small-body multi-moon systems (#607, added 2026-07-16) ---
    # GM = G*mass, G = 6.67430e-20 km^3 kg^-1 s^-2 (matches core/constants.py's
    # _G_KM3_KG_S2, so every conversion below uses the SAME constant this
    # project already uses elsewhere -- not a separately-chosen G).
    #
    # (87) Sylvia system GM: mass (1.44+/-0.01)e19 kg, DYNAMICALLY measured
    # from the Romulus/Remus mutual orbits (Kepler III), NOT a density
    # assumption -- Vernazza, Carry et al. (2021), A&A 654, A56, "VLT/SPHERE
    # imaging survey of the largest main-belt asteroids: Final results and
    # synthesis." GM = 6.67430e-20 * 1.44e19 = 0.9611 km^3/s^2.
    # CAVEAT: Sylvia's own shape is markedly non-spherical (triaxial
    # 374x248x194 km, Carry et al. 2021 ADAM model; shape-derived J2~0.024,
    # Berthier et al. 2014, Icarus 239, 118, arXiv:1407.1292) vs a near-zero
    # J2 implied by the moons' simple Keplerian orbit fits -- a documented
    # tension (non-homogeneous internal mass distribution). Both moons orbit
    # only ~5-10 Sylvia mean-radii out (mean radius ~136 km from the
    # volume-equivalent diameter ~271-274 km); point-mass is used here as a
    # first-pass approximation, NOT validated against a full multipole model.
    "Sylvia": 0.9611,
    # (130) Elektra system GM: mass (6.606 +0.007/-0.013)e18 kg, DYNAMICALLY
    # fit from all 3 moons' mutual orbits -- Fuksa, Broz, Hanus, Ferrais,
    # Fatka & Vernazza (2023), A&A 677, A189, doi:10.1051/0004-6361/202346386
    # (ADAM shape model: 60 lightcurves + 46 AO images + 2 occultations).
    # GM = 6.67430e-20 * 6.606e18 = 0.4409 km^3/s^2.
    # CAVEAT: Elektra's own shape is SEVERELY non-spherical (ellipsoid
    # 262x205x164 km, volume-equivalent diameter 201 km; shape-derived
    # J2~0.16-0.18, same source) -- an order of magnitude larger than
    # Sylvia's or any planet's J2. Point-mass is used here as a first-pass
    # approximation only (same discipline as every other system in this
    # registry); the moons orbit 5-13x Elektra's ~100 km mean radius out
    # (comparable to or better than the Uranus-Miranda ratio, ~5.1x, which
    # this registry already treats as point-mass), so this is not an
    # unprecedented approximation, but it is NOT validated against Fuksa et
    # al.'s own multipole+mutual-perturbation fit (which was required to
    # match the moons' real orbits precisely).
    "Elektra": 0.4409,
    # (45) Eugenia system GM: mass (5.69+/-0.12)e18 kg, DYNAMICALLY fit from
    # the 2-moon mutual orbit -- Beauvalet & Marchis (2014), Icarus 241, 13,
    # "A Dynamical Solution of the Triple Asteroid System (45) Eugenia"
    # follow-up (refines Marchis et al. 2010, Icarus 210, 635,
    # arXiv:1008.2164). GM = 6.67430e-20 * 5.69e18 = 0.3800 km^3/s^2.
    # CAVEAT: this is the SYSTEM mass (Eugenia + both moons); the moons'
    # individual masses are explicitly stated by Beauvalet & Marchis (2014)
    # to be too small to be constrained by the astrometry, so they are a
    # negligible fraction of this GM (consistent with every other primary in
    # this registry, where PRIMARIES values are system GMs).
    "Eugenia": 0.3800,
    # (216) Kleopatra GM: mass (2.97+/-0.32)e18 kg, from the 2-moon mutual
    # orbit multipole fit -- Marchis & Yang et al. (2021), A&A 653, A57,
    # "(216) Kleopatra, a low density critically rotating M-type asteroid"
    # (arXiv:2108.07207); supersedes the older Descamps et al. (2011,
    # arXiv:1011.5263) mass (4.64e18 kg), ~36% higher, now revised down.
    # GM = 6.67430e-20 * 2.97e18 = 0.19825 km^3/s^2.
    # CAVEAT (the #607 dogbone check): Kleopatra's shape is famously
    # bilobate/"dog-bone" (~270x94x78 km, Shepard et al. 2018, Icarus 311,
    # 197; Broz et al. 2021 companion multipole paper). The moons orbit at
    # only 1.8-2.4x Kleopatra's own half-length (~135 km) -- Broz et al.
    # (2021) needed a multipole expansion to l=10 to match their orbits
    # precisely. Point-mass is used here as a coarse first-pass screening
    # approximation ONLY, per this task's explicit scope; a real J2/C22
    # correction is NOT implemented (flagged as a known limitation, not
    # silently assumed valid).
    "Kleopatra": 0.19825,
}


def _sat(
    name: str,
    primary: str,
    mu_km3_s2: float,
    radius_eq_km: float,
    sma_km: float,
    safe_alt_km: float,
    *,
    retrograde: bool = False,
) -> SatelliteData:
    """Build a SatelliteData with mean motion DERIVED from sma + primary mu.

    ``retrograde`` (#599) defaults to False so every pre-existing call site
    is byte-for-byte unaffected; only Triton passes ``retrograde=True``.
    """
    return SatelliteData(
        name=name,
        code=name,  # full-name scheme
        primary=primary,
        mu_km3_s2=mu_km3_s2,
        radius_eq_km=radius_eq_km,
        sma_km=sma_km,
        mean_motion_deg_day=mean_motion_deg_day_about(sma_km, mu_primary=PRIMARIES[primary]),
        safe_alt_km=safe_alt_km,
        retrograde=retrograde,
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
    # safe_alt 200 km: SOURCED Galilean flyby floor (Campagnola 2014 Europa-tour,
    # docs/notes/2026-06-17-digest-campagnola-2014.md line 201; #428). Europa/Ganymede
    # share the 100 km Campagnola floor (already at 100); Callisto's is explicitly 200.
    "Callisto": _sat("Callisto", "Jupiter", 7179.289, 2410.3, 1882700.0, 200.0),
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
    # safe_alt 50 km for the four tour flyby bodies (Ariel/Umbriel/Titania/Oberon):
    # SOURCED design floor, Heaton-Longuski 2003 Table 4 minimum flyby altitude
    # (docs/notes/2026-06-23-digest-stone-miner-voyager2-uranus-neptune.md; #429).
    # Was an unsourced 100 km convention; lowered to the published design minimum.
    # Miranda stays 100 km (convention): it was a Voyager-2 mass-pass only, not a
    # Heaton-Longuski tour flyby body, so no published design floor exists for it.
    # Miranda: GM 4.3 km^3/s^2, mean R 235.8 km, a 129846 km.
    "Miranda": _sat("Miranda", "Uranus", 4.3, 235.8, 129846.0, 100.0),
    # Ariel: GM 83.5 km^3/s^2, mean R 578.9 km, a 190929 km.
    "Ariel": _sat("Ariel", "Uranus", 83.5, 578.9, 190929.0, 50.0),
    # Umbriel: GM 85.1 km^3/s^2, mean R 584.7 km, a 265986 km.
    "Umbriel": _sat("Umbriel", "Uranus", 85.1, 584.7, 265986.0, 50.0),
    # Titania: GM 226.9 km^3/s^2, mean R 788.9 km, a 436298 km.
    "Titania": _sat("Titania", "Uranus", 226.9, 788.9, 436298.0, 50.0),
    # Oberon: GM 205.3 km^3/s^2, mean R 761.4 km, a 583511 km.
    "Oberon": _sat("Oberon", "Uranus", 205.3, 761.4, 583511.0, 50.0),
    # --- Neptune system (added 2026-06-14) ---
    # JPL SSD satellite physical (phys_par, ref NEP097/NEP101) + mean elements,
    # accessed 2026-06-14. Triton is large (GM ~ 1428) but RETROGRADE + inclined
    # (a ~ 354800 km on a near-circular retrograde orbit) -> hostile flyby
    # geometry; carried for completeness. Proteus added (regular, prograde).
    # Nereid OMITTED: JPL SSD lists its GM as 0.0 (mass not determined), so per
    # the sourcing discipline we omit it rather than guess from a size estimate.
    # Triton: GM 1428.495 km^3/s^2, mean R 1352.6 km, a 354800 km.
    # retrograde=True (#599): the SAME JPL SSD satellite mean-elements source
    # already cited above for Triton's sma/orbital fit (ref NEP097/NEP101)
    # lists Triton's orbital inclination as ~156.885 deg to Neptune's Laplace
    # plane -- i > 90 deg is the standard convention for retrograde orbital
    # motion (it moves opposite the sense of Neptune's rotation and opposite
    # every other regular Neptunian moon, incl. Proteus). This is also the
    # modern-fit reference for the Neptune satellite system: Jacobson, R.A.
    # (2009), "The Orbits of the Neptunian Satellites and the Orientation of
    # the Pole of Neptune", AJ 137, 4322. Triton's retrograde sense itself was
    # first identified by Lassell (1846) from the reversed motion of its
    # apparent orbit; this is long-settled, not a recent-fit artifact.
    "Triton": _sat("Triton", "Neptune", 1428.49546, 1352.6, 354800.0, 100.0, retrograde=True),
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
    # --- Small-body multi-moon systems (#607, added 2026-07-16) ---
    # safe_alt_km = 10 km for every moon below: the SAME engineering-default
    # convention already used for Phobos/Deimos/Amalthea/Hyperion/Nix/Hydra
    # (mean radius < ~140 km) -- NOT a sourced design floor, a stated
    # convention (see the module-level CAVEAT note above SATELLITES). Every
    # GM below is DERIVED (mass * G), not copied from a paper's own "GM"
    # figure -- no source publishes these small-body moons' GM directly.
    #
    # (87) Sylvia's moons -- Fang, Margot & Rojo (2012), AJ 144, 70,
    # "Orbits, masses, and evolution of main belt triple (87) Sylvia"
    # (arXiv:1206.5755): individual moon masses from mutual-perturbation
    # dynamics (a REAL measurement, not a density assumption), though with
    # large (order-100%) asymmetric uncertainties -- the weakest-sourced
    # numbers in this addition. sma from the more recent Vernazza, Carry et
    # al. (2021) VLT/SPHERE fit (A&A 654, A56).
    # Romulus (outer): mass 9.32e14 kg -> GM=6.220e-5 km^3/s^2; mean radius
    # 11.55 km (diameter 23.1+/-0.7 km, Berthier et al. 2014, Icarus 239,
    # 118, itself notably elongated 38.0x14.0 km -- a single mean radius is
    # a simplification); a=1340.6+/-0.4 km, e~0 (near-circular).
    "Romulus": _sat("Romulus", "Sylvia", 6.220e-5, 11.55, 1340.6, 10.0),
    # Remus (inner): mass 7.33e14 kg -> GM=4.893e-5 km^3/s^2; mean radius
    # 3.5 km (diameter ~7+/-2 km, Marchis et al. 2005, Nature 436, 822 --
    # NO updated size measurement found post-2005, the weakest-sourced size
    # in this addition); a=694.2+/-0.1 km, e~0.005 (near-circular).
    "Remus": _sat("Remus", "Sylvia", 4.893e-5, 3.5, 694.2, 10.0),
    #
    # (130) Elektra's moons -- NONE has a published individual GM (Fuksa et
    # al. 2023 give only rough system-mass-fraction estimates, not measured
    # GMs); masses below are ASSUMED-DENSITY estimates (Elektra's own fitted
    # bulk density, 1.536 g/cm^3, same source) from each moon's photometric
    # (albedo-assumed) diameter -- flagged explicitly, matching the
    # #549 Didymos-Dimorphos precedent for an assumed-density-derived mass.
    # No official IAU names exist for any of the 3 -- "Beta"/"Gamma"/"Delta"
    # are the informal press/paper labels (Fuksa et al. 2023), NOT adopted
    # designations; prefixed here to avoid any future registry collision.
    # ElektraBeta (S/2003 (130) 1, outermost, largest, discovered Merline et
    # al. 2003, IAUC 8183): diameter 6.0+/-0.6 km -> mean radius 3.0 km,
    # assumed-density mass ~1.737e14 kg -> GM=1.160e-5 km^3/s^2; a=1297.58+/-
    # 0.54 km, P=5.287d (Kepler-III self-check at this a/mu gives 5.12 d,
    # ~3.3% off -- larger than every other moon in this addition, plausibly
    # real given Elektra's severe J2 and this moon's own e=0.0835,
    # Fuksa et al. 2023).
    "ElektraBeta": _sat("ElektraBeta", "Elektra", 1.160e-5, 3.0, 1297.58, 10.0),
    # ElektraGamma (S/2014 (130) 1, discovered Yang et al., reported Hanus et
    # al. 2017, arXiv:1611.03632): diameter 2.0+/-0.4 km -> mean radius
    # 1.0 km, assumed-density mass ~6.434e12 kg -> GM=4.295e-7 km^3/s^2;
    # a=501+/-7 km (discovery value; Kepler-III self-check against the
    # Fuksa et al. 2023 revised P~1.212 d gives 1.228 d, ~1.3% agreement).
    "ElektraGamma": _sat("ElektraGamma", "Elektra", 4.295e-7, 1.0, 501.0, 10.0),
    # ElektraDelta (S/2014 (130) 2, discovered Berdeu, Langlois & Vachier
    # 2022, A&A 658, L4; discovery orbit found DYNAMICALLY UNSTABLE against
    # Elektra's real gravity field by Valvano et al. 2023, MNRAS 522, 6196,
    # arXiv:2304.14967): diameter 1.6+/-0.4 km -> mean radius 0.8 km,
    # assumed-density mass ~3.294e12 kg -> GM=2.199e-7 km^3/s^2. sma is
    # BACK-DERIVED via Kepler III (a=(GM_Elektra*P^2/(4*pi^2))^(1/3)) from
    # Fuksa et al. (2023)'s own revised period P=1.642112 d -- their revised
    # semi-major axis could not be extracted from the paper (table not
    # machine-readable via the fetch tool used); this is NOT an independent
    # cross-check, it is definitionally self-consistent by construction.
    # Sanity check: 608.1 km falls between ElektraGamma (501 km) and
    # ElektraBeta (1297.58 km), consistent with Fuksa et al.'s own statement
    # that the revised period reorders Delta as the middle moon by period.
    "ElektraDelta": _sat("ElektraDelta", "Elektra", 2.199e-7, 0.8, 608.1, 10.0),
    #
    # (45) Eugenia's moons -- Beauvalet & Marchis (2014) state both moons'
    # individual masses are too small to be constrained by the astrometry;
    # masses below are ASSUMED-DENSITY estimates (Eugenia's own fitted bulk
    # density, 1.69 g/cm^3, same source) from each moon's diameter.
    # PetitPrince (Merline et al. 1999, Nature 401, 565 -- the original
    # discovery, S/1998 (45) 1): diameter ~7+/-2 km -> mean radius 3.5 km,
    # assumed-density mass ~3.035e14 kg -> GM=2.026e-5 km^3/s^2;
    # a=1164.4 km (Beauvalet & Marchis 2014 refined fit), e~0.002-0.01
    # (near-circular).
    "PetitPrince": _sat("PetitPrince", "Eugenia", 2.026e-5, 3.5, 1164.4, 10.0),
    # EugeniaS2 (S/2004 (45) 1, Marchis et al. 2007, IAUC 8817 -- CONFIRMED
    # per JPL SBDB's satellite list, not a retracted candidate; re-fit by
    # Marchis et al. 2010, Icarus 210, 635, arXiv:1008.2164): diameter
    # ~5+/-1 km -> mean radius 2.5 km, assumed-density mass ~1.106e14 kg ->
    # GM=7.383e-6 km^3/s^2; a=610.7 km (midpoint of the 610.6-610.8 km
    # range across cited fits), e~0.07-0.11 (the least circular orbit in
    # this addition besides ElektraBeta/Delta).
    "EugeniaS2": _sat("EugeniaS2", "Eugenia", 7.383e-6, 2.5, 610.7, 10.0),
    #
    # (216) Kleopatra's moons -- Marchis & Brox et al. (2021) companion
    # multipole paper, A&A, "An advanced multipole model for (216)
    # Kleopatra triple system": neither moon's GM is independently measured
    # (both too small to perturb each other or Kleopatra detectably); masses
    # below are ASSUMED-DENSITY estimates from that paper's own diameter +
    # density-model figures, flagged explicitly there as model-dependent.
    # AlexHelios (inner): diameter ~6.9 km -> mean radius 3.45 km, assumed
    # mass ~4e14 kg (assumed density ~2300 kg/m^3) -> GM=2.670e-5 km^3/s^2;
    # a=499 km, P=1.822359+/-0.004156 d, e~0 (circular, equatorial,
    # prograde) -- Kepler-III self-check gives 1.820 d, matches to <0.2%.
    "AlexHelios": _sat("AlexHelios", "Kleopatra", 2.670e-5, 3.45, 499.0, 10.0),
    # CleoSelene (outer): diameter ~8.9 km -> mean radius 4.45 km, assumed
    # mass ~6e14 kg (assumed density ~1600 kg/m^3) -> GM=4.005e-5 km^3/s^2;
    # a=655 km, P=2.745820+/-0.004820 d, e~0 (circular, equatorial,
    # prograde) -- Kepler-III self-check gives 2.738 d, matches to <0.3%.
    "CleoSelene": _sat("CleoSelene", "Kleopatra", 4.005e-5, 4.45, 655.0, 10.0),
}
