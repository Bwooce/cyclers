"""BCR4BP per-system constants registry (#334 Part A).

Sourced (``primary``-``secondary`` GM + SMA) parameters for building a
:class:`cyclerfinder.core.bcr4bp.BCR4BPSystem` for any Sun-primary-secondary
bicircular triple. Used by ``scripts/scan_334_bcr4bp_system_swap.py`` to
extend the #303/#304/#313 sweep ("does the SEM BCR4BP family pattern carry
over to other moon systems?") to the rest of the planetary-moon catalogue,
and to fit the geometric scaling rule

    Δx0_target / Δx0_SEM  ≈  (mu_sun_target / mu_sun_SEM)
                              x (a_sun_SEM / a_sun_target)^k

predicted by the #326 structural finding (commit ``c1896ef``).

Each entry's ``mu`` / ``mu_sun`` / ``a_sun_nondim`` / ``omega_sun_nondim``
is COMPUTED from sourced physical inputs (NASA JPL SSD GM tables / IAU
exact AU / heliocentric SMA from Standish & Williams 1992 used elsewhere
in the codebase) per the formulas in the module docstring; the ``sources``
field documents the provenance per entry.

The Sun-Earth-Moon entry SEM_ANDREU pins the four exact Rosales-Jorba 2023
Table 3 constants used throughout #292/#303/#304 (so the registry's SEM
record matches the published / built-in :func:`bcr4bp.andreu_default` exactly
to floating-point precision). The "derived" SEM record SEM_DERIVED uses the
same registry formulas as every other system and is kept for cross-check:
the two SEM records agree to ~4 decimal places in mu_sun / a_sun (model
difference between Andreu's exact Earth-Moon barycenter SMA and the
Standish & Williams J2000 Earth SMA / observed Moon SMA used elsewhere
in the codebase).

Formulas (mass-ratio + length-ratio bicircular triple, see Simo-Jorba-
Gomez 1995 / Andreu 1998 / Gimeno-Jorba 2018):

  mu          = GM_secondary / GM_primary_system
  mu_sun      = GM_sun / GM_primary_system
  a_sun_nondim = a_primary_around_sun / a_secondary_around_primary
  TU_seconds   = sqrt(a_secondary**3 / GM_primary_system)
  n_primary    = sqrt(GM_sun / a_primary_around_sun**3)
  omega_sun    = 1.0 - n_primary * TU_seconds

Discipline
----------
  * Sourced constants only -- every value cites a JPL SSD / DE440 / IAU
    upstream the codebase already uses elsewhere
    (``feedback_golden_tests_sourced_only``).
  * No catalogue writeback. No novelty claims.
  * READ-ONLY on ``bcr4bp.py`` / ``bcr4bp_genome.py`` / ``bcr4bp_continuation.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

__all__ = [
    "REGISTRY",
    "SEM_ANDREU",
    "SEM_DERIVED",
    "SJE_DERIVED",
    "SJI_DERIVED",
    "SMP_DERIVED",
    "SNT_DERIVED",
    "SPC_DERIVED",
    "SSE_DERIVED",
    "SST_DERIVED",
    "BCR4BPSystemConstants",
    "build_bcr4bp_system",
    "derive_from_sources",
]


# ---------------------------------------------------------------------------
# Dataclass.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BCR4BPSystemConstants:
    """Sourced constants for a Sun-primary-secondary BCR4BP triple.

    Attributes
    ----------
    name :
        Human-readable identifier, e.g. ``"Sun-Saturn-Titan"``.
    primary :
        Primary body name (e.g. ``"Saturn"``); a key in
        :data:`cyclerfinder.core.satellites.PRIMARIES` for GM lookup.
    secondary :
        Secondary body name (e.g. ``"Titan"``); a key in
        :data:`cyclerfinder.core.satellites.SATELLITES` (or, for the
        Earth-Moon case, a member-specific override -- the Moon entry exists
        in SATELLITES).
    mu :
        Secondary / primary-system mass ratio (BCR4BP definition).
    mu_sun :
        Sun / primary-system mass ratio (BCR4BP definition).
    a_sun_nondim :
        Heliocentric primary SMA in primary-secondary distance units.
    omega_sun_nondim :
        Sun synodic angular frequency in the primary-secondary rotating frame
        (rad / TU).
    tu_seconds :
        Primary-secondary nondim time unit in seconds (Kepler III).
    l_km :
        Primary-secondary distance in km (the BCR4BP length unit).
    a_primary_au :
        Heliocentric SMA of the primary, in AU (provenance / display).
    sources :
        Tuple of provenance strings: JPL SSD GM tables, IAU AU, Standish &
        Williams 1992 heliocentric SMAs, etc. Every numerical field above
        traces to one or more of these upstreams.
    """

    name: str
    primary: str
    secondary: str
    mu: float
    mu_sun: float
    a_sun_nondim: float
    omega_sun_nondim: float
    tu_seconds: float
    l_km: float
    a_primary_au: float
    sources: tuple[str, ...]


# ---------------------------------------------------------------------------
# Internal helpers.
# ---------------------------------------------------------------------------


# Heliocentric SMAs of the primaries (AU). Drawn from ``core.constants.PLANETS``
# where available; Earth/Mars/Jupiter/Saturn/Uranus/Neptune are present there.
# Pluto's heliocentric SMA is taken from ``core.constants._PLUTO_SMA_AU`` (in
# PLANETS as a dwarf-planet entry) via the same Standish & Williams 1992 table.
_PRIMARY_SMA_AU: dict[str, float] = {
    "Earth": PLANETS["E"].sma_au,
    "Mars": PLANETS["M"].sma_au,
    "Jupiter": PLANETS["J"].sma_au,
    "Saturn": PLANETS["S"].sma_au,
    "Uranus": PLANETS["U"].sma_au,
    "Neptune": PLANETS["N"].sma_au,
    "Pluto": PLANETS["Pl"].sma_au,
}


_PROVENANCE_GM = (
    "JPL SSD gm_de440 planetary-system / satellite GM tables "
    "(http://ssd.jpl.nasa.gov/, accessed 2026-06-07)"
)
_PROVENANCE_SUN = (
    "MU_SUN_KM3_S2 = 1.32712440018e11 km^3/s^2 (IAU 2015 nominal, cyclerfinder.core.constants)"
)
_PROVENANCE_AU = "AU_KM = 1.49597870700e8 km (IAU 2012 Resolution B2, exact)"
_PROVENANCE_PSMA = (
    "Heliocentric primary SMAs from Standish & Williams 1992 "
    "'Approximate Positions of the Planets' (J2000), as cached in "
    "cyclerfinder.core.constants.PLANETS"
)
_PROVENANCE_SECSMA = (
    "Secondary SMAs and GMs from cyclerfinder.core.satellites.SATELLITES "
    "(JPL SSD satellite physical/orbital tables, accessed 2026-06-07/14)"
)


def derive_from_sources(
    *,
    name: str,
    primary: str,
    secondary: str,
    extra_sources: tuple[str, ...] = (),
) -> BCR4BPSystemConstants:
    """Build a :class:`BCR4BPSystemConstants` from the codebase's sourced inputs.

    Pulls primary GM from :data:`cyclerfinder.core.satellites.PRIMARIES`,
    secondary GM + SMA from :data:`cyclerfinder.core.satellites.SATELLITES`,
    primary heliocentric SMA from :data:`cyclerfinder.core.constants.PLANETS`,
    and computes the four BCR4BP nondim constants per the docstring formulas.
    The returned record is "derived" (formulas) -- the SEM Andreu record is a
    sourced override that matches Rosales-Jorba 2023 Table 3 exactly.
    """
    if primary not in PRIMARIES:
        raise KeyError(f"unknown primary {primary!r}; not in PRIMARIES")
    if secondary not in SATELLITES:
        raise KeyError(f"unknown secondary {secondary!r}; not in SATELLITES")
    sat = SATELLITES[secondary]
    if sat.primary != primary:
        raise ValueError(
            f"{secondary} orbits {sat.primary!r}, not {primary!r}; registry primary mismatch"
        )
    if primary not in _PRIMARY_SMA_AU:
        raise KeyError(f"no heliocentric SMA known for primary {primary!r}")

    gm_primary_sys = PRIMARIES[primary]
    gm_secondary = sat.mu_km3_s2
    l_km = sat.sma_km
    a_primary_au = _PRIMARY_SMA_AU[primary]
    a_primary_km = a_primary_au * AU_KM

    mu = gm_secondary / gm_primary_sys
    mu_sun = MU_SUN_KM3_S2 / gm_primary_sys
    a_sun_nondim = a_primary_km / l_km
    tu_seconds = math.sqrt(l_km**3 / gm_primary_sys)
    n_primary = math.sqrt(MU_SUN_KM3_S2 / a_primary_km**3)
    omega_sun_nondim = 1.0 - n_primary * tu_seconds

    sources = (
        _PROVENANCE_GM,
        _PROVENANCE_SUN,
        _PROVENANCE_AU,
        _PROVENANCE_PSMA,
        _PROVENANCE_SECSMA,
        *extra_sources,
    )
    return BCR4BPSystemConstants(
        name=name,
        primary=primary,
        secondary=secondary,
        mu=mu,
        mu_sun=mu_sun,
        a_sun_nondim=a_sun_nondim,
        omega_sun_nondim=omega_sun_nondim,
        tu_seconds=tu_seconds,
        l_km=l_km,
        a_primary_au=a_primary_au,
        sources=sources,
    )


def build_bcr4bp_system(
    consts: BCR4BPSystemConstants, *, mu_sun_override: float | None = None
) -> bcr4bp.BCR4BPSystem:
    """Build a :class:`bcr4bp.BCR4BPSystem` from a registry entry.

    Pass ``mu_sun_override=0.0`` to construct the CR3BP-limit anchor for the
    sourced-golden round-trip test, or any intermediate value for a specific
    point in a mu_sun continuation.
    """
    mu_sun = consts.mu_sun if mu_sun_override is None else float(mu_sun_override)
    return bcr4bp.BCR4BPSystem(
        mu=consts.mu,
        mu_sun=mu_sun,
        a_sun_nondim=consts.a_sun_nondim,
        omega_sun_nondim=consts.omega_sun_nondim,
    )


# ---------------------------------------------------------------------------
# Sourced registry entries.
# ---------------------------------------------------------------------------


# SEM Andreu / Rosales-Jorba 2023 Table 3 -- the canonical baseline used by
# #292/#303/#304. These are PUBLISHED constants pinned exactly so the SEM
# registry record matches bcr4bp.andreu_default() to floating-point identity.
SEM_ANDREU: BCR4BPSystemConstants = BCR4BPSystemConstants(
    name="Sun-Earth-Moon (Andreu / Rosales-Jorba 2023 Table 3)",
    primary="Earth",
    secondary="Moon",
    mu=0.012150581600000,
    mu_sun=328900.5423094043,
    a_sun_nondim=388.8111430233511,
    omega_sun_nondim=0.925195985520347,
    tu_seconds=375190.0,  # ~ EM mean lunar sidereal frequency, matches scripts/run_303
    l_km=384400.0,  # JPL SSD Moon SMA
    a_primary_au=1.00000261,
    sources=(
        "Rosales-Jorba 2023 Table 3 BCR4BP parameters (mu_S, a_S, omega_S)",
        "Gimeno-Jorba 2018 Table 3 (Andreu 1998 published constants)",
        _PROVENANCE_GM,
        _PROVENANCE_SUN,
        _PROVENANCE_AU,
    ),
)


# Derived records: each computed by ``derive_from_sources``. The dynamic
# constants ARE sourced -- every formula input is a JPL/IAU value present
# elsewhere in the codebase -- but the four BCR4BP nondim values are OUR
# computation rather than the digest-published ones (which only exist for SEM).
SEM_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Earth-Moon (derived from codebase sources, registry self-consistency)",
    primary="Earth",
    secondary="Moon",
    extra_sources=(
        "Cross-check vs SEM_ANDREU: agreement to ~4 decimal places in mu_sun/a_sun "
        "reflects Andreu's exact Earth-Moon-barycenter SMA vs the J2000 Standish-Williams "
        "Earth SMA and observed Moon SMA used elsewhere in the codebase.",
    ),
)


# Sun-Jupiter-Europa: covered by #313, included for regression cross-check
# against scan_313_sun_jupiter_europa.jsonl.
SJE_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Jupiter-Europa",
    primary="Jupiter",
    secondary="Europa",
    extra_sources=(
        "Regression cross-check against scan_313_sun_jupiter_europa.jsonl (#313 Part B)",
    ),
)

# Sun-Jupiter-Io: also covered by #313, regression cross-check.
SJI_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Jupiter-Io",
    primary="Jupiter",
    secondary="Io",
    extra_sources=("Regression cross-check against scan_313_sun_jupiter_io.jsonl (#313 Part B)",),
)


# Sun-Saturn-Titan: new. Titan is the largest Saturnian moon (GM ~ 8978),
# moderate mu ~ 2.37e-4 -- the closest analogue to the EM mass ratio among
# the moon systems in the codebase.
SST_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Saturn-Titan",
    primary="Saturn",
    secondary="Titan",
)


# Sun-Saturn-Enceladus: small inner Saturnian moon (GM ~ 7.21), tiny mu ~ 1.9e-7.
SSE_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Saturn-Enceladus",
    primary="Saturn",
    secondary="Enceladus",
)


# Sun-Mars-Phobos: Phobos is a tiny irregular (GM ~ 7e-4), very low mu ~ 1.65e-8.
# Carried here for scope -- the screen will show whether the L1 Lyapunov family
# even exists at such low mu, and Phobos sits deep in Mars's Hill sphere so the
# Sun perturbation should be very weak (a_sun ~ 24300 -- huge).
SMP_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Mars-Phobos",
    primary="Mars",
    secondary="Phobos",
)


# Sun-Neptune-Triton: Triton is large (GM ~ 1428), mu ~ 2.1e-4 -- comparable to EM,
# but the orbit is RETROGRADE + inclined (hostile real-system geometry). The
# planar prograde BCR4BP doesn't honour the retrograde inclination -- the family
# we compute is the planar idealisation. Carried for the regime scan.
SNT_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Neptune-Triton",
    primary="Neptune",
    secondary="Triton",
    extra_sources=(
        "CAVEAT: real Triton orbit is retrograde + inclined; the planar BCR4BP "
        "model carried here is an idealisation that does not reflect real-system "
        "dynamics. Useful only for the geometric scaling-rule fit.",
    ),
)


# Sun-Pluto-Charon: large mass ratio (mu ~ 0.108) -- a genuine binary, well above
# the small-mu regime where the L1 Lyapunov family is well-defined. Carried for
# the registry but expected to fail the CR3BP-limit anchor if mu is too far
# from the L1 Lyapunov / small-mu basin.
SPC_DERIVED: BCR4BPSystemConstants = derive_from_sources(
    name="Sun-Pluto-Charon",
    primary="Pluto",
    secondary="Charon",
    extra_sources=(
        "CAVEAT: mu ~ 0.108 is a genuine binary (not small-mu); the L1 Lyapunov "
        "perpendicular-crossing seed strategy may not converge. Carried for "
        "scope; structural-failure mode is acceptable per the orbit-closure "
        "discipline.",
    ),
)


REGISTRY: tuple[BCR4BPSystemConstants, ...] = (
    SEM_ANDREU,
    SJE_DERIVED,
    SJI_DERIVED,
    SST_DERIVED,
    SSE_DERIVED,
    SMP_DERIVED,
    SNT_DERIVED,
    SPC_DERIVED,
)
"""Tuple of registry entries actually swept by ``scripts/scan_334_*``.

The SEM Andreu record is the geometric-scaling-rule anchor (the published
Δx0 ~ 1.055e-4 is associated with these EXACT constants). The SEM derived
record is constructed but kept out of the sweep tuple -- it's only used for
the self-consistency cross-check.
"""
