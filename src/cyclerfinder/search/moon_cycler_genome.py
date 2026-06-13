"""Repeated-moon multi-revolution cycler genome (Track A, #254).

A SEARCH genome for **repeated-moon multi-revolution cyclers** — the Liang
Callisto-Ganymede-Europa (CGE) class: trajectories that orbit a planet doing
several resonant revolutions between *repeated* encounters with the *same*
moons in a fixed sequence (Liang, Yang, Li, Bai & Qin, "Callisto-Ganymede-
Europa Triple Cyclers", *JGCD* Engineering Note, 2024, DOI 10.2514/1.G008387).

The project's prior genome (zero-rev, single-encounter, no repeated body)
structurally cannot express these, which is exactly why the Jovian region read
EMPTY in the Phase-6 sweeps. This module is the design of
``docs/notes/2026-06-14-repeated-moon-multirev-genome-design.md``, built in the
four bite-sized steps the design calls out:

1. **Moon-system registry + Tisserand/V_inf graph** for a planet's moons,
   reusing the body-agnostic Tisserand machinery
   (:mod:`cyclerfinder.search.tisserand`, ``mu=PRIMARIES[planet]``) and the
   satellite registry (:mod:`cyclerfinder.core.satellites`).
2. **Genome representation + decision vector** — moon sequence, per-leg
   ``(p:q resonance, n_rev)``, per-encounter flyby — with encode/decode
   round-trip and a validity predicate.
3. **Repeated-sequence periodicity corrector** — close the repeated moon
   sequence to periodicity in the planet-centric model, one canonical
   residual (km/s).
4. **REPRODUCE-BEFORE-SEARCH gate** — recover Liang's published CGE members
   from their sourced (V_inf, ToF, phase) data before any search is trusted.

DISCIPLINE: the gate goldens are Liang's PUBLISHED values, never our own
output (the existing same-model reproduction
:mod:`cyclerfinder.search.cge_scaffold` is the trusted reproduction route; this
genome's periodicity residual is checked to AGREE with that route, and the gate
asserts against the published Tables 3/5/7). No catalogue writeback: a closed
search candidate is SILVER, not a validated cycler.

This module deliberately does NOT run the open-ended search (design step 5 =
enumerate-and-close many candidates); that long combinatorial compute belongs
in the #253 discovery-campaign daemon, not a one-shot agent. The bounded
reproduce-before-search gate IS the deliverable here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.tisserand import linkable, tisserand_to_vinf, vinf_to_tisserand

# ---------------------------------------------------------------------------
# Step 1 — moon-system registry + Tisserand / V_inf graph
# ---------------------------------------------------------------------------

DAY_S: float = 86400.0

# Galilean moons of Jupiter, inner-to-outer (the Liang CGE system).
JUPITER_MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")


@dataclass(frozen=True)
class MoonSystem:
    """A planet + its registered moons, the registry side of the genome.

    All physical data is resolved from the in-repo registries
    (:data:`cyclerfinder.core.satellites.PRIMARIES` /
    :data:`~cyclerfinder.core.satellites.SATELLITES`) — never hand-copied. The
    Tisserand/V_inf graph below is centre-aware via ``mu`` = the primary GM, so
    the SAME body-agnostic Tisserand code (#75) serves heliocentric planets and
    Jovicentric/Saturnian moons alike.
    """

    planet: str
    moons: tuple[str, ...]

    @property
    def mu_primary(self) -> float:
        """Primary (planet-system) GM in km^3/s^2 (registry value)."""
        return PRIMARIES[self.planet]

    def sma_km(self, moon: str) -> float:
        """Semi-major axis of ``moon`` about the primary (km, registry)."""
        sat = SATELLITES[moon]
        if sat.primary != self.planet:
            raise ValueError(f"{moon} orbits {sat.primary}, not {self.planet}")
        return sat.sma_km

    def period_days(self, moon: str) -> float:
        """Sidereal orbital period of ``moon`` about the primary (days)."""
        a = self.sma_km(moon)
        return 2.0 * math.pi * math.sqrt(a**3 / self.mu_primary) / DAY_S

    def circular_speed_kms(self, moon: str) -> float:
        """Circular orbital speed of ``moon`` about the primary (km/s)."""
        return math.sqrt(self.mu_primary / self.sma_km(moon))


def jupiter_system() -> MoonSystem:
    """The Jupiter / Galilean-moon system (the Liang CGE validation target)."""
    return MoonSystem(planet="Jupiter", moons=JUPITER_MOONS)


def moon_vinf_to_tisserand(system: MoonSystem, moon: str, vinf_kms: float) -> float:
    """Tisserand parameter at ``moon`` for planet-centric V_inf (km/s).

    Thin centre-aware wrapper over :func:`cyclerfinder.search.tisserand.
    vinf_to_tisserand` with ``mu`` = the primary GM and ``a_p`` = the moon SMA.
    ``T = 3 - V_inf^2 a_moon / mu_primary``.
    """
    return vinf_to_tisserand(moon, vinf_kms, mu=system.mu_primary)


def moon_tisserand_to_vinf(system: MoonSystem, moon: str, t_p: float) -> float:
    """Inverse of :func:`moon_vinf_to_tisserand` (planet-centric V_inf, km/s)."""
    return tisserand_to_vinf(moon, t_p, mu=system.mu_primary)


def _moon_a_range_au(system: MoonSystem) -> tuple[float, float]:
    """``(a_min, a_max)`` in AU bracketing the system's moon SMAs.

    The body-agnostic Tisserand contour code works in AU internally
    (``a_p_au = a_p_km / AU_KM``); for a moon system the spacecraft semi-major
    axis lives on the moon SMA scale, so the AU range must be moon-scale (a few
    1e-3 AU for the Galileans), NOT the default heliocentric ``(0.3, 5.0)``.
    The window spans 0.4x the innermost to 2.5x the outermost moon SMA so a
    resonant spacecraft arc between any two moons is inside it.
    """
    from cyclerfinder.core.constants import AU_KM

    smas_au = [system.sma_km(m) / AU_KM for m in system.moons]
    return (0.4 * min(smas_au), 2.5 * max(smas_au))


def moon_linkable(
    system: MoonSystem,
    moon_a: str,
    moon_b: str,
    vinf_kms: float,
    *,
    n_points: int = 200,
) -> bool:
    """Do ``moon_a`` and ``moon_b`` share a constant-V_inf spacecraft orbit?

    The Tisserand-graph edge predicate for a moon tour: True iff a single
    planet-centric spacecraft ``(a, e)`` is reachable from both moons at the
    common ``vinf_kms`` (so a ballistic leg between them at that energy is
    energetically possible). Reuses the body-agnostic
    :func:`cyclerfinder.search.tisserand.linkable` with the primary GM and the
    moon-scale ``a_range_au``.
    """
    return linkable(
        moon_a,
        moon_b,
        vinf_kms,
        a_range_au=_moon_a_range_au(system),
        n_points=n_points,
        mu=system.mu_primary,
    )


def vinf_graph_edges(
    system: MoonSystem,
    vinf_kms: float,
    *,
    n_points: int = 200,
) -> dict[tuple[str, str], bool]:
    """All ordered moon-pair link flags at a common ``vinf_kms``.

    The Tisserand/V_inf graph at one energy level: ``{(m_i, m_j): linkable}``
    over distinct moon pairs (self-pairs are the resonant-return legs handled
    by the resonance machinery, not by the cross-moon link predicate, so they
    are omitted). Symmetric by construction; both orderings are stored for
    convenient lookup.
    """
    edges: dict[tuple[str, str], bool] = {}
    moons = system.moons
    for i, m_i in enumerate(moons):
        for m_j in moons[i + 1 :]:
            ok = moon_linkable(system, m_i, m_j, vinf_kms, n_points=n_points)
            edges[(m_i, m_j)] = ok
            edges[(m_j, m_i)] = ok
    return edges
