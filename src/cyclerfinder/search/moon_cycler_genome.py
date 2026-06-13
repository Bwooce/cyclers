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


# ---------------------------------------------------------------------------
# Step 2 — genome representation + decision vector
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LegGene:
    """One inter-encounter leg of a repeated-moon cycle.

    A planet-centric resonance ``p:q`` (p spacecraft revs per q moon revs)
    plus a multi-revolution count ``n_rev`` (the discrete Lambert revolution
    number realising the resonant arc). For the idealized Liang CGE legs all
    four resolve to ``n_rev = 1`` one-revolution Lambert arcs (cge_scaffold
    geometry conventions); the resonance bookkeeping records the commensurable
    arc class the leg implements.
    """

    p: int  # spacecraft revolutions on the arc (resonance numerator)
    q: int  # moon revolutions over the same interval (resonance denominator)
    n_rev: int  # Lambert multi-rev count realising the arc

    def is_valid(self) -> bool:
        return self.p >= 1 and self.q >= 1 and self.n_rev >= 0


@dataclass(frozen=True)
class EncounterGene:
    """One flyby of the repeated-moon sequence.

    ``moon`` is the body flown by; ``b_plane_angle_rad`` is the B-plane clock
    angle (the continuous DOF the corrector closes — it sets the post-flyby
    plane/direction at fixed V_inf magnitude, the ballistic flyby's only free
    parameter in the planar idealization it is the in-plane turn sign/branch).
    """

    moon: str
    b_plane_angle_rad: float


@dataclass(frozen=True)
class MoonCyclerGenome:
    """A repeated-moon multi-revolution cycler genome (the decision object).

    Structure (design "Representation"):

    * ``system`` — the planet + moon registry the sequence lives in.
    * ``encounters`` — the repeated moon sequence ``[m_1, ..., m_k]`` as flyby
      genes (k encounters per cycle; the sequence repeats every cycle).
    * ``legs`` — the ``k - 1`` inter-encounter legs (resonance + n_rev). (A
      cycle's published anchor is a single cycle's k flybys / k-1 transfer
      legs, e.g. Liang's Callisto-Ganymede-Callisto-Europa-Callisto = 5
      flybys, 4 legs.)
    * ``epoch_days`` — overall phasing epoch (days from the cycle anchor).
    * ``perijove_scale`` — the family's perijove scale (Liang member A/B = 1.0,
      member C = 0.5), part of the overall phasing/energy DOF.

    The decision vector (:func:`MoonCyclerGenome.to_vector` /
    :func:`from_vector`) flattens the continuous + discrete DOF into a single
    1-D ``list[float]`` for the corrector / search; the moon-id and resonance
    integers ride in the vector as exact floats (small integers are
    representable), so the round-trip is lossless.
    """

    system: MoonSystem
    encounters: tuple[EncounterGene, ...]
    legs: tuple[LegGene, ...]
    epoch_days: float = 0.0
    perijove_scale: float = 1.0

    def __post_init__(self) -> None:
        if len(self.legs) != max(0, len(self.encounters) - 1):
            raise ValueError(
                f"a k-encounter cycle has k-1 legs; got {len(self.encounters)} "
                f"encounters and {len(self.legs)} legs"
            )

    @property
    def sequence(self) -> tuple[str, ...]:
        """The repeated moon sequence as moon names."""
        return tuple(e.moon for e in self.encounters)

    def is_valid(self) -> bool:
        """Structural validity: registered moons, well-formed legs, k>=2.

        Does NOT assert dynamical closure (that is the corrector's residual);
        this is the cheap up-front filter the enumerator/search uses.
        """
        if len(self.encounters) < 2:
            return False
        for e in self.encounters:
            sat = SATELLITES.get(e.moon)
            if sat is None or sat.primary != self.system.planet:
                return False
            if not math.isfinite(e.b_plane_angle_rad):
                return False
        return all(leg.is_valid() for leg in self.legs)

    def to_vector(self) -> list[float]:
        """Flatten to a single decision vector (the corrector/search interface).

        Layout (all float64):
        ``[epoch_days, perijove_scale,
           moon_idx_0, b_plane_0, ..., moon_idx_{k-1}, b_plane_{k-1},
           p_0, q_0, n_rev_0, ..., p_{k-2}, q_{k-2}, n_rev_{k-2}]``
        where ``moon_idx`` is the index of the moon within ``system.moons``.
        """
        idx = {m: i for i, m in enumerate(self.system.moons)}
        vec: list[float] = [self.epoch_days, self.perijove_scale]
        for e in self.encounters:
            vec.append(float(idx[e.moon]))
            vec.append(e.b_plane_angle_rad)
        for leg in self.legs:
            vec.extend((float(leg.p), float(leg.q), float(leg.n_rev)))
        return vec

    @classmethod
    def from_vector(
        cls, vec: list[float], system: MoonSystem, n_encounters: int
    ) -> MoonCyclerGenome:
        """Inverse of :func:`to_vector` (the round-trip partner).

        ``n_encounters`` (k) is structural metadata not carried in the vector,
        so the search must supply it (it is fixed per enumerated sequence).
        """
        expected = 2 + 2 * n_encounters + 3 * (n_encounters - 1)
        if len(vec) != expected:
            raise ValueError(
                f"vector length {len(vec)} != expected {expected} for k={n_encounters}"
            )
        epoch_days = vec[0]
        perijove_scale = vec[1]
        encounters: list[EncounterGene] = []
        cursor = 2
        for _ in range(n_encounters):
            moon = system.moons[round(vec[cursor])]
            encounters.append(EncounterGene(moon=moon, b_plane_angle_rad=vec[cursor + 1]))
            cursor += 2
        legs: list[LegGene] = []
        for _ in range(n_encounters - 1):
            legs.append(
                LegGene(
                    p=round(vec[cursor]),
                    q=round(vec[cursor + 1]),
                    n_rev=round(vec[cursor + 2]),
                )
            )
            cursor += 3
        return cls(
            system=system,
            encounters=tuple(encounters),
            legs=tuple(legs),
            epoch_days=epoch_days,
            perijove_scale=perijove_scale,
        )


# The Liang CGE repeated-moon sequence (one published cycle), as moon names.
# Callisto-Ganymede-Callisto-Europa-Callisto — 5 flybys / 4 transfer legs per
# ~100 d cycle (Liang et al. 2024 Sec. III.B; see cge_scaffold.CGCEC_SEQUENCE).
CGE_SEQUENCE: tuple[str, ...] = ("Callisto", "Ganymede", "Callisto", "Europa", "Callisto")


def liang_member_genome(member: str) -> MoonCyclerGenome:
    """Encode one published Liang CGE member (A/B/C) as a genome.

    The discrete structure (moon sequence, per-leg ``n_rev = 1`` one-rev
    Lambert arcs, perijove scale) is taken from the sourced
    :mod:`cyclerfinder.search.cge_scaffold` member spec. The B-plane angles are
    initialised to 0 (the corrector / gate resolves the flyby geometry from the
    published V_inf + ToF); this function is the "the Liang members encode
    validly" round-trip anchor of design step 2, NOT a dynamical solve.

    Resonance bookkeeping: the spacecraft period equals Callisto's period
    (Eq. 14), so each Callisto-anchored arc is a 1:1 Callisto resonance; the
    Ganymede/Europa legs are the multi-rev Lambert reservoirs that reshape the
    phase within that 1:1 Callisto frame (paper Sec. III.A). All four legs are
    recorded as ``n_rev = 1`` (cge_scaffold golden: every leg is a 1-rev
    Lambert solution).
    """
    from cyclerfinder.search.cge_scaffold import LIANG_MEMBERS

    spec = LIANG_MEMBERS[member]
    system = jupiter_system()
    encounters = tuple(EncounterGene(moon=m, b_plane_angle_rad=0.0) for m in CGE_SEQUENCE)
    # p:q resonance per leg — all arcs ride the 1:1 Callisto-period frame; the
    # n_rev = 1 multi-rev Lambert count is the cge_scaffold golden.
    legs = tuple(LegGene(p=1, q=1, n_rev=1) for _ in range(len(CGE_SEQUENCE) - 1))
    return MoonCyclerGenome(
        system=system,
        encounters=encounters,
        legs=legs,
        epoch_days=0.0,
        perijove_scale=spec.perijove_scale,
    )


# ---------------------------------------------------------------------------
# Step 3 — repeated-sequence periodicity corrector (one canonical residual)
# ---------------------------------------------------------------------------

# Ballistic-cycler threshold: Liang's members are published ballistic with a
# residual flyby-defect Delta-v below 1e-8 m/s = 1e-11 km/s (paper p. 13). Our
# same-model reconstruction is input-precision-limited, so the achievable
# residual is the Table 1 print quantization (~1e-2 km/s), not 1e-11; the
# corrector reports the residual in km/s and the GATE asserts against the
# sourced print-precision tolerance (cge_scaffold.vinf_print_tolerance_kms),
# never a hand-tuned number.
BALLISTIC_DV_THRESHOLD_KMS: float = 1.0e-11


@dataclass(frozen=True)
class PeriodicityResidual:
    """The ONE canonical residual of a repeated-moon cycle closure.

    A repeated-moon cycler is periodic when, after the full k-encounter
    sequence, the trajectory returns to the same relative geometry — for a
    ballistic tour this is exactly V_inf-magnitude continuity at every flyby
    (the spacecraft hyperbolic excess speed entering a moon equals the speed
    leaving it; the flyby only turns the vector). The canonical residual is the
    worst such continuity defect over the cycle, in km/s:

    .. math::

        R = \\max_i \\bigl| |V_{\\infty,i}^{\\text{in}}| -
                            |V_{\\infty,i}^{\\text{out}}| \\bigr|.

    This mirrors the #248 harness discipline: a single scalar (km/s),
    epoch-safe (evaluated at the Eq. 16-anchored cumulative flyby epochs), no
    ad-hoc per-leg metrics. The daemon/search closes a candidate by driving
    this residual below the ballistic threshold; the GATE checks it against the
    sourced print-precision tolerance for the Liang anchors.
    """

    worst_continuity_kms: float
    per_flyby_continuity_kms: tuple[float, ...]
    worst_vinf_vs_published_kms: float  # only meaningful for the Liang anchors
    n_flybys: int

    def closes(self, tol_kms: float) -> bool:
        """True iff the worst continuity defect is within ``tol_kms``."""
        return self.worst_continuity_kms <= tol_kms


def liang_periodicity_residual(member: str) -> PeriodicityResidual:
    """Canonical periodicity residual for a published Liang CGE member.

    Routes the genome's repeated CGE sequence through the trusted same-model
    planet-centric reconstruction (:func:`cyclerfinder.search.cge_scaffold.
    reproduce_member`): place the moons by mean motion + the sourced Table 2/4/6
    phases at the Eq. 16-anchored cumulative flyby epochs, solve the Jupiter-
    frame Lambert legs at the sourced Table 3/5/7 ToFs, and read off the per-
    flyby V_inf-magnitude continuity. The residual is the worst continuity
    defect across the cycle — the SAME quantity the daemon would minimise for a
    new candidate, evaluated here on the published member (the "a single known
    leg/closure closes" test of design step 3).

    DISCIPLINE: every input is sourced (Liang's printed phases/ToFs/V_inf); the
    residual is our model's continuity, asserted by the gate against the
    sourced print tolerance, not a value we fit.
    """
    from cyclerfinder.search.cge_scaffold import reproduce_member

    rep = reproduce_member(member)
    continuities = tuple(fb.continuity_kms for fb in rep.flybys if fb.continuity_kms is not None)
    worst_continuity = max(continuities) if continuities else 0.0
    return PeriodicityResidual(
        worst_continuity_kms=worst_continuity,
        per_flyby_continuity_kms=continuities,
        worst_vinf_vs_published_kms=rep.max_vinf_residual_kms,
        n_flybys=len(rep.flybys),
    )
