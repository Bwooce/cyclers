"""Composed two-moon Keplerian-map screen (#688, Stage A of #686's CCR4BP plan).

This composes the existing, positive-controlled Ross-Scheeres (2007) exterior
Keplerian map (`keplerian_map.py`, `#500`) into an ALTERNATING two-map system
-- a Jupiter-Europa map and a Jupiter-Ganymede map used in alternation -- and
searches the composed map for periodic itineraries: sequences that alternate
between Europa-map and Ganymede-map segments and (approximately) close on
themselves, phase-locked to the moons' real ~2:1 Laplace commensurability.

SCREEN-GRADE ONLY (carried verbatim from `#686`'s strategy note):
    Each single-moon map ignores the OTHER moon's gravity while it is "in
    control".  Anything found here is SEED GEOMETRY for the real Stage-B CCR4BP
    build, not a registry-grade result on its own.  A negative here is NOT
    registry-grade (do not stamp `empty_regions.jsonl`); a positive is a SEED,
    not a result.

Physics of the composition
--------------------------
The RS07 map is the EXTERIOR periapsis map: it advances a spacecraft whose
periapsis sits just OUTSIDE a moon's circular orbit by one periapsis passage,
applying the moon's gravitational energy kick.  Validity requires the
spacecraft periapsis to be outside the moon, ``a_norm * (1 - e) >= 1`` in that
moon's length units (moon SMA = 1).

To compose the Europa and Ganymede maps, the physical (Jupiter-centred) orbit
is the shared invariant; only the NORMALISATION changes between segments.  At a
patch the physical semimajor axis ``a_phys = a_norm * SMA_moon`` and the
eccentricity ``e`` are preserved and re-expressed in the target moon's units,
and the target-frame Tisserand constant ``C_J`` is recomputed from the target
``a_norm`` and the preserved ``e`` (the Tisserand parameter w.r.t. each moon
genuinely differs for the same physical orbit).  Absolute time is tracked so
each moon's inertial phase -- and hence the periapsis-encounter geometry -- is
self-consistent across the whole itinerary (the per-encounter self-consistency
discipline: an analytic energy match is NOT sufficient; the moon must actually
be near the spacecraft periapsis for the kick to be physical).

References
----------
Ross & Scheeres, SIADS 6(3):576-596, 2007 (the single-moon map; `#500` PC0).
Kumar, Anderson, de la Llave, Gunter, AAS 21-651 / arXiv:2109.14815 (the CCR4BP
    Europa/Ganymede resonant-orbit target this screen pre-screens for Stage B).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.genome.keplerian_map import (
    KeplerianMap,
    eccentricity_from_tisserand,
)

_TWO_PI: float = 2.0 * math.pi
_GM_JUPITER: float = PRIMARIES["Jupiter"]  # km^3/s^2 (JPL SSD gm_de440 system GM)


# ---------------------------------------------------------------------------
# Per-moon configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MoonMapConfig:
    """Physical + map parameters for one Jupiter-moon Keplerian map.

    All physical constants trace to `core.satellites` (JPL SSD), never to
    values computed by this project's own code.
    """

    name: str
    mu: float  # perturbing-moon mass ratio m_moon/(m_primary+m_moon)
    sma_km: float  # moon circular-orbit radius about Jupiter [km]
    gm_moon_km3_s2: float

    @property
    def period_s(self) -> float:
        """Moon orbital period about Jupiter [s]."""
        return _TWO_PI * math.sqrt(self.sma_km**3 / _GM_JUPITER)

    @property
    def mean_motion_rad_s(self) -> float:
        """Moon mean motion [rad/s]."""
        return math.sqrt(_GM_JUPITER / self.sma_km**3)

    @property
    def orbital_speed_km_s(self) -> float:
        """Moon circular orbital speed [km/s] (for physical dv conversions)."""
        return math.sqrt(_GM_JUPITER / self.sma_km)


def moon_config(name: str) -> MoonMapConfig:
    """Build a MoonMapConfig for a Jupiter moon from `core.satellites`."""
    sat = SATELLITES[name]
    if sat.primary != "Jupiter":
        raise ValueError(f"{name} does not orbit Jupiter")
    mu = sat.mu_km3_s2 / (_GM_JUPITER + sat.mu_km3_s2)
    return MoonMapConfig(
        name=name,
        mu=mu,
        sma_km=sat.sma_km,
        gm_moon_km3_s2=sat.mu_km3_s2,
    )


# ---------------------------------------------------------------------------
# Resonance / geometry helpers (moon-normalised units, SMA_moon = 1)
# ---------------------------------------------------------------------------


def resonance_semimajor(p_sc: int, q_moon: int) -> float:
    """Semimajor axis (moon units) of a spacecraft:moon = p:q resonance.

    ``p_sc`` spacecraft revolutions per ``q_moon`` moon revolutions =>
    ``T_sc / T_moon = q_moon / p_sc`` => ``a = (q_moon / p_sc)^(2/3)``.
    Exterior (a > 1) requires ``q_moon > p_sc``.
    """
    return float((q_moon / p_sc) ** (2.0 / 3.0))


def tisserand_cj(a_norm: float, e: float) -> float:
    """Planar Tisserand constant C_J w.r.t. the moon for (a_norm, e).

    ``C_J = 1/a_norm + 2*sqrt(a_norm*(1-e^2))`` (RS07 eq. 2.5, planar i=0).
    Inverse of `eccentricity_from_tisserand`.
    """
    return 1.0 / a_norm + 2.0 * math.sqrt(a_norm * (1.0 - e * e))


def periapsis_norm(a_norm: float, e: float) -> float:
    """Periapsis radius in moon units."""
    return a_norm * (1.0 - e)


def apoapsis_norm(a_norm: float, e: float) -> float:
    """Apoapsis radius in moon units."""
    return a_norm * (1.0 + e)


# ---------------------------------------------------------------------------
# Composed-map state
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ComposedState:
    """Physical (frame-independent) spacecraft state carried across patches.

    ``a_phys_km`` and ``e`` fully specify the Keplerian orbit's size/shape;
    ``varpi_rad`` is the inertial argument of periapsis (Jupiter frame);
    ``t_s`` is absolute time (for moon-phase / encounter self-consistency).
    """

    a_phys_km: float
    e: float
    varpi_rad: float
    t_s: float

    @property
    def peri_km(self) -> float:
        return self.a_phys_km * (1.0 - self.e)

    @property
    def apo_km(self) -> float:
        return self.a_phys_km * (1.0 + self.e)


@dataclass(frozen=True)
class EncounterRecord:
    """One periapsis-passage step of the composed map, with a self-consistency flag."""

    moon: str
    a_phys_km: float
    e: float
    r_peri_norm: float  # periapsis in the ACTIVE moon's units (encounter geometry)
    encounter_valid: bool  # r_peri_norm in the exterior-encounter shell
    kicked: bool  # a real (self-consistent) kick was applied


class ComposedMoonMap:
    """Alternating Europa/Ganymede exterior-Keplerian-map propagator.

    Parameters
    ----------
    europa, ganymede : MoonMapConfig for each moon.
    c_j_ref : reference Tisserand constant used to size each single-moon kick
        table (RS07 running example uses C_J ~ 3).
    a_ref_norm : reference semimajor axis (moon units) for each kick table.
    shell_hi : upper bound on ``r_peri_norm`` for a periapsis passage to count
        as a real close encounter of the active moon (the exterior-encounter
        shell is ``[1.0, shell_hi]``; below 1.0 the exterior map is invalid).
    """

    def __init__(
        self,
        europa: MoonMapConfig,
        ganymede: MoonMapConfig,
        *,
        c_j_ref: float = 3.0,
        a_ref_norm: float = 1.35,
        shell_hi: float = 1.30,
        n_grid: int = 81,
        n_quad: int = 250,
    ) -> None:
        self.cfg: dict[str, MoonMapConfig] = {
            europa.name: europa,
            ganymede.name: ganymede,
        }
        self.c_j_ref = c_j_ref
        self.shell_hi = shell_hi
        k_ref = -0.5 / a_ref_norm
        self.maps: dict[str, KeplerianMap] = {
            name: KeplerianMap(
                mu=cfg.mu,
                C_J=c_j_ref,
                K_ref=k_ref,
                n_grid=n_grid,
                n_quad=n_quad,
            )
            for name, cfg in self.cfg.items()
        }

    # -- frame conversions --------------------------------------------------

    def to_map_state(self, state: ComposedState, moon: str) -> tuple[float, float, float]:
        """Express a physical ComposedState in ``moon``'s map units.

        Returns ``(omega, K, C_J)`` where omega is the rotating-frame periapsis
        angle at time ``state.t_s`` (omega = varpi - theta_moon(t)).
        """
        cfg = self.cfg[moon]
        a_norm = state.a_phys_km / cfg.sma_km
        k = -0.5 / a_norm
        c_j = tisserand_cj(a_norm, state.e)
        theta_moon = cfg.mean_motion_rad_s * state.t_s
        omega = ((state.varpi_rad - theta_moon + math.pi) % _TWO_PI) - math.pi
        return omega, k, c_j

    def from_map_state(
        self, omega: float, k: float, c_j: float, moon: str, t_s: float
    ) -> ComposedState:
        """Rebuild a physical ComposedState from ``moon``-frame map values."""
        cfg = self.cfg[moon]
        a_norm = -0.5 / k
        a_phys = a_norm * cfg.sma_km
        e = eccentricity_from_tisserand(k, c_j)
        theta_moon = cfg.mean_motion_rad_s * t_s
        varpi = ((omega + theta_moon + math.pi) % _TWO_PI) - math.pi
        return ComposedState(a_phys_km=a_phys, e=e, varpi_rad=varpi, t_s=t_s)

    # -- single periapsis passage under one moon ---------------------------

    def step(self, state: ComposedState, moon: str) -> tuple[ComposedState, EncounterRecord]:
        """Advance one periapsis passage under ``moon``.

        The kick is applied ONLY if the periapsis lies in the exterior-encounter
        shell of the active moon (``1.0 <= r_peri_norm <= shell_hi``); otherwise
        the passage is a pure Keplerian coast (time + phase advance, no kick),
        because the moon is not actually near the spacecraft periapsis.
        """
        omega, k, c_j = self.to_map_state(state, moon)
        a_norm = -0.5 / k
        r_peri_norm = periapsis_norm(a_norm, state.e)
        encounter_valid = 1.0 <= r_peri_norm <= self.shell_hi

        # Physical spacecraft period sets the true time advance.
        p_sc_s = _TWO_PI * math.sqrt(state.a_phys_km**3 / _GM_JUPITER)

        if encounter_valid and not math.isnan(state.e):
            omega_new, k_new = self.maps[moon].step(omega, k, u=0.0)
            kicked = abs(k_new - k) > 0.0
            new = self.from_map_state(omega_new, k_new, c_j, moon, state.t_s + p_sc_s)
        else:
            # Pure coast: energy unchanged, inertial periapsis fixed, clock advances.
            kicked = False
            new = ComposedState(
                a_phys_km=state.a_phys_km,
                e=state.e,
                varpi_rad=state.varpi_rad,
                t_s=state.t_s + p_sc_s,
            )

        rec = EncounterRecord(
            moon=moon,
            a_phys_km=state.a_phys_km,
            e=state.e,
            r_peri_norm=r_peri_norm,
            encounter_valid=encounter_valid,
            kicked=kicked,
        )
        return new, rec

    def run_segment(
        self, state: ComposedState, moon: str, n_steps: int
    ) -> tuple[ComposedState, list[EncounterRecord]]:
        """Propagate ``n_steps`` periapsis passages under one moon."""
        recs: list[EncounterRecord] = []
        cur = state
        for _ in range(n_steps):
            cur, rec = self.step(cur, moon)
            recs.append(rec)
            if math.isnan(cur.e):
                break
        return cur, recs

    def run_itinerary(
        self, state: ComposedState, itinerary: list[tuple[str, int]]
    ) -> tuple[ComposedState, list[EncounterRecord]]:
        """Propagate an alternating itinerary ``[(moon, n_steps), ...]``."""
        recs: list[EncounterRecord] = []
        cur = state
        for moon, n in itinerary:
            cur, seg = self.run_segment(cur, moon, n)
            recs.extend(seg)
        return cur, recs


def closure_residual(a: ComposedState, b: ComposedState) -> dict[str, float]:
    """Physical-element closure residual between two composed states.

    Returns fractional a-residual, absolute e-residual, and the inertial
    periapsis-direction residual [rad] (wrapped to [0, pi]).
    """
    da_rel = abs(a.a_phys_km - b.a_phys_km) / a.a_phys_km
    de = abs(a.e - b.e)
    dvarpi = abs(((a.varpi_rad - b.varpi_rad + math.pi) % _TWO_PI) - math.pi)
    return {"da_rel": da_rel, "de": de, "dvarpi_rad": dvarpi}
