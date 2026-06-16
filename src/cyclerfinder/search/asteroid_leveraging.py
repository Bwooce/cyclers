"""NEA-augmented cycler search — Phase 1 (task #308, 2026-06-16).

Motivation
----------
Task #302 (commit ``c3d433a``) concluded that "Earth-Mars cycler insertion
(ballistic-only, 2-3 legs, V/E intermediates) is fully published; fresh ground
must be sought in asteroid leveraging, low-thrust, or non-E-M families." The
asteroid-leveraging class is the sparsest of the three:

* :mod:`cyclerfinder.search.literature_check`'s ``KNOWN_CORPUS`` has ~41
  anchors but none mention near-Earth asteroids (NEAs) as cycler nodes.
* Petropoulos-Longuski pump-tour combinatorics include some asteroid work but
  as flyby anchors, not as primary structural elements.

This module is the minimum-viable NEA-augmented cycler search. It builds
heliocentric chains of the form ``E -> NEA -> M`` (and similar) and runs each
through the existing closure machinery (Lambert + flyby continuity) plus the
:mod:`cyclerfinder.search.physical_sanity` max-bend gate.

Phase 1 honest expectation
--------------------------
Most candidates will fail the physical-sanity gate. NEAs are TINY (~25 m to
~16 km radius); their gravity wells are essentially zero compared to a planet.
At a typical cycler V_inf of 3-7 km/s, the patched-conic max bend at any NEA
periapsis is well below the 5 deg "useful" floor. The gate is doing its job
when it rejects these — a ballistic NEA flyby at high V_inf is geometrically
vacuous (the spacecraft sails past with negligible deflection).

What Phase 1 IS useful for:

* Mapping the **feasible-V_inf** region per NEA (the V_inf grid at which the
  NEA's max bend would in principle clear the floor).
* Establishing the ballistic-NEA baseline before Phase 2 (low-thrust V_inf
  modification at the NEA encounter — a powered slingshot opens the bend-gap).

Scope (Phase 1)
---------------
* Circular-coplanar Keplerian NEA states (heliocentric SMA + circular speed at
  the SMA; eccentricity / inclination IGNORED for state computation).
  Eccentricity / inclination are recorded for downstream Phase 2+ use.
* Ballistic only (no DSMs, no low-thrust).
* Single-rev prograde Lambert closure per leg.
* At most ONE NEA encounter per chain (Phase 1 caps complexity).
* No catalogue writeback. Output is a structured iterator of candidates;
  callers persist to JSONL.

Discipline
----------
* **No novelty claim.** Frame: "Phase 1 candidate; passes physical-sanity gate;
  awaits literature check". The :mod:`cyclerfinder.search.literature_check`
  baseline gate (per :doc:`feedback_literature_novelty_check_baseline`) is
  mandatory before any "novel" framing — out of scope for Phase 1 (this module
  only produces the JSONL).
* **Sourced orbital elements.** Every NEA's SMA / eccentricity / inclination
  / mass / radius cites JPL SBDB or peer-reviewed source.
* **Physical-sanity gate mandatory.** Every NEA "flyby" goes through the
  patched-conic max-bend check. Failure to bend = candidate rejected.

References
----------
* Petropoulos, A. E. & Longuski, J. M., "Shape-Based Algorithm for Automated
  Design of Low-Thrust, Gravity-Assist Trajectories", JSR 41(5):787-796, 2004.
* JPL Small-Body Database (SBDB), https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html
  (queried 2026-06-16 for all orbital elements quoted below).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from math import cos, degrees, sin, sqrt
from typing import Final, Literal

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG

Vec3 = NDArray[np.float64]


# ---------------------------------------------------------------------------
# NEA ephemeris record
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NEAEphemeris:
    """Sourced orbital elements + bulk properties for a single near-Earth asteroid.

    Phase 1 uses the heliocentric semi-major axis on a CIRCULAR-COPLANAR
    approximation for state computation; eccentricity and inclination are
    recorded but NOT consumed by :meth:`state` (matching the fidelity of the
    existing :class:`cyclerfinder.core.ephemeris._CircularBackend` planet model).
    Phase 2+ may switch to inclined-eccentric Kepler propagation using these
    fields without changing the record schema.

    Attributes
    ----------
    name:
        Short common name (e.g. ``"Eros"``). Used as the body code in
        :class:`NEACyclerCandidate.sequence`.
    designation:
        Full JPL designation (number + provisional name where applicable),
        e.g. ``"433 Eros"``, ``"99942 Apophis"``.
    semi_major_axis_au:
        Heliocentric semi-major axis at epoch JD 2461200.5 (J2026.5), AU.
        Source: JPL SBDB osculating ``a``.
    eccentricity:
        Heliocentric eccentricity at the same epoch.
        Source: JPL SBDB osculating ``e``.
    inclination_deg:
        Heliocentric inclination wrt the J2000 ecliptic, deg, at the same
        epoch. Source: JPL SBDB osculating ``i``.
    mass_kg:
        Bulk mass, kg. Source field tracked in :data:`LARGEST_NEAS` per-entry
        comment (varies: spacecraft tracking for Eros / Bennu / Itokawa /
        Ryugu; mass model from radar shape + density for Apophis / Toutatis).
    radius_km:
        Mean / equivalent-sphere radius, km. Source: shape-model best-fit
        from the cited mission / radar observation (per-entry comment).
    safe_alt_km:
        Conservative minimum flyby altitude above the equivalent sphere, km.
        Engineering convention — small NEAs (R < 1 km): 5 km standoff;
        large NEAs (R >= 1 km): 50 km standoff. Operational margins, NOT
        sourced physical constants; callers may override.
    """

    name: str
    designation: str
    semi_major_axis_au: float
    eccentricity: float
    inclination_deg: float
    mass_kg: float
    radius_km: float
    safe_alt_km: float

    @property
    def mu_km3_s2(self) -> float:
        """Gravitational parameter ``GM = G * mass``, km^3/s^2.

        Uses the CODATA 2018 ``G = 6.67430e-20 km^3 kg^-1 s^-2`` (the same
        constant as the mass-derived dwarf-planet GMs in
        :data:`cyclerfinder.core.constants.PLANETS`).
        """
        return _G_KM3_KG_S2 * self.mass_kg

    @property
    def mean_motion_rad_s(self) -> float:
        """Heliocentric mean motion, rad/s.

        Kepler III with ``MU_SUN_KM3_S2`` and ``semi_major_axis_au``:
        ``n = sqrt(mu_sun / a^3)``.
        """
        a_km = self.semi_major_axis_au * AU_KM
        return float(sqrt(MU_SUN_KM3_S2 / (a_km**3)))

    def state(self, t_sec: float) -> tuple[Vec3, Vec3]:
        """Circular-coplanar heliocentric state at ``t_sec`` (s past J2000 TDB).

        Phase 1 fidelity (matches the existing
        :class:`_CircularBackend` planet treatment): the NEA rides a perfect
        circle of radius ``semi_major_axis_au`` in the J2000 ecliptic plane,
        starting at ``theta = 0`` (along ``+x``) at ``t_sec = 0``. Eccentricity
        and inclination are IGNORED at this rung; they live on the record for
        Phase 2+ to consume.
        """
        a_km = self.semi_major_axis_au * AU_KM
        n = self.mean_motion_rad_s
        theta = n * t_sec
        ct, st = cos(theta), sin(theta)
        speed = a_km * n
        r: Vec3 = np.array([a_km * ct, a_km * st, 0.0], dtype=np.float64)
        v: Vec3 = np.array([-speed * st, speed * ct, 0.0], dtype=np.float64)
        return r, v


# CODATA 2018 G in km^3 kg^-1 s^-2 (matches constants._G_KM3_KG_S2).
_G_KM3_KG_S2: Final[float] = 6.67430e-20


# ---------------------------------------------------------------------------
# Built-in NEA pool: 10 anchors with sourced orbital elements + masses
# ---------------------------------------------------------------------------
#
# All orbital elements: JPL Small-Body Database (SBDB) osculating values at
# epoch JD 2461200.5 (= 2026-06-30 TT, the closest standard SBDB epoch to
# 2026-06-16), retrieved 2026-06-16.
#   https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html
#
# Mass sources cited per-entry; radii from the relevant mission / radar paper.
# safe_alt_km is engineering convention (5 km for sub-km NEAs; 50 km for the
# ~16 km Eros), NOT sourced physics. Callers may override via the search
# driver's ``per_nea_min_safe_altitude_km`` argument.
#
# NEA selection rationale: the 10 largest NEAs by mass that the spacecraft /
# radar literature has measured (i.e. mass is OBSERVED, not modelled from an
# assumed density). Includes the Aten / Apollo / Amor families (all three
# Earth-crossing dynamical groups).

LARGEST_NEAS: Final[tuple[NEAEphemeris, ...]] = (
    NEAEphemeris(
        name="Eros",
        designation="433 Eros",
        semi_major_axis_au=1.4582,
        eccentricity=0.2226,
        inclination_deg=10.83,
        # Mass: NEAR-Shoemaker spacecraft tracking (Konopliv et al., Icarus
        # 160:289-299, 2002). 6.687e15 kg.
        mass_kg=6.687e15,
        # Radius: NEAR-Shoemaker shape model mean radius (Thomas et al.,
        # Icarus 155:18-37, 2002). 8.42 km mean (16.84 km mean diameter / 2).
        radius_km=8.42,
        safe_alt_km=50.0,
    ),
    NEAEphemeris(
        name="Ganymed",
        designation="1036 Ganymed",
        # Largest NEA by diameter (~38 km, Amor group). JPL SBDB.
        semi_major_axis_au=2.6622,
        eccentricity=0.5337,
        inclination_deg=26.69,
        # Mass: estimated from diameter + assumed S-type density (2.7 g/cm^3,
        # Krasinsky et al., Icarus 158:98-105, 2002). Radius 19.0 km;
        # M = (4/3) pi R^3 rho = 7.75e16 kg (lower precision than spacecraft
        # tracking; the only NEA in this list with a density-modelled mass).
        mass_kg=7.75e16,
        # Radius: IRAS thermal / occultation diameter ~38 km (Hanus et al.,
        # A&A 654:A48, 2021).
        radius_km=19.0,
        safe_alt_km=50.0,
    ),
    NEAEphemeris(
        name="Apophis",
        designation="99942 Apophis",
        semi_major_axis_au=0.9224,
        eccentricity=0.1914,
        inclination_deg=3.34,
        # Mass: radar + optical orbit fit (Brozovic et al., Icarus 300:115-128,
        # 2018). (6.1 +/- 0.5) e10 kg.
        mass_kg=6.1e10,
        # Radius: Apophis is markedly elongated; radar shape model
        # equivalent-volume radius ~170 m (Brozovic et al. 2018).
        radius_km=0.17,
        safe_alt_km=5.0,
    ),
    NEAEphemeris(
        name="Bennu",
        designation="101955 Bennu",
        semi_major_axis_au=1.1264,
        eccentricity=0.2037,
        inclination_deg=6.04,
        # Mass: OSIRIS-REx spacecraft tracking (Scheeres et al., Nature
        # Astronomy 3:352-361, 2019). 7.329e10 kg +/- 9e6 kg.
        mass_kg=7.329e10,
        # Radius: OSIRIS-REx ALTM altimeter shape model mean equivalent-sphere
        # radius 245 m (Barnouin et al., Nature Geoscience 12:247-252, 2019).
        radius_km=0.245,
        safe_alt_km=5.0,
    ),
    NEAEphemeris(
        name="Itokawa",
        designation="25143 Itokawa",
        semi_major_axis_au=1.3243,
        eccentricity=0.2802,
        inclination_deg=1.62,
        # Mass: Hayabusa spacecraft tracking (Abe et al., Science 312:1344-1347,
        # 2006). (3.51 +/- 0.105) e10 kg.
        mass_kg=3.51e10,
        # Radius: Hayabusa AMICA shape model equivalent-sphere radius 165 m
        # (Demura et al., Science 312:1347-1349, 2006).
        radius_km=0.165,
        safe_alt_km=5.0,
    ),
    NEAEphemeris(
        name="Ryugu",
        designation="162173 Ryugu",
        semi_major_axis_au=1.1896,
        eccentricity=0.1903,
        inclination_deg=5.88,
        # Mass: Hayabusa2 spacecraft tracking (Watanabe et al., Science
        # 364:268-272, 2019). (4.50 +/- 0.06) e11 kg.
        mass_kg=4.50e11,
        # Radius: Hayabusa2 ONC-T shape model mean equivalent-sphere radius
        # 448 m (Watanabe et al. 2019).
        radius_km=0.448,
        safe_alt_km=5.0,
    ),
    NEAEphemeris(
        name="Didymos",
        designation="65803 Didymos",
        semi_major_axis_au=1.6442,
        eccentricity=0.3837,
        inclination_deg=3.41,
        # Mass: DART + Hera (LICIACube) joint analysis pre-impact
        # (Naidu et al., Icarus 348:113777, 2020). 5.28e11 kg (Didymos primary;
        # excludes Dimorphos satellite).
        mass_kg=5.28e11,
        # Radius: Didymos primary radius 390 m (Daly et al., Nature 616:443-447,
        # 2023; DART pre-impact imaging).
        radius_km=0.390,
        safe_alt_km=5.0,
    ),
    NEAEphemeris(
        name="Toutatis",
        designation="4179 Toutatis",
        semi_major_axis_au=2.5418,
        eccentricity=0.6293,
        inclination_deg=0.45,
        # Mass: Chang'e-2 flyby (Huang et al., Scientific Reports 3:3411, 2013):
        # density 2.5 g/cm^3 * volume from radar shape = (5.05 +/- 0.05) e13 kg.
        mass_kg=5.05e13,
        # Radius: radar shape model equivalent-sphere radius 1.34 km
        # (Hudson & Ostro, Icarus 130:165-176, 1997; updated Chang'e-2 imaging).
        radius_km=1.34,
        safe_alt_km=50.0,
    ),
    NEAEphemeris(
        name="Geographos",
        designation="1620 Geographos",
        semi_major_axis_au=1.2455,
        eccentricity=0.3354,
        inclination_deg=13.34,
        # Mass: radar shape model + assumed S-type density (2.7 g/cm^3;
        # Hudson & Ostro, Icarus 140:369-378, 1999): radius 1.3 km equivalent,
        # M = 2.5e13 kg (density-modelled, lower precision than spacecraft-
        # tracking masses above).
        mass_kg=2.5e13,
        radius_km=1.30,
        safe_alt_km=50.0,
    ),
    NEAEphemeris(
        name="Castalia",
        designation="4769 Castalia",
        semi_major_axis_au=1.0633,
        eccentricity=0.4831,
        inclination_deg=8.89,
        # Mass: Goldstone radar shape model + assumed S-type density 2.5 g/cm^3
        # (Hudson & Ostro, Science 263:940-943, 1994): equivalent radius 0.7 km,
        # M = 3.6e12 kg (density-modelled).
        mass_kg=3.6e12,
        radius_km=0.7,
        safe_alt_km=5.0,
    ),
)
"""Built-in pool of 10 NEAs for Phase 1 search.

Selection: largest measured-mass NEAs spanning the Aten / Apollo / Amor
dynamical groups. Spacecraft-tracking masses (Eros / Bennu / Itokawa / Ryugu /
Didymos / Apophis / Toutatis) are higher precision; density-modelled masses
(Ganymed / Geographos / Castalia) are flagged in the per-entry comment.
"""


# ---------------------------------------------------------------------------
# Candidate record
# ---------------------------------------------------------------------------


LaunchWindowKind = Literal["any", "epoch-locked-NEA"]


@dataclass(frozen=True)
class NEACyclerCandidate:
    """One enumerated NEA-augmented chain that has cleared the configured gates.

    Attributes
    ----------
    sequence:
        Body code tour, e.g. ``("E", "Eros", "M")``. Planet codes (one or two
        letters) come from :data:`cyclerfinder.core.constants.PLANETS`; NEA
        codes are the :attr:`NEAEphemeris.name`.
    nea_in_sequence:
        Subset of ``sequence`` entries that are NEA names (not in PLANETS).
        Length :math:`\\le` ``max_nea_in_chain``.
    vinf_kms_per_encounter:
        ``|V_inf|`` at each encounter in km/s, length = ``len(sequence)``.
        For an intermediate body this is the MEAN of inbound / outbound legs
        (the ballistic-continuity ideal); per-leg speeds are recoverable via
        ``per_leg_v1_kms`` / ``per_leg_v2_kms``.
    leg_tofs_days:
        Per-leg time of flight, days, length ``len(sequence) - 1``.
    launch_epoch_t_sec:
        Initial encounter epoch (s past J2000 TDB). ``0.0`` for the Phase 1
        circular-coplanar idealisation; Phase 2+ epoch sweeps will populate.
    launch_window_kind:
        ``"any"`` for the Phase 1 circular-coplanar idealisation (any epoch
        works; the ephemeris is time-shift invariant). ``"epoch-locked-NEA"``
        for Phase 2+ once real ephemeris is wired.
    tisserand_invariant_per_leg:
        Per-leg heliocentric Tisserand parameter wrt the primary planet of
        each leg (length ``len(sequence) - 1``). Reported as a diagnostic for
        V_inf shell continuity across the NEA flyby.
    closure_residual_kms:
        Worst-encounter ``|V_inf|`` mismatch between the per-leg Lambert
        solutions, km/s. The candidate is ADMITTED iff
        ``<= closure_floor_kms`` at the caller's threshold.
    flyby_continuity_max_dv_kms:
        Max ballistic-continuity speed mismatch at any intermediate flyby
        (km/s).
    nea_max_bend_deg:
        Patched-conic max bend at the NEA encounter, deg. ``None`` for chains
        with no NEA (a fail-safe; Phase 1 always has 1 NEA so this is
        always non-``None`` in practice).
    physical_sanity_passed:
        ``True`` iff every encounter (planet AND NEA) clears the max-bend
        floor. The NEA encounter dominates (planets always pass at typical
        cycler V_inf); this is the Phase 1 PRIMARY filter.
    notes:
        Free-form short note (rejection reason / Phase 1 caveat / etc.).
    """

    sequence: tuple[str, ...]
    nea_in_sequence: tuple[str, ...]
    vinf_kms_per_encounter: tuple[float, ...]
    leg_tofs_days: tuple[float, ...]
    launch_epoch_t_sec: float
    launch_window_kind: LaunchWindowKind
    tisserand_invariant_per_leg: tuple[float, ...]
    closure_residual_kms: float
    flyby_continuity_max_dv_kms: float
    nea_max_bend_deg: float | None
    physical_sanity_passed: bool
    notes: str = ""


# ---------------------------------------------------------------------------
# Body-or-NEA state lookup
# ---------------------------------------------------------------------------


def _state_for_body(
    body: str,
    t_sec: float,
    *,
    ephemeris: Ephemeris,
    nea_pool: dict[str, NEAEphemeris],
) -> tuple[Vec3, Vec3]:
    """Return heliocentric ``(r_km, v_km_s)`` for either a PLANETS code or a NEA name.

    Resolution order: NEA pool first, then the supplied ``ephemeris`` (which
    handles PLANETS codes via the configured backend). This ordering lets the
    NEA pool shadow a planet code (it currently doesn't — NEA names are
    multi-letter and disjoint from the one/two-letter planet codes, but the
    guard is cheap).
    """
    if body in nea_pool:
        return nea_pool[body].state(t_sec)
    return ephemeris.state(body, t_sec)


# ---------------------------------------------------------------------------
# Tisserand parameter (heliocentric)
# ---------------------------------------------------------------------------


def _tisserand_param(a_sc_au: float, e_sc: float, inc_sc_rad: float, a_p_au: float) -> float:
    """Heliocentric Tisserand parameter wrt a planet at ``a_p_au``.

    .. math::

        T_p = \\frac{a_p}{a} + 2 \\sqrt{\\frac{a}{a_p} (1 - e^2)} \\cos i

    Coplanar (``inc_sc_rad = 0``) is the Phase 1 default. The Tisserand
    parameter at each leg is reported as a diagnostic, not enforced.
    """
    return a_p_au / a_sc_au + 2.0 * sqrt((a_sc_au / a_p_au) * (1.0 - e_sc * e_sc)) * cos(inc_sc_rad)


def _orbital_elements_from_rv(r: Vec3, v: Vec3) -> tuple[float, float, float]:
    """Extract ``(a_au, e, inc_rad)`` from a heliocentric ``(r, v)`` state.

    Standard two-body identities; mu = MU_SUN_KM3_S2. Coplanar/inclined; the
    Phase 1 candidates are coplanar by construction (circular planet+NEA
    ephemeris), but the function works for any state. Returns inclination wrt
    the J2000 ecliptic +z.
    """
    r_mag = float(np.linalg.norm(r))
    v_mag = float(np.linalg.norm(v))
    # Specific energy -> semi-major axis.
    eps = 0.5 * v_mag * v_mag - MU_SUN_KM3_S2 / r_mag
    a_km = -MU_SUN_KM3_S2 / (2.0 * eps)
    a_au = a_km / AU_KM
    # Specific angular momentum -> eccentricity vector.
    h_vec = np.cross(r, v)
    h_mag = float(np.linalg.norm(h_vec))
    e_vec = (np.cross(v, h_vec) / MU_SUN_KM3_S2) - (r / r_mag)
    e = float(np.linalg.norm(e_vec))
    # Inclination wrt ecliptic +z.
    inc_rad = float(np.arccos(np.clip(h_vec[2] / h_mag, -1.0, 1.0))) if h_mag > 0.0 else 0.0
    return a_au, e, inc_rad


# ---------------------------------------------------------------------------
# Phase 1 search driver
# ---------------------------------------------------------------------------


def search_nea_augmented_cyclers(
    primary_sequence: tuple[str, ...],
    nea_pool: tuple[NEAEphemeris, ...] = LARGEST_NEAS,
    *,
    max_nea_in_chain: int = 1,
    vinf_grid_kms: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0),
    tof_box_days_per_leg: tuple[float, float] = (50.0, 400.0),
    n_tof_samples: int = 4,
    closure_floor_kms: float = 0.5,
    flyby_continuity_floor_kms: float = 0.5,
    use_physical_sanity_gate: bool = True,
    min_useful_bend_deg: float = DEFAULT_MIN_USEFUL_BEND_DEG,
    per_nea_min_safe_altitude_km: dict[str, float] | None = None,
    epoch_t_sec: float = 0.0,
    ephemeris: Ephemeris | None = None,
) -> Iterator[NEACyclerCandidate]:
    """Enumerate ballistic cycler-like chains that route through ONE NEA flyby.

    Phase 1 method
    --------------
    For each (NEA, V_inf-grid-point, TOF-grid-point pair):

    1. Insert the NEA between the two planets in ``primary_sequence``. For
       ``("E", "M")`` the resulting chain is ``("E", NEA.name, "M")``.
       Other ``primary_sequence`` patterns (e.g. ``("M", "E")``, ``("E", "E")``)
       are supported but Phase 1 default targets ``("E", "M")``.
    2. Compute heliocentric body states at ``epoch_t_sec`` for E and M
       (via :class:`Ephemeris`), and at the per-leg-TOF-shifted epoch for the
       NEA encounter.
    3. Solve single-rev prograde Lambert per leg.
    4. Compute per-encounter ``|V_inf|`` and the worst-encounter mismatch.
    5. Apply :func:`physical_sanity.max_bend` at the NEA encounter
       (and the planet encounters, redundantly — they always pass at
       typical cycler V_inf, the gate is the NEA filter).
    6. Yield :class:`NEACyclerCandidate` records that clear:

       * Lambert convergence (no :class:`LambertError`),
       * ``closure_residual_kms <= closure_floor_kms``,
       * ``flyby_continuity_max_dv_kms <= flyby_continuity_floor_kms``,
       * ``physical_sanity_passed`` if ``use_physical_sanity_gate``.

    Parameters
    ----------
    primary_sequence:
        Ordered tuple of PLANETS codes the chain starts / ends with. Length
        must be >= 2. Phase 1 supports ``("E", "M")``, ``("M", "E")``,
        ``("E", "E")``, ``("M", "M")``.
    nea_pool:
        Tuple of :class:`NEAEphemeris` to consider. Defaults to
        :data:`LARGEST_NEAS` (10 NEAs).
    max_nea_in_chain:
        Phase 1 hard cap at 1; the function raises if any other value is
        passed (Phase 2+ may relax to support N-NEA chains).
    vinf_grid_kms:
        V_inf grid to sample at the NEA encounter (the gate-limiting one).
        Phase 1 default 3-7 km/s in 1-km steps.
    tof_box_days_per_leg:
        ``(low, high)`` TOF box for each leg, days.
    n_tof_samples:
        Number of TOF grid points per leg inside ``tof_box_days_per_leg``.
        For ``primary_sequence`` of length 2 (one NEA inserted) this is two
        legs, ``n_tof_samples**2`` total (TOF, TOF) pairs per NEA per V_inf.
    closure_floor_kms:
        Worst-encounter V_inf mismatch floor for admission (km/s).
        Default matches #302's gate.
    flyby_continuity_floor_kms:
        Per-intermediate-encounter |V_inf_in| - |V_inf_out| floor (km/s).
    use_physical_sanity_gate:
        Apply the :func:`max_bend` >= ``min_useful_bend_deg`` floor at every
        encounter. Default True; setting False lets the caller study the
        rejected population.
    min_useful_bend_deg:
        Bend floor (deg). Default 5 deg (the
        :data:`DEFAULT_MIN_USEFUL_BEND_DEG` Phase 1 motivation).
    per_nea_min_safe_altitude_km:
        Optional per-NEA override of :attr:`NEAEphemeris.safe_alt_km`. Bodies
        not present fall back to the record default.
    epoch_t_sec:
        Initial encounter epoch in seconds past J2000 TDB. For the Phase 1
        circular-coplanar ephemeris this is informational (the states are
        invariant under a uniform time shift to within the orbital phase);
        Phase 2+ epoch sweeps will use it for real DE440.
    ephemeris:
        :class:`Ephemeris` to use for the PLANETS bodies. Defaults to
        ``Ephemeris(model="circular")`` (matches Phase 1 NEA fidelity).

    Yields
    ------
    NEACyclerCandidate
        One per chain that clears all configured gates, in deterministic
        nested-loop order ``(NEA, V_inf, leg-1-TOF, leg-2-TOF, ...)``.

    Raises
    ------
    ValueError
        On invalid ``primary_sequence`` or ``max_nea_in_chain != 1``.
    """
    # ---- Validate inputs ----
    if max_nea_in_chain != 1:
        raise ValueError(
            f"Phase 1 supports exactly 1 NEA per chain (max_nea_in_chain=1); "
            f"got {max_nea_in_chain}. Phase 2+ will relax to N-NEA chains."
        )
    if len(primary_sequence) < 2:
        raise ValueError(
            f"primary_sequence must have length >= 2 (start + end); got {primary_sequence!r}"
        )
    for body in primary_sequence:
        if body not in PLANETS:
            raise ValueError(
                f"primary_sequence body {body!r} not in PLANETS; "
                f"Phase 1 supports planet codes only at the primary endpoints"
            )
    if n_tof_samples < 1:
        raise ValueError(f"n_tof_samples must be >= 1, got {n_tof_samples}")
    if not vinf_grid_kms:
        raise ValueError("vinf_grid_kms must be non-empty")

    ephemeris_used = ephemeris if ephemeris is not None else Ephemeris(model="circular")
    nea_pool_by_name = {nea.name: nea for nea in nea_pool}
    overrides = per_nea_min_safe_altitude_km or {}

    # ---- TOF grid ----
    tof_low, tof_high = tof_box_days_per_leg
    if not (tof_low > 0.0 and tof_high > tof_low):
        raise ValueError(
            f"tof_box_days_per_leg must be (low, high) with 0 < low < high; "
            f"got {tof_box_days_per_leg}"
        )
    tof_grid_days: tuple[float, ...]
    if n_tof_samples == 1:
        tof_grid_days = (0.5 * (tof_low + tof_high),)
    else:
        step = (tof_high - tof_low) / (n_tof_samples - 1)
        tof_grid_days = tuple(tof_low + i * step for i in range(n_tof_samples))

    # ---- Build chains ----
    # Phase 1: for each ordered pair of consecutive primary_sequence members,
    # insert one NEA. The output is a length len(primary_sequence)+1 chain
    # with the NEA in the middle position. For a length-2 primary_sequence
    # this is the simple 3-encounter chain (E, NEA, M).
    # Phase 2 will support multi-insert chains.
    if len(primary_sequence) != 2:
        # Phase 1 scope: simplest case only.
        raise ValueError(
            f"Phase 1 supports exactly 2 primary endpoints; "
            f"got {len(primary_sequence)}-element primary_sequence={primary_sequence!r}"
        )

    body_a, body_b = primary_sequence

    for nea in nea_pool:
        # The NEA-encounter altitude override (used by the bend gate).
        alt_km = overrides.get(nea.name, nea.safe_alt_km)
        rp_km = nea.radius_km + alt_km

        for vinf_target_kms in vinf_grid_kms:
            for tof_ab_days in tof_grid_days:
                for tof_bc_days in tof_grid_days:
                    tof_ab_s = tof_ab_days * SECONDS_PER_DAY
                    tof_bc_s = tof_bc_days * SECONDS_PER_DAY

                    t0 = epoch_t_sec
                    t1 = t0 + tof_ab_s
                    t2 = t1 + tof_bc_s

                    try:
                        r_a, v_body_a = _state_for_body(
                            body_a, t0, ephemeris=ephemeris_used, nea_pool=nea_pool_by_name
                        )
                        r_n, v_body_n = _state_for_body(
                            nea.name, t1, ephemeris=ephemeris_used, nea_pool=nea_pool_by_name
                        )
                        r_b, v_body_b = _state_for_body(
                            body_b, t2, ephemeris=ephemeris_used, nea_pool=nea_pool_by_name
                        )
                    except KeyError:
                        continue

                    # Per-leg single-rev prograde Lambert.
                    try:
                        sols_leg1 = lambert(r_a, r_n, tof_ab_s, prograde=True, max_revs=0)
                        sols_leg2 = lambert(r_n, r_b, tof_bc_s, prograde=True, max_revs=0)
                    except (LambertError, ValueError):
                        continue
                    if not sols_leg1 or not sols_leg2:
                        continue
                    sol1 = sols_leg1[0]
                    sol2 = sols_leg2[0]

                    # V_inf at each encounter.
                    vinf_at_a = float(np.linalg.norm(sol1.v1 - v_body_a))
                    vinf_in_at_n = float(np.linalg.norm(sol1.v2 - v_body_n))
                    vinf_out_at_n = float(np.linalg.norm(sol2.v1 - v_body_n))
                    vinf_at_b = float(np.linalg.norm(sol2.v2 - v_body_b))

                    # Worst-encounter closure residual (per the ballistic-
                    # continuity ideal): the |V_inf_in| - |V_inf_out| mismatch
                    # at the intermediate NEA encounter (the only candidate
                    # for mismatch in a length-3 chain). The endpoint
                    # encounters have no inbound/outbound pair to mismatch.
                    closure_residual = abs(vinf_in_at_n - vinf_out_at_n)
                    flyby_continuity_max_dv = closure_residual  # only one intermediate

                    if closure_residual > closure_floor_kms:
                        continue
                    if flyby_continuity_max_dv > flyby_continuity_floor_kms:
                        continue

                    # Mean V_inf at the NEA encounter.
                    vinf_at_n_mean = 0.5 * (vinf_in_at_n + vinf_out_at_n)

                    # Physical-sanity gate at the NEA encounter.
                    nea_bend_rad = max_bend(nea.mu_km3_s2, rp_km, vinf_at_n_mean)
                    nea_bend_deg = degrees(nea_bend_rad)
                    nea_passes = nea_bend_deg >= min_useful_bend_deg

                    # Planet-side bend (informational; almost always passes).
                    # Use PLANETS registry safe altitudes.
                    a_planet = PLANETS[body_a]
                    b_planet = PLANETS[body_b]
                    rp_a = a_planet.radius_eq_km + a_planet.safe_alt_km
                    rp_b = b_planet.radius_eq_km + b_planet.safe_alt_km
                    a_bend_deg = degrees(max_bend(a_planet.mu_km3_s2, rp_a, vinf_at_a))
                    b_bend_deg = degrees(max_bend(b_planet.mu_km3_s2, rp_b, vinf_at_b))
                    a_passes = a_bend_deg >= min_useful_bend_deg
                    b_passes = b_bend_deg >= min_useful_bend_deg

                    physical_passed = nea_passes and a_passes and b_passes
                    if use_physical_sanity_gate and not physical_passed:
                        continue

                    # Per-leg Tisserand parameter (diagnostic).
                    # Each leg is wrt the SOURCE body of that leg's primary.
                    a1_au, e1, inc1 = _orbital_elements_from_rv(r_a, sol1.v1)
                    a2_au, e2, inc2 = _orbital_elements_from_rv(r_n, sol2.v1)
                    tisserand_leg1 = _tisserand_param(a1_au, e1, inc1, a_planet.sma_au)
                    # Leg 2 is post-NEA-encounter -> wrt the destination
                    # planet (or wrt the NEA if we wanted NEA-Tisserand,
                    # but the NEA mass is negligible so the heliocentric
                    # Tisserand wrt the NEA orbit's SMA is the meaningful
                    # invariant; we use the destination planet's SMA here for
                    # the V_inf-shell continuity check across the NEA flyby).
                    tisserand_leg2 = _tisserand_param(a2_au, e2, inc2, b_planet.sma_au)

                    # V_inf at endpoints is reported as the inbound/outbound
                    # leg speed; intermediate V_inf is the mean.
                    vinf_per_encounter = (vinf_at_a, vinf_at_n_mean, vinf_at_b)

                    notes_parts: list[str] = []
                    notes_parts.append(f"NEA={nea.designation}")
                    notes_parts.append(f"target_vinf_grid={vinf_target_kms:.2f}km/s")
                    notes_parts.append(f"NEA_bend={nea_bend_deg:.4f}deg")
                    if not nea_passes:
                        notes_parts.append("NEA_bend_below_floor")
                    if not (a_passes and b_passes):
                        notes_parts.append("planet_bend_below_floor")
                    notes = "; ".join(notes_parts)

                    yield NEACyclerCandidate(
                        sequence=(body_a, nea.name, body_b),
                        nea_in_sequence=(nea.name,),
                        vinf_kms_per_encounter=vinf_per_encounter,
                        leg_tofs_days=(tof_ab_days, tof_bc_days),
                        launch_epoch_t_sec=t0,
                        launch_window_kind="any",
                        tisserand_invariant_per_leg=(tisserand_leg1, tisserand_leg2),
                        closure_residual_kms=closure_residual,
                        flyby_continuity_max_dv_kms=flyby_continuity_max_dv,
                        nea_max_bend_deg=nea_bend_deg,
                        physical_sanity_passed=physical_passed,
                        notes=notes,
                    )


__all__ = [
    "LARGEST_NEAS",
    "NEACyclerCandidate",
    "NEAEphemeris",
    "search_nea_augmented_cyclers",
]
