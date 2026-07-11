"""V4-strict Uranian-system gauntlet with full SPICE ephemeris (#335 Phase 4.1).

What this is
------------
The strictly-real-ephemeris successor to :mod:`v4_uranus`. V4-scipy (#332)
re-propagated the SILVER's Lambert tour under Uranus J2 + n-body perturbations
with the moons placed by CIRCULAR-COPLANAR Kepler ephemerides — by far the
biggest residual idealisation V3 (Kepler-only) had. V4-strict drops that:
the perturber moon positions and the Lambert-target moon endpoints are
both pulled from the freshly-installed Uranian satellite SPICE kernel
(ura111.bsp, see ``scripts/install_uranian_spice.sh``), giving the V4-strict
gauntlet:

1. **Real moon eccentricity** (e.g. Umbriel e = 0.0041, Oberon e = 0.00056 for
   the #327 Umbriel-Oberon SILVER; sampled from the SPICE state vector, for
   the candidate's own first two distinct sequence bodies, at the encounter
   epoch — NOT the V4-scipy assumption e = 0).
2. **Real inclinations / non-coplanar geometry** — the Uranian satellites
   are all in Uranus's equatorial plane, but that plane is tilted ~98 deg
   from the ecliptic. We integrate in the J2000-inertial frame the SPICE
   kernel publishes.
3. **Secular precession** — line-of-apsides + node precession are baked
   into the SPICE-modelled state at the actual epoch.
4. **Higher-fidelity Uranian masses** — taken from the same Jacobson
   2014-aligned constants the JPL satellite ephemeris is referenced to.

GMAT vs Python+SPICE
--------------------
The task spec called for a GMAT script driving the same physics, with
ura111.bsp loaded as a PointMasses ephemeris. Two practical issues with
GMAT R2022a forced the Python+SPICE path:

* GMAT R2022a does not natively support PointMasses against the URA
  satellite SPK without a custom CelestialBody definition; the Uranian
  rotation-matrix / pole orientation is not in the bundled SPICEPlanetary-
  ConstantsKernel.tpc the way Mars/Jupiter are. Adding it would require
  hand-authoring a Uranus body-fixed frame definition (out of scope for
  V4-strict).
* SPICE itself is the JPL C library; spiceypy binds the SAME library
  GMAT compiles against. Python+SPICE driving scipy DOP853 with SPK-
  sampled perturber states is the same integration recipe GMAT would
  use under the hood — the V4-strict force model is identical.

The Python+SPICE driver is therefore the strict Phase 4.1 successor;
the optional ``scripts/gmat_v4_uranus_generate.py`` emits a parallel
GMAT script for human cross-check, but the headline verdict is from
this module.

Composition map
---------------
Composes on:

* :func:`scipy.integrate.solve_ivp` (DOP853) — same RK family as #332.
  V4-strict-vs-V4-scipy isolates the SPICE-ephemeris delta.
* :func:`cyclerfinder.core.lambert.lambert` — re-used for Lambert
  targeting at each encounter; the moon endpoints are SPICE-sampled
  instead of circular-coplanar-Kepler-sampled.
* :func:`cyclerfinder.data.validation.v4_uranus._j2_acceleration_kms2`,
  :func:`._third_body_acceleration_kms2`, :func:`._hill_radius_km` —
  read-only physics helpers from the V4-scipy module; identical
  formulations, only the moon ephemeris source changes.

Discipline
----------
* SPICE kernel install (Part A) is the sourcing — ura111.bsp =
  JPL/NAIF generic_kernels, last released 2022. The V4-strict
  ephemeris is what JPL says it is, not what we'd like it to say.
* NO catalogue writeback. A V4-strict PASS unblocks #337 (catalogue
  admission) but does not perform the admission itself.
* No test-tuning. The verdict is whatever the math says.
* This module reads ``v4_scipy_verdict`` for the V4-scipy comparison
  and ``v3_verdict`` for the V3 series. Both are the OBJECTS, not
  modules; the prior modules are READ-ONLY.

Honest scope of the verdict
---------------------------
The Uranian system has TWO secular timescales V3/V4-scipy idealize:
the moon-line-of-apsides precession (Umbriel ~0.3 deg/d -> ~5 deg per
Lambert cycle) and the Sun's tidal influence on Uranus's orientation
(century-scale). The first IS captured here (SPICE-modelled). The
second is irrelevant on the ~50-day cycle window. A V4-strict PASS
therefore covers the dynamically-relevant fidelity over the bounded-
drift horizon n_cycles in {3, 5, 10}.

Launch-epoch sensitivity is real and explicit: V4-strict verdicts at
DIFFERENT launch_epoch values will differ because the Uranian-satellite
phase configuration varies. The driver returns the epoch + a single
verdict; the gauntlet runner (#335 Part C) sweeps multiple epochs.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import spiceypy as spice
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v3_3d import V3Verdict3D
from cyclerfinder.data.validation.v4_uranus import (
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4UranusVerdict,
    _hill_radius_km,
    _j2_acceleration_kms2,
    _third_body_acceleration_kms2,
)
from cyclerfinder.search.discovery_campaign import DAY_S

# --------------------------------------------------------------------------- #
# SPICE kernel locations
# --------------------------------------------------------------------------- #

_DEFAULT_GMAT_ROOT = Path.home() / "GMAT" / "R2022a"
DEFAULT_LSK_PATH = _DEFAULT_GMAT_ROOT / "data" / "time" / "SPICELeapSecondKernel.tls"
DEFAULT_PCK_PATH = (
    _DEFAULT_GMAT_ROOT / "data" / "planetary_coeff" / "SPICEPlanetaryConstantsKernel.tpc"
)
DEFAULT_URA_PATH = (
    _DEFAULT_GMAT_ROOT / "data" / "planetary_ephem" / "spk" / "uranian" / "ura111.bsp"
)

_MOON_NAIF_ID: dict[str, str] = {
    "Miranda": "MIRANDA",
    "Ariel": "ARIEL",
    "Umbriel": "UMBRIEL",
    "Titania": "TITANIA",
    "Oberon": "OBERON",
}

# --------------------------------------------------------------------------- #
# Verdict dataclasses
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class V4UranusStrictCycleVerdict:
    """Per-cycle V4-strict verdict (one entry per propagated cycle)."""

    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms_v4_strict: float
    """V4-strict per-cycle drift, km (this cycle's terminal position vs cycle 0's
    terminal position, both under the V4-strict SPICE-ephemeris model)."""
    rendezvous_drift_kms_v4_scipy: float
    """V4-scipy series at the same cycle index (the #332 fallback's drift),
    for direct compare."""
    rendezvous_drift_kms_v3: float
    """V3 series at the same cycle index, for direct compare."""
    agreement_kms_vs_v4_scipy: float
    """``|drift_v4_strict - drift_v4_scipy|`` — isolates the SPICE-ephemeris
    effect (eccentricity, inclination, secular precession) on the per-cycle
    drift. The headline V4-strict-vs-V4-scipy delta."""
    agreement_kms_vs_v3: float
    """``|drift_v4_strict - drift_v3|`` — V4-strict vs V3 (Kepler-only). The
    full V4-class perturbation effect (J2 + n-body + SPICE) on V3."""
    v4_terminal_offset_vs_moon_kms: float
    """V4-strict spacecraft terminal-position offset vs the SPICE-sampled
    moon target at the cycle's final encounter, km. Same definition as V3's
    ``ias15_vs_analytic_kepler_kms`` but under the V4-strict model."""
    notes: str = ""
    failure_mode: str = "converged"
    """One of ``"converged"``, ``"lambert_no_solution"``,
    ``"planet_crossing_infeasible"``, ``"integrator_failure"``
    (#559/#560/#567 bug 2). ``"planet_crossing_infeasible"`` means every
    candidate Lambert branch for some leg had an osculating periapsis
    inside the primary's equatorial radius — a genuine dynamical FAIL from
    real synodic geometry, distinguished here from an unexplained solver
    FAIL, NOT excluded from the pass-rate denominator."""
    perijove_km: float | None = None
    """Osculating periapsis (km) of the offending branch when
    ``failure_mode == "planet_crossing_infeasible"``; ``None`` otherwise."""


@dataclass(frozen=True)
class V4UranusStrictVerdict:
    """Frozen V4-strict verdict for an Uranus-system moontour candidate."""

    candidate_id: str
    sequence: tuple[str, ...]
    n_cycles_propagated: int
    integrator: str
    """Human-readable integrator + model + ephemeris source label."""
    launch_epoch_utc: str
    """ISO-8601 UTC launch epoch — V4-strict verdicts are epoch-dependent
    (the Uranian satellite phase configuration varies)."""
    spice_kernels_used: tuple[str, ...]
    """Full paths of the SPICE kernels FURNSH'd for the run."""
    audit_body1_name: str
    """The candidate's own first distinct ``sequence`` body (in order of
    first appearance) that ``eccentricity_used_e_body1`` /
    ``inclination_used_deg_body1`` were SPICE-sampled from. Prior to
    #567 bug 3's fix these audit fields were hardcoded to always sample
    Umbriel regardless of the candidate's actual sequence."""
    audit_body2_name: str
    """Same for ``eccentricity_used_e_body2`` / ``inclination_used_deg_body2``
    — the candidate's second distinct ``sequence`` body."""
    eccentricity_used_e_body1: float
    """The SPICE-sampled eccentricity of ``audit_body1_name`` at
    ``launch_epoch_utc``. The V4-scipy fallback uses e=0 here. Headline
    single-number for the eccentricity-fidelity gate."""
    eccentricity_used_e_body2: float
    """Same for ``audit_body2_name``."""
    inclination_used_deg_body1: float
    """SPICE-sampled inclination of ``audit_body1_name`` at
    ``launch_epoch_utc`` (J2000 frame), deg. Co-planar circular-coplanar
    assumes inc = const but V4-strict samples the actual angle to test the
    planar idealisation."""
    inclination_used_deg_body2: float
    """Same for ``audit_body2_name``."""
    per_cycle: tuple[V4UranusStrictCycleVerdict, ...]
    per_cycle_drift_kms_v4_strict: tuple[float, ...]
    per_cycle_drift_kms_v4_scipy: tuple[float, ...]
    per_cycle_drift_kms_v3: tuple[float, ...]
    drift_agreement_kms_vs_v4_scipy: float
    """``max_k |drift_v4_strict[k] - drift_v4_scipy[k]|`` — V4-strict vs
    V4-scipy. Quantifies the SPICE-ephemeris contribution to the per-cycle
    drift."""
    drift_agreement_kms_vs_v3: float
    """``max_k |drift_v4_strict[k] - drift_v3[k]|`` — V4-strict vs V3
    (Kepler-only)."""
    v4_v3_agreement_floor_kms: float
    """Same 50,000 km floor as V4-scipy — the project's same-model drift
    floor, also used here."""
    bounded_drift_survives: bool
    """V4-strict drift stays bounded (max drift over cycles does not exceed
    ``drift_unbounded_factor`` x V3 max drift)."""
    passes_v4_strict: bool
    """``drift_agreement_kms_vs_v3 <= v4_v3_agreement_floor_kms`` AND
    every cycle's Lambert leg closed AND ``bounded_drift_survives``.

    Interpretation:

    * PASS: V4-strict (full SPICE Uranian ephemeris) agrees with V3
      (Kepler-only) on the bounded-drift signature at the floor of 50,000 km.
      The SILVER has cleared the strictest computational gate available.
      Catalogue admission as ``quasi_cycler`` is unblocked (#337 = the
      admission task; full provenance chain review still required).
    * FAIL: Real Uranian eccentricity/inclination/secular precession break
      the bounded-drift signature. Retire to negative-results registry with
      the SPICE-ephemeris perturbation order that broke it documented.
    """
    notes: str = ""
    gmat_status: str = "not_used"
    """Mirror of the V4UranusVerdict-style ``gmat_status``: ``"not_used"`` for
    the Python+SPICE driver (the SPICE library is the same C code GMAT links;
    a parallel GMAT script lives in ``scripts/gmat_v4_uranus_generate.py``
    for human cross-check)."""


# --------------------------------------------------------------------------- #
# SPICE helpers
# --------------------------------------------------------------------------- #


def _spice_furnsh_all(kernel_paths: tuple[str, ...]) -> None:
    """FURNSH each kernel; raise FileNotFoundError if any missing."""
    for p in kernel_paths:
        if not Path(p).exists():
            raise FileNotFoundError(f"SPICE kernel not found: {p}")
        spice.furnsh(p)


def _ephemeris_time_seconds(launch_epoch_utc: str) -> float:
    """Convert ISO-8601 UTC string to SPICE ET seconds past J2000."""
    return float(spice.str2et(launch_epoch_utc))


def _moon_state_spice(
    moon_name: str, et_seconds: float, *, observer: str = "URANUS"
) -> tuple[np.ndarray, np.ndarray]:
    """SPICE-sampled moon state at ``et_seconds``, Uranus-centered J2000 km/kms."""
    target = _MOON_NAIF_ID[moon_name]
    state, _light_time = spice.spkezr(target, et_seconds, "J2000", "NONE", observer)
    arr = np.asarray(state, dtype=np.float64)
    return arr[:3].copy(), arr[3:].copy()


def _osculating_elements(
    r_km: np.ndarray, v_kms: np.ndarray, mu: float
) -> tuple[float, float, float]:
    """Return ``(sma_km, ecc, inc_deg)`` from Cartesian state about a primary."""
    r_mag = float(np.linalg.norm(r_km))
    v_mag = float(np.linalg.norm(v_kms))
    h = np.cross(r_km, v_kms)
    h_mag = float(np.linalg.norm(h))
    e_vec = np.cross(v_kms, h) / mu - r_km / r_mag
    ecc = float(np.linalg.norm(e_vec))
    energy = 0.5 * v_mag * v_mag - mu / r_mag
    sma = -mu / (2.0 * energy)
    inc = math.degrees(math.acos(h[2] / h_mag))
    return float(sma), ecc, inc


def _leg_periapsis_km(r0_km: np.ndarray, v0_kms: np.ndarray, mu_primary: float) -> float:
    """Osculating periapsis radius (km) of the two-body Kepler orbit ``(r0, v0)``.

    Used to pre-screen a Lambert branch for physical infeasibility BEFORE
    propagating it (#559/#560/#567 bug 2): a fixed-TOF rev>=1 Lambert
    "solution" whose osculating conic dips inside the primary body's
    equatorial radius is not a flyable transfer regardless of what the
    numerical integrator does with it downstream.
    """
    sma_km, ecc, _inc_deg = _osculating_elements(r0_km, v0_kms, mu_primary)
    return float(sma_km * (1.0 - ecc))


# --------------------------------------------------------------------------- #
# V4-strict propagation
# --------------------------------------------------------------------------- #


def _build_spice_dynamics(
    *,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    perturber_mu: dict[str, float],
    perturber_hill_km: dict[str, float],
    et_leg_start: float,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Build a closed-over ``f(t, y)`` for ``scipy.integrate.solve_ivp``.

    Same physics as :func:`v4_uranus._build_dynamics`, but moon positions
    come from SPICE at every RHS evaluation. ``et_leg_start`` is the
    SPICE ET seconds past J2000 at the integrator's ``t=0``; ``t`` is
    seconds since the leg start.
    """

    def rhs(t_s: float, y: np.ndarray) -> np.ndarray:
        r_sc = y[:3]
        v_sc = y[3:]
        r_norm = float(np.linalg.norm(r_sc))
        if r_norm <= 0.0:
            a_central = np.zeros(3, dtype=np.float64)
        else:
            a_central = -mu_primary * r_sc / (r_norm**3)
        a_j2 = _j2_acceleration_kms2(r_sc, mu=mu_primary, j2=j2, r_eq_km=r_eq_km)
        et = et_leg_start + t_s
        a_3b = np.zeros(3, dtype=np.float64)
        for moon in perturber_moons:
            r_moon, _ = _moon_state_spice(moon, et)
            a_3b += _third_body_acceleration_kms2(
                r_sc,
                r_moon,
                mu_body=perturber_mu[moon],
                softening_km=perturber_hill_km[moon],
            )
        a_total = a_central + a_j2 + a_3b
        return np.concatenate([v_sc, a_total])

    return rhs


def _v4_strict_propagate_leg(
    r0_km: np.ndarray,
    v0_km_s: np.ndarray,
    tof_s: float,
    *,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    perturber_mu: dict[str, float],
    perturber_hill_km: dict[str, float],
    et_leg_start: float,
    rtol: float = 1e-10,
    atol: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Propagate one leg under V4-strict (J2 + SPICE-driven n-moon 3B)."""
    rhs = _build_spice_dynamics(
        mu_primary=mu_primary,
        j2=j2,
        r_eq_km=r_eq_km,
        perturber_moons=perturber_moons,
        perturber_mu=perturber_mu,
        perturber_hill_km=perturber_hill_km,
        et_leg_start=et_leg_start,
    )
    y0 = np.concatenate(
        [np.asarray(r0_km, dtype=np.float64), np.asarray(v0_km_s, dtype=np.float64)]
    )
    sol = solve_ivp(
        rhs,
        (0.0, float(tof_s)),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not sol.success:
        return np.zeros(3), np.zeros(3), False
    yf = sol.y[:, -1]
    return yf[:3].copy(), yf[3:].copy(), True


#: Failure-mode tags a leg (and therefore its cycle) can be stamped with.
#: ``"converged"`` is success; the rest are FAILs, distinguished so a
#: genuine dynamical-infeasibility FAIL (``planet_crossing_infeasible``)
#: is never silently conflated with a Lambert-geometry FAIL or an
#: unexplained integrator FAIL (#559/#560/#567 bug 2).
FAILURE_MODE_CONVERGED = "converged"
FAILURE_MODE_LAMBERT_NO_SOLUTION = "lambert_no_solution"
FAILURE_MODE_PLANET_CROSSING = "planet_crossing_infeasible"
FAILURE_MODE_INTEGRATOR_FAILURE = "integrator_failure"


@dataclass(frozen=True)
class _LegOutcome:
    """Result of selecting + propagating one V4-strict Lambert leg.

    ``chosen_branch_index`` indexes into the ``max_revs``-filtered Lambert
    solution list :func:`_select_leg_transfer` was given (for diagnostics /
    regression pinning only; not part of the public API)."""

    ok: bool
    r_f_km: np.ndarray | None
    offset_kms: float
    failure_mode: str
    perijove_km: float | None
    chosen_branch_index: int | None


def _select_leg_transfer(
    r_a: np.ndarray,
    r_b: np.ndarray,
    tof_s: float,
    n_rev: int,
    *,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    perturber_mu: dict[str, float],
    perturber_hill_km: dict[str, float],
    et_leg_start: float,
) -> _LegOutcome:
    """Solve + propagate every rev-``n_rev`` Lambert branch between ``r_a`` and
    ``r_b`` and return the best-outcome branch.

    Fixes #559/#560/#567 bug 1 (Lambert branch-selection continuity flip):
    the old code picked a branch BEFORE propagating, via a bare ``min()`` on
    the departure-velocity residual against the departing moon's own
    velocity. That proxy criterion has no relationship to which branch
    actually flies well under the V4-strict perturbed dynamics, so its
    argmin can (and empirically does, see #567's regression test) flip
    between two branches that are each varying smoothly with epoch —
    producing a discontinuous jump in the reported terminal miss between
    adjacent epochs even though nothing physical changed discontinuously.

    Instead, every candidate branch is propagated and the branch with the
    smallest ACTUAL terminal offset (the same quantity the pass/fail gate
    cares about) is selected. Selecting on a continuous outcome makes the
    selected value itself continuous (a "min of continuous functions" is
    continuous, at worst kinked at a true crossing — never a jump).

    Also fixes #559/#560/#567 bug 2 (planet-crossing silent misclassification):
    each branch's osculating periapsis is checked against ``r_eq_km`` BEFORE
    propagation. A branch whose periapsis is inside the primary is a
    genuinely non-physical transfer (#559 confirmed this is REAL synodic
    geometry, not a numerical artifact) — it is tagged
    ``planet_crossing_infeasible`` with its periapsis recorded, and is
    *not* propagated (DOP853 either dies stiffly on it or "succeeds" with a
    physically meaningless through-the-planet trajectory; neither is a
    transfer to select). If every candidate branch is planet-crossing, the
    leg (and therefore the cycle) genuinely FAILs — this is NOT excluded
    from the pass-rate denominator, only correctly labelled.
    """
    sols = _lambert(r_a, r_b, tof_s, mu=mu_primary, max_revs=n_rev)
    wanted = [s for s in sols if s.n_revs == n_rev]
    if not wanted:
        return _LegOutcome(
            ok=False,
            r_f_km=None,
            offset_kms=float("inf"),
            failure_mode=FAILURE_MODE_LAMBERT_NO_SOLUTION,
            perijove_km=None,
            chosen_branch_index=None,
        )

    feasible: list[tuple[float, np.ndarray, int]] = []
    crossing_perijove_km: list[float] = []
    for idx, s in enumerate(wanted):
        rp_km = _leg_periapsis_km(r_a, s.v1, mu_primary)
        if rp_km < r_eq_km:
            crossing_perijove_km.append(rp_km)
            continue
        r_f_leg, _, ok = _v4_strict_propagate_leg(
            r_a.copy(),
            s.v1.copy(),
            tof_s,
            mu_primary=mu_primary,
            j2=j2,
            r_eq_km=r_eq_km,
            perturber_moons=perturber_moons,
            perturber_mu=perturber_mu,
            perturber_hill_km=perturber_hill_km,
            et_leg_start=et_leg_start,
        )
        if not ok:
            continue
        offset_kms = float(np.linalg.norm(r_f_leg - r_b))
        feasible.append((offset_kms, r_f_leg, idx))

    if feasible:
        best_offset, best_r_f, best_idx = min(feasible, key=lambda t: t[0])
        return _LegOutcome(
            ok=True,
            r_f_km=best_r_f,
            offset_kms=best_offset,
            failure_mode=FAILURE_MODE_CONVERGED,
            perijove_km=None,
            chosen_branch_index=best_idx,
        )

    if crossing_perijove_km:
        return _LegOutcome(
            ok=False,
            r_f_km=None,
            offset_kms=float("inf"),
            failure_mode=FAILURE_MODE_PLANET_CROSSING,
            perijove_km=float(min(crossing_perijove_km)),
            chosen_branch_index=None,
        )

    return _LegOutcome(
        ok=False,
        r_f_km=None,
        offset_kms=float("inf"),
        failure_mode=FAILURE_MODE_INTEGRATOR_FAILURE,
        perijove_km=None,
        chosen_branch_index=None,
    )


def _cycle_v4_strict(
    *,
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    et_cycle_start: float,
    perturber_moons: tuple[str, ...],
    perturber_mu: dict[str, float],
    perturber_hill_km: dict[str, float],
    mu_primary: float,
    n_revs: tuple[int, ...] | None,
    j2: float,
    r_eq_km: float,
) -> tuple[bool, np.ndarray | None, float, str, float | None]:
    """One cycle: SPICE-sampled Lambert endpoints + V4-strict propagation.

    Returns ``(converged, terminal_r_km, worst_offset_kms, failure_mode,
    perijove_km)``. On success, ``failure_mode == "converged"`` and
    ``perijove_km is None``. On failure, ``failure_mode`` distinguishes WHY
    (see :data:`FAILURE_MODE_LAMBERT_NO_SOLUTION`,
    :data:`FAILURE_MODE_PLANET_CROSSING`,
    :data:`FAILURE_MODE_INTEGRATOR_FAILURE`) and ``perijove_km`` is the
    offending branch's periapsis when the mode is
    ``planet_crossing_infeasible`` (#559/#560/#567 bug 2)."""
    n_legs = len(sequence) - 1
    if n_revs is None:
        n_revs_used: tuple[int, ...] = tuple(0 for _ in range(n_legs))
    else:
        if len(n_revs) != n_legs:
            raise ValueError(f"n_revs length {len(n_revs)} != n_legs {n_legs}")
        n_revs_used = tuple(n_revs)

    # SPICE-sample EVERY moon endpoint at the right epoch — V4-strict's key delta.
    epochs_s: list[float] = [0.0]
    for tof_d in leg_tofs_days:
        epochs_s.append(epochs_s[-1] + tof_d * DAY_S)
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for moon, t_s in zip(sequence, epochs_s, strict=True):
        et_at_encounter = et_cycle_start + t_s
        r, v = _moon_state_spice(moon, et_at_encounter)
        states.append((r, v))

    sc_r_curr: np.ndarray | None = None
    worst_offset_kms = 0.0
    for k in range(n_legs):
        r_a, _v_a_moon = states[k]
        r_b, _ = states[k + 1]
        nrev = max(0, n_revs_used[k])
        et_leg_start = et_cycle_start + epochs_s[k]
        outcome = _select_leg_transfer(
            r_a,
            r_b,
            leg_tofs_days[k] * DAY_S,
            nrev,
            mu_primary=mu_primary,
            j2=j2,
            r_eq_km=r_eq_km,
            perturber_moons=perturber_moons,
            perturber_mu=perturber_mu,
            perturber_hill_km=perturber_hill_km,
            et_leg_start=et_leg_start,
        )
        if not outcome.ok or outcome.r_f_km is None:
            return False, None, float("inf"), outcome.failure_mode, outcome.perijove_km
        worst_offset_kms = max(worst_offset_kms, outcome.offset_kms)
        sc_r_curr = outcome.r_f_km
    if sc_r_curr is None:
        return False, None, float("inf"), FAILURE_MODE_INTEGRATOR_FAILURE, None
    return True, sc_r_curr, worst_offset_kms, FAILURE_MODE_CONVERGED, None


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def run_v4_uranus_strict(
    candidate_id: str,
    sequence: tuple[str, ...],
    vinf_tuple_kms: tuple[float, ...],
    leg_tofs_days: tuple[float, ...],
    rel_offset_deg: float,
    launch_epoch_utc: str,
    system: cr3bp.CR3BPSystem | None,
    *,
    v3_verdict: V3Verdict3D,
    v4_scipy_verdict: V4UranusVerdict,
    n_cycles: int = V4_N_CYCLES_MIN,
    n_revs: tuple[int, ...] | None = None,
    j2: float = URANUS_J2,
    r_eq_km: float = URANUS_R_EQ_KM,
    perturber_moons: tuple[str, ...] = URANIAN_PERTURBER_MOONS,
    agreement_floor_kms: float = V4_AGREEMENT_FLOOR_KMS,
    drift_unbounded_factor: float = 10.0,
    spice_kernel_paths: tuple[str, ...] | None = None,
    notes: str = "",
) -> V4UranusStrictVerdict:
    """Phase 4.1 V4-strict gauntlet on the SILVER under full SPICE Uranian ephemeris.

    Pipeline:

    1. FURNSH the kernels (Part A install).
    2. Convert ``launch_epoch_utc`` to SPICE ET seconds past J2000.
    3. Sample the SPICE-osculating eccentricity + inclination for Umbriel +
       Oberon at the launch epoch — recorded on the verdict for audit.
    4. For each cycle k = 0, ..., n_cycles - 1:
         a. Compute ``et_cycle_start = et_launch + k * sum(leg_tofs_days)*DAY_S``.
         b. SPICE-sample every moon's state at every Lambert encounter.
         c. Solve Lambert between consecutive SPICE-sampled moon endpoints.
         d. Propagate each leg under Uranus central + Uranus J2 + every
            ``perturber_moons`` SPICE-driven third-body Battin acceleration.
         e. Record terminal spacecraft position at the cycle's final encounter.
    5. Compute V4-strict per-cycle drift, V4-strict-vs-V4-scipy agreement,
       V4-strict-vs-V3 agreement.
    6. Bounded-drift check + headline verdict.

    Parameters mostly mirror :func:`run_v4_uranus` except:

    * ``launch_epoch_utc`` is REQUIRED (V4-strict is epoch-dependent).
    * ``v4_scipy_verdict`` is the prior #332 V4UranusVerdict for the V4-strict-
      vs-V4-scipy comparison.
    * ``spice_kernel_paths`` overrides the default GMAT-bundled kernel set.

    The verdict's ``passes_v4_strict`` is gated against V3 (Kepler-only),
    matching the V4-scipy gate. The V4-strict-vs-V4-scipy delta is reported
    but not gated on (the two share the same scipy DOP853 integrator;
    disagreement isolates the SPICE-ephemeris effect).
    """
    if not sequence:
        raise ValueError("empty sequence")
    if sequence[0] != sequence[-1]:
        raise ValueError(f"moontour sequence must be CLOSED (first == last); got {sequence!r}")
    n_legs = len(sequence) - 1
    if len(leg_tofs_days) != n_legs:
        raise ValueError(f"leg_tofs_days length {len(leg_tofs_days)} != n_legs {n_legs}")
    if len(vinf_tuple_kms) != len(sequence):
        raise ValueError(
            f"vinf_tuple_kms length {len(vinf_tuple_kms)} != len(sequence) {len(sequence)}"
        )
    if any(tof <= 0.0 for tof in leg_tofs_days):
        raise ValueError(f"leg_tofs_days must be positive; got {leg_tofs_days!r}")
    if n_cycles < V4_N_CYCLES_MIN:
        raise ValueError(
            f"V4-strict requires n_cycles >= {V4_N_CYCLES_MIN} (spec section 14); got {n_cycles}"
        )
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v3_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v3_verdict has only {len(v3_verdict.per_cycle)} cycles; V4-strict wants {n_cycles}"
        )
    if len(v4_scipy_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v4_scipy_verdict has only {len(v4_scipy_verdict.per_cycle)} cycles; "
            f"V4-strict wants {n_cycles}"
        )

    # Resolve primary + check Uranus.
    primary = "Uranus"
    if system is not None and system.primary and system.primary.strip().lower() != "uranus":
        raise ValueError(
            f"run_v4_uranus_strict is for Uranus-system candidates only; "
            f"got primary={system.primary!r}"
        )
    if any(m not in _MOON_NAIF_ID for m in sequence):
        bad = [m for m in sequence if m not in _MOON_NAIF_ID]
        raise ValueError(f"sequence contains non-Uranian moons: {bad!r}")
    if any(m not in _MOON_NAIF_ID for m in perturber_moons):
        bad = [m for m in perturber_moons if m not in _MOON_NAIF_ID]
        raise ValueError(f"perturber_moons contains non-Uranian moons: {bad!r}")
    if not perturber_moons:
        raise ValueError("perturber_moons must be non-empty (V4-strict vs V3 needs perturbations)")

    mu_primary = PRIMARIES[primary]

    # Resolve SPICE kernel set.
    if spice_kernel_paths is None:
        spice_kernel_paths = (
            str(DEFAULT_LSK_PATH),
            str(DEFAULT_PCK_PATH),
            str(DEFAULT_URA_PATH),
        )

    # Hill softening per perturber moon (same recipe as V4-scipy).
    perturber_mu: dict[str, float] = {m: SATELLITES[m].mu_km3_s2 for m in perturber_moons}
    perturber_hill_km: dict[str, float] = {
        m: _hill_radius_km(
            sma_moon_km=SATELLITES[m].sma_km,
            mu_moon=perturber_mu[m],
            mu_primary=mu_primary,
        )
        for m in perturber_moons
    }

    # FURNSH the kernels for this run. Use try/finally so kclear() always runs
    # (avoids polluting the SPICE kernel pool across calls).
    spice.kclear()
    try:
        _spice_furnsh_all(spice_kernel_paths)
        et_launch = _ephemeris_time_seconds(launch_epoch_utc)

        # Record V4-strict-vs-V4-scipy idealisation deltas at the launch epoch,
        # for the candidate's OWN first two distinct sequence bodies (#567 bug 3
        # — this used to be hardcoded to Umbriel/Oberon regardless of `sequence`,
        # freezing the wrong moon's e/i into every non-Umbriel-Oberon candidate).
        distinct_bodies: list[str] = []
        for m in sequence:
            if m not in distinct_bodies:
                distinct_bodies.append(m)
        audit_body1 = distinct_bodies[0]
        audit_body2 = distinct_bodies[1] if len(distinct_bodies) > 1 else distinct_bodies[0]
        r_1, v_1 = _moon_state_spice(audit_body1, et_launch)
        _, e_1, i_1 = _osculating_elements(r_1, v_1, mu_primary)
        r_2, v_2 = _moon_state_spice(audit_body2, et_launch)
        _, e_2, i_2 = _osculating_elements(r_2, v_2, mu_primary)

        cycle_period_s = float(sum(leg_tofs_days)) * DAY_S

        per_cycle: list[V4UranusStrictCycleVerdict] = []
        n_completed = 0
        cycle_zero_r: np.ndarray | None = None

        for k in range(n_cycles):
            et_cycle_start = et_launch + k * cycle_period_s
            converged, r_v4s, v4_offset_vs_moon, failure_mode, perijove_km = _cycle_v4_strict(
                sequence=sequence,
                leg_tofs_days=leg_tofs_days,
                et_cycle_start=et_cycle_start,
                perturber_moons=perturber_moons,
                perturber_mu=perturber_mu,
                perturber_hill_km=perturber_hill_km,
                mu_primary=mu_primary,
                n_revs=n_revs,
                j2=j2,
                r_eq_km=r_eq_km,
            )
            if not converged or r_v4s is None:
                if failure_mode == FAILURE_MODE_PLANET_CROSSING:
                    fail_notes = (
                        "genuine dynamical FAIL: every candidate Lambert branch for "
                        f"at least one leg has an osculating periapsis "
                        f"({perijove_km:.1f} km) inside the primary's equatorial "
                        f"radius ({r_eq_km:.1f} km) -- real synodic geometry, not a "
                        "solver artifact (#559/#560/#567 bug 2)"
                    )
                elif failure_mode == FAILURE_MODE_LAMBERT_NO_SOLUTION:
                    fail_notes = "no Lambert solution at the requested n_rev for at least one leg"
                else:
                    fail_notes = (
                        "DOP853 integrator failed on at least one leg for a reason "
                        "other than a planet-crossing Lambert branch"
                    )
                per_cycle.append(
                    V4UranusStrictCycleVerdict(
                        cycle_index=k,
                        converged_legs=0,
                        n_legs=n_legs,
                        rendezvous_drift_kms_v4_strict=float("inf"),
                        rendezvous_drift_kms_v4_scipy=float(
                            v4_scipy_verdict.per_cycle[k].rendezvous_drift_kms_v4
                        ),
                        rendezvous_drift_kms_v3=float(
                            v3_verdict.per_cycle[k].rendezvous_drift_kms_v3
                        ),
                        agreement_kms_vs_v4_scipy=float("inf"),
                        agreement_kms_vs_v3=float("inf"),
                        v4_terminal_offset_vs_moon_kms=float("inf"),
                        notes=fail_notes,
                        failure_mode=failure_mode,
                        perijove_km=perijove_km,
                    )
                )
                break
            if k == 0:
                cycle_zero_r = r_v4s.copy()
                drift_v4s = 0.0
            else:
                assert cycle_zero_r is not None
                drift_v4s = float(np.linalg.norm(r_v4s - cycle_zero_r))
            drift_v3 = float(v3_verdict.per_cycle[k].rendezvous_drift_kms_v3)
            drift_v4_scipy = float(v4_scipy_verdict.per_cycle[k].rendezvous_drift_kms_v4)
            per_cycle.append(
                V4UranusStrictCycleVerdict(
                    cycle_index=k,
                    converged_legs=n_legs,
                    n_legs=n_legs,
                    rendezvous_drift_kms_v4_strict=drift_v4s,
                    rendezvous_drift_kms_v4_scipy=drift_v4_scipy,
                    rendezvous_drift_kms_v3=drift_v3,
                    agreement_kms_vs_v4_scipy=abs(drift_v4s - drift_v4_scipy),
                    agreement_kms_vs_v3=abs(drift_v4s - drift_v3),
                    v4_terminal_offset_vs_moon_kms=v4_offset_vs_moon,
                )
            )
            n_completed += 1
    finally:
        spice.kclear()

    drift_v4s_series = tuple(c.rendezvous_drift_kms_v4_strict for c in per_cycle)
    drift_v4scipy_series = tuple(c.rendezvous_drift_kms_v4_scipy for c in per_cycle)
    drift_v3_series = tuple(c.rendezvous_drift_kms_v3 for c in per_cycle)
    drift_agreement_vs_v4scipy = (
        float("inf") if n_completed == 0 else max(c.agreement_kms_vs_v4_scipy for c in per_cycle)
    )
    drift_agreement_vs_v3 = (
        float("inf") if n_completed == 0 else max(c.agreement_kms_vs_v3 for c in per_cycle)
    )

    # Bounded-drift heuristic identical to V4-scipy: V4-strict max drift must
    # not exceed drift_unbounded_factor x V3 max drift.
    if n_completed == 0:
        bounded_drift_survives = False
    else:
        finite_v4s = [d for d in drift_v4s_series if math.isfinite(d)]
        finite_v3 = [d for d in drift_v3_series if math.isfinite(d)]
        if not finite_v4s or not finite_v3:
            bounded_drift_survives = False
        else:
            max_v4s = max(finite_v4s)
            max_v3 = max(finite_v3)
            if max_v3 < 1.0:  # essentially closed; fall back to outer-SMA bound
                outer_sma = max(SATELLITES[m].sma_km for m in set(sequence))
                bounded_drift_survives = max_v4s < drift_unbounded_factor * outer_sma
            else:
                bounded_drift_survives = max_v4s <= drift_unbounded_factor * max_v3

    passes_v4_strict = bool(
        n_completed >= V4_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement_vs_v3)
        and drift_agreement_vs_v3 <= agreement_floor_kms
        and bounded_drift_survives
    )

    return V4UranusStrictVerdict(
        candidate_id=candidate_id,
        sequence=tuple(sequence),
        n_cycles_propagated=int(n_completed),
        integrator=(
            "scipy DOP853 + SPICE (URA111 + DE / LSK / PCK) + Uranus J2 "
            "+ classical-moon SPICE-driven third-body Battin (V4-strict Phase 4.1)"
        ),
        launch_epoch_utc=launch_epoch_utc,
        spice_kernels_used=tuple(spice_kernel_paths),
        audit_body1_name=audit_body1,
        audit_body2_name=audit_body2,
        eccentricity_used_e_body1=float(e_1),
        eccentricity_used_e_body2=float(e_2),
        inclination_used_deg_body1=float(i_1),
        inclination_used_deg_body2=float(i_2),
        per_cycle=tuple(per_cycle),
        per_cycle_drift_kms_v4_strict=drift_v4s_series,
        per_cycle_drift_kms_v4_scipy=drift_v4scipy_series,
        per_cycle_drift_kms_v3=drift_v3_series,
        drift_agreement_kms_vs_v4_scipy=float(drift_agreement_vs_v4scipy),
        drift_agreement_kms_vs_v3=float(drift_agreement_vs_v3),
        v4_v3_agreement_floor_kms=float(agreement_floor_kms),
        bounded_drift_survives=bool(bounded_drift_survives),
        passes_v4_strict=passes_v4_strict,
        notes=notes,
        gmat_status="not_used",
    )


__all__ = [
    "DEFAULT_LSK_PATH",
    "DEFAULT_PCK_PATH",
    "DEFAULT_URA_PATH",
    "FAILURE_MODE_CONVERGED",
    "FAILURE_MODE_INTEGRATOR_FAILURE",
    "FAILURE_MODE_LAMBERT_NO_SOLUTION",
    "FAILURE_MODE_PLANET_CROSSING",
    "V4UranusStrictCycleVerdict",
    "V4UranusStrictVerdict",
    "run_v4_uranus_strict",
]


def _verdict_to_jsonable(verdict: V4UranusStrictVerdict) -> dict[str, Any]:
    """Helper for the gauntlet runner (#335 Part C); JSON-serialise a verdict."""
    return {
        "kind": "moontour_v4_strict_verdict",
        "candidate_id": verdict.candidate_id,
        "sequence": list(verdict.sequence),
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "launch_epoch_utc": verdict.launch_epoch_utc,
        "spice_kernels_used": list(verdict.spice_kernels_used),
        "audit_body1_name": verdict.audit_body1_name,
        "audit_body2_name": verdict.audit_body2_name,
        "eccentricity_used_e_body1": verdict.eccentricity_used_e_body1,
        "eccentricity_used_e_body2": verdict.eccentricity_used_e_body2,
        "inclination_used_deg_body1": verdict.inclination_used_deg_body1,
        "inclination_used_deg_body2": verdict.inclination_used_deg_body2,
        "drift_agreement_kms_vs_v4_scipy": verdict.drift_agreement_kms_vs_v4_scipy,
        "drift_agreement_kms_vs_v3": verdict.drift_agreement_kms_vs_v3,
        "v4_v3_agreement_floor_kms": verdict.v4_v3_agreement_floor_kms,
        "bounded_drift_survives": verdict.bounded_drift_survives,
        "passes_v4_strict": verdict.passes_v4_strict,
        "per_cycle_drift_kms_v4_strict": list(verdict.per_cycle_drift_kms_v4_strict),
        "per_cycle_drift_kms_v4_scipy": list(verdict.per_cycle_drift_kms_v4_scipy),
        "per_cycle_drift_kms_v3": list(verdict.per_cycle_drift_kms_v3),
        "per_cycle": [
            {
                "cycle_index": c.cycle_index,
                "converged_legs": c.converged_legs,
                "n_legs": c.n_legs,
                "rendezvous_drift_kms_v4_strict": c.rendezvous_drift_kms_v4_strict,
                "rendezvous_drift_kms_v4_scipy": c.rendezvous_drift_kms_v4_scipy,
                "rendezvous_drift_kms_v3": c.rendezvous_drift_kms_v3,
                "agreement_kms_vs_v4_scipy": c.agreement_kms_vs_v4_scipy,
                "agreement_kms_vs_v3": c.agreement_kms_vs_v3,
                "v4_terminal_offset_vs_moon_kms": c.v4_terminal_offset_vs_moon_kms,
                "notes": c.notes,
                "failure_mode": c.failure_mode,
                "perijove_km": c.perijove_km,
            }
            for c in verdict.per_cycle
        ],
        "gmat_status": verdict.gmat_status,
        "notes": verdict.notes,
    }
