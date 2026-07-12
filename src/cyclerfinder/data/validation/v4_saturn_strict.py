"""V4-strict Saturn-system gauntlet with full SPICE ephemeris (#574 Stage B).

The strictly-real-ephemeris successor to :mod:`v4_saturn`, mirroring
:mod:`cyclerfinder.data.validation.v4_uranus_strict`'s role exactly: V4-scipy
re-propagated the candidate's Lambert tour under Saturn J2 + n-body perturbations with
the TOUR moons (Titan, Iapetus) placed by the idealized eccentric-3D Kepler model of
:mod:`cyclerfinder.genome.titan_iapetus_corrector` (fixed inclination/eccentricity/
periapsis-argument, no secular precession). V4-strict drops that idealization entirely:
EVERY moon's position (tour + all 6 non-tour perturbers) is pulled from the real SAT441
Saturnian-satellite SPICE ephemeris (see ``docs/OUTSTANDING.md`` #574 Stage-B kernel-fetch
note; NAIF IDs 606=Titan, 608=Iapetus, confirmed present with coverage 1749-2250),
giving:

1. **Real Titan/Iapetus eccentricity + inclination AT THE ACTUAL EPOCH** -- not the fixed
   ``ECC_TITAN``/``ECC_IAPETUS``/``INCLINATION_DEG`` constants the idealized corrector
   uses.
2. **Secular apsidal/node precession** -- baked into the SPICE-modelled state (the
   corrector fixes argument of periapsis at 0 for both moons; real Titan/Iapetus
   periapses precess).
3. **All 8 Saturnian moons' real (non-circular, non-coplanar) states** as third-body
   perturbers, not the mixed idealized-tour/circular-coplanar-non-tour model V4-scipy
   uses.

#567 fixes inherited from the start (mandatory per #574 Stage B scope)
--------------------------------------------------------------------------
This module is written FRESH (not copy-pasted from a pre-#567 state) with BOTH #567
robustness fixes built in from the outset, exactly mirroring
:mod:`v4_uranus_strict`'s post-#567 ``_select_leg_transfer``:

* **Continuous Lambert branch selection** (#559/#560/#567 bug 1): every candidate
  Lambert branch at the requested ``n_rev`` is PROPAGATED and the branch with the
  smallest ACTUAL terminal offset is selected -- never a pre-propagation
  departure-velocity-match proxy (which #567 proved produces discontinuous branch
  flips between adjacent epochs even though nothing physical changes discontinuously).
* **Tagged, not excluded, planet-crossing infeasibility** (#559/#560/#567 bug 2): a
  branch whose osculating periapsis dips inside Saturn's equatorial radius
  (:data:`cyclerfinder.data.validation.v4_saturn.SATURN_R_EQ_KM`) is a genuinely
  non-physical transfer -- tagged ``FAILURE_MODE_PLANET_CROSSING`` with its periapsis
  recorded and NOT propagated, but still counted as a real FAIL in the pass-rate
  denominator (never silently excluded).

Composition map
---------------
Reuses, VERBATIM, read-only:

* :func:`cyclerfinder.data.validation.v4_saturn._j2_acceleration_kms2`-equivalent
  (imported from :mod:`v4_uranus`, itself body-agnostic) via :mod:`v4_saturn`'s own
  re-exports (``_hill_radius_km``, ``_j2_acceleration_kms2``, ``_third_body_acceleration_kms2``).
* :func:`cyclerfinder.core.lambert.lambert` for Lambert targeting.
* :mod:`spiceypy` for SAT441 SPK reads.

Discipline
----------
* NO catalogue writeback. A V4-strict PASS does not perform admission.
* Framing (mandatory): a PASS is quasi-cycler-class evidence only, same standing as
  #312's own Uranian family -- NOT a ballistic-cycler finding and NOT a novelty claim.
* Launch-epoch sensitivity is real and explicit -- this driver returns a SINGLE verdict
  at a SINGLE epoch; a gauntlet runner sweeps multiple epochs (out of this dispatch's
  scope per #574 Stage B item 3's "single representative epoch" instruction).
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

from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v3_saturn_3d import V3Saturn3DVerdict
from cyclerfinder.data.validation.v4_saturn import (
    SATURN_J2,
    SATURN_PERTURBER_MOONS,
    SATURN_R_EQ_KM,
    V4_SATURN_AGREEMENT_FLOOR_KMS,
    V4_SATURN_N_CYCLES_MIN,
    V4SaturnVerdict,
)
from cyclerfinder.data.validation.v4_uranus import (
    _hill_radius_km,
    _j2_acceleration_kms2,
    _third_body_acceleration_kms2,
)
from cyclerfinder.genome.titan_iapetus_corrector import (
    ANCHOR,
    FLYBY,
    SEQUENCE,
    TitanIapetusClosureParams,
    leg_tof_days,
)
from cyclerfinder.search.discovery_campaign import DAY_S
from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel, ensure_sat441_kernel

_MOON_SPICE_NAME: dict[str, str] = {
    "Mimas": "MIMAS",
    "Enceladus": "ENCELADUS",
    "Tethys": "TETHYS",
    "Dione": "DIONE",
    "Rhea": "RHEA",
    "Titan": "TITAN",
    "Hyperion": "HYPERION",
    "Iapetus": "IAPETUS",
}

# --------------------------------------------------------------------------- #
# Verdict dataclasses
# --------------------------------------------------------------------------- #

FAILURE_MODE_CONVERGED = "converged"
FAILURE_MODE_LAMBERT_NO_SOLUTION = "lambert_no_solution"
FAILURE_MODE_PLANET_CROSSING = "planet_crossing_infeasible"
FAILURE_MODE_INTEGRATOR_FAILURE = "integrator_failure"


@dataclass(frozen=True)
class V4SaturnStrictCycleVerdict:
    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms_v4_strict: float
    rendezvous_drift_kms_v4_scipy: float
    rendezvous_drift_kms_v3: float
    agreement_kms_vs_v4_scipy: float
    agreement_kms_vs_v3: float
    v4_terminal_offset_vs_moon_kms: float
    notes: str = ""
    failure_mode: str = FAILURE_MODE_CONVERGED
    perijove_km: float | None = None
    """Osculating periapsis (km) of the offending branch when
    ``failure_mode == FAILURE_MODE_PLANET_CROSSING``; ``None`` otherwise. Named
    "perijove_km" for field-name parity with :mod:`v4_uranus_strict` (the quantity is
    the periapsis distance to SATURN here, not Jupiter -- kept for schema consistency
    across the project's V4-strict drivers)."""


@dataclass(frozen=True)
class V4SaturnStrictVerdict:
    candidate_id: str
    params: TitanIapetusClosureParams
    n_cycles_propagated: int
    integrator: str
    launch_epoch_utc: str
    spice_kernels_used: tuple[str, ...]
    eccentricity_used_e_titan: float
    """SPICE-sampled Titan eccentricity at ``launch_epoch_utc`` (vs the corrector's
    fixed :data:`cyclerfinder.genome.titan_iapetus_corrector.ECC_TITAN`)."""
    eccentricity_used_e_iapetus: float
    inclination_used_deg_iapetus: float
    """SPICE-sampled Iapetus inclination (to Saturn's equatorial/J2000 reference) at
    ``launch_epoch_utc``, deg (vs the corrector's fixed
    :data:`cyclerfinder.genome.titan_iapetus_corrector.INCLINATION_DEG`)."""
    per_cycle: tuple[V4SaturnStrictCycleVerdict, ...]
    per_cycle_drift_kms_v4_strict: tuple[float, ...]
    per_cycle_drift_kms_v4_scipy: tuple[float, ...]
    per_cycle_drift_kms_v3: tuple[float, ...]
    drift_agreement_kms_vs_v4_scipy: float
    drift_agreement_kms_vs_v3: float
    v4_v3_agreement_floor_kms: float
    bounded_drift_survives: bool
    passes_v4_strict: bool
    notes: str = ""


# --------------------------------------------------------------------------- #
# SPICE helpers
# --------------------------------------------------------------------------- #


def _spice_furnsh_all(kernel_paths: tuple[str, ...]) -> None:
    for p in kernel_paths:
        if not Path(p).exists():
            raise FileNotFoundError(f"SPICE kernel not found: {p}")
        spice.furnsh(p)


def _moon_state_spice(moon_name: str, et_seconds: float) -> tuple[np.ndarray, np.ndarray]:
    """SPICE-sampled moon state at ``et_seconds``, Saturn-centered J2000 km/km-s."""
    target = _MOON_SPICE_NAME[moon_name]
    state, _lt = spice.spkezr(target, et_seconds, "J2000", "NONE", "SATURN")
    arr = np.asarray(state, dtype=np.float64)
    return arr[:3].copy(), arr[3:].copy()


def _osculating_elements(
    r_km: np.ndarray, v_kms: np.ndarray, mu: float
) -> tuple[float, float, float]:
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
    sma_km, ecc, _inc = _osculating_elements(r0_km, v0_kms, mu_primary)
    return float(sma_km * (1.0 - ecc))


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
    def rhs(t_s: float, y: np.ndarray) -> np.ndarray:
        r_sc = y[:3]
        v_sc = y[3:]
        r_norm = float(np.linalg.norm(r_sc))
        a_central = (
            -mu_primary * r_sc / (r_norm**3) if r_norm > 0.0 else np.zeros(3, dtype=np.float64)
        )
        a_j2 = _j2_acceleration_kms2(r_sc, mu=mu_primary, j2=j2, r_eq_km=r_eq_km)
        et = et_leg_start + t_s
        a_3b = np.zeros(3, dtype=np.float64)
        for moon in perturber_moons:
            r_moon, _ = _moon_state_spice(moon, et)
            a_3b += _third_body_acceleration_kms2(
                r_sc, r_moon, mu_body=perturber_mu[moon], softening_km=perturber_hill_km[moon]
            )
        return np.concatenate([v_sc, a_central + a_j2 + a_3b])

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
    sol = solve_ivp(rhs, (0.0, float(tof_s)), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        return np.zeros(3), np.zeros(3), False
    yf = sol.y[:, -1]
    return yf[:3].copy(), yf[3:].copy(), True


@dataclass(frozen=True)
class _LegOutcome:
    ok: bool
    r_f_km: np.ndarray | None
    offset_kms: float
    failure_mode: str
    perijove_km: float | None


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
    """Solve + propagate every rev-``n_rev`` Lambert branch; select by actual terminal
    offset (#567 fix 1) and pre-screen planet-crossing branches (#567 fix 2). See the
    module docstring's "#567 fixes inherited from the start" section -- ported directly
    from :func:`cyclerfinder.data.validation.v4_uranus_strict._select_leg_transfer`,
    generalized to Saturn's ``r_eq_km``."""
    sols = _lambert(r_a, r_b, tof_s, mu=mu_primary, max_revs=n_rev)
    wanted = [s for s in sols if s.n_revs == n_rev]
    if not wanted:
        return _LegOutcome(False, None, float("inf"), FAILURE_MODE_LAMBERT_NO_SOLUTION, None)

    feasible: list[tuple[float, np.ndarray]] = []
    crossing_perijove_km: list[float] = []
    for s in wanted:
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
        feasible.append((float(np.linalg.norm(r_f_leg - r_b)), r_f_leg))

    if feasible:
        best_offset, best_r_f = min(feasible, key=lambda t: t[0])
        return _LegOutcome(True, best_r_f, best_offset, FAILURE_MODE_CONVERGED, None)

    if crossing_perijove_km:
        return _LegOutcome(
            False,
            None,
            float("inf"),
            FAILURE_MODE_PLANET_CROSSING,
            float(min(crossing_perijove_km)),
        )

    return _LegOutcome(False, None, float("inf"), FAILURE_MODE_INTEGRATOR_FAILURE, None)


def _cycle_v4_strict(
    *,
    tof_days: float,
    et_cycle_start: float,
    perturber_moons: tuple[str, ...],
    perturber_mu: dict[str, float],
    perturber_hill_km: dict[str, float],
    mu_primary: float,
    n_revs: tuple[int, int],
    j2: float,
    r_eq_km: float,
) -> tuple[bool, np.ndarray | None, float, str, float | None]:
    epochs_s = (0.0, tof_days * DAY_S, 2.0 * tof_days * DAY_S)
    states = [
        _moon_state_spice(moon, et_cycle_start + t_s)
        for moon, t_s in zip(SEQUENCE, epochs_s, strict=True)
    ]

    sc_r_curr: np.ndarray | None = None
    worst_offset_kms = 0.0
    for leg_idx, n_rev in enumerate(n_revs):
        r_a, _v_a = states[leg_idx]
        r_b, _v_b = states[leg_idx + 1]
        et_leg_start = et_cycle_start + epochs_s[leg_idx]
        outcome = _select_leg_transfer(
            r_a,
            r_b,
            tof_days * DAY_S,
            max(0, n_rev),
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


def run_v4_saturn_strict(
    candidate_id: str,
    params: TitanIapetusClosureParams,
    launch_epoch_utc: str,
    *,
    mu_primary: float,
    v3_verdict: V3Saturn3DVerdict,
    v4_scipy_verdict: V4SaturnVerdict,
    n_cycles: int = V4_SATURN_N_CYCLES_MIN,
    j2: float = SATURN_J2,
    r_eq_km: float = SATURN_R_EQ_KM,
    perturber_moons: tuple[str, ...] = SATURN_PERTURBER_MOONS,
    agreement_floor_kms: float = V4_SATURN_AGREEMENT_FLOOR_KMS,
    drift_unbounded_factor: float = 10.0,
    spice_kernel_paths: tuple[str, ...] | None = None,
    notes: str = "",
) -> V4SaturnStrictVerdict:
    """V4-strict gauntlet on a Titan-Iapetus candidate under full SAT441 SPICE ephemeris.

    Mirrors :func:`cyclerfinder.data.validation.v4_uranus_strict.run_v4_uranus_strict`'s
    pipeline exactly (see that function's docstring for the step-by-step description);
    the deltas are Saturn's own J2/R_eq/perturber set and the SAT441 kernel.
    """
    if n_cycles < V4_SATURN_N_CYCLES_MIN:
        raise ValueError(f"V4-strict-Saturn requires n_cycles >= {V4_SATURN_N_CYCLES_MIN}")
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v3_verdict.per_cycle) < n_cycles:
        raise ValueError(f"v3_verdict has only {len(v3_verdict.per_cycle)} cycles")
    if len(v4_scipy_verdict.per_cycle) < n_cycles:
        raise ValueError(f"v4_scipy_verdict has only {len(v4_scipy_verdict.per_cycle)} cycles")
    if not perturber_moons:
        raise ValueError("perturber_moons must be non-empty")

    if spice_kernel_paths is None:
        spice_kernel_paths = (ensure_leapseconds_kernel(), ensure_sat441_kernel())

    perturber_mu: dict[str, float] = {m: SATELLITES[m].mu_km3_s2 for m in perturber_moons}
    perturber_hill_km: dict[str, float] = {
        m: _hill_radius_km(
            sma_moon_km=SATELLITES[m].sma_km, mu_moon=perturber_mu[m], mu_primary=mu_primary
        )
        for m in perturber_moons
    }

    tof_days = leg_tof_days(params.tof_scale)
    n_legs = len(SEQUENCE) - 1

    spice.kclear()
    try:
        _spice_furnsh_all(spice_kernel_paths)
        et_launch = float(spice.str2et(launch_epoch_utc))

        r_titan, v_titan = _moon_state_spice(ANCHOR, et_launch)
        _, e_titan, _i_titan = _osculating_elements(r_titan, v_titan, mu_primary)
        r_iap, v_iap = _moon_state_spice(FLYBY, et_launch)
        _, e_iap, i_iap = _osculating_elements(r_iap, v_iap, mu_primary)

        cycle_period_s = 2.0 * tof_days * DAY_S

        per_cycle: list[V4SaturnStrictCycleVerdict] = []
        n_completed = 0
        cycle_zero_r: np.ndarray | None = None

        for k in range(n_cycles):
            et_cycle_start = et_launch + k * cycle_period_s
            converged, r_v4s, v4_offset_vs_moon, failure_mode, perijove_km = _cycle_v4_strict(
                tof_days=tof_days,
                et_cycle_start=et_cycle_start,
                perturber_moons=perturber_moons,
                perturber_mu=perturber_mu,
                perturber_hill_km=perturber_hill_km,
                mu_primary=mu_primary,
                n_revs=params.n_rev,
                j2=j2,
                r_eq_km=r_eq_km,
            )
            if not converged or r_v4s is None:
                if failure_mode == FAILURE_MODE_PLANET_CROSSING:
                    fail_notes = (
                        "genuine dynamical FAIL: every candidate Lambert branch for at "
                        f"least one leg has an osculating periapsis ({perijove_km:.1f} km) "
                        f"inside Saturn's equatorial radius ({r_eq_km:.1f} km) -- real "
                        "synodic geometry, not a solver artifact (#567-inherited fix)"
                    )
                elif failure_mode == FAILURE_MODE_LAMBERT_NO_SOLUTION:
                    fail_notes = "no Lambert solution at the requested n_rev for at least one leg"
                else:
                    fail_notes = (
                        "DOP853 integrator failed on at least one leg for a reason other "
                        "than a planet-crossing Lambert branch"
                    )
                per_cycle.append(
                    V4SaturnStrictCycleVerdict(
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
                V4SaturnStrictCycleVerdict(
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
            if max_v3 < 1.0:
                outer_sma = max(SATELLITES[m].sma_km for m in (ANCHOR, FLYBY))
                bounded_drift_survives = max_v4s < drift_unbounded_factor * outer_sma
            else:
                bounded_drift_survives = max_v4s <= drift_unbounded_factor * max_v3

    passes_v4_strict = bool(
        n_completed >= V4_SATURN_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement_vs_v3)
        and drift_agreement_vs_v3 <= agreement_floor_kms
        and bounded_drift_survives
    )

    return V4SaturnStrictVerdict(
        candidate_id=candidate_id,
        params=params,
        n_cycles_propagated=int(n_completed),
        integrator=(
            "scipy DOP853 + SPICE (SAT441 + LSK) + Saturn J2 + 8-moon SPICE-driven "
            "third-body Battin (V4-strict, #567 fixes inherited from the start)"
        ),
        launch_epoch_utc=launch_epoch_utc,
        spice_kernels_used=tuple(spice_kernel_paths),
        eccentricity_used_e_titan=float(e_titan),
        eccentricity_used_e_iapetus=float(e_iap),
        inclination_used_deg_iapetus=float(i_iap),
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
    )


def verdict_to_jsonable(verdict: V4SaturnStrictVerdict) -> dict[str, Any]:
    """JSON-serialise a verdict for gauntlet-runner output."""
    return {
        "kind": "titan_iapetus_v4_strict_verdict",
        "candidate_id": verdict.candidate_id,
        "params": {
            "omega_deg": verdict.params.omega_deg,
            "tof_scale": verdict.params.tof_scale,
            "n_rev": list(verdict.params.n_rev),
            "m0_titan_deg": verdict.params.m0_titan_deg,
            "m0_iapetus_deg": verdict.params.m0_iapetus_deg,
        },
        "n_cycles_propagated": verdict.n_cycles_propagated,
        "integrator": verdict.integrator,
        "launch_epoch_utc": verdict.launch_epoch_utc,
        "eccentricity_used_e_titan": verdict.eccentricity_used_e_titan,
        "eccentricity_used_e_iapetus": verdict.eccentricity_used_e_iapetus,
        "inclination_used_deg_iapetus": verdict.inclination_used_deg_iapetus,
        "drift_agreement_kms_vs_v4_scipy": verdict.drift_agreement_kms_vs_v4_scipy,
        "drift_agreement_kms_vs_v3": verdict.drift_agreement_kms_vs_v3,
        "v4_v3_agreement_floor_kms": verdict.v4_v3_agreement_floor_kms,
        "bounded_drift_survives": verdict.bounded_drift_survives,
        "passes_v4_strict": verdict.passes_v4_strict,
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
                "failure_mode": c.failure_mode,
                "perijove_km": c.perijove_km,
                "notes": c.notes,
            }
            for c in verdict.per_cycle
        ],
        "notes": verdict.notes,
    }


__all__ = [
    "FAILURE_MODE_CONVERGED",
    "FAILURE_MODE_INTEGRATOR_FAILURE",
    "FAILURE_MODE_LAMBERT_NO_SOLUTION",
    "FAILURE_MODE_PLANET_CROSSING",
    "V4SaturnStrictCycleVerdict",
    "V4SaturnStrictVerdict",
    "run_v4_saturn_strict",
    "verdict_to_jsonable",
]

# Re-export for callers.
SATURN_MU = PRIMARIES["Saturn"]
