"""V4 higher-fidelity Saturn-system gauntlet (J2 + n-body scipy fallback), #574 Stage B.

Same role as :mod:`cyclerfinder.data.validation.v4_uranus`: adds the dominant V4-class
perturbations the V2/V3 stack omits -- Saturn's J2 zonal harmonic and the classical
regular Saturnian moons as third-body perturbers -- on top of the SAME scipy DOP853
integrator family V2/V3 already compose on (V4 changes the MODEL, not the integrator;
V3 already isolated the integrator-architecture question).

What's reused verbatim from :mod:`v4_uranus`
----------------------------------------------
The physics helpers are fully body-agnostic (they take ``mu``/``j2``/``r_eq_km`` as
parameters, not hardcoded Uranus values), so they are imported and reused directly, NOT
reimplemented:

* :func:`cyclerfinder.data.validation.v4_uranus._j2_acceleration_kms2` (Vallado Eq. 8-37)
* :func:`cyclerfinder.data.validation.v4_uranus._third_body_acceleration_kms2`
  (Battin's third-body formulation)
* :func:`cyclerfinder.data.validation.v4_uranus._hill_radius_km`

What's Saturn-specific
-----------------------
* :data:`SATURN_J2` / :data:`SATURN_R_EQ_KM` (sourced below).
* :data:`SATURN_PERTURBER_MOONS` -- the classical regular Saturnian moons carried in
  ``core/satellites.py`` (Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Iapetus,
  Hyperion).
* The Lambert-targeting AND the tour-moon (Titan, Iapetus) third-body positions use
  :mod:`cyclerfinder.genome.titan_iapetus_corrector`'s eccentric 3D Kepler states (the
  "same defining model" as V2/V3 -- see :mod:`v2_saturn_3d`'s module docstring for why
  the circular-coplanar ``_moon_state`` is NOT appropriate here). The remaining
  non-tour perturber moons (Mimas, Enceladus, Tethys, Dione, Rhea, Hyperion) use the
  circular-coplanar ``_moon_state`` at a deterministic zero phase -- mirroring
  :mod:`v4_uranus`'s own documented treatment of Uranus's non-tour perturbers (their
  real ephemeris is "part of the fidelity gap this whole V4 fallback documents";
  randomizing their phase would introduce an unmotivated free parameter).

Discipline
----------
* NO catalogue writeback. V4-strict + literature check still gate.
* Framing (mandatory): a PASS is quasi-cycler-class evidence only.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v3_saturn_3d import V3Saturn3DVerdict
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
    cycle_period_days,
    iapetus_state,
    leg_tof_days,
    titan_state,
)
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day, _moon_state

# --------------------------------------------------------------------------- #
# SOURCED Saturn constants
# --------------------------------------------------------------------------- #

SATURN_J2: Final[float] = 16290.573e-6
"""Saturn J2 zonal harmonic coefficient.

Source: Iess, L. et al. (2019), "Measurement and implications of Saturn's gravity
field and ring mass", Science 364(6445), 1052-1056 (Cassini Grand Finale radio
science), Table 1: J2 = (16290.573 +/- 0.028) x 10^-6."""

SATURN_R_EQ_KM: Final[float] = 60268.0
"""Saturn equatorial radius (km). Source: IAU/JPL nominal value (also the Iess et al.
2019 reference radius for the zonal-harmonic expansion)."""

SATURN_PERTURBER_MOONS: Final[tuple[str, ...]] = (
    "Mimas",
    "Enceladus",
    "Tethys",
    "Dione",
    "Rhea",
    "Titan",
    "Iapetus",
    "Hyperion",
)
"""Every classical Saturnian moon carried in ``core/satellites.py`` -- all 8 are used
as V4 third-body perturbers (Titan/Iapetus are also the Lambert-tour moons; the other 6
are non-tour third-body perturbers, same role as Miranda/Ariel/Titania in the Uranian
V4)."""

V4_SATURN_AGREEMENT_FLOOR_KMS: Final[float] = 50_000.0
"""Same-model drift floor, km. Mirrors :data:`v4_uranus.V4_AGREEMENT_FLOOR_KMS`."""

V4_SATURN_N_CYCLES_MIN: Final[int] = 3


@dataclass(frozen=True)
class V4SaturnCycleVerdict:
    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms_v4: float
    rendezvous_drift_kms_v3: float
    agreement_kms: float
    v4_terminal_offset_vs_moon_kms: float
    notes: str = ""


@dataclass(frozen=True)
class V4SaturnVerdict:
    candidate_id: str
    params: TitanIapetusClosureParams
    n_cycles_propagated: int
    integrator: str
    per_cycle: tuple[V4SaturnCycleVerdict, ...]
    per_cycle_drift_kms_v4: tuple[float, ...]
    per_cycle_drift_kms_v3: tuple[float, ...]
    drift_agreement_kms: float
    v4_v3_agreement_floor_kms: float
    bounded_drift_survives: bool
    passes_v4: bool
    notes: str = ""


def _perturber_state(
    moon: str,
    params: TitanIapetusClosureParams,
    t_days: float,
    *,
    mu_primary: float,
    non_tour_consts: dict[str, tuple[float, float]],
) -> np.ndarray:
    """Position of ``moon`` at global time ``t_days`` (Saturn frame, km).

    Titan/Iapetus use the corrector's eccentric 3D model (the tour moons' own defining
    model); every other Saturnian moon uses circular-coplanar ``_moon_state`` at a fixed
    zero phase (see module docstring)."""
    if moon == ANCHOR:
        r, _ = titan_state(params, t_days)
        return r
    if moon == FLYBY:
        r, _ = iapetus_state(params, t_days)
        return r
    sma, n_rad_day = non_tour_consts[moon]
    r, _ = _moon_state(0.0, n_rad_day, t_days, sma, mu_primary)
    return r


def _build_dynamics(
    *,
    params: TitanIapetusClosureParams,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    non_tour_consts: dict[str, tuple[float, float]],
    t_cycle_offset_days: float,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Closed-over ``f(t, y)`` for ``scipy.integrate.solve_ivp`` (see module docstring)."""
    perturber_mu = {m: SATELLITES[m].mu_km3_s2 for m in perturber_moons}
    perturber_hill_km = {
        m: _hill_radius_km(
            sma_moon_km=SATELLITES[m].sma_km, mu_moon=perturber_mu[m], mu_primary=mu_primary
        )
        for m in perturber_moons
    }

    def rhs(t_s: float, y: np.ndarray) -> np.ndarray:
        r_sc = y[:3]
        v_sc = y[3:]
        r_norm = float(np.linalg.norm(r_sc))
        a_central = (
            -mu_primary * r_sc / (r_norm**3) if r_norm > 0.0 else np.zeros(3, dtype=np.float64)
        )
        a_j2 = _j2_acceleration_kms2(r_sc, mu=mu_primary, j2=j2, r_eq_km=r_eq_km)
        t_days = t_cycle_offset_days + t_s / DAY_S
        a_3b = np.zeros(3, dtype=np.float64)
        for moon in perturber_moons:
            r_moon = _perturber_state(
                moon, params, t_days, mu_primary=mu_primary, non_tour_consts=non_tour_consts
            )
            a_3b += _third_body_acceleration_kms2(
                r_sc, r_moon, mu_body=perturber_mu[moon], softening_km=perturber_hill_km[moon]
            )
        return np.concatenate([v_sc, a_central + a_j2 + a_3b])

    return rhs


def _v4_propagate_leg(
    r0_km: np.ndarray,
    v0_km_s: np.ndarray,
    tof_s: float,
    *,
    params: TitanIapetusClosureParams,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    non_tour_consts: dict[str, tuple[float, float]],
    t_cycle_offset_days: float,
    rtol: float = 1e-10,
    atol: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray, bool]:
    rhs = _build_dynamics(
        params=params,
        mu_primary=mu_primary,
        j2=j2,
        r_eq_km=r_eq_km,
        perturber_moons=perturber_moons,
        non_tour_consts=non_tour_consts,
        t_cycle_offset_days=t_cycle_offset_days,
    )
    y0 = np.concatenate(
        [np.asarray(r0_km, dtype=np.float64), np.asarray(v0_km_s, dtype=np.float64)]
    )
    sol = solve_ivp(rhs, (0.0, float(tof_s)), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        return np.zeros(3), np.zeros(3), False
    yf = sol.y[:, -1]
    return yf[:3].copy(), yf[3:].copy(), True


def _cycle_v4(
    params: TitanIapetusClosureParams,
    *,
    tof_days: float,
    mu_primary: float,
    perturber_moons: tuple[str, ...],
    non_tour_consts: dict[str, tuple[float, float]],
    j2: float,
    r_eq_km: float,
    t_cycle_offset_days: float,
) -> tuple[bool, np.ndarray | None, float]:
    n0, n1 = params.n_rev
    r0, v0 = titan_state(params, t_cycle_offset_days + 0.0)
    r1, v1 = iapetus_state(params, t_cycle_offset_days + tof_days)
    r2, _v2 = titan_state(params, t_cycle_offset_days + 2.0 * tof_days)
    tof_s = tof_days * DAY_S

    sc_r_curr: np.ndarray | None = None
    worst_offset_kms = 0.0
    for leg_start_days, r_a, v_a, r_b, n_rev in (
        (t_cycle_offset_days + 0.0, r0, v0, r1, n0),
        (t_cycle_offset_days + tof_days, r1, v1, r2, n1),
    ):
        sols = _lambert(r_a, r_b, tof_s, mu=mu_primary, max_revs=max(0, n_rev))
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return False, None, float("inf")
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        r_f_leg, _, ok = _v4_propagate_leg(
            r_a.copy(),
            best.v1.copy(),
            tof_s,
            params=params,
            mu_primary=mu_primary,
            j2=j2,
            r_eq_km=r_eq_km,
            perturber_moons=perturber_moons,
            non_tour_consts=non_tour_consts,
            t_cycle_offset_days=leg_start_days,
        )
        if not ok:
            return False, None, float("inf")
        worst_offset_kms = max(worst_offset_kms, float(np.linalg.norm(r_f_leg - r_b)))
        sc_r_curr = r_f_leg
    if sc_r_curr is None:
        return False, None, float("inf")
    return True, sc_r_curr, worst_offset_kms


def run_v4_saturn(
    candidate_id: str,
    params: TitanIapetusClosureParams,
    *,
    mu_primary: float,
    v3_verdict: V3Saturn3DVerdict,
    n_cycles: int = V4_SATURN_N_CYCLES_MIN,
    j2: float = SATURN_J2,
    r_eq_km: float = SATURN_R_EQ_KM,
    perturber_moons: tuple[str, ...] = SATURN_PERTURBER_MOONS,
    agreement_floor_kms: float = V4_SATURN_AGREEMENT_FLOOR_KMS,
    drift_unbounded_factor: float = 10.0,
    notes: str = "",
) -> V4SaturnVerdict:
    """Run V4 for the Titan-Iapetus 3D-eccentric family: J2 + n-body scipy fallback."""
    if n_cycles < V4_SATURN_N_CYCLES_MIN:
        raise ValueError(f"V4-Saturn requires n_cycles >= {V4_SATURN_N_CYCLES_MIN}")
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v3_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v3_verdict has only {len(v3_verdict.per_cycle)} cycles but V4 wants {n_cycles}"
        )
    if not perturber_moons:
        raise ValueError("perturber_moons must be non-empty")

    non_tour_moons = tuple(m for m in perturber_moons if m not in (ANCHOR, FLYBY))
    non_tour_consts = {
        m: (SATELLITES[m].sma_km, _mean_motion_rad_day(mu_primary, SATELLITES[m].sma_km))
        for m in non_tour_moons
    }

    tof_days = leg_tof_days(params.tof_scale)
    period_days = cycle_period_days(params.tof_scale)
    n_legs = len(SEQUENCE) - 1

    per_cycle: list[V4SaturnCycleVerdict] = []
    n_completed = 0
    cycle_zero_r_v4: np.ndarray | None = None

    for k in range(n_cycles):
        t_offset_days = k * period_days
        converged, r_v4, v4_offset_vs_moon = _cycle_v4(
            params,
            tof_days=tof_days,
            mu_primary=mu_primary,
            perturber_moons=perturber_moons,
            non_tour_consts=non_tour_consts,
            j2=j2,
            r_eq_km=r_eq_km,
            t_cycle_offset_days=t_offset_days,
        )
        if not converged or r_v4 is None:
            per_cycle.append(
                V4SaturnCycleVerdict(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms_v4=float("inf"),
                    rendezvous_drift_kms_v3=float(v3_verdict.per_cycle[k].rendezvous_drift_kms_v3),
                    agreement_kms=float("inf"),
                    v4_terminal_offset_vs_moon_kms=float("inf"),
                    notes="Lambert / DOP853 failed at least one leg in this cycle",
                )
            )
            break
        if k == 0:
            cycle_zero_r_v4 = r_v4.copy()
            drift_v4 = 0.0
        else:
            assert cycle_zero_r_v4 is not None
            drift_v4 = float(np.linalg.norm(r_v4 - cycle_zero_r_v4))
        drift_v3 = float(v3_verdict.per_cycle[k].rendezvous_drift_kms_v3)
        per_cycle.append(
            V4SaturnCycleVerdict(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms_v4=drift_v4,
                rendezvous_drift_kms_v3=drift_v3,
                agreement_kms=abs(drift_v4 - drift_v3),
                v4_terminal_offset_vs_moon_kms=v4_offset_vs_moon,
            )
        )
        n_completed += 1

    drift_v4_series = tuple(c.rendezvous_drift_kms_v4 for c in per_cycle)
    drift_v3_series = tuple(c.rendezvous_drift_kms_v3 for c in per_cycle)
    drift_agreement = float("inf") if n_completed == 0 else max(c.agreement_kms for c in per_cycle)

    if n_completed == 0:
        bounded_drift_survives = False
    else:
        finite_v4 = [d for d in drift_v4_series if math.isfinite(d)]
        finite_v3 = [d for d in drift_v3_series if math.isfinite(d)]
        if not finite_v4 or not finite_v3:
            bounded_drift_survives = False
        else:
            max_v4 = max(finite_v4)
            max_v3 = max(finite_v3)
            if max_v3 < 1.0:
                outer_sma = max(SATELLITES[m].sma_km for m in (ANCHOR, FLYBY))
                bounded_drift_survives = max_v4 < drift_unbounded_factor * outer_sma
            else:
                bounded_drift_survives = max_v4 <= drift_unbounded_factor * max_v3

    passes_v4 = bool(
        n_completed >= V4_SATURN_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement)
        and drift_agreement <= agreement_floor_kms
        and bounded_drift_survives
    )

    return V4SaturnVerdict(
        candidate_id=candidate_id,
        params=params,
        n_cycles_propagated=int(n_completed),
        integrator="scipy DOP853 (Saturn J2 + 8-moon third-body fallback; no SPICE)",
        per_cycle=tuple(per_cycle),
        per_cycle_drift_kms_v4=drift_v4_series,
        per_cycle_drift_kms_v3=drift_v3_series,
        drift_agreement_kms=float(drift_agreement),
        v4_v3_agreement_floor_kms=float(agreement_floor_kms),
        bounded_drift_survives=bool(bounded_drift_survives),
        passes_v4=passes_v4,
        notes=notes,
    )


__all__ = [
    "SATURN_J2",
    "SATURN_PERTURBER_MOONS",
    "SATURN_R_EQ_KM",
    "V4_SATURN_AGREEMENT_FLOOR_KMS",
    "V4_SATURN_N_CYCLES_MIN",
    "V4SaturnCycleVerdict",
    "V4SaturnVerdict",
    "run_v4_saturn",
]

# PRIMARIES re-exported for callers that want the Saturn GM without importing
# core.satellites directly.
SATURN_MU = PRIMARIES["Saturn"]
