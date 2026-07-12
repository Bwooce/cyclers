"""V2 long-span bounded-drift gauntlet for the Titan-Iapetus 3D eccentric-Keplerian
closure family (#574 Stage B).

Why a Saturn-specific V2, not the generic :mod:`v2_moontour`
--------------------------------------------------------------
:func:`cyclerfinder.data.validation.v2_moontour.run_v2_moontour` (the Uranian gauntlet's
V2) re-solves Lambert legs using the circular-coplanar
:func:`cyclerfinder.search.discovery_campaign._moon_state` -- that IS the Uranian
candidates' own defining model (the Uranian moons' real eccentricity is small enough,
e<=0.004, that the discovery genome itself never modeled it). For Titan-Iapetus the
whole point of #574 is that the candidates' defining model is NOT circular-coplanar:
Iapetus carries a real ~15.5 deg inclination and BOTH moons carry non-negligible
eccentricity (Titan ~0.0288, Iapetus ~0.028 -- 7-25x the Uranian moons'). Feeding these
candidates through the circular-coplanar V2 driver would silently drop the exact
fidelity axis #574 Stage A was built to test, defeating the point of the gate. This
module is therefore the "same model, re-solve over cycles" V2 gate for the Saturn
Titan-Iapetus family, built on
:mod:`cyclerfinder.genome.titan_iapetus_corrector`'s eccentric 3D Kepler propagator
instead of the circular-coplanar one.

Spec reference
--------------
Spec section 14 V2-ballistic -- ">= 3 continuous laps; BOUNDED drift in the dynamic
rotating frame (tolerant of geometric breathing), evaluated in the row's defining
model". Same discipline as :mod:`v2_moontour`, generalized to the 3D eccentric model:

* Cycle 0: the candidate's own converged closure (Titan t=0 -> Iapetus t=tof -> Titan
  t=2*tof).
* Cycle k (k=1, 2, ...): the SAME Lambert legs re-solved with both moons advanced to
  ``t = k * cycle_period_days + (leg-local time)``, Kepler-propagated (mean motion) from
  the SAME epoch mean anomalies ``m0_titan_deg`` / ``m0_iapetus_deg`` the candidate was
  found at (C1 discipline -- no free per-encounter re-phasing, mirroring
  :mod:`cyclerfinder.genome.titan_iapetus_corrector`'s own contract).
* "Drift" = position offset of cycle k's final Titan encounter vs cycle 0's final Titan
  encounter (km) -- exactly :mod:`v2_moontour`'s ``rendezvous_drift_kms`` definition,
  generalized from a circular-coplanar to an eccentric-3D orbit.

PASS criterion: identical form to :mod:`v2_moontour` -- ``n_cycles >= n_cycles_min``
complete AND every cycle's V_inf-continuity residual <= ``closure_floor_kms`` AND the
max inter-cycle rendezvous drift <= ``drift_floor_kms`` (same 50,000 km / 0.05 km/s
defaults as the Uranian gate, per :data:`V2_SATURN_DRIFT_FLOOR_KMS` /
:data:`V2_SATURN_CLOSURE_FLOOR_KMS`).

Discipline
----------
* NO catalogue writeback. V3 + V4 + V4-strict + literature check still gate.
* The 3-cycle minimum is spec section 14, not test-tunable.
* Framing (mandatory, carried from #571-#574): any PASS here is quasi-cycler-class
  evidence about our own idealized (eccentric-Keplerian, non-real-ephemeris) search
  space -- same standing as #312's own Uranian family -- NOT a ballistic-cycler finding
  and NOT a novelty claim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.genome.titan_iapetus_corrector import (
    SEQUENCE,
    TitanIapetusClosureParams,
    cycle_period_days,
    iapetus_state,
    leg_tof_days,
    titan_state,
)
from cyclerfinder.search.discovery_campaign import DAY_S

V2_SATURN_N_CYCLES_MIN: Final[int] = 3
"""Spec section 14 V2-ballistic minimum: >= 3 continuous laps."""

V2_SATURN_DRIFT_FLOOR_KMS: Final[float] = 50_000.0
"""Same-model inter-cycle rendezvous drift floor, km. Mirrors
:data:`cyclerfinder.data.validation.v2_moontour.V2_MOONTOUR_DRIFT_FLOOR_KMS`."""

V2_SATURN_CLOSURE_FLOOR_KMS: Final[float] = 0.05
"""Per-cycle V_inf-continuity residual floor. Matches the project-wide
:data:`cyclerfinder.genome.titan_iapetus_corrector.GATE_RESIDUAL_KMS`."""


@dataclass(frozen=True)
class V2SaturnCycleVerdict:
    """Per-cycle verdict for the Saturn 3D-eccentric V2 gauntlet."""

    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms: float
    closure_residual_kms: float
    notes: str = ""


@dataclass(frozen=True)
class V2Saturn3DVerdict:
    """Frozen V2 verdict for a Titan-Iapetus 3D-eccentric closure candidate."""

    candidate_id: str
    params: TitanIapetusClosureParams
    n_cycles_requested: int
    n_cycles_completed: int
    per_cycle: tuple[V2SaturnCycleVerdict, ...]
    max_drift_kms: float
    max_closure_residual_kms: float
    drift_floor_kms: float
    closure_floor_kms: float
    n_cycles_min: int
    passes_v2: bool
    notes: str = ""


def _cycle_residual(
    params: TitanIapetusClosureParams,
    *,
    tof_days: float,
    mu: float,
    t_cycle_offset_days: float,
) -> tuple[bool, float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Re-solve both legs of one cycle at the given global time offset.

    Mirrors :func:`v2_moontour._cycle_residual` exactly in spirit, with
    ``_moon_state`` calls replaced by the corrector's eccentric 3D Kepler states.
    """
    n0, n1 = params.n_rev
    r0, v0 = titan_state(params, t_cycle_offset_days + 0.0)
    r1, v1 = iapetus_state(params, t_cycle_offset_days + tof_days)
    r2, v2 = titan_state(params, t_cycle_offset_days + 2.0 * tof_days)

    tof_s = tof_days * DAY_S

    sols0 = _lambert(r0, r1, tof_s, mu=mu, max_revs=max(0, n0))
    cands0 = [s for s in sols0 if s.n_revs == n0]
    sols1 = _lambert(r1, r2, tof_s, mu=mu, max_revs=max(0, n1))
    cands1 = [s for s in sols1 if s.n_revs == n1]
    if not cands0 or not cands1:
        return False, math.inf, (r0, r1, r2)
    best0 = min(cands0, key=lambda s: float(np.linalg.norm(s.v1 - v0)))
    best1 = min(cands1, key=lambda s: float(np.linalg.norm(s.v1 - v1)))
    vinf0_out = float(np.linalg.norm(best0.v1 - v0))
    vinf1_in = float(np.linalg.norm(best0.v2 - v1))
    vinf1_out = float(np.linalg.norm(best1.v1 - v1))
    vinf2_in = float(np.linalg.norm(best1.v2 - v2))

    r_mid = abs(vinf1_in - vinf1_out)
    r_wrap = abs(vinf0_out - vinf2_in)
    residual = max(r_mid, r_wrap)
    return True, residual, (r0, r1, r2)


def run_v2_saturn_3d(
    candidate_id: str,
    params: TitanIapetusClosureParams,
    *,
    mu: float,
    n_cycles: int = V2_SATURN_N_CYCLES_MIN,
    drift_floor_kms: float = V2_SATURN_DRIFT_FLOOR_KMS,
    closure_floor_kms: float = V2_SATURN_CLOSURE_FLOOR_KMS,
    notes: str = "",
) -> V2Saturn3DVerdict:
    """Run V2 for a Titan-Iapetus 3D-eccentric closure: re-solve legs over ``n_cycles``.

    Parameters
    ----------
    candidate_id:
        Identifier carried into the verdict.
    params:
        The candidate's converged closure (from
        :mod:`cyclerfinder.genome.titan_iapetus_corrector`).
    mu:
        Saturn's GM (km^3/s^2) -- pass ``PRIMARIES["Saturn"]``.
    n_cycles:
        Cycles to attempt. Must be >= :data:`V2_SATURN_N_CYCLES_MIN`.
    drift_floor_kms, closure_floor_kms:
        Bars per the module docstring.
    notes:
        Free-form audit note.
    """
    if n_cycles < V2_SATURN_N_CYCLES_MIN:
        raise ValueError(
            f"V2-Saturn-3D requires n_cycles >= {V2_SATURN_N_CYCLES_MIN}; got {n_cycles}"
        )
    if drift_floor_kms <= 0.0:
        raise ValueError(f"drift_floor_kms must be > 0; got {drift_floor_kms}")
    if closure_floor_kms <= 0.0:
        raise ValueError(f"closure_floor_kms must be > 0; got {closure_floor_kms}")

    tof_days = leg_tof_days(params.tof_scale)
    period_days = cycle_period_days(params.tof_scale)
    n_legs = len(SEQUENCE) - 1

    per_cycle: list[V2SaturnCycleVerdict] = []
    cycle_zero_final_pos_km: np.ndarray | None = None
    n_completed = 0
    max_drift_kms = 0.0
    max_closure = 0.0

    for k in range(n_cycles):
        t_offset_days = k * period_days
        converged, residual, states = _cycle_residual(
            params, tof_days=tof_days, mu=mu, t_cycle_offset_days=t_offset_days
        )
        if not converged:
            per_cycle.append(
                V2SaturnCycleVerdict(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms=float("inf"),
                    closure_residual_kms=float("inf"),
                    notes="Lambert failed at least one leg in this cycle",
                )
            )
            break
        final_pos_km = states[-1]
        if k == 0:
            cycle_zero_final_pos_km = final_pos_km.copy()
            drift_kms = 0.0
        else:
            assert cycle_zero_final_pos_km is not None
            drift_kms = float(np.linalg.norm(final_pos_km - cycle_zero_final_pos_km))
            max_drift_kms = max(max_drift_kms, drift_kms)
        max_closure = max(max_closure, residual)
        per_cycle.append(
            V2SaturnCycleVerdict(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms=drift_kms,
                closure_residual_kms=residual,
            )
        )
        n_completed += 1

    passes_v2 = bool(
        n_completed >= V2_SATURN_N_CYCLES_MIN
        and max_drift_kms <= drift_floor_kms
        and max_closure <= closure_floor_kms
    )

    return V2Saturn3DVerdict(
        candidate_id=candidate_id,
        params=params,
        n_cycles_requested=int(n_cycles),
        n_cycles_completed=int(n_completed),
        per_cycle=tuple(per_cycle),
        max_drift_kms=float(max_drift_kms),
        max_closure_residual_kms=float(max_closure),
        drift_floor_kms=float(drift_floor_kms),
        closure_floor_kms=float(closure_floor_kms),
        n_cycles_min=V2_SATURN_N_CYCLES_MIN,
        passes_v2=passes_v2,
        notes=notes,
    )


__all__ = [
    "V2_SATURN_CLOSURE_FLOOR_KMS",
    "V2_SATURN_DRIFT_FLOOR_KMS",
    "V2_SATURN_N_CYCLES_MIN",
    "V2Saturn3DVerdict",
    "V2SaturnCycleVerdict",
    "run_v2_saturn_3d",
]
