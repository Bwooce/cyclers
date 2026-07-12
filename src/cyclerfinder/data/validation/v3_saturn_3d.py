"""V3 independent-integrator cross-check for the Titan-Iapetus 3D-eccentric family
(#574 Stage B).

Same discipline as :mod:`cyclerfinder.data.validation.v3_3d` (the Uranian V3): re-solve
the same Lambert legs as V2, then propagate the spacecraft NUMERICALLY under REBOUND
IAS15 (an architecturally distinct integrator family from the project's usual scipy
DOP853), and compare the IAS15 terminal position to the analytic target. If V3 and V2
AGREE within a tight floor, V2's bounded-drift signature is a real dynamical property of
the model, not an artifact of the DOP853-based Lambert+propagate stack.

The ONLY difference from :mod:`v3_3d` is the analytic state generator: this module uses
:mod:`cyclerfinder.genome.titan_iapetus_corrector`'s eccentric 3D Kepler states (Titan
in-plane, Iapetus inclined+eccentric) instead of the circular-coplanar ``_moon_state`` --
matching :mod:`v2_saturn_3d`'s own model swap, for the same reason (see that module's
docstring). The actual IAS15 leg-propagation helper
(:func:`cyclerfinder.data.validation.v3_3d._ias15_propagate_planet_frame`) is reused
VERBATIM -- it is bare two-body Kepler propagation about the primary and is already fully
system-agnostic (only needs ``mu_primary``); no Saturn-specific reimplementation is
warranted or provided.

Discipline
----------
* NO catalogue writeback. V4 + V4-strict + literature check still gate.
* Framing (mandatory): a PASS is quasi-cycler-class evidence only (see
  :mod:`v2_saturn_3d`'s module docstring for the full framing statement).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.data.validation.v2_saturn_3d import V2Saturn3DVerdict
from cyclerfinder.data.validation.v3_3d import _ias15_propagate_planet_frame
from cyclerfinder.genome.titan_iapetus_corrector import (
    SEQUENCE,
    TitanIapetusClosureParams,
    cycle_period_days,
    iapetus_state,
    leg_tof_days,
    titan_state,
)
from cyclerfinder.search.discovery_campaign import DAY_S

V3_SATURN_AGREEMENT_FLOOR_KMS: Final[float] = 100.0
"""Default V3-vs-V2 agreement floor, km. Mirrors
:data:`cyclerfinder.data.validation.v3_3d.V3_AGREEMENT_FLOOR_KMS`."""

V3_SATURN_N_CYCLES_MIN: Final[int] = 3


@dataclass(frozen=True)
class V3SaturnCycleVerdict:
    cycle_index: int
    converged_legs: int
    n_legs: int
    rendezvous_drift_kms_v3: float
    rendezvous_drift_kms_v2: float
    agreement_kms: float
    ias15_vs_analytic_kepler_kms: float
    notes: str = ""


@dataclass(frozen=True)
class V3Saturn3DVerdict:
    candidate_id: str
    params: TitanIapetusClosureParams
    n_cycles_propagated: int
    integrator: str
    per_cycle: tuple[V3SaturnCycleVerdict, ...]
    per_cycle_drift_kms_v3: tuple[float, ...]
    per_cycle_drift_kms_v2: tuple[float, ...]
    drift_agreement_kms: float
    v3_v2_agreement_floor_kms: float
    passes_v3: bool
    notes: str = ""


def _cycle_v3(
    params: TitanIapetusClosureParams,
    *,
    tof_days: float,
    mu: float,
    t_cycle_offset_days: float,
    ias15_epsilon: float,
) -> tuple[bool, np.ndarray | None, np.ndarray | None, str, float]:
    """One cycle: analytic (eccentric-3D) Lambert targeting + IAS15 leg propagation."""
    n0, n1 = params.n_rev
    r0, v0 = titan_state(params, t_cycle_offset_days + 0.0)
    r1, v1 = iapetus_state(params, t_cycle_offset_days + tof_days)
    r2, _v2 = titan_state(params, t_cycle_offset_days + 2.0 * tof_days)
    tof_s = tof_days * DAY_S

    integrator_label = ""
    worst_ias15_vs_analytic_kms = 0.0
    sc_r_curr: np.ndarray | None = None
    for r_a, v_a, r_b, n_rev in ((r0, v0, r1, n0), (r1, v1, r2, n1)):
        sols = _lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return False, None, None, integrator_label or "REBOUND IAS15", float("inf")
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        r_f_leg, _, label = _ias15_propagate_planet_frame(
            r_a.copy(), best.v1.copy(), tof_s, mu, epsilon=ias15_epsilon
        )
        if not integrator_label:
            integrator_label = label
        leg_offset_kms = float(np.linalg.norm(r_f_leg - r_b))
        worst_ias15_vs_analytic_kms = max(worst_ias15_vs_analytic_kms, leg_offset_kms)
        sc_r_curr = r_f_leg
    if sc_r_curr is None:
        return False, None, None, integrator_label or "REBOUND IAS15", float("inf")
    return True, sc_r_curr, r2.copy(), integrator_label, worst_ias15_vs_analytic_kms


def run_v3_saturn_3d(
    candidate_id: str,
    params: TitanIapetusClosureParams,
    *,
    mu: float,
    v2_verdict: V2Saturn3DVerdict,
    n_cycles: int = V3_SATURN_N_CYCLES_MIN,
    ias15_epsilon: float = 1e-12,
    agreement_floor_kms: float = V3_SATURN_AGREEMENT_FLOOR_KMS,
    notes: str = "",
) -> V3Saturn3DVerdict:
    """Run V3: re-propagate cycle terminal positions with IAS15 and compare to V2."""
    if n_cycles < V3_SATURN_N_CYCLES_MIN:
        raise ValueError(f"V3-Saturn-3D requires n_cycles >= {V3_SATURN_N_CYCLES_MIN}")
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v2_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v2_verdict has only {len(v2_verdict.per_cycle)} cycles but V3 wants {n_cycles}"
        )

    tof_days = leg_tof_days(params.tof_scale)
    period_days = cycle_period_days(params.tof_scale)
    n_legs = len(SEQUENCE) - 1

    per_cycle: list[V3SaturnCycleVerdict] = []
    n_completed = 0
    integrator_label_used = ""
    cycle_zero_r_v3: np.ndarray | None = None

    for k in range(n_cycles):
        t_offset_days = k * period_days
        converged, r_v3, _r_v2, ilabel, ias15_vs_analytic = _cycle_v3(
            params,
            tof_days=tof_days,
            mu=mu,
            t_cycle_offset_days=t_offset_days,
            ias15_epsilon=ias15_epsilon,
        )
        if not integrator_label_used:
            integrator_label_used = ilabel
        if not converged or r_v3 is None:
            per_cycle.append(
                V3SaturnCycleVerdict(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms_v3=float("inf"),
                    rendezvous_drift_kms_v2=float(v2_verdict.per_cycle[k].rendezvous_drift_kms),
                    agreement_kms=float("inf"),
                    ias15_vs_analytic_kepler_kms=float("inf"),
                    notes="Lambert / IAS15 failed at least one leg in this cycle",
                )
            )
            break
        if k == 0:
            cycle_zero_r_v3 = r_v3.copy()
            drift_v3 = 0.0
        else:
            assert cycle_zero_r_v3 is not None
            drift_v3 = float(np.linalg.norm(r_v3 - cycle_zero_r_v3))
        drift_v2 = float(v2_verdict.per_cycle[k].rendezvous_drift_kms)
        agreement = abs(drift_v3 - drift_v2)
        per_cycle.append(
            V3SaturnCycleVerdict(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms_v3=drift_v3,
                rendezvous_drift_kms_v2=drift_v2,
                agreement_kms=agreement,
                ias15_vs_analytic_kepler_kms=ias15_vs_analytic,
            )
        )
        n_completed += 1

    drift_v3_series = tuple(c.rendezvous_drift_kms_v3 for c in per_cycle)
    drift_v2_series = tuple(c.rendezvous_drift_kms_v2 for c in per_cycle)
    drift_agreement = float("inf") if n_completed == 0 else max(c.agreement_kms for c in per_cycle)

    passes_v3 = bool(
        n_completed >= V3_SATURN_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement)
        and drift_agreement <= agreement_floor_kms
    )

    return V3Saturn3DVerdict(
        candidate_id=candidate_id,
        params=params,
        n_cycles_propagated=int(n_completed),
        integrator=integrator_label_used or "REBOUND IAS15",
        per_cycle=tuple(per_cycle),
        per_cycle_drift_kms_v3=drift_v3_series,
        per_cycle_drift_kms_v2=drift_v2_series,
        drift_agreement_kms=float(drift_agreement),
        v3_v2_agreement_floor_kms=float(agreement_floor_kms),
        passes_v3=passes_v3,
        notes=notes,
    )


__all__ = [
    "V3_SATURN_AGREEMENT_FLOOR_KMS",
    "V3_SATURN_N_CYCLES_MIN",
    "V3Saturn3DVerdict",
    "V3SaturnCycleVerdict",
    "run_v3_saturn_3d",
]
