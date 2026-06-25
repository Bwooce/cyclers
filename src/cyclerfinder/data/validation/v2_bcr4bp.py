"""V2 long-span BCR4BP multi-lap drift gauntlet (#305 Part C — Option A).

Spec reinterpretation for BCR4BP periodic orbits (commensurate-only)
-------------------------------------------------------------------
Spec §14 V2 is ">=3 continuous laps; bounded drift in the defining rotating
frame". For a BCR4BP periodic orbit the load-bearing subtlety is the
two-frequency / quasi-periodicity problem (design draft §5): the model has two
incommensurate angular rates — the Earth-Moon synodic rate (normalized to 1)
and the Sun synodic rate ``omega_sun ~ 0.925196 rad/TU``. An orbit closes
strictly (``X(T) = X(0)`` at the SAME Sun phase) only if its period is
*Sun-commensurate* (``omega_sun*T = 2*pi*n``). For a NON-commensurate orbit,
after one nominal period the Sun is at a different phase, so ``X(T) != X(0)``
even for a perfect quasi-periodic orbit — the measured "drift" is the
Sun-phase mismatch, NOT instability. Reporting it as a V2 drift is a false
negative.

This module implements **Option A (commensurate-only)** from the design draft:

  * with ``theta_sun0`` fixed at the orbit's epoch, propagate ``k*T`` for
    ``k=1..n_cycles`` in the coherent BCR4BP and measure the position drift
    ``||X(kT)_pos - X(0)_pos||`` (km);
  * gate the max over k against a labelled drift floor;
  * if ``require_commensurate`` and the orbit is non-commensurate (its
    ``sun_phase_drift`` exceeds the V0 labelled convention), FAIL with reason
    ``non_commensurate_no_strict_period`` rather than reporting a meaningless
    drift.

The general (Option C, stroboscopic) path that admits quasi-periodic members is
a deferred follow-up (design draft §5/§6 step 5).

Floor rationale (labelled judgment call)
-----------------------------------------
``V2_BCR4BP_DRIFT_FLOOR_KMS = 50_000 km`` mirrors the 3D / moontour V2 drift
floor (:data:`cyclerfinder.data.validation.v2_3d.V2_DRIFT_FLOOR_KMS`) as the
starting judgment call. For a genuinely Sun-commensurate orbit that closed at
V1 (< 1 m/s), the position drift over 3 laps is many orders below this; a
breaking-up orbit (the Sun unbinds it) blows past it. The floor is a LABELLED
convention, not a sourced physical constant, exactly like the QP V2 floor.

Discipline
----------
* READ-ONLY on the BCR4BP genome (wrap, never re-solve).
* The commensurability gate is honest: a non-commensurate orbit FAILS with an
  explicit reason; it is NOT silently passed or silently failed on a bogus
  drift number.
* NO catalogue writeback.

References
----------
* spec §14 V2; design draft #305 §5 (two-frequency problem, Option A).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v0_bcr4bp import V0_BCR4BP_PHASE_DRIFT_CONVENTION
from cyclerfinder.data.validation.v1_bcr4bp import SEM_L_KM
from cyclerfinder.genome.bcr4bp_genome import BCR4BPPeriodicOrbit

V2_BCR4BP_N_CYCLES_MIN: Final[int] = 3
"""V2 minimum cycle count. Mirrors spec §14 ">=3 continuous laps". For a
BCR4BP periodic orbit a "lap" is the full period T (Sun-commensurate)."""

V2_BCR4BP_DRIFT_FLOOR_KMS: Final[float] = 50_000.0
"""V2 position-drift floor (km) over k laps. LABELLED judgment-call convention,
mirroring v2_3d's 50_000 km; recalibrate empirically. Not a sourced constant."""

_NON_COMMENSURATE_REASON: Final[str] = "non_commensurate_no_strict_period"
"""Verdict reason when an orbit has no strict multi-lap period (quasi-periodic)."""


@dataclass(frozen=True)
class V2VerdictBCR4BP:
    """Frozen V2 verdict for a BCR4BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id :
        Identifier carried for the audit trail.
    n_cycles_requested :
        Requested lap count (must be ``>= V2_BCR4BP_N_CYCLES_MIN``).
    n_cycles_propagated :
        Laps actually completed before the verdict formed. Equals
        ``n_cycles_requested`` unless the integrator failed or the orbit was
        rejected as non-commensurate (then 0).
    per_cycle_drift_km :
        Per-lap position drift ``||X(kT)_pos - X(0)_pos||`` (km), one element
        per completed lap. V3 chains on ``len(per_cycle_drift_km) >= n_cycles``.
    max_drift_km :
        Max of ``per_cycle_drift_km`` — the headline drift, gated against
        ``drift_floor_km``. ``inf`` if no lap completed.
    drift_floor_km :
        The labelled floor this verdict was held against.
    n_cycles_min :
        Spec minimum (3).
    require_commensurate :
        Whether the verdict required Sun-commensurability.
    is_commensurate :
        Whether the orbit's ``sun_phase_drift`` is below the V0 labelled
        convention (so a strict ``k*T`` lap is meaningful).
    converged_each_cycle :
        Whether every requested lap propagated successfully.
    passes_v2_bcr4bp :
        ``is_commensurate (if required) AND converged_each_cycle AND
        n_cycles_propagated >= n_cycles_min AND max_drift_km <= drift_floor_km``.
        Headline boolean.
    reason :
        Free-form reason string. ``non_commensurate_no_strict_period`` when the
        orbit was rejected for not being Sun-commensurate; ``""`` otherwise.
    notes :
        Free-form audit string.
    """

    candidate_id: str
    n_cycles_requested: int
    n_cycles_propagated: int
    per_cycle_drift_km: tuple[float, ...]
    max_drift_km: float
    drift_floor_km: float
    n_cycles_min: int
    require_commensurate: bool
    is_commensurate: bool
    converged_each_cycle: bool
    passes_v2_bcr4bp: bool
    reason: str
    notes: str = ""


def run_v2_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    n_cycles: int = V2_BCR4BP_N_CYCLES_MIN,
    drift_floor_km: float = V2_BCR4BP_DRIFT_FLOOR_KMS,
    require_commensurate: bool = True,
    phase_drift_convention: float = V0_BCR4BP_PHASE_DRIFT_CONVENTION,
    l_km: float = SEM_L_KM,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V2VerdictBCR4BP:
    """Run the V2 long-span BCR4BP multi-lap drift gauntlet (Option A).

    Pipeline:
      1. Determine commensurability: ``sun_phase_drift <= phase_drift_convention``.
         If ``require_commensurate`` and not commensurate, FAIL immediately with
         reason ``non_commensurate_no_strict_period`` (no drift propagation).
      2. With ``theta_sun0`` fixed at the orbit's epoch, propagate ``k*T`` for
         ``k=1..n_cycles`` from the IC in ONE shot per k (exposing instability
         amplification, not re-seeding from the previous lap).
      3. Record the position drift (km) per lap; gate the max against the floor.

    Parameters
    ----------
    candidate_id :
        Identifier carried into the verdict.
    orbit :
        A :class:`BCR4BPPeriodicOrbit` from the genome.
    n_cycles :
        Number of laps. Must be ``>= V2_BCR4BP_N_CYCLES_MIN``.
    drift_floor_km :
        Position-drift floor (default :data:`V2_BCR4BP_DRIFT_FLOOR_KMS`).
    require_commensurate :
        If True (default), a non-commensurate orbit FAILS with the explicit
        reason. If False, the drift is still measured (useful for a probe that
        wants the raw quasi-periodic recurrence number) but the verdict marks
        ``is_commensurate=False``.
    phase_drift_convention :
        Commensurability threshold (default :data:`V0_BCR4BP_PHASE_DRIFT_CONVENTION`).
    l_km :
        Length unit (km) for the drift conversion (default :data:`SEM_L_KM`).
    rtol, atol :
        Integrator tolerances.
    notes :
        Free-form audit note.

    Returns
    -------
    V2VerdictBCR4BP
        ``passes_v2_bcr4bp`` is the headline boolean.

    Notes
    -----
    A V2 PASS does NOT admit to the catalogue. Option A is commensurate-only;
    quasi-periodic members are honestly rejected (the stroboscopic Option C
    upgrade that admits them is a deferred follow-up).
    """
    if not isinstance(orbit, BCR4BPPeriodicOrbit):
        raise TypeError(f"orbit must be a BCR4BPPeriodicOrbit; got {type(orbit).__name__}")
    if n_cycles < V2_BCR4BP_N_CYCLES_MIN:
        raise ValueError(f"V2_bcr4bp requires n_cycles >= {V2_BCR4BP_N_CYCLES_MIN}; got {n_cycles}")
    if drift_floor_km <= 0.0:
        raise ValueError(f"drift_floor_km must be > 0; got {drift_floor_km}")
    if l_km <= 0.0:
        raise ValueError(f"l_km must be > 0; got {l_km}")

    is_commensurate = bool(orbit.sun_phase_drift <= phase_drift_convention)

    if require_commensurate and not is_commensurate:
        return V2VerdictBCR4BP(
            candidate_id=candidate_id,
            n_cycles_requested=int(n_cycles),
            n_cycles_propagated=0,
            per_cycle_drift_km=(),
            max_drift_km=float("inf"),
            drift_floor_km=float(drift_floor_km),
            n_cycles_min=V2_BCR4BP_N_CYCLES_MIN,
            require_commensurate=True,
            is_commensurate=False,
            converged_each_cycle=False,
            passes_v2_bcr4bp=False,
            reason=_NON_COMMENSURATE_REASON,
            notes=notes,
        )

    state0 = np.asarray(orbit.state_initial, dtype=np.float64)
    pos0 = state0[:3]
    period = float(orbit.period_nondim)

    per_cycle: list[float] = []
    converged = True
    n_done = 0
    for k in range(1, n_cycles + 1):
        try:
            arc = bcr4bp.propagate_bcr4bp(
                orbit.system, state0, k * period, with_stm=False, rtol=rtol, atol=atol
            )
        except RuntimeError:
            converged = False
            break
        pos_k = arc.state_f[:3]
        if not np.all(np.isfinite(pos_k)):
            converged = False
            break
        drift_km = float(np.linalg.norm(pos_k - pos0)) * l_km
        per_cycle.append(drift_km)
        n_done += 1

    max_drift = max(per_cycle) if per_cycle else float("inf")
    converged_each_cycle = bool(converged and n_done == n_cycles)
    passes = bool(
        is_commensurate
        and converged_each_cycle
        and n_done >= V2_BCR4BP_N_CYCLES_MIN
        and max_drift <= drift_floor_km
    )

    return V2VerdictBCR4BP(
        candidate_id=candidate_id,
        n_cycles_requested=int(n_cycles),
        n_cycles_propagated=int(n_done),
        per_cycle_drift_km=tuple(per_cycle),
        max_drift_km=float(max_drift),
        drift_floor_km=float(drift_floor_km),
        n_cycles_min=V2_BCR4BP_N_CYCLES_MIN,
        require_commensurate=bool(require_commensurate),
        is_commensurate=is_commensurate,
        converged_each_cycle=converged_each_cycle,
        passes_v2_bcr4bp=passes,
        reason="" if is_commensurate else _NON_COMMENSURATE_REASON,
        notes=notes,
    )


__all__ = [
    "V2_BCR4BP_DRIFT_FLOOR_KMS",
    "V2_BCR4BP_N_CYCLES_MIN",
    "V2VerdictBCR4BP",
    "run_v2_bcr4bp",
]
