"""V2 long-span bounded-drift 3D gauntlet (#306 Phase 1 Part B).

Spec reference
--------------
* §14 V2-ballistic — "≥3 continuous laps; **bounded** drift in the dynamic
  rotating frame (tolerant of geometric breathing), evaluated **in the row's
  defining model**".

For a CR3BP periodic 3D orbit, V2 means: propagate the IC for ``n_cycles``
consecutive periods in the SAME CR3BP model (Radau at ``rtol=atol=1e-12``),
and assert that at each return-to-period the state stays within a bounded
distance of the original IC. The max position drift across all cycle
boundaries is the headline number.

The drift floor (50,000 km) is the same-model V2-ballistic bar this project
already uses elsewhere (sourced from
:data:`cyclerfinder.verify.propagate.DRIFT_TOLERANCE_KM`). It is generous
enough to absorb hyperbolic-instability amplification of a corrector-clean
IC over 3 cycles AND tight enough to reject an orbit that has effectively
escaped its family. For real-ephemeris bridge calls, callers may override
via ``drift_floor_kms`` to the spec §14 V2-real number (200,000 km).

Composition map
---------------
This module composes directly on
:func:`cyclerfinder.core.cr3bp.propagate` (the 3D CR3BP propagator, 6D
state). Per the spec it does NOT need a separate maintenance step — V2
is the *unmaintained* periodicity gate; a powered cycler is a different
gate (V2-powered, see ``v2_powered.py``).

What V2 is NOT
--------------
* V1 (same-model closure of a single period) — see :mod:`v1_3d`.
* V2-powered (per-cycle retargeted) — see
  :mod:`cyclerfinder.verify.v2_powered`.
* V3 (real-ephemeris) — see :mod:`cyclerfinder.verify.real_closure`.
* V4 (HFEM real-eph) / V5 (mission quality).

Discipline
----------
* NO catalogue writeback. V2 passes alone do not admit.
* The 3-cycle minimum is spec §14, NOT test-tunable.
* The drift floor is a module constant; an override is for caller
  reflexivity (e.g. a real-ephemeris V2-bridge call) only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp

V2_N_CYCLES_MIN: Final[int] = 3
"""Spec §14 V2-ballistic minimum: ``>= 3`` continuous laps.

Module constant; NOT test-tunable (spec-fixed)."""

V2_DRIFT_FLOOR_KMS: Final[float] = 50_000.0
"""Default same-model drift floor in km.

Sourced from :data:`cyclerfinder.verify.propagate.DRIFT_TOLERANCE_KM` —
the V2-ballistic same-model bar this project already uses for E-M class
cyclers (~0.02 deg of geometric breathing per lap at Mars's ~1.5 AU =
~50,000 km / 0.02 deg). Tight enough to reject propagator regressions
and frame-transform errors; loose enough to absorb the natural
hyperbolic-instability amplification of a corrector-clean IC over 3
cycles (for the Braik-Ross C11a family, Floquet abs ~300 means a
1e-9-nondim closure amplifies to ~1e3 nondim ~= 4e5 km by cycle 3 —
which is precisely why §14 V2 says BOUNDED, not VANISHING).

Callers running V2 as the real-ephemeris bridge may override with the
spec §14 V2-real number
(:data:`cyclerfinder.verify.real_closure.REAL_DRIFT_TOLERANCE_KM` =
200,000 km), which absorbs real-ephemeris breathing on top of the
hyperbolic amplification."""


@dataclass(frozen=True)
class V2Verdict3D:
    """Frozen V2 verdict for a 3D CR3BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id:
        Identifier of the candidate (carried for the audit trail).
    n_cycles_requested:
        The ``n_cycles`` requested by the caller (>= :data:`V2_N_CYCLES_MIN`).
    n_cycles_propagated:
        Cycles the integrator actually completed before the verdict was
        formed. Equals ``n_cycles_requested`` unless the propagator failed
        mid-flight, in which case it's the count of successful returns.
    per_cycle_drift_kms:
        Cumulative position drift at each cycle boundary, km. Length equals
        ``n_cycles_propagated``. Element ``k`` is ``||X_pos((k+1)*T) -
        X_pos(0)||`` in km (NOT the cycle-to-cycle delta — the cumulative
        drift from the original IC, matching :attr:`max_drift_kms`).
    max_drift_kms:
        The max of ``per_cycle_drift_kms``. The headline drift number.
    per_cycle_velocity_drift_kms:
        Cumulative velocity drift at each cycle boundary, km/s. Diagnostic;
        does not gate the verdict (V2 is a POSITION drift gate per spec —
        velocity drift is reported for honesty).
    drift_floor_kms:
        The bar this verdict was held against. Stored for audit.
    n_cycles_min:
        The spec §14 minimum cycles (3). Stored for audit.
    converged_at_each_return:
        Whether the propagator successfully completed every requested
        return. Basic precondition for a credible V2 pass.
    passes_v2:
        ``converged_at_each_return AND n_cycles_propagated >= n_cycles_min
        AND max_drift_kms <= drift_floor_kms``. The headline boolean.
    notes:
        Free-form audit string.
    """

    candidate_id: str
    n_cycles_requested: int
    n_cycles_propagated: int
    per_cycle_drift_kms: tuple[float, ...]
    max_drift_kms: float
    per_cycle_velocity_drift_kms: tuple[float, ...]
    drift_floor_kms: float
    n_cycles_min: int
    converged_at_each_return: bool
    passes_v2: bool
    notes: str = ""


def run_v2_3d(
    candidate_id: str,
    state0: NDArray[np.float64],
    period_nondim: float,
    system: cr3bp.CR3BPSystem,
    *,
    n_cycles: int = V2_N_CYCLES_MIN,
    drift_floor_kms: float = V2_DRIFT_FLOOR_KMS,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V2Verdict3D:
    """Run the V2 long-span bounded-drift 3D gauntlet on a CR3BP IC.

    Pipeline:
      1. Propagate the IC for ``n_cycles`` consecutive periods using the
         6D CR3BP propagator (DOP853 at the requested tolerances, the
         integrator's stiff-tolerant default; same integrator/tolerances
         the corrector uses inside the Newton loop).
      2. At each cycle boundary record the cumulative position + velocity
         delta from the original IC.
      3. Convert position drift to km via the system's length unit; assert
         the max stays within ``drift_floor_kms``.

    Parameters
    ----------
    candidate_id:
        Identifier carried into the verdict.
    state0:
        6-vector IC in the CR3BP rotating frame, nondim.
    period_nondim:
        Full nondim period of the orbit (typically the corrector-refined
        period from V1's :class:`~v1_3d.V1Verdict3D`).
    system:
        CR3BP system. Used for the EOM (``system.mu``) AND the km
        conversion (``system.l_km``).
    n_cycles:
        Number of consecutive cycles to propagate. Must be
        ``>= V2_N_CYCLES_MIN`` (spec §14). Default 3.
    drift_floor_kms:
        Bar against which ``max_drift_kms`` is gated. Default
        :data:`V2_DRIFT_FLOOR_KMS` (50,000 km, same-model). A
        real-ephemeris bridge call may pass the V2-real floor instead.
    rtol, atol:
        Propagator tolerances. Defaults at 1e-12 match V1 / the corrector
        for like-for-like comparison.
    notes:
        Free-form audit note.

    Returns
    -------
    V2Verdict3D
        ``passes_v2`` is the headline.

    Notes
    -----
    The propagation is sequential cycle-by-cycle (not one shot of length
    ``n_cycles * period``) so that an intermediate failure surfaces in
    ``n_cycles_propagated`` rather than poisoning the whole verdict.

    A V2 PASS does NOT admit to the catalogue. V3+V4 follow.
    """
    state0_arr = np.asarray(state0, dtype=np.float64)
    if state0_arr.shape != (6,):
        raise ValueError(f"state0 must be 6D; got shape {state0_arr.shape}")
    if period_nondim <= 0.0:
        raise ValueError(f"period_nondim must be > 0; got {period_nondim}")
    if n_cycles < V2_N_CYCLES_MIN:
        raise ValueError(f"V2 requires n_cycles >= {V2_N_CYCLES_MIN} (spec §14); got {n_cycles}")
    if drift_floor_kms <= 0.0:
        raise ValueError(f"drift_floor_kms must be > 0; got {drift_floor_kms}")
    if system.l_km <= 0.0 or system.t_s <= 0.0:
        raise ValueError(
            f"invalid CR3BP system for V2 km conversion: l_km={system.l_km} t_s={system.t_s}"
        )

    l_km = float(system.l_km)
    v_unit_kms = l_km / float(system.t_s)

    pos_drifts_kms: list[float] = []
    vel_drifts_kms: list[float] = []
    current = state0_arr.copy()
    converged = True
    n_done = 0

    for _ in range(n_cycles):
        try:
            arc = cr3bp.propagate(
                system,
                current,
                period_nondim,
                with_stm=False,
                rtol=rtol,
                atol=atol,
            )
        except RuntimeError:
            converged = False
            break
        n_done += 1
        current = np.asarray(arc.state_f, dtype=np.float64)
        # Cumulative drift from the ORIGINAL IC, in km.
        pos_delta_nondim = float(np.linalg.norm(current[:3] - state0_arr[:3]))
        vel_delta_nondim = float(np.linalg.norm(current[3:] - state0_arr[3:]))
        pos_drifts_kms.append(pos_delta_nondim * l_km)
        vel_drifts_kms.append(vel_delta_nondim * v_unit_kms)

    max_drift = max(pos_drifts_kms) if pos_drifts_kms else float("inf")

    passes_v2 = bool(converged and n_done >= V2_N_CYCLES_MIN and max_drift <= drift_floor_kms)

    return V2Verdict3D(
        candidate_id=candidate_id,
        n_cycles_requested=int(n_cycles),
        n_cycles_propagated=int(n_done),
        per_cycle_drift_kms=tuple(pos_drifts_kms),
        max_drift_kms=float(max_drift),
        per_cycle_velocity_drift_kms=tuple(vel_drifts_kms),
        drift_floor_kms=float(drift_floor_kms),
        n_cycles_min=V2_N_CYCLES_MIN,
        converged_at_each_return=bool(converged and n_done == n_cycles),
        passes_v2=passes_v2,
        notes=notes,
    )


__all__ = [
    "V2_DRIFT_FLOOR_KMS",
    "V2_N_CYCLES_MIN",
    "V2Verdict3D",
    "run_v2_3d",
]
