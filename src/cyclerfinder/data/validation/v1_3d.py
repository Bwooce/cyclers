"""V1 same-model 3D gauntlet (#306 Phase 1 Part A).

Spec reference
--------------
* §14 V1 — "every leg re-solved with **lamberthub izzo + gooding**,
  agreement < 1e-3 m/s; full trajectory re-propagated with the **Kepler**
  propagator (not the Lambert that built it), planet positions met < tol."

For full-3D CR3BP periodic orbits (no Lambert legs to re-solve), V1 means
**re-close the IC under the 3D corrector AND independently re-propagate**
with a different integrator (Radau at rtol=atol=1e-12). Both must hold.

The km/s floor (1e-3) is the spec §14 V1 number; the corrector gates against
its nondim floor (default ``1e-6`` for the independent Radau cross-check,
matching :func:`correct_general_periodic_3d`'s ``independent_tol``).

Composition map
---------------
This module is a thin layer over the #291 Phase 1
:func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
(commit ``6da8ee1``), which already:

* runs single-shooting Newton on the STM to land the orbit
  (``corrector_residual``);
* independently re-propagates the corrected IC under Radau at the same
  tolerances and reports ``independent_closure_residual``;
* emits both numbers on the returned :class:`Periodic3DOrbit`.

V1_3D just wraps both gates into a single verdict, converts the nondim
residuals to km/s via the system's velocity unit ``l_km / t_s``, and asserts
the spec §14 1 m/s floor.

What V1 is NOT
--------------
V1 is a SAME-MODEL gate — it asserts the IC is genuinely periodic under the
*same CR3BP model* the search produced it in. V1 is silent on:

* fidelity persistence under a model widening (that is Axis B / §14 V3);
* whether the orbit is sourced (that is Axis C / §14 V3 corroboration);
* independent-tool falsification panels (Axis D);
* long-span periodicity over ≥3 cycles (that is V2 — see :mod:`v2_3d`);
* real-ephemeris realisation (V3-ballistic / V3-powered);
* HFEM real-eph closure (V4); mission quality (V5).

Discipline
----------
* NO catalogue writeback. A V1 pass DOES NOT admit to ``catalogue.yaml`` —
  that requires V3+V4 minimum plus the #328 lit-check verdict.
* The V1 floor (1 m/s) is spec §14, NOT test-tunable.
* The independent re-propagation is the whole point of V1 — single-integrator
  artefacts are exactly what V1 must catch.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)

V1_FLOOR_KMS: Final[float] = 1.0e-3
"""Spec §14 V1 floor: same-model closure must close to < 1 m/s = 1e-3 km/s.

Module constant; NOT test-tunable (spec-fixed)."""

V1_FLOOR_NONDIM_DEFAULT: Final[float] = 1.0e-6
"""Default nondim corrector-independent floor (matches
:func:`correct_general_periodic_3d` ``independent_tol`` default). The
corrector-side floor is checked in nondim; the km/s floor is the spec V1 bar
applied after conversion."""


@dataclass(frozen=True)
class V1Verdict3D:
    """Frozen V1 verdict for a 3D CR3BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id:
        Identifier of the candidate (carried for the audit trail).
    closure_residual_nondim:
        L2 norm of the corrector's masked closure residual (DOP853 inside the
        Newton loop). Reported nondim.
    independent_closure_nondim:
        L2 norm of the full 6D state difference ``X(T) - X(0)`` from an
        independent Radau re-propagation at ``rtol=atol=1e-12``. The whole
        point of V1: a different integrator gates against single-integrator
        artefacts.
    closure_residual_kms:
        ``closure_residual_nondim`` converted to km/s via the system's
        velocity unit ``l_km / t_s``. Reported for human + spec-floor compare.
    independent_closure_kms:
        ``independent_closure_nondim`` converted to km/s the same way.
    v1_floor_kms:
        The spec §14 V1 bar (1e-3 km/s). Stored for audit so a later reader
        sees exactly what the verdict was held against.
    independent_floor_nondim:
        The nondim floor the independent re-propagation was held against
        (default :data:`V1_FLOOR_NONDIM_DEFAULT`).
    converged_corrector:
        Whether the corrector itself converged (``corrector_residual`` <
        ``tol``) — basic precondition.
    converged_independent:
        Whether the independent re-propagation closed below
        ``independent_floor_nondim``.
    passes_v1:
        ``converged_corrector AND converged_independent AND
        independent_closure_kms <= v1_floor_kms``. The headline boolean.
    n_iter:
        Corrector Newton iterations consumed.
    degenerate_planar:
        ``True`` iff the corrected IC collapsed to the planar manifold
        (``|z0|`` and ``|zdot0|`` both below the corrector's ``planar_eps``).
        A genuinely 3D V1 verdict must have ``degenerate_planar=False``.
        Reported for diagnostics; does NOT veto the verdict (a planar
        member is still validly closed — caller chooses the family).
    notes:
        Free-form audit string.
    """

    candidate_id: str
    closure_residual_nondim: float
    independent_closure_nondim: float
    closure_residual_kms: float
    independent_closure_kms: float
    v1_floor_kms: float
    independent_floor_nondim: float
    converged_corrector: bool
    converged_independent: bool
    passes_v1: bool
    n_iter: int
    degenerate_planar: bool
    notes: str = ""


def _velocity_unit_kms(system: cr3bp.CR3BPSystem) -> float:
    """Return the CR3BP velocity unit in km/s: ``l_km / t_s``.

    Raised on a malformed system to fail loudly rather than silently divide
    by zero (a fabricated system would otherwise give a meaningless verdict).
    """
    if system.l_km <= 0.0 or system.t_s <= 0.0:
        raise ValueError(
            f"invalid CR3BP system for V1 km/s conversion: l_km={system.l_km} t_s={system.t_s}"
        )
    return float(system.l_km) / float(system.t_s)


def run_v1_3d(
    candidate_id: str,
    state0: NDArray[np.float64],
    period_nondim: float,
    system: cr3bp.CR3BPSystem,
    *,
    free_vars: tuple[int, ...] = FREE_VARS_FULL_ASYMMETRIC,
    residual_indices: tuple[int, ...] = RESIDUAL_FULL_STATE_AT_T,
    is_half_period_residual: bool = False,
    corrector_tol: float = 1e-10,
    independent_floor_nondim: float = V1_FLOOR_NONDIM_DEFAULT,
    v1_floor_kms: float = V1_FLOOR_KMS,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    max_iter: int = 60,
    notes: str = "",
) -> V1Verdict3D:
    """Run the V1 same-model 3D gauntlet on a CR3BP periodic-orbit IC.

    Pipeline:
      1. Re-close the IC under the #291 Phase 1 3D corrector
         (single-shooting Newton on the STM, DOP853).
      2. Independently re-propagate the corrected IC under Radau at the same
         rtol/atol. (The corrector already does this; we just wrap the
         numbers into the verdict.)
      3. Convert both nondim residuals to km/s via the system's velocity unit.
      4. Assert the spec §14 V1 floor: ``independent_closure_kms <=
         v1_floor_kms``.

    Parameters
    ----------
    candidate_id:
        Identifier carried into the verdict for the audit trail.
    state0:
        6-vector IC ``(x, y, z, xdot, ydot, zdot)``, nondim CR3BP rotating frame.
    period_nondim:
        Initial guess for the full nondim period (passed straight to the
        corrector; the period is in the default free-var bundle so the
        corrector may refine it).
    system:
        CR3BP system. Used for the EOM (via ``system.mu``) AND to convert
        the nondim residuals into km/s.
    free_vars, residual_indices, is_half_period_residual:
        Passed through to :func:`correct_general_periodic_3d`. Defaults to
        full asymmetric closure (the 7-unknown 6-residual mode).
    corrector_tol:
        Corrector Newton convergence tolerance (nondim L2 of the residual).
        Default 1e-10 matches the corrector's own default.
    independent_floor_nondim:
        Nondim bar for the independent re-propagation. Default 1e-6 matches
        the corrector's ``independent_tol`` default; 1e-6 nondim ~384 m in
        Earth-Moon units, ~580 m in Uranus-Oberon.
    v1_floor_kms:
        Spec §14 V1 bar in km/s. Default 1e-3 (1 m/s).
    rtol, atol:
        Integrator tolerances for BOTH the corrector and the independent
        re-propagation (the corrector applies the same tolerances to its
        Radau cross-check).
    max_iter:
        Corrector Newton iteration cap.
    notes:
        Free-form audit note carried into the verdict.

    Returns
    -------
    V1Verdict3D
        The frozen V1 verdict. ``passes_v1`` is the headline boolean.

    Notes
    -----
    A V1 PASS does NOT admit to the catalogue. V1 is the cheapest gate in the
    §14 ladder; V2 + V3 + V4 (+ V5) and Axis C (#328 lit-check) must all
    follow before any writeback.

    For genuinely 3D orbits the caller should screen on
    ``verdict.degenerate_planar=False`` AFTER the verdict; a corrector that
    collapsed to the planar manifold is mathematically valid (the manifold
    is invariant) but is NOT the 3D candidate that was searched for.
    """
    state0_arr = np.asarray(state0, dtype=np.float64)
    if state0_arr.shape != (6,):
        raise ValueError(f"state0 must be 6D; got shape {state0_arr.shape}")
    if period_nondim <= 0.0:
        raise ValueError(f"period_nondim must be > 0; got {period_nondim}")
    if v1_floor_kms <= 0.0:
        raise ValueError(f"v1_floor_kms must be > 0; got {v1_floor_kms}")
    if independent_floor_nondim <= 0.0:
        raise ValueError(f"independent_floor_nondim must be > 0; got {independent_floor_nondim}")

    v_unit_kms = _velocity_unit_kms(system)

    orbit: Periodic3DOrbit = correct_general_periodic_3d(
        system,
        state0_arr,
        period_nondim,
        free_vars=free_vars,
        residual_indices=residual_indices,
        is_half_period_residual=is_half_period_residual,
        tol=corrector_tol,
        max_iter=max_iter,
        rtol=rtol,
        atol=atol,
        independent_tol=independent_floor_nondim,
    )

    closure_kms = float(orbit.corrector_residual) * v_unit_kms
    independent_kms = float(orbit.independent_closure_residual) * v_unit_kms

    converged_corrector = bool(orbit.corrector_residual <= corrector_tol)
    converged_independent = bool(orbit.independent_closure_residual <= independent_floor_nondim)
    # The km/s spec floor applies to the INDEPENDENT closure (the same-model
    # re-propagation residual — what a downstream consumer of the IC would
    # observe). The corrector residual is reported for diagnostics but is
    # an internal Newton number, not the spec quantity.
    passes_kms_floor = bool(independent_kms <= v1_floor_kms)
    passes_v1 = bool(converged_corrector and converged_independent and passes_kms_floor)

    return V1Verdict3D(
        candidate_id=candidate_id,
        closure_residual_nondim=float(orbit.corrector_residual),
        independent_closure_nondim=float(orbit.independent_closure_residual),
        closure_residual_kms=closure_kms,
        independent_closure_kms=independent_kms,
        v1_floor_kms=v1_floor_kms,
        independent_floor_nondim=independent_floor_nondim,
        converged_corrector=converged_corrector,
        converged_independent=converged_independent,
        passes_v1=passes_v1,
        n_iter=int(orbit.n_iter),
        degenerate_planar=bool(orbit.degenerate_planar),
        notes=notes,
    )


__all__ = [
    "V1_FLOOR_KMS",
    "V1_FLOOR_NONDIM_DEFAULT",
    "V1Verdict3D",
    "run_v1_3d",
]
