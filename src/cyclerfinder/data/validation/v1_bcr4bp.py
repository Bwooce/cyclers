"""V1 same-model BCR4BP gauntlet (#305 Part A — RECOMMENDED FIRST BUILD).

Spec reinterpretation for BCR4BP periodic orbits
------------------------------------------------
The strict-periodic V1 in :mod:`v1_3d` asserts ``||X(T) - X(0)|| < 1 m/s`` as
the same-model closure gate, re-propagated with an integrator *other than* the
one that built the orbit. For a BCR4BP periodic orbit this carries over almost
verbatim: the model is non-autonomous (the Sun phase advances with time), so
the orbit closes on its *full* period under the same fixed ``theta_S0`` epoch
the corrector used. V1-BCR4BP re-propagates the converged IC under **Radau**
(independent of the corrector's DOP853) over the closure horizon (full period,
or half-period for a symmetric perpendicular-crossing residual) and gates the
masked closure residual against a nondim floor AND its km/s conversion.

Why this is the smallest, highest-confidence first build
--------------------------------------------------------
:func:`cyclerfinder.genome.bcr4bp_genome.correct_bcr4bp_periodic` already runs
an independent Radau full-period closure check at convergence and stores it on
``BCR4BPPeriodicOrbit.independent_closure_residual``. V1-BCR4BP *promotes and
freezes* that check into a standalone verdict with a sourced floor, AND
recomputes the masked closure FRESHLY under Radau (the QP-lane "fresh
confirmation" discipline) rather than re-quoting the corrector's number, so the
audit trail records an independently-evaluated result.

The two-frequency caveat (does NOT bite here)
---------------------------------------------
A BCR4BP orbit is strictly periodic only when its period is Sun-commensurate
(``omega_sun * T = 2*pi*n``). V1 validates a *single closed period* under a
fixed Sun epoch, so the quasi-periodicity issue (which bites V2 onward — see
:mod:`v2_bcr4bp`) does not block this tier. V1 reports
``sun_phase_drift`` for the audit trail but does not gate on it.

Floor rationale (sourced)
-------------------------
* ``V1_BCR4BP_FLOOR_NONDIM = 1e-6`` matches periodic V1
  (:data:`cyclerfinder.data.validation.v1_3d.V1_FLOOR_NONDIM_DEFAULT`) and the
  corrector's own ``independent_tol`` default.
* ``v1_floor_kms = 1e-3`` is the spec §14 V1 bar (1 m/s), reused verbatim from
  :data:`cyclerfinder.data.validation.v1_3d.V1_FLOOR_KMS`.

km/s conversion
---------------
:class:`cyclerfinder.core.bcr4bp.BCR4BPSystem` carries no length/time unit (it
is dimensionless). The conversion uses the registry's sourced Sun-Earth-Moon
units (``l_km = 384400`` km, ``tu_seconds = 375190`` s from
:data:`cyclerfinder.genome.bcr4bp_systems.SEM_ANDREU`), passed as parameters so
a non-SEM BCR4BP triple can supply its own units.

Discipline
----------
* READ-ONLY on the BCR4BP genome (wrap, never re-solve / re-correct here).
* The independent re-propagation is the whole point of V1 — single-integrator
  artefacts are exactly what V1 must catch.
* NO catalogue writeback. A V1 pass alone does NOT admit a BCR4BP family
  member; it is a *known-reproduction* candidate (Andreu lineage), flagged for
  human review, never self-admitted.

References
----------
* Rosales, Jorba et al. (2023), CeMDA 135:15 (BCR4BP parameters, POL1/POL2).
* spec §14 V1 (the 1 m/s floor and the independent-integrator requirement).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.genome.bcr4bp_genome import BCR4BPPeriodicOrbit

# Sourced Sun-Earth-Moon nondim units (Rosales-Jorba 2023 Table 3 / JPL SSD),
# mirroring cyclerfinder.genome.bcr4bp_systems.SEM_ANDREU. Used only for the
# km/s conversion of the closure residual; the gauntlet logic is dimensionless.
SEM_L_KM: Final[float] = 384400.0
"""Sun-Earth-Moon length unit (km): the Moon SMA (JPL SSD)."""

SEM_TU_SECONDS: Final[float] = 375190.0
"""Sun-Earth-Moon time unit (s): matches SEM_ANDREU.tu_seconds / run_303."""

V1_BCR4BP_FLOOR_NONDIM: Final[float] = 1.0e-6
"""V1 same-model nondim closure floor for BCR4BP periodic orbits.

Matches periodic V1's 1e-6 (:data:`v1_3d.V1_FLOOR_NONDIM_DEFAULT`) and the
corrector's ``independent_tol`` default. 1e-6 nondim ~ 0.38 km in EM units."""

V1_BCR4BP_FLOOR_KMS: Final[float] = 1.0e-3
"""Spec §14 V1 floor: same-model closure must close to < 1 m/s = 1e-3 km/s.

Reused verbatim from :data:`v1_3d.V1_FLOOR_KMS`; spec-fixed, NOT test-tunable."""


@dataclass(frozen=True)
class V1VerdictBCR4BP:
    """Frozen V1 verdict for a BCR4BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id :
        Identifier carried for the audit trail.
    corrector_residual_nondim :
        The corrector's own masked closure residual (DOP853 inside the Newton
        loop), echoed from the orbit. Reported nondim; a basic precondition.
    independent_closure_nondim :
        L2 norm of the masked closure residual ``X(T_event) - X(0)`` from a
        FRESH Radau re-propagation (rtol=atol=1e-12) at the SAME fixed Sun
        epoch (``theta_sun0``) the corrector used. The whole point of V1: a
        different integrator gates against single-integrator artefacts.
    corrector_residual_kms / independent_closure_kms :
        The two residuals converted to km/s via ``l_km / tu_seconds``.
    v1_floor_nondim / v1_floor_kms :
        The floors this verdict was held against (default
        :data:`V1_BCR4BP_FLOOR_NONDIM` / :data:`V1_BCR4BP_FLOOR_KMS`).
    sun_commensurate_n / sun_phase_drift :
        BCR4BP-specific bookkeeping echoed from the orbit. Recorded for the
        audit trail (V1 does NOT gate on commensurability — that is V2+).
    converged_corrector :
        Whether the corrector itself converged (precondition).
    converged_independent :
        Whether the FRESH Radau re-propagation closed below the nondim floor.
    passes_v1_bcr4bp :
        ``converged_corrector AND converged_independent AND
        independent_closure_kms <= v1_floor_kms``. Headline boolean.
    n_iter :
        Corrector Newton iterations consumed (echoed).
    is_half_period_residual :
        Whether the closure horizon is T/2 (symmetric perpendicular crossing)
        or T (full closure). Echoed so a reader knows what was propagated.
    notes :
        Free-form audit string.
    """

    candidate_id: str
    corrector_residual_nondim: float
    independent_closure_nondim: float
    corrector_residual_kms: float
    independent_closure_kms: float
    v1_floor_nondim: float
    v1_floor_kms: float
    sun_commensurate_n: int
    sun_phase_drift: float
    converged_corrector: bool
    converged_independent: bool
    passes_v1_bcr4bp: bool
    n_iter: int
    is_half_period_residual: bool
    notes: str = ""


def _velocity_unit_kms(l_km: float, tu_seconds: float) -> float:
    """Return the BCR4BP velocity unit in km/s: ``l_km / tu_seconds``."""
    if l_km <= 0.0 or tu_seconds <= 0.0:
        raise ValueError(
            f"invalid BCR4BP units for V1 km/s conversion: l_km={l_km} tu_seconds={tu_seconds}"
        )
    return float(l_km) / float(tu_seconds)


def _independent_radau_closure(
    orbit: BCR4BPPeriodicOrbit,
    *,
    rtol: float,
    atol: float,
) -> float:
    """Fresh Radau re-propagation of the converged IC over the closure horizon.

    Returns the L2 norm of the masked closure residual (over the orbit's
    ``residual_indices``) at ``T_event`` — full period or half-period per
    ``is_half_period_residual`` — propagated under Radau at the SAME fixed Sun
    epoch the corrector used. Returns ``inf`` on integrator failure.

    Distinct from ``orbit.independent_closure_residual`` (which is a FULL 6D
    closure at the full period) by recomputing the *masked* residual at the
    corrector's own closure horizon, so the V1 number is like-for-like with
    the corrector's convergence predicate but produced by a fresh integration.
    """
    t_event = 0.5 * orbit.period_nondim if orbit.is_half_period_residual else orbit.period_nondim
    try:
        sol = solve_ivp(
            bcr4bp.bcr4bp_eom,
            (0.0, t_event),
            np.asarray(orbit.state_initial, dtype=np.float64),
            args=(orbit.system,),
            method="Radau",
            rtol=max(rtol, 1e-12),
            atol=max(atol, 1e-12),
        )
    except (RuntimeError, ValueError):
        return float("inf")
    if not sol.success:
        return float("inf")
    diff = sol.y[:, -1] - np.asarray(orbit.state_initial, dtype=np.float64)
    masked = diff[list(orbit.residual_indices)]
    if not np.all(np.isfinite(masked)):
        return float("inf")
    return float(np.linalg.norm(masked))


def run_v1_bcr4bp(
    candidate_id: str,
    orbit: BCR4BPPeriodicOrbit,
    *,
    l_km: float = SEM_L_KM,
    tu_seconds: float = SEM_TU_SECONDS,
    corrector_floor_nondim: float = 1e-10,
    independent_floor_nondim: float = V1_BCR4BP_FLOOR_NONDIM,
    v1_floor_kms: float = V1_BCR4BP_FLOOR_KMS,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    notes: str = "",
) -> V1VerdictBCR4BP:
    """Run the V1 same-model BCR4BP gauntlet on a converged periodic orbit.

    Pipeline:
      1. Re-verify (never re-solve) the corrector residual against
         ``corrector_floor_nondim`` (basic precondition).
      2. FRESHLY re-propagate the corrected IC under Radau (independent of the
         corrector's DOP853) over the closure horizon, recomputing the masked
         closure residual.
      3. Convert both residuals to km/s via ``l_km / tu_seconds``.
      4. PASS iff the corrector converged AND the fresh Radau closure is below
         the nondim floor AND its km/s conversion is below the spec V1 floor.

    Parameters
    ----------
    candidate_id :
        Identifier carried into the verdict.
    orbit :
        A converged :class:`BCR4BPPeriodicOrbit` from the BCR4BP genome. V1
        re-verifies it; it does NOT re-correct.
    l_km, tu_seconds :
        BCR4BP nondim units for the km/s conversion. Default to the sourced
        Sun-Earth-Moon values (:data:`SEM_L_KM` / :data:`SEM_TU_SECONDS`).
    corrector_floor_nondim :
        Bar for the corrector's own residual (default 1e-10, the corrector tol).
    independent_floor_nondim :
        Bar for the fresh Radau closure (default :data:`V1_BCR4BP_FLOOR_NONDIM`).
    v1_floor_kms :
        Spec §14 V1 bar in km/s (default :data:`V1_BCR4BP_FLOOR_KMS`).
    rtol, atol :
        Integrator tolerances for the fresh Radau re-propagation.
    notes :
        Free-form audit note.

    Returns
    -------
    V1VerdictBCR4BP
        ``passes_v1_bcr4bp`` is the headline boolean.

    Notes
    -----
    A V1 PASS does NOT admit to the catalogue. For a published BCR4BP family
    (Andreu lineage) a passing member is a *known-reproduction* candidate to be
    flagged for human review with the proposed attribution, never self-admitted.
    """
    if not isinstance(orbit, BCR4BPPeriodicOrbit):
        raise TypeError(f"orbit must be a BCR4BPPeriodicOrbit; got {type(orbit).__name__}")
    if independent_floor_nondim <= 0.0:
        raise ValueError(f"independent_floor_nondim must be > 0; got {independent_floor_nondim}")
    if v1_floor_kms <= 0.0:
        raise ValueError(f"v1_floor_kms must be > 0; got {v1_floor_kms}")

    v_unit_kms = _velocity_unit_kms(l_km, tu_seconds)

    corrector_nondim = float(orbit.corrector_residual)
    independent_nondim = _independent_radau_closure(orbit, rtol=rtol, atol=atol)

    corrector_kms = corrector_nondim * v_unit_kms
    independent_kms = independent_nondim * v_unit_kms

    converged_corrector = bool(corrector_nondim <= corrector_floor_nondim)
    converged_independent = bool(independent_nondim <= independent_floor_nondim)
    passes_kms = bool(independent_kms <= v1_floor_kms)
    passes = bool(converged_corrector and converged_independent and passes_kms)

    return V1VerdictBCR4BP(
        candidate_id=candidate_id,
        corrector_residual_nondim=corrector_nondim,
        independent_closure_nondim=independent_nondim,
        corrector_residual_kms=corrector_kms,
        independent_closure_kms=independent_kms,
        v1_floor_nondim=float(independent_floor_nondim),
        v1_floor_kms=float(v1_floor_kms),
        sun_commensurate_n=int(orbit.sun_commensurate_n),
        sun_phase_drift=float(orbit.sun_phase_drift),
        converged_corrector=converged_corrector,
        converged_independent=converged_independent,
        passes_v1_bcr4bp=passes,
        n_iter=int(orbit.n_iter),
        is_half_period_residual=bool(orbit.is_half_period_residual),
        notes=notes,
    )


__all__ = [
    "SEM_L_KM",
    "SEM_TU_SECONDS",
    "V1_BCR4BP_FLOOR_KMS",
    "V1_BCR4BP_FLOOR_NONDIM",
    "V1VerdictBCR4BP",
    "run_v1_bcr4bp",
]
