r"""Sims-Flanagan low-thrust leg optimiser (Phase 3 of the v2 low-thrust scope).

Wires the pure leg model in :mod:`cyclerfinder.core.sims_flanagan` into the M5
optimiser pattern (scipy ``differential_evolution`` global pass + SLSQP local
polish, mirroring :mod:`cyclerfinder.search.optimize` and
:mod:`cyclerfinder.search.maintain`). The decision vector is the flattened
per-segment ``Delta V`` schedule (``3N`` for an ``N``-segment leg); the
match-point defect is the equality constraint driven to zero, and the
per-segment ``Delta V_max`` thrust capability is the inequality constraint.

Two-phase solve (Yam §1, recorded in ``docs/v2-future-references.md`` §1):

* **Phase 1** — :func:`solve_leg_min_dv` minimises the total ``Delta V``
  (``min Sum |Delta V_i|``, Yam Eq. 4) subject to the defect equality and the
  thrust-capability inequalities.
* **Phase 2** — :func:`solve_leg_max_mass` re-optimises *locally* from the
  Phase-1 solution to maximise the final spacecraft mass (Yam Eq. 6), with mass
  propagated by the rocket equation (Yam Eq. 5). Yam used SNOPT; we use SLSQP to
  match the existing stack — what Yam validates is the *transcription*, not the
  specific solver.

Golden discipline: there is no usable literature anchor for the leg model (the
Yam worked examples are out-of-scope Jupiter rendezvous; see the plan's Phase 5
section). The optimiser is validated on physics invariants only — a
ballistically-closed leg solves to ~zero ``Delta V``; a perturbed leg's defect
is driven below tolerance; the converged schedule respects the thrust bound; and
the max-mass phase never reduces the final mass below the min-ΔV phase.

Plan: ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md`` (Phase 3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import NonlinearConstraint, differential_evolution, minimize

from cyclerfinder.core.kepler import KeplerConvergenceError
from cyclerfinder.core.sims_flanagan import (
    SimsFlanaganError,
    SimsFlanaganLeg,
    final_mass,
    leg_feasible,
    match_point_defect,
    segment_dv_bounds,
)

# ---------------------------------------------------------------------------
# Module-level constants / tolerances
# ---------------------------------------------------------------------------

_DEFECT_POS_SCALE_KM: float = 1.0e6
"""Position-defect normaliser (km) so the equality constraint vector is
dimensionless and comparably scaled with the velocity block. 1e6 km ~ a few
times the Earth-Moon distance; the optimiser drives the raw defect well below
this, but the scaling keeps SLSQP's constraint Jacobian well-conditioned."""

_DEFECT_VEL_SCALE_KMS: float = 1.0
"""Velocity-defect normaliser (km/s). Heliocentric speeds are O(10) km/s, so a
unit scale keeps the velocity block O(1) without dominating the position one."""

_DE_MAXITER: int = 80
"""Scipy ``differential_evolution`` generation cap for the 3N-dim ΔV box."""

_DE_POPSIZE: int = 10
"""Scipy DE population multiplier."""

_DE_TOL: float = 1.0e-6
"""Scipy DE convergence tolerance (population stdev fraction)."""

_SLSQP_MAXITER: int = 300
"""Scipy SLSQP iteration cap."""

_SLSQP_FTOL: float = 1.0e-10
"""Scipy SLSQP objective tolerance — tight, since the objective (ΔV) is small."""

_FEASIBLE_POS_TOL_KM: float = 1.0
"""Position-defect tolerance (km) for the ``converged`` flag."""

_FEASIBLE_VEL_TOL_KMS: float = 1.0e-3
"""Velocity-defect tolerance (km/s) for the ``converged`` flag."""


@dataclass(frozen=True)
class LowThrustLegResult:
    r"""Outcome of one low-thrust leg solve.

    Attributes
    ----------
    dvs:
        The optimised per-segment ``Delta V`` schedule, ``(N, 3)`` km/s.
    total_dv_kms:
        Sum of per-segment ``Delta V`` magnitudes, km/s (the Phase-1
        objective).
    final_mass_kg:
        Spacecraft mass at the leg end under ``dvs``, kg (the Phase-2
        objective, negated for minimisation internally).
    converged:
        ``True`` iff the match-point defect at ``dvs`` is within the
        feasibility tolerances (position + velocity).
    defect_pos_km:
        Position-defect norm at ``dvs``, km.
    defect_vel_kms:
        Velocity-defect norm at ``dvs``, km/s.
    """

    dvs: NDArray[np.float64]
    total_dv_kms: float
    final_mass_kg: float
    converged: bool
    defect_pos_km: float
    defect_vel_kms: float


_DEFECT_PENALTY: float = 1.0e6
"""Large scaled-defect value returned when a trial schedule exhausts the
spacecraft mass (the rocket equation underflows ``final_mass`` to ~0 under an
absurdly large ΔV the optimiser may probe). Finite so SLSQP / DE see a smooth
penalty rather than an exception."""


def _defect_equality(leg: SimsFlanaganLeg, flat: NDArray[np.float64]) -> NDArray[np.float64]:
    """Scaled position+velocity defect 6-vector for the flattened schedule.

    The mass block is omitted from the equality constraint: the backward pass is
    seeded with the self-consistent :func:`final_mass`, so the mass defect is
    identically zero and would only add a redundant (always-satisfied) row.

    Returns a large finite penalty vector when the trial ΔV exhausts the
    spacecraft mass (``final_mass`` underflows so the backward propagation has no
    valid end mass) — the optimiser sees a smooth penalty, never an exception.
    """
    dvs = flat.reshape(leg.n_segments, 3)
    mf = final_mass(leg, dvs)
    if not np.isfinite(mf) or mf <= 0.0:
        return np.full(6, _DEFECT_PENALTY, dtype=np.float64)
    try:
        defect = match_point_defect(leg, dvs, mf)
    except (KeplerConvergenceError, SimsFlanaganError):
        # A giant trial ΔV can drive a coast arc to a near-parabolic state the
        # universal-variable Newton can't resolve; treat as a penalty point.
        return np.full(6, _DEFECT_PENALTY, dtype=np.float64)
    scaled = np.empty(6, dtype=np.float64)
    scaled[0:3] = defect[0:3] / _DEFECT_POS_SCALE_KM
    scaled[3:6] = defect[3:6] / _DEFECT_VEL_SCALE_KMS
    return scaled


def _thrust_slacks(leg: SimsFlanaganLeg, flat: NDArray[np.float64]) -> NDArray[np.float64]:
    """Per-segment thrust-capability slack ``bound_i - |Delta V_i|`` (km/s).

    Non-negative ⇒ within capability (SLSQP ``fun(x) >= 0`` convention).
    """
    dvs = flat.reshape(leg.n_segments, 3)
    bounds = segment_dv_bounds(leg, dvs)
    mags = np.linalg.norm(dvs, axis=1)
    return cast(NDArray[np.float64], bounds - mags)


def _total_dv(flat: NDArray[np.float64], leg: SimsFlanaganLeg) -> float:
    """Phase-1 objective: total ΔV (sum of per-segment magnitudes), km/s."""
    dvs = flat.reshape(leg.n_segments, 3)
    return float(np.sum(np.linalg.norm(dvs, axis=1)))


def _neg_final_mass(flat: NDArray[np.float64], leg: SimsFlanaganLeg) -> float:
    """Phase-2 objective: negative final mass (minimised ⇒ mass maximised)."""
    dvs = flat.reshape(leg.n_segments, 3)
    return -final_mass(leg, dvs)


def _dv_box(leg: SimsFlanaganLeg) -> list[tuple[float, float]]:
    """Per-component ΔV box bounds for DE / SLSQP.

    Each component is bounded by the largest per-segment capability (the bound
    grows as mass falls, so the *all-burns* profile gives the loosest valid
    box). A symmetric ``[-cap, cap]`` per component is a superset of the
    spherical ``|Delta V_i| <= cap`` feasible region; the nonlinear thrust
    inequality enforces the exact spherical bound.
    """
    zero = np.zeros((leg.n_segments, 3), dtype=np.float64)
    cap = float(np.max(segment_dv_bounds(leg, zero)))
    # When mass falls the cap rises; use a generous multiple so the box never
    # clips a feasible spherical-bound solution. The exact bound is enforced by
    # the nonlinear inequality.
    cap_box = cap * 2.0 if cap > 0.0 else 0.0
    return [(-cap_box, cap_box) for _ in range(3 * leg.n_segments)]


def _make_constraints(leg: SimsFlanaganLeg) -> list[dict[str, Any]]:
    """SLSQP constraint dicts: defect equality + per-segment thrust inequality."""
    constraints: list[dict[str, Any]] = [
        {"type": "eq", "fun": lambda flat: _defect_equality(leg, np.asarray(flat))},
        {"type": "ineq", "fun": lambda flat: _thrust_slacks(leg, np.asarray(flat))},
    ]
    return constraints


def _evaluate(leg: SimsFlanaganLeg, flat: NDArray[np.float64]) -> LowThrustLegResult:
    """Build a :class:`LowThrustLegResult` from a flattened schedule."""
    dvs = flat.reshape(leg.n_segments, 3).copy()
    mf = final_mass(leg, dvs)
    if not np.isfinite(mf) or mf <= 0.0:
        return LowThrustLegResult(
            dvs=dvs,
            total_dv_kms=_total_dv(flat, leg),
            final_mass_kg=0.0,
            converged=False,
            defect_pos_km=float("inf"),
            defect_vel_kms=float("inf"),
        )
    try:
        defect = match_point_defect(leg, dvs, mf)
    except (KeplerConvergenceError, SimsFlanaganError):
        return LowThrustLegResult(
            dvs=dvs,
            total_dv_kms=_total_dv(flat, leg),
            final_mass_kg=mf,
            converged=False,
            defect_pos_km=float("inf"),
            defect_vel_kms=float("inf"),
        )
    pos = float(np.linalg.norm(defect[0:3]))
    vel = float(np.linalg.norm(defect[3:6]))
    converged = pos <= _FEASIBLE_POS_TOL_KM and vel <= _FEASIBLE_VEL_TOL_KMS
    # Also require the thrust bound to hold (a feasible solve must respect it).
    if converged:
        converged = bool(np.all(_thrust_slacks(leg, flat) >= -1.0e-9))
    return LowThrustLegResult(
        dvs=dvs,
        total_dv_kms=_total_dv(flat, leg),
        final_mass_kg=final_mass(leg, dvs),
        converged=converged,
        defect_pos_km=pos,
        defect_vel_kms=vel,
    )


def _slsqp_polish(
    leg: SimsFlanaganLeg,
    x0: NDArray[np.float64],
    objective: Any,
    bounds: list[tuple[float, float]],
) -> NDArray[np.float64]:
    """One SLSQP local solve honouring the defect + thrust constraints."""
    constraints = _make_constraints(leg)
    minimize_any = cast(Any, minimize)
    result_any: Any = minimize_any(
        objective,
        x0,
        args=(leg,),
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": _SLSQP_MAXITER, "ftol": _SLSQP_FTOL},
    )
    return cast(NDArray[np.float64], np.asarray(result_any.x, dtype=np.float64))


def solve_leg_min_dv(
    leg: SimsFlanaganLeg,
    *,
    seed: int = 0,
    n_starts: int = 3,
    use_de: bool = True,
) -> LowThrustLegResult:
    r"""Phase 1: minimum-``Delta V`` low-thrust leg solve (Yam Eq. 4).

    Minimises ``Sum |Delta V_i|`` over the per-segment ``Delta V`` schedule
    subject to (a) the match-point defect equality (position + velocity → 0) and
    (b) the per-segment thrust-capability inequality. Uses an optional scipy
    ``differential_evolution`` global pass followed by SLSQP polishes from the
    zero schedule and ``n_starts - 1`` seeded random starts (the M5 pattern).

    A coast-only leg (``tmax_kn == 0``) admits only the zero schedule (every
    per-segment bound is zero); the solve returns it directly.

    Parameters
    ----------
    leg:
        The leg configuration.
    seed:
        RNG seed for DE and the random multi-starts (reproducible).
    n_starts:
        Number of SLSQP polish starts (the zero schedule plus ``n_starts - 1``
        random restarts), in addition to the optional DE pass.
    use_de:
        Run a DE global pass before the SLSQP polishes (defence-in-depth).

    Returns
    -------
    LowThrustLegResult
        The best feasible schedule found (lowest total ΔV among feasible
        candidates; if none is feasible, the lowest-defect candidate).
    """
    n = leg.n_segments
    dim = 3 * n
    if leg.tmax_kn == 0.0:
        # No thrust capability: only the zero schedule is admissible.
        return _evaluate(leg, np.zeros(dim, dtype=np.float64))

    bounds = _dv_box(leg)
    candidates: list[NDArray[np.float64]] = []

    if use_de:

        def de_constraint(flat: NDArray[np.float64]) -> NDArray[np.float64]:
            # Stack the (signed) defect and the thrust slacks; DE demands each
            # within [lb, ub]. The defect rows are equalities (lb = ub = 0).
            return np.concatenate([_defect_equality(leg, flat), _thrust_slacks(leg, flat)])

        lb = np.concatenate([np.zeros(6), np.zeros(n)])
        ub = np.concatenate([np.zeros(6), np.full(n, np.inf)])
        nlc = NonlinearConstraint(de_constraint, lb=lb, ub=ub)
        de_any = cast(Any, differential_evolution)
        try:
            de_res: Any = de_any(
                _total_dv,
                bounds,
                args=(leg,),
                constraints=nlc,
                seed=seed,
                polish=False,
                maxiter=_DE_MAXITER,
                popsize=_DE_POPSIZE,
                tol=_DE_TOL,
            )
            candidates.append(np.asarray(de_res.x, dtype=np.float64))
        except (ValueError, RuntimeError):
            pass

    # SLSQP starts: the zero schedule first (the ballistic-closure anchor), then
    # seeded random restarts inside a small fraction of the capability box.
    candidates.append(np.zeros(dim, dtype=np.float64))
    rng = np.random.default_rng(seed)
    cap = float(np.max(segment_dv_bounds(leg, np.zeros((n, 3)))))
    for _ in range(max(0, n_starts - 1)):
        candidates.append(rng.normal(scale=0.05 * cap, size=dim))

    results: list[LowThrustLegResult] = []
    for x0 in candidates:
        x_final = _slsqp_polish(leg, x0, _total_dv, bounds)
        results.append(_evaluate(leg, x_final))

    return _select_best(results, key="total_dv")


def solve_leg_max_mass(
    leg: SimsFlanaganLeg,
    phase1: LowThrustLegResult,
    *,
    seed: int = 0,
) -> LowThrustLegResult:
    r"""Phase 2: re-optimise the Phase-1 schedule to maximise final mass.

    Local SLSQP solve (Yam Eq. 6) seeded from the Phase-1 schedule, minimising
    ``-final_mass`` under the same defect equality + thrust inequality
    constraints. Because final mass falls monotonically with total ΔV (rocket
    equation), maximising mass is closely allied to minimising ΔV; this phase
    refines the Phase-1 optimum onto the mass objective without losing
    feasibility.

    Parameters
    ----------
    leg:
        The leg configuration.
    phase1:
        The Phase-1 (:func:`solve_leg_min_dv`) result used as the seed.
    seed:
        RNG seed (kept for signature symmetry; this phase is a single local
        polish from the Phase-1 schedule).

    Returns
    -------
    LowThrustLegResult
        The mass-refined schedule if feasible and at least as good as Phase 1;
        otherwise the Phase-1 result is returned unchanged (Phase 2 never
        degrades the solution).
    """
    if leg.tmax_kn == 0.0:
        return phase1
    bounds = _dv_box(leg)
    x0 = phase1.dvs.reshape(-1).copy()
    x_final = _slsqp_polish(leg, x0, _neg_final_mass, bounds)
    candidate = _evaluate(leg, x_final)
    # Phase 2 must not degrade: keep it only if feasible and mass >= Phase 1.
    if candidate.converged and candidate.final_mass_kg >= phase1.final_mass_kg - 1.0e-9:
        return candidate
    return phase1


def _select_best(results: list[LowThrustLegResult], *, key: str) -> LowThrustLegResult:
    """Pick the best result: feasible ones first, then by ``key`` ascending.

    ``key`` is ``"total_dv"`` (Phase 1) — lower is better. Infeasible
    candidates sort last, ranked among themselves by total defect so the
    fallback is the least-bad geometry.
    """
    if not results:
        raise SimsFlanaganError("_select_best requires a non-empty result list")

    def sort_key(r: LowThrustLegResult) -> tuple[int, float, float]:
        feasible_rank = 0 if r.converged else 1
        primary = r.total_dv_kms if r.converged else (r.defect_pos_km + r.defect_vel_kms)
        return (feasible_rank, primary, r.total_dv_kms)

    return min(results, key=sort_key)


def chain_feasible(
    legs: list[SimsFlanaganLeg],
    schedules: list[NDArray[np.float64]],
) -> bool:
    """Whether every leg in a chain is independently feasible.

    A thin convenience over :func:`cyclerfinder.core.sims_flanagan.leg_feasible`
    for the multi-leg case; the inter-leg flyby bend coupling is checked
    separately via
    :func:`cyclerfinder.core.sims_flanagan.flyby_bend_slacks`.
    """
    if len(legs) != len(schedules):
        raise SimsFlanaganError(
            f"need one schedule per leg: {len(legs)} legs, {len(schedules)} schedules"
        )
    return all(leg_feasible(leg, dvs) for leg, dvs in zip(legs, schedules, strict=True))


__all__ = [
    "LowThrustLegResult",
    "chain_feasible",
    "solve_leg_max_mass",
    "solve_leg_min_dv",
]
