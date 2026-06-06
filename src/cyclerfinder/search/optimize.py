"""Per-cell inner-timing optimisation (M5).

Spec references
---------------
* §6 — the top-level :func:`find_cyclers` discovery pipeline contract.
* §9 — published rediscovery anchors (5.65 km/s Earth, 3.05 km/s Mars on
  the 2-synodic E-M-E cycler) and the V∞ > 11 km/s degenerate-solution
  guard.
* §12(a) — two optimisation modes: *idealized* (strict closure in the
  circular-coplanar model) and *ephemeris* (finite-horizon TCM on a
  real ephemeris). The ephemeris-mode signature is locked here; the
  body is filled in M6.
* §12(d) — hard inequality constraints (V∞ ≤ cap, r_p ≥ r_p_min, bend
  achievable) on the optimiser. Not soft regularisers.
* §13.4 — structured inner search: free-return / resonance construction
  fixes most timing parameters; the residual continuous DOF is covered
  by a fixed, reproducible multi-start grid plus a local polish (and
  optional global wrapper). This is the single most important paragraph
  in M5's design — blind DE on a guessed sequence is the documented
  failure mode (§10).
* §13.8 — the ``Cell`` carries the catalogue's reproducibility identity;
  the ``seed`` argument here is the other half.

Design summary
--------------
M5 owns the **continuous timing layer** within a single :class:`Cell`
from M4. Discrete cell structure (sequence, period multiplier, per-leg
revolutions and branches) is the M4 enumerator's job; M5 fixes those
and searches the remaining encounter-time degrees of freedom.

The optimisation pipeline for ``optimise_cell_idealized``:

1. Compute the target heliocentric period from the cell's discrete
   structure (``period_k * T_syn(pair)``).
2. Generate ``n_starts`` initial guesses from a deterministic
   multi-start grid anchored on the free-return seed (interior
   encounter times equispaced across the period).
3. Local-polish every start with SLSQP, honouring the per-encounter
   hard inequality constraints (``V∞ ≤ vinf_cap`` and ``r_p ≥ rp_min``).
4. Optional global wrapper (scipy ``differential_evolution``) as
   defence-in-depth, followed by another SLSQP polish.
5. Select the best feasible result across all starts; infeasibles sort
   last via a ``+inf`` composite. Build the final cycler, score it via
   M4's :func:`~cyclerfinder.model.score.score`, and return the
   :class:`OptimisationResult`.

``find_cyclers`` is the spec §6 top-level wrapper: enumerate cells from
M4 (Tisserand-pruned), optimise each via :func:`optimise_cell_idealized`,
filter by hard-constraints, sort by composite score, return the top N.

Plan: ``docs/phases/m5-optimisation/plan.md``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import NonlinearConstraint, differential_evolution, minimize

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import flyby_dv_for
from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError
from cyclerfinder.model import Cycler, Score
from cyclerfinder.model.score import composite_score, score
from cyclerfinder.search.construct import construct_cycler
from cyclerfinder.search.resonance import (
    beat_period_days,
    multi_body_beat_days,
    synodic_period_days,
)
from cyclerfinder.search.sequence import Cell, feasible_cells

# ---------------------------------------------------------------------------
# Module-level constants / tolerances
# ---------------------------------------------------------------------------

_PATHOLOGICAL_OBJECTIVE_KMS: float = 1.0e6
"""Penalty value returned by :func:`_objective` when the underlying
cycler construction fails (e.g. Lambert pathology). Large enough to be
unambiguously worse than any feasible objective; finite so SLSQP and DE
don't choke on ``inf`` / ``nan``."""

_CONSTRAINT_VIOLATION_SLACK_KMS: float = -1.0
"""Negative slack returned from each per-encounter constraint when
``_build_cycler_from_x`` returns ``None`` (cycler couldn't be built).
SLSQP's convention is ``fun(x) >= 0`` ⇒ satisfied; any negative value
trips the constraint."""

_SLSQP_FEASIBILITY_TOL_KMS: float = 1.0e-6
"""Numerical slack allowed when re-checking constraints at SLSQP's
returned ``x_final``. SLSQP can land marginally on the wrong side of a
constraint boundary; this tolerance distinguishes that from a genuine
infeasibility."""

_DE_MAXITER: int = 50
"""Scipy ``differential_evolution`` generation cap. Small but enough to
cover a 1-6 dim parameter box; the multi-start grid + SLSQP polish
carry most of the global coverage per spec §13.4."""

_DE_POPSIZE: int = 8
"""Scipy DE population multiplier. With this and the maxiter above, the
total budget is ~400 evaluations — bounded per cell."""

_DE_TOL: float = 1.0e-4
"""Scipy DE convergence tolerance (population stdev fraction)."""

_SLSQP_MAXITER: int = 200
"""Scipy SLSQP iteration cap. SLSQP typically converges in 30-80 on this
problem size."""

_SLSQP_FTOL: float = 1.0e-6
"""Scipy SLSQP objective tolerance."""

_BOUNDS_INSET_FRAC: float = 0.01
"""Fractional inset of the per-encounter parameter bounds inside
``(0, T)``. ``t_i ∈ (0.01 T, 0.99 T)`` keeps SLSQP / DE strictly inside
the period and prevents the optimiser from collapsing two encounters."""


# ---------------------------------------------------------------------------
# Dataclasses (plan §3.4)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _StartRecord:
    """Diagnostic record for one start in the multi-start grid.

    Module-private. Surfaced on
    :attr:`OptimisationResult.optimiser_history` for debug / inspection;
    not part of the public API contract beyond ``len()`` and indexing.

    Attributes
    ----------
    start_index:
        Position in the multi-start grid (0-based). ``-1`` flags the DE
        wrapper's polished result; ``-2, -3, …`` flag caller-supplied
        warm starts (#52) — all distinguishable from the grid starts.
    x0:
        Initial parameter vector (the interior encounter times in
        seconds, length ``N - 2`` for an N-encounter cell).
    x_final:
        Converged parameter vector after SLSQP polish.
    objective_value:
        Objective at ``x_final``, km/s. Closure residual + sum of
        flyby ΔV.
    constraints_satisfied:
        ``True`` iff all per-encounter inequality constraints are
        satisfied at ``x_final`` to within
        :data:`_SLSQP_FEASIBILITY_TOL_KMS`.
    nit:
        SLSQP iteration count at exit.
    success:
        SLSQP's own convergence flag (``OptimizeResult.success``).
    """

    start_index: int
    x0: tuple[float, ...]
    x_final: tuple[float, ...]
    objective_value: float
    constraints_satisfied: bool
    nit: int
    success: bool


@dataclass(frozen=True)
class OptimisationResult:
    """The output of one cell's per-cell optimisation.

    All fields are immutable (frozen). The tuple-valued
    ``optimiser_history`` is intentionally a ``tuple`` rather than a
    ``list`` so the frozen invariant is structural, not just by-convention.

    Spec references: §12(a) (two-mode optimisation produces a result of
    this shape), §12(d) (``constraints_satisfied`` records the
    hard-inequality outcome), §13.8 (``cell`` and the seed used populate
    the catalogue record's ``reproducibility`` block).

    Attributes
    ----------
    cell:
        The M4 :class:`Cell` that was optimised. Carried so downstream
        consumers (M7 ledger, M8 reporter) can populate ``cell_id``
        directly from ``result.cell.id``.
    best_cycler:
        The optimised :class:`Cycler` from M3. Consumers can use any
        method on it (``maintenance_dv``, ``radial_span``, etc.).
    best_score:
        The :class:`Score` from M4 for the optimised cycler. Reads
        ``taxi_cost_kms``, ``max_vinf_kms``, etc. via attribute access.
    closure_residual_kms:
        Rotating-frame closure residual at the optimised geometry, km/s.
        Computed with ``omega = 2π / target_period_sec`` so consumers
        don't recompute the transform.
    optimiser_history:
        Per-start diagnostic records. Length is ``n_starts`` (or
        ``n_starts + 1`` when DE is enabled and produced a polished
        result). Treat as opaque diagnostic data.
    converged:
        Global success indicator — ``True`` iff at least one start
        produced a feasible converged result.
    constraints_satisfied:
        ``True`` iff every hard inequality constraint is satisfied at
        the selected best result. The "trustworthy" predicate for
        downstream consumers is ``converged ∧ constraints_satisfied``;
        either alone is insufficient.
    """

    cell: Cell
    best_cycler: Cycler
    best_score: Score
    closure_residual_kms: float
    optimiser_history: tuple[_StartRecord, ...]
    converged: bool
    constraints_satisfied: bool


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _target_period_sec(cell: Cell) -> float:
    """Resolve the target heliocentric period for ``cell``, seconds.

    Dispatch (M8 plan §2):

    * ``cell.period_basis`` set — use the catalogue's *anchor pair*:
      ``T_syn(*period_basis) * period_k``. ``period_k`` is the sourced
      catalogue value, never rewritten, so a Cell stays traceable to its
      YAML row (EMEEVE: ``T_syn(E,M) * 3 ~ 6.41 yr``).
    * ``period_basis is None`` and ``len(bodies) >= 3`` — fall back to the
      body set's *natural beat* via
      :func:`~cyclerfinder.search.resonance.multi_body_beat_days` /
      :func:`~cyclerfinder.search.resonance.beat_period_days`
      (``3*T_syn(E,M) ~ 4*T_syn(E,V) ~ 6.40 yr`` for ``["V","E","M"]``).
    * ``period_basis is None`` and ``len(bodies) == 2`` — the M5 native
      single-pair formula, preserved byte-for-byte.
    """
    if cell.period_basis is not None:
        a, b = cell.period_basis
        t_syn_days = synodic_period_days(a, b)
        return t_syn_days * cell.period_k * SECONDS_PER_DAY

    n = len(cell.bodies)
    if n < 2:
        raise ValueError(
            f"cell.bodies must have at least 2 entries; got {cell.bodies!r}",
        )
    if n == 2:
        t_syn_days = synodic_period_days(cell.bodies[0], cell.bodies[1])
        return t_syn_days * cell.period_k * SECONDS_PER_DAY
    bodies = list(cell.bodies)
    tuples = multi_body_beat_days(bodies)
    if not tuples:
        raise ValueError(
            f"no integer beat commensurability found for bodies={bodies!r} "
            f"within the resonance.multi_body_beat_days default tolerance; "
            f"this body set has no natural cycler beat in the searched k range",
        )
    beat_days = beat_period_days(bodies, tuples[0])
    return beat_days * cell.period_k * SECONDS_PER_DAY


def _free_return_seed(cell: Cell, target_period_sec: float) -> tuple[float, ...]:
    """Spec §13.4 structural anchor: equispaced encounter epochs across ``T``.

    For an ``N``-encounter cell, returns ``(0, T/(N-1), 2T/(N-1), ...,
    T)`` — the simplest free-return-style guess that places interior
    encounters at rational fractions of the period. The multi-start
    grid (:func:`_multi_start_grid`) perturbs the interior epochs
    around this anchor.

    For the M5 gate's 2-syn E-M-E cell (``N = 3``) this returns
    ``(0, T/2, T)``.

    Notes
    -----
    The first and last entries are pinned (``0`` and ``T``) and not
    parameters of the optimiser — they're returned for clarity so the
    caller can pass a full encounter-time vector to
    :func:`construct_cycler`. The free parameters are
    ``seed[1:-1]`` (the interior times).
    """
    n = len(cell.sequence)
    if n < 2:
        raise ValueError(f"cell.sequence must have >= 2 entries; got {cell.sequence!r}")
    return tuple(i * target_period_sec / (n - 1) for i in range(n))


def _multi_start_grid(
    cell: Cell,
    target_period_sec: float,
    n_starts: int,
    seed: int,
) -> list[tuple[float, ...]]:
    """Deterministic Latin-Hypercube multi-start over the interior epochs.

    Spec §13.4: "cover the remaining continuous DOF with a fixed,
    reproducible multi-start grid plus local polish, so coverage within
    a cell is systematic, not stochastic."

    Strategy (task #53):

    - Start ``0`` is exactly the free-return seed
      (:func:`_free_return_seed`) — the structural anchor SLSQP polishes
      first.
    - Starts ``1 …`` are a seeded Latin-Hypercube Sample (LHS) of the
      interior-epoch box ``[inset·T, (1-inset)·T]^{N-2}``. Each of the
      ``n_starts - 1`` samples sits in its own stratum along every
      dimension (one point per row and per column of the LHS grid), so
      coverage is stratified and the samples are mutually distinct — a
      strict improvement over the prior fixed ``±k·T`` perturbation
      table, which collapsed to ``≤ 4`` unique positive magnitudes for
      ``n_interior = 1`` cells (duplicate starts, zero minimum pairwise
      separation).

    Determinism: the LHS uses ``np.random.default_rng(seed)`` for both the
    in-stratum jitter and the per-dimension permutations, so the same
    ``seed`` always yields bitwise-identical start vectors (the optimiser's
    reproducibility tests rely on this).

    Returns
    -------
    list[tuple[float, ...]]
        Length ``n_starts``. Each entry is a tuple of length ``N - 2``
        (the interior epochs only — the pinned endpoints ``0`` and ``T``
        are added back by :func:`_build_cycler_from_x`).

    Notes
    -----
    For ``N = 2`` (a degenerate 1-leg "cell") the parameter vector is
    empty and every start is the same empty tuple. For ``N = 3`` (the
    M5 gate) there is one interior epoch and the LHS stratifies it across
    the interior bounds; each interior vector is sorted to keep the
    encounter epochs monotone (an unsorted start would make
    :func:`construct_cycler` raise immediately).
    """
    n = len(cell.sequence)
    if n < 2:
        raise ValueError(f"cell.sequence must have >= 2 entries; got {cell.sequence!r}")
    n_interior = n - 2
    free_return = _free_return_seed(cell, target_period_sec)
    interior_anchor = free_return[1:-1]

    starts: list[tuple[float, ...]] = []
    # Start 0: the free-return seed exactly.
    starts.append(tuple(interior_anchor))
    if n_starts <= 1 or n_interior == 0:
        return starts[:n_starts] if n_starts > 0 else []

    n_samples = n_starts - 1
    lo = _BOUNDS_INSET_FRAC * target_period_sec
    hi = (1.0 - _BOUNDS_INSET_FRAC) * target_period_sec
    span = hi - lo

    rng = np.random.default_rng(seed)
    # Standard LHS: for each dimension, place one jittered sample in each of
    # the n_samples equal strata, then permute the stratum->sample assignment
    # independently per dimension so the samples are not axis-aligned.
    # ``lhs_unit`` is (n_samples, n_interior) in the unit hypercube.
    lhs_unit = np.empty((n_samples, n_interior), dtype=np.float64)
    for d in range(n_interior):
        jitter = rng.random(n_samples)
        strata = (np.arange(n_samples) + jitter) / n_samples
        perm = rng.permutation(n_samples)
        lhs_unit[:, d] = strata[perm]

    for s in range(n_samples):
        sample = tuple(float(lo + span * lhs_unit[s, d]) for d in range(n_interior))
        # Keep interior epochs monotone (weak constraint; SLSQP refines).
        # Clip is defensive — LHS already lies inside [lo, hi].
        sample_clipped = tuple(_clip_interior(t, target_period_sec) for t in sorted(sample))
        starts.append(sample_clipped)
    return starts


def _clip_interior(t: float, target_period_sec: float) -> float:
    """Clip ``t`` to the strict interior ``(BOUNDS_INSET_FRAC * T,
    (1 - BOUNDS_INSET_FRAC) * T)``."""
    lo = _BOUNDS_INSET_FRAC * target_period_sec
    hi = (1.0 - _BOUNDS_INSET_FRAC) * target_period_sec
    return max(lo, min(hi, t))


def interior_epochs_from_leg_tofs(
    leg_tofs_days: Sequence[float],
    target_period_sec: float,
) -> tuple[float, ...]:
    """Map per-leg times-of-flight to clipped interior encounter epochs (seconds).

    Absolute-cumulative mapping (#52): encounter ``j`` sits at the
    cumulative sum of the first ``j`` leg ToFs. Encounter 0 is pinned at
    ``t = 0`` and the final encounter at ``t = T``, so only the interior
    cumulative sums are returned — ``cumsum(leg_tofs_days[:-1])`` — giving
    a vector of length ``len(leg_tofs_days) - 1`` (i.e. ``N - 2`` interior
    epochs for an ``N``-encounter cycle). Each epoch is clipped to the
    optimiser's strict interior via :func:`_clip_interior` and the result
    is sorted, so it satisfies the monotonic-epoch precondition that
    :func:`construct_cycler` enforces.

    This is the warm-start shape :func:`optimise_cell_idealized` expects:
    pass the result as one element of its ``warm_starts`` list. For the
    Aldrin classic cell (``legs = [146, 634]`` days) this yields a single
    interior epoch at ~146 days — the asymmetric Earth→Mars transit the
    equispaced free-return seed (``T/2``) misses.
    """
    cum = 0.0
    interior: list[float] = []
    for tof_days in leg_tofs_days[:-1]:  # drop final leg: last encounter pinned at T
        cum += float(tof_days) * SECONDS_PER_DAY
        interior.append(_clip_interior(cum, target_period_sec))
    return tuple(sorted(interior))


def _full_times_from_x(
    x: NDArray[np.float64],
    target_period_sec: float,
) -> list[float]:
    """Reconstruct the full encounter-time list from the free vector ``x``.

    The full list is ``[0.0, *x, target_period_sec]``. The free vector
    holds only the interior epochs; the endpoints are pinned by
    convention (encounter 0 at ``t=0``, encounter N-1 at ``t=T``).
    """
    return [0.0, *(float(xi) for xi in x), float(target_period_sec)]


def _build_cycler_from_x(
    x: NDArray[np.float64],
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
) -> Cycler | None:
    """Wrap :func:`construct_cycler` for an optimiser's parameter vector.

    On Lambert pathologies (non-monotonic times, ``ValueError`` from
    :func:`construct_cycler`), returns ``None`` so the optimiser's
    objective / constraint callbacks can convert into a finite penalty
    rather than propagating the exception into scipy's internals.

    Parameters
    ----------
    x:
        Length ``N - 2`` parameter vector (the interior encounter
        epochs in seconds).
    cell:
        Discrete structural specification of the cycler — supplies the
        sequence, per-leg revs, and per-leg branches.
    ephem:
        Planet-state provider.
    target_period_sec:
        Total cycler period; pins the last encounter to ``t = T``.

    Returns
    -------
    Cycler | None
        The constructed cycler on success, ``None`` on any
        :class:`ValueError` from :func:`construct_cycler`.
    """
    times = _full_times_from_x(x, target_period_sec)
    try:
        return construct_cycler(
            sequence=list(cell.sequence),
            encounter_times_sec=times,
            ephem=ephem,
            max_revs_per_leg=list(cell.per_leg_revs),
            branch_per_leg=list(cell.per_leg_branch),
        )
    except (ValueError, LambertConvergenceError, LambertGeometryError):
        # Lambert Newton can fail to converge in degenerate geometries the
        # optimiser may probe (very long ToF, near-180° transfer, large
        # multi-rev z values). Treat as a pathological point; the caller
        # converts ``None`` into ``_PATHOLOGICAL_OBJECTIVE_KMS`` so DE/SLSQP
        # see a finite penalty rather than an exception.
        return None


def _r_p_required(
    vin_vec: NDArray[np.float64],
    vout_vec: NDArray[np.float64],
    mu_planet: float,
) -> float:
    """Periapsis radius (km) needed for a single-pass ballistic bend
    from ``vin_vec`` to ``vout_vec`` at a planet of gravitational
    parameter ``mu_planet``.

    Derived from the bend formula
    ``sin(delta/2) = 1 / (1 + r_p * V_inf^2 / mu)``:

        ``r_p = (1 / sin(delta/2) - 1) * mu / V_inf^2``

    where ``delta`` is the angle between the V∞ vectors and ``V_inf``
    is their mean magnitude. The result is the *minimum* periapsis that
    can deliver the requested bend; if it exceeds the body's safe
    periapsis floor (:data:`SAFE_PERIHELION_KM`), the flyby is feasible.

    Returns
    -------
    float
        Periapsis radius in km. Returns a small positive number
        (effectively zero) when the bend exceeds ~180° (the formula
        diverges); the constraint then rejects the encounter as
        infeasible.
    """
    vin_mag = float(np.linalg.norm(vin_vec))
    vout_mag = float(np.linalg.norm(vout_vec))
    if vin_mag <= 0.0 or vout_mag <= 0.0:
        return 0.0
    v_mean = 0.5 * (vin_mag + vout_mag)
    cos_arg = float(np.dot(vin_vec, vout_vec)) / (vin_mag * vout_mag)
    cos_arg = max(-1.0, min(1.0, cos_arg))
    delta = math.acos(cos_arg)
    # sin(delta/2) in (0, 1] for delta in (0, pi]; near delta == 0 the
    # required r_p is huge (no bend needed → any flyby works); near
    # delta == pi the required r_p is tiny / negative (no flyby can
    # turn the velocity by 180°). Clamp the sine to a small floor so
    # the division stays finite.
    half = 0.5 * delta
    sin_half = math.sin(half)
    if sin_half < 1.0e-12:
        # Trivially feasible: delta ~ 0; return a very large r_p so the
        # constraint always passes.
        return 1.0e18
    return float((1.0 / sin_half - 1.0) * mu_planet / (v_mean * v_mean))


def _objective(
    x: NDArray[np.float64],
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    omega_rad_per_s: float,
) -> float:
    """Optimisation objective: closure residual + sum of flyby ΔV, km/s.

    Both terms are non-negative. On Lambert / construct failure returns
    :data:`_PATHOLOGICAL_OBJECTIVE_KMS` so the SLSQP / DE drivers don't
    propagate exceptions through scipy's internals.

    Parameters
    ----------
    x:
        Length ``N - 2`` parameter vector (interior encounter epochs).
    cell, ephem, target_period_sec:
        Forwarded to :func:`_build_cycler_from_x`.
    omega_rad_per_s:
        Rotating-frame angular rate used by
        :meth:`Cycler.closure_residual`. Typically
        ``2π / target_period_sec`` so the closure compares one period
        of the synodic frame.

    Returns
    -------
    float
        Composite km/s scalar suitable for direct minimisation.
    """
    cyc = _build_cycler_from_x(x, cell, ephem, target_period_sec)
    if cyc is None:
        return _PATHOLOGICAL_OBJECTIVE_KMS
    try:
        residual = cyc.closure_residual(omega_rad_per_s=omega_rad_per_s)
        dv = 0.0
        for enc in cyc.encounters:
            dv += flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out)
    except (ValueError, KeyError):
        return _PATHOLOGICAL_OBJECTIVE_KMS
    if not math.isfinite(residual) or not math.isfinite(dv):
        return _PATHOLOGICAL_OBJECTIVE_KMS
    return float(residual + dv)


def _per_encounter_slacks(
    x: NDArray[np.float64],
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    vinf_cap: float,
    rp_factors: dict[str, float],
) -> tuple[list[float], list[float]]:
    """Return ``(vinf_slacks, rp_slacks)`` per encounter.

    Each list has length ``len(cell.sequence)``. Positive slack ⇒
    constraint satisfied; negative ⇒ violated. On cycler-build failure
    returns ``[-1.0] * N`` for both lists (uniformly infeasible).

    Internal helper shared by :func:`_constraints` and the SLSQP
    post-hoc feasibility check.
    """
    n = len(cell.sequence)
    cyc = _build_cycler_from_x(x, cell, ephem, target_period_sec)
    if cyc is None:
        return [_CONSTRAINT_VIOLATION_SLACK_KMS] * n, [_CONSTRAINT_VIOLATION_SLACK_KMS] * n

    vinf_slacks: list[float] = []
    rp_slacks: list[float] = []
    for enc in cyc.encounters:
        vin_mag = float(np.linalg.norm(enc.vinf_in))
        vout_mag = float(np.linalg.norm(enc.vinf_out))
        vinf_slacks.append(vinf_cap - max(vin_mag, vout_mag))
        try:
            mu_planet = PLANETS[enc.body].mu_km3_s2
            rp_floor = rp_factors.get(enc.body, 1.0) * SAFE_PERIHELION_KM[enc.body]
        except KeyError:
            rp_slacks.append(_CONSTRAINT_VIOLATION_SLACK_KMS)
            continue
        rp_required = _r_p_required(enc.vinf_in, enc.vinf_out, mu_planet)
        rp_slacks.append(rp_required - rp_floor)
    return vinf_slacks, rp_slacks


def _constraints(
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    vinf_cap: float,
    rp_factors: dict[str, float],
) -> list[dict[str, Any]]:
    """Return a list of SLSQP-style ``{type: 'ineq', fun: callable}`` dicts.

    Two constraints per encounter:

    1. ``vinf_cap - max(||vinf_in||, ||vinf_out||) >= 0`` —
       per-encounter V∞ cap.
    2. ``r_p_required - rp_floor >= 0`` — per-encounter r_p floor.

    Per spec §12(d) these are **hard inequalities**, not soft
    regularisers; the objective :func:`_objective` carries no V∞ or
    r_p penalty terms. Per-encounter formulation (rather than aggregate
    max) keeps the constraints smooth.

    SLSQP convention: ``fun(x) >= 0`` ⇒ constraint satisfied.
    """
    n = len(cell.sequence)

    def make_vinf_fun(i: int) -> Callable[[NDArray[np.float64]], float]:
        def fun(x: NDArray[np.float64]) -> float:
            vinf_slacks, _ = _per_encounter_slacks(
                x,
                cell,
                ephem,
                target_period_sec,
                vinf_cap,
                rp_factors,
            )
            return vinf_slacks[i]

        return fun

    def make_rp_fun(i: int) -> Callable[[NDArray[np.float64]], float]:
        def fun(x: NDArray[np.float64]) -> float:
            _, rp_slacks = _per_encounter_slacks(
                x,
                cell,
                ephem,
                target_period_sec,
                vinf_cap,
                rp_factors,
            )
            return rp_slacks[i]

        return fun

    constraints: list[dict[str, Any]] = []
    for i in range(n):
        constraints.append({"type": "ineq", "fun": make_vinf_fun(i)})
        constraints.append({"type": "ineq", "fun": make_rp_fun(i)})
    return constraints


def _check_constraints_satisfied(
    x: NDArray[np.float64],
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    vinf_cap: float,
    rp_factors: dict[str, float],
) -> bool:
    """True iff every per-encounter constraint is satisfied at ``x``
    to within :data:`_SLSQP_FEASIBILITY_TOL_KMS`."""
    vinf_slacks, rp_slacks = _per_encounter_slacks(
        x,
        cell,
        ephem,
        target_period_sec,
        vinf_cap,
        rp_factors,
    )
    tol = _SLSQP_FEASIBILITY_TOL_KMS
    return all(s >= -tol for s in vinf_slacks) and all(s >= -tol for s in rp_slacks)


def _bounds_for(cell: Cell, target_period_sec: float) -> list[tuple[float, float]]:
    """Per-interior-epoch ``(lo, hi)`` bounds for SLSQP / DE.

    Each interior epoch is bounded to ``(BOUNDS_INSET_FRAC * T,
    (1 - BOUNDS_INSET_FRAC) * T)``. Length ``N - 2``.
    """
    n_interior = len(cell.sequence) - 2
    lo = _BOUNDS_INSET_FRAC * target_period_sec
    hi = (1.0 - _BOUNDS_INSET_FRAC) * target_period_sec
    return [(lo, hi) for _ in range(n_interior)]


def _polish(
    x0: tuple[float, ...],
    start_index: int,
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    omega_rad_per_s: float,
    vinf_cap: float,
    rp_factors: dict[str, float],
) -> _StartRecord:
    """Single SLSQP local solve from ``x0`` honouring hard constraints.

    Returns a :class:`_StartRecord` whose ``constraints_satisfied`` is
    re-evaluated at the returned ``x_final`` (SLSQP can land marginally
    outside the feasible region).

    Degenerate case (``N == 2``, no free DOF): the parameter vector is
    empty; SLSQP can't be called. We evaluate the objective and
    constraints at the empty vector and return the record directly.
    """
    n_interior = len(cell.sequence) - 2
    if n_interior == 0:
        # No free parameter; evaluate objective and constraints at the
        # empty vector and return a synthetic record.
        x_arr = np.asarray([], dtype=np.float64)
        obj = _objective(x_arr, cell, ephem, target_period_sec, omega_rad_per_s)
        feasible = _check_constraints_satisfied(
            x_arr,
            cell,
            ephem,
            target_period_sec,
            vinf_cap,
            rp_factors,
        )
        return _StartRecord(
            start_index=start_index,
            x0=(),
            x_final=(),
            objective_value=obj,
            constraints_satisfied=feasible,
            nit=0,
            success=True,
        )

    bounds = _bounds_for(cell, target_period_sec)
    constraints = _constraints(cell, ephem, target_period_sec, vinf_cap, rp_factors)
    x0_arr = np.asarray(x0, dtype=np.float64)
    # Clip start into bounds defensively (SLSQP misbehaves on starts
    # outside bounds with method='SLSQP' under certain scipy versions).
    for i, (lo, hi) in enumerate(bounds):
        x0_arr[i] = max(lo, min(hi, float(x0_arr[i])))

    args_tuple: tuple[object, ...] = (
        cell,
        ephem,
        target_period_sec,
        omega_rad_per_s,
    )
    # scipy's `minimize` overloads in the type stubs don't admit the
    # exact (ndarray, Cell, Ephemeris, float, float) -> float signature
    # of `_objective`, nor the SLSQP constraint-list-of-dicts shape under
    # strict mode. Both are documented scipy runtime conventions; cast
    # the entry point itself to keep runtime behaviour intact.
    minimize_any = cast(Any, minimize)
    result_any: Any = minimize_any(
        _objective,
        x0_arr,
        args=args_tuple,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": _SLSQP_MAXITER, "ftol": _SLSQP_FTOL},
    )
    # scipy.optimize.OptimizeResult is duck-typed; mypy sees it as Any.
    x_final_arr = cast(NDArray[np.float64], np.asarray(result_any.x, dtype=np.float64))
    x_final = tuple(float(xi) for xi in x_final_arr)
    obj_final = float(result_any.fun)
    feasible = _check_constraints_satisfied(
        x_final_arr,
        cell,
        ephem,
        target_period_sec,
        vinf_cap,
        rp_factors,
    )
    nit = int(result_any.nit) if hasattr(result_any, "nit") else 0
    success = bool(result_any.success)
    return _StartRecord(
        start_index=start_index,
        x0=tuple(float(xi) for xi in x0_arr),
        x_final=x_final,
        objective_value=obj_final,
        constraints_satisfied=feasible,
        nit=nit,
        success=success,
    )


def _de_pass(
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
    omega_rad_per_s: float,
    vinf_cap: float,
    rp_factors: dict[str, float],
    seed: int,
) -> tuple[float, ...] | None:
    """One scipy ``differential_evolution`` global pass.

    Returns the DE solution as a parameter tuple, or ``None`` if the
    parameter space is empty (``N == 2``, no free DOF) or the DE call
    fails. The caller (typically :func:`optimise_cell_idealized`) then
    feeds the result into :func:`_polish` for a final SLSQP polish.

    DE budget is bounded by :data:`_DE_MAXITER` and :data:`_DE_POPSIZE`;
    the multi-start grid + SLSQP carries most of the global coverage
    per spec §13.4, so DE is the defence-in-depth layer rather than the
    primary mechanism.
    """
    n_interior = len(cell.sequence) - 2
    if n_interior == 0:
        return None
    bounds = _bounds_for(cell, target_period_sec)

    # Nonlinear constraint for DE. NonlinearConstraint wants a callable
    # that returns an array of constraint function values, with
    # element-wise lower/upper bounds. We bundle all ``2 * N`` slacks
    # into one vector and demand each ``>= 0``.
    def constraint_fun(x: NDArray[np.float64]) -> NDArray[np.float64]:
        vinf_slacks, rp_slacks = _per_encounter_slacks(
            x,
            cell,
            ephem,
            target_period_sec,
            vinf_cap,
            rp_factors,
        )
        return np.asarray(vinf_slacks + rp_slacks, dtype=np.float64)

    n_slacks = 2 * len(cell.sequence)
    nlc = NonlinearConstraint(
        constraint_fun,
        lb=np.zeros(n_slacks),
        ub=np.full(n_slacks, np.inf),
    )

    args_tuple: tuple[object, ...] = (
        cell,
        ephem,
        target_period_sec,
        omega_rad_per_s,
    )
    # scipy's `differential_evolution` overloads in the type stubs are
    # likewise tighter than the runtime; cast through Any keeps runtime
    # behaviour intact while letting strict mypy pass.
    de_any = cast(Any, differential_evolution)
    try:
        result_any: Any = de_any(
            _objective,
            bounds,
            args=args_tuple,
            constraints=nlc,
            seed=seed,
            polish=False,
            maxiter=_DE_MAXITER,
            popsize=_DE_POPSIZE,
            tol=_DE_TOL,
            # Serial DE: workers=-1 was tested but the args_tuple (Cell,
            # Ephemeris) pickling overhead per evaluation dominated the
            # parallel speedup. The right parallelism point is the outer
            # 5x multi-start grid, not inside DE — future enhancement.
        )
    except (ValueError, RuntimeError):
        return None
    x_de = cast(NDArray[np.float64], np.asarray(result_any.x, dtype=np.float64))
    return tuple(float(xi) for xi in x_de)


def _composite_with_constraints(record: _StartRecord) -> float:
    """Sortable composite: infeasibles sort last via ``+inf``."""
    if not record.constraints_satisfied:
        return math.inf
    return record.objective_value


def _select_best(records: list[_StartRecord]) -> _StartRecord:
    """Pick the best record across all starts.

    Sort key: ``_composite_with_constraints`` ascending; ties broken by
    ``start_index`` so the result is deterministic across runs.
    """
    if not records:
        raise ValueError("_select_best requires a non-empty list")

    def key(r: _StartRecord) -> tuple[float, int]:
        return (_composite_with_constraints(r), r.start_index)

    return min(records, key=key)


def _sentinel_cycler(
    cell: Cell,
    ephem: Ephemeris,
    target_period_sec: float,
) -> Cycler:
    """Build a fallback cycler from the free-return seed.

    Used when the multi-start polish *and* the optional DE pass both
    fail to produce a buildable cycler. The cycler returned here is
    not feasible (its closure residual is whatever the free-return
    geometry produces); it exists only so :class:`OptimisationResult`
    carries a valid object reference for downstream consumers, with
    ``converged=False`` / ``constraints_satisfied=False`` flagging the
    failure.

    If even the free-return geometry can't be built (sequence pathology
    at the cell level), this function will raise — that's a genuine
    upstream bug (M4 yielded a structurally invalid cell) and should
    not be swallowed.
    """
    free_return = _free_return_seed(cell, target_period_sec)
    interior_arr: NDArray[np.float64] = np.asarray(free_return[1:-1], dtype=np.float64)
    cyc = _build_cycler_from_x(interior_arr, cell, ephem, target_period_sec)
    if cyc is None:
        raise RuntimeError(
            f"could not build sentinel cycler from free-return seed for "
            f"cell {cell.id!r}; the cell's structural specification appears "
            f"unbuildable. This is upstream of the optimisation layer.",
        )
    return cyc


# ---------------------------------------------------------------------------
# Public API: optimisation
# ---------------------------------------------------------------------------


def optimise_cell_idealized(
    cell: Cell,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    n_starts: int = 5,
    seed: int = 0,
    use_de: bool = True,
    rp_factors: dict[str, float] | None = None,
    target_period_sec: float | None = None,
    warm_starts: Sequence[Sequence[float]] | None = None,
) -> OptimisationResult:
    """Optimise encounter timings within ``cell`` (spec §12(a) idealized).

    Performs a multi-start SLSQP search over the free interior
    encounter epochs (spec §13.4 structured inner search), optionally
    wrapped by a scipy ``differential_evolution`` global pass for
    defence-in-depth. The objective is ``closure_residual +
    Σ flyby_dv`` in km/s, both non-negative; the constraints are
    per-encounter hard inequalities (``V∞ ≤ vinf_cap`` and
    ``r_p ≥ rp_min``) per spec §12(d).

    Parameters
    ----------
    cell:
        Discrete structural cell from M4. The optimiser does **not**
        change ``cell.sequence``, ``cell.period_k``, ``cell.per_leg_revs``,
        or ``cell.per_leg_branch`` — only the continuous timing within.
    ephem:
        Planet-state provider. M5 closes in
        ``Ephemeris(model="circular")``; passing an astropy backend
        will still run but the closure semantics are then approximate
        (spec §12(c)).
    vinf_cap:
        Hard ceiling on V∞ at every encounter, km/s.
    n_starts:
        Multi-start grid size. ``5`` is the spec §13.4 default; for the
        higher-dimensional VEM problem (M8) this may need to grow.
    seed:
        Reproducibility seed. Two runs with the same seed and same
        cell / ephem produce bitwise-identical results.
    use_de:
        If ``True`` (default), run a scipy ``differential_evolution``
        global pass after the multi-start polish and merge its
        polished result into the candidate set. Disable for the M8
        bulk-enumeration case where the per-cell budget is tight.
    rp_factors:
        Per-body multipliers on
        :data:`~cyclerfinder.core.constants.SAFE_PERIHELION_KM` for the
        r_p floor constraint. Default ``1.0`` everywhere.
    target_period_sec:
        Total cycler period, seconds. ``None`` (default) ⇒ derived from
        ``cell`` via :func:`_target_period_sec`.
    warm_starts:
        Optional caller-supplied interior-epoch seeds (#52). Each entry
        is a vector of ``N - 2`` interior encounter epochs in seconds —
        the same shape :func:`_multi_start_grid` produces and
        :func:`interior_epochs_from_leg_tofs` returns. Warm starts are
        polished through the same SLSQP path as the grid, *before* it,
        and ranked into the same candidate set; they let a caller seed
        the search from a known geometry (e.g. a catalogue entry's leg
        ToFs) whose interior timing the equispaced free-return seed
        misses. Each is defensively clipped to the strict interior and
        sorted. A warm start of the wrong length raises ``ValueError``.
        ``None`` (default) leaves the start set bitwise-unchanged.

    Returns
    -------
    OptimisationResult
        Frozen record bundling the best cycler, its score, the closure
        residual, and the per-start diagnostic history.

    Notes
    -----
    The "trustworthy" predicate for downstream consumers is
    ``result.converged and result.constraints_satisfied``; either alone
    is insufficient. SLSQP can converge to a local minimum outside the
    feasible region (the spec §10 degenerate-solution basin), in which
    case ``constraints_satisfied=False`` is the correct rejection
    signal even though ``converged=True``.
    """
    if rp_factors is None:
        rp_factors = {}
    if target_period_sec is None:
        target_period_sec = _target_period_sec(cell)
    omega_rad_per_s = 2.0 * math.pi / target_period_sec

    starts = _multi_start_grid(cell, target_period_sec, n_starts, seed)
    records: list[_StartRecord] = []

    # Caller-supplied warm starts (#52) are polished first. They share
    # the SLSQP path and candidate ranking with the grid; their negative
    # idx (DE uses -1) marks them in the diagnostic history.
    n_interior = len(cell.sequence) - 2
    for i, ws in enumerate(warm_starts or ()):
        if len(ws) != n_interior:
            raise ValueError(
                f"warm start {i} has {len(ws)} interior epochs; "
                f"cell {cell.id!r} expects {n_interior}",
            )
        x0_warm = tuple(float(_clip_interior(t, target_period_sec)) for t in sorted(ws))
        records.append(
            _polish(
                x0_warm,
                -(2 + i),
                cell,
                ephem,
                target_period_sec,
                omega_rad_per_s,
                vinf_cap,
                rp_factors,
            ),
        )

    for idx, x0 in enumerate(starts):
        record = _polish(
            x0,
            idx,
            cell,
            ephem,
            target_period_sec,
            omega_rad_per_s,
            vinf_cap,
            rp_factors,
        )
        records.append(record)

    if use_de:
        x_de = _de_pass(
            cell,
            ephem,
            target_period_sec,
            omega_rad_per_s,
            vinf_cap,
            rp_factors,
            seed,
        )
        if x_de is not None:
            de_record = _polish(
                x_de,
                -1,
                cell,
                ephem,
                target_period_sec,
                omega_rad_per_s,
                vinf_cap,
                rp_factors,
            )
            records.append(de_record)

    best = _select_best(records)

    best_x_arr: NDArray[np.float64] = np.asarray(best.x_final, dtype=np.float64)
    final_cyc = _build_cycler_from_x(best_x_arr, cell, ephem, target_period_sec)
    if final_cyc is None:
        # Fall back to the seed geometry; flag as failed.
        final_cyc = _sentinel_cycler(cell, ephem, target_period_sec)
        converged = False
        constraints_satisfied = False
        residual = float("inf")
    else:
        residual = final_cyc.closure_residual(omega_rad_per_s=omega_rad_per_s)
        converged = best.success
        constraints_satisfied = best.constraints_satisfied

    final_score = score(
        final_cyc,
        ephem,
        vinf_cap=vinf_cap,
        target_period_sec=target_period_sec,
        rp_factors=rp_factors,
    )

    return OptimisationResult(
        cell=cell,
        best_cycler=final_cyc,
        best_score=final_score,
        closure_residual_kms=residual,
        optimiser_history=tuple(records),
        converged=converged,
        constraints_satisfied=constraints_satisfied,
    )


def _multirev_min_tof_days(body_a: str, body_b: str, n_revs: int) -> float:
    """Conservative lower-bound ToF (days) for an ``n_revs``-revolution leg
    between ``body_a`` and ``body_b``.

    A multi-revolution Lambert leg is infeasible below its physical minimum
    time-of-flight. Rather than solving for the exact minimum-energy
    ``t_min(n_revs)`` (which depends on the actual departure / arrival
    positions and is unavailable at the seed-construction stage), this uses the
    Hohmann transfer between the two bodies' circular orbits as a conservative
    proxy: the minimum-energy semi-major axis is the mean of the two orbital
    radii, so a single direct half-revolution takes ``pi * sqrt(a^3 / mu)`` and
    ``n`` full revolutions add ``2*pi*n`` more of the transfer period. The
    ``n_revs`` leg must therefore last at least one such transfer-period span:

        t_min(n) ~= (2*n + 1) * pi * sqrt(a_h^3 / mu_sun)

    where ``a_h = (r_a + r_b) / 2``. ``n_revs == 0`` returns ``0.0`` (no floor;
    direct legs keep the historic ``0.1 * share`` lower bound).
    """
    if n_revs <= 0:
        return 0.0
    r_a = PLANETS[body_a].sma_au * AU_KM
    r_b = PLANETS[body_b].sma_au * AU_KM
    a_h = 0.5 * (r_a + r_b)
    t_half = math.pi * math.sqrt(a_h**3 / MU_SUN_KM3_S2)
    return (2 * n_revs + 1) * t_half / SECONDS_PER_DAY


def _ephemeris_tof_seed_and_bounds(
    cell: Cell, target_period_sec: float
) -> tuple[list[float], list[tuple[float, float]]]:
    """Equispaced per-leg ToF seed (days) and per-leg bounds for the
    real-ephemeris optimiser, derived from the cell's period.

    Mirrors the interior-epoch bounds logic of the idealised optimiser
    (:func:`_free_return_seed` / :func:`_multi_start_grid`): each leg
    starts at ``T/(N-1)`` and may range over ``[0.1 * share, 0.9 * T]``
    so the optimiser can redistribute time between legs while keeping
    every ToF strictly positive (landmine #2: ``construct_cycler``
    requires strictly increasing epochs).

    Multi-revolution legs (``cell.per_leg_revs[j] >= 1``) raise their lower
    bound to a Hohmann-proxy minimum ToF floor (see
    :func:`_multirev_min_tof_days`) so the optimiser is not seeded below the
    leg's physical ``t_min`` — a seed there would trip ``LambertConvergenceError``
    inside the objective on every trial. The seed itself is lifted onto the
    floor when the equispaced share falls below it.
    """
    n_legs = len(cell.sequence) - 1
    if n_legs < 1:
        raise ValueError(f"cell.sequence must have >= 2 entries; got {cell.sequence!r}")
    period_days = target_period_sec / SECONDS_PER_DAY
    share = period_days / n_legs
    hi = 0.9 * period_days  # a single leg may absorb most of the period
    base_lo = 0.1 * share
    seed: list[float] = []
    bounds: list[tuple[float, float]] = []
    for j in range(n_legs):
        revs = cell.per_leg_revs[j] if j < len(cell.per_leg_revs) else 0
        floor = _multirev_min_tof_days(cell.sequence[j], cell.sequence[j + 1], revs)
        lo = max(base_lo, floor)
        # Keep the lower bound below the upper one even for aggressive floors.
        lo = min(lo, 0.95 * hi)
        seed.append(max(share, lo))
        bounds.append((lo, hi))
    return seed, bounds


def _resolve_t0_multi_seed(
    cell: Cell,
    seed_days: list[float],
    priority_date: datetime,
    ephem: Ephemeris,
    vinf_targets_kms: dict[str, float],
    target_period_sec: float,
    n_candidates: int = 5,
) -> float | None:
    """Resolve a real launch epoch via the STAGE 3 multi-seed pool.

    Builds the primary :class:`PhaseSignature` from ``seed_days``, fans it into
    asymmetric leg-duration perturbation seeds (so a strongly asymmetric family
    phase-matches its own basin rather than a symmetric degenerate one), and
    returns the lowest-V_inf-mismatch window's J2000-relative seconds, or
    ``None`` if no window beats the cap.

    Provenance: the epoch is COMPUTED; ``vinf_targets_kms`` carries the sourced
    V_inf anchors used only as match targets.

    Parameters
    ----------
    cell:
        The cell being optimised; its ``sequence`` defines the body chain and
        per-encounter V_inf target lookup.
    seed_days:
        Primary per-leg ToF seed in days (length ``len(cell.sequence) - 1``).
    priority_date:
        Centre of the ±10 yr launch-window search.
    ephem:
        Ephemeris backend (``"astropy"`` for real dates).
    vinf_targets_kms:
        Per-body V_inf targets (km/s); must cover every body in
        ``cell.sequence``.
    target_period_sec:
        The cell's period in seconds; the conservation reference for the
        leg-duration perturbations.
    n_candidates:
        Per-seed window count requested before merging/ranking.
    """
    from datetime import UTC, timedelta

    from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
    from cyclerfinder.search.phase_match import (
        PhaseSignature,
        _dt_to_t_sec,
        find_candidate_windows,
        leg_duration_seeds,
    )

    try:
        vinf_target_per_enc = tuple(float(vinf_targets_kms[body]) for body in cell.sequence)
    except KeyError as exc:
        raise ValueError(
            f"vinf_targets_kms is missing body {exc.args[0]!r} required by "
            f"cell.sequence {cell.sequence!r}",
        ) from exc

    # Thread the cell's per-leg revs/branch into the phase-match Lambert so
    # same-body multi-rev legs (e.g. S1L1's Earth->Earth L1 loop) are scored on
    # the correct revolution/branch rather than mis-scored single-rev. Stay
    # empty (single-rev default) when every leg is 0-rev — byte-identical to the
    # pre-multi-rev epoch resolver for Aldrin and other direct-leg cells.
    n_legs = len(cell.sequence) - 1
    cell_leg_revs = tuple(cell.per_leg_revs[:n_legs])
    cell_leg_branches = tuple(cell.per_leg_branch[:n_legs])
    if all(r == 0 for r in cell_leg_revs):
        cell_leg_revs = ()
        cell_leg_branches = ()

    primary = PhaseSignature(
        bodies=tuple(cell.sequence),
        leg_durations_s=tuple(s * SECONDS_PER_DAY for s in seed_days),
        vinf_target_kms=vinf_target_per_enc,
        leg_revs=cell_leg_revs,
        leg_branches=cell_leg_branches,
    )
    seeds = leg_duration_seeds(
        bodies=primary.bodies,
        primary_leg_durations_s=primary.leg_durations_s,
        vinf_target_kms=primary.vinf_target_kms,
        period_s=target_period_sec,
        leg_revs=cell_leg_revs,
        leg_branches=cell_leg_branches,
    )

    delta = timedelta(days=10.0 * DAYS_PER_JULIAN_YEAR)
    # Ensure tz-aware bounds for the J2000-relative conversion downstream.
    if priority_date.tzinfo is None:
        priority_date = priority_date.replace(tzinfo=UTC)
    windows = find_candidate_windows(
        seeds,
        ephem,
        (priority_date - delta, priority_date + delta),
        n=n_candidates * len(seeds),
        mismatch_cap_kms=20.0,
    )
    if not windows:
        return None
    return _dt_to_t_sec(windows[0].departure_date)


def optimise_cell_ephemeris(
    cell: Cell,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    priority_date_iso: str | None = None,
    vinf_targets_kms: dict[str, float] | None = None,
    n_laps: int = 5,
    n_starts: int = 5,
    seed: int = 0,
    rp_factors: dict[str, float] | None = None,
    tof_seed_days: Sequence[float] | None = None,
    mode: str = "maintenance",
    scan_epochs: int = 1,
    scan_window_years: float | None = None,
    scan_max_workers: int | None = None,
) -> OptimisationResult:
    """Spec §12(a) ephemeris-mode optimisation over the general engine.

    Optimises a :class:`Cell` to a closed periodic cycler on the *real*
    ephemeris (typically ``Ephemeris("astropy")``) by deriving the
    inputs of the body-agnostic
    :func:`~cyclerfinder.search.maintain.optimise_maintenance_dv` from
    the cell plus a resolved launch epoch, then mapping its
    :class:`~cyclerfinder.search.maintain.MaintenanceOptimResult` back
    onto the same :class:`OptimisationResult` shape that
    :func:`optimise_cell_idealized` produces.

    Closure in the rotating frame is physically unreachable on a real
    ephemeris (documented ``bvp.py:39-51``); the optimisation *objective*
    is the summed flyby turn-deficit maintenance ΔV (carried in
    ``best_score.total_maintenance_dv_kms``), and ``closure_residual_kms``
    carries the same maintenance-ΔV magnitude as a feasibility proxy
    rather than an exact rotating-frame residual.

    Parameters
    ----------
    cell:
        The discrete cell to optimise. Its first and last
        ``sequence`` entries must match (closed loop); open-sequence
        ephemeris cyclers are out of scope (landmine #1).
    ephem:
        Planet-state provider — use ``Ephemeris("astropy")`` for real
        DE440 states.
    vinf_cap:
        Hard ceiling on V∞ at every encounter, km/s (spec §12(d)).
    priority_date_iso:
        ISO-8601 ``"YYYY-MM-DD"`` literature/priority epoch centring the
        real launch-window search. ``None`` ⇒ no epoch can be resolved
        and a non-converged "no real window" result is returned.
    vinf_targets_kms:
        Per-body V∞ targets (km/s) used to phase-match a real launch
        window. Required for epoch resolution (blind discovery without
        targets is out of scope); ``None`` ⇒ non-converged result.
    n_laps:
        Reserved for an optional multi-lap drift diagnostic; not consumed
        by the maintenance-ΔV objective. Kept for API compatibility.
    n_starts, seed:
        Multi-start count and RNG seed threaded into
        :func:`~cyclerfinder.search.maintain.optimise_maintenance_dv`.
    rp_factors:
        Per-body ``SAFE_PERIHELION_KM`` multipliers, threaded into the
        final :func:`~cyclerfinder.model.score.score` call.
    tof_seed_days:
        Optional per-leg time-of-flight seed (days), length
        ``len(cell.sequence) - 1``. Overrides the default equispaced seed
        for both the epoch-resolution leg signature and the maintenance
        guesses, so a strongly asymmetric family (e.g. S1L1's ~154 d
        outbound + long return) phase-matches to its own launch epoch
        rather than a symmetric degenerate basin. Bounds widen to
        ``[0.05, 0.95] * period`` per leg when supplied.
    mode:
        ``"maintenance"`` (default, byte-identical to pre-M-ED): the
        summed flyby turn-deficit maintenance-ΔV objective, answering the
        M7 TCM question; ``closure_residual_kms`` carries the ΔV proxy.
        ``"ballistic"`` (M-ED): runs the N-arc ballistic differential
        corrector (:func:`~cyclerfinder.search.correct.ballistic_correct`)
        — V∞-magnitude continuity at every node driven to zero — and
        reports a **real** closure residual (spec §2.2). Any other value
        raises ``ValueError``.
    scan_epochs:
        Ballistic mode only. ``1`` (default) runs a single start from the
        resolved/priority epoch (byte-identical to pre-scan behaviour).
        ``> 1`` drives the **scan rung** (task #110): an epoch grid of
        ``scan_epochs`` launch seeds across ``scan_window_years`` is run in
        parallel via :mod:`cyclerfinder.search.scan`, and the best closed
        solution (lowest residual, then lowest max-V∞) is returned. This is
        the density lever that family selection across epochs needs (spec §3.4;
        the prototype's ``main()`` scan loop).
    scan_window_years:
        Width of the epoch-scan window (years), centred on the resolved epoch.
        ``None`` ⇒ one full target period (the sourced repeat period). Ignored
        unless ``scan_epochs > 1``.
    scan_max_workers:
        Worker process count for the scan rung. ``None`` ⇒ ``os.cpu_count()``.

    Returns
    -------
    OptimisationResult
        Same shape as :func:`optimise_cell_idealized`. ``converged`` and
        ``constraints_satisfied`` follow the maintenance solve and the
        V∞-cap check; an unresolved real launch window yields a
        non-converged, non-feasible result (not an exception).

    Raises
    ------
    ValueError
        If ``cell.sequence`` is not a closed loop.
    """
    from cyclerfinder.search.maintain import optimise_maintenance_dv
    from cyclerfinder.verify.real_closure import _parse_priority_date

    del n_laps  # reserved for a future multi-lap drift diagnostic
    if rp_factors is None:
        rp_factors = {}

    if mode not in {"maintenance", "ballistic"}:
        raise ValueError(
            f"optimise_cell_ephemeris mode must be 'maintenance' or 'ballistic'; got {mode!r}."
        )

    # Landmine #1: the maintenance chain assumes a closed loop.
    if cell.sequence[0] != cell.sequence[-1]:
        raise ValueError(
            f"optimise_cell_ephemeris requires a closed sequence "
            f"(sequence[0] == sequence[-1]); got {cell.sequence!r}. "
            f"Open-sequence ephemeris cyclers are out of scope.",
        )

    if mode == "ballistic":
        return _optimise_cell_ephemeris_ballistic(
            cell,
            ephem,
            vinf_cap=vinf_cap,
            priority_date_iso=priority_date_iso,
            vinf_targets_kms=vinf_targets_kms,
            n_starts=n_starts,
            rp_factors=rp_factors,
            tof_seed_days=tof_seed_days,
            scan_epochs=scan_epochs,
            scan_window_years=scan_window_years,
            scan_max_workers=scan_max_workers,
        )

    target_period_sec = _target_period_sec(cell)
    seed_days, bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)

    # A caller-supplied asymmetric ToF seed overrides the equispaced one for
    # BOTH the phase-match leg signature (epoch resolution) and the maintenance
    # guesses. This matters for families whose legs are strongly asymmetric
    # (e.g. S1L1's ~154 d outbound + long return): the equispaced seed phase-
    # matches a symmetric leg signature and resolves a launch epoch in a
    # degenerate high-V_inf basin, never reaching the target family.
    if tof_seed_days is not None:
        n_legs = len(cell.sequence) - 1
        if len(tof_seed_days) != n_legs:
            raise ValueError(
                f"tof_seed_days has {len(tof_seed_days)} entries; "
                f"cell.sequence implies {n_legs} legs.",
            )
        seed_days = [float(t) for t in tof_seed_days]
        period_days = target_period_sec / 86400.0
        # Widen bounds to bracket any reasonable asymmetric seed while keeping
        # every leg strictly positive and shorter than the full period.
        bounds = [(0.05 * period_days, 0.95 * period_days)] * n_legs

    # --- Resolve a real launch epoch by phase-matching against the sourced
    # V∞ anchors. Both the priority date and the V∞ targets are required;
    # without them the cell is unanchored on the real ephemeris. STAGE 3: the
    # multi-seed resolver fans the (possibly asymmetric) seed ToFs into a
    # perturbation grid and ranks candidate windows by V∞ mismatch, so an
    # asymmetric family lands its own basin instead of a symmetric degenerate.
    priority = _parse_priority_date(priority_date_iso)
    t0_sec: float | None = None
    if priority is not None and vinf_targets_kms is not None:
        t0_sec = _resolve_t0_multi_seed(
            cell,
            seed_days,
            priority,
            ephem,
            vinf_targets_kms,
            target_period_sec,
            n_candidates=n_starts,
        )

    if t0_sec is None:
        # Honest "no real window" outcome (landmine #3): surface, don't crash.
        sentinel = _sentinel_cycler(cell, ephem, target_period_sec)
        sentinel_score = score(
            sentinel,
            ephem,
            vinf_cap=vinf_cap,
            target_period_sec=target_period_sec,
            rp_factors=rp_factors,
        )
        return OptimisationResult(
            cell=cell,
            best_cycler=sentinel,
            best_score=sentinel_score,
            closure_residual_kms=float("inf"),
            optimiser_history=(),
            converged=False,
            constraints_satisfied=False,
        )

    # --- General real-ephemeris maintenance solve. Thread the cell's per-leg
    # revolution / branch metadata so multi-rev legs (e.g. an Earth-to-Earth
    # resonant interval) are Lambert-solved with the right topology rather than
    # silently flattened to direct legs.
    maint = optimise_maintenance_dv(
        list(cell.sequence),
        ephem,
        t0_guess_sec=t0_sec,
        tof_days_guesses=seed_days,
        tof_bounds_days=bounds,
        per_leg_revs=cell.per_leg_revs,
        per_leg_branch=cell.per_leg_branch,
        synodic_pair=(cell.bodies[0], cell.bodies[1]),
        closure_body=cell.sequence[0],
        n_starts=n_starts,
        seed=seed,
    )

    best_cycler = maint.cycler
    best_score = score(
        best_cycler,
        ephem,
        vinf_cap=vinf_cap,
        target_period_sec=target_period_sec,
        rp_factors=rp_factors,
    )

    # V∞ cap is the spec §12(d) hard inequality; check every encounter.
    vinf_ok = all(
        max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        <= vinf_cap
        for enc in best_cycler.encounters
    )
    constraints_satisfied = bool(maint.converged and vinf_ok)

    # Closure is unreachable on a real ephemeris; carry the maintenance-ΔV
    # magnitude as the feasibility proxy in lieu of an exact residual.
    closure_residual_kms = float(maint.maintenance_dv_kms)

    return OptimisationResult(
        cell=cell,
        best_cycler=best_cycler,
        best_score=best_score,
        closure_residual_kms=closure_residual_kms,
        optimiser_history=(),
        converged=bool(maint.converged),
        constraints_satisfied=constraints_satisfied,
    )


def _optimise_cell_ephemeris_ballistic(
    cell: Cell,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    priority_date_iso: str | None,
    vinf_targets_kms: dict[str, float] | None,
    n_starts: int,
    rp_factors: dict[str, float],
    tof_seed_days: Sequence[float] | None,
    scan_epochs: int = 1,
    scan_window_years: float | None = None,
    scan_max_workers: int | None = None,
) -> OptimisationResult:
    """M-ED ballistic mode (spec §2.2): run the N-arc ballistic differential
    corrector and report a real V∞-continuity closure residual.

    Reuses the maintenance path's epoch resolution
    (:func:`_resolve_t0_multi_seed`) and seed/bounds helpers verbatim, then
    calls :func:`~cyclerfinder.search.correct.ballistic_correct` instead of the
    maintenance-ΔV solver. The corrected ``[t0, *tofs]`` is rebuilt into a
    :class:`Cycler` so downstream consumers read genuine per-encounter V∞,
    and ``closure_residual_kms`` carries the real residual (not the ΔV proxy).
    """
    from cyclerfinder.search.correct import ballistic_correct
    from cyclerfinder.search.maintain import _build_chain
    from cyclerfinder.verify.real_closure import _parse_priority_date

    target_period_sec = _target_period_sec(cell)
    seed_days, _bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)
    n_legs = len(cell.sequence) - 1
    if tof_seed_days is not None:
        if len(tof_seed_days) != n_legs:
            raise ValueError(
                f"tof_seed_days has {len(tof_seed_days)} entries; "
                f"cell.sequence implies {n_legs} legs.",
            )
        seed_days = [float(t) for t in tof_seed_days]

    priority = _parse_priority_date(priority_date_iso)
    t0_sec: float | None = None
    if priority is not None and vinf_targets_kms is not None:
        t0_sec = _resolve_t0_multi_seed(
            cell,
            seed_days,
            priority,
            ephem,
            vinf_targets_kms,
            target_period_sec,
            n_candidates=n_starts,
        )
    if t0_sec is None and priority is not None:
        # The V∞ phase-match resolver cannot seat multi-rev / strongly
        # asymmetric topologies (e.g. S1L1's E-M-E-E with a multi-rev E-E
        # loop) — find_candidate_windows returns no window. The corrector,
        # unlike the maintenance solver, drives a root-find from a *direct*
        # epoch seed (prototype scripts/correct_s1l1_twoarc.py seeds t0 from a
        # fixed date, never the phase-match resolver). Fall back to the
        # priority date as the direct t0 seed so the corrector still runs;
        # family selection across epochs is the seeding ladder's job (spec §3).
        from datetime import UTC

        from cyclerfinder.search.phase_match import _dt_to_t_sec

        if priority.tzinfo is None:
            priority = priority.replace(tzinfo=UTC)
        t0_sec = _dt_to_t_sec(priority)

    if t0_sec is None:
        # Honest "no real window" outcome — surface, don't crash.
        sentinel = _sentinel_cycler(cell, ephem, target_period_sec)
        sentinel_score = score(
            sentinel,
            ephem,
            vinf_cap=vinf_cap,
            target_period_sec=target_period_sec,
            rp_factors=rp_factors,
        )
        return OptimisationResult(
            cell=cell,
            best_cycler=sentinel,
            best_score=sentinel_score,
            closure_residual_kms=float("inf"),
            optimiser_history=(),
            converged=False,
            constraints_satisfied=False,
        )

    # Eliminate the longest seed leg as the period slack leg (spec §2.1(a)).
    slack_leg = int(np.argmax(seed_days)) if seed_days else 0
    free_tof = [t for i, t in enumerate(seed_days) if i != slack_leg]

    if scan_epochs > 1:
        # --- Scan rung (task #110, spec §3.4): a parallel epoch grid across one
        # period (or scan_window_years) from the resolved epoch. Family selection
        # across launch epochs is the density lever the single start lacks.
        corr = _ballistic_scan_rung(
            cell,
            ephem,
            t0_sec=float(t0_sec),
            free_tof=free_tof,
            target_period_sec=target_period_sec,
            slack_leg=slack_leg,
            vinf_cap=vinf_cap,
            rp_factors=rp_factors,
            scan_epochs=scan_epochs,
            scan_window_years=scan_window_years,
            scan_max_workers=scan_max_workers,
        )
    else:
        corr = ballistic_correct(
            sequence=tuple(cell.sequence),
            per_leg_revs=tuple(cell.per_leg_revs),
            per_leg_branch=tuple(cell.per_leg_branch),
            t0_seed_sec=float(t0_sec),
            tof_seed_days=free_tof,
            period_sec=target_period_sec,
            ephem=ephem,
            vinf_cap=vinf_cap,
            rp_factors=rp_factors or None,
            slack_leg=slack_leg,
        )

    # Rebuild the closed cycler at the corrected geometry so consumers read
    # genuine per-encounter V∞ from the same construct path as every other mode.
    x = np.array([corr.t0_sec, *corr.tof_days], dtype=np.float64)
    best_cycler = _build_chain(
        x,
        cell.sequence,
        ephem,
        per_leg_revs=cell.per_leg_revs,
        per_leg_branch=cell.per_leg_branch,
    )
    if best_cycler is None:
        best_cycler = _sentinel_cycler(cell, ephem, target_period_sec)

    best_score = score(
        best_cycler,
        ephem,
        vinf_cap=vinf_cap,
        target_period_sec=target_period_sec,
        rp_factors=rp_factors,
    )

    return OptimisationResult(
        cell=cell,
        best_cycler=best_cycler,
        best_score=best_score,
        closure_residual_kms=float(corr.max_residual_kms),
        optimiser_history=(),
        converged=bool(corr.converged),
        constraints_satisfied=bool(corr.constraints_satisfied),
    )


def _ballistic_scan_rung(
    cell: Cell,
    ephem: Ephemeris,
    *,
    t0_sec: float,
    free_tof: Sequence[float],
    target_period_sec: float,
    slack_leg: int,
    vinf_cap: float,
    rp_factors: dict[str, float],
    scan_epochs: int,
    scan_window_years: float | None,
    scan_max_workers: int | None,
) -> Any:
    """Scan rung (task #110): parallel epoch grid -> best closed corrector result.

    Builds ``scan_epochs`` launch-epoch seeds across one period (or
    ``scan_window_years``) centred on ``t0_sec``, runs the corrector at each in
    parallel via :func:`cyclerfinder.search.scan.scan_parallel`, and returns the
    best result: lowest residual among closed solutions, ties broken by lowest
    max-V∞; if none close, the lowest-residual result overall. The cell's own
    topology is used for every grid point (epoch is the swept axis here; the
    descriptor/topology axis is the caller's via separate cells).
    """
    from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
    from cyclerfinder.search.scan import build_epoch_branch_grid, scan_parallel

    window_days = (
        scan_window_years * DAYS_PER_JULIAN_YEAR
        if scan_window_years is not None
        else target_period_sec / SECONDS_PER_DAY
    )
    half = 0.5 * window_days * SECONDS_PER_DAY
    if scan_epochs == 1:
        t0_seeds = [t0_sec]
    else:
        step = (2.0 * half) / (scan_epochs - 1)
        t0_seeds = [t0_sec - half + i * step for i in range(scan_epochs)]

    grid = build_epoch_branch_grid(
        sequence=tuple(cell.sequence),
        period_sec=target_period_sec,
        vinf_cap=vinf_cap,
        t0_seeds_sec=t0_seeds,
        branch_topologies=[(tuple(cell.per_leg_revs), tuple(cell.per_leg_branch))],
        tof_seed_days=free_tof,
        slack_leg=slack_leg,
        rp_factors=rp_factors or None,
    )
    results = scan_parallel(grid, ephem_model=ephem.model, max_workers=scan_max_workers)
    if not results:
        # No grid -> fall back to a single direct start.
        from cyclerfinder.search.correct import ballistic_correct

        return ballistic_correct(
            sequence=tuple(cell.sequence),
            per_leg_revs=tuple(cell.per_leg_revs),
            per_leg_branch=tuple(cell.per_leg_branch),
            t0_seed_sec=t0_sec,
            tof_seed_days=free_tof,
            period_sec=target_period_sec,
            ephem=ephem,
            vinf_cap=vinf_cap,
            rp_factors=rp_factors or None,
            slack_leg=slack_leg,
        )
    closed = [r for r in results if r.closed]
    pool = closed if closed else results
    best = min(
        pool,
        key=lambda r: (
            r.max_residual_kms,
            max(r.result.vinf_per_encounter_kms)
            if r.result.vinf_per_encounter_kms
            else float("inf"),
        ),
    )
    return best.result


# ---------------------------------------------------------------------------
# Public API: top-level discovery (spec §6)
# ---------------------------------------------------------------------------


def find_cyclers(
    bodies: tuple[str, ...],
    k_synodic: int,
    vinf_cap: float,
    *,
    n_keep: int = 20,
    ephem: Ephemeris | None = None,
    l_max: int = 4,
    n_max: int = 0,
    branch_set: tuple[str, ...] = ("single",),
    n_starts: int = 5,
    seed: int = 0,
    use_de: bool = True,
    rp_factors: dict[str, float] | None = None,
) -> list[OptimisationResult]:
    """Spec §6 top-level discovery interface — the v1 pipeline end-to-end.

    Pipeline (composes M4 + M5):

    1. ``ephem = ephem or Ephemeris("circular")``
    2. ``cells = list(feasible_cells(bodies, l_max, k_synodic, n_max,
       vinf_cap, ephem, branch_set))``
    3. ``results = [optimise_cell_idealized(c, ...) for c in cells]``
    4. Filter to feasible results
       (``constraints_satisfied AND hard_constraints_pass``).
    5. Sort by composite score ascending (best first).
    6. Return the top ``n_keep``.

    Parameters
    ----------
    bodies:
        Canonical body codes for the enumerator. The 2-body case is the
        M5 native case; ``len(bodies) >= 3`` is M8 territory (the
        per-cell ``target_period_sec`` will still resolve via the
        single-pair formula in :func:`_target_period_sec`).
    k_synodic:
        Used both as ``k_max`` for the M4 enumerator AND as the
        per-cell ``period_k`` resolved by :func:`_target_period_sec`.
    vinf_cap:
        Hard ceiling on V∞ at every encounter, km/s. Forwarded to both
        the Tisserand pruning and the per-cell optimiser.
    n_keep:
        Maximum number of results returned.
    ephem:
        Planet-state provider. Defaults to ``Ephemeris("circular")``;
        callers passing ``Ephemeris("astropy")`` get circular-style
        closure semantics because ``optimise_cell_idealized`` is the
        only mode wired in M5.
    l_max:
        Maximum sequence length passed to
        :func:`~cyclerfinder.search.sequence.feasible_cells`. Default
        ``4`` matches the M4 gate test cap.
    n_max:
        Maximum heliocentric revolutions per leg. Default ``0`` (direct
        legs only); ``>= 1`` requires expanding ``branch_set``.
    branch_set:
        Allowed Lambert branches per leg. Default ``("single",)``.
    n_starts, seed, use_de, rp_factors:
        Forwarded to :func:`optimise_cell_idealized`.

    Returns
    -------
    list[OptimisationResult]
        Sorted ascending by composite score; length at most ``n_keep``.
        Every entry satisfies the hard inequalities.

    Notes
    -----
    Single-process per spec §13.6 deferred to M7. Single-pair only —
    multi-body beats (the §3.4 VEM 6.4-yr case) require the §13.4
    multi-pair resonance computation which M8's VEM campaign is the
    right place to add. The ``ephem`` parameter is accepted for forward
    compatibility but the M5 contract calls
    :func:`optimise_cell_idealized` regardless of ``ephem.model``;
    M6/M7's extensions will add a ``mode=`` parameter to switch backends.
    """
    if ephem is None:
        ephem = Ephemeris(model="circular")

    cells = list(
        feasible_cells(
            tuple(bodies),
            l_max=l_max,
            k_max=k_synodic,
            n_max=n_max,
            vinf_cap=vinf_cap,
            ephem=ephem,
            branch_set=branch_set,
        ),
    )

    results: list[OptimisationResult] = []
    for cell in cells:
        try:
            r = optimise_cell_idealized(
                cell,
                ephem,
                vinf_cap=vinf_cap,
                n_starts=n_starts,
                seed=seed,
                use_de=use_de,
                rp_factors=rp_factors,
            )
        except (ValueError, RuntimeError):
            # An optimiser-internal failure on one cell should not block
            # the rest of the enumeration; record as a skip.
            continue
        results.append(r)

    feasible = [
        r for r in results if r.constraints_satisfied and r.best_score.hard_constraints_pass
    ]
    feasible.sort(key=lambda r: composite_score(r.best_score))
    return feasible[:n_keep]


__all__ = [
    "OptimisationResult",
    "find_cyclers",
    "interior_epochs_from_leg_tofs",
    "optimise_cell_ephemeris",
    "optimise_cell_idealized",
]
