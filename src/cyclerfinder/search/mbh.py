"""Monotonic Basin Hopping (MBH) wrapper around the existing correctors (task #145).

Why this exists
---------------
#110/#120/#122 proved that our corrector failures are basin-*selection* problems,
not grid-density problems: the local solver converges, but to a degenerate /
off-family basin. A denser epoch x branch grid does not cure that. MBH is the
trajectory community's standard answer to exactly this disease (Englander &
Conway; the EMTG architecture): instead of gridding seeds, hop from the best
incumbent by a random *perturbation*, re-run the local corrector, and accept the
new point only if it strictly improves the objective. The randomness is used to
escape the current funnel; the local solver does the refining.

The honest framing (survey Thread 1, ``docs/notes/2026-06-07-external-algorithms-
survey.md``): MBH cures selection *within a fixed transcription*. If a target
basin does not exist for the single-ellipse topology, MBH confirms that faster
than a denser grid -- itself a positive (negative-science) result.

Perturbation distribution -- SPEC CAVEAT
----------------------------------------
The canonical perturbation-distribution paper (Englander & Englander 2014,
"Tuning Monotonic Basin Hopping", ISSFD24 S7-3) prescribes long-tailed
(Cauchy / Pareto) per-gene perturbations -- the long tail is what lets the hop
jump *between* basins rather than only jittering within one. That paper is NOT
yet acquired (it is on the #116 acquisition list). Until it is, this module
implements a *documented, sensible default*: per-gene relative perturbation drawn
from a configurable distribution (``"gaussian"``, ``"uniform"``, or ``"cauchy"``)
with a configurable scale. ``"cauchy"`` is provided as the long-tailed option and
is the closest available stand-in for the Englander spec, but its TUNING (scale
schedule, per-gene scaling) is NOT sourced and must not be claimed as such. The
default is ``"cauchy"`` because the long tail is the property that makes MBH hop;
the Gaussian/uniform options are offered for the mechanics gate and for
ablation. Refinement to the exact Englander 2014 tuning is a follow-up once the
paper is in hand.

Design
------
* :func:`mbh` is generic: it takes a closure ``objective_and_solve(x, rng) ->
  MBHStep`` that runs ONE local solve from a seed and reports the resulting point
  + scalar objective + feasibility, plus a starting vector ``x0``. The inner step
  is an existing corrector call (see the adapters below); ``mbh`` never imports a
  corrector itself, keeping it reusable.
* Monotonic acceptance: a hop is accepted iff it is feasible AND its objective is
  strictly below the incumbent by more than ``accept_tol`` (so floating noise
  does not churn the incumbent). The "hop" comes from the perturbation, never
  from accepting an uphill move.
* Full audit trail on :class:`MBHResult`: every hop's objective, accept/reject,
  the rng seed, the best-objective history, and the stall counter -- this repo's
  culture is auditable solver provenance.
* Determinism: ``rng_seed`` is REQUIRED and the only source of randomness is a
  local ``numpy.random.Generator``; no global RNG state is touched.

Adapters
--------
* :func:`make_ballistic_step` wraps
  :func:`cyclerfinder.search.correct.ballistic_correct` (genome ``x = [t0_sec,
  *free_tof_days]``; objective = ``max_residual_kms``).
* :func:`make_free_return_step` wraps
  :func:`cyclerfinder.search.free_return.free_return_correct` (genome
  ``x = [a_au, e, t0_sec]``; objective = ``max_residual_kms``). The #137
  radial-crossing genome; pair with
  :func:`cyclerfinder.search.free_return.seed_ae_from_aphelion_transit` to build
  the descriptor seed.

Pure: depends only on numpy + the existing correctors (wrapped, never edited).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct
from cyclerfinder.search.free_return import free_return_correct

# Distributions the perturbation kernel understands. "cauchy" is the long-tailed
# default (closest available stand-in for the Englander 2014 Cauchy/Pareto spec,
# which is NOT yet acquired -- see module docstring).
PERTURBATIONS = ("cauchy", "gaussian", "uniform")


@dataclass(frozen=True)
class MBHStep:
    """Result of ONE local solve inside an MBH hop.

    Attributes
    ----------
    x:
        The decision vector at the converged local solution (the corrector may
        move the seed; this is where it landed).
    objective:
        Scalar to MINIMISE (lower is better). For both adapters this is the
        corrector's ``max_residual_kms``.
    feasible:
        Whether the local solution is acceptable as an incumbent at all (e.g. the
        corrector ``converged``). An infeasible step can never be accepted, no
        matter how low its objective.
    info:
        Free-form audit payload from the underlying corrector (the full result
        object's salient fields), carried onto the accepted incumbent.
    """

    x: np.ndarray
    objective: float
    feasible: bool
    info: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MBHResult:
    """Outcome + full audit trail of an :func:`mbh` run.

    Attributes
    ----------
    best_x, best_objective, best_info:
        The best accepted incumbent and its corrector payload.
    best_feasible:
        Whether ``best_x`` is feasible. False means NO feasible point was ever
        found (``best_objective`` is then the lowest *infeasible* objective seen,
        for diagnostics) -- the caller MUST check this before trusting the result.
    hops_attempted, hops_accepted:
        Counts over the whole run (the seed solve is hop 0 and is counted as an
        attempt; it is "accepted" iff it produced the first feasible incumbent).
    objective_history:
        Per-hop objective of the *candidate* (length == ``hops_attempted``), so
        the search trace is fully reconstructable.
    accept_history:
        Per-hop accept (True) / reject (False) decision, aligned with
        ``objective_history``.
    best_history:
        Per-hop best-incumbent objective AFTER that hop's decision (the monotone
        non-increasing envelope -- the headline curve).
    final_stall:
        Consecutive non-improving hops at the end of the run (== ``stop_after_
        stall`` if the run stopped on the stall criterion).
    stopped_on_stall:
        Whether the run terminated early on the stall criterion (vs exhausting
        ``n_hops``).
    rng_seed:
        The seed used -- echoed for reproducibility (the run is a pure function
        of (objective_and_solve, x0, params, rng_seed)).
    """

    best_x: np.ndarray
    best_objective: float
    best_feasible: bool
    best_info: dict[str, Any]
    hops_attempted: int
    hops_accepted: int
    objective_history: tuple[float, ...]
    accept_history: tuple[bool, ...]
    best_history: tuple[float, ...]
    final_stall: int
    stopped_on_stall: bool
    rng_seed: int


def _perturb(
    x: np.ndarray,
    rng: np.random.Generator,
    *,
    distribution: str,
    scale: float | Sequence[float] | np.ndarray | None,
    absolute_scale: np.ndarray | None = None,
) -> np.ndarray:
    """Perturb ``x`` componentwise.

    Two sizing modes (a gene uses ABSOLUTE if ``absolute_scale`` gives it a
    finite, positive value, else RELATIVE):

    * RELATIVE (``scale``): step is a fraction of each gene's magnitude (with a
      unit floor so a near-zero gene still moves). One ``scale`` then stays
      meaningful across genes whose magnitudes differ wildly.
    * ABSOLUTE (``absolute_scale``): step is in the gene's OWN physical units.
      This is essential for a gene like ``t0`` (~1e8 seconds): a relative
      fraction there is enormous, so a basin that is only a few days wide is
      unreachable. Sizing ``t0`` absolutely (e.g. a few days) lets the hop
      actually land in such a basin.

    ``scale`` / ``absolute_scale`` may be scalar or per-gene. Use ``0`` (relative)
    or ``NaN``/``0`` (absolute) on a gene to FREEZE it (no perturbation).
    """
    if scale is None:
        scale_arr = np.zeros(x.shape, dtype=np.float64)
    else:
        scale_arr = np.asarray(scale, dtype=np.float64)
        if scale_arr.ndim == 0:
            scale_arr = np.full(x.shape, float(scale_arr))
    if scale_arr.shape != x.shape:
        raise ValueError(f"scale shape {scale_arr.shape} != x shape {x.shape}")

    # Relative step size: fraction of |gene|, floored so a ~0 gene still moves.
    step = scale_arr * np.maximum(np.abs(x), 1.0)

    if absolute_scale is not None:
        abs_arr = np.asarray(absolute_scale, dtype=np.float64)
        if abs_arr.ndim == 0:
            abs_arr = np.full(x.shape, float(abs_arr))
        if abs_arr.shape != x.shape:
            raise ValueError(f"absolute_scale shape {abs_arr.shape} != x shape {x.shape}")
        # A finite positive absolute scale OVERRIDES the relative step for that gene.
        use_abs = np.isfinite(abs_arr) & (abs_arr > 0.0)
        step = np.where(use_abs, abs_arr, step)

    if distribution == "gaussian":
        r = rng.standard_normal(x.shape)
    elif distribution == "uniform":
        r = rng.uniform(-1.0, 1.0, size=x.shape)
    elif distribution == "cauchy":
        # Standard Cauchy via ratio of normals (long-tailed -- the property that
        # lets a hop escape the current basin). NOT the sourced Englander 2014
        # tuning; see the module docstring spec caveat.
        r = rng.standard_cauchy(x.shape)
    else:
        raise ValueError(
            f"unknown perturbation distribution {distribution!r}; expected one of {PERTURBATIONS}"
        )
    return np.asarray(x + step * r, dtype=np.float64)


def mbh(
    objective_and_solve: Callable[[np.ndarray, np.random.Generator], MBHStep],
    x0: Sequence[float] | np.ndarray,
    *,
    n_hops: int,
    perturbation: str = "cauchy",
    perturbation_scale: float | Sequence[float] | None = 0.05,
    perturbation_absolute_scale: Sequence[float] | None = None,
    rng_seed: int,
    accept_tol: float = 1e-9,
    stop_after_stall: int | None = None,
) -> MBHResult:
    """Generic Monotonic Basin Hopping over an arbitrary local-solve closure.

    Parameters
    ----------
    objective_and_solve:
        ``(x_seed, rng) -> MBHStep``. Runs ONE local solve from ``x_seed`` and
        reports where it landed, its scalar objective (lower better), and whether
        it is feasible. The ``rng`` is passed through for solvers that want
        internal randomness; the two provided adapters are deterministic and
        ignore it.
    x0:
        Initial seed decision vector (hop 0 solves from here, unperturbed).
    n_hops:
        Number of PERTURBED hops after the initial seed solve (so the corrector
        runs at most ``n_hops + 1`` times).
    perturbation:
        Perturbation distribution -- one of :data:`PERTURBATIONS`. Default
        ``"cauchy"`` (long-tailed; see the module docstring spec caveat).
    perturbation_scale:
        Relative per-gene step fraction (scalar or per-gene). Default 0.05 (5%).
        ``None`` (or ``0`` on a gene) freezes the relative step for that gene.
    perturbation_absolute_scale:
        Optional per-gene ABSOLUTE step size, in each gene's own units. A finite
        positive entry OVERRIDES the relative step for that gene -- needed when a
        gene's magnitude is huge (e.g. ``t0`` in seconds) so that a few-day basin
        is reachable. ``None`` (default) uses the relative step everywhere.
    rng_seed:
        REQUIRED. Seeds a local ``numpy.random.Generator``; no global RNG state
        is used, so the run is fully reproducible.
    accept_tol:
        A hop is accepted only if it is feasible AND improves the incumbent
        objective by strictly more than this (guards against floating-noise
        churn). Default 1e-9.
    stop_after_stall:
        If set, stop early once this many consecutive PERTURBED hops fail to
        improve the incumbent. ``None`` (default) runs the full ``n_hops``.

    Returns
    -------
    MBHResult
        Best incumbent + the full audit trail. Check ``best_feasible`` before
        trusting ``best_x`` (False => no feasible point was found).
    """
    if perturbation not in PERTURBATIONS:
        raise ValueError(
            f"unknown perturbation distribution {perturbation!r}; expected one of {PERTURBATIONS}"
        )
    if n_hops < 0:
        raise ValueError("n_hops must be >= 0")
    rng = np.random.default_rng(rng_seed)
    x0_arr = np.asarray(x0, dtype=np.float64)
    abs_scale = (
        np.asarray(perturbation_absolute_scale, dtype=np.float64)
        if perturbation_absolute_scale is not None
        else None
    )

    objective_history: list[float] = []
    accept_history: list[bool] = []
    best_history: list[float] = []

    # Hop 0: solve from the unperturbed seed. This establishes the first
    # incumbent (feasible or not). It counts as an attempt; it is "accepted" iff
    # it yields the first feasible incumbent.
    seed_step = objective_and_solve(x0_arr, rng)
    best_x = np.asarray(seed_step.x, dtype=np.float64)
    best_obj = float(seed_step.objective)
    best_feasible = bool(seed_step.feasible)
    best_info = dict(seed_step.info)
    objective_history.append(best_obj)
    accept_history.append(best_feasible)
    best_history.append(best_obj if best_feasible else float("inf"))

    stall = 0
    stopped_on_stall = False
    for _ in range(n_hops):
        # Perturb the BEST incumbent's location (the classic MBH "hop from the
        # incumbent" -- the perturbation, not an uphill accept, is the escape).
        x_seed = _perturb(
            best_x,
            rng,
            distribution=perturbation,
            scale=perturbation_scale,
            absolute_scale=abs_scale,
        )
        step = objective_and_solve(x_seed, rng)
        cand_obj = float(step.objective)
        objective_history.append(cand_obj)

        # Monotonic acceptance: must be feasible AND strictly improve.
        # If we have no feasible incumbent yet, the first feasible candidate is
        # accepted regardless of the (infeasible) incumbent's objective.
        improves = step.feasible and (not best_feasible or cand_obj < best_obj - accept_tol)
        if improves:
            best_x = np.asarray(step.x, dtype=np.float64)
            best_obj = cand_obj
            best_feasible = True
            best_info = dict(step.info)
            accept_history.append(True)
            stall = 0
        else:
            accept_history.append(False)
            stall += 1

        best_history.append(best_obj if best_feasible else float("inf"))

        if stop_after_stall is not None and stall >= stop_after_stall:
            stopped_on_stall = True
            break

    hops_attempted = len(objective_history)
    hops_accepted = sum(accept_history)
    return MBHResult(
        best_x=best_x,
        best_objective=best_obj,
        best_feasible=best_feasible,
        best_info=best_info,
        hops_attempted=hops_attempted,
        hops_accepted=hops_accepted,
        objective_history=tuple(objective_history),
        accept_history=tuple(accept_history),
        best_history=tuple(best_history),
        final_stall=stall,
        stopped_on_stall=stopped_on_stall,
        rng_seed=int(rng_seed),
    )


# ---------------------------------------------------------------------------
# Concrete adapters: wrap the existing correctors as objective_and_solve closures.
# These import and CALL the correctors; they never modify them.
# ---------------------------------------------------------------------------


def make_ballistic_step(
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_sec: float,
    ephem: Ephemeris,
    vinf_cap: float,
    slack_leg: int | None = None,
    tol_kms: float = 0.1,
    residual_mode: str = "magnitude",
    rp_factors: dict[str, float] | None = None,
) -> Callable[[np.ndarray, np.random.Generator], MBHStep]:
    """Adapter: :func:`cyclerfinder.search.correct.ballistic_correct` as an
    MBH local-solve closure.

    Genome ``x = [t0_sec, *free_tof_days]`` (the same free-variable layout
    ``ballistic_correct`` consumes). Objective = ``max_residual_kms``; feasible =
    ``converged``. The corrector is deterministic, so the ``rng`` argument is
    accepted (for the generic signature) and ignored.
    """

    def step(x: np.ndarray, rng: np.random.Generator) -> MBHStep:
        t0 = float(x[0])
        free_tof = [float(v) for v in x[1:]]
        r = ballistic_correct(
            sequence=sequence,
            per_leg_revs=per_leg_revs,
            per_leg_branch=per_leg_branch,
            t0_seed_sec=t0,
            tof_seed_days=free_tof,
            period_sec=period_sec,
            ephem=ephem,
            vinf_cap=vinf_cap,
            slack_leg=slack_leg,
            tol_kms=tol_kms,
            residual_mode=residual_mode,
            rp_factors=rp_factors,
        )
        # The corrector returns the FULL (slack-reinserted) ToFs; rebuild the
        # free-leg layout so the landed x is in the same coordinates as the seed.
        eff_slack = slack_leg
        if eff_slack is None:
            eff_slack = int(np.argmax(free_tof)) if free_tof else 0
        full = list(r.tof_days)
        landed_free = [full[i] for i in range(len(full)) if i != eff_slack]
        landed_x = np.array([r.t0_sec, *landed_free], dtype=np.float64)
        return MBHStep(
            x=landed_x,
            objective=float(r.max_residual_kms),
            feasible=bool(r.converged),
            info={
                "vinf_per_encounter_kms": tuple(r.vinf_per_encounter_kms),
                "bend_feasible": bool(r.bend_feasible),
                "vinf_cap_ok": bool(r.vinf_cap_ok),
                "tof_days": tuple(r.tof_days),
                "hyperbolic_impossible": bool(r.hyperbolic_impossible),
            },
        )

    return step


def make_free_return_step(
    *,
    period_sec: float,
    ephem: Ephemeris,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float | None = None,
    tol_kms: float = 0.1,
) -> Callable[[np.ndarray, np.random.Generator], MBHStep]:
    """Adapter: :func:`cyclerfinder.search.free_return.free_return_correct` as an
    MBH local-solve closure (the #137 radial-crossing genome).

    Genome ``x = [a_au, e, t0_sec]``. Objective = ``max_residual_kms``; feasible =
    ``converged``. The emerged per-body V_inf is carried in ``info`` (it is the
    EVIDENCE the free-return gate compares against the sourced anchor; it is never
    used as the objective, which would impose it). Deterministic: ``rng`` ignored.
    """

    def step(x: np.ndarray, rng: np.random.Generator) -> MBHStep:
        a_au = float(x[0])
        e = float(x[1])
        t0 = float(x[2])
        kwargs: dict[str, Any] = {"bodies": bodies, "tol_kms": tol_kms}
        if mu is not None:
            kwargs["mu"] = mu
        r = free_return_correct(
            t0_seed_sec=t0,
            a_seed_au=a_au,
            e_seed=e,
            period_sec=period_sec,
            ephem=ephem,
            **kwargs,
        )
        landed_x = np.array([r.a_au, r.e, r.t0_sec], dtype=np.float64)
        return MBHStep(
            x=landed_x,
            objective=float(r.max_residual_kms),
            feasible=bool(r.converged),
            info={
                "a_au": float(r.a_au),
                "e": float(r.e),
                "vinf_kms": dict(r.vinf_kms),
                "transfer_tof_days": float(r.transfer_tof_days),
                "ee_interval_days": float(r.ee_interval_days),
                "solver_nfev": int(r.solver_nfev),
            },
        )

    return step
