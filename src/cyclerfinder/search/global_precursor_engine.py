"""#430 Unified global MGA-DSM precursor search engine.

Replaces the precursor matcher's local-optimiser path with a global
differential_evolution search over (launch epoch, per-leg TOFs, per-leg DSM),
seeded with eccentric-body Tisserand-Poincaré candidates, ranking survivors by
dv_band / total ΔV. See docs/superpowers/specs/2026-06-23-global-precursor-engine-design.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import Bounds, differential_evolution

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import (
    DSMSpec,
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.tisserand_mga_window import (
    MGAChainCandidate,
    _add_days_utc,
    find_mga_chains,
)
from cyclerfinder.verify.dv_band_acceptance import classify_dv_band


def eccentric_tp_linkable_radius_au(body: str, t_sec: float, ephemeris: Ephemeris) -> float:
    """Body's ACTUAL heliocentric radius (AU) at ``t_sec``, for the eccentric
    Tisserand-Poincaré graph (Campagnola-Russell 2009 Part B). The T-P contour
    drawn at the real encounter radius shifts/widens the linkable set vs the
    mean-``a`` circular form. Reduces to ``sma_au`` on the circular backend."""
    r_km, _v = ephemeris.state(body, t_sec)
    return float(np.linalg.norm(np.asarray(r_km, dtype=np.float64))) / AU_KM


def eccentric_tp_seeds(
    *,
    first_body: str,
    seed_vinf_kms: float,
    launch_window: tuple[str, str],
    ephemeris: Ephemeris,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    max_legs: int = 3,
    vinf_grid_kms: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    epoch_step_days: float = 60.0,
    vinf_terminal_tol_kms: float = 0.8,
) -> list[MGAChainCandidate]:
    """Enumerate Earth-launched MGA chains, filtered to those terminating at
    ``first_body`` with terminal V_inf within ``vinf_terminal_tol_kms`` of
    ``seed_vinf_kms``. Returns the DE init population (MGAChainCandidate list).
    Eccentric-body (real-radius) re-screening via
    :func:`eccentric_tp_linkable_radius_au` is a planned follow-on consumer (it
    is not applied here — the enumerator order is preserved).

    Notes
    -----
    The real :func:`find_mga_chains` signature (tisserand_mga_window.py:615)
    takes ``launch_window`` and ``planet_set`` as positional args, is pure
    geometry (it does NOT accept an ``ephemeris`` argument — the ephemeris is
    consumed later by the Phase-1 closure/validation functions), and returns a
    lazy ``Iterator[MGAChainCandidate]``. The ``ephemeris`` parameter here is
    retained for the eccentric-radius linkability (see
    :func:`eccentric_tp_linkable_radius_au`) and for downstream DE seeding; it
    is not forwarded into the pure-geometry enumerator.
    """
    candidates = find_mga_chains(
        launch_window,
        tuple(dict.fromkeys((first_body, *intermediate_bodies))),
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
    )
    out: list[MGAChainCandidate] = []
    for c in candidates:
        if c.sequence[-1] != first_body:
            continue
        if abs(c.vinf_tuple_kms[-1] - seed_vinf_kms) > vinf_terminal_tol_kms:
            continue
        out.append(c)
    return out


_COST_FLOOR = 1.0e6
_DSM_EPS_KMS = 1.0e-4  # below this magnitude a leg is treated as ballistic (no DSMSpec)
_FRACTION_CLAMP = 1.0e-3  # keep fraction_along_leg strictly inside (0, 1)


@dataclass(frozen=True)
class DecisionEval:
    """Result of evaluating one DE decision vector against the real-DE440 oracle."""

    closure: EpochLockedClosure
    total_dsm_dv_kms: float
    per_leg_dsm_kms: tuple[float, ...]
    feasible: bool


def evaluate_decision_vector(
    x: list[float] | np.ndarray,
    *,
    sequence: tuple[str, ...],
    seed_launch_epoch_utc: str,
    vinf_expected_kms: tuple[float, ...],
    ephemeris: Ephemeris,
    inserts_into: str,
    max_revs: int = 2,
) -> DecisionEval:
    """Map a flat decision vector to a closed :class:`EpochLockedTrajectory`.

    Decision-vector layout for an ``N``-leg chain (``len(sequence) == N + 1``)::

        x = [epoch_offset_days,
             tof_1 .. tof_N,
             (eta_i, dvx_i, dvy_i, dvz_i) * N]

    i.e. length ``1 + N + 4 * N``. ``epoch_offset_days`` is added to
    ``seed_launch_epoch_utc``; each per-leg ``(eta, dv)`` block places one
    optional DSM (only legs whose ``|dv| >= _DSM_EPS_KMS`` get a
    :class:`DSMSpec`, so an all-zero DSM block is a ballistic no-op). The
    trajectory is closed against the reused real-ephemeris oracle
    (:func:`close_epoch_locked`).

    Non-convergence / infeasible geometry (non-positive TOF, constructor or
    closure failure) returns a :class:`DecisionEval` with ``_COST_FLOOR``
    residuals and ``feasible=False``.
    """
    xv = np.asarray(x, dtype=np.float64)
    n_legs = len(sequence) - 1
    epoch_offset = float(xv[0])
    tofs = tuple(float(t) for t in xv[1 : 1 + n_legs])
    if any(t <= 0.0 for t in tofs):
        return _infeasible(sequence, vinf_expected_kms, seed_launch_epoch_utc, inserts_into)
    dsm_block = xv[1 + n_legs :]
    dsm_specs: list[DSMSpec] = []
    per_leg_dsm: list[float] = []
    for i in range(n_legs):
        eta, dvx, dvy, dvz = (float(v) for v in dsm_block[4 * i : 4 * i + 4])
        mag = float(np.linalg.norm((dvx, dvy, dvz)))
        per_leg_dsm.append(mag)
        if mag < _DSM_EPS_KMS:
            continue
        frac = min(1.0 - _FRACTION_CLAMP, max(_FRACTION_CLAMP, eta))
        dsm_specs.append(
            DSMSpec(leg_index=i, fraction_along_leg=frac, delta_v_kms=(dvx, dvy, dvz)),
        )
    launch = _add_days_utc(seed_launch_epoch_utc, epoch_offset)
    end_utc = _add_days_utc(launch, sum(tofs))
    try:
        traj = EpochLockedTrajectory(
            sequence=sequence,
            leg_tofs_days=tofs,
            vinf_kms_at_encounters=vinf_expected_kms,
            launch_epoch_utc=launch,
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc=launch,
            validity_window_end_utc=end_utc,
            inserts_into=inserts_into,
            dsm_specs=tuple(dsm_specs) if dsm_specs else None,
        )
        closure = close_epoch_locked(
            traj,
            ephemeris,
            closure_tol_kms=1.0e6,
            flyby_continuity_tol_kms=1.0e6,
            independent_cross_check=False,
            independent_tol_kms=1.0e6,
            max_revs=max_revs,
        )
    except Exception:
        return _infeasible(sequence, vinf_expected_kms, seed_launch_epoch_utc, inserts_into)
    return DecisionEval(
        closure=closure,
        total_dsm_dv_kms=float(sum(per_leg_dsm)),
        per_leg_dsm_kms=tuple(per_leg_dsm),
        feasible=True,
    )


def _infeasible(
    sequence: tuple[str, ...],
    vinf_expected_kms: tuple[float, ...],
    launch: str,
    inserts_into: str,
) -> DecisionEval:
    """Sentinel :class:`DecisionEval` for an infeasible vector: ``_COST_FLOOR`` residuals.

    Assumes ``len(vinf_expected_kms) == len(sequence)`` (the same caller
    precondition the feasible path enforces); a mismatch raises here from the
    sentinel's own ``EpochLockedTrajectory`` constructor rather than being
    silently absorbed, which is the intended behaviour for a caller bug.
    """
    n = len(sequence)
    end = _add_days_utc(launch, 1.0)
    traj = EpochLockedTrajectory(
        sequence=sequence,
        leg_tofs_days=tuple(1.0 for _ in range(n - 1)),
        vinf_kms_at_encounters=vinf_expected_kms,
        launch_epoch_utc=launch,
        orbit_class="precursor_mga",
        n_returns=1,
        validity_window_start_utc=launch,
        validity_window_end_utc=end,
        inserts_into=inserts_into,
    )
    closure = EpochLockedClosure(
        trajectory=traj,
        closure_residual_kms=_COST_FLOOR,
        flyby_continuity_max_dv_kms=_COST_FLOOR,
        per_leg_lambert_solutions=(),
        per_encounter_vinf_kms=tuple(0.0 for _ in range(n)),
        independent_check_residual_kms=None,
        converged=False,
        dsm_delta_v_kms_per_leg=(),
    )
    return DecisionEval(
        closure=closure,
        total_dsm_dv_kms=_COST_FLOOR,
        per_leg_dsm_kms=(),
        feasible=False,
    )


def decision_cost(ev: DecisionEval, *, w_cont: float = 1.0, w_dsm: float = 0.5) -> float:
    """Scalar objective: closure + w_cont*continuity + w_dsm*total_DSM (km/s).
    w_dsm > 0 makes a ballistic solution always score below an equal powered one."""
    return (
        ev.closure.closure_residual_kms
        + w_cont * ev.closure.flyby_continuity_max_dv_kms
        + w_dsm * ev.total_dsm_dv_kms
    )


def rank_band(total_dsm_dv_kms: float) -> str:
    """dv_band for a candidate's total DSM ΔV, via the sourced classifier.

    NB: ``classify_dv_band`` is a *7-cycle maintenance*-ΔV classifier; here we
    reinterpret its sourced bin edges (1/10/300 m/s) for a one-shot *insertion*
    DSM total. At the default ``n_cycles`` no pro-rata scaling occurs, so the
    bins apply as-is — we reuse the sourced thresholds rather than inventing new
    ones, accepting that the budget being binned is insertion not maintenance."""
    return classify_dv_band(total_dsm_dv_kms * 1000.0)


def _bounds_for(
    n_legs: int,
    *,
    epoch_half_width_days: float,
    tof_box_days_per_leg: tuple[float, float],
    dsm_max_kms: float,
) -> Bounds:
    """Decision-box bounds matching :func:`evaluate_decision_vector`'s layout:
    ``[epoch_offset, tof_1..tof_N, (eta, dvx, dvy, dvz) * N]``."""
    lo = [-epoch_half_width_days]
    hi = [epoch_half_width_days]
    lo += [tof_box_days_per_leg[0]] * n_legs
    hi += [tof_box_days_per_leg[1]] * n_legs
    for _ in range(n_legs):  # (eta, dvx, dvy, dvz) per leg
        lo += [0.0, -dsm_max_kms, -dsm_max_kms, -dsm_max_kms]
        hi += [1.0, dsm_max_kms, dsm_max_kms, dsm_max_kms]
    return Bounds(lo, hi)


def search_sequence(
    *,
    sequence: tuple[str, ...],
    seed_launch_epoch_utc: str,
    vinf_expected_kms: tuple[float, ...],
    ephemeris: Ephemeris,
    inserts_into: str,
    epoch_half_width_days: float = 120.0,
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    dsm_max_kms: float = 2.0,
    max_revs: int = 2,
    popsize: int = 15,
    maxiter: int = 60,
    seed: int = 0,
    w_cont: float = 1.0,
    w_dsm: float = 0.5,
) -> DecisionEval:
    """Global differential_evolution over the decision box for ONE sequence.
    Returns the best DecisionEval."""
    n_legs = len(sequence) - 1
    bounds = _bounds_for(
        n_legs,
        epoch_half_width_days=epoch_half_width_days,
        tof_box_days_per_leg=tof_box_days_per_leg,
        dsm_max_kms=dsm_max_kms,
    )

    def _obj(x: np.ndarray) -> float:
        ev = evaluate_decision_vector(
            x,
            sequence=sequence,
            seed_launch_epoch_utc=seed_launch_epoch_utc,
            vinf_expected_kms=vinf_expected_kms,
            ephemeris=ephemeris,
            inserts_into=inserts_into,
            max_revs=max_revs,
        )
        return decision_cost(ev, w_cont=w_cont, w_dsm=w_dsm)

    # workers=1: this repo has documented BLAS reduction-order noise (see the
    # tests/nbody flake fixes); single-worker keeps the DE bit-reproducible
    # under a fixed seed so the determinism test holds.
    result = differential_evolution(
        _obj,
        bounds,
        popsize=popsize,
        maxiter=maxiter,
        seed=seed,
        polish=True,
        tol=1e-8,
        workers=1,
    )
    return evaluate_decision_vector(
        result.x,
        sequence=sequence,
        seed_launch_epoch_utc=seed_launch_epoch_utc,
        vinf_expected_kms=vinf_expected_kms,
        ephemeris=ephemeris,
        inserts_into=inserts_into,
        max_revs=max_revs,
    )


@dataclass(frozen=True)
class PrecursorSurvivor:
    sequence: tuple[str, ...]
    eval: DecisionEval
    dv_band: str
    cost: float


def search_precursors(
    *,
    cycler_id: str,
    first_body: str,
    seed_vinf_kms: float,
    launch_window: tuple[str, str],
    ephemeris: Ephemeris,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    max_legs: int = 3,
    vinf_grid_kms: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    epoch_step_days: float = 60.0,
    dsm_max_kms: float = 2.0,
    max_revs: int = 2,
    popsize: int = 15,
    maxiter: int = 60,
    seed: int = 0,
) -> list[PrecursorSurvivor]:
    """Enumerate sequences (eccentric-T-P seeds), run a global DE per distinct
    sequence, rank survivors by total ΔV ascending (ballistic first)."""
    seeds = eccentric_tp_seeds(
        first_body=first_body,
        seed_vinf_kms=seed_vinf_kms,
        launch_window=launch_window,
        ephemeris=ephemeris,
        intermediate_bodies=intermediate_bodies,
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
    )
    survivors: list[PrecursorSurvivor] = []
    seen: set[tuple[str, ...]] = set()
    for s in seeds:
        if s.sequence in seen:
            continue
        seen.add(s.sequence)
        ev = search_sequence(
            sequence=s.sequence,
            seed_launch_epoch_utc=s.launch_epoch_utc,
            vinf_expected_kms=s.vinf_tuple_kms,
            ephemeris=ephemeris,
            inserts_into=cycler_id,
            tof_box_days_per_leg=tof_box_days_per_leg,
            dsm_max_kms=dsm_max_kms,
            max_revs=max_revs,
            popsize=popsize,
            maxiter=maxiter,
            seed=seed,
        )
        survivors.append(
            PrecursorSurvivor(
                sequence=s.sequence,
                eval=ev,
                dv_band=rank_band(ev.total_dsm_dv_kms),
                cost=decision_cost(ev),
            )
        )
    survivors.sort(key=lambda s: s.eval.total_dsm_dv_kms)
    return survivors
