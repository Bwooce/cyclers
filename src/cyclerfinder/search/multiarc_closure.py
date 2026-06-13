"""Multi-arc closure harness (#248) — clean, reliable, epoch-safe.

Answers, for a multi-arc E-E-M-M / E-E-E-M-M catalogue row, the single question
the throwaway ``scratch_*`` scripts could not answer reliably: *does the row
converge below the 0.1 km/s gate on the Lambert optimiser lane, robustly?* A
converged row is a **SILVER candidate -> gauntlet**, never a catalogue writeback
(the V0-V5 gauntlet still governs).

This module is the production replacement for the scratch experiments. It designs
out the four infrastructure failures that produced the false "0.163 km/s" signal
and the three agent hangs (see ``docs/notes/2026-06-14-multiarc-closure-harness-
design.md`` and the ``2026-06-13-multiarc-seed-basin-fix.md`` companion):

1. **ONE canonical metric.** Everything reports
   :attr:`cyclerfinder.search.dsm_leg.DsmChainResult.max_residual_kms` /
   ``.converged`` against ``tol_kms=0.1``. No ad-hoc "res", no ``dv_total``
   masquerading as residual (the 0.163-vs-2.06 disagreement was a metric mismatch).
2. **Epoch-range-SAFE evaluation.** :func:`safe_chain_residual` wraps the corrector
   so a ``jplephem`` ``OutOfRangeError`` (or any ephemeris-range error — DE440
   covers 1549-2650) returns a large PENALTY residual instead of crashing the run.
   The optimiser is then free to explore wild epochs.
3. **Correct resonant-leg seeding.** :func:`resonant_return_seeds` seeds the transit
   leg via :func:`cyclerfinder.search.self_seeding.joint_epoch_tof_close` (the
   structured-epoch seed) and ENUMERATES the discrete small-N resonant returns of
   each resonant leg's body (descriptor resonant period x {1,2,...,N} plus the
   descriptor value), NOT a +-% scale of one value (which never reaches the right
   return and gives identical junk).
4. **Multi-start.** :func:`close_multiarc_row` drives the structured-epoch x
   discrete-resonant-return seed set through repeated epoch-safe canonical solves
   (the MBH wrapper is available for basin hops within a seed), sampling basins
   systematically instead of trusting one charged seed.

Pure-ish: depends only on numpy + the existing seed/corrector machinery, which it
CALLS and never edits.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import dsm_leg
from cyclerfinder.search.dsm_descriptor_seed import seed_dsm_chain_from_descriptor
from cyclerfinder.search.self_seeding import FamilyAnchors, joint_epoch_tof_close

DAY_S = 86400.0
YEAR_DAYS = 365.25

# Penalty residual returned (instead of crashing) when a decision vector drives the
# corrector's epoch search outside the ephemeris kernel's coverage. Large enough
# that the optimiser/grid always rejects it, finite so it never poisons a min().
PENALTY_RESIDUAL_KMS: float = 1.0e6

# Default number of discrete resonant returns enumerated per resonant leg
# (descriptor period x {1, .., N}). Small-N: the resonant loop arcs of a multi-arc
# cycler are a few body-periods, not dozens.
_DEFAULT_N_RESONANT: int = 4


# ---------------------------------------------------------------------------
# Resolving the jplephem out-of-range exception robustly. astropy's DE440 reader
# raises jplephem's OutOfRangeError for epochs outside kernel coverage; import it
# defensively (the circular backend never needs jplephem) and fall back to a
# never-matching sentinel so the except-tuple is always well-formed.
# ---------------------------------------------------------------------------
try:  # pragma: no cover - import shape depends on the optional astropy/jplephem dep
    from jplephem.exceptions import OutOfRangeError as _OutOfRangeError
except Exception:  # pragma: no cover - jplephem absent (circular-only environments)

    class _OutOfRangeError(Exception):  # type: ignore[no-redef]
        """Sentinel used when jplephem is unavailable (never raised by our code)."""


# Errors that mean "this decision vector wandered outside the ephemeris range" and
# must become a penalty, not a crash. OutOfRangeError is the documented culprit
# (astropy's DE440 reader for epochs past the kernel coverage); ValueError covers
# the astropy/erfa range guards that surface the same out-of-range condition by a
# different type. dsm_chain_correct already absorbs every per-leg LambertError /
# KeplerError / ValueError internally and returns an INFEASIBLE result, so a
# ValueError that escapes the corrector is an ephemeris-range surface, not a
# corrector logic bug — safe to penalise here (design point 2: "any ephemeris-range
# error returns a penalty").
_EPHEMERIS_RANGE_ERRORS: tuple[type[BaseException], ...] = (_OutOfRangeError, ValueError)


def safe_chain_residual(
    x: NDArray[np.float64],
    *,
    sequence: tuple[str, ...],
    ephem: Ephemeris,
    bounds: dsm_leg.DsmBounds | None = None,
    tol_kms: float = 0.1,
    gradient: str = "lambert",
    charge_flyby_continuity: bool = True,
    max_nfev: int = 200,
) -> tuple[float, bool, dsm_leg.DsmChainResult | None]:
    """Epoch-safe canonical-metric wrapper around :func:`dsm_leg.dsm_chain_correct`.

    Runs the corrector from decision vector ``x`` and returns the CANONICAL triple
    ``(max_residual_kms, converged, result_or_None)`` — the single metric the gate
    is decided on everywhere in this harness (design point 1).

    If the corrector drives the ephemeris epoch search outside the kernel's coverage
    (``jplephem`` ``OutOfRangeError`` or an equivalent range ``ValueError`` — DE440
    covers 1549-2650), this returns ``(PENALTY_RESIDUAL_KMS, False, None)`` instead
    of letting the exception abort the whole run (design point 2). The optimiser/grid
    is then free to explore wild epochs; an out-of-range seed simply scores badly.

    Parameters
    ----------
    x:
        The chain decision vector (see :func:`dsm_leg.dsm_chain_decision_vector`).
    sequence, ephem, bounds, tol_kms, gradient, charge_flyby_continuity, max_nfev:
        Passed straight through to :func:`dsm_leg.dsm_chain_correct`. ``gradient``
        defaults to ``"lambert"`` (the convergence reference lane, #245);
        ``charge_flyby_continuity`` defaults to ``True`` (the only mode that rewards
        the bend-feasible low-V_inf basin — the descriptor-seed lane's mode).

    Returns
    -------
    (max_residual_kms, converged, result_or_None)
        The canonical metric, the converged flag, and the full ``DsmChainResult``
        (``None`` only on an out-of-range penalty).
    """
    try:
        result = dsm_leg.dsm_chain_correct(
            np.asarray(x, dtype=np.float64),
            sequence=sequence,
            ephem=ephem,
            bounds=bounds,
            tol_kms=tol_kms,
            charge_flyby_continuity=charge_flyby_continuity,
            gradient=gradient,
            max_nfev=max_nfev,
        )
    except _EPHEMERIS_RANGE_ERRORS:
        # The corrector's epoch search left the ephemeris kernel's coverage. Treat
        # the seed as infeasible with a large finite penalty rather than crashing
        # the campaign (the unguarded crash that aborted the inline probe).
        return (PENALTY_RESIDUAL_KMS, False, None)

    return (float(result.max_residual_kms), bool(result.converged), result)


# ---------------------------------------------------------------------------
# Resonant-return seed enumeration (design point 3).
# ---------------------------------------------------------------------------


def _body_period_days(body: str) -> float:
    """Sidereal orbital period of ``body`` in days (Kepler from the SMA)."""
    a_km = PLANETS[body].sma_au * AU_KM
    return float(2.0 * np.pi * np.sqrt(a_km**3 / MU_SUN_KM3_S2) / DAY_S)


def _transit_leg_index(sequence: tuple[str, ...]) -> int:
    """Index of the unique E<->M transit leg in the sequence."""
    for i in range(len(sequence) - 1):
        if {sequence[i], sequence[i + 1]} == {"E", "M"}:
            return i
    raise ValueError(f"no E<->M transit leg in {sequence}")


def _resonant_leg_indices(sequence: tuple[str, ...]) -> tuple[int, ...]:
    """Indices of the same-body (resonant return) legs, e.g. E-E or M-M."""
    return tuple(i for i in range(len(sequence) - 1) if sequence[i] == sequence[i + 1])


def resonant_return_tof_grid(
    body: str,
    descriptor_tof_days: float,
    *,
    n_resonant: int = _DEFAULT_N_RESONANT,
) -> tuple[float, ...]:
    """Discrete resonant-return ToF candidates for one resonant leg (design point 3).

    The integer multiples of the body's sidereal period (``period x {1, .., N}``)
    PLUS the row's descriptor value itself, de-duplicated and sorted. The integer
    multiples are the physical small-N resonant returns of the body; the descriptor
    value is the row's own tabulated arc duration. This is the discrete ENUMERATION
    the design mandates in place of a +-% scale of one value (which never reaches the
    right return and gives identical-across-seeds junk).
    """
    period = _body_period_days(body)
    cands = {round(float(k) * period, 6) for k in range(1, n_resonant + 1)}
    cands.add(round(float(descriptor_tof_days), 6))
    return tuple(sorted(cands))


@dataclass(frozen=True)
class MultiArcSeedSpec:
    """A fully-built multi-arc seed: decision vector + bounds + provenance.

    Attributes
    ----------
    x0:
        The chain decision vector (see :func:`dsm_leg.dsm_chain_decision_vector`).
    bounds:
        The box bounds for this seed (epoch recentred on the structured seed, ToF
        box widened to admit the seeded resonant returns).
    sequence:
        The body sequence the corrector runs over.
    transit_tof_days:
        The structured-epoch transit-leg ToF used (the working part).
    resonant_tof_days_per_leg:
        The resonant-return ToF assigned to each resonant leg, by leg index.
    """

    x0: NDArray[np.float64]
    bounds: dsm_leg.DsmBounds
    sequence: tuple[str, ...]
    transit_tof_days: float
    resonant_tof_days_per_leg: dict[int, float]


def resonant_return_seeds(
    row: dict[str, Any],
    eph: Ephemeris,
    t_center_sec: float,
    *,
    n_resonant: int = _DEFAULT_N_RESONANT,
    epoch_halfwidth_days: float = 400.0,
) -> Iterator[MultiArcSeedSpec]:
    """Yield multi-arc seed specs = structured transit epoch x discrete resonant ToFs.

    The transit leg is seeded from the structured self-seed
    (:func:`cyclerfinder.search.self_seeding.joint_epoch_tof_close`, the working
    pattern reused verbatim from ``scratch_structured_seed_experiment.structured_
    seed``); each RESONANT leg's ToF is crossed over its discrete small-N resonant
    returns (:func:`resonant_return_tof_grid`). The Cartesian product of the
    per-resonant-leg grids is the multi-start seed set (design points 3 + 4).

    Yields nothing (a clean empty iterator) for rows without a per-arc descriptor or
    without a structured transit-leg seed — the same clean-negative discipline the
    rest of the closure stack uses.
    """
    base = seed_dsm_chain_from_descriptor(row)
    if base is None:
        return
    seq = base.sequence
    n_legs = len(seq) - 1

    vinf = {e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])}
    if "E" not in vinf or "M" not in vinf:
        return

    sig_list = sorted(
        {float(t) for t in (row.get("invariants") or {}).get("transit_times_days") or []}
    )
    if not sig_list:
        return
    sig = sig_list[0]

    anchors = FamilyAnchors(vinf_e=vinf["E"], vinf_m=vinf["M"], vinf_band_kms=1.5)
    ss = joint_epoch_tof_close(eph, anchors, t_center_sec, sig, max_revs=2)
    if ss is None:
        return

    try:
        t_idx = _transit_leg_index(seq)
    except ValueError:
        return
    transit_tof_days = float(ss.tof_g_days)
    resonant_idxs = _resonant_leg_indices(seq)

    # Descriptor resonant period (free_return_arcs[0] = g = the resonant return),
    # in days — the per-leg descriptor value blended into each leg's grid.
    fra = row.get("free_return_arcs") or []
    g_tofs_yr = [a.get("tof_years") for a in fra if a.get("tof_years") is not None]
    if not g_tofs_yr:
        return
    descriptor_tof_days = float(g_tofs_yr[0]) * YEAR_DAYS

    # Per-resonant-leg discrete grids; non-resonant non-transit legs (if any) keep
    # the descriptor value as a single fixed candidate.
    per_leg_grid: dict[int, tuple[float, ...]] = {}
    for i in range(n_legs):
        if i == t_idx:
            continue
        if i in resonant_idxs:
            per_leg_grid[i] = resonant_return_tof_grid(
                seq[i], descriptor_tof_days, n_resonant=n_resonant
            )
        else:
            per_leg_grid[i] = (descriptor_tof_days,)

    # Cartesian product over the non-transit legs (preserves leg ordering).
    other_legs = sorted(per_leg_grid)
    grids = [per_leg_grid[i] for i in other_legs]

    bl_base = base.bounds.lower.copy()
    bu_base = base.bounds.upper.copy()

    def _combos(idx: int, acc: dict[int, float]) -> Iterator[dict[int, float]]:
        if idx == len(other_legs):
            yield dict(acc)
            return
        leg = other_legs[idx]
        for val in grids[idx]:
            acc[leg] = val
            yield from _combos(idx + 1, acc)
        acc.pop(leg, None)

    for combo in _combos(0, {}):
        tof_seed: list[float] = []
        for i in range(n_legs):
            tof_seed.append(transit_tof_days if i == t_idx else combo[i])

        # t0 so the transit leg DEPARTS at the structured epoch.
        t0 = ss.t_depart_sec - sum(tof_seed[:t_idx]) * DAY_S
        eta_seed = tuple(0.5 for _ in range(n_legs))
        x0 = dsm_leg.dsm_chain_decision_vector(
            t0_sec=t0,
            vinf_out0_kms=vinf["E"],
            alpha0=0.0,
            beta0=0.0,
            tof_days_per_leg=tuple(tof_seed),
            eta_per_leg=eta_seed,
            alpha_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
            beta_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
        )

        bl = bl_base.copy()
        bu = bu_base.copy()
        bl[0] = t0 - epoch_halfwidth_days * DAY_S
        bu[0] = t0 + epoch_halfwidth_days * DAY_S
        for i in range(n_legs):
            ui = 4 + i
            bl[ui] = min(bl[ui], 30.0)
            bu[ui] = max(bu[ui], tof_seed[i] * 1.5)
        bounds = dsm_leg.DsmBounds(lower=bl, upper=bu)

        yield MultiArcSeedSpec(
            x0=x0,
            bounds=bounds,
            sequence=seq,
            transit_tof_days=transit_tof_days,
            resonant_tof_days_per_leg={i: combo[i] for i in resonant_idxs},
        )


# ---------------------------------------------------------------------------
# The driver: multi-start over the seed set (design point 4).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ClosureReport:
    """Outcome of a multi-arc closure run over a row (the canonical answer).

    Attributes
    ----------
    row_id:
        The catalogue row id.
    sequence:
        The body sequence closed (empty tuple if the row had no seed).
    best_max_residual_kms:
        The lowest canonical ``max_residual_kms`` seen across all starts (the
        headline number; ``PENALTY_RESIDUAL_KMS`` if every start was out-of-range,
        ``inf`` if no start ran at all).
    converged:
        Whether the best start crossed ``tol_kms`` (the gate decision).
    best_seed:
        The seed spec that achieved ``best_max_residual_kms`` (``None`` if no seed
        was available for the row).
    best_result:
        The full ``DsmChainResult`` of the best start (``None`` on penalty / no seed).
    n_starts_run:
        How many starts were actually solved.
    n_seeds_available:
        How many seeds the enumerator produced (the multi-start budget ceiling).
    gradient, tol_kms:
        The lane and gate the run used (echoed for provenance).
    """

    row_id: str
    sequence: tuple[str, ...]
    best_max_residual_kms: float
    converged: bool
    best_seed: MultiArcSeedSpec | None
    best_result: dsm_leg.DsmChainResult | None
    n_starts_run: int
    n_seeds_available: int
    gradient: str = "lambert"
    tol_kms: float = 0.1
    per_start_residual_kms: tuple[float, ...] = field(default_factory=tuple)


def close_multiarc_row(
    row: dict[str, Any],
    eph: Ephemeris,
    *,
    n_starts: int,
    gradient: str = "lambert",
    tol_kms: float = 0.1,
    t_center_sec: float | None = None,
    n_resonant: int = _DEFAULT_N_RESONANT,
    charge_flyby_continuity: bool = True,
    max_nfev: int = 200,
) -> ClosureReport:
    """Drive the multi-arc seed set through the epoch-safe canonical corrector.

    Enumerates the structured-epoch x discrete-resonant-return seed set
    (:func:`resonant_return_seeds`), solves the first ``n_starts`` of them through
    :func:`safe_chain_residual` (the canonical, epoch-safe metric), and returns the
    best ``(max_residual_kms, converged, best_seed)`` as a :class:`ClosureReport`
    (design point 4). Each start is one bounded local solve from a distinct discrete
    seed — basin sampling over the discrete seed set, which is exactly what the MBH
    wrapper hops within; the coordinator can later widen ``n_starts`` and add MBH
    hops per start for the full campaign.

    Never raises on an out-of-range epoch (the wrapper penalises it). A row with no
    descriptor / no structured seed returns a ``converged=False`` report with an
    empty sequence and ``best_max_residual_kms = inf`` (a clean negative).

    Parameters
    ----------
    row:
        Raw catalogue row dict.
    eph:
        Ephemeris instance (``Ephemeris("astropy")`` for real DE440).
    n_starts:
        Maximum number of distinct seeds to solve (the multi-start budget). The
        enumerator may yield fewer; the report records both.
    gradient:
        Corrector lane (``"lambert"`` default — the convergence reference, #245).
    tol_kms:
        Convergence gate (km/s).
    t_center_sec:
        Centre epoch for the structured transit-leg search. ``None`` (default) uses
        2027-01-01 (the era the #181 ToF-fix closer used; the structured seed is
        insensitive to the exact centre — see the run-1 note).
    n_resonant:
        Number of discrete resonant returns enumerated per resonant leg.
    charge_flyby_continuity, max_nfev:
        Passed through to the corrector (charged mode is the bend-feasible lane).
    """
    row_id = str(row.get("id", "<unknown>"))
    if t_center_sec is None:
        # 2027-01-01 in t_sec (TDB seconds since J2000); ~852055200 s. Computed from
        # the day count to avoid a datetime dependency here.
        t_center_sec = (27.0 * YEAR_DAYS - 0.5) * DAY_S

    seeds = list(resonant_return_seeds(row, eph, t_center_sec, n_resonant=n_resonant))
    n_available = len(seeds)
    if not seeds:
        return ClosureReport(
            row_id=row_id,
            sequence=(),
            best_max_residual_kms=float("inf"),
            converged=False,
            best_seed=None,
            best_result=None,
            n_starts_run=0,
            n_seeds_available=0,
            gradient=gradient,
            tol_kms=tol_kms,
        )

    best_res = float("inf")
    best_conv = False
    best_seed: MultiArcSeedSpec | None = None
    best_result: dsm_leg.DsmChainResult | None = None
    per_start: list[float] = []

    to_run = seeds[: max(0, n_starts)]
    for seed in to_run:
        res_kms, conv, result = safe_chain_residual(
            seed.x0,
            sequence=seed.sequence,
            ephem=eph,
            bounds=seed.bounds,
            tol_kms=tol_kms,
            gradient=gradient,
            charge_flyby_continuity=charge_flyby_continuity,
            max_nfev=max_nfev,
        )
        per_start.append(res_kms)
        # Prefer a converged start; otherwise track the lowest residual seen.
        better = (conv and not best_conv) or (conv == best_conv and res_kms < best_res)
        if best_seed is None or better:
            best_res = res_kms
            best_conv = conv
            best_seed = seed
            best_result = result

    return ClosureReport(
        row_id=row_id,
        sequence=best_seed.sequence if best_seed is not None else (),
        best_max_residual_kms=best_res,
        converged=best_conv,
        best_seed=best_seed,
        best_result=best_result,
        n_starts_run=len(to_run),
        n_seeds_available=n_available,
        gradient=gradient,
        tol_kms=tol_kms,
        per_start_residual_kms=tuple(per_start),
    )
