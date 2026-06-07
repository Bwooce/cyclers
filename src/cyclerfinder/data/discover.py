"""Ledger-backed discovery runner (M7) — spec §13.6 / §13.8 / §14.

``discover`` is the M7 composition layer: it walks the M4 enumerator's
feasible cells, runs the M5 idealised optimiser on each, reduces the
result to a canonical signature, matches it against the seed catalogue,
auto-validates it through the cheapest-first V0→V2 gauntlet, and records
the outcome to the JSONL ledger. It yields
``(OptimisationResult, MatchResult, level)`` per solved cell so a
downstream reporter (M8) can present finds with their match outcome and
validation level.

Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md`` §3.5.

Ledger interaction
------------------
The plan pseudocode sketches a ``claim`` / ``record`` two-phase write.
The shipped :class:`~cyclerfinder.data.ledger.Ledger` is append-only with
first-write-wins on reload, and :meth:`Ledger.record` refuses to
double-record a ``cell_id``. So the runner uses a single terminal
``record`` per cell guarded by :meth:`Ledger.has`; a concurrent writer
that loses the race surfaces as :class:`LedgerError`, which the runner
swallows (first-writer-wins) and skips the cell.
"""

from __future__ import annotations

import contextlib
import itertools
import socket
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import (
    DAYS_PER_JULIAN_YEAR,
    MU_SUN_KM3_S2,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import (
    Catalog,
    MatchResult,
    canonical_signature,
    load_catalog,
    match,
)
from cyclerfinder.data.ledger import Ledger, LedgerEntry, LedgerError
from cyclerfinder.search.optimize import (
    OptimisationResult,
    optimise_cell_ephemeris,
    optimise_cell_idealized,
)
from cyclerfinder.search.sequence import Cell, feasible_cells
from cyclerfinder.verify.crosscheck import crosscheck_cycler
from cyclerfinder.verify.propagate import verify_long_term_stability

_V2_N_LAPS: int = 3
"""Laps propagated by the V2 auto-gate (matches the M6a 3-lap gate)."""

_DEFAULT_CATALOGUE_PATH: Path = Path(__file__).resolve().parents[3] / "data" / "catalogue.yaml"
"""Repo-root ``data/catalogue.yaml`` — the default source for the free-return
descriptor discover path (task #106)."""


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _auto_validate(
    result: OptimisationResult,
    catalog: Catalog,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    enable_v3: bool = False,
    priority_date_iso: str | None = None,
    vinf_targets_kms: dict[str, float] | None = None,
) -> str:
    """Run the cheapest-first V0→V2(/V3) gauntlet; return the highest level passed.

    V0 is the optimiser's ``constraints_satisfied`` (the caller only
    reaches here when it is ``True``). V1 is the lamberthub
    cross-check; V2 is M6a's multi-lap closure. The gauntlet stops at
    the first failing gate, so ``level`` is the highest *consecutively*
    passed gate per spec §14.

    V1 / V2 are wrapped defensively: a solver error (e.g. a degenerate
    leg the cross-check cannot re-solve) degrades the level rather than
    aborting the whole run.
    """
    del catalog  # reserved for future catalogue-aware gates
    cycler = result.best_cycler
    level = "V0"

    try:
        crosscheck = crosscheck_cycler(cycler, ephem)
    except Exception:
        return level
    if not (crosscheck and all(r.passed for r in crosscheck)):
        return level
    level = "V1"

    try:
        report = verify_long_term_stability(cycler, _V2_N_LAPS, ephem)
    except Exception:
        return level
    if not report.stable:
        return level
    level = "V2"

    if enable_v3:
        # V3 is the M-ED ballistic-closure gate (spec §0 finding 3, task #109):
        # run optimise_cell_ephemeris in ballistic mode and promote to V3 only
        # when the N-arc V_inf-continuity closure is genuinely satisfied
        # (constraints_satisfied = converged AND bend-feasible AND V_inf-cap),
        # NOT on an exception (the prior branch assumed the M5 stub still raised
        # NotImplementedError; it no longer does). A solver error degrades the
        # level rather than aborting the run, consistent with V1/V2.
        from cyclerfinder.search.optimize import optimise_cell_ephemeris

        try:
            v3_result = optimise_cell_ephemeris(
                result.cell,
                ephem,
                vinf_cap=vinf_cap,
                mode="ballistic",
                priority_date_iso=priority_date_iso,
                vinf_targets_kms=vinf_targets_kms,
            )
        except Exception:
            return level
        if v3_result.constraints_satisfied:
            level = "V3"

    return level


def discover(
    bodies: tuple[str, ...],
    k_synodic: int,
    vinf_cap: float,
    ledger_path: str,
    *,
    ephem: Ephemeris | None = None,
    l_max: int = 4,
    n_max: int = 0,
    branch_set: tuple[str, ...] = ("single",),
    max_cells: int | None = None,
    n_starts: int = 5,
    seed: int = 0,
    use_de: bool = True,
    rp_factors: dict[str, float] | None = None,
    optimiser: str = "idealized",
    priority_date_iso: str | None = None,
    vinf_targets_kms: dict[str, float] | None = None,
    host: str | None = None,
    enable_v3: bool = False,
    finder_version: str = "0.7.0",
) -> Iterator[tuple[OptimisationResult, MatchResult, str]]:
    """Run the ledger-backed discovery pipeline over the feasible cells.

    Yields ``(result, match_result, level)`` for every *solved* cell
    (one whose optimiser converged with constraints satisfied).
    Non-solved cells (``searched`` / ``failed``) are recorded in the
    ledger but not yielded. Cells already present in the ledger are
    skipped, so an interrupted run resumes without repeating work.

    ``n_max`` / ``branch_set`` widen the enumeration to multi-revolution
    legs (threaded straight into :func:`feasible_cells`); ``max_cells``
    bounds the feasible-cell stream so a sweep can be capped for cost. The
    defaults (``n_max=0``, ``branch_set=("single",)``, ``max_cells=None``)
    keep the single-rev behaviour byte-identical for existing callers.

    ``optimiser`` selects the per-cell solver: ``"idealized"`` (default,
    circular-coplanar :func:`optimise_cell_idealized`) or ``"ephemeris"``
    (real-DE440 :func:`optimise_cell_ephemeris`). Ephemeris mode phase-
    matches a launch epoch from ``priority_date_iso`` + ``vinf_targets_kms``
    (both required to resolve a real window — without them every cell
    returns a non-converged result and is recorded but not yielded), and
    signs results into the ``"analytic-ephemeris"`` signature pool rather
    than ``"circular-coplanar"``.
    """
    del finder_version  # provenance is the operator-driven writeback's concern
    ephem = ephem or Ephemeris(model="circular")
    ledger = Ledger(ledger_path)
    catalog = load_catalog()
    host = host or socket.gethostname()

    cell_stream: Iterator[Cell] = feasible_cells(
        bodies,
        l_max=l_max,
        k_max=k_synodic,
        n_max=n_max,
        vinf_cap=vinf_cap,
        ephem=ephem,
        branch_set=branch_set,
    )
    if max_cells is not None:
        cell_stream = itertools.islice(cell_stream, max_cells)

    for cell in cell_stream:
        if ledger.has(cell.id):
            continue

        try:
            if optimiser == "ephemeris":
                # Real-ephemeris (DE440) mode: phase-matches a launch epoch from
                # the supplied priority date + V-infinity targets, then minimises
                # maintenance dV. Without targets it returns a non-converged
                # result (recorded as "searched"), so a blind sweep is a no-op.
                result = optimise_cell_ephemeris(
                    cell,
                    ephem,
                    vinf_cap=vinf_cap,
                    priority_date_iso=priority_date_iso,
                    vinf_targets_kms=vinf_targets_kms,
                    n_starts=n_starts,
                    seed=seed,
                    rp_factors=rp_factors,
                )
            else:
                result = optimise_cell_idealized(
                    cell,
                    ephem,
                    vinf_cap=vinf_cap,
                    n_starts=n_starts,
                    seed=seed,
                    use_de=use_de,
                    rp_factors=rp_factors,
                )
        except (ValueError, RuntimeError):
            _record(
                ledger,
                LedgerEntry(cell.id, "failed", 0, None, (), None, _now(), host),
            )
            continue

        if not result.constraints_satisfied:
            _record(
                ledger,
                LedgerEntry(cell.id, "searched", 0, None, (), "V0", _now(), host),
            )
            continue

        model_assumption = "analytic-ephemeris" if optimiser == "ephemeris" else "circular-coplanar"
        signature = canonical_signature(result.best_cycler, model_assumption=model_assumption)
        match_result = match(signature, catalog)
        level = _auto_validate(
            result,
            catalog,
            ephem,
            vinf_cap=vinf_cap,
            enable_v3=enable_v3,
            priority_date_iso=priority_date_iso,
            vinf_targets_kms=vinf_targets_kms,
        )
        _record(
            ledger,
            LedgerEntry(
                cell.id,
                "solved",
                1,
                result.closure_residual_kms,
                (signature.hash,),
                level,
                _now(),
                host,
            ),
        )
        yield result, match_result, level


def _record(ledger: Ledger, entry: LedgerEntry) -> None:
    """Record ``entry``, swallowing a lost-race double-record (first wins)."""
    with contextlib.suppress(LedgerError):
        ledger.record(entry)


# ---------------------------------------------------------------------------
# Free-return descriptor discover path (task #106, the SnLm re-scope).
#
# The MULTI_ENCOUNTER_SEQUENCE catalogue rows that carry a ``free_return_arcs[]``
# descriptor (the McConaghy/Russell-12 SnLm chains) are NOT reachable through the
# 2-encounter idealised-optimiser sweep above. They ARE reachable through the
# #137 free-return (radial-crossing) genome: seed (a, e) from the SOURCED aphelion
# + transit, close the single ellipse, and the per-body V_inf EMERGES for
# non-circular comparison against the sourced anchor.
#
# This is a real path that REUSES the existing physics
# (cyclerfinder.search.free_return + .free_return_v1) — it does NOT reimplement
# the corrector. It is a sibling of discover(): same role (walk catalogue-eligible
# structures, attempt closure, classify), different genome. Like-for-like circular
# closure of a circular-coplanar source — NOT a real-ephemeris V3 result.
# ---------------------------------------------------------------------------

_DAY_S = SECONDS_PER_DAY

# Free-return match tolerance vs the row's SOURCED V_inf anchor (km/s). Mirrors
# the campaign default (scripts/campaign_russell12.py:TOL_VINF_KMS) and the §14
# V_inf-continuity ceiling — NOT loosened per row (golden discipline).
_FR_VINF_TOL_KMS = 0.5
_FR_TOL_KMS = 0.1  # corrector residual closure floor
# Dense t0 phase floor: deep-aphelion high-e rows have a narrow residual basin a
# coarse grid steps over (#137 Part 3). Pure residual evals (no Lambert), cheap.
_FR_PHASE_EPOCHS = 4096


@dataclass(frozen=True)
class FreeReturnDiscovery:
    """One descriptor row's free-return discover outcome (task #106).

    Attributes
    ----------
    row_id:
        Catalogue id of the descriptor-bearing SnLm row.
    outcome:
        ``"CLOSE-AND-MATCH"`` (closed AND emerged V_inf within tolerance of the
        sourced anchor, symmetric single-ellipse row),
        ``"CLOSE-MATCH-SYMMETRIC-ONLY"`` (V_inf matches but the row's transit legs
        are asymmetric, so a single symmetric ellipse only reproduces it
        partially — reported honestly, not over-claimed),
        ``"CLOSE-OFF-ANCHOR"`` (closed but off the sourced anchor), ``"NO-CLOSE"``
        (corrector did not reach the residual floor), or ``"NO-SEED"`` (missing
        sourced aphelion/transit to seed from).
    closed:
        The corrector reached the residual floor.
    vinf_match:
        Every sourced per-body V_inf matched the emerged V_inf within tolerance.
    v1_passed:
        §14 V1 mechanics (paths a+c) + V_inf-continuity passed on the closed,
        reconstructed E->M->E arc (``None`` when not run / not closed).
    max_residual_kms:
        Corrector closure residual (km/s).
    derived_vinf_kms:
        Per-body EMERGED V_inf (evidence; never imposed).
    sourced_vinf_kms:
        Per-body SOURCED V_inf anchor (the EXPECTED side — golden).
    """

    row_id: str
    outcome: str
    closed: bool
    vinf_match: bool
    v1_passed: bool | None
    max_residual_kms: float
    derived_vinf_kms: dict[str, float]
    sourced_vinf_kms: dict[str, float]


def _descriptor_rows(catalogue_path: Path) -> list[dict[str, Any]]:
    """Catalogue rows carrying a ``free_return_arcs[]`` descriptor (sorted by id).

    N-agnostic: driven by descriptor presence, not a hard-coded id list — exactly
    the DESCRIPTOR_CLOSABLE sub-classification (tests/_catalogue_loader.py).
    """
    rows = yaml.safe_load(catalogue_path.read_text())
    out = [r for r in rows if r.get("free_return_arcs")]
    out.sort(key=lambda r: r["id"])
    return out


def _best_phase_t0(
    a_seed: float, e_seed: float, period_sec: float, ephem: Ephemeris, mu: float
) -> tuple[float, float]:
    """Scan t0 over one period; return ``(best_t0_sec, best_residual_kms)`` at the
    SOURCED ``(a, e)``. Mirrors the campaign/probe best-phase selection."""
    from cyclerfinder.search.free_return import _residuals as _fr_residuals

    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, _FR_PHASE_EPOCHS, endpoint=False):
        t0 = float(frac) * period_sec
        try:
            res = _fr_residuals(
                np.array([a_seed, e_seed, t0]),
                period_days=period_sec / _DAY_S,
                ephem=ephem,
                bodies=("E", "M"),
                mu=mu,
            )
        except Exception:
            continue
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0, best_res


def discover_free_return(
    *,
    ephem: Ephemeris | None = None,
    catalogue_path: str | Path | None = None,
    run_v1: bool = True,
    mu: float = MU_SUN_KM3_S2,
) -> Iterator[FreeReturnDiscovery]:
    """Attempt free-return closure for every descriptor-bearing SnLm row.

    The sibling discover path for the DESCRIPTOR_CLOSABLE bucket (task #106): for
    each catalogue row carrying a ``free_return_arcs[]`` descriptor it derives the
    seed ``(a, e)`` from the SOURCED aphelion + transit
    (:func:`~cyclerfinder.search.free_return.seed_ae_from_aphelion_transit`), scans
    t0 for the best phase, closes the single radial-crossing ellipse
    (:func:`~cyclerfinder.search.free_return.free_return_correct`), and compares the
    EMERGED per-body V_inf against the row's SOURCED anchor — yielding the row OUT
    of the "unreachable" bucket with an honest outcome.

    Reuses the #137 physics verbatim (the corrector + the §14 V1 mechanics); it
    does NOT reimplement closure. Closure is circular-coplanar like-for-like
    (reproducing a circular-coplanar source), NOT a real-ephemeris V3 result.

    GOLDEN DISCIPLINE: the EXPECTED side of every match is the row's SOURCED V_inf
    anchor; the emerged V_inf (evidence) is never imposed.
    """
    from cyclerfinder.search.free_return import (
        free_return_correct,
        seed_ae_from_aphelion_transit,
    )

    ephem = ephem or Ephemeris(model="circular")
    cat_path = Path(catalogue_path) if catalogue_path else _DEFAULT_CATALOGUE_PATH

    for row in _descriptor_rows(cat_path):
        rid = row["id"]
        aphelion = (row.get("orbit_elements") or {}).get("aphelion_au")
        transit = (row.get("invariants") or {}).get("transit_times_days")
        sourced_vinf = {
            e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])
        }
        if aphelion is None or not transit:
            yield FreeReturnDiscovery(
                row_id=rid,
                outcome="NO-SEED",
                closed=False,
                vinf_match=False,
                v1_passed=None,
                max_residual_kms=float("inf"),
                derived_vinf_kms={},
                sourced_vinf_kms=sourced_vinf,
            )
            continue

        period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * _DAY_S
        # A single symmetric ellipse cannot represent a different-per-leg transit;
        # the campaign flags these CLOSE-MATCH-SYMMETRIC-ONLY (honest, not
        # over-claimed) even when the emerged V_inf matches.
        asymmetric = len(transit) >= 2 and abs(float(transit[0]) - float(transit[1])) > 1.0
        a_seed, e_seed = seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]), mu=mu)
        best_t0, _ = _best_phase_t0(a_seed, e_seed, period_sec, ephem, mu)
        sol = free_return_correct(
            t0_seed_sec=best_t0,
            a_seed_au=a_seed,
            e_seed=e_seed,
            period_sec=period_sec,
            ephem=ephem,
            mu=mu,
            tol_kms=_FR_TOL_KMS,
        )

        derived = {k: float(v) for k, v in sol.vinf_kms.items()}
        vinf_match = bool(sourced_vinf) and all(
            (body in derived) and abs(derived[body] - src) <= _FR_VINF_TOL_KMS
            for body, src in sourced_vinf.items()
        )

        v1_passed: bool | None = None
        if run_v1 and sol.converged:
            from cyclerfinder.search.free_return_v1 import free_return_v1_mechanics

            v1_passed = free_return_v1_mechanics(sol, ephem, period_sec, mu=mu).v1_passed

        if not sol.converged:
            outcome = "NO-CLOSE"
        elif vinf_match and asymmetric:
            outcome = "CLOSE-MATCH-SYMMETRIC-ONLY"
        elif vinf_match:
            outcome = "CLOSE-AND-MATCH"
        else:
            outcome = "CLOSE-OFF-ANCHOR"

        yield FreeReturnDiscovery(
            row_id=rid,
            outcome=outcome,
            closed=bool(sol.converged),
            vinf_match=vinf_match,
            v1_passed=v1_passed,
            max_residual_kms=float(sol.max_residual_kms),
            derived_vinf_kms=derived,
            sourced_vinf_kms=sourced_vinf,
        )


__all__ = ["FreeReturnDiscovery", "discover", "discover_free_return"]
