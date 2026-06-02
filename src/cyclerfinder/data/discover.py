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
from datetime import UTC, datetime

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import (
    Catalog,
    MatchResult,
    canonical_signature,
    load_catalog,
    match,
)
from cyclerfinder.data.ledger import Ledger, LedgerEntry, LedgerError
from cyclerfinder.search.optimize import OptimisationResult, optimise_cell_idealized
from cyclerfinder.search.sequence import Cell, feasible_cells
from cyclerfinder.verify.crosscheck import crosscheck_cycler
from cyclerfinder.verify.propagate import verify_long_term_stability

_V2_N_LAPS: int = 3
"""Laps propagated by the V2 auto-gate (matches the M6a 3-lap gate)."""


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _auto_validate(
    result: OptimisationResult,
    catalog: Catalog,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    enable_v3: bool = False,
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
        # V3 (ephemeris-mode TCM) needs M6b's optimise_cell_ephemeris,
        # which is an M5 stub raising NotImplementedError. Guarded off
        # by default; M6b flips the flag once it lands.
        from cyclerfinder.search.optimize import optimise_cell_ephemeris

        optimise_cell_ephemeris(result.cell, ephem, vinf_cap=vinf_cap)  # raises until M6b lands
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

        signature = canonical_signature(result.best_cycler, model_assumption="circular-coplanar")
        match_result = match(signature, catalog)
        level = _auto_validate(result, catalog, ephem, vinf_cap=vinf_cap, enable_v3=enable_v3)
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


__all__ = ["discover"]
