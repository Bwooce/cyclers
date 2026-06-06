"""M-ED seeding ladder (spec §3): descriptor -> sourced anchor -> coplanar ->
epoch scan.

Resolves the per-leg ``(n_revs, branch)`` topology + asymmetric ToF seed a cell
needs before the N-arc ballistic corrector (search/correct.py) runs. Each rung is
used only when every higher-priority rung is absent for that row, so a row with a
sourced Russell descriptor uses it; a Jones VEM row with no descriptor but sourced
transit ToFs uses the anchor rung; an unanchored row falls back to a coplanar warm
start and finally a bare equispaced epoch-scan seed.

The single-start corrector lands the degenerate high-V_inf basin (the Phase 1
finding; project memory project_s1l1_realeph_closure_blocker.md); family selection
across epochs/branches lives in the scan rung's epoch window (spec §3.4), the same
mechanism the prototype's main() loop uses (scripts/correct_s1l1_twoarc.py).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.descriptor import parse_free_return_arcs
from cyclerfinder.search.sequence import Cell

SeedSource = Literal["descriptor", "anchor", "coplanar", "scan"]


@dataclass(frozen=True)
class SeedPlan:
    """A resolved seed for one cell: per-leg topology + ToF seed + the rung used.

    ``tof_seed_days`` is the full per-leg seed (length ``len(cell.sequence) - 1``)
    EXCEPT for the descriptor rung, where it carries one seed per Earth-Earth
    return arc (the corrector reconstructs the transfer-leg ToFs from the period).
    ``source`` records which ladder rung produced the plan.
    """

    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    tof_seed_days: tuple[float, ...]
    source: SeedSource


def _equispaced_seed(cell: Cell) -> tuple[float, ...]:
    """Bare equispaced per-leg seed (days) over the cell's resolved period."""
    from cyclerfinder.search.optimize import (
        _ephemeris_tof_seed_and_bounds,
        _target_period_sec,
    )

    target_period_sec = _target_period_sec(cell)
    seed_days, _bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)
    return tuple(float(s) for s in seed_days)


def _coplanar_seed(cell: Cell, ephem: Ephemeris) -> tuple[float, ...] | None:
    """Per-leg ToF seed (days) from a coplanar warm start (spec §3.3).

    Runs the circular-coplanar idealized optimiser and reads its converged leg
    ToFs (consecutive encounter-epoch diffs). Returns ``None`` if the warm start
    does not converge so the ladder can fall through to the scan rung. The leg
    ToFs are a *seed only* (refined by the corrector), so this is not a golden-
    discipline violation: no self-computed value is asserted as an EXPECTED.
    """
    from cyclerfinder.search.optimize import optimise_cell_idealized

    try:
        result = optimise_cell_idealized(cell, ephem, vinf_cap=99.0, seed=0)
    except (ValueError, RuntimeError):
        return None
    if not result.converged:
        return None
    encs = result.best_cycler.encounters
    if len(encs) < 2:
        return None
    tofs = [(encs[i + 1].t - encs[i].t) / 86400.0 for i in range(len(encs) - 1)]
    if any(t <= 0 for t in tofs):
        return None
    return tuple(float(t) for t in tofs)


def resolve_seed(
    cell: Cell,
    *,
    free_return_arcs: Sequence[Mapping[str, Any]] | None = None,
    anchor_tofs: Sequence[float] | None = None,
    anchor_vinf: Mapping[str, float] | None = None,
    coplanar_tofs: Sequence[float] | None = None,
    ephem: Ephemeris | None = None,
) -> SeedPlan:
    """Resolve a cell's seed via the ladder (spec §3): descriptor -> anchor ->
    coplanar -> scan. Each rung fires only when every higher rung is absent.

    Rung 1 (descriptor): a Russell ``free_return_arcs[]`` list -> per-leg topology
    + per-arc ToF seeds via :func:`parse_free_return_arcs`.
    Rung 2 (anchor): sourced transit ToFs (``anchor_tofs``) + per-body V_inf
    targets (``anchor_vinf``) -> the asymmetric seed; topology from the cell.
    Rung 3 (coplanar): a coplanar warm start (explicit ``coplanar_tofs`` or, when
    an ``ephem`` is supplied, the idealized optimiser's converged leg ToFs).
    Rung 4 (scan, last resort): an equispaced seed; family selection across epochs
    is the corrector's ±10 yr scan window (spec §3.4).
    """
    revs = tuple(cell.per_leg_revs)
    branches = tuple(cell.per_leg_branch)

    # Rung 1 — descriptor.
    if free_return_arcs:
        d_revs, d_branches, d_seeds = parse_free_return_arcs(free_return_arcs)
        return SeedPlan(
            per_leg_revs=d_revs,
            per_leg_branch=d_branches,
            tof_seed_days=d_seeds,
            source="descriptor",
        )

    # Rung 2 — sourced anchor.
    if anchor_tofs is not None and anchor_vinf is not None:
        return SeedPlan(
            per_leg_revs=revs,
            per_leg_branch=branches,
            tof_seed_days=tuple(float(t) for t in anchor_tofs),
            source="anchor",
        )

    # Rung 3 — coplanar warm start.
    if coplanar_tofs is not None:
        return SeedPlan(
            per_leg_revs=revs,
            per_leg_branch=branches,
            tof_seed_days=tuple(float(t) for t in coplanar_tofs),
            source="coplanar",
        )
    if ephem is not None:
        warm = _coplanar_seed(cell, ephem)
        if warm is not None:
            return SeedPlan(
                per_leg_revs=revs,
                per_leg_branch=branches,
                tof_seed_days=warm,
                source="coplanar",
            )

    # Rung 4 — epoch scan (last resort).
    return SeedPlan(
        per_leg_revs=revs,
        per_leg_branch=branches,
        tof_seed_days=_equispaced_seed(cell),
        source="scan",
    )
