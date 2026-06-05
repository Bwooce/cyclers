"""Structural cell enumeration and Tisserand pruning (M4).

Spec references
---------------
* §13.1 — the structural ``cell`` definition (the atomic unit of search).
* §13.2 — iterative-deepening frontier under ``(L_max, k_max, N_max)`` caps.
* §13.3 — Tisserand pruning predicate (the M4 binding gate).
* §13.8 — deterministic, sortable ``Cell.id`` format and the deepening
  loop sketch.

Design split (binding M4 boundary)
----------------------------------
M4 owns the **discrete combinatorial layer**: enumerate cells, prune them
by Tisserand energetic feasibility, and yield the survivors. M4 does
**not** build cyclers (M3 already does, against an explicit schedule),
does **not** search encounter times within a cell (M5), and does **not**
persist a work queue / ledger (M7). The ``deepening_frontier`` generator
ships in M4 as a single-process, in-memory-deduped iterator suitable for
M5 development; M7 replaces (not extends) it with a ledger-backed work
queue.

Plan: ``docs/phases/m4-enumeration-scoring/plan.md`` §3.1.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterator
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from cyclerfinder.search.tisserand import linkable, linkable_3d

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris


# ---------------------------------------------------------------------------
# Branch alphabet (spec §13.8)
# ---------------------------------------------------------------------------

_BRANCH_LETTER: dict[str, str] = {"single": "s", "low": "l", "high": "h"}
"""Letters used by :attr:`Cell.id` per spec §13.8. ``single`` ≡ direct
(0-rev) leg; ``low``/``high`` are the two multi-rev Lambert branches."""

_MULTI_REV_BRANCHES: frozenset[str] = frozenset({"low", "high"})
"""Branch labels valid for a leg with ``n_revs >= 1``."""


# ---------------------------------------------------------------------------
# Cell — the atomic unit of search (spec §13.1, §13.8)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Cell:
    """Atomic unit of search per spec §13.1.

    Discrete structural specification of a cycler candidate. Continuous
    DOF (encounter epochs, phases) are bounded sub-problems *inside* a
    cell, optimised in M5 — not held on the dataclass.

    Attributes
    ----------
    bodies:
        Canonical body set (e.g. ``("V","E","M")``). The caller supplies
        the canonical order; M4 does not sort. Stored verbatim into the
        ``Cell.id`` prefix.
    sequence:
        Ordered flyby sequence (e.g. ``("E","V","M","E","M")``).
        Adjacent entries must differ; the enumerator enforces this at
        generation time.
    period_k:
        Period in synodic multiples, ``k >= 1``. Combined with the
        per-pair synodic period this gives the cell's target heliocentric
        period; the actual target is M5's responsibility to compute.
    per_leg_revs:
        Heliocentric revolutions per leg, length ``len(sequence) - 1``.
        ``0`` denotes a direct leg (single Lambert branch); ``>= 1``
        requires the leg's branch to be ``"low"`` or ``"high"``.
    per_leg_branch:
        Per-leg Lambert branch in ``{"single","low","high"}``, length
        ``len(sequence) - 1``. The enumerator enforces the
        ``revs == 0 -> branch == "single"`` and
        ``revs >= 1 -> branch in {"low","high"}`` invariant.
    period_basis:
        Optional catalogue anchor pair (e.g. ``("E","M")``) for a
        ≥3-body cell, so the cell echoes its YAML row verbatim while
        ``period_k`` stays the sourced value (no silent rewrite, plan
        §2). ``None`` for 2-body cells — the native single-pair period
        path needs no explicit basis. When set, ``Cell.id`` gains a
        ``|p<AB>`` token.

    Notes
    -----
    * No canonicalisation of ``sequence`` happens at the cell layer; spec
      §16.2's lexicographically-minimal rotation is M7's catalogue
      concern, kept out of M4 deliberately so the enumerator does not
      yield the same structure under multiple rotations.
    * Tuples are hashable, so ``Cell`` is hashable and usable as a
      ``dict`` key — needed by M7's ledger; M4 pays nothing to provide it.
    """

    bodies: tuple[str, ...]
    sequence: tuple[str, ...]
    period_k: int
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    period_basis: tuple[str, str] | None = None

    @property
    def id(self) -> str:
        """Deterministic sortable identifier per spec §13.8.

        Format: ``{bodyset}|{sequence}|k{K}|r{revs}|b{branches}``.

        The branch alphabet is ``{"single":"s", "low":"l", "high":"h"}``.

        Examples
        --------
        ``VEM|E-V-M-E-M-E|k3|r00101|blllll`` — the spec §13.8 worked
        example (all multi-rev legs taking the ``low`` branch).

        ``EM|E-M-E|k2|r00|bss`` — a 2-synodic E-M-E cell with two direct
        legs (the M4 native case).
        """
        bodyset = "".join(self.bodies)
        sequence = "-".join(self.sequence)
        revs = "".join(str(r) for r in self.per_leg_revs)
        branches = "".join(_BRANCH_LETTER[b] for b in self.per_leg_branch)
        base = f"{bodyset}|{sequence}|k{self.period_k}|r{revs}|b{branches}"
        if self.period_basis is not None:
            base += f"|p{''.join(self.period_basis)}"
        return base


# ---------------------------------------------------------------------------
# enumerate_cells — the unpruned combinatorial generator (spec §13.2)
# ---------------------------------------------------------------------------


def _branch_choices_for(revs: int, branch_set: tuple[str, ...]) -> tuple[str, ...]:
    """Branch labels valid for a leg with ``revs`` heliocentric revolutions.

    ``revs == 0`` accepts only ``"single"`` (the only Lambert branch that
    exists for a direct leg). ``revs >= 1`` accepts whatever of
    ``{"low","high"}`` the caller exposed in ``branch_set``.
    """
    if revs == 0:
        return ("single",) if "single" in branch_set else ()
    return tuple(b for b in branch_set if b in _MULTI_REV_BRANCHES)


def enumerate_cells(
    body_set: tuple[str, ...],
    l_max: int,
    k_max: int,
    n_max: int,
    branch_set: tuple[str, ...] = ("single",),
) -> Iterator[Cell]:
    """Yield every combinatorial cell under the given caps. No pruning.

    Per spec §13.2 the cell set is finite under fixed caps; this
    generator walks it in lexicographically-stable order:

    .. code-block:: text

        for L in range(2, l_max + 1):              # min length 2 = one leg
          for k in range(1, k_max + 1):
            for sequence in body_set^L:            # adjacency-distinct
              for per_leg_revs in {0..n_max}^(L-1):
                for per_leg_branch in branch_set^(L-1):
                  yield Cell(...)

    The enumerator enforces two invariants at generation time, rather
    than as post-filters, to keep the yielded count honest:

    1. **Adjacency distinctness** — ``sequence[i] != sequence[i+1]``.
       A flyby of body X immediately followed by another flyby of X with
       no intervening heliocentric arc is degenerate.
    2. **Rev / branch consistency** — a 0-rev leg pins ``branch="single"``;
       a ``>=1``-rev leg cycles through ``branch_set ∩ {"low","high"}``.
       If ``branch_set`` does not include ``"low"``/``"high"`` but
       ``n_max >= 1``, the multi-rev rev tuples yield nothing —
       documented behaviour, not a bug.

    Parameters
    ----------
    body_set:
        Canonical body codes the sequence can visit. ``Cell.bodies`` is
        set to this tuple verbatim.
    l_max:
        Maximum encounters in the sequence (``>= 2``).
    k_max:
        Maximum period in synodic multiples (``>= 1``).
    n_max:
        Maximum heliocentric revolutions per leg (``>= 0``). ``n_max == 0``
        restricts to direct legs (the M1/M3 only-supported Lambert
        regime); ``n_max >= 1`` requires the caller to also expand
        ``branch_set`` to include ``"low"``/``"high"`` for any cells to
        be yielded at the multi-rev rev tuples.
    branch_set:
        Allowed branches per leg. Default ``("single",)`` matches the M1
        Lambert solver's current single-branch capability.

    Yields
    ------
    Cell
        Every combinatorial cell under the caps. Caller responsible for
        pruning (e.g. :func:`tisserand_feasible`).

    Notes
    -----
    The yielded count is bounded above by
    ``sum_{L=2..l_max} k_max · |body_set| · (|body_set|-1)^(L-1)
    · (effective branch product)``. For
    ``body_set=("E","M"), l_max=4, k_max=2, n_max=0,
    branch_set=("single",)`` the count is exactly 12 (plan §4.4).

    The generator is **lazy** (``Iterator[Cell]``, not ``list[Cell]``):
    at higher caps with multi-rev branches the count explodes
    combinatorially, so callers should consume the iterator with
    ``itertools.islice`` rather than materialising it.
    """
    if l_max < 2:
        return
    if k_max < 1:
        return
    if n_max < 0:
        return
    bodies = tuple(body_set)
    for length in range(2, l_max + 1):
        # Adjacency-distinct sequences of the requested length.
        for seq in itertools.product(bodies, repeat=length):
            if any(seq[i] == seq[i + 1] for i in range(length - 1)):
                continue
            n_legs = length - 1
            for k in range(1, k_max + 1):
                # Per-leg revs.
                for revs in itertools.product(range(n_max + 1), repeat=n_legs):
                    # For each leg, choose the valid branches given its revs.
                    per_leg_choices = [_branch_choices_for(r, branch_set) for r in revs]
                    if any(len(choices) == 0 for choices in per_leg_choices):
                        # No valid branch for some leg under the supplied
                        # branch_set; documented silent skip.
                        continue
                    for branch_combo in itertools.product(*per_leg_choices):
                        yield Cell(
                            bodies=bodies,
                            sequence=tuple(seq),
                            period_k=k,
                            per_leg_revs=tuple(revs),
                            per_leg_branch=tuple(branch_combo),
                        )


# ---------------------------------------------------------------------------
# tisserand_feasible — the M4 pruning gate (spec §13.3)
# ---------------------------------------------------------------------------

_TISSERAND_VINF_SAMPLES: int = 24
"""Number of V∞ samples in ``(0.5, vinf_cap]`` used by
:func:`tisserand_feasible`. Coarser than
:func:`cyclerfinder.search.tisserand.linkable_region` would give but
sufficient for the gate: a real linkable pair has a *band* of linkable
V∞, not a discrete value. Tests pin this resolution."""

_TISSERAND_VINF_FLOOR_KMS: float = 0.5
"""Lower bound of the V∞ sampling grid (km/s). Avoids the always-false
asymptote near V∞ → 0 where every pair's contour collapses to its own
body's circular orbit."""


def tisserand_feasible(
    cell: Cell,
    vinf_cap: float,
    ephem: Ephemeris | None = None,
) -> bool:
    """Tisserand pruning per spec §13.3 — the M4 binding gate.

    Returns ``True`` iff for every consecutive body pair
    ``(cell.sequence[i], cell.sequence[i+1])`` there exists a V∞ in
    ``(0.5, vinf_cap]`` at which the two bodies are
    :func:`~cyclerfinder.search.tisserand.linkable` — i.e. their
    constant-V∞ contours intersect in ``(a, e)`` space, so a spacecraft
    orbit of fixed ``(a, e)`` is reachable from both bodies at that V∞.

    The vast majority of cells die here; compute on the survivors is
    spent only on energetically viable structures (spec §13.3).

    Parameters
    ----------
    cell:
        The cell to test.
    vinf_cap:
        Common V∞ ceiling, km/s. Each consecutive pair must be linkable
        somewhere in ``(0.5, vinf_cap]``.
    ephem:
        When ``None`` (the default) the coplanar :func:`~cyclerfinder.search
        .tisserand.linkable` predicate is used — byte-identical to the M4
        behaviour. When an :class:`~cyclerfinder.core.ephemeris.Ephemeris` is
        supplied the 3-D :func:`~cyclerfinder.search.tisserand.linkable_3d`
        predicate is consulted for each pair, with a coplanar-``True``
        short-circuit: a pair that is already coplanar-linkable needs no 3-D
        scan (coplanar ``True`` ⇒ 3-D ``True``), so the expensive 3-D scan
        runs only on the pairs the coplanar predicate rejects. The ephemeris
        object itself is not sampled — planet inclinations/eccentricities are
        epoch-independent mean elements read from :data:`PLANETS`.

    Returns
    -------
    bool
        ``True`` iff every consecutive pair has at least one linkable
        V∞ in the sampling band. **Never raises** — failures (NaN,
        no-bracket, unknown body) return ``False``, mirroring the
        :func:`~cyclerfinder.search.tisserand.linkable` contract so the
        enumerator can call this in a tight loop without try/except
        scaffolding.

    Notes
    -----
    Sampling resolution is fixed at :data:`_TISSERAND_VINF_SAMPLES`
    points evenly spaced in ``(0.5, vinf_cap]``. With ``ephem=None`` the
    coplanar (i=0) :func:`linkable` predicate is used; with an ephemeris the
    3-D :func:`linkable_3d` predicate extends it (see the ``ephem`` note).
    """
    use_3d = ephem is not None
    try:
        if vinf_cap <= _TISSERAND_VINF_FLOOR_KMS:
            return False
        if len(cell.sequence) < 2:
            return False
        grid = np.linspace(_TISSERAND_VINF_FLOOR_KMS, vinf_cap, _TISSERAND_VINF_SAMPLES)
        for i in range(len(cell.sequence) - 1):
            body_a = cell.sequence[i]
            body_b = cell.sequence[i + 1]
            found = False
            for vinf in grid:
                v = float(vinf)
                # Coplanar predicate first: a coplanar-True is also 3-D-True,
                # so it short-circuits the expensive 3-D scan. Only the
                # coplanar-False samples fall through to linkable_3d (and only
                # when an ephemeris was supplied).
                if linkable(body_a, body_b, v):
                    found = True
                    break
                if use_3d and linkable_3d(body_a, body_b, v):
                    found = True
                    break
            if not found:
                return False
        return True
    except Exception:
        # Mirror tisserand.linkable's never-raise contract — the
        # enumerator must be able to call this from a tight loop without
        # try/except scaffolding.
        return False


# ---------------------------------------------------------------------------
# feasible_cells — convenience composition
# ---------------------------------------------------------------------------


def feasible_cells(
    body_set: tuple[str, ...],
    l_max: int,
    k_max: int,
    n_max: int,
    vinf_cap: float,
    ephem: Ephemeris | None = None,
    branch_set: tuple[str, ...] = ("single",),
) -> Iterator[Cell]:
    """Cells from :func:`enumerate_cells` filtered by :func:`tisserand_feasible`.

    The standard M4 entry point for any caller wanting the Tisserand-
    pruned cell stream. The unpruned :func:`enumerate_cells` is exposed
    mostly for tests and for M7's ledger (which counts pruned-vs-searched
    separately).
    """
    return (
        cell
        for cell in enumerate_cells(body_set, l_max, k_max, n_max, branch_set)
        if tisserand_feasible(cell, vinf_cap, ephem)
    )


# ---------------------------------------------------------------------------
# deepening_frontier — the spec §13.2 iterative-deepening generator
# ---------------------------------------------------------------------------


def deepening_frontier(
    body_set: tuple[str, ...],
    ephem: Ephemeris | None = None,
    *,
    vinf_cap: float,
    l_step: int = 1,
    k_step: int = 1,
    n_step: int = 1,
    l_initial: int = 3,
    k_initial: int = 1,
    n_initial: int = 0,
    branch_set: tuple[str, ...] = ("single",),
    max_tiers: int | None = None,
) -> Iterator[Cell]:
    """Iterative-deepening frontier per spec §13.2.

    Yields cells in monotonically increasing complexity by raising the
    caps stepwise. After exhausting all feasible cells at
    ``(l_initial, k_initial, n_initial)``, the generator raises the caps
    by ``(l_step, k_step, n_step)`` and continues, yielding only the
    *newly-added* cells at each tier (those whose id was not yielded in
    a prior tier).

    In-memory deduplication via a ``set`` of cell ids. **Single-process,
    in-memory only** — M7's ledger replaces (not extends) this with a
    persistent, parallel-safe work queue. The contract here is simply:
    eventually visits all cells under monotonic cap growth, no repeats,
    no fancy ordering.

    Parameters
    ----------
    body_set:
        Canonical body codes for :func:`enumerate_cells`.
    ephem:
        Passed through to :func:`tisserand_feasible` (unused in M4).
    vinf_cap:
        Common V∞ ceiling, km/s.
    l_step, k_step, n_step:
        Per-tier cap increments. ``>= 1`` to guarantee progress (an
        all-zero step would loop forever yielding nothing new).
    l_initial, k_initial, n_initial:
        Starting caps for the first tier.
    branch_set:
        Allowed Lambert branches.
    max_tiers:
        Optional hard stop on the number of tiers. ``None`` (default) ⇒
        infinite frontier (caller controls termination with
        :func:`itertools.islice` or similar). M4 tests pass an explicit
        cap to keep iteration finite.

    Yields
    ------
    Cell
        Each cell at most once, ordered by tier and then by
        :func:`enumerate_cells`'s lexicographic walk within a tier.
    """
    if l_step < 1 or k_step < 1 or n_step < 1:
        raise ValueError(
            f"step sizes must be >= 1 to guarantee monotonic frontier growth; "
            f"got l_step={l_step}, k_step={k_step}, n_step={n_step}"
        )
    seen: set[str] = set()
    tier = 0
    l_cap = l_initial
    k_cap = k_initial
    n_cap = n_initial
    while max_tiers is None or tier < max_tiers:
        for cell in feasible_cells(
            body_set,
            l_max=l_cap,
            k_max=k_cap,
            n_max=n_cap,
            vinf_cap=vinf_cap,
            ephem=ephem,
            branch_set=branch_set,
        ):
            cid = cell.id
            if cid in seen:
                continue
            seen.add(cid)
            yield cell
        tier += 1
        l_cap += l_step
        k_cap += k_step
        n_cap += n_step


__all__ = [
    "Cell",
    "deepening_frontier",
    "enumerate_cells",
    "feasible_cells",
    "tisserand_feasible",
]
