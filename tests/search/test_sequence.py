"""Tests for :mod:`cyclerfinder.search.sequence` — the M4 enumeration layer.

Spec gate (plan §4 / spec §8): two Tisserand-pruning tests at the Aldrin
neighbourhood, the enumeration-count derivation at ``(E,M), L=4, k=2``,
and the spec §13.8 ``Cell.id`` worked example.

Plan: ``docs/phases/m4-enumeration-scoring/plan.md`` §4.1.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterator
from dataclasses import FrozenInstanceError

import pytest

from cyclerfinder.search.sequence import (
    Cell,
    deepening_frontier,
    enumerate_cells,
    feasible_cells,
    tisserand_feasible,
)

# ---------------------------------------------------------------------------
# Tolerances and shared anchors (module-level per house style)
# ---------------------------------------------------------------------------

# Aldrin neighbourhood V∞ — M2 confirmed linkable("E","M", 5.5) is True.
ALDRIN_VINF_KMS: float = 5.5
LOW_VINF_CAP_KMS: float = 2.0  # below Aldrin family → E-M not linkable
HIGH_VINF_CAP_KMS: float = 8.0  # above Aldrin family → E-M linkable


# ---------------------------------------------------------------------------
# Cell dataclass
# ---------------------------------------------------------------------------


def test_cell_is_frozen() -> None:
    """Cell is a frozen dataclass — assignment raises ``FrozenInstanceError``."""
    c = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    with pytest.raises(FrozenInstanceError):
        c.period_k = 2  # type: ignore[misc]


def test_cell_id_format_spec_worked_example() -> None:
    """Spec §13.8 worked example: ``VEM|E-V-M-E-M-E|k3|r00101|blllll``.

    Five legs all multi-rev on the ``low`` branch (revs 0 use ``"single"``
    which would map to ``s``; revs 1 use ``"low"`` → ``l``). The worked
    example uses ``low`` everywhere except the 0-rev legs which take
    ``single``. Per the plan §3.1.1 branch-alphabet mapping
    ``{"single":"s","low":"l","high":"h"}``, the cell with the rev
    tuple ``(0,0,1,0,1)`` and the matching branch tuple
    ``("single","single","low","single","low")`` resolves to
    ``VEM|E-V-M-E-M-E|k3|r00101|bsslsl``. (The spec's ``blllll`` example
    corresponds to a cell with all five legs multi-rev on the low branch,
    revs e.g. ``(1,1,1,1,1)`` and branches ``("low",)*5``.)
    """
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=("E", "V", "M", "E", "M", "E"),
        period_k=3,
        per_leg_revs=(0, 0, 1, 0, 1),
        per_leg_branch=("single", "single", "low", "single", "low"),
    )
    assert cell.id == "VEM|E-V-M-E-M-E|k3|r00101|bsslsl"


def test_cell_id_format_all_low_branches() -> None:
    """The spec's literal ``blllll`` example: all five legs on the low branch."""
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=("E", "V", "M", "E", "M", "E"),
        period_k=3,
        per_leg_revs=(1, 1, 1, 1, 1),
        per_leg_branch=("low", "low", "low", "low", "low"),
    )
    assert cell.id == "VEM|E-V-M-E-M-E|k3|r11111|blllll"


def test_cell_id_format_em_2syn_direct() -> None:
    """The 2-synodic E-M-E direct case — the M4 native cell."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert cell.id == "EM|E-M-E|k2|r00|bss"


def test_cell_hashable_in_dict() -> None:
    """Cells are usable as dict keys (M7 ledger dependency)."""
    c1 = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    c2 = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    d: dict[Cell, str] = {c1: "x"}
    assert d[c2] == "x"  # equal cells map to the same key
    assert hash(c1) == hash(c2)


def test_cell_id_uniqueness_in_em_l4_k2() -> None:
    """All 12 cells in the gate enumeration have distinct ids."""
    cells = list(enumerate_cells(("E", "M"), l_max=4, k_max=2, n_max=0))
    ids = {c.id for c in cells}
    assert len(ids) == len(cells), f"id collisions: {len(cells) - len(ids)}"


# ---------------------------------------------------------------------------
# enumerate_cells — counts and invariants (plan §4.4 gate)
# ---------------------------------------------------------------------------


def test_enumeration_count_em_l2_k1() -> None:
    """Plan §4.4: ``(E,M), L=2, k=1, N=0`` → 2 cells (``E-M`` and ``M-E``)."""
    cells = list(enumerate_cells(("E", "M"), l_max=2, k_max=1, n_max=0, branch_set=("single",)))
    assert len(cells) == 2


def test_enumeration_count_em_l4_k2_gate() -> None:
    """**Gate (plan §4.4):** ``(E,M), L=4, k=2, N=0, single`` → 12 cells.

    Per the §4.4 derivation:
      L=2: 2 sequences * 2 k * 1 rev * 1 branch = 4
      L=3: 2 sequences * 2 k * 1 * 1            = 4
      L=4: 2 sequences * 2 k * 1 * 1            = 4
      total                                     = 12
    """
    cells = list(enumerate_cells(("E", "M"), l_max=4, k_max=2, n_max=0, branch_set=("single",)))
    assert len(cells) == 12


def test_enumeration_excludes_consecutive_same_body() -> None:
    """No yielded cell has ``sequence[i] == sequence[i+1]``."""
    for cell in enumerate_cells(("V", "E", "M"), l_max=5, k_max=1, n_max=0):
        for i in range(len(cell.sequence) - 1):
            assert cell.sequence[i] != cell.sequence[i + 1], cell.id


def test_enumeration_iterator_not_list() -> None:
    """``enumerate_cells`` returns a lazy iterator, not a materialised list."""
    it = enumerate_cells(("E", "M"), l_max=3, k_max=1, n_max=0)
    assert isinstance(it, Iterator)
    assert not isinstance(it, list)


def test_enumeration_multirev_requires_branch() -> None:
    """At ``n_max >= 1`` with default ``branch_set=("single",)`` only direct
    legs (rev=0) yield cells — multi-rev rev tuples have no valid branch
    and are silently skipped (documented behaviour).
    """
    cells_single_only = list(
        enumerate_cells(("E", "M"), l_max=3, k_max=1, n_max=1, branch_set=("single",))
    )
    cells_with_low = list(
        enumerate_cells(("E", "M"), l_max=3, k_max=1, n_max=1, branch_set=("single", "low"))
    )
    assert len(cells_with_low) > len(cells_single_only), (
        f"{len(cells_with_low)=} should exceed {len(cells_single_only)=}"
    )


# ---------------------------------------------------------------------------
# Tisserand pruning gate (plan §4.1 binding gate)
# ---------------------------------------------------------------------------


def test_tisserand_pruning_rejects_low_vinf_em_gate() -> None:
    """**Gate:** E-M at ``vinf_cap=2.0`` km/s → False (well below Aldrin).

    The Aldrin family lives at V∞ ≈ 5.5 km/s; at vinf_cap=2.0 the
    contours of Earth and Mars do not intersect.
    """
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    assert tisserand_feasible(cell, vinf_cap=LOW_VINF_CAP_KMS) is False


def test_tisserand_pruning_accepts_aldrin_vinf_em_gate() -> None:
    """**Gate:** E-M at ``vinf_cap=8.0`` km/s → True (covers Aldrin band)."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    assert tisserand_feasible(cell, vinf_cap=HIGH_VINF_CAP_KMS) is True


def test_tisserand_pruning_propagates_through_sequence() -> None:
    """A 3-encounter cell at low vinf_cap fails because some adjacent pair fails.

    Use ``vinf_cap=1.0`` — far below any planetary pair's linkable band —
    so an E-M-V cell must be False because both (E,M) and (M,V) are
    infeasible at this energy.
    """
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=("E", "M", "V"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert tisserand_feasible(cell, vinf_cap=1.0) is False


@pytest.mark.parametrize("vinf_cap", [0.0, -1.0, 1.0e6])
def test_tisserand_feasible_never_raises(vinf_cap: float) -> None:
    """Per :func:`tisserand_feasible`'s contract: never raises."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M"),
        period_k=1,
        per_leg_revs=(0,),
        per_leg_branch=("single",),
    )
    result = tisserand_feasible(cell, vinf_cap=vinf_cap)
    assert isinstance(result, bool)


# ---------------------------------------------------------------------------
# feasible_cells composition
# ---------------------------------------------------------------------------


def test_feasible_cells_subset_of_enumerate() -> None:
    """``feasible_cells(..., vinf_cap=8.0)`` ids ⊆ ``enumerate_cells(...)`` ids."""
    all_ids = {c.id for c in enumerate_cells(("E", "M"), l_max=3, k_max=1, n_max=0)}
    feasible_ids = {
        c.id
        for c in feasible_cells(("E", "M"), l_max=3, k_max=1, n_max=0, vinf_cap=HIGH_VINF_CAP_KMS)
    }
    assert feasible_ids <= all_ids
    assert len(feasible_ids) > 0  # some cells must survive at 8 km/s


def test_feasible_cells_strict_subset_at_low_cap() -> None:
    """At ``vinf_cap=1.0`` the feasible set is strictly smaller than the
    unpruned set — the gate is doing work.
    """
    all_ids = {c.id for c in enumerate_cells(("E", "M"), l_max=3, k_max=1, n_max=0)}
    feasible_ids = {
        c.id for c in feasible_cells(("E", "M"), l_max=3, k_max=1, n_max=0, vinf_cap=1.0)
    }
    assert feasible_ids < all_ids


# ---------------------------------------------------------------------------
# deepening_frontier
# ---------------------------------------------------------------------------


def test_deepening_frontier_yields_in_increasing_complexity() -> None:
    """Sequence length is monotonically non-decreasing when only ``L`` is
    raised across tiers (with ``k`` and ``N`` held large enough to be
    fully exhausted in tier 0).

    Plan §3.1.5: "yields cells in monotonically increasing complexity by
    raising the caps stepwise". The strict-length monotonicity only holds
    when the cap being raised tier-to-tier is ``L``; if ``k`` or ``N``
    also bumps, length-``L`` cells with the new ``k``/``N`` are
    legitimately new at the later tier. We test the strict case here by
    pinning ``k`` and ``N`` to their initial values via huge steps that
    we never actually reach (max_tiers caps the loop first).
    """
    # Initial: L=2, k=2, n=0 — exhausts all small-L cells in tier 0.
    # Steps: l=1 raises L tier-by-tier; k=999, n=999 effectively never
    # bump within max_tiers.
    cells = list(
        itertools.islice(
            deepening_frontier(
                ("E", "M"),
                vinf_cap=HIGH_VINF_CAP_KMS,
                l_initial=2,
                k_initial=2,
                n_initial=0,
                l_step=1,
                k_step=999,
                n_step=999,
                max_tiers=5,
            ),
            20,
        )
    )
    assert len(cells) > 0
    lengths = [len(c.sequence) for c in cells]
    for i in range(len(lengths) - 1):
        assert lengths[i] <= lengths[i + 1], f"non-monotone at {i}: {lengths}"


def test_deepening_frontier_no_repeats() -> None:
    """First 50 cells have distinct ids — the in-memory dedup works."""
    cells = list(
        itertools.islice(
            deepening_frontier(
                ("E", "M"),
                vinf_cap=HIGH_VINF_CAP_KMS,
                l_initial=2,
                k_initial=1,
                n_initial=0,
                max_tiers=10,
            ),
            50,
        )
    )
    ids = [c.id for c in cells]
    assert len(set(ids)) == len(ids), f"duplicate ids in frontier: {len(ids) - len(set(ids))}"


def test_deepening_frontier_step_validation() -> None:
    """All-zero step sizes are rejected — they would loop forever yielding nothing new."""
    with pytest.raises(ValueError, match="step sizes must be >= 1"):
        list(
            itertools.islice(
                deepening_frontier(
                    ("E", "M"),
                    vinf_cap=HIGH_VINF_CAP_KMS,
                    l_step=0,
                    max_tiers=1,
                ),
                1,
            )
        )
