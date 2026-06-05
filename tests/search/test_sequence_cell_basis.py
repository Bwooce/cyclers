"""M8: Cell.period_basis carries the catalogue's anchor pair without
mutating period_k (plan §2 design decision)."""

from __future__ import annotations

from cyclerfinder.search.sequence import Cell


def _emeeve_basis_cell() -> Cell:
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=3,  # sourced anchor-pair k, NOT rewritten
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=("E", "M"),
    )


def test_period_basis_defaults_to_none() -> None:
    """A 2-body cell built the old way has period_basis is None."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert cell.period_basis is None


def test_period_basis_preserved_and_period_k_untouched() -> None:
    cell = _emeeve_basis_cell()
    assert cell.period_basis == ("E", "M")
    assert cell.period_k == 3  # traceable to catalogue.yaml:1782


def test_id_unchanged_when_basis_none() -> None:
    """Existing 2-body cell ids stay byte-identical (no ledger churn)."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert cell.id == "EM|E-M-E|k2|r00|bss"


def test_id_appends_basis_token_when_set() -> None:
    """When period_basis is set, the id gains a |p<AB> token so a basis-bearing
    cell is distinguishable (and the YAML-traceable k3 stays visible)."""
    cell = _emeeve_basis_cell()
    assert cell.id == "VEM|E-M-E-E-V-E|k3|r00000|bsssss|pEM"
