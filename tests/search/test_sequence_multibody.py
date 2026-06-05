"""M8 Tisserand feasibility for >=3-body sequences (spec §13.3).

Honest scope: tisserand_feasible is COPLANAR (i=0) only — it inherits
tisserand.linkable's restriction (sequence.py:308). These tests pin the
coplanar baseline so the M-3D inclination lift is a reviewed diff, not a
silent behaviour change. They do NOT assert real-geometry feasibility.
"""

from __future__ import annotations

from cyclerfinder.search.sequence import Cell, tisserand_feasible


def _emeeve_cell() -> Cell:
    seq = ("E", "M", "E", "E", "V", "E")  # vem-emeeve-3syn sequence_canonical
    # NOTE: adjacency E-E at index 3 is the same-body loop leg; enumerate_cells
    # forbids it but a catalogue-derived cell can carry it. tisserand_feasible
    # bypasses same-body pairs (Task 3.2) because linkable(X,X,..) is trivially
    # and meaninglessly True; loop legs are validated by M-L/M-ED, not Tisserand.
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )


def _vem_simple_cell() -> Cell:
    seq = ("E", "V", "M", "E")  # no same-body adjacency
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )


def test_vem_simple_sequence_coplanar_linkable_at_7kms() -> None:
    """E-V-M-E: each consecutive pair (E-V, V-M, M-E) is coplanar-linkable
    somewhere in (0.5, 7.0] km/s (spec §13.3).

    Pinned baseline (observed): True — the production function already
    supports distinct-pair VEM sequences in the coplanar model.
    """
    assert tisserand_feasible(_vem_simple_cell(), vinf_cap=7.0) is True


def test_emeeve_loop_leg_bypass_returns_true() -> None:
    """An EMEEVE-style cell with an E-E loop leg returns True under the
    same-body bypass — the loop leg does NOT falsely trigger rejection.

    Finding (plan §3 Task 3.2 step 1): linkable(E, E, vinf) is trivially True
    (a contour equals itself) and physically meaningless, so we skip same-body
    pairs rather than consult linkable. Every distinct pair (E-M, M-E, E-V,
    V-E) is coplanar-linkable in (0.5, 7.0], so the cell is feasible.
    """
    assert tisserand_feasible(_emeeve_cell(), vinf_cap=7.0) is True


def test_distinct_pair_feasibility_unchanged() -> None:
    """The bypass does not alter distinct-pair behaviour: a plain E-V-M-E cell
    still passes exactly as before."""
    assert tisserand_feasible(_vem_simple_cell(), vinf_cap=7.0) is True
