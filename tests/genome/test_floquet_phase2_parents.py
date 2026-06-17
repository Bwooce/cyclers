"""#347 Phase 2 P2.2 — Sweep parents enumeration tests.

Gates:

  1. The Phase 2 parent inventory has at least 4 entries (the four Braik-Ross
     Table 2 cyclers C11a, C11b, C21, C32 — the family-discovery anchors)
     and is capped at 8 (per the Phase 2 mandate "Cap initial sweep at ~6-8
     (k1, k2) pairs to keep scope tractable").
  2. Every parent traces to a published source (sourced_period_days in the
     Braik-Ross Table 2 range; sourced_sigma_d derivable from the published
     sigma in TU⁻¹).
  3. Every parent's label is recognised by the existing recovery pipeline
     :func:`cyclerfinder.search.reachable_representatives.recover_all_cyclers_braik_ross`
     (modulo the `_down` suffix for the inverse-direction variants).
  4. Every parent's cj_window encloses the jacobi_anchor.
  5. Every parent's (k1, k2) is in the set Braik-Ross Table 2 documents
     ((1, 1), (2, 1), (3, 2)).
  6. The (3, 2) C32 entry's window matches Phase 1's P1.2 walk window
     (3.1294, 3.1544) — the control case.
"""

from __future__ import annotations

from cyclerfinder.genome.floquet_phase2_parents import (
    PHASE2_SWEEP_PARENTS,
    SweepParent,
)
from cyclerfinder.search.reachable_representatives import (
    C_J_BRAIK_ROSS,
    C_J_C21,
    SOURCED_PERIODS_DAYS,
)


def test_phase2_parent_inventory_size_in_4_to_8() -> None:
    """The inventory has 4-8 entries per the Phase 2 mandate."""
    assert 4 <= len(PHASE2_SWEEP_PARENTS) <= 8, (
        f"PHASE2_SWEEP_PARENTS has {len(PHASE2_SWEEP_PARENTS)} entries; "
        f"expected 4-8 per the Phase 2 mandate"
    )


def test_phase2_parents_cover_known_braik_ross_cyclers() -> None:
    """All four Braik-Ross Table 2 cycler labels (C11a, C11b, C21, C32) are present."""
    labels = {p.label.removesuffix("_down") for p in PHASE2_SWEEP_PARENTS}
    expected = {"C11a", "C11b", "C21", "C32"}
    missing = expected - labels
    assert not missing, f"Phase 2 parent inventory missing Braik-Ross cyclers: {missing}"


def test_phase2_parents_sourced_periods_match_table2() -> None:
    """Every parent's sourced_period_days matches the recovery pipeline's Table 2 value."""
    for parent in PHASE2_SWEEP_PARENTS:
        base_label = parent.label.removesuffix("_down")
        sourced = SOURCED_PERIODS_DAYS[base_label]
        assert parent.sourced_period_days == sourced, (
            f"parent {parent.label} sourced_period_days {parent.sourced_period_days} "
            f"!= SOURCED_PERIODS_DAYS[{base_label}] {sourced}"
        )


def test_phase2_parents_jacobi_anchor_inside_cj_window() -> None:
    """Every parent's jacobi_anchor lies within its cj_window."""
    for parent in PHASE2_SWEEP_PARENTS:
        c_min, c_max = parent.cj_window
        assert c_min <= parent.jacobi_anchor <= c_max, (
            f"parent {parent.label} jacobi_anchor {parent.jacobi_anchor} "
            f"outside cj_window ({c_min}, {c_max})"
        )


def test_phase2_parents_topology_subset_of_braik_ross_table2() -> None:
    """Every (k1, k2) is in the Braik-Ross Table 2 cycler set: (1, 1), (2, 1), (3, 2)."""
    allowed = {(1, 1), (2, 1), (3, 2)}
    for parent in PHASE2_SWEEP_PARENTS:
        assert (parent.k1, parent.k2) in allowed, (
            f"parent {parent.label} topology ({parent.k1}, {parent.k2}) "
            f"not in Braik-Ross Table 2 set {allowed}"
        )


def test_phase2_c21_uses_unrounded_jacobi() -> None:
    """The (2, 1) C21 family must sit at the unrounded C_J_C21 (NOT the rounded 3.1294)."""
    c21_parents = [p for p in PHASE2_SWEEP_PARENTS if p.label.removesuffix("_down") == "C21"]
    assert c21_parents, "expected at least one C21 entry"
    for p in c21_parents:
        assert p.jacobi_anchor == C_J_C21, (
            f"C21 parent {p.label} sits at jacobi_anchor={p.jacobi_anchor} != C_J_C21={C_J_C21} — "
            f"per #262 the (2, 1) family only exists at the unrounded value"
        )


def test_phase2_c11a_c11b_c32_use_literal_braik_ross_jacobi() -> None:
    """C11a/C11b/C32 sit at the literal Braik-Ross C_J=3.1294, per #262 reproduction."""
    for p in PHASE2_SWEEP_PARENTS:
        base = p.label.removesuffix("_down")
        if base in {"C11a", "C11b", "C32"}:
            assert p.jacobi_anchor == C_J_BRAIK_ROSS, (
                f"parent {p.label} (base {base}) sits at jacobi_anchor={p.jacobi_anchor} "
                f"!= C_J_BRAIK_ROSS={C_J_BRAIK_ROSS}"
            )


def test_phase2_c32_anchor_window_matches_phase1_p1_2() -> None:
    """The C32 upward-direction entry's cj_window matches Phase 1 P1.2's walk window.

    Phase 1 walked C ∈ [3.1294, 3.1544] in 250 steps at dC=1e-4; the Phase 2
    re-sweep on (3, 2) must use the same parameters to reproduce the
    saddle-center at C* ∈ (3.14170, 3.14180).
    """
    c32_up = next(p for p in PHASE2_SWEEP_PARENTS if p.label == "C32")
    assert c32_up.cj_window == (3.1294, 3.1544), (
        f"C32 cj_window {c32_up.cj_window} != Phase 1 P1.2's (3.1294, 3.1544)"
    )
    assert c32_up.dc == 1e-4
    assert c32_up.n_steps == 250


def test_phase2_parents_are_frozen_dataclass() -> None:
    """SweepParent is frozen — the inventory is immutable across tests."""
    import pytest

    p = PHASE2_SWEEP_PARENTS[0]
    with pytest.raises((AttributeError, Exception)):
        p.label = "modified"  # type: ignore[misc]


def test_phase2_parent_dc_positive() -> None:
    """All dc values are positive (direction is encoded in the label suffix _down).

    The driver picks `direction=-1` for `_down` labels and `direction=+1` otherwise.
    """
    for parent in PHASE2_SWEEP_PARENTS:
        assert parent.dc > 0, f"parent {parent.label} dc={parent.dc} <= 0"


def test_phase2_parent_n_steps_positive() -> None:
    """All n_steps are positive."""
    for parent in PHASE2_SWEEP_PARENTS:
        assert parent.n_steps > 0, f"parent {parent.label} n_steps={parent.n_steps} <= 0"


def test_phase2_sweep_parent_dataclass_contract() -> None:
    """Smoke test: SweepParent can be constructed with the documented fields."""
    p = SweepParent(
        label="test",
        k1=1,
        k2=1,
        jacobi_anchor=3.1294,
        cj_window=(3.0, 3.2),
        dc=1e-4,
        n_steps=10,
        sourced_period_days=42.0,
        sourced_sigma_d=0.5,
        notes="smoke",
    )
    assert p.label == "test"
    assert p.k1 == 1 and p.k2 == 1
