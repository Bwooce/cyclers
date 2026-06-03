"""Phase 5: period.basis tests for n-body (VEM) cycler rows.

Golden-discipline: period.basis is populated ONLY where the source explicitly
states the synodic commensurabilities. For rows where the source does not give
per-pair resonances, a data_gaps[] entry documents the gap; basis is absent.

2-body rows must be unchanged (no period.basis added).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_schema_invariants

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"


@pytest.fixture(scope="module")
def catalogue_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


@pytest.fixture(scope="module")
def by_id(catalogue_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {row["id"]: row for row in catalogue_rows}


# ---------------------------------------------------------------------------
# vem-emeeve-3syn — 1 T_syn cycle (3 E-M ≈ 4 E-V).
# Source: Jones et al. 2017 (AAS 17-577) p.3:
#   "T_syn is approximately 6.4 years. This is about three Earth-Mars synodic
#    periods."  The same page states the three planets align every ~32 yr = 5 T_syn,
#    confirming T_syn = ~6.4 yr.  The 4xE-V commensurability follows from
#    4 x 1.599 yr ~= 6.396 yr ~= T_syn (stated in the existing period.note).
# ---------------------------------------------------------------------------


def test_vem_emeeve_3syn_has_period_basis(by_id: dict[str, dict[str, Any]]) -> None:
    """vem-emeeve-3syn must carry period.basis documenting the E-M and E-V resonances."""
    row = by_id["vem-emeeve-3syn"]
    period = row.get("period", {})
    assert "basis" in period, "vem-emeeve-3syn: period.basis is missing"
    assert isinstance(period["basis"], list), "vem-emeeve-3syn: period.basis must be a list"
    assert len(period["basis"]) >= 2, "vem-emeeve-3syn: period.basis must have >=2 entries"


def test_vem_emeeve_3syn_basis_em(by_id: dict[str, dict[str, Any]]) -> None:
    """vem-emeeve-3syn period.basis includes {pair: E-M, k: 3}.
    Source: Jones 2017 p.3 'approximately 6.4 years ... about three Earth-Mars
    synodic periods'."""
    period = by_id["vem-emeeve-3syn"]["period"]
    em_entries = [b for b in period["basis"] if b.get("pair") == "E-M"]
    assert em_entries, "vem-emeeve-3syn: no E-M entry in period.basis"
    assert em_entries[0]["k"] == 3, (
        f"vem-emeeve-3syn: E-M basis k should be 3, got {em_entries[0]['k']}"
    )


def test_vem_emeeve_3syn_basis_ev(by_id: dict[str, dict[str, Any]]) -> None:
    """vem-emeeve-3syn period.basis includes {pair: E-V, k: 4}.
    Source: existing period.note and Jones 2017 (4 x 1.599 yr ~= T_syn ~= 6.4 yr)."""
    period = by_id["vem-emeeve-3syn"]["period"]
    ev_entries = [b for b in period["basis"] if b.get("pair") == "E-V"]
    assert ev_entries, "vem-emeeve-3syn: no E-V entry in period.basis"
    assert ev_entries[0]["k"] == 4, (
        f"vem-emeeve-3syn: E-V basis k should be 4, got {ev_entries[0]['k']}"
    )


# ---------------------------------------------------------------------------
# jones-2017-vem-triple-family — 2 T_syn (2xE-M-synodic-based) family.
# Jones 2017 p.3 states the paper "is restricted to families with 1 or 2
# synodic periods in a cycle".  The paper's own abstract says "two synodic
# period Earth-Mars-Venus triple cyclers" but defines T_syn = 6.4 yr = 3 E-M.
# However the existing catalogue entry interprets "two synodic" as k=2 E-M
# (2 x 2.135 = 4.27 yr), and the paper does NOT explicitly tabulate the
# per-pair basis decomposition for this entry.  Per golden-discipline, if the
# source does not directly state the basis decomposition, we gap it.
# ---------------------------------------------------------------------------


def test_jones_vem_triple_has_data_gap_for_period_basis(
    by_id: dict[str, dict[str, Any]],
) -> None:
    """jones-2017-vem-triple-family: period.basis is absent and a data_gaps
    entry documents the gap (the paper does not give an explicit resonance basis
    for this 2-synodic family entry as modelled in this catalogue row)."""
    row = by_id["jones-2017-vem-triple-family"]
    period = row.get("period", {})
    # Must NOT have a basis (because we cannot source it for this entry's interpretation)
    assert "basis" not in period or period.get("basis") is None, (
        "jones-2017-vem-triple-family: period.basis should be absent for this entry; "
        "the source does not give an explicit per-pair resonance basis for the "
        "2-synodic (4.27 yr) family as interpreted in this row."
    )
    # Must have a data_gaps entry
    data_gaps = row.get("data_gaps", [])
    gap_paths = {d.get("path", "") for d in data_gaps}
    assert "period.basis" in gap_paths, (
        "jones-2017-vem-triple-family: no data_gaps entry for 'period.basis'. "
        "Add one explaining why the basis is not set for this row."
    )


# ---------------------------------------------------------------------------
# 2-body rows must not have period.basis added
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "entry_id",
    [
        "s1l1-2syn-em-cpom",
        "aldrin-classic-em-k1-outbound",
        "aldrin-classic-em-k1-inbound",
    ],
)
def test_two_body_rows_unchanged_no_basis(
    by_id: dict[str, dict[str, Any]],
    entry_id: str,
) -> None:
    """2-body E-M cycler rows must not carry period.basis (unchanged by Phase 5)."""
    row = by_id[entry_id]
    period = row.get("period", {})
    assert "basis" not in period or period.get("basis") is None, (
        f"{entry_id}: 2-body row should not have period.basis added"
    )


# ---------------------------------------------------------------------------
# validate_schema_invariants accepts both forms (with and without basis)
# ---------------------------------------------------------------------------


def test_schema_invariants_accept_basis_form() -> None:
    """validate_schema_invariants returns [] for a row with a valid period.basis."""
    row = {
        "id": "test-vem-basis",
        "period": {
            "years": 6.41,
            "basis": [
                {"pair": "E-M", "k": 3},
                {"pair": "E-V", "k": 4},
            ],
        },
    }
    errs = validate_schema_invariants([row])
    assert errs == [], f"Unexpected invariant errors: {errs}"


def test_schema_invariants_accept_no_basis_form() -> None:
    """validate_schema_invariants returns [] for a row without period.basis."""
    row = {
        "id": "test-2body",
        "period": {"pair": "E-M", "k": 2, "years": 4.27},
    }
    errs = validate_schema_invariants([row])
    assert errs == [], f"Unexpected invariant errors: {errs}"
