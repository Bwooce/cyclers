"""Phase 4: CR3BP identity block tests for the 6 non-keplerian catalogue rows.

Golden-discipline: every asserted value traces to a published source cited inline.
Missing dynamical identity (jacobi_constant, period_nd, stability_index, etc.)
must appear as data_gaps[] entries, never as guessed values.

Mass-ratio physical constants (Earth-Moon ~0.01215, Sun-Jupiter ~9.54e-4,
Sun-Saturn ~2.86e-4) are standard physical constants and may be golden-asserted
with a citing comment.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

NON_KEPLERIAN_IDS = [
    "arenstorf-em-figure8-1963",
    "genova-aldrin-2015-em-3petal-cycler",
    "wittal-2022-em-cycler-family",
    "hernandez-2017-jovian-ieg-triple-family",
    "russell-strange-2009-jovian-multimoon-family",
    "russell-strange-2009-saturnian-multimoon-family",
]


@pytest.fixture(scope="module")
def catalogue_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


@pytest.fixture(scope="module")
def non_keplerian_by_id(catalogue_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {row["id"]: row for row in catalogue_rows if row.get("id") in NON_KEPLERIAN_IDS}


# ---------------------------------------------------------------------------
# All 6 rows must have a cr3bp block
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("entry_id", NON_KEPLERIAN_IDS)
def test_non_keplerian_row_has_cr3bp_block(
    non_keplerian_by_id: dict[str, dict[str, Any]],
    entry_id: str,
) -> None:
    """Every non-keplerian row must carry an orbit_elements.cr3bp sub-dict."""
    row = non_keplerian_by_id[entry_id]
    oe = row.get("orbit_elements", {})
    assert isinstance(oe, dict), f"{entry_id}: orbit_elements must be a dict"
    assert "cr3bp" in oe, f"{entry_id}: orbit_elements.cr3bp block is missing"
    assert isinstance(oe["cr3bp"], dict), f"{entry_id}: orbit_elements.cr3bp must be a dict"


# ---------------------------------------------------------------------------
# Earth-Moon rows: mass_ratio (physical constant, not guessed)
# ---------------------------------------------------------------------------
# Source: Genova & Aldrin 2015 (AAS-15) p.1 introduction states
# "periodic orbits in the Earth-Moon system (μ = 0.01215)".
# The Cook 2020 explainer (corroborating source for the Arenstorf entry)
# cites μ = 0.012277471.  Both are published values for the Earth-Moon CR3BP
# mass ratio; 0.01215 is the rounded value used in Genova 2015.
EM_MASS_RATIO = pytest.approx(0.01215, rel=1e-3)


def test_arenstorf_mass_ratio(non_keplerian_by_id: dict[str, dict[str, Any]]) -> None:
    """Arenstorf row carries mass_ratio = 0.01215 (sourced: Genova 2015 p.1)."""
    cr3bp = non_keplerian_by_id["arenstorf-em-figure8-1963"]["orbit_elements"]["cr3bp"]
    assert "mass_ratio" in cr3bp, "arenstorf: mass_ratio missing from cr3bp block"
    assert cr3bp["mass_ratio"] == EM_MASS_RATIO


def test_arenstorf_family_name(non_keplerian_by_id: dict[str, dict[str, Any]]) -> None:
    """Arenstorf row carries family = 'figure-8' (sourced: Genova 2015 p.1)."""
    cr3bp = non_keplerian_by_id["arenstorf-em-figure8-1963"]["orbit_elements"]["cr3bp"]
    assert "family" in cr3bp, "arenstorf: family missing from cr3bp block"
    assert cr3bp["family"] == "figure-8"


def test_genova_mass_ratio(non_keplerian_by_id: dict[str, dict[str, Any]]) -> None:
    """Genova/Aldrin row carries mass_ratio = 0.01215 (sourced: Genova 2015 p.1)."""
    cr3bp = non_keplerian_by_id["genova-aldrin-2015-em-3petal-cycler"]["orbit_elements"]["cr3bp"]
    assert "mass_ratio" in cr3bp, "genova: mass_ratio missing from cr3bp block"
    assert cr3bp["mass_ratio"] == EM_MASS_RATIO


def test_genova_family_name(non_keplerian_by_id: dict[str, dict[str, Any]]) -> None:
    """Genova/Aldrin row carries family = '3:1-lunar-resonance' (sourced: Genova 2015 p.2)."""
    cr3bp = non_keplerian_by_id["genova-aldrin-2015-em-3petal-cycler"]["orbit_elements"]["cr3bp"]
    assert "family" in cr3bp, "genova: family missing from cr3bp block"
    assert cr3bp["family"] == "3:1-lunar-resonance"


def test_wittal_mass_ratio(non_keplerian_by_id: dict[str, dict[str, Any]]) -> None:
    """Wittal row carries mass_ratio = 0.01215 (Earth-Moon standard constant)."""
    cr3bp = non_keplerian_by_id["wittal-2022-em-cycler-family"]["orbit_elements"]["cr3bp"]
    assert "mass_ratio" in cr3bp, "wittal: mass_ratio missing from cr3bp block"
    assert cr3bp["mass_ratio"] == EM_MASS_RATIO


# ---------------------------------------------------------------------------
# Jovian / Saturnian rows: model is circular-coplanar, not CR3BP;
# cr3bp block exists but is empty (all values gapped).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "entry_id",
    [
        "hernandez-2017-jovian-ieg-triple-family",
        "russell-strange-2009-jovian-multimoon-family",
        "russell-strange-2009-saturnian-multimoon-family",
    ],
)
def test_jovian_saturnian_cr3bp_block_exists_but_no_mass_ratio(
    non_keplerian_by_id: dict[str, dict[str, Any]],
    entry_id: str,
) -> None:
    """Jovian/Saturnian rows have a cr3bp block but no mass_ratio
    (these use circular-coplanar patched-conic, not CR3BP; no published
    CR3BP identity parameters were found in accessible sources)."""
    cr3bp = non_keplerian_by_id[entry_id]["orbit_elements"]["cr3bp"]
    # The block must exist (tested above), but mass_ratio should NOT be set
    # (no single Sun-moon mass ratio governs these multi-moon patched-conic cyclers)
    assert cr3bp.get("mass_ratio") is None, (
        f"{entry_id}: mass_ratio should be null/absent for "
        "circular-coplanar Jovian/Saturnian cyclers"
    )


# ---------------------------------------------------------------------------
# Gapped dynamical identity: jacobi_constant, period_nd, stability_index
# must be absent or null AND accompanied by data_gaps[] entries for each.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("entry_id", NON_KEPLERIAN_IDS)
def test_gapped_dynamical_identity_has_data_gaps_entries(
    non_keplerian_by_id: dict[str, dict[str, Any]],
    entry_id: str,
) -> None:
    """If jacobi_constant / period_nd / stability_index are absent or null in cr3bp,
    there must be corresponding data_gaps[] entries documenting the gap."""
    row = non_keplerian_by_id[entry_id]
    cr3bp = row["orbit_elements"]["cr3bp"]
    data_gaps = row.get("data_gaps", [])
    gap_paths = {d.get("path", "") for d in data_gaps}

    for field in ("jacobi_constant", "period_nd", "stability_index"):
        val = cr3bp.get(field)
        if val is None:
            # The gap must be documented
            expected_path = f"orbit_elements.cr3bp.{field}"
            assert expected_path in gap_paths, (
                f"{entry_id}: cr3bp.{field} is absent/null but no data_gaps "
                f"entry for path '{expected_path}' was found. "
                "Add a data_gaps entry or populate the value from a source."
            )
