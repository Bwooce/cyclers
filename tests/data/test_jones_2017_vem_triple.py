"""Data-integrity (V0) gate for the Jones 2017 VEM triple-cycler catalogue rows.

Asserts that the sourced anchors (per-encounter V∞ multiset, 12.8-yr repeat
period, transit-leg ToFs, sequence/class) in the two newly-ingested member
rows match the values read directly from:

  Jones, Hernandez & Jesick (2017). "Low Excess Speed Triple Cyclers of Venus,
  Earth, and Mars." AAS/AIAA Astrodynamics Specialist Conference, AAS Paper
  17-577 (NTRS 20190028464). No DOI (AAS conference paper).

  - jones-2017-vem-emevve-outbound  ← Table 2 (p.10) + p.9 transit-leg text
  - jones-2017-vem-meevem-inbound   ← Table 3 (p.10) + p.9 transit-leg text

This is a traceability test: it guards against the catalogue rows drifting
from their sourced numbers. The expected values below are the ONLY legitimate
source for this test — they were read directly from the paper's tables and
prose and must NOT be replaced with values computed by this codebase
(golden-test discipline).

It also pins the correction made 2026-06-05 to the repeat-period
interpretation: "two synodic period" in this paper means 2x the 6.4-yr VEM
synodic = 12.8 yr (p.9: "the repeat period T is 12.8 years"), NOT 2x the
2.135-yr E-M synodic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

_REPO = Path(__file__).resolve().parent.parent.parent
_CATALOGUE = _REPO / "data" / "catalogue.yaml"

_OUTBOUND_ID = "jones-2017-vem-emevve-outbound"
_INBOUND_ID = "jones-2017-vem-meevem-inbound"
_FAMILY_ID = "jones-2017-vem-triple-family"
_EMEEVE_ID = "vem-emeeve-3syn"

# ---------------------------------------------------------------------------
# SOURCED EXPECTED values — AAS 17-577 Tables 2 & 3 (p.10) + p.9 text.
#
# Table 2 (EMEVVE outbound) "Excess speed, km/sec" column, encounter order:
_EMEVVE_VINF_KMS: list[float] = [4.72, 2.50, 5.81, 7.00, 7.00, 4.21, 2.79, 5.04, 4.49, 4.49, 5.67]
# Table 2 "Periapsis altitude, km" column, encounter order:
_EMEVVE_ALT_KM: list[float] = [100, 4164, 3814, 684, 1985, 1998, 1754, 3213, 3319, 836, 100]
# p.9: outbound transit legs are "309 and 259 days" (E→M 309, M→E 259):
_EMEVVE_TRANSIT_DAYS: list[int] = [309, 259]

# Table 3 (MEEVEM inbound) "Excess speed, km/sec" column, encounter order:
_MEEVEM_VINF_KMS: list[float] = [3.85, 3.48, 3.42, 5.16, 3.84, 3.12, 2.98, 2.92, 4.31, 4.97, 2.42]
# Table 3 "Periapsis altitude, km" column, encounter order:
_MEEVEM_ALT_KM: list[float] = [100, 7831, 967, 29777, 2545, 249, 3719, 2224, 12349, 7484, 100]
# p.9: inbound transit legs are "268 and 223 days" (M→E 268, E→M 223):
_MEEVEM_TRANSIT_DAYS: list[int] = [268, 223]

# p.9: "Recall that the repeat period T is 12.8 years."
_REPEAT_PERIOD_YEARS: float = 12.8
# ---------------------------------------------------------------------------


def _load_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(_CATALOGUE.read_text())  # type: ignore[no-any-return]


def _row(row_id: str) -> dict[str, Any]:
    matches = [r for r in _load_rows() if r["id"] == row_id]
    assert matches, f"catalogue row {row_id!r} not found"
    return matches[0]


# --- structural / class anchors --------------------------------------------


@pytest.mark.parametrize("row_id", [_OUTBOUND_ID, _INBOUND_ID])
def test_member_rows_are_multi_arc(row_id: str) -> None:
    """Both new member rows are tagged multi-arc (6-flyby VEM triple cyclers)."""
    row = _row(row_id)
    assert row["cycler_class"] == "multi-arc", (
        f"{row_id}: expected cycler_class=multi-arc, got {row['cycler_class']!r}"
    )


@pytest.mark.parametrize("row_id", [_OUTBOUND_ID, _INBOUND_ID])
def test_member_rows_are_analytic_ephemeris(row_id: str) -> None:
    """Model assumption is analytic-ephemeris (real-ephemeris patched conic, p.5)."""
    row = _row(row_id)
    assert row["model_assumption"] == "analytic-ephemeris", (
        f"{row_id}: expected analytic-ephemeris, got {row['model_assumption']!r}"
    )


@pytest.mark.parametrize("row_id", [_OUTBOUND_ID, _INBOUND_ID])
def test_member_rows_orbit_elements_null(row_id: str) -> None:
    """Multi-arc invariant: a_au/e null (AAS 17-577 publishes no orbital elements)."""
    row = _row(row_id)
    oe = row.get("orbit_elements") or {}
    assert oe.get("a_au") is None, f"{row_id}: a_au must be null for multi-arc entry"
    assert oe.get("e") is None, f"{row_id}: e must be null for multi-arc entry"


def test_outbound_sequence_and_sense() -> None:
    """EMEVVE outbound: sequence E-M-E-V-V-E, sense outbound (Table 1/2)."""
    row = _row(_OUTBOUND_ID)
    assert row["sequence_canonical"] == "E-M-E-V-V-E", (
        f"expected E-M-E-V-V-E, got {row['sequence_canonical']!r}"
    )
    assert row["sense"] == "outbound"


def test_inbound_sequence_and_sense() -> None:
    """MEEVEM inbound: sequence M-E-E-V-E-M, sense inbound (Table 1/3)."""
    row = _row(_INBOUND_ID)
    assert row["sequence_canonical"] == "M-E-E-V-E-M", (
        f"expected M-E-E-V-E-M, got {row['sequence_canonical']!r}"
    )
    assert row["sense"] == "inbound"


# --- period anchor (the 12.8-yr correction) --------------------------------


@pytest.mark.parametrize("row_id", [_OUTBOUND_ID, _INBOUND_ID])
def test_repeat_period_12_8_years(row_id: str) -> None:
    """Repeat period T = 12.8 yr — SOURCED: AAS 17-577 p.9."""
    period = _row(row_id).get("period") or {}
    assert period.get("years") == pytest.approx(_REPEAT_PERIOD_YEARS, abs=0.01), (
        f"{row_id}: expected period.years = 12.8, got {period.get('years')!r}"
    )


# --- V∞ multiset anchors ----------------------------------------------------


def test_outbound_vinf_multiset_matches_table2() -> None:
    """EMEVVE row V∞ multiset matches AAS 17-577 Table 2 (p.10)."""
    row = _row(_OUTBOUND_ID)
    got = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"])
    exp = sorted(_EMEVVE_VINF_KMS)
    assert len(got) == len(exp), f"expected {len(exp)} V∞ entries, got {len(got)}"
    for g, e in zip(got, exp, strict=True):
        assert g == pytest.approx(e, abs=0.001), f"Table 2 V∞ mismatch: {g} vs {e}"


def test_inbound_vinf_multiset_matches_table3() -> None:
    """MEEVEM row V∞ multiset matches AAS 17-577 Table 3 (p.10)."""
    row = _row(_INBOUND_ID)
    got = sorted(v["vinf_kms"] for v in row["vinf_kms_at_encounters"])
    exp = sorted(_MEEVEM_VINF_KMS)
    assert len(got) == len(exp), f"expected {len(exp)} V∞ entries, got {len(got)}"
    for g, e in zip(got, exp, strict=True):
        assert g == pytest.approx(e, abs=0.001), f"Table 3 V∞ mismatch: {g} vs {e}"


# --- transit-leg ToF anchors (p.9 prose) -----------------------------------


def test_outbound_transit_times_309_259() -> None:
    """EMEVVE transit-leg ToFs = [309, 259] days — SOURCED: AAS 17-577 p.9."""
    inv = _row(_OUTBOUND_ID).get("invariants") or {}
    assert inv.get("transit_times_days") == _EMEVVE_TRANSIT_DAYS, (
        f"expected {_EMEVVE_TRANSIT_DAYS}, got {inv.get('transit_times_days')!r}"
    )


def test_inbound_transit_times_268_223() -> None:
    """MEEVEM transit-leg ToFs = [268, 223] days — SOURCED: AAS 17-577 p.9."""
    inv = _row(_INBOUND_ID).get("invariants") or {}
    assert inv.get("transit_times_days") == _MEEVEM_TRANSIT_DAYS, (
        f"expected {_MEEVEM_TRANSIT_DAYS}, got {inv.get('transit_times_days')!r}"
    )


# --- periapsis-altitude anchors --------------------------------------------


def test_outbound_altitudes_match_table2() -> None:
    """EMEVVE flyby periapsis altitudes match AAS 17-577 Table 2 (p.10)."""
    fm = _row(_OUTBOUND_ID)["flyby_mechanics"]
    got = [m["min_altitude_km"] for m in fm]
    assert got == _EMEVVE_ALT_KM, f"Table 2 altitude mismatch: {got} vs {_EMEVVE_ALT_KM}"


def test_inbound_altitudes_match_table3() -> None:
    """MEEVEM flyby periapsis altitudes match AAS 17-577 Table 3 (p.10)."""
    fm = _row(_INBOUND_ID)["flyby_mechanics"]
    got = [m["min_altitude_km"] for m in fm]
    assert got == _MEEVEM_ALT_KM, f"Table 3 altitude mismatch: {got} vs {_MEEVEM_ALT_KM}"


# --- segment-ToF sum cross-check (internal consistency, not a golden anchor)-


@pytest.mark.parametrize(
    ("row_id", "n_segments"),
    [(_OUTBOUND_ID, 10), (_INBOUND_ID, 10)],
)
def test_member_rows_have_ten_segments(row_id: str, n_segments: int) -> None:
    """Two repeat periods over 11 encounters ⇒ 10 inter-encounter segments."""
    segs = (_row(row_id).get("trajectory") or {}).get("segments") or []
    assert len(segs) == n_segments, f"{row_id}: expected {n_segments} segments, got {len(segs)}"


# --- placeholder-row corrections (Task 2) ----------------------------------


def test_family_seed_period_corrected_to_12_8() -> None:
    """jones-2017-vem-triple-family period corrected 4.27 → 12.8 yr (p.9)."""
    period = _row(_FAMILY_ID).get("period") or {}
    assert period.get("years") == pytest.approx(12.8, abs=0.01), (
        f"family seed period.years should be 12.8, got {period.get('years')!r}"
    )


def test_emeeve_3syn_records_no_feasible_conflict() -> None:
    """vem-emeeve-3syn keeps a data_gap recording the empty 6.4-yr family (p.8)."""
    row = _row(_EMEEVE_ID)
    gaps = row.get("data_gaps") or []
    paths = {g.get("path") for g in gaps}
    assert "period.feasibility" in paths, (
        "vem-emeeve-3syn must carry a data_gap at period.feasibility recording "
        "that Jones 2017 found NO feasible 6.4-yr (1-synodic) cyclers (p.8)"
    )
    # The 6.4-yr beat value is retained (archetype), not deleted.
    assert (row.get("period") or {}).get("years") == pytest.approx(6.41, abs=0.01)
