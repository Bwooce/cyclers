"""Tests for the M6b catalogue loader (test infrastructure).

Per plan §4.5: the loader must filter to circular-coplanar, impulsive
(ballistic or powered), Sun-primary entries; reject cr3bp and
analytic-ephemeris entries; and contain every id in
:data:`M6B_REGRESSION_IDS`.
"""

from __future__ import annotations

from tests.data._catalogue_loader_m6b import (
    M6B_REGRESSION_IDS,
    load_m6b_entries,
)


def test_loader_filters_v1_pass_circular_coplanar_ballistic_sun_only() -> None:
    """Plan §4.5: returned set is bounded; every entry's fields obey the filter."""
    entries = load_m6b_entries()
    # Bounded count: catalogue ~ 302 entries; M6b scope drops the cr3bp,
    # analytic-ephemeris, and non-Sun rows. Plan §3.2 expects ~180-215 at
    # original publish; ceiling raised 250->270 after #367 wave 2 admitted
    # 11 Rogers 2015 Table 3 circular-coplanar precursor_mga rows (all
    # M6b-scope eligible).
    assert 100 <= len(entries) <= 270, f"unexpected M6b-scope entry count: {len(entries)}"
    for entry in entries:
        assert entry.get("model_assumption") in (None, "circular-coplanar")
        assert entry.get("trajectory_regime") in (None, "ballistic", "powered")
        assert entry.get("primary") in (None, "Sun")


def test_loader_excludes_cr3bp_entries() -> None:
    entries = load_m6b_entries()
    for entry in entries:
        assert entry.get("model_assumption") != "cr3bp", entry.get("id")


def test_loader_excludes_analytic_ephemeris_entries() -> None:
    entries = load_m6b_entries()
    for entry in entries:
        assert entry.get("model_assumption") != "analytic-ephemeris", entry.get("id")


def test_loader_excludes_non_sun_primaries() -> None:
    entries = load_m6b_entries()
    for entry in entries:
        primary = entry.get("primary")
        assert primary in (None, "Sun"), f"{entry.get('id')}: primary={primary!r}"


def test_m6b_regression_ids_all_in_loader() -> None:
    entries = load_m6b_entries()
    ids = {entry["id"] for entry in entries}
    for rid in M6B_REGRESSION_IDS:
        assert rid in ids, f"regression id {rid!r} missing from loader output"


def test_loader_returns_aldrin_outbound() -> None:
    """Aldrin outbound is the M6b binding gate's fixture; must be loadable."""
    entries = load_m6b_entries()
    aldrin = next((e for e in entries if e["id"] == "aldrin-classic-em-k1-outbound"), None)
    assert aldrin is not None
    assert aldrin["bodies"] == ["E", "M"]
    assert aldrin["model_assumption"] == "circular-coplanar"
    # Powered (reclassified 2026-06-01, task #70): the 1L1 Mars flyby cannot
    # supply the required 84 deg turn (max 72 deg) even in the idealized model.
    assert aldrin["trajectory_regime"] == "powered"
