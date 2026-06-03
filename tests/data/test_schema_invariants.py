"""Task 1.3: Semantic/cross-field schema invariant validation tests.

Tests for validate_schema_invariants() which enforces rules JSON Schema
cannot express: cross-field semantics (multi-arc/non-keplerian must not
have top-level a/e, non-keplerian implies non-Sun primary, period.basis
items have pair+k).

The live catalogue ratchet asserts the current data/catalogue.yaml
passes — this must remain green through Phase 2.
"""

from __future__ import annotations

from pathlib import Path

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_schema_invariants

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"


# ---------------------------------------------------------------------------
# Multi-arc invariant: no top-level a_au / e
# ---------------------------------------------------------------------------


def test_multi_arc_must_not_have_top_level_a() -> None:
    """multi-arc row with top-level orbit_elements.a_au triggers an error."""
    bad = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"a_au": 1.6}}
    errs = validate_schema_invariants([bad])
    assert any("a_au" in m for m in errs), f"Expected a_au error, got: {errs}"


def test_multi_arc_must_not_have_top_level_e() -> None:
    """multi-arc row with top-level orbit_elements.e triggers an error."""
    bad = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"e": 0.3}}
    errs = validate_schema_invariants([bad])
    assert any("e" in m for m in errs), f"Expected e error, got: {errs}"


def test_multi_arc_without_ae_is_clean() -> None:
    """multi-arc row without a_au/e at orbit_elements level is clean."""
    ok = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"center": "Sun"}}
    assert validate_schema_invariants([ok]) == []


def test_non_keplerian_must_not_have_top_level_a() -> None:
    """non-keplerian row with top-level orbit_elements.a_au triggers an error."""
    bad = {
        "id": "y",
        "cycler_class": "non-keplerian",
        "orbit_elements": {"a_au": 1.0},
        "primary": "Earth",
    }
    errs = validate_schema_invariants([bad])
    assert any("a_au" in m for m in errs), f"Expected a_au error, got: {errs}"


# ---------------------------------------------------------------------------
# Non-keplerian invariant: primary must not be "Sun"
# ---------------------------------------------------------------------------


def test_non_keplerian_implies_non_sun_primary() -> None:
    """non-keplerian row with primary=Sun triggers an error."""
    bad = {"id": "y", "cycler_class": "non-keplerian", "primary": "Sun"}
    errs = validate_schema_invariants([bad])
    assert any("primary" in m for m in errs), f"Expected primary error, got: {errs}"


def test_non_keplerian_non_sun_primary_is_clean() -> None:
    """non-keplerian row with primary=Earth (non-Sun) is clean (no primary error)."""
    ok = {"id": "y", "cycler_class": "non-keplerian", "primary": "Earth"}
    errs = validate_schema_invariants([ok])
    # Should have no primary-related error
    assert not any("primary" in m for m in errs), f"Unexpected primary error: {errs}"


# ---------------------------------------------------------------------------
# period.basis items must have pair + k
# ---------------------------------------------------------------------------


def test_period_basis_items_must_have_pair_and_k() -> None:
    """period.basis item missing 'pair' triggers an error."""
    bad = {
        "id": "z",
        "period": {"basis": [{"k": 3}, {"pair": "E-M", "k": 4}]},
    }
    errs = validate_schema_invariants([bad])
    assert any("pair" in m or "basis" in m for m in errs), f"Expected basis error, got: {errs}"


def test_period_basis_items_must_have_k() -> None:
    """period.basis item missing 'k' triggers an error."""
    bad = {
        "id": "z",
        "period": {"basis": [{"pair": "E-M"}, {"pair": "E-V", "k": 4}]},
    }
    errs = validate_schema_invariants([bad])
    assert any("k" in m or "basis" in m for m in errs), f"Expected basis error, got: {errs}"


def test_valid_period_basis_is_clean() -> None:
    """period.basis items with pair+k are clean."""
    ok = {
        "id": "z",
        "period": {"basis": [{"pair": "E-M", "k": 3}, {"pair": "E-V", "k": 4}]},
    }
    assert validate_schema_invariants([ok]) == []


# ---------------------------------------------------------------------------
# Single-ellipse default: no restrictions
# ---------------------------------------------------------------------------


def test_single_ellipse_with_ae_is_clean() -> None:
    """single-ellipse row with a_au+e is clean (normal v3 case)."""
    ok = {
        "id": "w",
        "cycler_class": "single-ellipse",
        "orbit_elements": {"a_au": 1.5, "e": 0.3},
        "primary": "Sun",
    }
    assert validate_schema_invariants([ok]) == []


def test_no_cycler_class_with_ae_is_clean() -> None:
    """Row with no cycler_class key (defaults single-ellipse) + a/e is clean."""
    ok = {"id": "w", "orbit_elements": {"a_au": 1.5, "e": 0.3}}
    assert validate_schema_invariants([ok]) == []


# ---------------------------------------------------------------------------
# Multiple rows, multiple errors
# ---------------------------------------------------------------------------


def test_multiple_errors_returned() -> None:
    """All errors from all rows are collected (no early exit)."""
    from typing import Any

    rows: list[dict[str, Any]] = [
        {"id": "a", "cycler_class": "multi-arc", "orbit_elements": {"a_au": 1.5}},
        {"id": "b", "cycler_class": "non-keplerian", "primary": "Sun"},
    ]
    errs = validate_schema_invariants(rows)
    assert len(errs) >= 2, f"Expected >=2 errors, got: {errs}"


# ---------------------------------------------------------------------------
# Live catalogue ratchet
# ---------------------------------------------------------------------------


def test_current_catalogue_passes_invariants() -> None:
    """Live data/catalogue.yaml must pass all schema invariants.

    This is a ratchet: it must remain green through Phase 2 (data tagging)
    and beyond. If it fails after cycler_class is backfilled, the tagging
    must be fixed, not this test.
    """
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    errs = validate_schema_invariants(rows)
    assert errs == [], "Live catalogue has schema invariant violations:\n" + "\n".join(errs)
