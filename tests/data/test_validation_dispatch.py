"""Phase 6.1: anchors_for() dispatch helper and construction guard tests.

Asserts that:
* ``anchors_for`` returns the correct anchor categories for each
  ``cycler_class`` (spec §16.7.5).
* multi-arc entries include ``invariants`` and ``cr3bp`` is False;
  ``a_e`` is False (single-ellipse constructor must not be applied).
* single-ellipse entries include ``a_e``; ``invariants`` and ``cr3bp``
  are False.
* non-keplerian entries include ``cr3bp``; ``a_e``, ``invariants``, and
  ``period`` are False.
* The loader's ``CONSTRUCTIBLE`` set never includes a non-single-ellipse
  entry — this is the gauntlet guard: the single-ellipse constructor
  (``_build_cell_from_entry``) is only ever called on entries for which
  ``anchors_for``'s ``a_e`` is True.
"""

from __future__ import annotations

from typing import Any

from cyclerfinder.data.validate import anchors_for
from tests._catalogue_loader import ExclusionReason, classify_catalogue, classify_row

# ---------------------------------------------------------------------------
# Unit tests for anchors_for() with synthetic rows
# ---------------------------------------------------------------------------


def test_single_ellipse_includes_a_e() -> None:
    """single-ellipse → a_e anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "single-ellipse"}
    a = anchors_for(row)
    assert a["a_e"] is True


def test_single_ellipse_excludes_invariants() -> None:
    """single-ellipse → invariants anchor is False."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "single-ellipse"}
    a = anchors_for(row)
    assert a["invariants"] is False


def test_single_ellipse_excludes_cr3bp() -> None:
    """single-ellipse → cr3bp anchor is False."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "single-ellipse"}
    a = anchors_for(row)
    assert a["cr3bp"] is False


def test_single_ellipse_includes_period() -> None:
    """single-ellipse → period anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "single-ellipse"}
    a = anchors_for(row)
    assert a["period"] is True


def test_single_ellipse_includes_vinf() -> None:
    """single-ellipse → vinf anchor is True (universal)."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "single-ellipse"}
    a = anchors_for(row)
    assert a["vinf"] is True


def test_multi_arc_includes_invariants() -> None:
    """multi-arc → invariants anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "multi-arc"}
    a = anchors_for(row)
    assert a["invariants"] is True


def test_multi_arc_excludes_a_e() -> None:
    """multi-arc → a_e anchor is False (no single conic; constructor must not run)."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "multi-arc"}
    a = anchors_for(row)
    assert a["a_e"] is False


def test_multi_arc_excludes_cr3bp() -> None:
    """multi-arc → cr3bp anchor is False."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "multi-arc"}
    a = anchors_for(row)
    assert a["cr3bp"] is False


def test_multi_arc_includes_period() -> None:
    """multi-arc → period anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "multi-arc"}
    a = anchors_for(row)
    assert a["period"] is True


def test_multi_arc_includes_vinf() -> None:
    """multi-arc → vinf anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "multi-arc"}
    a = anchors_for(row)
    assert a["vinf"] is True


def test_non_keplerian_includes_cr3bp() -> None:
    """non-keplerian → cr3bp anchor is True."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "non-keplerian"}
    a = anchors_for(row)
    assert a["cr3bp"] is True


def test_non_keplerian_excludes_a_e() -> None:
    """non-keplerian → a_e anchor is False."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "non-keplerian"}
    a = anchors_for(row)
    assert a["a_e"] is False


def test_non_keplerian_excludes_invariants() -> None:
    """non-keplerian → invariants anchor is False."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "non-keplerian"}
    a = anchors_for(row)
    assert a["invariants"] is False


def test_non_keplerian_excludes_period() -> None:
    """non-keplerian → period anchor is False (CR3BP period is non-dimensional)."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "non-keplerian"}
    a = anchors_for(row)
    assert a["period"] is False


def test_non_keplerian_includes_vinf() -> None:
    """non-keplerian → vinf anchor is True (universal; may be null in data)."""
    row: dict[str, Any] = {"id": "x", "cycler_class": "non-keplerian"}
    a = anchors_for(row)
    assert a["vinf"] is True


def test_default_class_is_single_ellipse() -> None:
    """Row with no cycler_class key defaults to single-ellipse behaviour."""
    row: dict[str, Any] = {"id": "x"}
    a = anchors_for(row)
    assert a["a_e"] is True
    assert a["invariants"] is False
    assert a["cr3bp"] is False
    assert a["period"] is True


# ---------------------------------------------------------------------------
# Live catalogue spot-checks — specific known rows
# ---------------------------------------------------------------------------


def _row_by_id(rows: list[dict[str, Any]], rid: str) -> dict[str, Any]:
    for r in rows:
        if r.get("id") == rid:
            return r
    raise KeyError(f"Row {rid!r} not found in catalogue")


def test_mcconaghy_2006_multi_arc_dispatch(catalogue_rows: list[dict[str, Any]]) -> None:
    """mcconaghy-2006-em-k2 is multi-arc → invariants True, a_e False."""
    rows = catalogue_rows
    row = _row_by_id(rows, "mcconaghy-2006-em-k2")
    a = anchors_for(row)
    assert a["invariants"] is True
    assert a["a_e"] is False
    assert a["cr3bp"] is False


def test_aldrin_single_ellipse_dispatch(catalogue_rows: list[dict[str, Any]]) -> None:
    """aldrin-classic-em-k1-outbound is single-ellipse → a_e True, invariants False."""
    rows = catalogue_rows
    row = _row_by_id(rows, "aldrin-classic-em-k1-outbound")
    a = anchors_for(row)
    assert a["a_e"] is True
    assert a["invariants"] is False
    assert a["cr3bp"] is False


def test_arenstorf_non_keplerian_dispatch(catalogue_rows: list[dict[str, Any]]) -> None:
    """arenstorf-em-figure8-1963 is non-keplerian → cr3bp True, a_e False, period False."""
    rows = catalogue_rows
    row = _row_by_id(rows, "arenstorf-em-figure8-1963")
    a = anchors_for(row)
    assert a["cr3bp"] is True
    assert a["a_e"] is False
    assert a["period"] is False
    assert a["invariants"] is False


# ---------------------------------------------------------------------------
# Construction guard: no CONSTRUCTIBLE entry may be non-single-ellipse
#
# The v1 gauntlet calls _build_cell_from_entry (the single-ellipse
# resonant-constructor path) on every CONSTRUCTIBLE entry.  If a
# multi-arc or non-keplerian entry were ever classified as CONSTRUCTIBLE
# by the loader, it would be fed to the wrong constructor.  This test
# asserts that invariant: every CONSTRUCTIBLE entry must have
# anchors_for(row)["a_e"] == True.
# ---------------------------------------------------------------------------


def test_all_constructible_entries_are_single_ellipse(
    catalogue_rows: list[dict[str, Any]],
) -> None:
    """Every CONSTRUCTIBLE entry must have cycler_class == single-ellipse
    (anchors_for a_e == True) — the v1 single-ellipse constructor must
    only be called for single-ellipse entries.

    If this fails, a non-single-ellipse row leaked into the gauntlet and
    would be fed to the wrong constructor.  Fix the loader's guard, not
    this test.
    """
    rows = catalogue_rows
    row_by_id = {r["id"]: r for r in rows}
    violations: list[str] = []
    for cid, reason in classify_catalogue():
        if reason is not ExclusionReason.CONSTRUCTIBLE:
            continue
        row = row_by_id.get(cid, {})
        if not anchors_for(row)["a_e"]:
            cls = row.get("cycler_class", "single-ellipse")
            violations.append(f"{cid}: cycler_class={cls!r}, but classified CONSTRUCTIBLE")
    assert violations == [], (
        "Non-single-ellipse entries classified as CONSTRUCTIBLE — "
        "the single-ellipse constructor would be misapplied:\n" + "\n".join(violations)
    )


def test_classify_row_excludes_multi_arc_by_sequence(
    catalogue_rows: list[dict[str, Any]],
) -> None:
    """classify_row excludes multi-arc entries via MULTI_ENCOUNTER_SEQUENCE
    (they always have >2-encounter sequences), never via cycler_class directly.

    Confirms the structural guard is in place: even without an explicit
    cycler_class check, no multi-arc row can slip into CONSTRUCTIBLE because
    its sequence is always a multi-encounter (e.g. E-E-M-M) form.
    """
    rows = catalogue_rows
    multi_arc_constructible: list[str] = []
    for row in rows:
        if row.get("cycler_class") != "multi-arc":
            continue
        reason, _ = classify_row(row)
        if reason is ExclusionReason.CONSTRUCTIBLE:
            multi_arc_constructible.append(row["id"])
    assert multi_arc_constructible == [], (
        f"multi-arc entries classified as CONSTRUCTIBLE (constructor would be misapplied): "
        f"{multi_arc_constructible}"
    )
