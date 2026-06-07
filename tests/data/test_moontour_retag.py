"""Tier-1 Phase 2: Jovian patched-conic rows route as multi-arc, not CR3BP
(plan Phase 2 Task 2.1). The gauntlet anchors_for dispatch keys off cycler_class
(validate.py anchors_for); a non-keplerian tag demands the cr3bp identity triple
from a patched-conic row, which is the misroute the design (§7) flags."""

from __future__ import annotations

import pytest
import yaml

from cyclerfinder.data.validate import anchors_for
from tests._catalogue_loader import CATALOGUE_PATH

_JOVIAN = (
    "hernandez-2017-jovian-ieg-triple-family",
    "russell-strange-2009-jovian-multimoon-family",
)


def _row(entry_id: str) -> dict:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row["id"] == entry_id:
            return row
    raise AssertionError(f"catalogue row {entry_id!r} not found")


@pytest.mark.parametrize("entry_id", _JOVIAN)
def test_jovian_rows_are_multi_arc(entry_id: str) -> None:
    assert _row(entry_id)["cycler_class"] == "multi-arc"


@pytest.mark.parametrize("entry_id", _JOVIAN)
def test_jovian_dispatch_wants_invariants_not_cr3bp(entry_id: str) -> None:
    a = anchors_for(_row(entry_id))
    assert a["invariants"] is True
    assert a["cr3bp"] is False


def test_saturnian_row_stays_non_keplerian_with_titan_split_note() -> None:
    row = _row("russell-strange-2009-saturnian-multimoon-family")
    # NOT silently re-tagged: midsize members are genuinely CR3BP (Tier 2).
    assert row["cycler_class"] == "non-keplerian"
    notes = row.get("notes") or ""
    assert "Titan" in notes
    # The honest split must name the Tier-1/Tier-2 boundary.
    assert "Tier-1" in notes or "patched-conic" in notes
    assert "Tier-2" in notes or "CR3BP" in notes
