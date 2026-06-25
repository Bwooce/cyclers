"""#450 Task 7: DA/HOTM capability-subsumption stamp on the re-opened negatives.

Once the global multi-rev enumeration lane ships, the four EM/Saturn empty-region
entries whose verbatim re-open key was "a DA/HOTM-class GLOBAL multi-rev discovery
lane" carry a ``reverification`` stamp recording the subsuming method
``da-hotm-enumeration-v1`` and the Png'-recovery outcome. The entries are NOT
deleted (capability-subsumption rule: record the subsumption, don't erase the
prior-method negative).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_REGISTRY = Path(__file__).resolve().parents[2] / "data" / "empty_regions.jsonl"

_SUBSUMED_REGIONS = {
    "cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13",
    "cr3bp-continuation-saturn-mimas-lyapunov-2026-06-12",
    "cr3bp-continuation-saturn-enceladus-lyapunov-2026-06-12",
    "cr3bp-continuation-saturn-tethys-lyapunov-2026-06-12",
}
_METHOD_TAG = "da-hotm-enumeration-v1"


def _rows() -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for line in _REGISTRY.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            rows[r["region_id"]] = r
    return rows


def test_subsumed_regions_carry_da_hotm_stamp() -> None:
    rows = _rows()
    for region_id in _SUBSUMED_REGIONS:
        assert region_id in rows, region_id
        reverifs = rows[region_id].get("reverification", [])
        stamps = [r for r in reverifs if r.get("method") == _METHOD_TAG]
        assert stamps, (region_id, "no da-hotm-enumeration-v1 reverification stamp")
        stamp = stamps[0]
        assert stamp.get("git_sha"), region_id
        assert "result" in stamp, region_id


def test_em_band_records_png_recovery_outcome() -> None:
    """The EM band's stamp records the Png' recovery (re-opened -> Png' as a PO)."""
    rows = _rows()
    em = rows["cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13"]
    stamp = next(r for r in em["reverification"] if r.get("method") == _METHOD_TAG)
    outcome = (stamp.get("result", "") + " " + stamp.get("reason", "")).lower()
    assert "png" in outcome, stamp
    # Re-opened and Png' recovered as a periodic orbit (NOT a cycler/catalogue row).
    assert "po" in outcome or "periodic" in outcome, stamp


def test_entries_not_deleted() -> None:
    """Capability-subsumption: the prior-method negative entries still exist."""
    rows = _rows()
    for region_id in _SUBSUMED_REGIONS:
        assert region_id in rows
        # The original prior-method verdict is preserved (not overwritten).
        assert rows[region_id].get("verdict")
