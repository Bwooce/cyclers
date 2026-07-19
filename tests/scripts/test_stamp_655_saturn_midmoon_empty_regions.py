"""#655 -- regression coverage for the 3 clean Saturn mid-moon-pair empty-region
stamps (Dione-Rhea, Tethys-Dione, Enceladus-Tethys).

Confirms the 3 committed ``data/empty_regions.jsonl`` entries this dispatch
appended are present, pass :func:`validate_empty_region`, and are bounded on
the actual evaluated-candidate counts from the committed #655 enumeration
outputs -- not re-deriving new numbers, just checking round-trip fidelity.
"""

from __future__ import annotations

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    load_empty_regions_list,
    validate_empty_region,
)

_EXPECTED_IDS = {
    "saturn-dione-rhea-symmetric-closure-empty-655",
    "saturn-tethys-dione-symmetric-closure-empty-655",
    "saturn-enceladus-tethys-symmetric-closure-empty-655",
}


def test_all_3_655_stamps_present_and_valid() -> None:
    reports = {r.region_id: r for r in load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)}
    for region_id in _EXPECTED_IDS:
        assert region_id in reports, region_id
        report = reports[region_id]
        validate_empty_region(report)  # must not raise
        assert report.verdict == "EMPTY"
        assert report.result["n_all_gates_passed"] == 0
        assert report.search_extent["points_total"] > 0


def test_rhea_titan_not_stamped() -> None:
    """The companion Rhea-Titan pair (3 gate-passing base closures) must NOT
    be auto-stamped here -- that decision is held for the coordinating
    session per #655's explicit dispatch discipline."""
    reports = {r.region_id: r for r in load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)}
    assert "saturn-rhea-titan-symmetric-closure-empty-655" not in reports
