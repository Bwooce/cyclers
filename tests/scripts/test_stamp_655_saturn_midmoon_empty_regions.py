"""#655 -- regression coverage for all 4 Saturn mid-moon-pair empty-region
stamps (Dione-Rhea, Tethys-Dione, Enceladus-Tethys, Rhea-Titan).

Confirms the 4 committed ``data/empty_regions.jsonl`` entries are present,
pass :func:`validate_empty_region`, and are bounded on the actual
evaluated-candidate counts from the committed #655 enumeration/probe
outputs -- not re-deriving new numbers, just checking round-trip fidelity.

The first 3 pairs are a clean base-gate negative (0/N candidates pass at
all). Rhea-Titan is different: 3/512 candidates DO pass the base gate, but
0/3 survive the inclination-extension + multi-cycle repeat check (the same
qualitative failure #575 found for Titan-Iapetus). The #655 dispatch agent
deliberately did NOT stamp Rhea-Titan itself, holding that disposition
decision for the coordinating session; the coordinator decided (following
direct #575 precedent -- a repeat-check negative is a clean
method-conditional close, not a novelty claim, so no separate
literature-check/adjudication gate applies) to stamp it the same day.
"""

from __future__ import annotations

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    load_empty_regions_list,
    validate_empty_region,
)

_EXPECTED_BASE_GATE_EMPTY_IDS = {
    "saturn-dione-rhea-symmetric-closure-empty-655",
    "saturn-tethys-dione-symmetric-closure-empty-655",
    "saturn-enceladus-tethys-symmetric-closure-empty-655",
}

_RHEA_TITAN_ID = "saturn-rhea-titan-symmetric-closure-empty-655"


def test_all_3_655_base_gate_stamps_present_and_valid() -> None:
    reports = {r.region_id: r for r in load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)}
    for region_id in _EXPECTED_BASE_GATE_EMPTY_IDS:
        assert region_id in reports, region_id
        report = reports[region_id]
        validate_empty_region(report)  # must not raise
        assert report.verdict == "EMPTY"
        assert report.result["n_all_gates_passed"] == 0
        assert report.search_extent["points_total"] > 0


def test_rhea_titan_stamped_by_coordinator_decision() -> None:
    """Rhea-Titan's 3 base-gate survivors all failed the repeat-check
    extension (0/3 survive), so the coordinator stamped it as a clean
    method-conditional negative -- same disposition as #575's own closure,
    no adjudication needed since nothing novel is being claimed."""
    reports = {r.region_id: r for r in load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)}
    assert _RHEA_TITAN_ID in reports
    report = reports[_RHEA_TITAN_ID]
    validate_empty_region(report)  # must not raise
    assert report.result["n_base_gate_survivors"] == 3
    assert report.result["n_repeat_check_survivors"] == 0
    assert report.search_extent["points_total"] > 0
