"""Phase 6 Phase 4: empty-region report is bounded + reproducible (Tasks 4.0/4.1)."""

from __future__ import annotations

from pathlib import Path

import pytest

from cyclerfinder.data.empty_regions import (
    EmptyRegionReport,
    append_empty_region,
    is_catalogue_source,
    load_empty_regions,
    validate_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability

_METHOD = MethodCapability(
    genome="single-ellipse free-return",
    corrector="ballistic_correct (no-leveraging)",
    capability_tags=frozenset({"ballistic", "patched-conic", "single-arc", "coplanar"}),
    git_sha="deadbeef",
)


def _valid_report(region_id: str = "x") -> EmptyRegionReport:
    return EmptyRegionReport(
        region_id=region_id,
        family="planet-centric moon system (Jupiter)",
        centre="Jupiter",
        topologies=({"sequence": ["Io", "Europa", "Ganymede", "Io"], "period_k": 1},),
        method_capability=_METHOD,
        search_extent={"points_total": 2816, "n_epochs": 256, "center": "Jupiter"},
        prune_gates=("vilm_dv_floor<=budget", "linkable(Jovicentric)", "max_bend_deg"),
        result={"closed": 0, "bend_feasible": 0, "best_max_vinf_kms": 10.4},
        verdict="EMPTY",
        interpretation="no bend-feasible closure below the V_inf floor",
        source_anchors="none populated in (circular-coplanar, Jupiter) bucket",
        run={"date": "2026-06-08", "git_sha": "deadbeef", "wall_s": 12.0},
    )


def test_empty_region_requires_bounded_search_extent() -> None:
    r = EmptyRegionReport(
        region_id="x",
        family="f",
        centre="Jupiter",
        topologies=(),
        method_capability=_METHOD,
        search_extent={"points_total": 0},
        prune_gates=("vilm",),
        result={},
        verdict="EMPTY",
        interpretation="",
        source_anchors="",
        run={},
    )
    with pytest.raises(ValueError):
        validate_empty_region(r)  # points_total == 0 -> unbounded -> invalid


def test_empty_region_requires_prune_gates() -> None:
    r = EmptyRegionReport(
        region_id="x",
        family="f",
        centre="Jupiter",
        topologies=(),
        method_capability=_METHOD,
        search_extent={"points_total": 2816},
        prune_gates=(),
        result={},
        verdict="EMPTY",
        interpretation="",
        source_anchors="",
        run={},
    )
    with pytest.raises(ValueError):
        validate_empty_region(r)  # no prune gates -> can't bound over-pruning


def test_empty_region_requires_method_capability() -> None:
    bad_method = MethodCapability(
        genome="x", corrector="y", capability_tags=frozenset(), git_sha="z"
    )
    r = EmptyRegionReport(
        region_id="x",
        family="f",
        centre="Jupiter",
        topologies=(),
        method_capability=bad_method,
        search_extent={"points_total": 2816},
        prune_gates=("vilm",),
        result={},
        verdict="EMPTY",
        interpretation="",
        source_anchors="",
        run={},
    )
    with pytest.raises(ValueError):
        validate_empty_region(r)  # empty capability tags -> unconditional "empty"


def test_valid_report_passes_validation() -> None:
    validate_empty_region(_valid_report())  # no raise


def test_is_catalogue_source_is_false() -> None:
    assert is_catalogue_source() is False


def test_append_and_load_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "empty_regions.jsonl"
    report = _valid_report(region_id="round-trip")
    append_empty_region(path, report)
    loaded = list(load_empty_regions(path))
    assert len(loaded) == 1
    got = loaded[0]
    assert got.region_id == "round-trip"
    assert got.centre == "Jupiter"
    assert got.method_capability == _METHOD
    assert got.prune_gates == report.prune_gates
    assert got.search_extent == report.search_extent
