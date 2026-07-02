"""Tests for the #521 phase-2 mandatory pre-flight search gate."""

from __future__ import annotations

from pathlib import Path

import pytest

from cyclerfinder.data.empty_regions import EmptyRegionReport
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import (
    LARGE_GRID_THRESHOLD,
    PreflightBlockedError,
    PreflightResult,
    preflight_search,
)

_WEAK_METHOD = MethodCapability(
    genome="single-ellipse free-return",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "coplanar"}),
    git_sha="abc1234",
)

_STRONG_METHOD = MethodCapability(
    genome="single-ellipse free-return",
    corrector="ballistic_correct",
    capability_tags=frozenset({"powered", "coplanar"}),  # powered subsumes ballistic
    git_sha="def5678",
)


def _outstanding_with(*task_nos: int) -> str:
    lines = ["# Outstanding\n", "\n", "**TASK ALLOCATIONS:**\n"]
    for n in task_nos:
        lines.append(f"- **#{n}** — Some task description.\n")
    return "".join(lines)


def _empty_region(region_id: str, method: MethodCapability) -> EmptyRegionReport:
    return EmptyRegionReport(
        region_id=region_id,
        family="test-family",
        centre="Earth",
        topologies=(),
        method_capability=method,
        search_extent={"points_total": 10},
        prune_gates=("test-gate",),
        result={},
        verdict="EMPTY",
        interpretation="test fixture",
        source_anchors="",
        run={},
    )


def _call(
    tmp_path: Path,
    *,
    task_no: int = 999,
    region_id: str = "test-region",
    method: MethodCapability = _WEAK_METHOD,
    script_name: str = "run_999_something.py",
    n_points: int = 10,
    timing_pilot_seconds_per_point: float | None = None,
    override_reason: str | None = None,
    allocated_task_nos: tuple[int, ...] = (999,),
    registry: tuple[EmptyRegionReport, ...] = (),
) -> PreflightResult:
    outstanding_path = tmp_path / "OUTSTANDING.md"
    outstanding_path.write_text(_outstanding_with(*allocated_task_nos))
    runlog_path = tmp_path / "runlogs" / "preflight_runlog.jsonl"
    return preflight_search(
        task_no=task_no,
        region_id=region_id,
        method=method,
        script_path=tmp_path / "scripts" / script_name,
        n_points=n_points,
        timing_pilot_seconds_per_point=timing_pilot_seconds_per_point,
        override_reason=override_reason,
        outstanding_path=outstanding_path,
        registry=registry,
        runlog_path=runlog_path,
    )


def test_allows_registered_task_small_grid_open_region(tmp_path: Path) -> None:
    result = _call(tmp_path)
    assert result.proceed is True
    assert result.warnings == ()


def test_blocks_unregistered_task_number(tmp_path: Path) -> None:
    with pytest.raises(PreflightBlockedError, match="not recorded"):
        _call(tmp_path, task_no=42, allocated_task_nos=(999,))


def test_blocks_filename_task_number_mismatch(tmp_path: Path) -> None:
    with pytest.raises(PreflightBlockedError, match="filename and declared task"):
        _call(
            tmp_path,
            task_no=888,
            script_name="run_777_something.py",
            allocated_task_nos=(888, 777),
        )


def test_allows_matching_filename_task_number(tmp_path: Path) -> None:
    result = _call(
        tmp_path,
        task_no=777,
        script_name="run_777_thing.py",
        allocated_task_nos=(777,),
    )
    assert result.proceed is True


def test_scripts_without_run_nnn_filename_skip_the_filename_check(tmp_path: Path) -> None:
    """A script not following the run_NNN_*.py convention isn't penalised."""
    result = _call(tmp_path, task_no=999, script_name="ad_hoc_helper.py")
    assert result.proceed is True


def test_blocks_subsumed_region(tmp_path: Path) -> None:
    prior = _empty_region("cross-system-region", _STRONG_METHOD)
    with pytest.raises(PreflightBlockedError, match="already covered"):
        _call(
            tmp_path,
            region_id="cross-system-region",
            method=_WEAK_METHOD,
            registry=(prior,),
        )


def test_allows_stronger_method_over_previously_swept_region(tmp_path: Path) -> None:
    prior = _empty_region("cross-system-region", _WEAK_METHOD)
    result = _call(
        tmp_path,
        region_id="cross-system-region",
        method=_STRONG_METHOD,
        registry=(prior,),
    )
    assert result.proceed is True


def test_allows_unrelated_region_regardless_of_registry(tmp_path: Path) -> None:
    prior = _empty_region("some-other-region", _STRONG_METHOD)
    result = _call(
        tmp_path,
        region_id="never-swept-region",
        method=_WEAK_METHOD,
        registry=(prior,),
    )
    assert result.proceed is True


def test_blocks_large_grid_without_timing_pilot(tmp_path: Path) -> None:
    with pytest.raises(PreflightBlockedError, match="LARGE_GRID_THRESHOLD"):
        _call(tmp_path, n_points=LARGE_GRID_THRESHOLD + 1)


def test_allows_large_grid_with_timing_pilot(tmp_path: Path) -> None:
    result = _call(
        tmp_path,
        n_points=LARGE_GRID_THRESHOLD + 1,
        timing_pilot_seconds_per_point=1.0,
    )
    assert result.proceed is True


def test_small_grid_does_not_require_timing_pilot(tmp_path: Path) -> None:
    result = _call(tmp_path, n_points=LARGE_GRID_THRESHOLD)
    assert result.proceed is True


def test_override_reason_downgrades_failure_to_warning(tmp_path: Path) -> None:
    result = _call(
        tmp_path,
        task_no=42,
        script_name="run_42_something.py",
        allocated_task_nos=(999,),
        override_reason="testing the override path",
    )
    assert result.proceed is True
    assert result.override_reason == "testing the override path"
    assert len(result.warnings) == 1
    assert "not recorded" in result.warnings[0]


def test_blocked_exception_lists_all_failing_checks(tmp_path: Path) -> None:
    with pytest.raises(PreflightBlockedError) as exc_info:
        _call(
            tmp_path,
            task_no=42,
            allocated_task_nos=(999,),
            n_points=LARGE_GRID_THRESHOLD + 1,
        )
    message = str(exc_info.value)
    assert "not recorded" in message
    assert "LARGE_GRID_THRESHOLD" in message


def test_runlog_records_every_invocation(tmp_path: Path) -> None:
    runlog_path = tmp_path / "runlogs" / "preflight_runlog.jsonl"
    outstanding_path = tmp_path / "OUTSTANDING.md"
    outstanding_path.write_text(_outstanding_with(999))

    preflight_search(
        task_no=999,
        region_id="region-a",
        method=_WEAK_METHOD,
        script_path=tmp_path / "scripts" / "run_999_a.py",
        n_points=10,
        outstanding_path=outstanding_path,
        registry=(),
        runlog_path=runlog_path,
    )
    with pytest.raises(PreflightBlockedError):
        preflight_search(
            task_no=42,
            region_id="region-b",
            method=_WEAK_METHOD,
            script_path=tmp_path / "scripts" / "run_999_b.py",
            n_points=10,
            outstanding_path=outstanding_path,
            registry=(),
            runlog_path=runlog_path,
        )

    import json

    lines = runlog_path.read_text().strip().splitlines()
    assert len(lines) == 2
    first, second = (json.loads(line) for line in lines)
    assert first["task_no"] == 999
    assert first["proceed"] is True
    assert second["task_no"] == 42
    assert second["proceed"] is False
    assert second["failures"]
