"""Tests for the isolated-singleton-flip sweep guard (task #559 follow-up)."""

from __future__ import annotations

import pytest

from cyclerfinder.data.sweep_diagnostics import (
    detect_isolated_singleton_anomalies,
    singleton_anomaly_summary,
)


def test_no_anomalies_in_all_pass() -> None:
    assert detect_isolated_singleton_anomalies([True] * 10) == []


def test_no_anomalies_in_smooth_boundary() -> None:
    # A genuine smooth transition (a run of True then a run of False) has
    # no isolated singleton -- every flip has a same-valued neighbor on one
    # side of the transition.
    values = [True, True, True, False, False, False]
    assert detect_isolated_singleton_anomalies(values) == []


def test_single_isolated_flip_detected() -> None:
    values = [True, True, False, True, True]
    anomalies = detect_isolated_singleton_anomalies(values)
    assert len(anomalies) == 1
    assert anomalies[0].index == 2
    assert anomalies[0].value is False


def test_multiple_isolated_flips_detected() -> None:
    values = [True, False, True, True, False, True, True]
    anomalies = detect_isolated_singleton_anomalies(values)
    assert [a.index for a in anomalies] == [1, 4]


def test_a_two_point_dip_is_not_flagged() -> None:
    # Two consecutive FAILs are a cluster, not an isolated singleton --
    # each disagrees with only one neighbor, agreeing with the other.
    values = [True, True, False, False, True, True]
    assert detect_isolated_singleton_anomalies(values) == []


def test_endpoints_never_flagged() -> None:
    # index 0 and the last index have no neighbor on one side, so they
    # can never be flagged regardless of value.
    values = [False, True, True, True]
    anomalies = detect_isolated_singleton_anomalies(values)
    assert all(a.index not in (0, len(values) - 1) for a in anomalies)


def test_labels_carried_through() -> None:
    values = [True, False, True]
    labels = ["2000-01-01", "2000-01-02", "2000-01-03"]
    anomalies = detect_isolated_singleton_anomalies(values, labels)
    assert anomalies[0].label == "2000-01-02"


def test_mismatched_labels_length_raises() -> None:
    with pytest.raises(ValueError, match="labels length"):
        detect_isolated_singleton_anomalies([True, False, True], ["only-one"])


def test_559_actual_pattern_reproduces_known_singleton_count() -> None:
    """Regression pin: replays the #559 finding (isolated FAIL spikes in
    the daily V4-strict sweep on catalogue row #312) against the real
    committed sweep output, so this guard is verified against the case
    that motivated it, not just synthetic examples.

    The 2000 and 2030 windows are DISJOINT (a 30-year gap between them) and
    must be analyzed separately -- concatenating them naively manufactures
    one spurious cross-boundary anomaly at the seam (2030-01-01, sandwiched
    between 2000-12-31 and 2030-01-02), confirmed empirically while writing
    this test. This test doubles as the regression pin for that gotcha,
    documented in this module's own docstring.
    """
    import json
    from pathlib import Path

    path = (
        Path(__file__).resolve().parents[2] / "data" / "silver_327_v4_strict_daily_sweep_559.jsonl"
    )
    if not path.exists():
        pytest.skip("#559 sweep output not present in this checkout")

    rows = []
    with path.open() as f:
        for line in f:
            row = json.loads(line)
            if row.get("kind") == "daily_sweep_row":
                rows.append(row)
    rows.sort(key=lambda r: r["launch_epoch_utc"])
    rows_2000 = [r for r in rows if r["launch_epoch_utc"].startswith("2000")]
    rows_2030 = [r for r in rows if r["launch_epoch_utc"].startswith("2030")]
    assert len(rows_2000) == 366  # leap year
    assert len(rows_2030) == 365

    anomalies_2000 = detect_isolated_singleton_anomalies(
        [bool(r["passes_v4_strict"]) for r in rows_2000],
        [r["launch_epoch_utc"] for r in rows_2000],
    )
    anomalies_2030 = detect_isolated_singleton_anomalies(
        [bool(r["passes_v4_strict"]) for r in rows_2030],
        [r["launch_epoch_utc"] for r in rows_2030],
    )
    assert len(anomalies_2000) == 28
    assert len(anomalies_2030) == 29

    # The naive concatenation (both windows as one sequence) manufactures
    # exactly one extra spurious anomaly at the 2000/2030 seam -- pin that
    # too, so the CAUTION note in the module docstring stays verified.
    all_values = [bool(r["passes_v4_strict"]) for r in rows]
    all_labels = [r["launch_epoch_utc"] for r in rows]
    naive_anomalies = detect_isolated_singleton_anomalies(all_values, all_labels)
    assert len(naive_anomalies) == len(anomalies_2000) + len(anomalies_2030) + 1

    summary = singleton_anomaly_summary(
        [bool(r["passes_v4_strict"]) for r in rows_2000],
        [r["launch_epoch_utc"] for r in rows_2000],
    )
    assert "28/366" in summary
    assert "SUSPECT NUMERICAL ARTIFACT" in summary
