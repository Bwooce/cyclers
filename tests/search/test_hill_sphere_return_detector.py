"""Tests for the #535 repeated-Hill-sphere-return detector.

POSITIVE CONTROL: a hand-constructed, closed-form synthetic distance-vs-time
series with a KNOWN dip pattern (exact dip times/depths/spacing) -- the
EXPECTED return count/timing is a property of the construction itself, not a
value this module computed and then asserted against itself.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.hill_sphere_return_detector import (
    find_admission_windows,
    find_returns,
)


def _dip_series(
    dip_times: list[float],
    *,
    dip_depth: float,
    baseline: float,
    dip_half_width: float,
    t_max: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """A synthetic 1D "distance" trajectory (embedded as (x, 0) 2D positions)
    with Gaussian dips of ``dip_depth`` centered at each of ``dip_times``.
    """
    t = np.arange(0.0, t_max, dt)
    dist = np.full_like(t, baseline)
    for tc in dip_times:
        dist = dist - (baseline - dip_depth) * np.exp(-0.5 * ((t - tc) / dip_half_width) ** 2)
    positions = np.stack([dist, np.zeros_like(dist)], axis=1)
    return t, positions


def test_find_returns_counts_well_separated_dips() -> None:
    """5 dips spaced 6 years apart, each well below r_hill=1.0, baseline well
    above it -- exactly 5 distinct returns expected, no merging.
    """
    dip_times = [0.0, 6.0, 12.0, 18.0, 24.0]
    t, pos = _dip_series(
        dip_times, dip_depth=0.3, baseline=5.0, dip_half_width=0.1, t_max=30.0, dt=0.01
    )
    returns = find_returns(t, pos, r_hill=1.0, min_separation=1.0)
    assert len(returns) == 5
    for r, expected_tc in zip(returns, dip_times, strict=True):
        assert abs(r.t_closest - expected_tc) < 0.05
        assert abs(r.closest_distance - 0.3) < 0.05


def test_find_returns_merges_closely_spaced_dips() -> None:
    """Two dips only 0.3 years apart (below the 1-year min_separation) must
    merge into ONE return, not two.
    """
    dip_times = [10.0, 10.3]
    t, pos = _dip_series(
        dip_times, dip_depth=0.3, baseline=5.0, dip_half_width=0.05, t_max=30.0, dt=0.005
    )
    returns = find_returns(t, pos, r_hill=1.0, min_separation=1.0)
    assert len(returns) == 1


def test_find_returns_does_not_merge_well_separated_dips() -> None:
    """Two dips 2 years apart (above the 1-year floor) must stay distinct."""
    dip_times = [10.0, 12.0]
    t, pos = _dip_series(
        dip_times, dip_depth=0.3, baseline=5.0, dip_half_width=0.05, t_max=30.0, dt=0.01
    )
    returns = find_returns(t, pos, r_hill=1.0, min_separation=1.0)
    assert len(returns) == 2


def test_admission_window_finds_a_qualifying_window() -> None:
    """5 returns spaced 6 years apart over 30 years -- a 15-year window
    starting near t=0 should contain exactly 3 returns (0, 6, 12), an
    admissible count under n_returns in [3,15].
    """
    dip_times = [0.0, 6.0, 12.0, 18.0, 24.0]
    t, pos = _dip_series(
        dip_times, dip_depth=0.3, baseline=5.0, dip_half_width=0.1, t_max=30.0, dt=0.01
    )
    returns = find_returns(t, pos, r_hill=1.0, min_separation=1.0)
    windows = find_admission_windows(
        returns,
        float(t[0]),
        float(t[-1]),
        window_lo=10.0,
        window_hi=15.0,
        n_returns_lo=3,
        n_returns_hi=15,
        geometry_factor=3.0,
    )
    assert len(windows) > 0
    counts = [len(w.returns) for w in windows]
    assert any(3 <= c <= 15 for c in counts)
    # Every reported window is genuinely admissible (the function's own contract).
    for w in windows:
        assert 3 <= len(w.returns) <= 15


def test_admission_window_geometry_check_flags_uneven_approaches() -> None:
    """If one return is much closer than the others (ratio > geometry_factor),
    geometry_ok must be False for any window containing both.
    """
    dip_times = [0.0, 6.0, 12.0]
    t = np.arange(0.0, 20.0, 0.01)
    dist = np.full_like(t, 5.0)
    depths = [0.3, 0.3, 0.03]  # last dip 10x closer than the others
    for tc, depth in zip(dip_times, depths, strict=True):
        dist = dist - (5.0 - depth) * np.exp(-0.5 * ((t - tc) / 0.1) ** 2)
    pos = np.stack([dist, np.zeros_like(dist)], axis=1)
    returns = find_returns(t, pos, r_hill=1.0, min_separation=1.0)
    assert len(returns) == 3
    windows = find_admission_windows(
        returns,
        float(t[0]),
        float(t[-1]),
        window_lo=15.0,
        window_hi=15.0,
        n_returns_lo=3,
        n_returns_hi=15,
        geometry_factor=3.0,
    )
    full_windows = [w for w in windows if len(w.returns) == 3]
    assert len(full_windows) > 0
    for w in full_windows:
        assert not w.geometry_ok
        assert w.geometry_ratio > 3.0
