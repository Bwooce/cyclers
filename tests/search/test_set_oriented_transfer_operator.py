"""Tests for :mod:`cyclerfinder.search.set_oriented_transfer_operator` (task #664).

Model-agnostic generic GAIO-style box/transfer-operator machinery -- tested
here with small synthetic maps (no CR3BP dependency); the Sun-Jupiter
Dellnitz-2005 positive control is tested separately in
``test_quasi_hilda_positive_control.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.search.set_oriented_transfer_operator import (
    BoxGrid,
    almost_invariance_ratio,
    almost_invariant_sets_spectral,
    build_transition_matrix,
    transport_probability,
)


def _grid_1d(n: int, lo: float = 0.0, hi: float = 1.0) -> BoxGrid:
    return BoxGrid(lower=np.array([lo]), upper=np.array([hi]), counts=(n,))


def _grid_2d(nx: int, ny: int) -> BoxGrid:
    return BoxGrid(lower=np.array([0.0, 0.0]), upper=np.array([1.0, 1.0]), counts=(nx, ny))


# ---------------------------------------------------------------------------
# BoxGrid indexing
# ---------------------------------------------------------------------------


def test_box_grid_rejects_bad_bounds() -> None:
    with pytest.raises(ValueError, match="upper must exceed lower"):
        BoxGrid(lower=np.array([1.0]), upper=np.array([0.0]), counts=(4,))


def test_box_grid_round_trip_flat_multi_index() -> None:
    grid = _grid_2d(5, 7)
    for flat in range(grid.n_boxes):
        multi = grid.flat_to_multi(flat)
        assert grid.multi_to_flat(multi) == flat


def test_box_grid_index_of_point_matches_bounds() -> None:
    grid = _grid_1d(10)
    for flat in range(grid.n_boxes):
        lo, hi = grid.box_bounds(flat)
        center = (lo + hi) / 2.0
        assert grid.index_of_point(center) == flat
    # Outside the domain -> None.
    assert grid.index_of_point(np.array([-0.5])) is None
    assert grid.index_of_point(np.array([1.5])) is None
    # Exactly on the upper bound is exclusive (half-open boxes).
    assert grid.index_of_point(np.array([1.0])) is None


def test_box_grid_sample_uniform_stays_within_box() -> None:
    grid = _grid_2d(4, 4)
    rng = np.random.default_rng(0)
    for flat in range(grid.n_boxes):
        pts = grid.sample_uniform(flat, 50, rng)
        lo, hi = grid.box_bounds(flat)
        assert np.all(pts >= lo)
        assert np.all(pts <= hi)


# ---------------------------------------------------------------------------
# Transition matrix construction
# ---------------------------------------------------------------------------


def test_identity_map_gives_near_identity_transition_matrix() -> None:
    """A map that returns its input unchanged should transition every box
    almost entirely to itself, with ~zero escape."""
    grid = _grid_1d(8)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(1)

    def identity_map(pt: np.ndarray) -> np.ndarray:
        return pt

    result = build_transition_matrix(grid, mask, identity_map, n_samples_per_box=20, rng=rng)
    diag = result.p.diagonal()
    assert np.all(diag > 0.99)
    assert np.allclose(result.escape_fraction, 0.0)
    # Off-diagonal mass should be ~0 (identity map never leaves its own box).
    off_diag_mass = result.p.sum() - diag.sum()
    assert off_diag_mass < 1e-9


def test_always_escaping_map_gives_full_escape_fraction() -> None:
    grid = _grid_1d(6)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(2)

    def escaping_map(pt: np.ndarray) -> np.ndarray | None:
        return None

    result = build_transition_matrix(grid, mask, escaping_map, n_samples_per_box=10, rng=rng)
    assert np.allclose(result.escape_fraction, 1.0)
    assert result.p.nnz == 0


def test_masked_out_boxes_are_never_sampled() -> None:
    grid = _grid_1d(6)
    mask = np.array([True, True, False, False, True, True])
    rng = np.random.default_rng(3)
    calls: list[np.ndarray] = []

    def recording_map(pt: np.ndarray) -> np.ndarray:
        calls.append(pt.copy())
        return pt

    result = build_transition_matrix(grid, mask, recording_map, n_samples_per_box=5, rng=rng)
    assert result.n_samples[2] == 0
    assert result.n_samples[3] == 0
    assert result.n_samples[0] == 5
    for pt in calls:
        idx = grid.index_of_point(pt)
        assert idx is not None
        assert mask[idx]


def test_shift_map_transitions_mass_to_neighbor_box() -> None:
    """A deterministic shift-by-one-box-width map should send ~all of box i's
    mass to box i+1 (except the last box, which escapes off the right edge)."""
    n = 10
    grid = _grid_1d(n)
    mask = np.ones(n, dtype=bool)
    rng = np.random.default_rng(4)
    width = 1.0 / n

    def shift_map(pt: np.ndarray) -> np.ndarray | None:
        new = pt + width
        if new[0] >= 1.0:
            return None
        return new

    result = build_transition_matrix(grid, mask, shift_map, n_samples_per_box=30, rng=rng)
    for i in range(n - 1):
        assert result.p[i, i + 1] > 0.95
    assert result.escape_fraction[n - 1] > 0.95


# ---------------------------------------------------------------------------
# Spectral almost-invariant set extraction
# ---------------------------------------------------------------------------


def test_two_disconnected_clusters_recovered_by_spectral_decomposition() -> None:
    """Two blocks of boxes that only transition among themselves (never
    across) is the textbook case: the top eigenvector should be near-
    constant-sign (invariant density) and the second eigenvector's sign
    should split the two blocks perfectly."""
    grid = _grid_1d(20)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(5)

    def block_map(pt: np.ndarray) -> np.ndarray:
        # Stay within [0, 0.5) or [0.5, 1.0) -- never cross the midline.
        idx = grid.index_of_point(pt)
        assert idx is not None
        if idx < 10:
            return rng.uniform(0.0, 0.5 - 1e-9, size=1)
        return rng.uniform(0.5, 1.0 - 1e-9, size=1)

    result = build_transition_matrix(grid, mask, block_map, n_samples_per_box=40, rng=rng)
    decomp = almost_invariant_sets_spectral(result, n_eigvecs=2)
    # Every box in [0,10) shares one label; every box in [10,20) shares the other.
    labels_low = decomp.labels[:10]
    labels_high = decomp.labels[10:]
    assert len(set(labels_low.tolist())) == 1
    assert len(set(labels_high.tolist())) == 1
    assert labels_low[0] != labels_high[0]
    # Leading eigenvalue should be close to 1 (this P is exactly row-stochastic,
    # closed system, no escape).
    assert abs(decomp.eigenvalues[0].real - 1.0) < 1e-6


def test_almost_invariance_ratio_high_for_recovered_cluster_low_for_random_mix() -> None:
    grid = _grid_1d(20)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(6)

    def block_map(pt: np.ndarray) -> np.ndarray:
        idx = grid.index_of_point(pt)
        assert idx is not None
        if idx < 10:
            return rng.uniform(0.0, 0.5 - 1e-9, size=1)
        return rng.uniform(0.5, 1.0 - 1e-9, size=1)

    result = build_transition_matrix(grid, mask, block_map, n_samples_per_box=40, rng=rng)
    ratio_true_cluster = almost_invariance_ratio(result, np.arange(10))
    ratio_random_mix = almost_invariance_ratio(
        result, np.array([0, 1, 2, 3, 4, 15, 16, 17, 18, 19])
    )
    assert ratio_true_cluster > 0.99
    assert ratio_random_mix < 0.6


def test_almost_invariance_ratio_empty_set_is_zero() -> None:
    grid = _grid_1d(4)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(7)
    result = build_transition_matrix(grid, mask, lambda pt: pt, n_samples_per_box=5, rng=rng)
    assert almost_invariance_ratio(result, np.array([], dtype=np.intp)) == 0.0


# ---------------------------------------------------------------------------
# Transport probability
# ---------------------------------------------------------------------------


def test_transport_probability_zero_when_disjoint_and_no_mixing() -> None:
    grid = _grid_1d(20)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(8)

    def block_map(pt: np.ndarray) -> np.ndarray:
        idx = grid.index_of_point(pt)
        assert idx is not None
        if idx < 10:
            return rng.uniform(0.0, 0.5 - 1e-9, size=1)
        return rng.uniform(0.5, 1.0 - 1e-9, size=1)

    result = build_transition_matrix(grid, mask, block_map, n_samples_per_box=40, rng=rng)
    source = np.arange(10)
    target = np.arange(10, 20)
    probs = transport_probability(result, source, target, n_iters=5)
    assert np.allclose(probs, 0.0)


def test_transport_probability_escapes_to_zero_with_leaky_map() -> None:
    """A map that always escapes should drain ALL probability out of the
    covered domain after one iterate -- transport to any target is 0 for
    n >= 1 (and the source/target overlap at n=0 is whatever it is)."""
    grid = _grid_1d(10)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(9)
    result = build_transition_matrix(grid, mask, lambda pt: None, n_samples_per_box=5, rng=rng)
    probs = transport_probability(result, np.arange(5), np.arange(5, 10), n_iters=3)
    assert np.allclose(probs[1:], 0.0)


def test_transport_probability_requires_nonempty_source() -> None:
    grid = _grid_1d(4)
    mask = np.ones(grid.n_boxes, dtype=bool)
    rng = np.random.default_rng(10)
    result = build_transition_matrix(grid, mask, lambda pt: pt, n_samples_per_box=3, rng=rng)
    with pytest.raises(ValueError, match="source_boxes is empty"):
        transport_probability(result, np.array([], dtype=np.intp), np.array([0]), n_iters=2)
