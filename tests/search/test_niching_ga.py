"""Deterministic-crowding niching GA tests (#581 stage 1).

The load-bearing property (Gurfil-Kasdin 2002 p. 5687; Mahfoud's Deterministic
Crowding): with two well-separated optima of EQUAL height in a toy landscape,
both niches survive to the final generation of a single run — the mechanism a
plain single-optimum evolutionary loop (every existing DE lane in this repo)
does not provide.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from cyclerfinder.search.niching_ga import (
    DeterministicCrowdingConfig,
    _Encoding,
    run_deterministic_crowding,
)


def _two_peak_fitness(x: np.ndarray) -> float:
    """Two equal Gaussian peaks at 0.25 and 0.75 in x[0]; x[1] is fixed."""
    return float(np.exp(-(((x[0] - 0.25) / 0.05) ** 2)) + np.exp(-(((x[0] - 0.75) / 0.05) ** 2)))


def test_encoding_bounds_and_fixed_variables() -> None:
    enc = _Encoding([(0.0, 1.0), (0.3, 0.3), (-2.0, 2.0)], bits=8)
    assert enc.free_idx == [0, 2]
    assert enc.n_bits == 16
    zeros = np.zeros((1, 16), dtype=bool)
    ones = np.ones((1, 16), dtype=bool)
    np.testing.assert_allclose(enc.decode(zeros)[0], [0.0, 0.3, -2.0])
    np.testing.assert_allclose(enc.decode(ones)[0], [1.0, 0.3, 2.0])


def test_two_separated_optima_both_survive() -> None:
    """Both equal-height niches retain a substantial subpopulation."""
    config = DeterministicCrowdingConfig(
        population_size=60, generations=60, bits_per_variable=16, seed=7
    )
    result = run_deterministic_crowding(_two_peak_fitness, [(0.0, 1.0), (0.5, 0.5)], config)
    x = result.phenotypes[:, 0]
    near_a = int(np.sum(np.abs(x - 0.25) < 0.05))
    near_b = int(np.sum(np.abs(x - 0.75) < 0.05))
    assert near_a >= 10, f"peak at 0.25 lost: {near_a} members"
    assert near_b >= 10, f"peak at 0.75 lost: {near_b} members"
    # Fixed variable stays pinned.
    np.testing.assert_allclose(result.phenotypes[:, 1], 0.5)


def test_checkpoint_resume_is_deterministic(tmp_path: Path) -> None:
    """One 30-generation call == 15 + 15 via checkpoint resume (same RNG path)."""
    config = DeterministicCrowdingConfig(
        population_size=20, generations=30, bits_per_variable=12, seed=3
    )
    bounds = [(0.0, 1.0)]
    full = run_deterministic_crowding(_two_peak_fitness, bounds, config)

    ckpt = tmp_path / "dc.npz"
    part1 = run_deterministic_crowding(
        _two_peak_fitness, bounds, config, checkpoint_path=ckpt, max_generations_this_call=15
    )
    assert part1.generations_run == 15
    part2 = run_deterministic_crowding(_two_peak_fitness, bounds, config, checkpoint_path=ckpt)
    assert part2.generations_run == 30
    np.testing.assert_array_equal(part2.phenotypes, full.phenotypes)
    np.testing.assert_array_equal(part2.fitness, full.fitness)
    assert len(part2.history) == 30


def test_child_replaces_only_closer_parent() -> None:
    """DC replacement is distance-paired, not fitness-paired.

    With mutation and crossover disabled, children are clones of their
    parents; the distance pairing must then match each child to its own
    parent (d=0) and the population must be exactly preserved — no member
    is displaced by a fitter but distant individual.
    """
    config = DeterministicCrowdingConfig(
        population_size=20,
        generations=25,
        bits_per_variable=10,
        p_crossover=0.0,
        p_mutation=0.0,
        seed=11,
    )
    result = run_deterministic_crowding(_two_peak_fitness, [(0.0, 1.0)], config)
    # Reproduce generation-0 population from the same seed and encoding.
    rng = np.random.default_rng(11)
    enc = _Encoding([(0.0, 1.0)], bits=10)
    genomes0 = rng.random((20, enc.n_bits)) < 0.5
    x0 = np.sort(enc.decode(genomes0)[:, 0])
    np.testing.assert_allclose(np.sort(result.phenotypes[:, 0]), x0)
