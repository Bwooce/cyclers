"""Tests for #582's isolated_3d_asymmetric_fitness (GA fitness + bounds).

Property/regression tests for the fitness FORMULA and the bounds-box
construction, not sourced-discovery goldens (that discipline applies to
CLAIMED trajectory results, not to unit-testing a scoring function's own
mathematical behaviour). Where the tests do reference #440's own converged
resonant members (via ``resonant_po_seed``/``all_mmr_seeds``), those in turn
trace to the GO-gate-validated Antoniadou & Libert 2018 a1 table -- see
``tests/search/test_er3bp_isolated_seeds.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.er3bp_isolated_seeds import MMR_SEMI_MAJOR_AXES, all_mmr_seeds
from cyclerfinder.search.isolated_3d_asymmetric_fitness import (
    DEFAULT_PRIMARY_EXCLUSION_RADIUS,
    IsolatedAsymmetricFitnessConfig,
    genome_to_state0,
    isolated_3d_asymmetric_fitness,
    mmr_bounds,
    mmr_t0,
)

MU = 0.001


def _system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=MU, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)


def _config_for(a1: float) -> IsolatedAsymmetricFitnessConfig:
    _bounds, x0_guess, ydot0_guess, t0 = mmr_bounds(a1, mu=MU)
    seed_state = np.array([x0_guess, 0.0, 0.0, 0.0, ydot0_guess, 0.0])
    jacobi_target = float(cr3bp.jacobi_constant(seed_state, MU))
    return IsolatedAsymmetricFitnessConfig(mu=MU, t0=t0, jacobi_target=jacobi_target)


def test_mmr_t0_matches_440_docstring_3_2() -> None:
    """#440's own docstring: T0 = 2*pi/(a1**-1.5 - 1), '4*pi for 3/2'."""
    t0 = mmr_t0(0.763143)
    assert t0 == pytest.approx(4.0 * math.pi, rel=1e-3)


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_genome_to_state0_pins_y0_zero(p: int, q: int, a1: float) -> None:
    genome = np.array([a1, 0.01, 0.0, 0.4, 0.0, mmr_t0(a1)])
    state0 = genome_to_state0(genome)
    assert state0[1] == 0.0
    assert state0.shape == (6,)


def test_death_penalty_short_period() -> None:
    """T <= T0/2 is a death penalty (excludes degenerate near-equilibrium loops)."""
    config = _config_for(0.763143)
    genome = np.array([0.75, 0.0, 0.0, 0.4, 0.0, config.t0 * 0.4])
    assert isolated_3d_asymmetric_fitness(genome, config=config) == 0.0
    # just above the floor is NOT death-penalized (only <= T0/2 is)
    genome_ok = np.array([0.75, 0.0, 0.0, 0.4, 0.0, config.t0 * 0.9])
    assert isolated_3d_asymmetric_fitness(genome_ok, config=config) > 0.0


def test_death_penalty_primary_collision() -> None:
    """A trajectory that grazes the primary mid-arc is death-penalized."""
    config = _config_for(0.763143)
    # start just inside the primary exclusion radius (r ~ 0 immediately)
    genome = np.array([-MU + DEFAULT_PRIMARY_EXCLUSION_RADIUS * 0.5, 0.0, 0.0, 0.4, 0.0, config.t0])
    assert isolated_3d_asymmetric_fitness(genome, config=config) == 0.0


def test_death_penalty_secondary_collision() -> None:
    config = _config_for(0.763143)
    genome = np.array([1.0 - MU, 0.0, 0.0, 0.0, 0.0, config.t0])
    assert isolated_3d_asymmetric_fitness(genome, config=config) == 0.0


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_fitness_bounded_zero_to_one(p: int, q: int, a1: float) -> None:
    config = _config_for(a1)
    bounds, *_ = mmr_bounds(a1, mu=MU)
    rng = np.random.default_rng(42)
    for _ in range(25):
        genome = np.array([rng.uniform(lo, hi) for lo, hi in bounds])
        f = isolated_3d_asymmetric_fitness(genome, config=config)
        assert 0.0 <= f <= 1.0, f"fitness {f} out of [0,1] for genome {genome}"


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_mmr_bounds_contains_known_440_seed(p: int, q: int, a1: float) -> None:
    """The GA bounds box actually contains the known circular member (#440).

    Caught a real bug while building #582: the initial ydot0_frac=0.15 box
    did NOT contain the 3:2 member's converged ydot0 (it lands ~20% off the
    analytic guess, #440's own "most eccentric member" case) -- silently
    excluding the answer the positive control needed to find. This is a
    regression guard against that exact failure mode recurring.
    """
    system = _system()
    seed = all_mmr_seeds(mu=MU)[MMR_SEMI_MAJOR_AXES.index((p, q, a1))]
    assert seed.converged
    bounds, *_ = mmr_bounds(a1, mu=MU)
    x0_lo, x0_hi = bounds[0]
    ydot0_lo, ydot0_hi = bounds[3]
    t_lo, t_hi = bounds[5]
    assert x0_lo <= seed.state0[0] <= x0_hi, f"{p}:{q} x0 {seed.state0[0]} outside {bounds[0]}"
    assert ydot0_lo <= seed.state0[4] <= ydot0_hi, (
        f"{p}:{q} ydot0 {seed.state0[4]} outside {bounds[3]}"
    )
    assert t_lo <= seed.period <= t_hi, f"{p}:{q} T {seed.period} outside {bounds[5]}"
    del system  # constructed for parity with other tests; unused here


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_fitness_near_maximal_at_known_440_seed(p: int, q: int, a1: float) -> None:
    """The fitness function scores the known #440 circular member highly.

    Necessary property for the GA to be able to find it: a near-zero
    periodicity defect at the seed's own (near-target) Jacobi constant must
    score near the fitness ceiling of 1.0.
    """
    system = _system()
    seed = all_mmr_seeds(mu=MU)[MMR_SEMI_MAJOR_AXES.index((p, q, a1))]
    config = _config_for(a1)
    genome = np.array(
        [
            seed.state0[0],
            seed.state0[2],
            seed.state0[3],
            seed.state0[4],
            seed.state0[5],
            seed.period,
        ]
    )
    f = isolated_3d_asymmetric_fitness(genome, config=config)
    assert f > 0.95, f"{p}:{q} fitness at known seed only {f}"
    del system


def test_fitness_lower_for_far_off_periodicity() -> None:
    """A grossly non-periodic state scores much lower than the known seed."""
    a1 = 0.763143
    config = _config_for(a1)
    system = _system()
    seed = all_mmr_seeds(mu=MU)[0]
    assert seed.p == 3 and seed.q == 2
    good_genome = np.array(
        [
            seed.state0[0],
            seed.state0[2],
            seed.state0[3],
            seed.state0[4],
            seed.state0[5],
            seed.period,
        ]
    )
    bad_genome = good_genome.copy()
    bad_genome[3] = 1.2  # xdot0: large, breaks periodicity badly
    f_good = isolated_3d_asymmetric_fitness(good_genome, config=config)
    f_bad = isolated_3d_asymmetric_fitness(bad_genome, config=config)
    assert f_good > f_bad
    del system


def test_mmr_bounds_raises_on_nonpositive_a1() -> None:
    with pytest.raises(ValueError):
        mmr_bounds(0.0)
    with pytest.raises(ValueError):
        mmr_t0(-1.0)


def test_isolated_3d_asymmetric_fitness_rejects_wrong_shape() -> None:
    config = _config_for(0.763143)
    with pytest.raises(ValueError):
        isolated_3d_asymmetric_fitness(np.array([1.0, 2.0, 3.0]), config=config)
