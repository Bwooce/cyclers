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
    mmr_a1_from_t0,
    mmr_bounds,
    mmr_t0,
)

MU = 0.001

# #585's Fable-reviewed 2-rung ladder (OUTSTANDING.md #585). Tests below use
# S_RUNG1 as the "default" representative width except where explicitly
# comparing both rungs.
S_RUNG1 = 0.15
S_RUNG2 = 0.30


def _system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=MU, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)


def _config_for(a1: float, *, s: float = S_RUNG1) -> IsolatedAsymmetricFitnessConfig:
    _bounds, x0_guess, ydot0_guess, t0 = mmr_bounds(a1, mu=MU, s=s)
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
    bounds, *_ = mmr_bounds(a1, mu=MU, s=S_RUNG1)
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
    bounds, *_ = mmr_bounds(a1, mu=MU, s=S_RUNG1)
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


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_mmr_a1_from_t0_inverts_mmr_t0(p: int, q: int, a1: float) -> None:
    """#585's drift-detection check needs an exact round-trip inverse."""
    t0 = mmr_t0(a1)
    assert mmr_a1_from_t0(t0) == pytest.approx(a1, rel=1e-12)


def test_mmr_a1_from_t0_raises_on_nonpositive_or_nonfinite() -> None:
    with pytest.raises(ValueError):
        mmr_a1_from_t0(0.0)
    with pytest.raises(ValueError):
        mmr_a1_from_t0(-1.0)
    with pytest.raises(ValueError):
        mmr_a1_from_t0(float("nan"))


def test_mmr_bounds_raises_on_nonpositive_a1() -> None:
    with pytest.raises(ValueError):
        mmr_bounds(0.0, s=S_RUNG1)
    with pytest.raises(ValueError):
        mmr_t0(-1.0)


def test_isolated_3d_asymmetric_fitness_rejects_wrong_shape() -> None:
    config = _config_for(0.763143)
    with pytest.raises(ValueError):
        isolated_3d_asymmetric_fitness(np.array([1.0, 2.0, 3.0]), config=config)


def test_mmr_bounds_requires_s_explicitly() -> None:
    """#585: ``s`` is keyword-only with NO default -- every caller must opt in.

    A flat 0.05 box (#582's original) has no single ``s`` that reproduces it
    identically across all 5 MMRs (that anisotropy is exactly what #585
    fixes), so silently defaulting ``s`` would either quietly resurrect the
    old flat box for one MMR only or silently change behaviour for callers
    who forgot the new parameter exists. Forcing a ``TypeError`` makes that
    an explicit, reviewed choice at every call site.
    """
    with pytest.raises(TypeError):
        mmr_bounds(0.763143, mu=MU)  # type: ignore[call-arg]


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
@pytest.mark.parametrize("s", [S_RUNG1, S_RUNG2], ids=("s=0.15", "s=0.30"))
def test_mmr_bounds_s_scaling_formula(p: int, q: int, a1: float, s: float) -> None:
    """#585's exact resonance-scaled formula: xdot0_abs=zdot0_abs=s*v_circ,
    z0_abs=max(0.05, s*a1)."""
    bounds, _x0_guess, _ydot0_guess, _t0 = mmr_bounds(a1, mu=MU, s=s)
    v_circ = math.sqrt((1.0 - MU) / a1)
    expected_xdot0_abs = s * v_circ
    expected_z0_abs = max(0.05, s * a1)
    z0_lo, z0_hi = bounds[1]
    xdot0_lo, xdot0_hi = bounds[2]
    zdot0_lo, zdot0_hi = bounds[4]
    assert z0_lo == pytest.approx(-expected_z0_abs)
    assert z0_hi == pytest.approx(expected_z0_abs)
    assert xdot0_lo == pytest.approx(-expected_xdot0_abs)
    assert xdot0_hi == pytest.approx(expected_xdot0_abs)
    assert zdot0_lo == pytest.approx(-expected_xdot0_abs)
    assert zdot0_hi == pytest.approx(expected_xdot0_abs)


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_mmr_bounds_z0_never_shrinks_below_582_floor(p: int, q: int, a1: float) -> None:
    """The ``max(0.05, s*a1)`` floor keeps z0_abs >= #582's stamped 0.05 box."""
    for s in (0.01, S_RUNG1, S_RUNG2):
        bounds, *_ = mmr_bounds(a1, mu=MU, s=s)
        z0_lo, z0_hi = bounds[1]
        assert z0_hi >= 0.05
        assert z0_lo <= -0.05


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_mmr_bounds_rung2_wider_than_rung1(p: int, q: int, a1: float) -> None:
    """Rung 2 (s=0.30) must be at least as wide as rung 1 (s=0.15) on every
    symmetry-breaking component, and x0/ydot0/T (unchanged fractions) must be
    IDENTICAL across rungs -- only the symmetry-breaking bounds move."""
    bounds1, x0_guess1, ydot01, t01 = mmr_bounds(a1, mu=MU, s=S_RUNG1)
    bounds2, x0_guess2, ydot02, t02 = mmr_bounds(a1, mu=MU, s=S_RUNG2)
    assert x0_guess1 == x0_guess2
    assert ydot01 == ydot02
    assert t01 == t02
    for idx in (0, 3, 5):  # x0, ydot0, T: unchanged across rungs
        assert bounds1[idx] == bounds2[idx]
    for idx in (1, 2, 4):  # z0, xdot0, zdot0: rung 2 >= rung 1
        lo1, hi1 = bounds1[idx]
        lo2, hi2 = bounds2[idx]
        assert hi2 >= hi1
        assert lo2 <= lo1
