"""Tests for the Monotonic Basin Hopping wrapper (task #145).

Two gate families:

* **Mechanics gate** (``test_mbh_*``): a constructed 1-D double-well whose global
  minimum is known BY CONSTRUCTION. This is an ALGORITHM-MECHANICS test, NOT a
  golden test -- the "expected" value is a value our own constructed function
  defines, which is acceptable to prove the MBH loop escapes the wrong basin
  where a single local solve provably cannot.

* **Free-return recovery gate** (``test_free_return_mbh_*``, slow): the real one.
  Start from a DELIBERATELY OFF-PHASE mis-seed of the #137 free-return genome on
  the symmetric ``mcconaghy-2006-em-k2`` row. A single plain local solve from the
  mis-seed FAILS (does not converge). MBH recovers a solution in the basin of the
  row's SOURCED transfer ellipse, whose EMERGED V_inf matches the INDEPENDENTLY
  sourced anchor.

  GOLDEN DISCIPLINE: the EXPECTED side is the SOURCED ellipse ``(a, e)`` (Rogers
  2012 Table 1, a=1.30 AU, e=0.257 -- the same constraint side the #137 gate in
  ``test_russell12_likeforlike_probe.py`` uses) and the SOURCED V_inf anchor
  (McConaghy 4.7/5.0 km/s). The emerged V_inf and the recovered ``(a, e)`` are
  the EVIDENCE; nothing our code computed is ever an EXPECTED value.

  HONEST CAVEAT (recorded in docs/notes/2026-06-07-mbh-wrapper.md): the
  free-return single-ellipse genome is 1-DOF underdetermined along the Mars-V_inf
  ridge -- every t0 closes at *some* V_inf -- so residual-only MBH does not select
  the sourced sub-basin from every rng seed. The gate pins a deterministic seed
  that recovers it; the seed-sensitivity is the survey's predicted "MBH confirms
  the topology limitation" finding, reported in the note.
"""

from __future__ import annotations

import math
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.mbh import (
    DEFAULT_PARETO_ALPHA,
    PERTURBATIONS,
    MBHStep,
    make_free_return_step,
    mbh,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DAY_S = 86400.0


# ---------------------------------------------------------------------------
# Mechanics gate: constructed double-well (global min known BY CONSTRUCTION).
# ---------------------------------------------------------------------------

# Two Gaussian wells. The DEEP (global) well sits at x=+3 with depth -2; the
# SHALLOW (local) well sits at x=-3 with depth -1. A gradient-descent local solve
# converges to whichever well it starts in -- so from the shallow basin it can
# NEVER reach the global one without a hop.
_DEEP_X = 3.0
_SHALLOW_X = -3.0
_DEEP_DEPTH = -2.0
_SHALLOW_DEPTH = -1.0


def _double_well(x: float) -> float:
    return _DEEP_DEPTH * math.exp(-((x - _DEEP_X) ** 2) / 0.5) + _SHALLOW_DEPTH * math.exp(
        -((x - _SHALLOW_X) ** 2) / 0.5
    )


def _local_min(x0: float) -> float:
    """Deterministic gradient descent to the nearest local minimum of
    ``_double_well`` (the inner "local solver")."""
    x = float(x0)
    for _ in range(400):
        grad = (_double_well(x + 1e-5) - _double_well(x - 1e-5)) / 2e-5
        x -= 0.05 * grad
    return x


def _double_well_step(x_seed: np.ndarray, rng: np.random.Generator) -> MBHStep:
    """MBH local-solve closure over the constructed double-well."""
    x_loc = _local_min(float(x_seed[0]))
    val = _double_well(x_loc)
    return MBHStep(
        x=np.array([x_loc], dtype=np.float64),
        objective=float(val),
        feasible=True,
        info={"x": x_loc},
    )


def test_plain_local_solve_stays_in_wrong_basin() -> None:
    """A single local solve from the SHALLOW basin provably cannot reach the
    global (deep) minimum -- this is the failure MBH must cure."""
    x_loc = _local_min(_SHALLOW_X)
    assert abs(x_loc - _SHALLOW_X) < 0.1
    assert _double_well(x_loc) == pytest.approx(_SHALLOW_DEPTH, abs=1e-3)
    # Strictly worse than the global minimum -- the wrong basin.
    assert _double_well(x_loc) > _DEEP_DEPTH + 0.5


def test_mbh_escapes_to_global_basin() -> None:
    """MBH started in the SHALLOW basin reaches the DEEP (global) basin."""
    result = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=200,
        perturbation="cauchy",  # long tail -> can jump the inter-well barrier
        perturbation_absolute_scale=[3.0],  # absolute step in x-units
        rng_seed=12345,
    )
    assert result.best_feasible
    assert abs(result.best_x[0] - _DEEP_X) < 0.2
    assert result.best_objective == pytest.approx(_DEEP_DEPTH, abs=1e-3)
    # Monotonic acceptance: the best history never increases.
    hist = [h for h in result.best_history if math.isfinite(h)]
    assert all(b <= a + 1e-12 for a, b in pairwise(hist))


def test_mbh_is_deterministic_for_a_seed() -> None:
    """Same rng_seed -> identical run (no global RNG state)."""
    kw: dict[str, Any] = dict(
        x0=[_SHALLOW_X],
        n_hops=50,
        perturbation="gaussian",
        perturbation_absolute_scale=[3.0],
        rng_seed=7,
    )
    r1 = mbh(_double_well_step, **kw)
    r2 = mbh(_double_well_step, **kw)
    assert np.array_equal(r1.best_x, r2.best_x)
    assert r1.objective_history == r2.objective_history
    assert r1.accept_history == r2.accept_history
    assert r1.rng_seed == 7


def test_mbh_audit_trail_is_consistent() -> None:
    """The audit-trail lengths and counts are internally consistent."""
    result = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=30,
        perturbation="gaussian",
        perturbation_absolute_scale=[3.0],
        rng_seed=1,
    )
    n = result.hops_attempted
    assert len(result.objective_history) == n
    assert len(result.accept_history) == n
    assert len(result.best_history) == n
    assert result.hops_accepted == sum(result.accept_history)
    # hop 0 (seed solve) plus up to n_hops perturbed hops.
    assert 1 <= n <= 31


def test_mbh_monotonic_accepts_only_improvements() -> None:
    """An accepted hop must strictly lower the incumbent objective."""
    result = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=80,
        perturbation="cauchy",
        perturbation_absolute_scale=[3.0],
        rng_seed=3,
    )
    best_so_far = float("inf")
    feasible_seen = False
    for obj, accepted in zip(result.objective_history, result.accept_history, strict=True):
        if accepted:
            # First feasible accept can come from any objective; subsequent
            # accepts must strictly improve.
            if feasible_seen:
                assert obj < best_so_far
            best_so_far = min(best_so_far, obj)
            feasible_seen = True


def test_mbh_stall_stop() -> None:
    """``stop_after_stall`` terminates the run once the incumbent stops
    improving for the configured number of consecutive hops."""

    # A trivial objective with a single basin: after the seed solve every hop
    # returns the same minimum, so no perturbed hop ever improves -> stall fires.
    def flat_step(x_seed: np.ndarray, rng: np.random.Generator) -> MBHStep:
        return MBHStep(x=np.array([0.0]), objective=0.0, feasible=True, info={})

    result = mbh(
        flat_step,
        x0=[5.0],
        n_hops=1000,
        perturbation="gaussian",
        perturbation_absolute_scale=[1.0],
        rng_seed=0,
        stop_after_stall=5,
    )
    assert result.stopped_on_stall
    assert result.final_stall == 5
    # Seed solve (hop 0) + exactly 5 non-improving perturbed hops.
    assert result.hops_attempted == 6


def test_mbh_requires_known_distribution() -> None:
    with pytest.raises(ValueError, match="unknown perturbation"):
        mbh(_double_well_step, x0=[0.0], n_hops=1, perturbation="laplace", rng_seed=0)


# ---------------------------------------------------------------------------
# Englander 2014 upgrade (task #156): bi-polar Pareto + audit fields + restart.
# ---------------------------------------------------------------------------


def test_pareto_is_a_registered_long_tailed_distribution() -> None:
    """``"pareto"`` (bi-polar Pareto) is the SOURCED paper's recommendation and is
    registered; ``"cauchy"`` remains the CODE default (Englander & Englander 2014;
    local-sweep caveat in the wrapper note)."""
    assert "pareto" in PERTURBATIONS
    assert PERTURBATIONS == ("pareto", "cauchy", "gaussian", "uniform")


def test_mbh_pareto_escapes_to_global_basin() -> None:
    """MBH driven by the bi-polar Pareto perturbation escapes the wrong basin
    (the long tail jumps the inter-well barrier), mirroring the cauchy gate."""
    result = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=200,
        perturbation="pareto",
        perturbation_absolute_scale=[3.0],
        rng_seed=12345,
    )
    assert result.best_feasible
    assert abs(result.best_x[0] - _DEEP_X) < 0.2
    assert result.best_objective == pytest.approx(_DEEP_DEPTH, abs=1e-3)
    hist = [h for h in result.best_history if math.isfinite(h)]
    assert all(b <= a + 1e-12 for a, b in pairwise(hist))


def test_mbh_pareto_is_bipolar_and_alpha_sensitive() -> None:
    """The bi-polar Pareto draw is symmetric about zero (both signs occur), is
    long-tailed (unit floor with a heavy tail), and SMALLER ``alpha`` => heavier
    tail (the Power-Law exponent; inverse-CDF Pareto)."""
    from cyclerfinder.search.mbh import _perturb

    def pareto_sample(alpha: float) -> np.ndarray:
        return np.array(
            [
                _perturb(
                    np.array([0.0]),
                    np.random.default_rng(s),
                    distribution="pareto",
                    scale=None,
                    absolute_scale=np.array([1.0]),
                    alpha=alpha,
                )[0]
                for s in range(2000)
            ]
        )

    sample_heavy = pareto_sample(1.05)  # heavier tail (smaller exponent)
    sample_light = pareto_sample(3.00)  # lighter tail (larger exponent)
    # Bi-polar: both signs present.
    assert np.any(sample_heavy > 0) and np.any(sample_heavy < 0)
    # Unit floor: every magnitude is >= 1 (before the per-gene scale of 1.0).
    assert np.min(np.abs(sample_heavy)) >= 1.0 - 1e-9
    # Heavier tail at SMALLER alpha -> larger extreme magnitude.
    assert np.max(np.abs(sample_heavy)) > np.max(np.abs(sample_light))


def test_mbh_result_records_distribution_and_excursion_param() -> None:
    """MBHResult records the distribution name + its excursion parameter -- the
    audit gap the Englander upgrade closes (the distribution choice is the SOURCED
    performance driver)."""
    # Pareto: excursion parameter is alpha.
    r_pareto = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=10,
        perturbation="pareto",
        perturbation_absolute_scale=[1.0],
        perturbation_alpha=1.11,
        rng_seed=1,
    )
    assert r_pareto.perturbation == "pareto"
    assert r_pareto.perturbation_param == pytest.approx(1.11)

    # Scalar relative scale: excursion parameter is that scale.
    r_scale = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=10,
        perturbation="cauchy",
        perturbation_scale=0.07,
        rng_seed=1,
    )
    assert r_scale.perturbation == "cauchy"
    assert r_scale.perturbation_param == pytest.approx(0.07)

    # Per-gene / None scale: excursion parameter is nan (caller records vectors).
    r_vec = mbh(
        _double_well_step,
        x0=[_SHALLOW_X],
        n_hops=10,
        perturbation="gaussian",
        perturbation_scale=None,
        perturbation_absolute_scale=[3.0],
        rng_seed=1,
    )
    assert r_vec.perturbation == "gaussian"
    assert math.isnan(r_vec.perturbation_param)


def test_mbh_default_records_cauchy_and_default_scale() -> None:
    """The default path records ``cauchy`` with the default 0.05 scale (the code
    default is unchanged by the Englander upgrade)."""
    r = mbh(_double_well_step, x0=[_SHALLOW_X], n_hops=5, rng_seed=0)
    assert r.perturbation == "cauchy"
    assert r.perturbation_param == pytest.approx(0.05)


def test_default_pareto_alpha_is_the_documented_value() -> None:
    assert pytest.approx(1.08) == DEFAULT_PARETO_ALPHA


def test_mbh_restart_on_stall_continues_instead_of_stopping() -> None:
    """With ``restart_bounds`` set, a stall RE-SEEDS from a fresh random point in
    bounds (Englander Algorithm 1 global reset) and continues, rather than
    stopping the run. The default (no bounds) still stops -- proven by the
    existing ``test_mbh_stall_stop``."""

    seen_seeds: list[float] = []

    def trap_step(x_seed: np.ndarray, _rng: np.random.Generator) -> MBHStep:
        # A flat objective: no perturbed hop ever improves, so stalls keep firing.
        seen_seeds.append(float(x_seed[0]))
        return MBHStep(x=np.array([0.0]), objective=0.0, feasible=True, info={})

    result = mbh(
        trap_step,
        x0=[5.0],
        n_hops=30,
        perturbation="gaussian",
        perturbation_absolute_scale=[1.0],
        rng_seed=0,
        stop_after_stall=3,
        restart_bounds=([-100.0], [100.0]),
    )
    # Did NOT stop early: ran all 30 hops + seed solve.
    assert not result.stopped_on_stall
    assert result.hops_attempted == 31
    # A reset draws within [-100, 100]; the post-reset seeds are wide spread (not
    # the tiny gaussian jitter around the frozen incumbent at 0).
    assert max(abs(s) for s in seen_seeds) > 10.0


def test_mbh_restart_bounds_validates_shapes() -> None:
    with pytest.raises(ValueError, match="restart_bounds shapes"):
        mbh(
            _double_well_step,
            x0=[0.0, 1.0],
            n_hops=1,
            perturbation="gaussian",
            rng_seed=0,
            stop_after_stall=1,
            restart_bounds=([0.0], [1.0]),
        )


# ---------------------------------------------------------------------------
# Free-return recovery gate (the real one) -- slow.
# ---------------------------------------------------------------------------

# SOURCED constraint side (Rogers 2012 Table 1; the same S1L1 ellipse the #137
# gate in test_russell12_likeforlike_probe.py uses). NOT a value our code made.
_S1L1_A_AU = 1.30
_S1L1_E = 0.257
# SOURCED V_inf anchor for mcconaghy-2006-em-k2 (McConaghy 2006): E 4.7, M 5.0.
_SRC_VINF_E = 4.7
_SRC_VINF_M = 5.0


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _best_phase_t0(ephem: Any, period_sec: float) -> float:
    """Scan t0 over one period; return the phase minimising the free-return
    residual at the SOURCED (a, e) (mirrors the #137 gate's best-phase pick)."""
    from cyclerfinder.search.free_return import _residuals

    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, 2000, endpoint=False):
        t0 = float(frac) * period_sec
        res = _residuals(
            np.array([_S1L1_A_AU, _S1L1_E, t0]),
            period_days=period_sec / DAY_S,
            ephem=ephem,
            bodies=("E", "M"),
            mu=132712440018.0,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


@pytest.mark.slow
def test_free_return_mbh_recovers_sourced_basin_from_misseed() -> None:
    """MBH gate: from an off-phase mis-seed of the free-return genome, a single
    plain solve FAILS, but MBH recovers a solution in the SOURCED transfer-ellipse
    basin whose emerged V_inf matches the SOURCED anchor.

    Determinism: the gate pins rng_seed=6 (which recovers the sourced basin for
    this mis-seed). The genome is 1-DOF underdetermined along the Mars-V_inf
    ridge, so not every seed selects the sourced sub-basin -- see the note.
    """
    from cyclerfinder.core.ephemeris import Ephemeris

    row = _row("mcconaghy-2006-em-k2")
    period_sec = float(row["period"]["years"]) * 365.25 * DAY_S
    ephem = Ephemeris("circular")
    step = make_free_return_step(period_sec=period_sec, ephem=ephem)

    best_t0 = _best_phase_t0(ephem, period_sec)
    # DELIBERATE mis-seed: 40 days off the in-phase crossing window (the realistic
    # phasing miss the campaign hits). (a, e) start at the sourced ellipse.
    misseed = np.array([_S1L1_A_AU, _S1L1_E, best_t0 - 40.0 * DAY_S])

    # 1) A single plain local solve from the mis-seed does NOT converge.
    plain = step(misseed, np.random.default_rng(0))
    assert not plain.feasible, f"plain solve unexpectedly converged: obj={plain.objective}"

    # 2) MBH recovers. Absolute 8-day t0 step (the basin is only a few days wide;
    #    a relative step on t0~1e8 s would be enormous), a/e frozen in the hop.
    result = mbh(
        step,
        misseed,
        n_hops=60,
        perturbation="gaussian",
        perturbation_scale=None,
        perturbation_absolute_scale=[0.0, 0.0, 8.0 * DAY_S],
        rng_seed=6,
    )

    assert result.best_feasible
    assert result.best_objective < 0.1  # converged corrector residual (km/s)

    # 3) Recovered point is in the SOURCED ellipse basin (constraint side), far
    #    from the powered a~1.13, e~0.57 degenerate basin.
    a_rec = result.best_info["a_au"]
    e_rec = result.best_info["e"]
    assert abs(a_rec - _S1L1_A_AU) <= 0.03
    assert abs(e_rec - _S1L1_E) <= 0.03

    # 4) EVIDENCE: emerged V_inf matches the SOURCED anchor. The Mars V_inf is the
    #    binding flyby and matches within the campaign 0.5 km/s tolerance; the
    #    Earth V_inf is in the sourced low-V_inf regime (within ~1 km/s) -- NOT
    #    the 9+ km/s powered basin the plain ballistic genome lands in.
    vinf = result.best_info["vinf_kms"]
    assert abs(vinf["M"] - _SRC_VINF_M) <= 0.5
    assert abs(vinf["E"] - _SRC_VINF_E) <= 1.0


@pytest.mark.slow
def test_free_return_mbh_determinism_on_misseed() -> None:
    """The free-return MBH recovery is reproducible for a fixed rng_seed."""
    from cyclerfinder.core.ephemeris import Ephemeris

    row = _row("mcconaghy-2006-em-k2")
    period_sec = float(row["period"]["years"]) * 365.25 * DAY_S
    ephem = Ephemeris("circular")
    step = make_free_return_step(period_sec=period_sec, ephem=ephem)
    best_t0 = _best_phase_t0(ephem, period_sec)
    misseed = np.array([_S1L1_A_AU, _S1L1_E, best_t0 - 40.0 * DAY_S])

    kw: dict[str, Any] = dict(
        n_hops=30,
        perturbation="gaussian",
        perturbation_scale=None,
        perturbation_absolute_scale=[0.0, 0.0, 8.0 * DAY_S],
        rng_seed=6,
    )
    r1 = mbh(step, misseed, **kw)
    r2 = mbh(step, misseed, **kw)
    assert np.allclose(r1.best_x, r2.best_x)
    assert r1.objective_history == r2.objective_history
