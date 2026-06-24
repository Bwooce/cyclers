from __future__ import annotations

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_discovery import (
    Er3bpSeed,
    continue_and_monitor,
    standard_family_seeds,
)


def test_standard_family_seeds_returns_usable_floor() -> None:
    seeds = standard_family_seeds(target_e=0.0549)
    assert len(seeds) >= 1
    for s in seeds:
        assert isinstance(s, Er3bpSeed)
        assert s.state0.shape == (6,)
        assert s.period_f > 0.0
        assert 0.0 < s.target_e < 1.0
        assert s.system.primary_name and s.system.secondary_name
        assert s.source


def test_continue_and_monitor_classifies_survival() -> None:
    seed = standard_family_seeds(target_e=0.01)[0]
    trace = continue_and_monitor(seed, n_steps=5)
    assert trace.outcome in {"survives", "dies", "bifurcates"}
    assert len(trace.steps) >= 1
    s0 = trace.steps[0]
    assert 0.0 <= s0.e <= 0.01 + 1e-9
    assert s0.stability_tag in {"stable", "unstable", "marginal", "unknown"}


def test_continue_and_monitor_dies_on_infeasible_seed() -> None:
    seed = Er3bpSeed(
        label="garbage",
        system=ER3BPSystem(mu=0.0121550, e=0.5, primary_name="E", secondary_name="M"),
        state0=np.array([5.0, 5.0, 0.0, 9.0, 9.0, 0.0]),
        period_f=2.0 * np.pi,
        is_half_period_residual=True,
        target_e=0.5,
        source="test-garbage",
    )
    trace = continue_and_monitor(seed, n_steps=3)
    assert trace.outcome == "dies"
    assert trace.e_max_reached < 0.5
