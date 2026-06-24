"""#436 Task 1: CR3BP-independent direct-e>0 seed grid."""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_direct_seeding import DirectEr3bpSeed, direct_e_seed_grid


def test_direct_e_seed_grid_corners() -> None:
    sys = ER3BPSystem(mu=0.012155, e=0.0549, primary_name="E", secondary_name="M")
    seeds = direct_e_seed_grid(sys, (0.1, 0.9), (-3.5, 3.5), n_x=3, n_ydot=3, period_f=2 * np.pi)
    assert len(seeds) == 9
    assert all(isinstance(s, DirectEr3bpSeed) for s in seeds)
    assert all(s.state0.shape == (6,) for s in seeds)
    assert all(s.state0[1] == 0 and s.state0[3] == 0 and s.state0[5] == 0 for s in seeds)
    assert all(s.target_e == 0.0549 for s in seeds)
    assert seeds[0].state0[0] == 0.1 and seeds[-1].state0[0] == 0.9
    assert "grid" in seeds[0].source.lower() and "broucke" not in seeds[0].source.lower()
