"""#436 Task 1: CR3BP-independent direct-e>0 seed grid."""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_direct_seeding import (
    DirectEr3bpSeed,
    classify_no_cr3bp_limit,
    converge_direct_seed,
    direct_e_seed_grid,
)


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


def test_broucke_classifies_cr3bp_continuous() -> None:
    # Broucke 1969 TR 32-1360 Table 12 Family 7P is CR3BP-continuous by
    # construction (it is parameterised by e down to e=0). The published anchor
    # IC (0.1520965, 3.1608994) is the e=0.0001 member; the e=0.0549 member of
    # the SAME family lives at (0.13892235, 3.24428595) (reached by continuing
    # the anchor along the family). ``converge_direct_seed`` is a pure
    # single-shot corrector (no homotopy), so the seed must be the genuine
    # e=0.0549 member rather than the far-off e=0.0001 anchor, which a single
    # Newton shot at e=0.0549 cannot reach. Either way this validates the
    # NEGATIVE (cr3bp_continuous) branch against a known family.
    sys = ER3BPSystem(mu=0.012155, e=0.0549, primary_name="E", secondary_name="M")
    seed = DirectEr3bpSeed(
        "broucke-check",
        sys,
        np.array([0.13892235, 0, 0, 0, 3.24428595, 0]),
        2 * np.pi,
        True,
        0.0549,
        "test",
    )
    orb = converge_direct_seed(seed)
    assert orb is not None
    res = classify_no_cr3bp_limit(orb, sys, n_steps=20)
    assert res["status"] == "cr3bp_continuous"
