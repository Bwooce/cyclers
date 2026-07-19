"""Task #656: Pluto-Charon higher-(k1,k2) topology extension — test suite.

`#504`/`#549` only ever swept `(k1,k2)` in {(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)}
at Pluto-Charon's own mass ratio (Ross-RT 2026's own Table I stops at (3,3)
too). `scripts/run_656_pc_higher_kk_sweep.py` extends the topology list to
the 9 higher topologies (4,1)-(5,5), reusing #504/#549's `sweep_family_grid`
machinery verbatim (no changes to `real_binary_kk_sweep.py` itself).

This module is a REGRESSION suite for the driver script, not a re-run of
the full 9-topology discovery grid search (each topology's ~256-point grid
takes ~2 minutes single-core — see docs/notes/scratch/656_pc_kk45_sweep_raw.txt
for the actual discovery-run results, summarized in `data/OUTSTANDING.md`'s
`#656` entry). It covers:

1. The topology-enumeration convention (k2<=k1, matching #549's own
   convention exactly — no invented pairs).
2. The per-topology grid-parameter scaling function.
3. A FAST positive-control regression of `_grid_seed_search` — the exact
   seed-finding step behind every one of the 9 new topologies' "no orbit
   found in grid" negatives (not just the anchor-seeded
   `sweep_32_positive_control`) — a tight grid seeded around the already
   -admitted (3,2) solution must still recover it via the grid-search seed
   step, proving that step (not just the anchor path) is trustworthy before
   trusting any of the 9 new grid-only topology negatives. Also documents a
   pre-existing, out-of-scope quirk found while writing this test: the
   downstream C-sweep stage's `hc=None` auto-redetection can walk a
   correctly-seeded orbit onto an unrelated branch (see the test's own
   docstring) — a property of `sweep_family_grid`/`c_sweep_find_nu_zero`
   shared with every topology #504/#549/#656 have ever swept, not something
   this task introduces or is in scope to fix.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from run_656_pc_higher_kk_sweep import (  # noqa: E402
    TOPOLOGIES,
    _grid_for_topology,
)

from cyclerfinder.search.binary_star_search import winding_topology  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    PC_MU,
    _grid_seed_search,
    make_pluto_charon_system,
)

# ---------------------------------------------------------------------------
# 1. Topology-enumeration convention
# ---------------------------------------------------------------------------


def test_656_topology_list_is_exactly_4x_and_5x_with_k2_le_k1() -> None:
    """(4,1)-(5,5) means the 9 pairs with k2<=k1, matching #549's own convention.

    Every anchor/grid target #504/#549 ever swept has k2<=k1 (e.g. (3,1),
    (3,2),(3,3) but never a k2>k1 pair like (1,3)) -- this extension must
    match that, not invent a new convention.
    """
    expected = {
        (4, 1),
        (4, 2),
        (4, 3),
        (4, 4),
        (5, 1),
        (5, 2),
        (5, 3),
        (5, 4),
        (5, 5),
    }
    assert set(TOPOLOGIES) == expected
    assert len(TOPOLOGIES) == 9
    for k1, k2 in TOPOLOGIES:
        assert k2 <= k1, f"({k1},{k2}) violates the k2<=k1 convention"
        assert k1 in (4, 5)


# ---------------------------------------------------------------------------
# 2. Per-topology grid-parameter scaling
# ---------------------------------------------------------------------------


def test_656_grid_for_topology_scales_with_k1_plus_k2() -> None:
    """hc_list and period_guess must grow with (k1+k2), not stay fixed.

    A one-size-fits-all grid (as #504/#549 used for (2,1)/(2,2)) would
    systematically miss higher-winding orbits that need a longer period
    guess and a higher half-crossings count.
    """
    x0_41, c_41, hc_41, period_41 = _grid_for_topology(4, 1)
    x0_55, c_55, hc_55, period_55 = _grid_for_topology(5, 5)

    assert period_55 > period_41, "period_guess must grow with k1+k2"
    assert max(hc_55) > max(hc_41), "hc_list must shift upward with k1+k2"
    assert min(hc_41) >= 1, "hc must stay positive"

    # x0/C ranges are shared across topologies (only hc/period scale).
    assert np.array_equal(x0_41, x0_55)
    assert np.array_equal(c_41, c_55)

    # C range must stay strictly below C_L1(PC) (a physically valid CR3BP
    # C-sweep upper bound), matching #504's own (2,1)/(2,2) convention.
    from cyclerfinder.search.pluto_charon_kk_sweep import _c_l1

    c_l1 = _c_l1(PC_MU)
    assert c_41.max() < c_l1


def test_656_all_nine_topologies_have_sane_grids() -> None:
    for k1, k2 in TOPOLOGIES:
        x0_grid, c_grid, hc_list, period_guess = _grid_for_topology(k1, k2)
        assert len(x0_grid) == 8
        assert len(c_grid) == 8
        assert len(hc_list) >= 1
        assert all(h >= 1 for h in hc_list)
        assert period_guess > 0


# ---------------------------------------------------------------------------
# 3. Fast positive-control regression of the grid-SEED-finding step itself
# ---------------------------------------------------------------------------


def test_656_grid_seed_search_recovers_admitted_pc_32_seed() -> None:
    """`_grid_seed_search` (the exact seed-finding step behind every one of
    the 9 new topologies' "no orbit found in grid" negatives) must recover
    the already-admitted PC (3,2) member from a tight grid seeded around its
    own known solution.

    #504/#549's grid-search path was previously only exercised on (2,1)/
    (2,2), which are BOTH clean negatives -- never proven capable of
    recovering a real positive from the grid path alone. This closes that
    gap cheaply (a 2x2x1 grid around the known solution, not the full
    8x8x4 discovery grid) before trusting any (4,1)-(5,5) "no orbit found in
    grid" negative as genuine rather than a grid-machinery blind spot.

    Deliberately tests `_grid_seed_search` directly rather than the full
    `sweep_family_grid` pipeline: while writing this test, the downstream
    C-sweep stage (`c_sweep_find_nu_zero`, called with `hc=None` --
    auto-redetecting the crossing count rather than holding the seed's own
    `hc` fixed) was found to walk this exact (3,2) seed onto an UNRELATED
    stable branch ((k1,k2)=(4,0)-family region, topology_ok=False) even
    though the seed step itself was correct. That is a pre-existing,
    already-documented property of this shared machinery (see
    `mu_step_to_system_tracking_c_l1`'s own docstring: "an auto-redetected
    crossing index can snap onto a different, unrelated branch") -- NOT
    something this task's scope (a thin driver reusing `sweep_family_grid`
    verbatim, per its own dispatch) is meant to fix. Isolating the seed step
    keeps this regression test honest about what it actually verifies.
    """
    pc = make_pluto_charon_system()
    # Tight grid bracketing the catalogued ross-rt-pc-cycler-32-2026 solution
    # (x0=-0.693198287043394, C=3.57951501972907, T=11.8334625170346 TU,
    # hc=6). The Newton corrector is local: seeded even ~0.01 nd away in x0
    # or ~0.5 TU away in the period guess, it converges to an ENTIRELY
    # different branch (empirically checked: a (6,0) retrograde-Pluto-only
    # family, not (3,2)) -- so this grid/period_guess must stay close to the
    # known solution to test the intended thing (does the grid-search CODE
    # correctly recover a known answer), not "does the corrector find (3,2)
    # from a blind seed" (that's what the real (4,1)-(5,5) discovery grid is
    # for, and it need not be tight).
    x0_grid = np.array([-0.6935, -0.6930])
    c_grid = np.array([3.5793, 3.5798])
    hc_list = (6,)
    seed = _grid_seed_search(pc, 3, 2, x0_grid, c_grid, hc_list, 11.8335, per_call_timeout=5)
    assert seed is not None, "grid seed search failed to recover the known PC (3,2) seed"
    assert seed.converged
    assert abs(seed.x0 - (-0.693198287043394)) < 5e-3
    assert abs(seed.jacobi - 3.57951501972907) < 1e-3
    assert abs(seed.period - 11.8334625170346) < 0.5

    state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0])
    topo = winding_topology(pc.mu, state0, seed.period)
    assert topo.k1 == 3 and topo.k2 == 2, f"seed topology ({topo.k1},{topo.k2}) != (3,2)"
    assert topo.prograde
    assert topo.reaches_secondary
