"""#436 Task 3: ER3BP fixed-period saddle-center branch-switcher.

Safety / golden test. #432/#435 found ZERO bifurcations in the e-continuation of
the Broucke 7P family (it stays hyperbolic — no elliptic<->hyperbolic Floquet
transition), so there is no genuine bifurcating parent to branch off. This test
therefore exercises the SAFE behaviour on a non-bifurcating parent: the
switcher must not raise, must return the documented type, and must NOT fabricate
a spurious distinct family on an orbit that does not actually bifurcate.
"""

from __future__ import annotations

import math

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_branching import branch_at_saddle_center_er3bp
from cyclerfinder.genome.er3bp_periodic import ER3BPPeriodicOrbit

# Broucke 1969 TR 32-1360 Table 12, Family 7P — the genuine e=0.0549 member of
# the family (same anchor used by Task 2's test_er3bp_direct_seeding). #432 showed
# this family stays hyperbolic: NO bifurcation along the e-continuation.
MU_BROUCKE = 0.012155
BROUCKE_E0549_STATE0 = np.array([0.13892235, 0.0, 0.0, 0.0, 3.24428595, 0.0])

# Convention: parent_period_f is the integration span the corrector expects, i.e.
# the HALF-period (pi) when the corrector runs in half-period-residual mode — the
# same span passed by er3bp_discovery.continue_and_monitor / converge_direct_seed.
BROUCKE_HALF_PERIOD_F = math.pi


def test_branch_switcher_on_non_bifurcating_parent_is_safe() -> None:
    """On a non-bifurcating parent the switcher must not invent a false branch."""
    sys = ER3BPSystem(mu=MU_BROUCKE, e=0.0549, primary_name="E", secondary_name="M")

    result = branch_at_saddle_center_er3bp(sys, BROUCKE_E0549_STATE0, BROUCKE_HALF_PERIOD_F)

    # (b) returns the documented tuple[ER3BPPeriodicOrbit | None, dict]
    assert isinstance(result, tuple)
    assert len(result) == 2
    orbit, info = result
    assert isinstance(info, dict)
    assert orbit is None or isinstance(orbit, ER3BPPeriodicOrbit)

    # (c) if it returns an orbit, it must NOT be a spurious distinct family — it
    # must be within a small tolerance of the parent (i.e. the switcher did not
    # fabricate a fake branch on an orbit that does not actually bifurcate).
    if orbit is not None:
        np.testing.assert_allclose(orbit.state0, BROUCKE_E0549_STATE0, atol=1e-3)
    else:
        assert "reason" in info


def test_branch_switcher_returns_reason_when_no_branch() -> None:
    """When no marginal eigenvector / perturbation converges, info carries a reason."""
    sys = ER3BPSystem(mu=MU_BROUCKE, e=0.0549, primary_name="E", secondary_name="M")

    orbit, info = branch_at_saddle_center_er3bp(sys, BROUCKE_E0549_STATE0, BROUCKE_HALF_PERIOD_F)

    if orbit is None:
        assert info.get("reason") in {
            "no marginal eigenvector",
            "no perturbation converged",
        }
    else:
        # A converged orbit must report the perturbation that produced it.
        assert "epsilon" in info and "sign" in info and "eigenvalue" in info
