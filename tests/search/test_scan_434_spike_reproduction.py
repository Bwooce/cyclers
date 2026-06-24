"""#434 spike-reproduction golden — the #287 3D Braik-Ross (1,1) family.

Re-runs the #287 spike's converged 3D C11a seed through the production
``continue_general_3d_family`` pseudo-arclength tracer and asserts the family
reproduces. This is a CLOSURE-SELF-CHECK golden (per the orbit-closure /
golden-sourced discipline): the expected side is the spike's published-on-disk
result (``data/spike_287.jsonl``) — the family locks near z0 = -0.2408 with
~157 members all closing to ~2.7e-9, and we assert against safe bounds
(>=75 members, closure < 1e-8) rather than an exact member count, because the
pseudo-arclength walk's fold-traversal makes the realized count step-dependent.

The seed is the #287 spike's ``z0_seed_5.0000e-02`` 3D lock (the moderate-z0
amplitude rung that locked onto the 3D branch, NOT a planar collapse).
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.search.cr3bp_3d_family_tracer import continue_general_3d_family
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
)

# #287 spike converged 3D C11a (1,1) Earth-Moon seed (data/spike_287.jsonl,
# case "z0_seed_5.0000e-02" — the moderate-z0 lock onto the 3D branch).
_SPIKE_X0 = -0.8116406668238195
_SPIKE_Z0 = -0.2408102083477011
_SPIKE_YDOT0 = -0.10629710963669947
_SPIKE_PERIOD = 10.204301970414399


@pytest.mark.slow
def test_spike_287_3d_family_reproduces() -> None:
    system = cr3bp_system("Earth", "Moon")
    seed_state = [_SPIKE_X0, 0.0, _SPIKE_Z0, 0.0, _SPIKE_YDOT0, 0.0]

    family = continue_general_3d_family(
        system,
        seed_state,
        _SPIKE_PERIOD,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=100,
        direction="both",
        monodromy_eval=False,
    )

    # Spike walked to ~201 members through folds; assert the safe lower bound.
    assert len(family.members) >= 75, f"expected >=75 converged members, got {len(family.members)}"

    # Every accepted member must close independently. Spike worst case ~2.7e-9;
    # 1e-8 is the safe bound.
    worst_closure = max(m.orbit.independent_closure_residual for m in family.members)
    assert worst_closure < 1e-8, (
        f"worst independent closure residual {worst_closure:.3e} exceeds 1e-8"
    )

    # The family stays genuinely 3D (no collapse back to the planar manifold).
    for m in family.members:
        assert not m.orbit.degenerate_planar, (
            f"member at step {m.step_index} collapsed to the planar manifold "
            f"(z0={m.orbit.state0[2]:.3e})"
        )
