"""N-body Phase C: multiple-shooting node/defect structure + conic seed map (plan Phase C).

Reads node V_inf VECTORS from correct._vinf_nodes (b{i}_in/_out), NOT from a
best_cycler attribute (design §3 drift: BallisticClosureResult carries only scalar
vinf_per_encounter_kms; the vectors live in _vinf_nodes).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.shooter import build_shooting_vector, defect_count


def test_defect_count_is_legs_times_state_dim() -> None:
    # E-M-E-V-V-E: 5 legs -> 5 full-state defects (6 components each).
    assert defect_count(n_encounters=6) == 5 * 6


def test_shooting_vector_packs_states_epochs_tofs() -> None:
    nodes = {f"b{i}": np.arange(6.0) for i in range(4)}
    epochs = [0.0, 1.0, 2.0, 3.0]
    tofs = [1.0, 1.0, 1.0]
    x = build_shooting_vector(nodes, epochs, tofs, slack_leg=2, period_days=4.0)
    # slack leg eliminated from the free vector.
    assert len(x) == 4 * 6 + 4 + (3 - 1)
