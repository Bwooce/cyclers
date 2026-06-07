"""N-body Phase B: node-impulse correction-dV convention (plan Phase B; design §3 Q3)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.nbody.correction_dv import node_impulse_correction_dv


def test_zero_correction_when_seed_equals_corrected() -> None:
    nodes = {"e0": np.zeros(3), "m1": np.array([1.0, 0.0, 0.0])}
    dv = node_impulse_correction_dv(nodes, nodes)
    assert dv.total_kms == pytest.approx(0.0)


def test_correction_is_sum_of_per_node_discontinuity_changes() -> None:
    seed = {"e0": np.zeros(3), "m1": np.zeros(3)}
    corr = {"e0": np.array([0.1, 0.0, 0.0]), "m1": np.array([0.0, 0.05, 0.0])}
    dv = node_impulse_correction_dv(seed, corr)
    assert dv.total_kms == pytest.approx(0.15)
    assert set(dv.per_node_kms) == {"e0", "m1"}
