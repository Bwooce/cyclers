"""The node-impulse correction-ΔV metric (design §3, Q3).

The SILVER-rung headline number. The **correction ΔV** is the sum of the
impulsive ΔV at the encounter nodes that the flybys could *not* absorb
gravitationally and that a real mission would burn to restore n-body periodicity
from the patched-conic seed. Concretely it is the per-node change in the
spacecraft velocity discontinuity between the raw conic seed and the
n-body-corrected solution.

**Comparable to ``maintain.py``'s per-synodic maintenance ΔV** (``maintain.py:24,
41``; the 2.9138 km/s the V4 consumer reproduces): the same physical quantity —
ΔV to keep the cycle periodic — with a *documented* mapping, NOT forced
numerically identical (design Q3). The accounting is **sign- and node-explicit**.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Vec3 = NDArray[np.float64]


@dataclass(frozen=True)
class CorrectionDV:
    """Frozen node-impulse correction-ΔV accounting (design §3, Q3)."""

    total_kms: float
    per_node_kms: dict[str, float]
    per_node_vector_kms: dict[str, Vec3]

    def compare_to_maintenance(self, maintenance_dv_kms: float) -> dict[str, float]:
        """Document the mapping to ``maintain.py``'s per-synodic maintenance ΔV.

        Comparable (same physical meaning) but not identical (design Q3): the
        maintenance ΔV is the per-synodic burn for the *powered* Aldrin cycler,
        whereas this is the per-period node-impulse to restore n-body periodicity
        of a *ballistic* seed. Returns both numbers + their ratio / difference so
        the gauntlet can compare like with like across consumers 1 and 3.
        """
        return {
            "rung_correction_dv_kms": self.total_kms,
            "maintenance_dv_kms": float(maintenance_dv_kms),
            "difference_kms": self.total_kms - float(maintenance_dv_kms),
            "ratio": (
                self.total_kms / float(maintenance_dv_kms) if maintenance_dv_kms else float("inf")
            ),
        }


def node_impulse_correction_dv(
    seed_nodes: Mapping[str, Vec3],
    corrected_nodes: Mapping[str, Vec3],
) -> CorrectionDV:
    """Sum the per-node velocity-discontinuity changes (seed -> corrected).

    For each shared node, the correction is ``|v_corrected - v_seed|`` (the
    impulsive ΔV applied at that node); the total is their sum. Nodes present in
    only one mapping contribute their own magnitude (a node the correction
    introduced / removed). Sign- and node-explicit: the per-node *vector* change
    is retained alongside its magnitude.
    """
    keys = set(seed_nodes) | set(corrected_nodes)
    per_node_kms: dict[str, float] = {}
    per_node_vector_kms: dict[str, Vec3] = {}
    total = 0.0
    for key in sorted(keys):
        v_seed = np.asarray(seed_nodes.get(key, np.zeros(3)), dtype=np.float64)
        v_corr = np.asarray(corrected_nodes.get(key, np.zeros(3)), dtype=np.float64)
        delta = v_corr - v_seed
        mag = float(np.linalg.norm(delta))
        per_node_vector_kms[key] = delta
        per_node_kms[key] = mag
        total += mag
    return CorrectionDV(
        total_kms=total,
        per_node_kms=per_node_kms,
        per_node_vector_kms=per_node_vector_kms,
    )


__all__ = ["CorrectionDV", "node_impulse_correction_dv"]
