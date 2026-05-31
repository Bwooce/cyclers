"""Cross-tool consistency: Lambert solution propagates back to its target.

For each of the three canonical legs, solve Lambert for ``(v1, v2)``, then
propagate ``(r1, v1)`` forward by ``tof`` with the universal-variable Kepler
propagator. The final position must agree with ``r2`` to within 1 km and the
final velocity with ``v2`` to within ~0.1 m/s. This is the "cross-tool"
self-consistency check separate from the lamberthub-only agreement.

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §4.2.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert

from .conftest import Leg


def _reclose(leg: Leg) -> None:
    sol = lambert(leg.r1, leg.r2, leg.tof)[0]
    r_end, v_end = propagate(leg.r1, sol.v1, leg.tof)
    pos_err_km = float(np.linalg.norm(r_end - leg.r2))
    vel_err_km_s = float(np.linalg.norm(v_end - sol.v2))
    assert pos_err_km < 1.0, f"{leg.name}: |r_end - r2| = {pos_err_km:.6f} km"
    # Plan §4.2: < 1e-4 km/s.
    assert vel_err_km_s < 1.0e-4, f"{leg.name}: |v_end - v2| = {vel_err_km_s:.6e} km/s"


def test_aldrin_leg_kepler_reclose(leg_aldrin: Leg) -> None:
    """Aldrin E->M leg: Lambert solution re-closes under Kepler propagation."""
    _reclose(leg_aldrin)


def test_short_arc_kepler_reclose(leg_short: Leg) -> None:
    """Short Earth-to-Earth arc: same."""
    _reclose(leg_short)


def test_long_arc_kepler_reclose(leg_long: Leg) -> None:
    """Long Earth-to-Mars arc: same."""
    _reclose(leg_long)


@pytest.mark.parametrize("scale", [0.5, 1.0, 1.5])
def test_stumpff_handles_circular_propagation_scale(scale: float) -> None:
    """Stumpff/Kepler stable across short to long arcs of a circular orbit."""
    from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2

    r0 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    v_circ = float(np.sqrt(MU_SUN_KM3_S2 / AU_KM))
    v0 = np.array([0.0, v_circ, 0.0], dtype=np.float64)
    period_s = 2.0 * np.pi * float(np.sqrt(AU_KM**3 / MU_SUN_KM3_S2))
    r1, _ = propagate(r0, v0, scale * period_s)
    # After scale * period the orbit lands on the predicted point.
    theta = 2.0 * np.pi * scale
    expected = np.array([AU_KM * np.cos(theta), AU_KM * np.sin(theta), 0.0])
    assert float(np.linalg.norm(r1 - expected)) < 1.0
