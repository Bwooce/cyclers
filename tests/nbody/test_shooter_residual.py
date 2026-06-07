"""N-body Phase C: defect residual is zero on a continuous fixture (plan Phase C).

NON-GOLDEN solver-health check: nodes sampled FROM a single n-body propagation
make the leg-continuity defects vanish by construction (a self-consistency check,
not a sourced anchor). The Jones-multiset gate is Task C.4.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import ShootingSeed, defect_residual  # noqa: E402


def _arc_consistent_seed(ephem: Ephemeris) -> ShootingSeed:
    """Build a 3-node seed whose two legs lie ON one continuous n-body arc.

    Sampling node states from the SAME propagation guarantees the leg-continuity
    defects are zero up to integrator accuracy. Flyby hinges are made trivially
    feasible (v_in == v_out -> no bend), and the wrap term is just measured, not
    asserted-zero (the synthetic arc is not periodic).
    """
    bodies = ("E", "M")
    prop = RestrictedNBody("rebound")
    # Seed near 1 AU on a mildly eccentric heliocentric arc.
    r0 = np.array([1.30e8, 0.0, 0.0])
    v0 = np.array([0.0, 26.0, 1.5])
    day = 86400.0
    t_nodes = [0.0, 120.0 * day, 260.0 * day]
    states: list[np.ndarray] = [np.concatenate([r0, v0])]
    cur_r, cur_v, cur_t = r0, v0, 0.0
    for t1 in t_nodes[1:]:
        arc = prop.propagate(
            cur_r, cur_v, t0_sec=cur_t, t1_sec=t1, bodies=bodies, accuracy=1e-11, ephem=ephem
        )
        states.append(np.concatenate([arc.r_km, arc.v_km_s]))
        cur_r, cur_v, cur_t = arc.r_km, arc.v_km_s, t1
    tofs = [120.0, 140.0]
    zero = np.zeros(3)
    return ShootingSeed(
        node_states=states,
        epochs=t_nodes,
        tofs=tofs,
        sequence=("E", "M", "E"),
        slack_leg=1,
        period_days=260.0,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )


@pytest.mark.slow
def test_defect_zero_for_self_consistent_nodes() -> None:
    ephem = Ephemeris("astropy")
    seed = _arc_consistent_seed(ephem)
    res = defect_residual(seed, ephem=ephem, bodies=("E", "M"), accuracy=1e-11)
    # Leg-continuity defects (first 12 components: 2 legs x 6) must vanish.
    leg_defects = res[:12]
    assert np.max(np.abs(leg_defects)) < 1e-3  # km / (km/s) defect floor
