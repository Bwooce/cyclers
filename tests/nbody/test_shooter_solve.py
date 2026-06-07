"""N-body Phase C: multiple-shooting solve drives defects down (plan Phase C).

NON-GOLDEN convergence check: asserts the solver reduces the defect from a
slightly-perturbed self-consistent seed (a solver-health test, not a sourced
rediscovery). The Jones-multiset gate is Task C.4.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import ShootingSeed, shoot  # noqa: E402


def _consistent_seed(ephem: Ephemeris) -> ShootingSeed:
    """A 3-node seed whose legs lie on one continuous n-body arc (defects ~0)."""
    bodies = ("E", "M")
    prop = RestrictedNBody("rebound")
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
    zero = np.zeros(3)
    return ShootingSeed(
        node_states=states,
        epochs=t_nodes,
        tofs=[120.0, 140.0],
        sequence=("E", "M", "E"),
        slack_leg=1,
        period_days=260.0,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )


def _perturbed_consistent_seed(ephem: Ephemeris) -> ShootingSeed:
    """Self-consistent seed with the interior node velocity nudged off the arc."""
    seed = _consistent_seed(ephem)
    states = [s.copy() for s in seed.node_states]
    states[1][3:] += np.array([0.02, -0.02, 0.01])  # ~30 m/s interior nudge
    return ShootingSeed(
        node_states=states,
        epochs=seed.epochs,
        tofs=seed.tofs,
        sequence=seed.sequence,
        slack_leg=seed.slack_leg,
        period_days=seed.period_days,
        vinf_in=seed.vinf_in,
        vinf_out=seed.vinf_out,
    )


@pytest.mark.slow
def test_shoot_reduces_defect_from_perturbed_seed() -> None:
    ephem = Ephemeris("astropy")
    # Bounded nfev: a solver-health check (defect goes DOWN), not a full converge.
    result = shoot(
        _perturbed_consistent_seed(ephem),
        ephem=ephem,
        bodies=("E", "M"),
        accuracy=1e-10,
        max_nfev=25,
    )
    assert result.defect_norm < result.seed_defect_norm
