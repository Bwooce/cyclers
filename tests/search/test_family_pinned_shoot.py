from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import ShootingSeed  # noqa: E402
from cyclerfinder.search.family_pinned_shoot import (  # noqa: E402
    FamilyPinnedResult,
    family_pinned_shoot,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    ephem = Ephemeris("circular")
    prop = RestrictedNBody("rebound")
    r0 = np.array([1.30e8, 0.0, 0.0])
    v0 = np.array([0.0, 26.0, 1.5])
    day = 86400.0
    t_nodes = [0.0, 120.0 * day, 260.0 * day]
    states: list[np.ndarray] = [np.concatenate([r0, v0])]
    cur_r, cur_v, cur_t = r0, v0, 0.0
    for t1 in t_nodes[1:]:
        arc = prop.propagate(cur_r, cur_v, t0_sec=cur_t, t1_sec=t1, bodies=(), accuracy=1e-11)
        states.append(np.concatenate([arc.r_km, arc.v_km_s]))
        cur_r, cur_v, cur_t = arc.r_km, arc.v_km_s, t1
    zero = np.zeros(3)
    seed = ShootingSeed(
        node_states=states,
        epochs=t_nodes,
        tofs=[120.0, 140.0],
        sequence=("E", "M", "E"),
        slack_leg=1,
        period_days=260.0,
        vinf_in=[zero, zero, zero],
        vinf_out=[zero, zero, zero],
    )
    return seed, ephem


def test_weight_ladder_must_end_at_zero() -> None:
    seed, ephem = _two_body_seed()
    with pytest.raises(ValueError, match="end at 0"):
        family_pinned_shoot(
            seed,
            ephem=ephem,
            bodies=(),
            vinf_anchors={"E": 6.0},
            weight_ladder=(10.0, 1.0),  # does not end at 0.0
        )


@pytest.mark.slow
def test_family_pinned_shoot_returns_result_and_trace() -> None:
    seed, ephem = _two_body_seed()
    res = family_pinned_shoot(
        seed,
        ephem=ephem,
        bodies=(),
        vinf_anchors={"E": 6.0, "M": 6.0},
        weight_ladder=(10.0, 1.0, 0.0),
        accuracy=1e-11,
        max_nfev=20,
    )
    assert isinstance(res, FamilyPinnedResult)
    # final rung is the unpenalized (weight=0) solve
    assert res.final_weight == 0.0
    assert res.final.sequence == seed.sequence
    # one trace entry per ladder rung
    assert [w for (w, _d, _v) in res.trace] == [10.0, 1.0, 0.0]
    # anchor_retention is finite
    assert np.isfinite(res.anchor_retention_kms)


@pytest.mark.slow
def test_family_pinned_shoot_penalty_pulls_during_ladder() -> None:
    """The penalty wired through the driver pulls wrap-node V∞ toward the anchor.

    Compared at the FIRST (most-penalized) ladder rung vs the plain unpenalized
    shoot: the penalized rung's wrap-node V∞ is closer to a deliberately-far
    target. (The final λv=0 rung then relaxes; the wrap-node minimizer is not
    unique at this nfev budget, so the unpenalized endpoint itself is not a tight
    invariant — that retention behaviour is what the batch characterises.)
    """
    from cyclerfinder.nbody.shooter import shoot

    seed, ephem = _two_body_seed()
    wrap = len(seed.sequence) - 1  # node 2, body E, has slack
    plain = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=20, jacobian="stm")
    target = plain.vinf_per_encounter_kms[wrap] + 5.0
    res = family_pinned_shoot(
        seed,
        ephem=ephem,
        bodies=(),
        vinf_anchors={"E": target},
        weight_ladder=(50.0, 0.0),
        accuracy=1e-11,
        max_nfev=40,
    )
    # trace[0] is the most-penalized rung (w=50): its wrap V∞ is pulled toward target.
    _w0, _d0, vinf0 = res.trace[0]
    assert abs(vinf0[wrap] - target) < abs(plain.vinf_per_encounter_kms[wrap] - target)
