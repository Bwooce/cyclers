from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootingSeed,
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _stm_jacobian,
    _x_to_states,
    defect_residual,
    shoot,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    """A 3-node Sun-only seed sampled from one continuous two-body arc.

    Mirrors the fixture in tests/nbody/test_shooter_stm_jacobian.py: cheap,
    deterministic, bodies=() so the penalty (which depends only on node velocity
    vs planet velocity from the ephemeris) is exercised without perturber cost.
    """
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


def test_vinf_penalty_off_is_identical() -> None:
    seed, ephem = _two_body_seed()
    base = defect_residual(seed, ephem=ephem, bodies=())
    none_set = defect_residual(seed, ephem=ephem, bodies=(), vinf_anchors=None, vinf_weight=0.0)
    np.testing.assert_array_equal(base, none_set)
    # anchors present but zero weight -> still identical (no rows appended)
    zero_w = defect_residual(
        seed, ephem=ephem, bodies=(), vinf_anchors={"E": 5.0, "M": 6.0}, vinf_weight=0.0
    )
    np.testing.assert_array_equal(base, zero_w)


def test_vinf_penalty_rows_values() -> None:
    seed, ephem = _two_body_seed()
    w = 4.0
    anchors = {"E": 5.0, "M": 6.0}
    base = defect_residual(seed, ephem=ephem, bodies=())
    res = defect_residual(seed, ephem=ephem, bodies=(), vinf_anchors=anchors, vinf_weight=w)
    pen = res[len(base) :]
    # sequence is E, M, E -> all three nodes carry an anchor -> 3 penalty rows
    assert pen.shape == (3,)
    sw = float(np.sqrt(w))
    expected = []
    for i, body in enumerate(seed.sequence):
        _, v_pl = ephem.state(body, seed.epochs[i])
        mag = float(np.linalg.norm(np.asarray(seed.node_states[i][3:]) - np.asarray(v_pl)))
        expected.append(sw * (mag - anchors[body]))
    np.testing.assert_allclose(pen, np.asarray(expected), rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_stm_jacobian_with_penalty_matches_fd() -> None:
    """The augmented STM Jacobian (incl. penalty rows) matches the FD oracle."""
    seed, ephem = _two_body_seed()
    bodies: tuple[str, ...] = ()
    anchors = {"E": 5.0, "M": 6.0}
    w = 4.0
    x0 = _states_to_x(seed.node_states)

    def resid(x: np.ndarray) -> np.ndarray:
        trial = _seed_with_states(seed, _x_to_states(x, len(seed.sequence)))
        return defect_residual(
            trial, ephem=ephem, bodies=bodies, accuracy=1e-11, vinf_anchors=anchors, vinf_weight=w
        )

    f0 = resid(x0)
    fd = _fd_jacobian(resid, x0, f0, column_eval=_serial_columns)
    stm = _stm_jacobian(
        seed, x0, ephem=ephem, bodies=bodies, accuracy=1e-11, vinf_anchors=anchors, vinf_weight=w
    )
    assert stm.shape == fd.shape
    rel = np.linalg.norm(stm - fd) / np.linalg.norm(fd)
    print(f"penalty STM Jacobian vs FD rel = {rel:.3e}")
    assert rel < 5e-3, f"penalty STM Jacobian vs FD rel={rel}"


def test_shoot_penalty_off_matches_plain() -> None:
    """vinf_weight=0 makes the penalty layer a true no-op.

    The penalty block is gated on ``vinf_anchors and vinf_weight > 0.0``
    (shooter.py block 4), so at ``vinf_weight=0.0`` both calls run the IDENTICAL
    residual/Jacobian path — bit-identical single-threaded. We assert *numerical*
    identity rather than exact ``==`` because under multithreaded BLAS (CI) the
    two least-squares solves carry reduction-order noise (~1e-10 rel on the
    objective). The tolerances below sit far above that noise yet far below any
    real penalty effect (which, if the weight=0 gate failed, would shift the
    solution by O(km/s)), so the no-op invariant is tested robustly.
    """
    seed, ephem = _two_body_seed()
    plain = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=20, jacobian="stm")
    off = shoot(
        seed,
        ephem=ephem,
        bodies=(),
        accuracy=1e-11,
        max_nfev=20,
        jacobian="stm",
        vinf_anchors={"E": 5.0, "M": 6.0},
        vinf_weight=0.0,
    )
    assert plain.defect_norm == pytest.approx(off.defect_norm, rel=1e-6)
    for a, b in zip(plain.corrected_states, off.corrected_states, strict=True):
        np.testing.assert_allclose(a, b, rtol=1e-4, atol=1e-2)


def test_shoot_penalty_rejects_parallel_fd() -> None:
    """A positive penalty with the parallel FD path is unsupported -> ValueError."""
    seed, ephem = _two_body_seed()
    with pytest.raises(ValueError, match="vinf penalty"):
        shoot(
            seed,
            ephem=ephem,
            bodies=(),
            accuracy=1e-11,
            max_nfev=5,
            n_jobs=4,
            vinf_anchors={"E": 5.0},
            vinf_weight=1.0,
        )


@pytest.mark.slow
def test_shoot_penalty_biases_vinf_toward_anchor() -> None:
    """A strong penalty pulls a node's corrected V∞ toward the anchor.

    Target the WRAP node (index n-1, body E): node 0's V∞ is locked by leg-0
    continuity + the wrap term (a continuous-arc seed has no slack there, so the
    penalty correctly cannot move it), but the wrap node trades leg-defect ↔ wrap
    ↔ penalty and so responds to the anchor pull.
    """
    seed, ephem = _two_body_seed()
    wrap = len(seed.sequence) - 1  # node 2, body "E"
    base = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=40, jacobian="stm")
    base_e = base.vinf_per_encounter_kms[wrap]
    target = base_e + 5.0  # pull the wrap-node V∞ 5 km/s away from natural
    pinned = shoot(
        seed,
        ephem=ephem,
        bodies=(),
        accuracy=1e-11,
        max_nfev=40,
        jacobian="stm",
        vinf_anchors={"E": target},
        vinf_weight=50.0,
    )
    # The penalized solve's wrap-node V∞ moves toward the target vs the baseline.
    assert abs(pinned.vinf_per_encounter_kms[wrap] - target) < abs(base_e - target)
