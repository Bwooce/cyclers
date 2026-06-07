"""N-body perf lever 1 (#159): parallel FD Jacobian == serial FD Jacobian.

Correctness gate for the ``n_jobs`` parallel-finite-difference path of
:func:`cyclerfinder.nbody.shooter.shoot`. The parallel Jacobian evaluates the
independent FD columns across a process pool; the arithmetic (step size, forward
difference) is identical to the serial reference, so the two Jacobians must agree
to working precision. This fixture is deliberately small and Sun-only (``bodies
= ()``, two-body) so it runs in well under a second and is fully deterministic.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootingSeed,
    _fd_jacobian,
    _parallel_columns_for_test,
    _serial_columns,
    _states_to_x,
    _x_to_states,
    defect_residual,
    shoot,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    """A 3-node Sun-only seed sampled from one continuous two-body arc.

    Cheap and deterministic: no rails perturbers (``bodies = ()``), so each leg
    propagation is a single REBOUND chunk. The circular ephemeris supplies the
    encounter-body states for the wrap term but is irrelevant to the leg defects.
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


def test_parallel_fd_jacobian_equals_serial() -> None:
    seed, ephem = _two_body_seed()
    n = len(seed.sequence)

    def residual_of_x(x: np.ndarray) -> np.ndarray:
        trial = ShootingSeed(
            node_states=_x_to_states(x, n),
            epochs=list(seed.epochs),
            tofs=list(seed.tofs),
            sequence=seed.sequence,
            slack_leg=seed.slack_leg,
            period_days=seed.period_days,
            vinf_in=list(seed.vinf_in),
            vinf_out=list(seed.vinf_out),
        )
        return defect_residual(trial, ephem=ephem, bodies=(), accuracy=1e-11)

    x0 = _states_to_x(seed.node_states)
    f0 = residual_of_x(x0)

    jac_serial = _fd_jacobian(residual_of_x, x0, f0, column_eval=_serial_columns)
    jac_parallel = _parallel_columns_for_test(seed, ephem, x0, f0, n_jobs=4)

    # Same FD arithmetic, different execution mode -> equal to working precision.
    assert jac_serial.shape == jac_parallel.shape
    np.testing.assert_allclose(jac_parallel, jac_serial, rtol=1e-9, atol=1e-12)


def test_shoot_n_jobs_reduces_defect() -> None:
    """A bounded solve with n_jobs>1 corrects the seed (solver health).

    The column-equality gate above proves the parallel Jacobian is numerically
    identical to the serial one of the SAME FD scheme. A full ``shoot`` with
    ``n_jobs>1`` swaps ``scipy``'s internal-FD ``lm`` path for our explicit-FD
    ``lm`` path (a different FD scheme from ``n_jobs == 1``, so the LM iterates
    legitimately differ at FP level); the invariant that holds end-to-end is that
    the parallel solve drives the defect below the seed, i.e. it is solving.
    """
    seed, ephem = _two_body_seed()
    # Nudge the interior node off the arc so there is something to correct.
    states = [s.copy() for s in seed.node_states]
    states[1][3:] += np.array([0.02, -0.02, 0.01])
    perturbed = ShootingSeed(
        node_states=states,
        epochs=seed.epochs,
        tofs=seed.tofs,
        sequence=seed.sequence,
        slack_leg=seed.slack_leg,
        period_days=seed.period_days,
        vinf_in=seed.vinf_in,
        vinf_out=seed.vinf_out,
    )
    parallel = shoot(perturbed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=15, n_jobs=4)
    assert parallel.defect_norm < parallel.seed_defect_norm
