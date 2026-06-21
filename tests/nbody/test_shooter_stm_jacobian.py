"""Analytic block-bidiagonal STM Jacobian == FD Jacobian (#388).

Parity gate for the ``jacobian="stm"`` path of :func:`cyclerfinder.nbody.shooter.shoot`.
The multiple-shooting residual is full-state continuity, so the Jacobian is exactly
known from per-leg state-transition matrices (one co-integrated variational
propagation per leg) — the lever that replaces the ``6*n_nodes+1`` finite-difference
re-propagations and makes the multi-year shoot tractable. FD stays the oracle: the
analytic STM Jacobian must agree with it. The fixture is the same cheap 3-node
Sun-only seed the parallel-FD test uses, so this runs in well under a second.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody import shooter  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootingSeed,
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _x_to_states,
    defect_residual,
    shoot,
)


def _two_body_seed() -> tuple[ShootingSeed, Ephemeris]:
    """A 3-node Sun-only seed sampled from one continuous two-body arc.

    Cheap and deterministic (``bodies = ()``), mirroring the parallel-FD fixture.
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


def _residual(
    seed: ShootingSeed, x: np.ndarray, ephem: Ephemeris, bodies: tuple[str, ...]
) -> np.ndarray:
    n = len(seed.sequence)
    trial = _seed_with_states(seed, _x_to_states(x, n))
    return defect_residual(trial, ephem=ephem, bodies=bodies, accuracy=1e-11)


@pytest.mark.slow
def test_stm_jacobian_matches_fd_jacobian() -> None:
    """The analytic STM Jacobian agrees with the FD oracle on the cheap fixture."""
    seed, ephem = _two_body_seed()
    bodies: tuple[str, ...] = ()
    x0 = _states_to_x(seed.node_states)
    f0 = _residual(seed, x0, ephem, bodies)
    fd = _fd_jacobian(
        lambda x: _residual(seed, x, ephem, bodies), x0, f0, column_eval=_serial_columns
    )
    stm_jac = shooter._stm_jacobian(seed, x0, ephem=ephem, bodies=bodies, accuracy=1e-11)
    assert stm_jac.shape == fd.shape
    rel = np.linalg.norm(stm_jac - fd) / np.linalg.norm(fd)
    print(f"STM Jacobian vs FD rel = {rel:.3e}")
    assert rel < 5e-3, f"STM Jacobian vs FD rel={rel}"


@pytest.mark.slow
def test_shoot_stm_reaches_same_fixed_point_as_fd() -> None:
    """``shoot(jacobian='stm')`` reaches the same fixed-point quality as ``'fd'``.

    The leg defects start ~0 (the seed is one continuous arc), so the LM solve
    minimises the structurally-irreducible periodicity wrap (this arc is not a
    closed cycler). With node 0 and node ``n-1`` both free the wrap minimiser is
    non-unique, so the exact corrected states differ between Jacobian modes — the
    invariant is the *defect-norm quality* of the fixed point: the analytic STM
    Jacobian drives the LM solve to the same minimum as the FD oracle (and at least
    as deep), which is what proves it steers the solver correctly.
    """
    seed, ephem = _two_body_seed()
    res_fd = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=40, jacobian="fd")
    res_stm = shoot(seed, ephem=ephem, bodies=(), accuracy=1e-11, max_nfev=40, jacobian="stm")
    print(
        f"seed={res_fd.seed_defect_norm:.3e} "
        f"defect_norm fd={res_fd.defect_norm:.3e} stm={res_stm.defect_norm:.3e}"
    )
    # Both reach the same fixed-point quality (within 50% of each other) ...
    rel = abs(res_stm.defect_norm - res_fd.defect_norm) / max(res_fd.defect_norm, 1e-12)
    assert rel < 0.5, (
        f"STM vs FD defect-norm quality differs: fd={res_fd.defect_norm} stm={res_stm.defect_norm}"
    )
    # ... and both meaningfully reduce the seed defect.
    assert res_stm.defect_norm < res_fd.seed_defect_norm
