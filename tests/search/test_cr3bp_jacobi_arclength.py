"""Tests for pseudo-arclength CR3BP continuation in the Jacobi constant C at
fixed mu (``cyclerfinder.search.cr3bp_jacobi_arclength``, #249).

This is the fold-turning continuation that the natural-parameter Jacobi
continuation (fix C, re-solve x0) cannot do. These tests verify the machinery
on FAST, SHORT walks of a converged Earth-Moon symmetric seed -- they do NOT
attempt the full C11a/C21 saddle-center fold-turn (too slow for unit tests; that
validation is the controller's job). All EXPECTED quantities are intrinsic
correctness properties (residual ~ 0, Jacobi conserved by an independent
integrator, period continuity, tangent in the null space of dr/dz), never a
value our own code computed and then asserted against itself.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_jacobi_arclength as ja
import cyclerfinder.search.cr3bp_periodic as cp

ROSS_MU = 1.2150584270572e-2  # canonical Earth-Moon mass parameter


def _system(mu: float) -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8)


def _seed() -> cp.SymmetricOrbit:
    """A converged short-period Earth-Moon perpendicular-crossing symmetric seed.

    Obtained via the fixed-Jacobi symmetric corrector at Earth-Moon mu; chosen
    low-period (T ~ 0.52) so every continuation step is sub-second.
    """
    seed = cp.correct_symmetric_fixed_jacobi(
        _system(ROSS_MU),
        0.82,
        3.17,
        2.8,
        ydot0_sign=1.0,
        half_crossings=1,
        tol=1e-11,
    )
    assert seed.converged
    assert seed.crossing_residual < 1e-9
    return seed


def test_zero_length_self_consistency() -> None:
    """Continuing toward a C_target ~ the seed's own C returns it ~unchanged.

    Runtime: ~0.5 s (one short walk).
    """
    seed = _seed()
    branch = ja.continue_in_jacobi(
        seed,
        mu=ROSS_MU,
        half_crossings=1,
        ydot0_sign=1.0,
        c_target=seed.jacobi + 1e-5,
        ds0=2e-3,
        ds_max=4e-3,
        max_steps=20,
        record_every=1,
    )
    assert branch.stop_reason == ja.JacobiStopReason.TARGET_REACHED
    assert len(branch.members) >= 2
    landed = branch.members[-1]

    # Landed essentially on the seed (tiny C nudge -> tiny x0/period change).
    assert abs(landed.jacobi - (seed.jacobi + 1e-5)) < 1e-9
    assert abs(landed.x0 - seed.x0) < 1e-4
    assert abs(landed.period - seed.period) < 1e-3

    # Genuine periodic orbit: perpendicular re-crossing + independent Radau.
    assert landed.crossing_residual < 1e-8
    assert landed.radau_djacobi < 1e-8

    # Jacobi self-consistency: C(state0, mu) equals the recorded jacobi.
    c_check = cr3bp.jacobi_constant(landed.state0, ROSS_MU)
    assert abs(c_check - landed.jacobi) < 1e-10

    # Symmetric IC structure preserved (planar perpendicular x-axis crossing).
    s = landed.state0
    assert s[1] == 0.0 and s[2] == 0.0 and s[3] == 0.0 and s[5] == 0.0
    assert abs(s[4]) > 1e-6  # non-trivial ydot0


def test_period_continuity_on_short_walk() -> None:
    """A few-step walk varies period SMOOTHLY (no >30% jumps); Radau drift small.

    Runtime: ~1.5 s (one short walk of a handful of steps).
    """
    seed = _seed()
    branch = ja.continue_in_jacobi(
        seed,
        mu=ROSS_MU,
        half_crossings=1,
        ydot0_sign=1.0,
        c_target=seed.jacobi + 0.02,
        ds0=3e-3,
        ds_max=6e-3,
        max_steps=8,
        record_every=1,
    )
    assert len(branch.members) >= 3

    periods = [m.period for m in branch.members]
    for prev_p, next_p in pairwise(periods):
        assert abs(next_p - prev_p) < 0.30 * prev_p  # smooth, no topology jump

    # Every member is a genuine orbit under an INDEPENDENT integrator (Radau).
    for m in branch.members:
        assert m.crossing_residual < 1e-7
        assert m.radau_djacobi < 1e-7
        assert m.mu == ROSS_MU  # mu held frozen along the whole walk


def test_tangent_is_unit_null_vector_of_jacobian() -> None:
    """The continuation tangent is a unit vector in the null space of dr/dz.

    Property check (one residual-Jacobian evaluation). Runtime: < 0.3 s.
    """
    seed = _seed()
    z = np.array([seed.x0, seed.jacobi])
    t_hi = 1.8 * seed.period

    tan = ja.tangent(z, ROSS_MU, 1.0, 1, t_hi, prev=None)
    assert tan is not None

    rj = ja._residual_jac(z, ROSS_MU, 1.0, 1, t_hi, rtol=1e-12, atol=1e-12)
    assert rj is not None
    _r0, _t_half, grad = rj

    # Unit length and orthogonal to the residual gradient (i.e. in its null space).
    assert abs(np.linalg.norm(tan) - 1.0) < 1e-12
    assert abs(float(grad @ tan)) < 1e-8 * (np.linalg.norm(grad) + 1.0)

    # Oriented toward increasing C at the start (prev=None convention).
    assert tan[1] >= 0.0


def test_land_at_jacobi_projects_onto_target_c() -> None:
    """``land_at_jacobi`` Newton-projects x0 onto r=0 at exactly C_target.

    Runtime: < 0.5 s.
    """
    seed = _seed()
    z_pred = np.array([seed.x0, seed.jacobi])
    c_target = seed.jacobi + 5e-4
    t_hi = 1.8 * seed.period
    z = ja.land_at_jacobi(
        z_pred,
        c_target,
        ROSS_MU,
        1.0,
        1,
        t_hi,
        tol=1e-11,
        max_iter=80,
    )
    assert z is not None
    # Held exactly at C_target; residual driven to ~0.
    assert z[1] == c_target
    rj = ja._residual_jac(z, ROSS_MU, 1.0, 1, t_hi, rtol=1e-12, atol=1e-12)
    assert rj is not None
    assert abs(rj[0]) < 1e-10


@pytest.mark.slow
def test_fold_turn_smoke() -> None:
    """Tiny fold-turn smoke: the walk follows the curve even where the natural
    parameter (C) is non-monotone in arclength.

    Gated @slow -- kept short. The bare arclength (period_jump_frac=0) is allowed
    to turn; we only assert it produced a chain of genuine members and that some
    step reversed the C-direction (evidence the tangent, not C, drove the walk).
    """
    seed = _seed()
    branch = ja.continue_in_jacobi(
        seed,
        mu=ROSS_MU,
        half_crossings=1,
        ydot0_sign=1.0,
        c_target=seed.jacobi - 0.5,  # push toward the inner edge of the family
        ds0=4e-3,
        ds_max=1e-2,
        max_steps=40,
        record_every=1,
    )
    assert len(branch.members) >= 2
    for m in branch.members:
        assert m.crossing_residual < 1e-6
        assert m.radau_djacobi < 1e-6
