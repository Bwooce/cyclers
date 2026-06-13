"""Tests for the mass-parameter (mu) continuation of CR3BP symmetric cyclers
(``cyclerfinder.search.mu_continuation``).

The acceptance anchor is the held Ross & Roberts-Tsoukkas 2025 (3,1) stable
Earth-Moon cycler at mu = 1.2150584270572e-2: a small step in mu must recover a
genuine periodic orbit (perpendicular re-crossing residual ~ machine precision,
independent-Radau full-period re-closure, Jacobi conserved) that varies
*continuously* from the seed -- the "a known small-mu step recovers a known
member" deliverable.

All EXPECTED quantities are intrinsic correctness properties (residual ~ 0,
Jacobi conserved, period continuous), never a value our own code computed and
then asserted against itself.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.mu_continuation as mc

ROSS_MU = 1.2150584270572e-2


def _system(mu: float) -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8)


def _seed_31() -> cp.SymmetricOrbit:
    # Ross Table 3 (3,1): C^stable, T^stable (sourced); x0 seed is derived.
    return cp.correct_symmetric_fixed_jacobi(
        _system(ROSS_MU),
        -0.3209891696,
        3.161784147013429,
        14.78849241668140,
        ydot0_sign=-1.0,
        half_crossings=3,
        tol=1e-11,
    )


def test_seed_is_genuine_periodic_orbit() -> None:
    """The starting member is a real perpendicular-crossing periodic orbit."""
    seed = _seed_31()
    assert seed.converged
    assert seed.crossing_residual < 1e-9
    # Jacobi enforced exactly by construction.
    assert abs(seed.jacobi - 3.161784147013429) < 1e-12


def test_small_mu_step_recovers_genuine_member() -> None:
    """A small mu step yields a genuine periodic orbit, continuous from the seed."""
    seed = _seed_31()
    branch = mc.continue_in_mu(
        seed,
        ROSS_MU,
        half_crossings=3,
        ydot0_sign=-1.0,
        mu_target=ROSS_MU + 0.01,
        label="(3,1)",
        record_every=1,
        ds0=4e-3,
        ds_max=8e-3,
    )
    assert branch.stop_reason == mc.MuStopReason.TARGET_REACHED
    assert len(branch.members) >= 2
    landed = branch.members[-1]

    # Reached (close to) the target mu and it MOVED from the seed mu.
    assert landed.mu > ROSS_MU + 5e-3
    assert abs(landed.mu - (ROSS_MU + 0.01)) < 5e-3

    # Genuine periodic orbit: perpendicular re-crossing + independent Radau.
    assert landed.crossing_residual < 1e-8
    assert landed.radau_djacobi < 1e-7

    # Jacobi self-consistency: C(state0, mu) equals the recorded jacobi.
    c_check = cr3bp.jacobi_constant(landed.state0, landed.mu)
    assert abs(c_check - landed.jacobi) < 1e-10

    # Symmetric IC structure preserved (planar perpendicular x-axis crossing).
    s = landed.state0
    assert s[1] == 0.0 and s[2] == 0.0 and s[3] == 0.0 and s[5] == 0.0
    assert abs(s[4]) > 1e-6  # non-trivial ydot0 (not an equilibrium)

    # Continuity: the family parameters drift smoothly (no topology jump).
    assert abs(landed.period - seed.period) < 0.3 * seed.period
    assert abs(landed.x0 - seed.x0) < 0.2


def test_independent_radau_closure_on_landed_member() -> None:
    """The landed member re-closes under a DIFFERENT integrator (Radau vs DOP853)."""
    seed = _seed_31()
    branch = mc.continue_in_mu(
        seed,
        ROSS_MU,
        half_crossings=3,
        ydot0_sign=-1.0,
        mu_target=ROSS_MU + 0.006,
        label="(3,1)",
        record_every=1,
        ds0=3e-3,
        ds_max=6e-3,
    )
    landed = branch.members[-1]
    system = _system(landed.mu)
    po = cp.PeriodicOrbit(
        state0=landed.state0,
        period=landed.period,
        jacobi=landed.jacobi,
        converged=True,
        closure_residual=landed.crossing_residual,
    )
    ok, dj = cp.crosscheck_periodic(system, po, closure_tol=1e-2, jacobi_tol=1e-7)
    assert ok
    assert dj < 1e-7


def test_c_scan_at_seed_mu_contains_a_stable_member() -> None:
    """The fixed-mu C-scan around the held (3,1) member finds the stable subfamily.

    The held member is itself the published nu=0 midpoint, so a scan about it
    must contain at least one |nu|<1 member -- exercises the inner stable-window
    search used for figure-matching at binary-star mu.
    """
    seed = _seed_31()
    members = mc.scan_c_family_at_mu(
        ROSS_MU,
        seed.x0,
        seed.jacobi,
        seed.period,
        half_crossings=3,
        ydot0_sign=-1.0,
        dc=5e-5,
        n_each=6,
    )
    assert members, "C-scan returned no converged members"
    assert any(m.stable for m in members), "no stable member near the published nu=0 midpoint"
    # Every returned member is a genuine periodic orbit.
    for m in members:
        assert m.crossing_residual < 1e-7
        assert np.isfinite(m.nu)
