"""Tests for the CR3BP natural-parameter (Jacobi) continuation driver.

Spec: ``docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md``
(Phase 1 testing block).

SOURCED-GOLDEN DISCIPLINE: the only EXPECTED numbers traceable to a publication
are Ross & Roberts-Tsoukkas 2025 (AAS 25-621): the mass ratio mu (p. 3) and the
(k1,k2) families' C^stable / T^stable / STABLE verdict (Table 3, p. 11). The
binding golden here is the spec's: *walking from a Ross seed reproduces the
published stable-window edge* -- i.e. the seed (the published nu=0 midpoint) is
linearly STABLE, and continuing in Jacobi walks OUT of that finite stable window
into instability (|nu| crosses 1). The stable-window EXISTENCE and the seed's
stability verdict are Ross's; the per-member nu values we compute are derived
diagnostics, not goldens.

The fault-injection tests (equilibrium / period-collapse) and the fold-stop test
use no published numbers -- they check the gauntlet/fold logic structurally.

Model: pure planar CR3BP (PCR3BP); mu = 1.2150584270572e-2 (Ross p. 3).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp

ROSS_MU = 1.2150584270572e-2


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=ROSS_MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
    )


def _seed_33() -> cp.SymmetricOrbit:
    # Ross (3,3): wide stable window (Delta_p_m 2041 km, 5 windows). The seed is
    # the published nu=0 midpoint (Table 3): C^stable, T^stable.
    return cp.correct_symmetric_fixed_jacobi(
        _em_system(),
        -0.3217380626,
        3.177224018696528,
        18.14546057589189,
        ydot0_sign=-1.0,
        half_crossings=5,
        tol=1e-10,
    )


def test_continuation_walks_a_branch_with_passing_gauntlet() -> None:
    # Walking from the (3,3) seed must produce an ordered branch of members, each
    # passing the full gauntlet (the independent-Radau dJ is recorded and small).
    sysm = _em_system()
    seed = _seed_33()
    assert seed.converged
    branch = cc.continue_family(
        sysm,
        seed,
        direction=-1,
        d_jacobi=1e-4,
        n_steps=8,
        min_jacobi=3.10,
        max_jacobi=3.1883,
        half_crossings=5,
        ydot0_sign=-1.0,
        seed_label="(3,3)",
    )
    # Seed + at least a few stepped members.
    assert len(branch.members) >= 4
    # Members are ordered in C by the step direction (decreasing here).
    cs = [m.jacobi for m in branch.members]
    assert all(cs[i + 1] < cs[i] for i in range(len(cs) - 1))
    # Every kept member passed the independent-Radau cross-check (small dJ) and
    # carries a Barden nu.
    for m in branch.members:
        assert m.radau_djacobi < 1e-8
        assert np.isfinite(m.nu)
        assert m.crossing_residual < 1e-6


def test_seed_is_stable_and_continuation_reproduces_published_stable_window_edge() -> None:
    # GOLDEN (spec): Ross's (3,3) Table-3 member is the nu=0 midpoint of a finite
    # STABLE window. The seed must be linearly STABLE (|nu|<1, Ross's verdict),
    # and walking in Jacobi must exit that finite window -- some continued member
    # becomes UNSTABLE (|nu| > 1). Reproducing the *edge* = the stable->unstable
    # transition the published window implies.
    sysm = _em_system()
    seed = _seed_33()
    seed_nu, _lam = cp.barden_stability(sysm, seed)
    assert abs(seed_nu) < 1.0, f"Ross (3,3) seed must be STABLE, got nu={seed_nu}"

    # Walk far enough (both directions) that the window edge is crossed.
    saw_stable = False
    saw_unstable = False
    for direction in (1, -1):
        branch = cc.continue_family(
            sysm,
            seed,
            direction=direction,
            d_jacobi=1e-4,
            n_steps=15,
            min_jacobi=3.10,
            max_jacobi=3.1883,
            half_crossings=5,
            ydot0_sign=-1.0,
            seed_label="(3,3)",
        )
        for m in branch.members:
            saw_stable = saw_stable or m.stable
            saw_unstable = saw_unstable or (not m.stable)
    assert saw_stable, "seed (the nu=0 midpoint) must register STABLE"
    assert saw_unstable, "continuation must walk OUT of the finite stable window"


def test_direction_must_be_plus_or_minus_one() -> None:
    sysm = _em_system()
    seed = _seed_33()
    with pytest.raises(ValueError, match="direction must be"):
        cc.continue_family(
            sysm,
            seed,
            direction=0,
            d_jacobi=1e-4,
            n_steps=1,
            min_jacobi=3.0,
            max_jacobi=3.2,
            half_crossings=5,
            ydot0_sign=-1.0,
        )


def test_stops_at_jacobi_bound() -> None:
    # A tight bound just past the seed forces an immediate JACOBI_BOUND stop.
    sysm = _em_system()
    seed = _seed_33()
    c_seed = seed.jacobi
    branch = cc.continue_family(
        sysm,
        seed,
        direction=1,
        d_jacobi=1e-4,
        n_steps=20,
        min_jacobi=3.0,
        max_jacobi=c_seed + 1.5e-4,  # room for ~1 step, then the bound
        half_crossings=5,
        ydot0_sign=-1.0,
    )
    assert branch.stop_reason is cc.StopReason.JACOBI_BOUND
    assert all(m.jacobi <= c_seed + 1.5e-4 for m in branch.members)


def _collinear_l1_x(mu: float) -> float:
    """x of L1 by Newton on dU/dx = 0 from a Hill-radius start (x-axis equilibrium)."""
    x = 1.0 - mu - (mu / 3.0) ** (1.0 / 3.0)
    for _ in range(100):
        s1 = x + mu
        s2 = x - 1.0 + mu
        f = x - (1.0 - mu) * s1 / abs(s1) ** 3 - mu * s2 / abs(s2) ** 3
        fp = 1.0 + 2.0 * (1.0 - mu) / abs(s1) ** 3 + 2.0 * mu / abs(s2) ** 3
        dx = f / fp
        x -= dx
        if abs(dx) < 1e-15:
            break
    return float(x)


def test_gauntlet_rejects_equilibrium() -> None:
    # Fault injection: a collinear libration point (on the x-axis, ydot0 = 0) is an
    # equilibrium -- it trivially "closes" for any period with max|v| ~ 0 and zero
    # amplitude. The corrector's symmetric IC is (x0, 0, 0, 0, ydot0, 0); putting
    # x0 = x_L1 and ydot0 = 0 reproduces the equilibrium the gate must reject.
    sysm = _em_system()
    x_l1 = _collinear_l1_x(sysm.mu)
    state0 = np.array([x_l1, 0.0, 0.0, 0.0, 0.0, 0.0])
    fake = cp.SymmetricOrbit(
        x0=x_l1,
        ydot0=0.0,
        jacobi=cr3bp.jacobi_constant(state0, sysm.mu),
        t_half=1.0,
        period=2.0,
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    ok, reason, member = cc._run_gauntlet(
        sysm,
        fake,
        period_floor=0.5,
        period_ceiling=10.0,
        max_speed_floor=cc.MAX_SPEED_FLOOR_ND,
        amplitude_floor=cc.AMPLITUDE_FLOOR_ND,
        jacobi_tol=cc.JACOBI_CONSERVATION_TOL,
        radau_closure_tol=1e-3,
        radau_jacobi_tol=1e-8,
        rtol=1e-12,
        atol=1e-12,
    )
    assert not ok
    assert reason == "equilibrium"
    assert member is None


def test_gauntlet_rejects_period_collapse() -> None:
    # Fault injection: a collapsed period (below the floor) must be rejected by
    # the period-bounds gate before any propagation cost.
    sysm = _em_system()
    seed = _seed_33()
    collapsed = cp.SymmetricOrbit(
        x0=seed.x0,
        ydot0=seed.ydot0,
        jacobi=seed.jacobi,
        t_half=1e-4,
        period=2e-4,  # far below period_floor
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    ok, reason, member = cc._run_gauntlet(
        sysm,
        collapsed,
        period_floor=0.5 * seed.period,
        period_ceiling=2.0 * seed.period,
        max_speed_floor=cc.MAX_SPEED_FLOOR_ND,
        amplitude_floor=cc.AMPLITUDE_FLOOR_ND,
        jacobi_tol=cc.JACOBI_CONSERVATION_TOL,
        radau_closure_tol=1e-3,
        radau_jacobi_tol=1e-8,
        rtol=1e-12,
        atol=1e-12,
    )
    assert not ok
    assert reason == "period_bounds"
    assert member is None


def test_gauntlet_passes_genuine_seed() -> None:
    # Positive control: the genuine Ross (3,3) seed passes every gate.
    sysm = _em_system()
    seed = _seed_33()
    ok, reason, member = cc._run_gauntlet(
        sysm,
        seed,
        period_floor=0.5 * seed.period,
        period_ceiling=2.0 * seed.period,
        max_speed_floor=cc.MAX_SPEED_FLOOR_ND,
        amplitude_floor=cc.AMPLITUDE_FLOOR_ND,
        jacobi_tol=cc.JACOBI_CONSERVATION_TOL,
        radau_closure_tol=1e-3,
        radau_jacobi_tol=1e-8,
        rtol=1e-12,
        atol=1e-12,
    )
    assert ok, f"genuine seed rejected: {reason}"
    assert member is not None
    assert member.stable  # Ross (3,3) is STABLE (Table 3)
    assert member.jacobi == pytest.approx(3.177224018696528, abs=1e-12)
