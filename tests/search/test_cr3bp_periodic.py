"""CR3BP periodic-orbit corrector (plan 2026-06-10, Phase 2)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# Sourced Arenstorf golden (μ, x0, vy0, period) — Arenstorf 1963 / Hairer et al.
# (Hairer, Nørsett, Wanner, "Solving ODEs I", p. 129, test problem B5)
MU = 0.012277471
X0, VY0, PERIOD = 0.994, -2.0015851063790825, 17.0652165601579625


def test_arenstorf_orbit_is_periodic() -> None:
    # The published Arenstorf IC already (very nearly) closes after one period —
    # the corrector confirms it and returns a tight closure residual.
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0])
    res = cp.correct_periodic(sysm, s0, PERIOD)
    assert res.converged
    assert res.closure_residual < 1e-6
    assert res.period == pytest.approx(PERIOD, rel=1e-2)


def test_non_periodic_guess_does_not_converge() -> None:
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([0.7, 0.0, 0.0, 0.0, -0.2, 0.0])
    res = cp.correct_periodic(sysm, s0, 3.0, max_iter=8)
    assert not res.converged


def test_cuevas2026_em_l1_southern_halo_seed_converges() -> None:
    # SOURCED SEED: Cuevas del Valle, Urrutxua & Solano-Lopez 2026 ("Fuel-optimal
    # Rendezvous in the CR3BP via MPC and Proximal Operators", CEAS EuroGNC 2026,
    # CEAS-GNC-2026-012), Sec. 7.3 p. 18 — the Scenario II target state, printed
    # at full precision. SEED, not golden: the paper prints no mu and no Jacobi
    # constant, and the only period is the CHASER's, "approximately 2.7549" nd
    # (4 s.f., nearby family member, opposite z sign).
    #
    # LABEL ADJUDICATION: the paper calls this an "L2 southern halo", but the
    # label is wrong — x0 = 0.8240 is Earth-side of L1 (x_L1 = 0.836915 for the
    # physical EM mu; the L2 family has no y = 0 crossing there), and the same
    # group's 2023 Aerospace paper (10.3390/aerospace10050393, Sec. 3.1 p. 12)
    # labels the nearly identical x = 0.82413 state an Earth-Moon L1 standard
    # halo. See docs/notes/2026-06-11-cuevas-del-valle-2023-floquet-mining.md
    # Sec. 3.3. Treated here as an Earth-Moon L1 SOUTHERN halo seed.
    sysm = cr3bp.cr3bp_system("Earth", "Moon")  # physical mu = 0.0121505844
    s0 = np.array([0.824024728136525, 0.0, -0.054501847320725, 0.0, 0.164671964079122, 0.0])
    res = cp.correct_periodic(sysm, s0, 2.7549)

    # Converged with a small correction off the printed seed (recorded run
    # 2026-06-13: residual 7.7e-13, |ds| = 2.3e-3 nd, T = 2.76021 nd,
    # C = 3.15169, monodromy |lambda|_max ~ 1.58e3 i.e. strongly unstable).
    assert res.converged
    assert res.closure_residual < 1e-10
    assert float(np.linalg.norm(res.state0 - s0)) < 5e-3
    # Period lands within 1% of the published ~2.7549 nd (the chaser's
    # approximate period on the adjacent northern member; 0.19% observed).
    assert res.period == pytest.approx(2.7549, rel=1e-2)

    # Independent-integrator crosscheck (Radau vs DOP853).
    ok, dj = cp.crosscheck_periodic(sysm, res)
    assert ok
    assert dj < 1e-10

    # Topology confirms the L1 relabel: the orbit straddles x_L1 = 0.836915
    # (root of dUbar/dx between Earth and Moon for this mu) and never nears the
    # Moon at x = 1 - mu, so it cannot be an L2-family member; the dominant z
    # excursion is southern (z_min = -0.0530 vs z_max = +0.0436 observed).
    x_l1 = 0.836915
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, res.period),
        res.state0,
        args=(sysm.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.sol is not None
    samples = sol.sol(np.linspace(0.0, res.period, 400))
    xs, zs = samples[0], samples[2]
    assert xs.min() < x_l1 < xs.max()
    assert xs.max() < 1.0 - sysm.mu - 0.05  # well Earth-side of the Moon/L2
    assert abs(zs.min()) > abs(zs.max())  # southern member
    # Strongly unstable, as expected for an EM L1 halo of this size: the
    # monodromy's largest eigenvalue is O(1e3) (nu = (lambda + 1/lambda)/2 >> 1).
    arc = cr3bp.propagate(sysm, res.state0, res.period, with_stm=True)
    assert arc.stm is not None
    lam_max = float(np.max(np.abs(np.linalg.eigvals(arc.stm))))
    assert lam_max > 100.0


def test_periodic_orbit_crosscheck_independent_integrator() -> None:
    # A corrected Arenstorf orbit re-propagated with a DIFFERENT integrator ("Radau")
    # stays closed and conserves Jacobi -- an independent confirmation.
    # Sourced golden: Arenstorf 1963 / Hairer et al. (same as above).
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    res = cp.correct_periodic(sysm, np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0]), PERIOD)
    ok, dc = cp.crosscheck_periodic(sysm, res)
    assert ok
    assert dc < 1e-8
