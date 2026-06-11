"""CR3BP dynamics core (plan 2026-06-10, Phase 1)."""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.satellites as sats


def test_jacobi_constant_value() -> None:
    # At a sample state the Jacobi constant matches the closed-form convention
    # C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2.
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0])
    x, y, z, vx, vy, vz = s
    r1 = math.hypot(x + mu, y, z)
    r2 = math.hypot(x - 1 + mu, y, z)
    expect = (x * x + y * y) + 2 * (1 - mu) / r1 + 2 * mu / r2 - (vx * vx + vy * vy + vz * vz)
    assert np.isclose(cr3bp.jacobi_constant(s, mu), expect)


def test_eom_shape_and_coriolis_sign() -> None:
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0])
    d = cr3bp.cr3bp_eom(0.0, s, mu)
    assert d.shape == (6,)
    # d[0:3] == velocity; ax includes +2*vy Coriolis term.
    assert np.allclose(d[0:3], s[3:6])


def test_jacobi_conserved_over_propagation() -> None:
    # Jacobi is conserved to ~1e-10 over a propagation (the integrator self-check).
    mu = 0.012277471
    s0 = np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0])
    c0 = cr3bp.jacobi_constant(s0, mu)
    arc = cr3bp.propagate(
        cr3bp.CR3BPSystem(mu=mu, primary="test", secondary="test", l_km=1.0, t_s=1.0), s0, 5.0
    )
    c1 = cr3bp.jacobi_constant(arc.state_f, mu)
    assert abs(c1 - c0) < 1e-9


def test_stm_matches_finite_difference() -> None:
    # The propagated 6x6 STM matches a finite-difference of the flow.
    mu = 0.012277471
    s0 = np.array([0.9, 0.0, 0.0, 0.0, -0.5, 0.0])
    sysm = cr3bp.CR3BPSystem(mu=mu, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    arc = cr3bp.propagate(sysm, s0, 1.0, with_stm=True)
    eps = 1e-6
    col0 = (
        cr3bp.propagate(sysm, s0 + np.array([eps, 0, 0, 0, 0, 0]), 1.0).state_f
        - cr3bp.propagate(sysm, s0 - np.array([eps, 0, 0, 0, 0, 0]), 1.0).state_f
    ) / (2 * eps)
    assert arc.stm is not None
    assert np.allclose(arc.stm[:, 0], col0, atol=1e-4)


def test_earth_moon_mu_physical() -> None:
    # cr3bp_system derives the physical Earth-Moon mass ratio from the registry
    # GMs (distinct from the Arenstorf test-problem 0.012277471).
    #
    # SOURCED expected value: Ross & Roberts-Tsoukkas 2025 (AAS 25-621) p.3
    # prints the Earth-Moon CR3BP mass parameter mu = 1.2150584270572e-2.
    # The registry-derived GM_Moon / GM_EarthSystem agrees to ~6e-9 relative;
    # the pre-fix formula GM_Moon / (GM_EarthSystem + GM_Moon) double-counted
    # the Moon (PRIMARIES["Earth"] is the Earth+Moon SYSTEM GM) and was off by
    # -1.2e-2 relative (mu ~0.0120047) — #212 Part A.
    assert "Moon" in sats.SATELLITES
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    assert sysm.mu == pytest.approx(1.2150584270572e-2, rel=1e-6)
    # Characteristic time uses the same G(m1+m2): t_s = sqrt(l^3 / GM_system).
    # The double-count made t_s ~0.6% too small (sqrt(1 + mu') ~ 1.00606).
    gm_system = sats.PRIMARIES["Earth"]
    assert sysm.t_s == pytest.approx(math.sqrt(sysm.l_km**3 / gm_system), rel=1e-12)


def test_saturnian_mu_double_count_was_negligible() -> None:
    # Quantify the SAME double-count effect (#212 Part A) for the Saturnian
    # pairs backing the 14 SILVER Lyapunov members. The relative shift between
    # the old formula gm2/(gm_sys+gm2) and the corrected gm2/gm_sys is exactly
    # mu itself: Mimas ~6.60e-8, Enceladus ~1.90e-7, Tethys ~1.09e-6 — orders
    # of magnitude below the members' reported precision, so the SILVER results
    # are unaffected (unlike Earth-Moon, where mu ~1.2e-2 made the shift -1.2%).
    expected_mu = {"Mimas": 6.599e-8, "Enceladus": 1.901e-7, "Tethys": 1.0864e-6}
    for moon, approx_mu in expected_mu.items():
        gm2 = sats.SATELLITES[moon].mu_km3_s2
        gm_sys = sats.PRIMARIES["Saturn"]
        sysm = cr3bp.cr3bp_system("Saturn", moon)
        assert sysm.mu == pytest.approx(approx_mu, rel=1e-3)
        rel_shift = abs(sysm.mu - gm2 / (gm_sys + gm2)) / sysm.mu
        assert rel_shift < 1.1e-6  # worst case (Tethys); EM was 1.2e-2


def test_propagate_raises_on_integrator_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    # propagate() must NOT silently return sol.y[:, -1] when solve_ivp fails
    # (it would be the state where the integrator gave up, not at time t).
    # A reliable REAL failure trigger was probed (start at the primary/secondary
    # singularity, radial collision plunges, sub-eps tolerances): DOP853 either
    # succeeds or grinds for minutes shrinking steps near the singularity rather
    # than promptly returning success=False, so no fast deterministic trigger is
    # constructible. The contract is pinned by monkeypatching solve_ivp instead.

    class _FailedSol:
        success = False
        message = "Required step size is less than spacing between numbers."
        t = np.array([0.0, 0.5])
        y = np.zeros((6, 2))

    def _fake_solve_ivp(*args: object, **kwargs: object) -> _FailedSol:
        return _FailedSol()

    monkeypatch.setattr(cr3bp, "solve_ivp", _fake_solve_ivp)
    sysm = cr3bp.CR3BPSystem(mu=0.012277471, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0])
    with pytest.raises(RuntimeError, match="step size"):
        cr3bp.propagate(sysm, s0, 1.0)
    with pytest.raises(RuntimeError, match="step size"):
        cr3bp.propagate(sysm, s0, 1.0, with_stm=True)
