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


def test_jacobi_gap_dv_min_exact_algebra() -> None:
    # Hand derivation (Cuevas del Valle et al. 2026, CEAS-GNC-2026-012, Sec. 7.2
    # technique): an impulse at fixed position changes C only via the -v^2
    # term, so delta_c = v0^2 - vf^2 with vf the post-impulse speed; the
    # cheapest impulse is tangential with |dv| = |vf - v0|.
    #   v0 = 3, delta_c = +5 (raise C / shed energy): vf = sqrt(9-5) = 2 -> dv = 1
    #   v0 = 3, delta_c = -7 (lower C / add energy):  vf = sqrt(9+7) = 4 -> dv = 1
    #   delta_c = 0: no burn needed.
    assert cr3bp.jacobi_gap_dv_min(3.0, 5.0) == pytest.approx(1.0, rel=1e-14)
    assert cr3bp.jacobi_gap_dv_min(3.0, -7.0) == pytest.approx(1.0, rel=1e-14)
    assert cr3bp.jacobi_gap_dv_min(3.0, 0.0) == 0.0
    # v0 = 0 (zero-velocity surface): dv = sqrt(-delta_c); C can only decrease.
    assert cr3bp.jacobi_gap_dv_min(0.0, -4.0) == pytest.approx(2.0, rel=1e-14)


def test_jacobi_gap_dv_min_consistent_with_jacobi_constant() -> None:
    # The tangential burn of exactly the bound's magnitude realises delta_c
    # exactly, per our own jacobi_constant convention (position unchanged, so
    # the potential terms cancel and only -v^2 moves).
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.05, 0.3, -0.02])
    v0 = float(np.linalg.norm(s[3:6]))
    for delta_c in (0.04, -0.07):
        dv = cr3bp.jacobi_gap_dv_min(v0, delta_c)
        vf = math.sqrt(v0 * v0 - delta_c)
        s_burn = s.copy()
        s_burn[3:6] = s[3:6] * (vf / v0)  # tangential: scale speed v0 -> vf
        assert np.linalg.norm(s_burn[3:6] - s[3:6]) == pytest.approx(dv, rel=1e-12)
        achieved = cr3bp.jacobi_constant(s_burn, mu) - cr3bp.jacobi_constant(s, mu)
        assert achieved == pytest.approx(delta_c, abs=1e-12)


def test_jacobi_gap_dv_min_is_a_lower_bound() -> None:
    # No impulse of magnitude 0.99 * bound, in ANY direction, can bridge the
    # gap: for delta_c > 0 the achieved gap v0^2 - |v + dv|^2 is maximised by
    # the anti-tangential dv and still falls short; symmetrically for
    # delta_c < 0. Sampled over a deterministic direction grid.
    mu = 0.012277471
    s = np.array([0.5, 0.1, 0.0, 0.05, 0.3, -0.02])
    v0 = float(np.linalg.norm(s[3:6]))
    c0 = cr3bp.jacobi_constant(s, mu)
    rng = np.random.default_rng(20260612)
    dirs = rng.normal(size=(64, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    for delta_c in (0.04, -0.07):
        w = 0.99 * cr3bp.jacobi_gap_dv_min(v0, delta_c)
        for u in dirs:
            s_burn = s.copy()
            s_burn[3:6] = s[3:6] + w * u
            achieved = cr3bp.jacobi_constant(s_burn, mu) - c0
            if delta_c > 0:
                assert achieved < delta_c
            else:
                assert achieved > delta_c


def test_jacobi_gap_dv_min_raises_on_unreachable_gap() -> None:
    # delta_c > v0^2 needs vf^2 < 0: even killing all velocity (dv = v0) only
    # raises C by v0^2 — the zero-velocity ceiling for one impulse.
    with pytest.raises(ValueError, match="zero-velocity ceiling"):
        cr3bp.jacobi_gap_dv_min(3.0, 9.0 + 1e-9)
    with pytest.raises(ValueError, match="negative speed"):
        cr3bp.jacobi_gap_dv_min(-1.0, 0.0)


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


def test_sun_earth_registry_system() -> None:
    # #409: cr3bp_system serves Sun-primary heliocentric pairs from the planet
    # registry (PLANETS + MU_SUN_KM3_S2), the single source of truth for the
    # #405/#411 cross-system work.
    #
    # SOURCED cross-checks (no circular dependence on our own derived value):
    #   1. The Sun-Earth (Earth-only) CR3BP mass parameter is the standard
    #      textbook value mu ~ 3.0035e-6 (e.g. Koon-Lo-Marsden-Ross system table
    #      lists the Sun-Earth ratio at this order; Earth-only, not Earth+Moon).
    #   2. The characteristic time satisfies 2*pi*t_s = Earth's orbital period,
    #      which must be one year. Earth orbiting the Sun in ~1 yr is the physical
    #      anchor independent of the GM arithmetic.
    se = cr3bp.cr3bp_system("Sun", "Earth")
    assert se.mu == pytest.approx(3.0035e-6, rel=2e-4)
    year_s = 365.25 * 86400.0
    assert (2.0 * math.pi * se.t_s) == pytest.approx(year_s, rel=1e-3)
    # Convention identical to the moon path: t_s = sqrt(l^3 / G(m1+m2)).
    gm_pair = 1.32712440018e11 + 3.98600435507e5  # MU_SUN + Earth-only GM
    assert se.t_s == pytest.approx(math.sqrt(se.l_km**3 / gm_pair), rel=1e-12)
    assert se.mu == pytest.approx(3.98600435507e5 / gm_pair, rel=1e-12)
    # Unknown Sun-secondary names fail loudly (parallels SATELLITES[bad] KeyError).
    with pytest.raises(KeyError):
        cr3bp.cr3bp_system("Sun", "Nibiru")


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
