"""Ross & Roberts-Tsoukkas 2025 (AAS 25-621) stable Earth-Moon cycler families.

Same-model golden gate for #212 Part B: the fixed-Jacobi symmetric corrector
(``correct_symmetric_fixed_jacobi``) + Barden half-period monodromy stability
(``barden_stability``) reproduce the paper's published stable (k1,k2) cycler
members.

SOURCED-GOLDEN DISCIPLINE: every EXPECTED value below is a number *printed* in
Ross & Roberts-Tsoukkas 2025 — the mass ratio mu (p. 3), and the C^stable /
T^stable / stability-verdict columns of Table 3 (p. 11). The recovered ``x0``
and ``ydot0`` are DERIVED quantities (the 1-D solve of §5), never goldens; they
are only checked for self-consistency (the orbit closes perpendicularly and the
Jacobi constant is preserved). The acceptance checks are:

  * Jacobi constant: enforced algebraically -> matches C^stable to ~machine eps;
  * period T = 2 * t_half: matches the published T^stable (limited by the
    finite width of the stable subfamily -- the exact nu=0 x0 sits a hair off
    the published 16-digit C, so T lands within ~3e-4 of T^stable for the two
    razor-thin windows, and to <1e-7 for the wide ones);
  * stability verdict: |nu| < 1 (linearly stable) AND nu ~ 0 (the Table-3
    member is the nu=0 midpoint of the largest stable subfamily, p. 11) -- the
    CR3BP lane's first |nu|<1 verdicts.

All five published families reproduce. (3,2)'s half-period perpendicular
crossing is the 6th x-axis crossing (``half_crossings=6``); seeded from the
same x0 ~ -0.321 region as (3,1)/(3,3) it lands T within ~2e-8 of T^stable
with nu ~ -0.012.

Model: pure planar CR3BP (PCR3BP); mu = 1.2150584270572e-2 (paper p. 3).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# Ross & Roberts-Tsoukkas 2025, AAS 25-621, p. 3.
ROSS_MU = 1.2150584270572e-2
ROSS_TU_DAYS = 27.321661 / (2.0 * math.pi)  # p. 3 -> 4.348377401631 d
ROSS_C1_L1 = 3.188341105401253  # p. 13, C(L1)
# Table 1 (p. 8) Earth-side critical Jacobi constants C^{u1}_{k1}.
ROSS_TABLE1_CU1 = {1: 3.151763728314920, 2: 3.129751730201047, 3: 3.188341092440989}


def _em_system() -> cr3bp.CR3BPSystem:
    # Paper's nd scales: a_m = 384,400 km; 1 TU = 27.321661 d / 2pi = 4.348377... d.
    return cr3bp.CR3BPSystem(
        mu=ROSS_MU,
        primary="Earth",
        secondary="Moon",
        l_km=384400.0,
        t_s=375699.8,
    )


# Per-family acceptance record (sourced EXPECTED on the left of each check).
# (label, C^stable, T^stable, x0_seed, ydot0_sign, half_crossings, T_tol)
# Source for C^stable/T^stable: Ross & Roberts-Tsoukkas 2025, Table 3, p. 11.
_FAMILIES = [
    ("(1,1)", 3.151175879508174, 10.29206921007976, -0.7682140805, -1.0, 3, 1e-5),
    ("(2,1)", 3.129389531088256, 19.44043166795154, 0.7237335857, 1.0, 4, 5e-4),
    ("(3,1)", 3.161784147013429, 14.78849241668140, -0.3209891696, -1.0, 3, 5e-4),
    ("(3,2)", 3.182762663084288, 17.90058010350006, -0.3210000000, -1.0, 6, 1e-5),
    ("(3,3)", 3.177224018696528, 18.14546057589189, -0.3217380626, -1.0, 5, 1e-5),
]


@pytest.mark.parametrize(
    ("label", "c_stable", "t_stable", "x0_seed", "sign", "n_half", "t_tol"),
    _FAMILIES,
    ids=[f[0] for f in _FAMILIES],
)
def test_ross_family_reproduced(
    label: str,
    c_stable: float,
    t_stable: float,
    x0_seed: float,
    sign: float,
    n_half: int,
    t_tol: float,
) -> None:
    sysm = _em_system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        sysm, x0_seed, c_stable, t_stable, ydot0_sign=sign, half_crossings=n_half, tol=1e-10
    )
    # The corrector converged to a perpendicular x-axis crossing.
    assert orbit.converged, f"{label}: corrector did not converge"

    # Jacobi constant is enforced -> matches the published C^stable exactly.
    assert orbit.jacobi == pytest.approx(c_stable, abs=1e-12)

    # Period matches the published T^stable (within the stable-window width).
    assert orbit.period == pytest.approx(t_stable, abs=t_tol)

    # Barden stability: the lane's first |nu| < 1 verdict, and the Table-3 member
    # is the nu=0 midpoint of the stable subfamily.
    nu, _lam = cp.barden_stability(sysm, orbit)
    assert abs(nu) < 1.0, f"{label}: expected linearly STABLE (|nu|<1), got nu={nu}"
    assert abs(nu) < 0.2, f"{label}: expected nu~0 midpoint, got nu={nu}"


@pytest.mark.parametrize(
    ("label", "c_stable", "t_stable", "x0_seed", "sign", "n_half", "t_tol"),
    _FAMILIES,
    ids=[f[0] for f in _FAMILIES],
)
def test_ross_family_independent_crosscheck(
    label: str,
    c_stable: float,
    t_stable: float,
    x0_seed: float,
    sign: float,
    n_half: int,
    t_tol: float,
) -> None:
    # Independent-integrator (Radau) confirmation that the recovered member is a
    # genuine periodic orbit: re-propagate one full period and require closure +
    # Jacobi conservation. (Degeneracy gate: a perpendicular crossing alone is
    # not enough; the full-period orbit must close.)
    sysm = _em_system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        sysm, x0_seed, c_stable, t_stable, ydot0_sign=sign, half_crossings=n_half, tol=1e-10
    )
    po = cp.PeriodicOrbit(
        state0=np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0]),
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        closure_residual=orbit.crossing_residual,
    )
    ok, dj = cp.crosscheck_periodic(sysm, po, closure_tol=1e-3, jacobi_tol=1e-8)
    assert ok, f"{label}: independent Radau cross-check failed (dj={dj:.2e})"


def test_ydot0_from_jacobi_inverts_jacobi_constant() -> None:
    # ydot0_from_jacobi(x0, C) must invert jacobi_constant for a perpendicular IC.
    mu = ROSS_MU
    x0, c = -0.7682140805, 3.151175879508174
    ydot0 = cp.ydot0_from_jacobi(x0, c, mu, sign=-1.0)
    state = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    assert cr3bp.jacobi_constant(state, mu) == pytest.approx(c, abs=1e-13)


def test_ydot0_from_jacobi_rejects_unattainable_jacobi() -> None:
    # A Jacobi constant above the -2*Ubar ceiling has no real perpendicular IC.
    with pytest.raises(ValueError, match="negative radicand"):
        cp.ydot0_from_jacobi(-0.768, 100.0, ROSS_MU)


def test_barden_monodromy_is_symplectic() -> None:
    # Structural check (independent of any published number): the 4x4 planar
    # Barden monodromy must be symplectic (det = 1).
    sysm = _em_system()
    _label, c, t, x0, sign, n_half, _tt = _FAMILIES[0]
    orbit = cp.correct_symmetric_fixed_jacobi(
        sysm, x0, c, t, ydot0_sign=sign, half_crossings=n_half
    )
    s0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    phi = cr3bp.propagate(sysm, s0, orbit.t_half, with_stm=True).stm
    assert phi is not None
    idx = [0, 1, 3, 4]
    g4 = np.diag([1.0, -1.0, -1.0, 1.0])
    mono = g4 @ np.linalg.inv(phi[np.ix_(idx, idx)]) @ g4 @ phi[np.ix_(idx, idx)]
    # det = 1 to integrator precision (the half-period STM over the ~10 TU (1,1)
    # orbit accumulates ~1e-7 symplecticity error at the DOP853 default tols).
    assert np.linalg.det(mono) == pytest.approx(1.0, abs=1e-6)


def test_tu_days_matches_published() -> None:
    # EXPECTED: Ross p. 3 prints 1 TU = 4.348377401631 d.
    assert pytest.approx(4.348377401631, abs=5e-13) == ROSS_TU_DAYS


def test_period_day_conversions_match_table4() -> None:
    # EXPECTED: Table 4 (p. 13) day column for the (k1,1) rows.
    expected_days = {
        (3.151175879508174, 10.29206921007976): 44.753800,
        (3.129389531088256, 19.44043166795154): 84.534335,
        (3.161784147013429, 14.78849241668140): 64.305944,
    }
    for _label, c, t, *_ in _FAMILIES:
        if (c, t) in expected_days:
            assert t * ROSS_TU_DAYS == pytest.approx(expected_days[(c, t)], abs=5e-4)


def test_all_c_stable_below_c1() -> None:
    # EXPECTED: C(L1) = 3.188341105401253 (p. 13); every C^stable sits below it.
    for _label, c, *_ in _FAMILIES:
        assert c < ROSS_C1_L1


def test_data_gap_c21_bound_inconsistency() -> None:
    # data_gap (do NOT silently resolve): Table 3 prints C_(2,1)=3.1297495000000
    # but Eq.8 + Table 1 give min{C^{u1}_2, C^{u2}_1}=3.129751730201047
    # (mining note section 4 item 1). Assert the ~2.2e-6 inconsistency persists.
    eq8_bound = min(ROSS_TABLE1_CU1[2], 3.1833333078762)  # C^{u2}_1 from Table 2
    table3_bound = 3.1297495000000
    assert abs(eq8_bound - table3_bound) == pytest.approx(2.235e-6, abs=1e-7)


def test_data_gap_c31_bound_inconsistency() -> None:
    # data_gap (do NOT silently resolve): Table 3 prints C_(3,1)=3.1833333078762
    # but Table 4's dC_(3,1)=1.272710e-2 implies C_(3,1)=C1-dC=3.175614005
    # (mining note section 4 item 2). Assert the ~7.7e-3 discrepancy persists.
    table3_bound = 3.1833333078762
    table4_implied = ROSS_C1_L1 - 1.272710e-2
    assert abs(table3_bound - table4_implied) == pytest.approx(7.72e-3, abs=5e-5)


def test_barden_stability_detects_unstable_member() -> None:
    # Sanity: the Arenstorf figure-eight (an UNSTABLE orbit) yields |nu| > 1.
    # This pins the verdict direction so a stable-only test cannot pass vacuously.
    # Sourced: Arenstorf 1963 / Hairer et al. (test-problem mu = 0.012277471).
    sysm = cr3bp.CR3BPSystem(mu=0.012277471, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    # The Arenstorf IC is itself a perpendicular x-axis crossing (x0, 0, 0, 0, vy0, 0)
    # with a known half period = T/2; build the SymmetricOrbit directly.
    orbit = cp.SymmetricOrbit(
        x0=0.994,
        ydot0=-2.0015851063790825,
        jacobi=cr3bp.jacobi_constant(
            np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0]), 0.012277471
        ),
        t_half=0.5 * 17.0652165601579625,
        period=17.0652165601579625,
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    nu, _lam = cp.barden_stability(sysm, orbit)
    assert abs(nu) > 1.0
