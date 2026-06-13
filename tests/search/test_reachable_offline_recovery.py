"""Offline (network-independent) member recovery at C_J=3.1294 (#247).

The free-(x0, t_half) corrector
:func:`cyclerfinder.search.reachable_representatives.correct_symmetric_free_period`
recovers Braik-Ross 2026 (arXiv:2605.31543) Table-2 representatives without the
JPL oracle. The seed (x0 region + velocity sign + target half-period) is the only
family input; the recovered period is a PREDICTION confirmed against the sourced
Table-2 period AND the sourced Floquet rate sigma. All EXPECTED values trace to
Braik-Ross Table 2; recovered ICs are derived quantities, never goldens.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.reachable_representatives as rr


def test_free_period_holds_jacobi_and_closes() -> None:
    """Structural: the recovered IC sits on C_J and re-closes (independent Radau)."""
    sysm = rr.braik_ross_system()
    o = correct = rr.correct_symmetric_free_period(
        sysm, 0.4485, rr.C_J_BRAIK_ROSS, 0.5 * 26.500 / rr.TU_DAYS, ydot0_sign=1.0
    )
    assert o.converged
    assert o.jacobi == pytest.approx(rr.C_J_BRAIK_ROSS, abs=1e-9)
    assert o.crossing_residual < 1e-9
    po = cp.PeriodicOrbit(
        state0=np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0]),
        period=o.period,
        jacobi=o.jacobi,
        converged=o.converged,
        closure_residual=o.crossing_residual,
    )
    ok, dj = cp.crosscheck_periodic(sysm, po, closure_tol=1e-3, jacobi_tol=1e-8)
    assert ok, f"independent Radau cross-check failed (dj={dj:.2e})"
    _ = correct


def test_free_period_separates_r31s_from_r21s() -> None:
    """The free-period corrector recovers R31-S (27.252 d), NOT R21-S (26.500 d).

    The 1-DOF fixed-crossing corrector collapses this seed region onto R21-S; the
    free-period corrector, seeded at the R31-S half-period, recovers the distinct
    27.25 d member -- the branch-separation capability that motivated it.
    """
    sysm = rr.braik_ross_system()
    o = rr.correct_symmetric_free_period(
        sysm, 0.357, rr.C_J_BRAIK_ROSS, 0.5 * 27.252 / rr.TU_DAYS, ydot0_sign=-1.0
    )
    assert o.converged
    assert o.period * rr.TU_DAYS == pytest.approx(27.252, abs=0.2)
    assert abs(o.period * rr.TU_DAYS - 26.500) > 0.5


def test_lagrange_collinear_points() -> None:
    """L1/L2 collinear-point x positions (structural; no network)."""
    mu = rr.ROSS_MU
    l1 = rr.lagrange_collinear_x(mu, "L1")
    l2 = rr.lagrange_collinear_x(mu, "L2")
    # Standard Earth-Moon collinear points (Szebehely): L1~0.8369, L2~1.1557.
    assert l1 == pytest.approx(0.8369, abs=1e-3)
    assert l2 == pytest.approx(1.1557, abs=1e-3)
    # dUbar/dx vanishes there.
    assert abs(cp._ubar_grad_x_at_axis(l1, mu)) < 1e-8
    assert abs(cp._ubar_grad_x_at_axis(l2, mu)) < 1e-8


@pytest.mark.slow
def test_offline_recovery_confirms_sourced_periods() -> None:
    """Each offline-seeded member recovers to its Table-2 period AND sigma (no JPL)."""
    sysm = rr.braik_ross_system()
    reps = rr.recover_offline_set(sysm)
    by = {r.label: r for r in reps}
    for r in reps:
        src_d, src_s = rr.SOURCED_TABLE2[r.label]
        arc = cr3bp.propagate(sysm, r.state0, r.period, with_stm=False)
        print(
            f"{r.label:6s} conv={r.converged!s:5s} T={r.period_days:8.3f} d "
            f"(sourced {src_d:7.3f}, sigma {src_s}) C={r.jacobi:.6f} "
            f"confirmed={r.confirmed}"
        )
        _ = arc
    # These five recover offline to <0.5 d of the sourced period with the right
    # stability character (sigma): LL1, LL2, DPO, R21-S, R31-S.
    for label in ("LL1", "LL2", "DPO", "R21-S", "R31-S"):
        r = by[label]
        assert r.confirmed, (
            f"{label}: recovered {r.period_days:.3f} d vs sourced "
            f"{r.sourced_period_days} d (converged={r.converged})"
        )
        assert r.jacobi == pytest.approx(rr.C_J_BRAIK_ROSS, abs=1e-9)


def test_sourced_table2_matches_sourced_periods_subset() -> None:
    """The full Table-2 dict is consistent with the legacy SOURCED_PERIODS_DAYS subset."""
    for label, days in rr.SOURCED_PERIODS_DAYS.items():
        assert rr.SOURCED_TABLE2[label][0] == pytest.approx(days, abs=1e-9)
    # sigma = 0 exactly for the stable resonants.
    for label in ("R21-S", "R31-S", "R52-S"):
        assert rr.SOURCED_TABLE2[label][1] == 0.0
    assert math.isclose(rr.SOURCED_TABLE2["LL1"][1], 2.4884)
