"""#896: the elliptic three-body model against the Sun-Mercury ME-Halo orbits of Peng, Bai & Xu.

Source: H. Peng, X. Bai and S. Xu, "Continuation of periodic orbits in the Sun-Mercury elliptic
restricted three-body problem", Communications in Nonlinear Science and Numerical Simulation 47
(2017) 1-15, doi:10.1016/j.cnsns.2016.11.005; filed in the private paper corpus as
peng-xu-2017-continuation-periodic-orbits-sun-mercury-ERTBP-cnsns-47-1.pdf.

Model (p. 2, Eqs. 1-4): pulsating synodic frame, true anomaly ``f`` as the independent variable,
``Omega~ = (x^2 + y^2)/2 + (1-mu)/r1 + mu/r2 + mu(1-mu)/2 - e cos f z^2 / 2`` divided by
``1 + e cos f``. ``core/er3bp.py`` uses the same potential (without the constant term) and the
same independent variable, so the paper's velocities (derivatives with respect to ``f``) are the
module's state. Constants: ``mu = 1.660e-7`` "approximately" (p. 3) and ``e = 0.2056`` (p. 1,
Figs. 2-5 and 17); these printed values are used.

Checked here:

* Table 1 (p. 4): six circular-problem halo orbits whose period is ``2 N pi / M`` (Eq. 5) for
  M:N = 5:2, 7:3, 9:4 around L1 and L2, states ``(x0, z0, ydot0)`` on the x-z plane. A symmetric
  corrector with the half period fixed at ``N pi / M`` and ``e = 0``, started from the printed
  state, must land within one unit of the last printed digit (measured: largest difference
  6.7e-7, in z0 of the 5:2 L1 row).
* Table 2 (p. 14): eight orbits of the elliptic problem with ``e = 0.2056`` and full period
  ``4 pi`` (M5N2), starting at ``f0 = 0`` or ``pi`` on the x-z plane. The same corrector (half
  period ``2 pi`` from ``f0``) must land within one unit of the last printed digit. Three rows
  do (L1 orbits 2 and 3, L2 orbit 3; largest difference 1.3e-6 in the five-decimal x0 of L2
  orbit 3). The other five do not converge to the printed orbit from the printed state with a
  damped single-shooting Newton corrector: L1 orbit 1, L1 orbit 4 and L2 orbit 1 stall (half
  period residual 0.04 to 0.11), and L2 orbits 2 and 4 converge to other orbits (x0 = 1.00514
  and 1.00505). With ``mu = 1.6601209e-7`` and ``e = 0.205630`` (scratch run, not kept) L1
  orbit 1 converges to within 2e-7 of the print, so the printed constants' rounding is a
  likely cause for that row; the others still fail. A four-segment multiple-shooting corrector
  without damping diverged for all five (scratch). These rows are held as strict expected
  failures.
* Table 3 (p. 14): monodromy eigenvalues over ``4 pi``. For the three reproduced orbits the
  largest eigenvalue does not match its row: L1 orbit 2 measured 7.346e4 against printed
  -5.7936e5, L1 orbit 3 measured 8.777e4 against 74,342, L2 orbit 3 measured 5.074e4 against
  43,176. The measured L1 values are within 1.2% of the printed values one row lower (74,342
  and 88,537), which suggests the L1 rows of Table 3 may be shifted by one relative to
  Table 2; this is not asserted. Held as strict expected failures. Structural checks that do
  hold: the monodromy is symplectic (determinant 1 to 2e-6) and its eigenvalues come in
  reciprocal pairs.

Controls: the 7:3 L1 halo corrected with the 5:2 period lands far from its print; an ``f0 = pi``
orbit's printed state propagated from ``f0 = 0`` (the sign of ``e`` reversed, Eq. 6) misses the
x-z plane by far more than from ``f0 = pi``.
"""

from __future__ import annotations

import math
from functools import cache

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.er3bp as er3bp

FloatArray = NDArray[np.float64]

MU_PRINTED = 1.660e-7  # p. 3
E_PRINTED = 0.2056  # p. 1 and figures

# Table 1 (p. 4): M, N, x0, z0, ydot0 (CRTBP, period 2 N pi / M)
TABLE1 = [
    ("5:2 L1", 5, 2, 0.997148, 0.004534, 0.005621),
    ("5:2 L2", 5, 2, 1.002834, 0.004548, -0.005592),
    ("7:3 L1", 7, 3, 0.996848, 0.004315, 0.005790),
    ("7:3 L2", 7, 3, 1.003134, 0.004335, -0.005765),
    ("9:4 L1", 9, 4, 0.996664, 0.004096, 0.005803),
    ("9:4 L2", 9, 4, 1.003318, 0.004125, -0.005784),
]

# Table 2 (p. 14): f0, x0, z0, ydot0, and the unit of the last printed digit of each.
TABLE2 = {
    "L1-1": (0.0, 0.997862, 0.00587395, 0.00198001, (1e-6, 1e-8, 1e-8)),
    "L1-2": (math.pi, 0.996917, 0.00389863, 0.00635667, (1e-6, 1e-8, 1e-8)),
    "L1-3": (math.pi, 0.99671, 0.00359628, 0.00724124, (1e-5, 1e-8, 1e-8)),
    "L1-4": (math.pi, 0.996201, 0.00296523, 0.00870505, (1e-6, 1e-8, 1e-8)),
    "L2-1": (0.0, 1.00211, 0.00588612, -0.00195444, (1e-5, 1e-8, 1e-8)),
    "L2-2": (math.pi, 1.00329, 0.00360261, -0.00722798, (1e-5, 1e-8, 1e-8)),
    "L2-3": (math.pi, 1.00308, 0.00390416, -0.00634580, (1e-5, 1e-8, 1e-8)),
    "L2-4": (math.pi, 1.00380, 0.00297788, -0.00868286, (1e-5, 1e-8, 1e-8)),
}

# Table 3 (p. 14): largest eigenvalue lambda1 of the monodromy over 4 pi.
TABLE3_LAMBDA1 = {
    "L1-1": 1.092e6,
    "L1-2": -5.7936e5,
    "L1-3": 74342.0,
    "L1-4": 88537.0,
    "L2-1": 1.042e6,
    "L2-2": -5.1771e5,
    "L2-3": 43176.0,
    "L2-4": 65906.0,
}

REPRODUCED = ("L1-2", "L1-3", "L2-3")
_NOT_REPRODUCED_REASON = (
    "#896: a damped single-shooting symmetric corrector started from the printed state at "
    "mu=1.660e-7, e=0.2056 does not converge to within one printed digit of the printed orbit"
)
_EIG_REASON = (
    "#896: largest monodromy eigenvalue of the reproduced orbit differs from its Table 3 row "
    "(L1 values match the next row down to 1.2 percent)"
)


def _propagate(f0: float, state: FloatArray, df: float, e: float) -> tuple[FloatArray, FloatArray]:
    y0 = np.concatenate([state, np.eye(6).ravel()])
    sol = solve_ivp(
        er3bp.er3bp_stm_eom,
        (f0, f0 + df),
        y0,
        args=(MU_PRINTED, e),
        method="DOP853",
        rtol=1e-13,
        atol=1e-14,
    )
    return sol.y[:6, -1], sol.y[6:, -1].reshape(6, 6)


def _correct(
    f0: float, guess: tuple[float, float, float], half: float, e: float, max_iter: int = 8
) -> tuple[FloatArray, float]:
    """Symmetric corrector: free ``(x0, z0, ydot0)``, target ``y = xdot = zdot = 0`` after
    ``half`` from ``f0``; Newton with residual-halving damping. Returns the state and the final
    residual norm."""
    g = np.array(guess, dtype=np.float64)

    def residual(v: FloatArray) -> tuple[FloatArray, FloatArray]:
        sf, phi = _propagate(f0, np.array([v[0], 0.0, v[1], 0.0, v[2], 0.0]), half, e)
        return sf[[1, 3, 5]], phi[np.ix_([1, 3, 5], [0, 2, 4])]

    r, jac = residual(g)
    for _ in range(max_iter):
        if float(np.linalg.norm(r)) < 1e-12:
            break
        step = np.linalg.solve(jac, -r)
        lam = 1.0
        while True:
            trial = g + lam * step
            r_trial, jac_trial = residual(trial)
            if float(np.linalg.norm(r_trial)) < float(np.linalg.norm(r)) or lam < 1e-4:
                break
            lam *= 0.5
        g, r, jac = trial, r_trial, jac_trial
    return g, float(np.linalg.norm(r))


@cache
def _table2_orbit(name: str) -> tuple[FloatArray, float]:
    f0, x0, z0, yd0, _ = TABLE2[name]
    return _correct(f0, (x0, z0, yd0), 2.0 * math.pi, E_PRINTED)


@cache
def _monodromy(name: str) -> FloatArray:
    f0 = TABLE2[name][0]
    g, _ = _table2_orbit(name)
    _, phi = _propagate(f0, np.array([g[0], 0.0, g[1], 0.0, g[2], 0.0]), 4.0 * math.pi, E_PRINTED)
    return phi


@pytest.mark.parametrize(("label", "m", "n", "x0", "z0", "yd0"), TABLE1)
def test_table1_circular_halo_has_the_resonant_period(
    label: str, m: int, n: int, x0: float, z0: float, yd0: float
) -> None:
    """With the half period fixed at ``N pi / M`` the corrected halo is the printed one."""
    g, res = _correct(0.0, (x0, z0, yd0), n * math.pi / m, 0.0)
    assert res < 1e-10, label
    np.testing.assert_allclose(g, [x0, z0, yd0], rtol=0.0, atol=1e-6, err_msg=label)


def test_table1_wrong_resonance_period_misses_print() -> None:
    """Control: the 7:3 L1 state corrected with the 5:2 half period lands far from the print."""
    _, _, _, x0, z0, yd0 = TABLE1[2]
    g, _ = _correct(0.0, (x0, z0, yd0), 2.0 * math.pi / 5.0, 0.0)
    assert float(np.max(np.abs(g - np.array([x0, z0, yd0])))) > 1e-4


@pytest.mark.parametrize(
    "name",
    [
        n
        if n in REPRODUCED
        else pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=_NOT_REPRODUCED_REASON))
        for n in TABLE2
    ],
)
def test_table2_orbit_closes_at_printed_state(name: str) -> None:
    """The corrector started from the printed state converges, and lands within one unit of
    the last printed digit of each component."""
    _f0, x0, z0, yd0, units = TABLE2[name]
    g, res = _table2_orbit(name)
    assert res < 1e-10, f"{name}: half-period residual {res:.2e}"
    for got, printed, unit in zip(g, (x0, z0, yd0), units, strict=True):
        assert abs(got - printed) <= unit, f"{name}: {got:.10f} vs printed {printed}"


@pytest.mark.parametrize("name", REPRODUCED)
def test_table2_monodromy_is_symplectic_with_reciprocal_pairs(name: str) -> None:
    phi = _monodromy(name)
    assert np.linalg.det(phi) == pytest.approx(1.0, abs=1e-5)
    lam = np.linalg.eigvals(phi)
    # every eigenvalue has a partner with product 1 (for a pair on the unit circle the partner
    # is its conjugate)
    for a in lam:
        assert min(abs(a * b - 1.0) for b in lam) < 1e-3


@pytest.mark.parametrize(
    "name",
    [pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=_EIG_REASON)) for n in REPRODUCED],
)
def test_table3_largest_eigenvalue(name: str) -> None:
    lam = np.linalg.eigvals(_monodromy(name))
    big = lam[np.argmax(np.abs(lam))]
    assert big.real == pytest.approx(TABLE3_LAMBDA1[name], rel=0.02)


def test_apoapsis_orbit_started_at_periapsis_misses_by_far_more() -> None:
    """Control (Eq. 6): L1 orbit 2 is an ``f0 = pi`` orbit. Its printed state, propagated for
    half a period from ``f0 = pi``, returns to the x-z plane up to the print rounding
    (measured residual 4.1e-4); started at ``f0 = 0`` (equivalent to reversing the sign of
    ``e``) it misses by 2.0e-2."""
    _, x0, z0, yd0, _ = TABLE2["L1-2"]
    s0 = np.array([x0, 0.0, z0, 0.0, yd0, 0.0])
    right, _ = _propagate(math.pi, s0, 2.0 * math.pi, E_PRINTED)
    wrong, _ = _propagate(0.0, s0, 2.0 * math.pi, E_PRINTED)
    r_right = float(np.linalg.norm(right[[1, 3, 5]]))
    r_wrong = float(np.linalg.norm(wrong[[1, 3, 5]]))
    assert r_wrong > 20.0 * r_right, (r_right, r_wrong)
