"""#896: circular-limit monodromy eigenvalues of the M5N2 halo against Peng & Xu (2015).

Source: H. Peng and S. Xu, "Stability of two groups of multi-revolution elliptic halo orbits in
the elliptic restricted three-body problem", Celestial Mechanics and Dynamical Astronomy 123
(2015) 279-303, doi:10.1007/s10569-015-9635-2; filed in the private paper corpus as
peng-xu-2015-stability-two-groups-multirev-elliptic-halo-orbits-ERTBP-cmda-123-279-doi-10.1007-s10569-015-9635-2.pdf.

The paper prints no initial conditions, only monodromy eigenvalues of L1 north-halo ME-Halo
orbits with M = 5, N = 2 over the elliptic-problem period T_E = 4 pi (Eq. 10). At ``e = 0`` the
orbit is the circular-problem halo of period ``T_C = 2 N pi / M = 4 pi / 5``, extracted from the
L1 halo family by bisection in the z amplitude (p. 285), traversed five times; p. 294 says the
``e = 0`` rows are "the halo orbits revolving M = 5 revolutions in the corresponding CRTBP
systems, so their true value of lambda2 (or 1/lambda2) should be one", and attributes the printed
1.0002 and 1.0000 + 0.0001i to integration error. So the ``e = 0`` rows are the eigenvalues of the
fifth power of the single-period monodromy of that halo, which ``core/cr3bp.py`` must reproduce.

Printed ``e = 0`` rows used (digits reread from the PDF):

* Table 1 (p. 294), periapsis group, mu = 0.009 (Fig. 10): 1.3112e+06, 7.6258e-07, 1.0002,
  0.9998, 0.9587 +- 0.2843i.
* Table 4 (p. 299), apoapsis group, mu = 0.009: 1.3112e+06, 7.6272e-07, 1.0000 +- 0.0001i,
  0.9587 +- 0.2843i.
* Table 2 (p. 295), periapsis group, mu = 0.015 (Fig. 11): 1.5966e+06, 6.2662e-07,
  1.0000 +- 0.0001i, 0.9021 +- 0.4316i.
* Table 5 (p. 300), apoapsis group, mu = 0.015: 1.5966e+06, 6.2636e-07, 1.0002, 0.9998,
  0.9021 +- 0.4316i.
* Table 3 (p. 299), apoapsis group, mu = 0.004: 9.3110e+05, 1.0740e-06, 1.0004, 0.9996,
  0.9021 +- 0.4316i.

Method here: a symmetric single-shooting corrector with the half period fixed at ``2 pi / 5``
(free ``x0, z0, ydot0`` on the x-z plane, target ``y = xdot = zdot = 0``), integrated with
``core.cr3bp.cr3bp_stm_eom``. The seed is the M5N2 halo of Neelakantan & Ramanan 2022, Table 8
(mu = 0.0122, e = 0.0554), used only as a starting guess; it is corrected at ``e = 0`` and
continued in mu to 0.009, 0.015 and 0.004. The single-period monodromy ``M1`` is diagonalised
and its eigenvalues raised to the fifth power (diagonalising ``M1**5`` directly loses the small
eigenvalue to conditioning near 1e12).

Measured (fifth powers of the single-period eigenvalues):

* mu = 0.009: x0 = 0.870119, z0 = 0.162204, ydot0 = 0.233997; 1.311204e6, 7.62658e-7,
  1 +- 3.6e-6 i, 0.958732 +- 0.284311i.
* mu = 0.015: x0 = 0.838551, z0 = 0.190396, ydot0 = 0.282355; 1.596597e6, 6.26332e-7,
  1 +- 3.4e-6 i, 0.902069 +- 0.431592i.
* mu = 0.004: x0 = 0.906241, z0 = 0.125478, ydot0 = 0.174558; 9.31098e5, 1.07400e-6,
  1 +- 4.9e-6 i, 0.989628 +- 0.143652i.

These initial states are consistent with the ``e = 0`` (green) ends of the characteristic curves
of Fig. 6 (p. 291), whose axes span x0 0.8 to 0.95, z0 0.1 to 0.25, y0' 0.1 to 0.3.

Findings, held as strict expected failures:

* The printed small eigenvalue 1/lambda1 does not match 1/lambda1 of the measured orbit (equal to
  the reciprocal of the measured lambda1, which itself matches every print to within 0.2 units
  of the last printed digit) to one unit of its last printed digit in Tables 1, 2, 4 and 5:
  printed 7.6258 and 7.6272 against 7.6266 (e-07), printed 6.2662 and 6.2636 against 6.2633
  (e-07), misses of 7.8, 6.2, 28.8 and 2.8 units of the last printed digit. The two
  printings of the same e = 0 orbit (Tables 1 and 4; 2 and 5) differ from each other by 1.4 and
  2.6 units of the fourth digit, so the paper's small eigenvalue carries about 1e-4 to 5e-4
  relative error. Table 3's 1.0740e-06 matches.
* Table 3's lambda3 column (0.9021 +- 0.4316i at e = 0) is identical in all twelve rows to the
  lambda3 column of Table 2 (mu = 0.015). The mu = 0.004 halo measured here has the unit pair
  0.989628 +- 0.143652i, while its lambda1 and 1/lambda1 match Table 3; the printed column is a
  copy of Table 2's.

Controls, each of which must miss the print by far more than the tolerance: the single-period
monodromy (largest eigenvalue 16.73), the halo of the other mass ratio, the halo of period
``6 pi / 7`` (M7N3) raised to the fifth power, and the other member of the L1 halo family that
has the same period ``4 pi / 5`` (beyond the minimum-period orbit, z0 = 0.351 at mu = 0.009,
whose fifth-power monodromy has a complex quartet of modulus 4.0e6 and no large real eigenvalue).

Not covered: the elliptic rows (e > 0) and the Earth-Moon numbers of Section 3.3 (mu = 0.0122,
e = 0.0554, lambda1 about 1.5427e6) need an elliptic multi-segment corrector the project lacks;
blocked on #912.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_stm_eom

FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]

T_C = 4.0 * math.pi / 5.0  # Eq. 10 with M = 5, N = 2

# Neelakantan & Ramanan 2022, Table 8, M5N2 halo (mu 0.0122, e 0.0554): corrector seed only.
_SEED_NR_M5N2 = (0.851666641652152, 0.183285539178136, 0.25828972225268)
_MU_SEED = 0.0122

# Corrector seeds for the other same-period member of the family (from a scratch fixed-z0
# continuation through the minimum-period orbit); seeds, not expected values.
_SEED_OTHER_MEMBER = {
    0.009: (0.9197660734568641, 0.35106144371003845, 0.08223632798047165),
    0.015: (0.8854176798669642, 0.4107106094751497, 0.11535532783625413),
}


@dataclass(frozen=True)
class Row:
    table: str
    mu: float
    lam1: float
    lam1_unit: float
    inv_lam1: float
    inv_lam1_unit: float
    unit_pair: complex  # printed lambda3 with positive imaginary part
    trivial_dev: float  # largest printed departure of the trivial pair from 1


ROWS = {
    "T1": Row("Table 1 p294", 0.009, 1.3112e6, 1e2, 7.6258e-07, 1e-11, 0.9587 + 0.2843j, 2e-4),
    "T4": Row("Table 4 p299", 0.009, 1.3112e6, 1e2, 7.6272e-07, 1e-11, 0.9587 + 0.2843j, 1e-4),
    "T2": Row("Table 2 p295", 0.015, 1.5966e6, 1e2, 6.2662e-07, 1e-11, 0.9021 + 0.4316j, 1e-4),
    "T5": Row("Table 5 p300", 0.015, 1.5966e6, 1e2, 6.2636e-07, 1e-11, 0.9021 + 0.4316j, 2e-4),
    "T3": Row("Table 3 p299", 0.004, 9.3110e5, 1e1, 1.0740e-06, 1e-10, 0.9021 + 0.4316j, 4e-4),
}

_INV_REASON = (
    "#896: printed 1/lambda1 (7.6258e-07 T1, 7.6272e-07 T4, 6.2662e-07 T2, 6.2636e-07 T5) is "
    "not within one unit of its last digit of 1/lambda1 of the reproduced halo (7.6266e-07, "
    "6.2633e-07), whose lambda1 matches the print; the paper's small eigenvalue is inaccurate"
)
_T3_LAMBDA3_REASON = (
    "#896: Table 3 lambda3 (0.9021 + 0.4316i) is Table 2's column copied; the mu = 0.004 halo "
    "that matches Table 3's lambda1 has unit pair 0.989628 + 0.143652i"
)


def _propagate(state: FloatArray, t: float, mu: float) -> tuple[FloatArray, FloatArray]:
    y0 = np.concatenate([state, np.eye(6).ravel()])
    sol = solve_ivp(
        cr3bp_stm_eom, (0.0, t), y0, args=(mu,), method="DOP853", rtol=1e-13, atol=1e-14
    )
    return sol.y[:6, -1], sol.y[6:, -1].reshape(6, 6)


def _state(v: FloatArray) -> FloatArray:
    return np.array([v[0], 0.0, v[1], 0.0, v[2], 0.0])


def _correct(guess: FloatArray, half: float, mu: float) -> tuple[FloatArray, float]:
    """Newton on ``(x0, z0, ydot0)`` for ``y = xdot = zdot = 0`` at fixed time ``half``."""
    g = np.array(guess, dtype=np.float64)
    res = math.inf
    for _ in range(25):
        sf, phi = _propagate(_state(g), half, mu)
        r = sf[[1, 3, 5]]
        res = float(np.linalg.norm(r))
        if res < 1e-12:
            break
        g = g - np.linalg.solve(phi[np.ix_([1, 3, 5], [0, 2, 4])], r)
    return g, res


def _continue(g: FloatArray, mu_from: float, mu_to: float, period: float) -> FloatArray:
    for mu in np.linspace(mu_from, mu_to, 9)[1:]:
        g, res = _correct(g, period / 2.0, float(mu))
        assert res < 1e-11, (mu, res)
    return g


@dataclass(frozen=True)
class Halo:
    ic: FloatArray  # x0, z0, ydot0
    single: ComplexArray  # eigenvalues of the single-period monodromy

    @property
    def fifth(self) -> ComplexArray:
        return self.single**5


def _halo(ic: FloatArray, mu: float, period: float = T_C) -> Halo:
    _, m1 = _propagate(_state(ic), period, mu)
    return Halo(ic, np.linalg.eigvals(m1).astype(np.complex128))


@pytest.fixture(scope="module")
def halos() -> dict[float, Halo]:
    g, res = _correct(np.array(_SEED_NR_M5N2), T_C / 2.0, _MU_SEED)
    assert res < 1e-11
    g09 = _continue(g, _MU_SEED, 0.009, T_C)
    g15 = _continue(g, _MU_SEED, 0.015, T_C)
    g04 = _continue(g09, 0.009, 0.004, T_C)
    return {0.009: _halo(g09, 0.009), 0.015: _halo(g15, 0.015), 0.004: _halo(g04, 0.004)}


def _largest(lam: ComplexArray) -> complex:
    return complex(lam[np.argmax(np.abs(lam))])


def _smallest(lam: ComplexArray) -> complex:
    return complex(lam[np.argmin(np.abs(lam))])


def _trivial_pair(lam: ComplexArray) -> ComplexArray:
    return lam[np.argsort(np.abs(lam - 1.0))[:2]]


def _unit_pair_upper(lam: ComplexArray) -> complex:
    on_circle = [z for z in lam if abs(abs(z) - 1.0) < 1e-3 and abs(z.imag) > 1e-2]
    assert len(on_circle) == 2, lam
    return complex(max(on_circle, key=lambda z: z.imag))


@pytest.mark.parametrize("key", list(ROWS))
def test_largest_eigenvalue_matches_print(halos: dict[float, Halo], key: str) -> None:
    """lambda1 of the fifth power of the single-period monodromy is the printed value to one
    unit of the last printed digit (measured 1.311204e6, 1.596597e6, 9.31098e5)."""
    row = ROWS[key]
    big = _largest(halos[row.mu].fifth)
    assert abs(big.imag) < 1e-6 * abs(big)
    assert abs(big.real - row.lam1) <= row.lam1_unit, (row.table, big)


@pytest.mark.parametrize("key", ["T1", "T4", "T2", "T5"])
def test_unit_circle_pair_matches_print(halos: dict[float, Halo], key: str) -> None:
    """The unit-circle pair matches the printed four decimals to one unit (measured
    0.958732 + 0.284311i and 0.902069 + 0.431592i)."""
    row = ROWS[key]
    z = _unit_pair_upper(halos[row.mu].fifth)
    assert abs(z.real - row.unit_pair.real) <= 1e-4, (row.table, z)
    assert abs(z.imag - row.unit_pair.imag) <= 1e-4, (row.table, z)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=_T3_LAMBDA3_REASON)
def test_table3_unit_circle_pair(halos: dict[float, Halo]) -> None:
    row = ROWS["T3"]
    z = _unit_pair_upper(halos[row.mu].fifth)
    assert abs(z.real - row.unit_pair.real) <= 1e-4, (row.table, z)
    assert abs(z.imag - row.unit_pair.imag) <= 1e-4, (row.table, z)


@pytest.mark.parametrize("key", list(ROWS))
def test_trivial_pair_is_one_within_printed_error(halos: dict[float, Halo], key: str) -> None:
    """p. 294: the true value is one; the printed departure (up to 4e-4) is the authors' error.
    Measured departures 3.4e-6 to 4.9e-6."""
    row = ROWS[key]
    for z in _trivial_pair(halos[row.mu].fifth):
        assert abs(z - 1.0) <= row.trivial_dev, (row.table, z)


@pytest.mark.parametrize(
    "key",
    [
        k
        if k == "T3"
        else pytest.param(
            k, marks=pytest.mark.xfail(strict=True, raises=AssertionError, reason=_INV_REASON)
        )
        for k in ROWS
    ],
)
def test_smallest_eigenvalue_matches_print(halos: dict[float, Halo], key: str) -> None:
    """1/lambda1 to one unit of the last printed digit (measured 7.62658e-7, 6.26332e-7,
    1.07400e-6)."""
    row = ROWS[key]
    small = _smallest(halos[row.mu].fifth)
    assert abs(small.real - row.inv_lam1) <= row.inv_lam1_unit, (row.table, small)


def test_symplectic_reciprocal_of_largest(halos: dict[float, Halo]) -> None:
    """Model property used above: the smallest eigenvalue is the reciprocal of the largest."""
    for h in halos.values():
        assert abs(_largest(h.fifth) * _smallest(h.fifth) - 1.0) < 1e-6


# ---- controls -------------------------------------------------------------------------------


@pytest.mark.parametrize("mu", [0.009, 0.015])
def test_control_single_period_monodromy_misses(halos: dict[float, Halo], mu: float) -> None:
    """The single-period monodromy (largest eigenvalue 16.73 and 17.40) is not the printed one."""
    big = _largest(halos[mu].single).real
    assert 10.0 < big < 30.0
    printed = ROWS["T1"].lam1 if mu == 0.009 else ROWS["T2"].lam1
    assert abs(big - printed) > 1e6


def test_control_wrong_mass_ratio_misses(halos: dict[float, Halo]) -> None:
    """Each mass ratio's halo misses the other table's lambda1 by more than 2.8e5, against a
    tolerance of 100."""
    for mu, other in ((0.009, "T2"), (0.015, "T1")):
        big = _largest(halos[mu].fifth).real
        assert abs(big - ROWS[other].lam1) > 1e3 * ROWS[other].lam1_unit
        z = _unit_pair_upper(halos[mu].fifth)
        assert abs(z - ROWS[other].unit_pair) > 0.1


def test_control_m7n3_period_misses(halos: dict[float, Halo]) -> None:
    """The mu = 0.009 halo of period 6 pi / 7 (M7N3), fifth power: largest eigenvalue 6.9e8."""
    g = halos[0.009].ic
    for period in np.linspace(T_C, 6.0 * math.pi / 7.0, 6)[1:]:
        g, res = _correct(g, float(period) / 2.0, 0.009)
        assert res < 1e-11
    h = _halo(g, 0.009, 6.0 * math.pi / 7.0)
    big = _largest(h.fifth).real
    assert abs(big - ROWS["T1"].lam1) > 1e3 * ROWS["T1"].lam1_unit, big


@pytest.mark.parametrize("mu", [0.009, 0.015])
def test_control_other_family_member_of_same_period_misses(
    halos: dict[float, Halo], mu: float
) -> None:
    """The L1 halo family has a second member of period 4 pi / 5 past the minimum-period orbit
    (z0 = 0.351 and 0.411, outside Fig. 6's z0 range). Its fifth-power monodromy has a complex
    quartet of modulus 4.0e6 (mu = 0.009) and 3.3e6 (mu = 0.015) and no real eigenvalue near the
    printed lambda1, so the print selects the low-amplitude member."""
    g, res = _correct(np.array(_SEED_OTHER_MEMBER[mu]), T_C / 2.0, mu)
    assert res < 1e-11
    assert abs(g[1] - halos[mu].ic[1]) > 0.1
    lam = _halo(g, mu).fifth
    printed = ROWS["T1"].lam1 if mu == 0.009 else ROWS["T2"].lam1
    big = _largest(lam)
    assert abs(big.imag) > 1e6, big
    real_large = [z for z in lam if abs(z.imag) < 1e-6 * abs(z) and abs(z) > 10.0]
    assert not real_large, lam
    assert abs(abs(big) - printed) > 1e3 * ROWS["T1"].lam1_unit
