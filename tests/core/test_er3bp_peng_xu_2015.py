"""Published positive controls for ``core/er3bp.py``: the Earth-Moon L1 M5N2 multi-revolution halo.

Sources (expected values are printed numbers only, never values this code computed):

* Peng & Xu, "Transfer to a multi-revolution elliptic halo orbit in Earth-Moon elliptic restricted
  three-body problem using stable manifold", Adv. Space Res. 55:1015-1027 (2015), DOI
  10.1016/j.asr.2014.11.013, Table 2 (p1020): mu = 0.0122, e = 0.0554, initial state
  [0.8516666, 0, 0.1832855, 0, 0.2582897, 0] at f0 = 0 (periapsis), period 4 pi (M = 5, N = 2).
  Filed in the private paper corpus as peng-xu-2015-transfer-multirev-elliptic-halo-earth-moon-
  ER3BP-stable-manifold-asr-55-1015-doi-10.1016-j.asr.2014.11.013.pdf.
* Neelakantan & Ramanan, J. Astrophys. Astron. 43:50 (2022), DOI 10.1007/s12036-022-09830-x,
  Table 8, row "M5N2 halo" (15 digits).
* Peng & Xu, Celest. Mech. Dyn. Astron. 123:279 (2015), DOI 10.1007/s10569-015-9635-2, p300: for
  this orbit the monodromy has real eigenvalues 1.5427e6 and 1.0086 (Earth-Moon, e = 0.0554).

The seven-digit Table 2 state is a rounding of the 15-digit row; with a multiplier of 1.5e6 that
rounding alone makes the printed state miss closure by about 4e-2, so it is tested by
re-correction (three unknowns, three half-period perpendicular-crossing residuals), not by closure.
Notes: ``docs/notes/2026-10-04-digest-peng-xu-2015-asr-transfer-multirev-elliptic-halo.md``.
"""

from __future__ import annotations

import math
from itertools import pairwise

import numpy as np
from scipy.optimize import fsolve

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp
from cyclerfinder.search.er3bp_periodic import monodromy_eigenstructure

_MU = 0.0122
_ECC = 0.0554
_T_E = 4.0 * math.pi  # N = 2 revolutions of the primaries (printed)
_SYS = ER3BPSystem(mu=_MU, e=_ECC, primary_name="Earth", secondary_name="Moon")

# Peng & Xu ASR Table 2 (7 digits) and Neelakantan & Ramanan Table 8 (15 digits).
_PX = np.array([0.8516666, 0.0, 0.1832855, 0.0, 0.2582897, 0.0])
_NR = np.array(
    [0.851666641652152, 0.0, 0.183285539178136, 0.0, 0.25828972225268, 0.0],
)
_TOL = 1e-13


def _prop(state: np.ndarray, f0: float, span: float, sys: ER3BPSystem = _SYS) -> np.ndarray:
    _f, hist, _stm = propagate_er3bp(state, (f0, f0 + span), sys, rtol=_TOL, atol=_TOL)
    return np.asarray(hist)


def _monodromy(state: np.ndarray, nseg: int = 4) -> np.ndarray:
    edges = np.linspace(0.0, _T_E, nseg + 1)
    s = state.copy()
    mono = np.eye(6)
    for a, b in pairwise(edges):
        _f, hist, stm = propagate_er3bp(s, (a, b), _SYS, rtol=_TOL, atol=_TOL, with_stm=True)
        s = hist[:, -1]
        mono = stm @ mono
    return mono


def _half_period_residual(u: np.ndarray) -> list[float]:
    """Perpendicular crossing of the x-z plane at f = T_E / 2 = 2 pi: y, x', z' vanish."""
    s0 = np.array([u[0], 0.0, u[1], 0.0, u[2], 0.0])
    h = _prop(s0, 0.0, 0.5 * _T_E)
    return [float(h[1, -1]), float(h[3, -1]), float(h[5, -1])]


def test_printed_seven_digit_state_recorrects_to_the_published_orbit() -> None:
    """Re-correction moves the Table 2 state by under 1e-7 and lands on the 15-digit row."""
    u0 = np.array([_PX[0], _PX[2], _PX[4]])
    u = fsolve(_half_period_residual, u0, xtol=1e-14)
    assert float(np.max(np.abs(u - u0))) < 1e-7  # measured 4.2e-8
    nr = np.array([_NR[0], _NR[2], _NR[4]])
    assert float(np.max(np.abs(u - nr))) < 1e-7  # measured 2e-13
    assert float(np.max(np.abs(_half_period_residual(u)))) < 1e-9


def test_printed_seven_digit_state_is_within_rounding_of_a_perpendicular_crossing() -> None:
    """Not a closure test: the printed digits are amplified by the 1.5e6 multiplier."""
    res = _half_period_residual(np.array([_PX[0], _PX[2], _PX[4]]))
    assert max(abs(r) for r in res) < 5e-4  # measured 1.8e-4


def test_fifteen_digit_row_closes_from_periapsis() -> None:
    final = _prop(_NR, 0.0, _T_E)[:, -1]
    assert float(np.linalg.norm(final - _NR)) < 1e-6  # measured 1.6e-8


def test_apoapsis_start_does_not_close() -> None:
    """Control: the periapsis-group orbit is not periodic when started at f0 = pi."""
    final = _prop(_NR, math.pi, _T_E)[:, -1]
    assert float(np.linalg.norm(final - _NR)) > 0.1  # measured 1.18


def test_monodromy_real_eigenvalues_match_printed_values() -> None:
    eigs = np.linalg.eigvals(_monodromy(_NR))
    real_gt1 = sorted(
        (float(z.real) for z in eigs if abs(z.imag) < 1e-9 and z.real > 1.0), reverse=True
    )
    assert len(real_gt1) == 2
    assert abs(real_gt1[0] - 1.5427e6) / 1.5427e6 < 1e-3  # computed 1.54273e6
    assert abs(real_gt1[1] - 1.0086) < 5e-4  # computed 1.00858


def test_monodromy_eigenstructure_positive_control() -> None:
    """The existing helper returns the printed largest eigenvalue and finds the unit-circle pair."""
    r, w = monodromy_eigenstructure(_monodromy(_NR))
    assert abs(r - 1.5427e6) / 1.5427e6 < 1e-3
    assert 0.0 < w < math.pi


def test_apoapsis_equals_periapsis_with_negative_eccentricity() -> None:
    """The equations depend on f only through e cos f, so (f0 = pi, e) is (f0 = 0, -e)."""
    neg = ER3BPSystem(mu=_MU, e=-_ECC, primary_name="Earth", secondary_name="Moon")
    a = _prop(_NR, math.pi, 2.0 * math.pi)[:, -1]
    b = _prop(_NR, 0.0, 2.0 * math.pi, neg)[:, -1]
    assert float(np.linalg.norm(a - b)) < 1e-12  # measured 6e-15
