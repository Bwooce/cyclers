"""#896: the Uranus-Oberon L1 and L2 Jacobi constants printed by Kumar & Anderson (2024).

Source: B. Kumar and R. L. Anderson, "A survey of Oberon mean motion resonant unstable orbit
properties and connections for Uranian tours", AAS/AIAA paper AAS 24-288 (2024). Filed in the
private paper corpus as anderson-kumar-2024-oberon-mmr-unstable-orbit-survey-aas-24-288.pdf. Digest:
``docs/notes/2026-07-27-728-anderson-kumar-2024-oberon-mmr-survey-digest.md``.

What the paper prints (page 7): "for the value of mu = 3.54326 x 10^-5 we use for the system,
C = 3.00454 at L1 and C = 3.00450 at L2". Page 2, eq. (2): H0 = (px^2 + py^2)/2 + px y - py x
- (1 - mu)/r1 - mu/r2 with Uranus at -mu, and "C = -2H is referred to as the
Jacobi constant". At rest in the rotating frame that is C = x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2,
the convention of ``core.cr3bp.jacobi_constant`` (no added constant mu(1 - mu)).

Measured (2026-10-04) with ``core.cr3bp`` at the printed mu: C(L1) = 3.0045498, C(L2) = 3.0045025.
Both printed values are the computed ones cut to five decimals; rounded, L1 would be 3.00455, so
the paper truncates (or computed L1 slightly differently; it does not say). The test asserts what
holds under either reading: each printed value is at most one unit of its last decimal below the
computed one, and L1 lies above L2. Control: with the constant mu(1 - mu) added (the convention of
some other papers, e.g. Font, Nunes & Simo) the values are 3.0045852 and 3.0045379, and the L1 value
is 4.5 units of the last decimal above the print.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from cyclerfinder.core.cr3bp import jacobi_constant

_MU = 3.54326e-5  # page 2 and page 7
_PRINTED_C_L1 = 3.00454  # page 7
_PRINTED_C_L2 = 3.00450  # page 7
_LAST_DECIMAL = 1e-5


def _collinear_jacobi(lo: float, hi: float, mu: float) -> float:
    def d_omega(x: float) -> float:
        return (
            x
            - (1.0 - mu) * (x + mu) / abs(x + mu) ** 3
            - mu * (x - 1.0 + mu) / abs(x - 1.0 + mu) ** 3
        )

    x = float(brentq(d_omega, lo, hi, xtol=1e-15))
    return jacobi_constant(np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0]), mu)


def _c_l1(mu: float) -> float:
    return _collinear_jacobi(0.9, 1.0 - mu - 1e-6, mu)


def _c_l2(mu: float) -> float:
    return _collinear_jacobi(1.0 - mu + 1e-6, 1.1, mu)


def test_printed_libration_jacobi_constants_within_one_last_decimal() -> None:
    c_l1, c_l2 = _c_l1(_MU), _c_l2(_MU)
    assert 0.0 <= c_l1 - _PRINTED_C_L1 < _LAST_DECIMAL
    assert 0.0 <= c_l2 - _PRINTED_C_L2 < _LAST_DECIMAL
    assert c_l1 > c_l2


def test_jacobi_with_added_constant_misses_the_print() -> None:
    """Control: the convention C' = C + mu(1 - mu) puts L1 outside one last decimal."""
    shift = _MU * (1.0 - _MU)
    assert _c_l1(_MU) + shift - _PRINTED_C_L1 > 4.0 * _LAST_DECIMAL
