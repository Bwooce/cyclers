"""Published positive control for ``core/er3bp.py``: multi-revolution orbits of the Earth-Moon
elliptic restricted three-body problem.

Source: Neelakantan & Ramanan, "Two-impulse transfer to multi-revolution halo orbits in the
Earth-Moon elliptic restricted three body problem framework", Journal of Astrophysics and Astronomy
43:50 (2022), DOI 10.1007/s12036-022-09830-x, Table 8 ("Initial conditions and period of target MR
orbits"), with mu = 0.0122 and e = 0.0554 as used in the paper. Filed in the private paper corpus
as neelakantan-ramanan-2022-two-impulse-transfer-multi-revolution-halo-orbits-earth-moon-ER3BP-
jaa-43-50-doi-10.1007-s12036-022-09830-x.pdf.

The expected side is the paper's printed initial conditions and periods. The assertion is a
property the model must have if it is the same model: started with the primaries at periapsis,
each printed state returns to itself after the printed period (in true anomaly) and crosses the
x-z plane perpendicularly at the half period. The same states do not close in the circular
problem (control).

Not reproduced: the table's fifth row, "M4N2 Lyapunov" (x0 = 0.804504659626012, ydot0 =
0.31826866733409, period 4 pi), does not close in this model from either apse (closure 4.5 and
2.0). Its sibling M2N1, the same kind of orbit over one revolution of the primaries, closes to
1e-9. Whether the printed digits are too few for an orbit that unstable over two revolutions, or
the row has a slip, is not known; it is left out and recorded here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.er3bp as er3bp

_MU = 0.0122
_ECC = 0.0554

# (label, x0, z0, ydot0, period in true anomaly); y0 = xdot0 = zdot0 = 0.
_TABLE_8 = [
    ("M2N1 Lyapunov", 0.804125156956177, 0.0, 0.31182413982453, 2.0 * math.pi),
    ("M3N1 halo", 0.875404052867664, 0.201620653468241, 0.21551063837955, 2.0 * math.pi),
    ("M4N2 halo", 0.895820897947402, 0.194415672884192, 0.34702042423622, 4.0 * math.pi),
    ("M5N2 halo", 0.851666641652152, 0.183285539178136, 0.25828972225268, 4.0 * math.pi),
]


@pytest.mark.parametrize(("label", "x0", "z0", "vy0", "period"), _TABLE_8)
def test_printed_orbit_closes_from_periapsis(
    label: str, x0: float, z0: float, vy0: float, period: float
) -> None:
    """Measured closures: 1.1e-9, 4.4e-11, 1.4e-11 and 2.9e-7; half-period residuals 5e-12,
    7e-11, 9e-13 and 4e-10."""
    state0 = np.array([x0, 0.0, z0, 0.0, vy0, 0.0])
    sol = solve_ivp(
        er3bp.er3bp_eom,
        (0.0, period),
        state0,
        args=(_MU, _ECC),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert sol.success, label
    assert float(np.linalg.norm(sol.y[:, -1] - state0)) < 1e-6, label
    assert sol.sol is not None
    half = sol.sol(0.5 * period)
    assert math.hypot(half[1], half[3]) < 1e-8, label


@pytest.mark.parametrize(("label", "x0", "z0", "vy0", "period"), _TABLE_8)
def test_printed_orbit_does_not_close_without_eccentricity_or_from_apoapsis(
    label: str, x0: float, z0: float, vy0: float, period: float
) -> None:
    """Controls: the closure above is a property of the elliptic model at the right phase."""
    state0 = np.array([x0, 0.0, z0, 0.0, vy0, 0.0])
    for start, ecc in ((0.0, 0.0), (math.pi, _ECC)):
        sol = solve_ivp(
            er3bp.er3bp_eom,
            (start, start + period),
            state0,
            args=(_MU, ecc),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        assert float(np.linalg.norm(sol.y[:, -1] - state0)) > 0.1, (label, start, ecc)
