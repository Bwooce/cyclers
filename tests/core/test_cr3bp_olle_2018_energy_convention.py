"""#933/#935 convention golden: Olle, Rodriguez & Soler 2018 printed energies against core.cr3bp.

Source: M. Olle, O. Rodriguez and J. Soler, "Ejection-collision orbits in the RTBP", Commun.
Nonlinear Sci. Numer. Simulat. 55:298-315 (2018), DOI 10.1016/j.cnsns.2017.07.013 (manuscript
pp. 15, 20, 25; digest docs/notes/2026-10-05-digest-rodriguez-del-rio-2021-thesis-part-a.md
section 8). The paper works in H = -C/2 with the mu(1 - mu) term kept in the Hamiltonian (eq. 6),
so H = -(C_core + mu (1 - mu)) / 2 where C_core is core.cr3bp.jacobi_constant. Printed values:

* H_L1(0.5) = -2.125
* H_L2(0.5) = -1.853398112043077
* H_L1(0.1) = -1.843476614939948

The paper's frame (P1 big at (-mu, 0), P2 at (1 - mu, 0)) is core.cr3bp's, so no transformation
applies. The collinear point is located here by an independent Newton iteration on the
pseudo-potential gradient (not by core code) so that the L point and the energy are both
independent of the function under test except for jacobi_constant itself.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import jacobi_constant


def _collinear(mu: float, which: str) -> float:
    """Newton on Omega_x(x, 0) = 0 (written out here, independent of src)."""

    def f(x: float) -> float:
        r1, r2 = abs(x + mu), abs(x - 1.0 + mu)
        return x - (1 - mu) * (x + mu) / r1**3 - mu * (x - 1 + mu) / r2**3

    lo, hi = {
        "L1": (1 - mu - 0.5, 1 - mu - 1e-9),
        "L2": (1 - mu + 1e-9, 1 - mu + 1.0),
    }[which]
    for _ in range(200):  # bisection to machine precision
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


@pytest.mark.parametrize(
    ("mu", "point", "h_printed"),
    [
        (0.5, "L1", -2.125),
        (0.5, "L2", -1.853398112043077),
        (0.1, "L1", -1.843476614939948),
    ],
)
def test_olle_2018_printed_h_equals_core_jacobi_convention(
    mu: float, point: str, h_printed: float
) -> None:
    x = _collinear(mu, point)
    state = np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0])
    h_core = -(jacobi_constant(state, mu) + mu * (1.0 - mu)) / 2.0
    assert h_core == pytest.approx(h_printed, abs=1e-12)
