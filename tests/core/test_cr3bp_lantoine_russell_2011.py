"""#932(a): published CR3BP orbits of Lantoine & Russell (2011), Tables 1, 2 and 4.

Source: G. Lantoine and R. P. Russell, "Near Ballistic Halo-to-Halo Transfers between Planetary
Moons", J. Astronaut. Sci. 58(3):335-363 (2011), DOI 10.1007/BF03321174. Mass ratios from Table 1
(p. 350), the two halo orbits from Table 2 (p. 352), the five planar resonant orbits from Table 4
(p. 353). Every expected value below is a printed number; nothing is taken from this project's own
output. The printed column "Y0 (DU/TU)" is the velocity ydot0 (the caption states y0 = 0).

How the bounds were chosen
--------------------------
* Jacobi constant: printed to 4 decimals, so the bound is the rounding half-unit, 5e-5.
* Closure after one printed period (norm of the 6-state difference, nondimensional): the initial
  conditions are printed to 11 to 13 digits but the mass ratios only to 5 significant digits
  (7.8037e-5, 2.5280e-5), i.e. an absolute uncertainty of 5e-10. A check with mu shifted by
  +-5e-10 moves the closure by about 4e-4 for the unstable halos and by 2e-6 to 4e-5 for the
  resonant orbits, so the bounds are set a little above that mass-ratio rounding effect: 1e-3 for
  the halos (unstable, 3 TU) and 2e-4 for the resonant orbits. They are not set from our own
  residual on the printed mu, which is smaller. The control (the other moon's mass ratio) must
  miss by more than 1e-2, over 40 times the largest bound.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core import cr3bp

MU_GANYMEDE = 7.8037e-5  # Table 1
MU_EUROPA = 2.5280e-5  # Table 1

# name, mu, state (x, y, z, vx, vy, vz), period (TU), printed C, closure bound
ORBITS = [
    (
        "halo1",
        MU_GANYMEDE,
        [0.9768297703815, 0, 0.0067575508137, 0, -0.0338705588350, 0],
        3.01368039319,
        3.0066,
        1e-3,
    ),
    (
        "halo2",
        MU_EUROPA,
        [1.0118043920085, 0, 0.0087547927135, 0, 0.0357066388227, 0],
        3.06456025428,
        3.0024,
        1e-3,
    ),
    (
        "res_3_4",
        MU_GANYMEDE,
        [0.9639250025000, 0, 0, 0, -0.037537693295765, 0],
        19.1527202833,
        3.0066,
        2e-4,
    ),
    (
        "res_9_7",
        MU_EUROPA,
        [1.022912512093, 0, 0, 0, 0.035866768601460, 0],
        56.8415853699,
        3.0024,
        2e-4,
    ),
    (
        "res_4_3",
        MU_EUROPA,
        [1.025860244947, 0, 0, 0, 0.038403138969070, 0],
        25.3393083838,
        3.0024,
        2e-4,
    ),
    (
        "res_11_8",
        MU_EUROPA,
        [1.028261885259, 0, 0, 0, 0.040854642490345, 0],
        69.268896450,
        3.0024,
        2e-4,
    ),
    (
        "res_7_5",
        MU_EUROPA,
        [1.029619747357, 0, 0, 0, 0.042564223924969, 0],
        44.117093502,
        3.0024,
        2e-4,
    ),
]
IDS = [o[0] for o in ORBITS]


def _system(mu: float) -> cr3bp.CR3BPSystem:
    # Only mu is used by propagate(); the dimensional fields are placeholders.
    return cr3bp.CR3BPSystem(mu=mu, primary="Jupiter", secondary="moon", l_km=1.0, t_s=1.0)


def _closure(mu: float, state: list[float], period: float) -> float:
    s0 = np.array(state, dtype=float)
    sf = cr3bp.propagate(_system(mu), s0, period).state_f
    return float(np.linalg.norm(sf - s0))


@pytest.mark.parametrize(("name", "mu", "state", "period", "c_printed", "bound"), ORBITS, ids=IDS)
def test_jacobi_constant_matches_printed_digits(
    name: str, mu: float, state: list[float], period: float, c_printed: float, bound: float
) -> None:
    c = cr3bp.jacobi_constant(np.array(state, dtype=float), mu)
    assert abs(c - c_printed) <= 5e-5, (name, c, c_printed)


@pytest.mark.parametrize(("name", "mu", "state", "period", "c_printed", "bound"), ORBITS, ids=IDS)
def test_orbit_closes_after_printed_period(
    name: str, mu: float, state: list[float], period: float, c_printed: float, bound: float
) -> None:
    assert _closure(mu, state, period) <= bound, name


def test_control_wrong_moon_mass_ratio_fails_to_close() -> None:
    # Ganymede's halo with Europa's mass ratio is a different orbit: it must not close.
    name, _mu, state, period, _c, bound = ORBITS[0]
    assert bound < 1e-2, name
    assert _closure(MU_EUROPA, state, period) > 1e-2, name
