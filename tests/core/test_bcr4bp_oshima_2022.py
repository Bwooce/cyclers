"""#896: Oshima (2022) Tables 2 to 5 as published controls for the bicircular model.

K. Oshima (2022), "Multiple families of synodic resonant periodic orbits in the bicircular
restricted four-body problem", Advances in Space Research 70:1325-1335,
doi:10.1016/j.asr.2022.06.009. Filed in the private paper corpus as
oshima-2022-multiple-families-synodic-resonant-periodic-orbits-bicircular-restricted-four-body-problem-asr-70-1325-doi-10.1016-j.asr.2022.06.009.pdf

Tables 2 to 5 (p1332) print, for four families of 1:1 synodic resonant three-dimensional
retrograde orbits about the Earth (z0-S0, z0-Spi, vz0-S0, vz0-Spi), the state
``(x, z, vx, vy, vz)`` and the Sun angle ``theta_S`` at the four crossings #1 to #4 of
``y = 0``, "the modulo operation express[ing] theta_S between 0 and 2 pi". Constants are Table 1
(p1326). The Sun is at ``a_S (cos theta_S, sin theta_S, 0)`` with ``theta_S = theta_S0 + omega_S
t`` (eqs. 4 and 5) and Table 1 prints ``omega_S = -0.925195985``: the angle falls, the Sun goes
clockwise, Earth at ``-mu``. This module has the Sun at ``theta_sun0 - omega_sun t`` with
``omega_sun > 0`` and the Earth at ``-mu``, so a printed ``theta_S`` is this module's Sun angle at
that instant without any change, and the time from crossing k to crossing k + 1 is the fall of
the Sun angle, ``((theta_k - theta_{k+1}) mod 2 pi) / |omega_S|``. The period is one Sun period,
``2 pi / |omega_S|`` (eq. 10 with N = 1).

Table 4 row 1's one-period closure is already in ``test_bcr4bp_sun_sense.py`` and is not repeated.

Measured (rtol = atol = 1e-13): the other fifteen rows close after one Sun period to 9.6e-10 to
2.7e-8; each printed crossing propagated to the next printed Sun angle lands on the next printed
row to 1.5e-9 to 6.7e-9, and the ``y = 0`` crossings fall at the times the Sun angles imply to
within 1.6e-9. With the Sun advancing counter-clockwise (this module before 2026-10-04) the rows
miss by 0.14 to 0.17, with the Sun a quarter turn off by 0.40 to 0.48, with the Sun off by 0.12
to 0.31, and with the mass ratio of ``andreu_default()`` (0.0121505816 against the printed
0.0121506683) by 1.4e-6 to 6.2e-6. A Sun half a turn off is a weak control (2.4e-4 to 2.5e-3),
because the Sun's tidal field is unchanged by that shift to first order in 1/a_S, and is not used.

Monodromy (Fig. 8, p1331, read off the figure): the z0 families have a real pair of modulus about
1.06 and 0.94 at epsilon = 1 ("weakly unstable"); the vz0 families have all six at modulus 1
("linearly stable"). Measured: 1.06045 and 0.94299 for z0; the remaining moduli (all six for
vz0, four for z0) within 2e-13 of 1.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp

FloatArray = NDArray[np.float64]

# Table 1 (p1326), "Parameters used for the CR3BP and BCR4BP (Topputo, 2013)".
_OMEGA_S = 0.925195985  # printed as -0.925195985; the sign is the module's convention
_SYSTEM = bcr4bp.BCR4BPSystem(
    mu=0.0121506683,
    mu_sun=328900.541,
    a_sun_nondim=388.811143,
    omega_sun_nondim=_OMEGA_S,
    theta_sun0=0.0,
)
_PERIOD = 2.0 * math.pi / _OMEGA_S
_TOL = 1e-13

# Tables 2 to 5 (p1332), rows #1 to #4: (x, z, vx, vy, vz, theta_S), y = 0.
_TABLES: dict[str, list[tuple[float, float, float, float, float, float]]] = {
    "table2_z0_S0": [
        (-1.107328855, 0.0, 0.0, 2.056256309, -0.156230339, 0.0),
        (1.111203054, -0.180245301, -0.000134973, -2.062045643, -0.000001655, 4.712253115),
        (-1.106981252, 0.0, 0.0, 2.056236796, 0.156277994, 3.141592654),
        (1.111203054, 0.180245301, 0.000134973, -2.062045643, -0.000001655, 1.570932192),
    ],
    "table3_z0_Spi": [
        (-1.106981252, 0.0, 0.0, 2.056236796, -0.156277994, 3.141592654),
        (1.111203054, -0.180245301, 0.000134973, -2.062045643, 0.000001655, 1.570932192),
        (-1.107328855, 0.0, 0.0, 2.056256309, 0.156230339, 0.0),
        (1.111203054, 0.180245301, -0.000134973, -2.062045643, 0.000001655, 4.712253115),
    ],
    "table4_vz0_S0": [
        (1.090174251, -0.204803847, 0.0, -2.061909684, 0.0, 0.0),
        (-1.120233045, 0.000079580, -0.000178477, 2.042532822, 0.177628656, 4.712573766),
        (1.090649738, 0.204909100, 0.0, -2.061914819, 0.0, 3.141592654),
        (-1.120233045, 0.000079580, 0.000178477, 2.042532822, -0.177628656, 1.570611541),
    ],
    "table5_vz0_Spi": [
        (1.090649738, -0.204909100, 0.0, -2.061914819, 0.0, 3.141592654),
        (-1.120233045, -0.000079580, 0.000178477, 2.042532822, 0.177628656, 1.570611541),
        (1.090174251, 0.204803847, 0.0, -2.061909684, 0.0, 0.0),
        (-1.120233045, -0.000079580, -0.000178477, 2.042532822, -0.177628656, 4.712573766),
    ],
}

# Every row except Table 4 #1, whose closure is tested in test_bcr4bp_sun_sense.py.
_ROWS = [
    pytest.param(table, k, id=f"{table}_row{k + 1}")
    for table in _TABLES
    for k in range(4)
    if (table, k) != ("table4_vz0_S0", 0)
]


def _row(table: str, k: int) -> tuple[FloatArray, float]:
    x, z, vx, vy, vz, theta = _TABLES[table][k % 4]
    return np.array([x, 0.0, z, vx, vy, vz]), theta


def _propagate(system: bcr4bp.BCR4BPSystem, state: FloatArray, t: float) -> FloatArray:
    out: FloatArray = bcr4bp.propagate_bcr4bp(system, state, t, rtol=_TOL, atol=_TOL).state_f
    return out


def _time_between(theta_from: float, theta_to: float) -> float:
    """Time for the Sun angle to fall from ``theta_from`` to ``theta_to`` (mod 2 pi)."""
    return ((theta_from - theta_to) % (2.0 * math.pi)) / _OMEGA_S


def _miss(system: bcr4bp.BCR4BPSystem, state: FloatArray) -> float:
    return float(np.linalg.norm(_propagate(system, state, _PERIOD) - state))


@pytest.mark.parametrize(("table", "k"), _ROWS)
def test_printed_row_closes_after_one_sun_period(table: str, k: int) -> None:
    """Started at its printed Sun angle, each printed crossing returns to itself after
    ``2 pi / |omega_S|``. Controls: the Sun advancing the wrong way, and the Sun a quarter turn
    from its printed angle, both miss by more than 0.1."""
    state, theta = _row(table, k)
    system = dataclasses.replace(_SYSTEM, theta_sun0=theta)
    assert _miss(system, state) < 1e-7
    wrong_sense = dataclasses.replace(system, omega_sun_nondim=-_OMEGA_S)
    assert _miss(wrong_sense, state) > 0.1
    quarter_turn = dataclasses.replace(system, theta_sun0=theta + 0.5 * math.pi)
    assert _miss(quarter_turn, state) > 0.1


@pytest.mark.parametrize(("table", "k"), _ROWS)
def test_printed_row_needs_the_printed_mass_ratio(table: str, k: int) -> None:
    """The closure resolves the seventh significant digit of mu: with ``andreu_default()``'s
    mass ratio and Oshima's other constants the rows miss by 1.4e-6 or more (measured), against
    at most 2.7e-8 with the printed one."""
    state, theta = _row(table, k)
    system = dataclasses.replace(_SYSTEM, theta_sun0=theta, mu=bcr4bp.andreu_default().mu)
    assert _miss(system, state) > 5e-7


@pytest.mark.parametrize(("table", "k"), [(t, k) for t in _TABLES for k in range(4)])
def test_consecutive_printed_crossings_are_one_trajectory(table: str, k: int) -> None:
    """Row k propagated for the time implied by the fall of the Sun angle lands on row k + 1
    (row #4 on row #1). The Sun-angle convention is load-bearing: reading the angle as rising
    instead (time ``(theta_{k+1} - theta_k) mod 2 pi``) would put rows #2 and #4 in swapped
    order and the arrival far from the printed row."""
    state, theta = _row(table, k)
    target, theta_next = _row(table, k + 1)
    system = dataclasses.replace(_SYSTEM, theta_sun0=theta)
    arrival = _propagate(system, state, _time_between(theta, theta_next))
    assert float(np.linalg.norm(arrival - target)) < 1e-7
    rising = _propagate(system, state, _time_between(theta_next, theta))
    assert float(np.linalg.norm(rising - target)) > 0.1


@pytest.mark.parametrize("table", list(_TABLES))
def test_y_crossings_fall_at_the_times_the_sun_angles_imply(table: str) -> None:
    """Over one period from row #1 the orbit crosses ``y = 0`` exactly three times in between,
    at the times given by the printed Sun angles of rows #2, #3 and #4, and at no other time."""
    state, theta = _row(table, 0)
    system = dataclasses.replace(_SYSTEM, theta_sun0=theta)

    def y_zero(t: float, s: FloatArray, _system: bcr4bp.BCR4BPSystem) -> float:
        return float(s[1])

    # solve_ivp passes ``args`` to the event too; scipy-stubs does not model that.
    sol = solve_ivp(  # type: ignore[call-overload]
        bcr4bp.bcr4bp_eom,
        (0.0, _PERIOD),
        state,
        args=(system,),
        method="DOP853",
        rtol=_TOL,
        atol=_TOL,
        events=y_zero,
    )
    crossings = [t for t in sol.t_events[0] if 1e-6 < t < _PERIOD - 1e-6]
    implied = [_time_between(theta, _row(table, j)[1]) for j in (1, 2, 3)]
    assert len(crossings) == 3
    assert np.max(np.abs(np.array(crossings) - np.array(implied))) < 1e-8


@pytest.mark.parametrize(
    ("table", "unstable"),
    [
        ("table2_z0_S0", True),
        ("table3_z0_Spi", True),
        ("table4_vz0_S0", False),
        ("table5_vz0_Spi", False),
    ],
)
def test_monodromy_moduli_match_fig_8(table: str, unstable: bool) -> None:
    """Fig. 8 (p1331), values read off the figure at epsilon = 1: the z0 families are weakly
    unstable with a real pair near 1.06 and 0.94 (tolerance 0.01, the figure's legibility);
    the vz0 families are linearly stable, all six moduli equal to 1."""
    state, theta = _row(table, 0)
    system = dataclasses.replace(_SYSTEM, theta_sun0=theta)
    arc = bcr4bp.propagate_bcr4bp(system, state, _PERIOD, with_stm=True, rtol=_TOL, atol=_TOL)
    assert arc.stm is not None
    moduli = np.sort(np.abs(np.linalg.eigvals(arc.stm)))
    if unstable:
        assert abs(moduli[-1] - 1.06) < 0.01
        assert abs(moduli[0] - 0.94) < 0.01
        assert np.all(np.abs(moduli[1:-1] - 1.0) < 1e-8)
    else:
        assert np.all(np.abs(moduli - 1.0) < 1e-8)
