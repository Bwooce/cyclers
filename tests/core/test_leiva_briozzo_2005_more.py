"""Leiva & Briozzo (2005): the printed three-body orbits, and the Earth distances, lunar distances
and stability parameters of the two quasi-bicircular orbits (#896).

Source: A. M. Leiva and C. B. Briozzo, "Fast periodic transfer orbits in the Sun-Earth-Moon
Quasi-Bicircular Problem", Celestial Mechanics and Dynamical Astronomy 91:357-372 (2005), DOI
10.1007/s10569-004-7818-3. Filed in the private paper corpus as
leiva-briozzo-2005-fast-periodic-transfer-orbits-sun-earth-moon-quasi-bicircular-problem-cmda-91-357-doi-10.1007-s10569-004-7818-3.pdf.
The paper has no numbered tables; the numbers are in the text of Sect. 4.1 (p365) and Sect. 4.3
(p368-369), read from the page images (signs included; the text layer drops minus signs).

Already tested in ``tests/core/test_qbcp.py``: the closure of the two quasi-bicircular orbits and
their closest approach to the Moon (loosely, 0.030 to 0.033 length units). This file adds the rest
of what is printed and cheap.

Frames: Sect. 4.1 and 4.2 use the Earth at -mu and the Moon at 1 - mu, the project's frame; Sect.
4.3 uses the Earth at +mu and the Moon at -1 + mu, so its x, y, xdot and ydot change sign here.
The paper prints mu only as "~0.012150" (p360); the module default 0.012150581623433623 is used.

Printed but not tested: the speeds at closest approach, 2023 and 2024 m/s at the Moon and 2443
and 2446 m/s at the Earth (p368-369). The velocity unit is printed only as "~1024 m/s" and the
frame of the speed is not stated. Measured rotating-frame speeds are 2023.6, 2021.9, 2443.2 and
2444.2 m/s with 384,400 km per sidereal month / 2 pi (1023.16 m/s), or 2025.3, 2023.6, 2445.2 and
2446.3 m/s with 1024 m/s; neither unit gives all four to the printed metre per second.
"""

from __future__ import annotations

import functools
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from cyclerfinder.core import cr3bp, qbcp

MU = qbcp.qbcp_default().mu
KM_PER_LENGTH_UNIT = 384400.0  # p360: "length units of ~384,400 km"
# The Sun period T_sun = 2 pi / n with n = 0.9251959855 (p362), printed as 6.7911939 (p363).
SUN_PERIOD = qbcp.qbcp_default().sun_period_tu


def _cr3bp_miss(state0: np.ndarray, duration: float) -> float:
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, duration),
        state0,
        args=(MU,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    return float(np.linalg.norm(sol.y[:, -1] - state0))


# Sect. 4.3 (p368), in the frame with the Earth at +mu: "the RTBP T_sun-periodic orbit with
# x0 = -1.02379270, y0 = 0, xdot0 = 0, ydot0 = 1.91110553". Sign-reversed for this frame.
_SEED = np.array([1.02379270, 0.0, 0.0, 0.0, -1.91110553, 0.0])


def test_seed_three_body_orbit_returns_after_one_sun_period() -> None:
    """The three-body orbit from which both quasi-bicircular orbits were continued returns to its
    printed state after one Sun period.

    Measured: 4.7e-5 in the four state components (mostly xdot); 6.3e-5 at mu = 0.0121505482
    and 3.2e-4 at mu = 0.01215. The paper gives no closure figure for this state and prints it to
    nine digits, so it is a good but not a precise periodic orbit. Control: after the period of
    the localisation orbit (6.372441, below) the same state misses by 1.8.
    """
    assert _cr3bp_miss(_SEED, SUN_PERIOD) < 1e-4
    assert _cr3bp_miss(_SEED, 6.372441) > 1.0


# Sect. 4.1 (p365), in the project's frame: "the PO with x0 = 1.107569 and ydot = -1.644251. This
# PO has a Jacobi constant h = -0.245294 and a period tau = 6.372441". (The paper calls h, the value
# of the Hamiltonian of its Eq. 1, the Jacobi constant; C = -2h.) It lies on the section
# y = 0, xdot = 0, ydot < 0.
_LOCALISATION = np.array([1.107569, 0.0, 0.0, 0.0, -1.644251, 0.0])


def test_localisation_orbit_energy() -> None:
    """The printed h follows from the printed x0 and ydot.

    Measured: h = -0.2452952, 1.2e-6 from the printed -0.245294 (one unit in the last digit).
    Budget from printing alone: x0 and ydot are rounded to 5e-7, which moves h by up to
    |dh/dx| 5e-7 + |ydot| 5e-7 = 0.53 x 5e-7 + 1.64 x 5e-7 = 1.1e-6, and h itself is rounded to
    5e-7; 1.6e-6 in all.
    """
    h = -0.5 * cr3bp.jacobi_constant(_LOCALISATION, MU)
    assert abs(h - (-0.245294)) < 1.6e-6


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#896: the printed localisation orbit (x0 = 1.107569, ydot0 = -1.644251) misses its "
        "start by 3.5e-2 after the printed period 6.372441 (y by -1.0e-2, xdot by -3.3e-2) and "
        "crosses y = 0 at the half period with xdot = 4.5e-3, so it is not the symmetric periodic "
        "orbit it is said to be; the symmetric periodic orbit at the printed h has x0 = 1.110654 "
        "and period 6.370242 (the same to 1e-5 at mu = 0.01215)"
    ),
)
def test_localisation_orbit_returns_after_its_period() -> None:
    """Printed numbers taken at face value: the orbit should return after tau = 6.372441. Seven
    printed digits and a stable orbit (the paper: elliptic below h1 = -0.034536) allow a miss
    near 1e-6; 1e-4 is generous."""
    assert _cr3bp_miss(_LOCALISATION, 6.372441) < 1e-4


# Sect. 4.3 (p368-369): the two quasi-bicircular orbits in the frame with the Earth at +mu, as
# (start time in Sun periods, x0, ydot0) with y0 = xdot0 = 0, and the printed minimum distance to
# the Earth's centre (km), minimum distance to the lunar surface (km) and stability parameters.
_ORBITS = {
    "orbit 1": (0.0, -1.01950751115, 1.97782573253, 134588.0, 10431.0, (5.496, 2.069)),
    "orbit 2": (0.5, -1.01940558303, 1.97615899365, 134381.0, 10400.0, (5.531, 2.070)),
}


@functools.cache
def _orbit(name: str) -> tuple[float, float, float, tuple[float, float]]:
    """Closure in position, minimum Earth and Moon centre distances (length units) and the two
    stability parameters s = lambda + 1/lambda, sorted by size, of the planar monodromy matrix,
    from one propagation over one Sun period."""
    start_fraction, x, ydot = _ORBITS[name][:3]
    system = qbcp.qbcp_default()
    t0 = start_fraction * SUN_PERIOD
    state0 = qbcp.state_pv_to_pm(np.array([-x, 0.0, 0.0, 0.0, -ydot, 0.0]), t0, system)
    sol = solve_ivp(
        qbcp.qbcp_stm_eom,
        (t0, t0 + SUN_PERIOD),
        np.concatenate([state0, np.eye(6).ravel()]),
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    assert sol.sol is not None
    dense = sol.sol
    end = sol.y[:, -1]
    closure = math.hypot(end[0] - state0[0], end[1] - state0[1])

    def minimum(centre: float) -> float:
        def distance(t: float) -> float:
            s = dense(t)
            return math.hypot(float(s[0]) - centre, float(s[1]))

        times = np.linspace(t0, t0 + SUN_PERIOD, 20001)
        samples = dense(times)
        values = np.hypot(samples[0] - centre, samples[1])
        i = int(np.argmin(values))
        lo, hi = times[max(i - 1, 0)], times[min(i + 1, len(times) - 1)]
        refined = minimize_scalar(
            distance, bounds=(lo, hi), method="bounded", options={"xatol": 1e-13}
        )
        return min(float(refined.fun), float(values[i]))

    # The monodromy matrix in canonical variables has the eigenvalues of the one in velocities,
    # since the change of variables is the same at t0 and t0 + T_sun.
    planar = [0, 1, 3, 4]
    monodromy = end[6:].reshape(6, 6)[np.ix_(planar, planar)]
    eigenvalues = sorted(np.linalg.eigvals(monodromy), key=abs, reverse=True)
    s1, s2 = (abs(complex(v + 1.0 / v)) for v in eigenvalues[:2])
    return closure, minimum(-MU), minimum(1.0 - MU), (s1, s2)


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(
            "orbit 1",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "#896: the minimum distance to the Earth's centre is 134,582.8 km along the "
                    "orbit integrated from the printed state against the printed 134,588 km; the "
                    "5.2 km exceeds the 1.2 km budget (134,583.0 km at mu = 0.0121505482)"
                ),
            ),
        ),
        pytest.param(
            "orbit 2",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "#896: the minimum distance to the Earth's centre is 134,377.7 km along the "
                    "orbit integrated from the printed state against the printed 134,381 km; the "
                    "3.3 km exceeds the 1.9 km budget (134,378.0 km at mu = 0.0121505482)"
                ),
            ),
        ),
    ],
)
def test_earth_distance(name: str) -> None:
    """The printed minimum distance to the Earth's centre (p368-369), against the minimum along the
    orbit times 384,400 km. Budget: half a kilometre for the printing plus the orbit's own
    closure miss in position (1.8e-6 and 3.6e-6, 0.7 and 1.4 km), the rule used for the 2008
    Table 5 test. Both printed values lie above the computed ones, by 5.2 and 3.3 km; the
    difference between the two orbits is 207 km printed and 205.1 km computed. The length unit is
    printed only as "~384,400 km"; a scale alone does not reconcile both values.
    """
    closure, earth, _, _ = _orbit(name)
    budget = 0.5 + KM_PER_LENGTH_UNIT * closure
    assert abs(earth * KM_PER_LENGTH_UNIT - _ORBITS[name][3]) <= budget


@pytest.mark.xfail(
    strict=True,
    reason=(
        "#896: the printed starting abscissas place orbit 1 39.2 km farther from the Moon's "
        "centre than orbit 2 (the closest approach is the starting point in both, an x-axis "
        "crossing), but the printed distances to the lunar surface differ by 31 km (10,431 and "
        "10,400 km). Orbit 1 alone matches: 10,432.0 km with a 1737.4 km radius"
    ),
)
def test_lunar_distance_difference() -> None:
    """The lunar radius is not printed, so test the difference of the two printed distances to
    the lunar surface, which does not depend on it. Budget: 1 km for the printed rounding of
    the two values plus the closure miss of each orbit in position (0.7 and 1.4 km), 3.1 km
    in all.

    Dividing each centre distance by alpha_6 at its start time (a pulsating scale) would give a
    difference of 32.1 km, near the printed 31, but then orbit 1 alone would need a lunar radius
    of 1637 km. The 7 km is therefore not explained by that scale factor either.
    """
    closure_1, _, moon_1, _ = _orbit("orbit 1")
    closure_2, _, moon_2, _ = _orbit("orbit 2")
    computed = (moon_1 - moon_2) * KM_PER_LENGTH_UNIT
    budget = 1.0 + KM_PER_LENGTH_UNIT * (closure_1 + closure_2)
    assert abs(computed - (10431.0 - 10400.0)) <= budget


def test_closest_lunar_approach_is_the_start() -> None:
    """Both orbits pass closest to the Moon at their printed starting point, the x-axis crossing
    behind the Moon (or, equivalently, at its return one period later, to within the closure
    miss), so the printed lunar distance is a property of the printed x0 alone."""
    for name, (_, x, *_rest) in _ORBITS.items():
        _, _, moon, _ = _orbit(name)
        assert abs(moon - abs(-x - (1.0 - MU))) < 1e-5


_S_XFAIL = (
    "#896: |s{k}| of {orbit} is {got} from the monodromy matrix of the orbit integrated from the "
    "printed state, against the printed {want}"
)


@pytest.mark.parametrize(
    ("name", "k"),
    [
        pytest.param(
            "orbit 1",
            0,
            marks=pytest.mark.xfail(
                strict=True,
                reason=_S_XFAIL.format(k=1, orbit="orbit 1", got=5.5036, want=5.496)
                + " (5.5024 at mu = 0.0121505482)",
            ),
        ),
        pytest.param(
            "orbit 1",
            1,
            marks=pytest.mark.xfail(
                strict=True,
                reason=_S_XFAIL.format(k=2, orbit="orbit 1", got=2.0632, want=2.069)
                + " (2.0558 at mu = 0.0121505482)",
            ),
        ),
        pytest.param(
            "orbit 2",
            0,
            marks=pytest.mark.xfail(
                strict=True,
                reason=_S_XFAIL.format(k=1, orbit="orbit 2", got=5.5401, want=5.531)
                + " (5.5389 at mu = 0.0121505482)",
            ),
        ),
        ("orbit 2", 1),
    ],
)
def test_stability_parameters(name: str, k: int) -> None:
    """The printed stability parameters |s1|, |s2| (p369), s = lambda + 1/lambda for the two
    eigenvalue pairs of the planar monodromy matrix over one Sun period, to the printed three
    decimals.

    Measured: 5.5036 and 2.0632 for orbit 1 (printed 5.496 and 2.069), 5.5401 and 2.0697 for
    orbit 2 (printed 5.531 and 2.070). The differences, 0.1 to 0.3 percent, are far larger than
    anything the orbits' closure misses (4e-6 at most) could produce; they also stay at the paper's mass
    ratio. The change from orbit 1 to orbit 2 agrees better: +0.0365 computed against +0.035
    printed for |s1|. s1 is negative and s2 positive, so both orbits are unstable through a pair
    of real eigenvalues (-5.32, -0.188) and (1.285, 0.778).
    """
    _, _, _, s_values = _orbit(name)
    assert abs(s_values[k] - _ORBITS[name][5][k]) <= 5e-4
