"""Leiva & Briozzo (2008), Tables 1, 3, 4 and 5, as tests of the three-body and quasi-bicircular
models (#896).

Source: A. M. Leiva and C. B. Briozzo, "Extension of fast periodic transfer orbits from the
Earth-Moon RTBP to the Sun-Earth-Moon Quasi-Bicircular Problem", Celestial Mechanics and
Dynamical Astronomy 101:225-245 (2008), DOI 10.1007/s10569-008-9134-9. Filed in the private paper
corpus as
leiva-briozzo-2008-extension-fast-periodic-transfer-orbits-earth-moon-rtbp-to-sun-earth-moon-qbcp-cmda-101-225-doi-10.1007-s10569-008-9134-9.pdf.
Table 1 is on p234, Table 2 on p238, Tables 3 and 4 on p240, Table 5 on p241. Every expected value
below is the paper's, typed digit by digit from the page and checked against the text layer.
Table 2 closure is already tested in ``tests/core/test_qbcp.py``; this file reuses the printed
Table 2 states only for their Table 5 distances.

The paper's frame has the Earth at +mu and the Moon at -1 + mu (p228); the project's frames are
rotated by pi, so x, y, xdot and ydot all change sign. The velocities are time derivatives of the
synodic coordinates (p229, p238).

Mass ratio. The paper prints mu only as "~0.0121505" (p228). The Table 1 section abscissa
x = -0.836915310 (note to Table 1) is the L1 abscissa for mu = 0.0121505482 (to 1.5e-10); for the
project's Earth-Moon mu = 0.012150581623433623 the L1 abscissa is 0.836915145, 1.6e-7 away. The
three-body orbits are tested at 0.0121505482. With the project's mu instead, every one of the 34
orbits closes worse, by a factor of 12 to 1100: the misses (position and velocity together) run
from 2.9e-5 to 3.4e-2, against 3.1e-8 to 4.6e-4 at the paper's value.

The quasi-bicircular checks use ``core.qbcp`` with its coefficient tables (Andreu 1998) and mu set
to the paper's value, which is how the paper describes its own model: Andreu's coefficients
(p229) with its own mass ratio. With that mu the Table 3 and 4 arcs return to within 7e-7 to 1.7e-4
in position (the Table 2 orbits to 5e-6 to 6e-5); with the module's default mu they return to
within 4.7e-6 to 2.1e-3 (Table 2: 5e-6 to 3e-4), and the Table 5 distances agree less well.
"""

from __future__ import annotations

import dataclasses
import functools
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from cyclerfinder.core import cr3bp, qbcp

MU_PAPER = 0.0121505482
MU_PROJECT = qbcp.qbcp_default().mu
SECTION_X = 0.836915310  # printed as -0.836915310 in the paper's frame
KM_PER_LENGTH_UNIT = 384400.0  # p228: "length units of ~384400 km"
# p228: n = -Omega = 0.925195985520347; one Sun period T = 2 pi / n (printed as 6.7911939).
SUN_PERIOD = 2.0 * math.pi / 0.925195985520347

# Table 1 (p234): number, h, y, ydot, tau / T_sun as (p, q). Section x = -0.836915310, xdot < 0.
_TABLE_1 = [
    ("053_1", -1.56587846, 0.0746237099, 0.0382328315, 5, 2),
    ("053_2", -1.58753537, 0.0348725952, -0.0149583173, 5, 2),
    ("077_1", -1.58840213, -0.0041344629, 0.0619882192, 5, 2),
    ("077_2", -1.57183324, -0.0446025645, 0.125546118, 5, 2),
    ("084_1", -1.57183324, 0.0485238196, 0.0516615071, 5, 2),
    ("084_2", -1.58840213, 0.0357836953, -0.0150779150, 5, 2),
    ("180A_1", -1.59005198, 0.00520342002, 0.0479318298, 5, 2),
    ("180A_2", -1.57583831, -0.0283283340, 0.117872065, 5, 2),
    ("357", -1.52791268, 0.129037155, 0.115396082, 5, 2),
    ("146A", -1.59170073, 0.00684686737, 0.0389091502, 3, 1),
    ("018A", -1.58890359, -0.0348199799, 0.0299219506, 7, 2),
    ("144A", -1.58963300, -0.0414498069, -0.00372442160, 7, 2),
    ("147", -1.59196333, 0.00534538143, 0.0491647641, 7, 2),
    ("157A_1", -1.59096752, 0.0111493986, 0.0260247945, 7, 2),
    ("157A_2", -1.58960218, 0.0207084632, 0.00687076626, 7, 2),
    ("187A", -1.59288010, -0.0159741427, -0.0215813955, 7, 2),
    ("187B_1", -1.59003721, 0.0220349307, -0.0202753029, 7, 2),
    ("187B_2", -1.59003721, -0.0359084333, -0.0350860129, 7, 2),
    ("188A_1", -1.59064811, 0.0209945181, -0.0103532292, 7, 2),
    ("188A_2", -1.59316603, 0.00892212868, 0.00616109122, 7, 2),
    ("194", -1.59203408, -0.0208044831, -0.0327233527, 7, 2),
    ("244", -1.59203408, 0.000196827658, 0.0549913768, 7, 2),
    ("305", -1.59029796, -0.0347974421, 0.0281011828, 7, 2),
    ("013", -1.58740571, -0.0399746624, -0.0441622383, 4, 1),
    ("020", -1.58710323, -0.0369266173, -0.00906635816, 4, 1),
    ("021", -1.58710323, -0.0446655750, -0.0470435962, 4, 1),
    ("171_1", -1.58930426, 0.0275869044, -0.0153042817, 4, 1),
    ("171_2", -1.59262227, 0.00719782489, 0.0385102976, 4, 1),
    ("081A", -1.58708764, -0.0494577087, -0.0323811597, 9, 2),
    ("172", -1.58857380, 0.0325459456, -0.0264123993, 9, 2),
    ("032B_1", -1.58703219, -0.0280775876, 0.00719954961, 5, 1),
    ("032B_2", -1.58703219, -0.0191844708, -0.00487293371, 5, 1),
    ("287", -1.59413574, 0.00344758201, -0.000982857301, 5, 1),
    ("238", -1.59331122, -0.00774529571, -0.0196143210, 11, 2),
]

# Table 2 (p238): the eleven QBCP periodic orbits. Number, period in Sun periods, t_i, x, xdot,
# y, ydot (paper frame). Same data as ``_LEIVA_BRIOZZO_2008_TABLE_2`` in test_qbcp.py.
_TABLE_2 = [
    ("146A_t3", 3, 3.32657957, -0.833881068, -0.0658016572, -0.00176277248, 0.0366759999),
    ("146A_t4", 3, 6.72217651, -0.833912081, -0.0658968551, -0.00207184538, 0.0369412534),
    ("013_t3", 4, 1.92708674, -0.841058432, -0.0710601802, -0.0415648661, -0.0231934953),
    ("013_t4", 4, 5.32268367, -0.841255581, -0.0709901229, -0.0417347404, -0.0224675395),
    ("020_t1", 4, 1.37929365, -0.840861762, -0.0890884586, -0.0359687313, 0.0101908036),
    ("020_t2", 4, 4.77489058, -0.841778236, -0.0889071284, -0.0358040150, 0.0138070651),
    ("171_2_t3", 4, 1.73158585, -0.837975852, -0.0330903804, 0.0107192583, 0.0366443351),
    ("171_2_t4", 4, 5.12718279, -0.838292090, -0.0322885877, 0.0111104984, 0.0373817317),
    ("032B_1_t4", 5, 6.35357835, -0.824049581, -0.112839541, -0.0276354646, -0.0392670699),
    ("053d_2_t3", 5, 3.28799799, -0.833040875, -0.103358636, 0.0365770176, -0.0240383643),
    ("053d_2_t4", 5, 6.68359492, -0.832339057, -0.102739753, 0.0372962615, -0.0265393587),
]

# Tables 3 and 4 (p240): the periodic arcs. Number, arc length in Sun periods (tau = 5/2 or 7/2
# T_sun), t_i, x, xdot, y, ydot (paper frame).
#
# Row 187A_t1 of Table 4 is left out. It is printed with xdot = -0.0371102305 and
# y = -0.0371102305, the same ten digits in two independent variables; the page image and the
# text layer both show it that way (checked 2026-10-04). Taken as printed, the state misses its
# return by 1.25 in position, so at least one of the two numbers is a misprint. Measured
# (2026-10-04, paper's mu): keeping the printed xdot and solving for y alone gives
# y = -0.0175001 and a return to 1.4e-5 (four state components); keeping the printed y and solving
# for xdot alone leaves a miss of 8.6e-2. So the printed y is the misprint, and the true value is
# close to 187A_t2's y = -0.0175746564. Not used as a control, since the digits are not known.
# Its Table 5 row is left out for the same reason.
_TABLES_3_4 = [
    ("053_1_t3", 2.5, 2.94372263, -0.823291283, -0.204967189, 0.0729037633, 0.0157937233),
    ("053_1_t4", 2.5, 6.33931957, -0.822940134, -0.205554219, 0.0728395206, 0.0160284110),
    ("053_2_t3", 2.5, 3.28799799, -0.832638313, -0.102970156, 0.0369785750, -0.0257859571),
    ("053_2_t4", 2.5, 6.68359492, -0.832737737, -0.103123506, 0.0368880646, -0.0248040354),
    ("077_1_t1", 2.5, 1.13012471, -0.807338411, -0.141817858, -0.00903646943, -0.00511122361),
    ("077_1_t2", 2.5, 4.82781242, -0.839269824, -0.0853997094, 0.000509872176, 0.0631047327),
    ("084_2_t3", 2.5, 3.34190762, -0.837184965, -0.0936535576, 0.0340469056, -0.0151707186),
    ("084_2_t4", 2.5, 6.73750455, -0.837104209, -0.0935333914, 0.0341170567, -0.0146415891),
    ("180A_1_t1", 2.5, 1.51986327, -0.838273181, -0.0700359357, 0.0103226185, 0.0420907706),
    ("180A_1_t2", 2.5, 4.91546020, -0.837710299, -0.0705053126, 0.00978465371, 0.0414271551),
    ("180A_2_t1", 2.5, 1.17849132, -0.845704702, -0.126198028, -0.0202900217, 0.146995322),
    ("357_t3", 2.5, 2.82792318, -0.816461641, -0.296940242, 0.129114141, 0.0853142743),
    ("357_t4", 2.5, 6.22352012, -0.815847321, -0.298659169, 0.128445946, 0.0850832278),
    ("018A_t1", 3.5, 1.36400740, -0.838283891, -0.0626291135, -0.0327661089, 0.0319939392),
    ("018A_t2", 3.5, 4.75960434, -0.837560836, -0.0634742122, -0.0331814869, 0.0289125536),
    ("144A_t3", 3.5, 1.89115633, -0.837604853, -0.0458798745, -0.0426290387, -0.00335455397),
    ("144A_t4", 3.5, 5.28675326, -0.837485334, -0.0461119947, -0.0426245640, -0.00379001677),
    ("147_t1", 3.5, 1.65992167, -0.838588024, -0.0358277014, 0.00979166411, 0.0511566376),
    ("147_t2", 3.5, 5.05551861, -0.838334830, -0.0364152988, 0.00928545481, 0.0511342220),
    ("187A_t2", 3.5, 4.12681402, -0.838137701, -0.0374160755, -0.0175746564, -0.0143312619),
    ("188A_1_t3", 3.5, 3.38816739, -0.839117139, -0.0908157085, 0.0174832264, -0.00186802548),
    ("188A_1_t4", 3.5, 6.78376433, -0.838350309, -0.0899598374, 0.0184230777, -0.00458069786),
    ("188A_2_t1", 3.5, 0.157104258, -0.834216258, -0.0478467375, 0.00272691115, 0.00189692582),
    ("188A_2_t2", 3.5, 3.55270119, -0.834027198, -0.0473041023, 0.00229340183, 0.00144173781),
]

# Table 5 (p241): minimal distances d_E to the Earth and d_M to the Moon, in km. 187A_t1 omitted
# (see above).
_TABLE_5 = {
    "053_1_t3": (111191, 18702),
    "053_1_t4": (111038, 18728),
    "053_2_t3": (121705, 4593),
    "053_2_t4": (121558, 4610),
    "077_1_t1": (118266, 3110),
    "077_1_t2": (118325, 3150),
    "084_2_t3": (119167, 4683),
    "084_2_t4": (119052, 4705),
    "180A_1_t1": (118458, 7371),
    "180A_1_t2": (118471, 7365),
    "180A_2_t1": (106956, 23966),
    "357_t3": (96931, 14713),
    "357_t4": (96755, 14716),
    "146A_t3": (121237, 7311),
    "146A_t4": (121361, 7346),
    "018A_t1": (133112, 12440),
    "018A_t2": (133120, 12448),
    "144A_t3": (137294, 10432),
    "144A_t4": (137295, 10442),
    "147_t1": (118788, 3701),
    "147_t2": (118786, 3696),
    "187A_t2": (126668, 6684),
    "188A_1_t3": (122518, 5908),
    "188A_1_t4": (122431, 5816),
    "188A_2_t1": (124181, 6473),
    "188A_2_t2": (124332, 6519),
    "013_t3": (137125, 2729),
    "013_t4": (137126, 2733),
    "020_t1": (134724, 3253),
    "020_t2": (134709, 3250),
    "171_2_t3": (119714, 4199),
    "171_2_t4": (119714, 4197),
    "032B_1_t4": (129841, 725),
    "053d_2_t3": (121518, 4542),
    "053d_2_t4": (121518, 4542),
}

_ROWS = {row[0]: row for row in _TABLE_2 + _TABLES_3_4}


# --- Table 1: three-body orbits -----------------------------------------------------------------


def _l1_abscissa(mu: float) -> float:
    def d_omega(x: float) -> float:
        return (
            x
            - (1.0 - mu) * (x + mu) / abs(x + mu) ** 3
            - mu * (x - 1.0 + mu) / abs(x - 1.0 + mu) ** 3
        )

    return float(brentq(d_omega, 0.5, 0.95, xtol=1e-15, rtol=1e-15))


def test_table_1_section_is_at_l1_for_the_paper_mass_ratio() -> None:
    """The section abscissa printed under Table 1 is the L1 abscissa at mu = 0.0121505482.

    Measured: 0.8369153098 at mu = 0.0121505482 (1.5e-10 from the printed 0.836915310) and
    0.8369151454 at the project's mu (1.6e-7 away). The paper does not call the section L1.
    """
    assert abs(_l1_abscissa(MU_PAPER) - SECTION_X) < 5e-10
    assert abs(_l1_abscissa(MU_PROJECT) - SECTION_X) > 1e-7


def _table_1_closure(row: tuple[str, float, float, float, int, int], mu: float) -> float:
    _, h, y, ydot, p, q = row
    x, y, ydot = SECTION_X, -y, -ydot
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    # Eq. 1 with C = -2h: v^2 = 2h + x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2; xdot > 0 in this frame.
    xdot = math.sqrt(2.0 * h + x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - ydot**2)
    state0 = np.array([x, y, 0.0, xdot, ydot, 0.0])
    assert math.isclose(cr3bp.jacobi_constant(state0, mu), -2.0 * h, rel_tol=1e-14)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, p / q * SUN_PERIOD),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    return float(np.linalg.norm(sol.y[:, -1] - state0))


@pytest.mark.parametrize("row", _TABLE_1, ids=lambda r: r[0])
def test_table_1_orbits_close_in_the_three_body_model(
    row: tuple[str, float, float, float, int, int],
) -> None:
    """Each of the 34 printed three-body orbits returns to its section state after its printed
    period (p/q) T_sun, with xdot solved from the printed h.

    Expected size of the miss: h, y and ydot are printed to 9 or 10 significant figures, so
    rounding alone moves the state by up to about 1e-9 (xdot by about 1e-8 through h), and these
    orbits are unstable over periods of 17 to 37 time units, multiplying that by 1e3 to 1e5.
    Measured at mu = 0.0121505482: 3.1e-8 (053_2) to 4.6e-4 (147) in the four state
    components together. The control is the mass ratio: at the project's mu every orbit misses
    by more (2.9e-5 to 3.4e-2), so the closure resolves mu at the 3e-8 level.
    """
    closure = _table_1_closure(row, MU_PAPER)
    assert closure < 1e-3
    assert _table_1_closure(row, MU_PROJECT) > closure


# --- Tables 2 to 5: quasi-bicircular orbits and periodic arcs -----------------------------------


def _paper_system() -> qbcp.QBCPSystem:
    return dataclasses.replace(qbcp.qbcp_default(), mu=MU_PAPER)


@dataclasses.dataclass(frozen=True)
class _Run:
    closure_position: float
    closure_velocity: float
    earth_distance: float  # minimum along the orbit or arc, centre to centre, length units
    moon_distance: float
    alpha6_at_earth_min: float
    alpha6_at_moon_min: float


def _closure(
    pv0: np.ndarray, t0: float, duration: float, system: qbcp.QBCPSystem, *, dense: bool = False
) -> tuple[float, float, object]:
    state0 = qbcp.state_pv_to_pm(pv0, t0, system)
    sol = solve_ivp(
        qbcp.qbcp_eom,
        (t0, t0 + duration),
        state0,
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=dense,
    )
    # Coordinates and velocities are compared, not the canonical momenta: for an arc the
    # end time differs from the start by a half-integer number of Sun periods, where the
    # momentum-velocity relation (alpha_1, alpha_2, alpha_3) is not the same.
    pv1 = qbcp.state_pm_to_pv(sol.y[:, -1], t0 + duration, system)
    return (
        math.hypot(pv1[0] - pv0[0], pv1[1] - pv0[1]),
        math.hypot(pv1[3] - pv0[3], pv1[4] - pv0[4]),
        sol.sol,
    )


def _minimum_distance(dense: object, t0: float, t1: float, centre: float) -> tuple[float, float]:
    def distance(t: float) -> float:
        s = dense(t)  # type: ignore[operator]
        return math.hypot(float(s[0]) - centre, float(s[1]))

    times = np.linspace(t0, t1, 20001)
    samples = dense(times)  # type: ignore[operator]
    values = np.hypot(samples[0] - centre, samples[1])
    i = int(np.argmin(values))
    lo, hi = times[max(i - 1, 0)], times[min(i + 1, len(times) - 1)]
    refined = minimize_scalar(distance, bounds=(lo, hi), method="bounded", options={"xatol": 1e-13})
    if float(refined.fun) < float(values[i]):
        return float(refined.fun), float(refined.x)
    return float(values[i]), float(times[i])


@functools.cache
def _run(name: str) -> _Run:
    """One propagation of a printed Table 2 orbit or Table 3/4 arc at the paper's mass ratio."""
    _, n_periods, t_i, x, xdot, y, ydot = _ROWS[name]
    system = _paper_system()
    duration = n_periods * system.sun_period_tu
    pv0 = np.array([-x, -y, 0.0, -xdot, -ydot, 0.0])
    dpos, dvel, dense = _closure(pv0, t_i, duration, system, dense=True)
    d_e, t_e = _minimum_distance(dense, t_i, t_i + duration, -system.mu)
    d_m, t_m = _minimum_distance(dense, t_i, t_i + duration, 1.0 - system.mu)
    return _Run(
        closure_position=dpos,
        closure_velocity=dvel,
        earth_distance=d_e,
        moon_distance=d_m,
        alpha6_at_earth_min=float(qbcp.evaluate_alphas(t_e, system)[6]),
        alpha6_at_moon_min=float(qbcp.evaluate_alphas(t_m, system)[6]),
    )


@pytest.mark.parametrize("name", [row[0] for row in _TABLE_2])
def test_table_2_orbits_close_at_the_paper_mass_ratio(name: str) -> None:
    """With the paper's mass ratio the Table 2 orbits return closer than with the module's
    default (test_qbcp.py): measured 5.0e-6 to 6.0e-5 in position against 5e-6 to 3e-4."""
    assert _run(name).closure_position < 1e-4


@pytest.mark.parametrize("name", [row[0] for row in _TABLES_3_4])
def test_tables_3_4_periodic_arcs_return_after_tau(name: str) -> None:
    """The paper's definition (p239): "These arcs are periodic in the sense that after a time
    tau equal to the period of the RTBP PO, their coordinates and velocities return to their
    initial values, but being q = 2 the solar phase at t = tau is pi instead of zero."

    So each printed state, started at its printed time t_i (the clock of ``core.qbcp``, the
    same reading as the Table 2 test in test_qbcp.py), must return in position and velocity
    after tau = 5/2 or 7/2 Sun periods. Measured at the paper's mu: position 7.8e-7 (053_2_t3)
    to 1.7e-4 (187A_t2), velocity 1.8e-6 to 4.3e-4. The bound is that of the Table 2 test.

    Control: the same state started at t_i + T_sun/2, where the Sun is at the opposite phase,
    misses (next test). Started a quarter period late, the misses are 0.1 to 4.
    """
    run = _run(name)
    assert run.closure_position < 1e-3
    assert run.closure_velocity < 1e-3


@pytest.mark.parametrize("name", [row[0] for row in _TABLES_3_4])
def test_tables_3_4_arcs_do_not_return_with_the_sun_at_the_opposite_phase(name: str) -> None:
    """An arc is a fixed point of the map from Sun phase 0 to Sun phase pi only, not a periodic
    orbit: started half a Sun period later (Sun at phase pi, the end condition of the arc), the
    same printed state misses by at least a hundred times more. Measured at the paper's mu:
    misses of 2.7e-3 (077_1_t1) to 1.0 (147_t1) in the four state components, 500 to 12,800
    times the miss from t_i.

    Why the arcs nearly close at all: the dominant part of the Sun's tidal field is quadrupolar
    and repeats every half Sun period, so the maps from phase 0 and from phase pi differ only by
    the odd harmonics (the Sun's octupole and the primaries' odd-harmonic motion). The reversing
    symmetry (x, y, t) -> (x, -y, -t) plays no part: the t_i are not symmetric epochs.
    """
    _, n_periods, t_i, x, xdot, y, ydot = _ROWS[name]
    system = _paper_system()
    pv0 = np.array([-x, -y, 0.0, -xdot, -ydot, 0.0])
    dpos, dvel, _ = _closure(
        pv0, t_i + 0.5 * system.sun_period_tu, n_periods * system.sun_period_tu, system
    )
    run = _run(name)
    assert math.hypot(dpos, dvel) > 100.0 * math.hypot(run.closure_position, run.closure_velocity)


def _distance_budget_km(name: str) -> float:
    # Chosen after seeing the data (2026-10-04), and stated here so the reader can judge it: half
    # a kilometre for the printed rounding, plus the distance by which the integrated orbit itself
    # departs from the paper's, taken as its closure miss in position. The budget exceeds 10 km for
    # 053_1, 357, 144A, 147, 187A_t2, 020 and 053d_2 (14.8 to 63.9 km), where the check is weak.
    return 0.5 + KM_PER_LENGTH_UNIT * _run(name).closure_position


_SAMPLED_MINIMUM = (
    "#896: the printed d_M exceeds the minimum distance to the Moon's centre along the "
    "integrated orbit by {miss} km, more than the {budget} km budget; the printed d_M is larger "
    "than the computed one in every row where they differ by more than 0.6 km, which is what a "
    "minimum taken over the paper's integration steps instead of the continuous orbit would give"
)
_MOON_XFAIL: dict[str, str] = {
    "053_2_t3": _SAMPLED_MINIMUM.format(miss=3.3, budget=0.8),
    "053_2_t4": _SAMPLED_MINIMUM.format(miss=1.4, budget=0.8),
    "077_1_t2": _SAMPLED_MINIMUM.format(miss=7.5, budget=1.7),
    "084_2_t4": _SAMPLED_MINIMUM.format(miss=1.9, budget=1.8),
    "013_t4": _SAMPLED_MINIMUM.format(miss=5.8, budget=2.4),
    "171_2_t3": _SAMPLED_MINIMUM.format(miss=3.7, budget=3.5),
    "032B_1_t4": (
        "#896: the printed d_M is 725 km; the orbit integrated from the printed state passes "
        "483.4 km from the Moon's centre (inside the Moon; the paper says some orbits are Moon "
        "colliders), at 7.3075 time units after t_i, unchanged to 1e-6 km at rtol 3e-14, "
        "atol 1e-15 with the periapsis located by an event on r.v = 0. The orbit closes to "
        "1.4e-5, so the 242 km is not an integration or reproduction error; the paper's minimum "
        "was probably taken over its integration steps, which near a 483 km pass are coarse"
    ),
}


def _params(xfails: dict[str, str]) -> list[object]:
    return [
        pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=xfails[n]))
        if n in xfails
        else n
        for n in _TABLE_5
    ]


@pytest.mark.parametrize("name", list(_TABLE_5))
def test_table_5_earth_distances(name: str) -> None:
    """Minimum distance to the Earth's centre along each orbit (n T_sun from t_i) or arc (tau
    from t_i), times 384,400 km, against the printed d_E (Table 5, p241).

    The paper does not say whether d_E and d_M are to the centre or the surface, or in which
    length scale. Measured: centre distances in plain length units of 384,400 km reproduce all
    35 printed d_E within the budget, 26 of them within 1.4 km; the largest differences are
    -8.3 km (147_t1) and +6.8 km (147_t2) on a 17 km budget. See the control test below.
    """
    run = _run(name)
    computed = run.earth_distance * KM_PER_LENGTH_UNIT
    assert abs(computed - _TABLE_5[name][0]) <= _distance_budget_km(name)


@pytest.mark.parametrize("name", _params(_MOON_XFAIL))
def test_table_5_moon_distances(name: str) -> None:
    """Minimum distance to the Moon's centre, as above, against the printed d_M.

    Measured: 28 of 35 rows within the budget. In every row where computed and printed differ
    by more than 0.6 km, the printed value is the larger (by 0.66 to 7.5 km, and 242 km for
    032B_1_t4); the positive differences are all 0.54 km or less. Seven rows exceed the budget
    and are held as strict expected failures.
    """
    run = _run(name)
    computed = run.moon_distance * KM_PER_LENGTH_UNIT
    assert abs(computed - _TABLE_5[name][1]) <= _distance_budget_km(name)


def test_table_5_distances_are_centre_distances_in_unscaled_units() -> None:
    """Control for the two distance tests: the other readings fail the same budget.

    Distance to the lunar surface (centre distance less 1737.4 km) fails every d_M row, and the
    centre distance divided by alpha_6 at the time of closest approach (a pulsating length
    scale) fails 34 of the 35 d_E rows (by 2.4 to 1015 km); the exception is 187A_t2, whose
    budget is 64 km.
    """
    surface_ok = 0
    scaled_ok = 0
    for name, (d_e, d_m) in _TABLE_5.items():
        run = _run(name)
        budget = _distance_budget_km(name)
        surface = run.moon_distance * KM_PER_LENGTH_UNIT - 1737.4
        surface_ok += abs(surface - d_m) <= budget
        scaled = run.earth_distance * KM_PER_LENGTH_UNIT / run.alpha6_at_earth_min
        scaled_ok += abs(scaled - d_e) <= budget
    assert surface_ok == 0
    assert scaled_ok <= 1
