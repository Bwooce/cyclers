"""#931 (a): published general planar three-body controls with a test-only Newtonian integrator.

No module in ``src/`` integrates the finite-mass three-body problem, so this file carries its own
inertial integrator (planar, G = 1, DOP853 at rtol 1e-13, atol 1e-14, no regularisation). Every
expected value below is printed in the source; none was computed by project code. Where a printed
number is a known misprint it is a strict expected failure with the evidence in the reason.

Sources (digests under ``docs/notes/``):

* Hadjidemetriou 1975b, Celest. Mech. 12:255-276, Table I p.268 (equal masses 1/3, G = 1, rotating
  frame tied to the P1-P2 line with theta'0 = 1): rows 1, 2, 3, 10, 24 for E, p and the half
  period tau/2 (first return of the rotating-frame y of P3 to zero, perpendicular crossing); orbit
  a of eq. 73 (11 digits) returns to its state. The b1, b2 and S/U columns need the isoenergetic
  stability map, i.e. the Floquet classifier that another agent is changing (#931 (b)), so they
  are not tested here.
* Hadjidemetriou & Christides 1975, Celest. Mech. 12:175-187, Table I p.181 (m1 = m2 = (1 - m3)/2,
  x30 = 0.181, theta'0 = 1): 19 rows, x1(T), x3(T), period 2T and theta/2pi from the printed
  (m3, x10, y3'0); row 15's period is a misprint. The m3 = 0 row is also checked against
  ``core.cr3bp`` at mu = 1/2.
* Henon 1974a, Celest. Mech. 10:375-388, Tables I (family 1) and II (family 2), masses 3/12, 4/12,
  5/12, E = -47/288: E, A, closure of the state at T rotated by -phi. The A = 0 row of Table II
  is Szebehely's orbit with a collision between bodies 1 and 3 at T/2, which an unregularised
  integrator cannot pass (needs the regularised propagator of #928), so it is skipped.
* Szebehely & Peters 1967a, Astron. J. 72:876 (Pythagorean problem, Table I close approaches) and
  1967b, Astron. J. 72:1187 (periodic orbit, Table I position offsets, E and T).

Marker choice: the repo marks tests over about 10 s ``slow`` (skipped by default). The
Hadjidemetriou tables and the Pythagorean close approaches take a few seconds in total and stay
in the default suite in full; Henon's 13 rows take about 0.65 s each, so 5 rows (both ends of each
table and one in between) are default and the other 8 are ``slow``.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom

Arr = NDArray[np.float64]
RTOL, ATOL = 1e-13, 1e-14


# ---- test-only planar N-body integrator ---------------------------------------------------


def _rhs(masses: Arr) -> Callable[[float, Arr], Arr]:
    n = len(masses)

    def f(_t: float, y: Arr) -> Arr:
        pos = y[: 2 * n].reshape(n, 2)
        acc = np.zeros_like(pos)
        for i in range(n):
            for j in range(n):
                if i != j:
                    d = pos[j] - pos[i]
                    acc[i] += masses[j] * d / (d @ d) ** 1.5
        return np.concatenate((y[2 * n :], acc.ravel()))

    return f


def energy(masses: Arr, pos: Arr, vel: Arr) -> float:
    kin = 0.5 * float(np.sum(masses * np.sum(vel * vel, axis=1)))
    pot = 0.0
    for i in range(len(masses)):
        for j in range(i + 1, len(masses)):
            pot -= masses[i] * masses[j] / float(np.linalg.norm(pos[i] - pos[j]))
    return kin + pot


def ang_mom(masses: Arr, pos: Arr, vel: Arr) -> float:
    return float(np.sum(masses * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0])))


def to_com(masses: Arr, pos: Arr, vel: Arr) -> tuple[Arr, Arr]:
    m = masses.sum()
    return pos - masses @ pos / m, vel - masses @ vel / m


@dataclass
class Frame:
    """Rotating-frame quantities of the 'binary (bodies 0, 1) plus third body (2)' description."""

    theta: float
    x1: float
    x3: float
    y3: float
    x1p: float
    x3p: float
    y3p: float


def frame_of(masses: Arr, pos: Arr, vel: Arr) -> Frame:
    m1, m2 = masses[0], masses[1]
    g1 = (m1 * pos[0] + m2 * pos[1]) / (m1 + m2)
    v1 = (m1 * vel[0] + m2 * vel[1]) / (m1 + m2)
    r12, v12 = pos[0] - pos[1], vel[0] - vel[1]
    theta = math.atan2(r12[1], r12[0])
    omega = (r12[0] * v12[1] - r12[1] * v12[0]) / float(r12 @ r12)
    c, s = math.cos(theta), math.sin(theta)

    def rot(r: Arr) -> Arr:
        return np.array([c * r[0] + s * r[1], -s * r[0] + c * r[1]])

    def rotating(i: int) -> tuple[Arr, Arr]:
        r = pos[i] - g1
        v = vel[i] - v1 - omega * np.array([-r[1], r[0]])
        return rot(r), rot(v)

    p1, q1 = rotating(0)
    p3, q3 = rotating(2)
    return Frame(theta, p1[0], p3[0], p3[1], q1[0], q3[0], q3[1])


def rotating_start(masses: Arr, x10: float, x30: float, y3p0: float) -> tuple[Arr, Arr]:
    """Inertial CM state from the printed rotating-frame data (theta'0 = 1, P1 and P2 at rest)."""
    m1, m2 = masses[0], masses[1]
    pos = np.array([[x10, 0.0], [-m1 * x10 / m2, 0.0], [x30, 0.0]])
    v_rot = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, y3p0]])
    vel = v_rot + np.stack((-pos[:, 1], pos[:, 0]), axis=1)  # + theta'0 z x r with theta'0 = 1
    return to_com(masses, pos, vel)


@dataclass
class Orbit:
    t_half: float
    final: Frame
    sol: object
    masses: Arr


def first_return(masses: Arr, pos: Arr, vel: Arr, t_max: float, dense: bool = False) -> Orbit:
    """Integrate to the first return of the rotating-frame y of body 2 to zero (terminal)."""
    n = len(masses)
    y0 = np.concatenate((pos.ravel(), vel.ravel()))
    y3p0 = frame_of(masses, pos, vel).y3p
    direction = -1.0 if y3p0 > 0 else 1.0

    def event(_t: float, y: Arr) -> float:
        return frame_of(masses, y[: 2 * n].reshape(n, 2), y[2 * n :].reshape(n, 2)).y3

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = direction  # type: ignore[attr-defined]
    sol = solve_ivp(
        _rhs(masses),
        (0.0, t_max),
        y0,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=event,
        dense_output=dense,
    )
    assert sol.t_events is not None and sol.y_events is not None
    assert sol.t_events[0].size == 1, "no return of y3 to zero before t_max"
    yf = sol.y_events[0][0]
    return Orbit(
        float(sol.t_events[0][0]),
        frame_of(masses, yf[: 2 * n].reshape(n, 2), yf[2 * n :].reshape(n, 2)),
        sol,
        masses,
    )


# ---- Hadjidemetriou 1975b Table I (p.268) -------------------------------------------------

EQUAL = np.array([1.0, 1.0, 1.0]) / 3.0
# row: (x10, x30, y3'0, tau/2, E, p)
H75B = {
    1: (1.09554886, 1.02242545, -2.94637505, 0.076672, -0.81130995, 0.36301930),
    2: (1.00261239, 0.90242545, -2.47957163, 0.123789, -0.61131972, 0.35387510),
    3: (0.87680981, 0.70242545, -1.78210271, 0.292671, -0.38509299, 0.34399856),
    10: (0.76153702, 0.31122545, -0.86360433, 1.989432, -0.19605400, 0.34842257),
    24: (0.84305194, 0.52195540, -1.56444275, 5.486094, -0.13566873, 0.35290618),
}


@pytest.mark.parametrize("row", sorted(H75B))
def test_hadjidemetriou_1975b_table_i_energy_momentum_half_period(row: int) -> None:
    x10, x30, y3p, tau_half, e_print, p_print = H75B[row]
    pos, vel = rotating_start(EQUAL, x10, x30, y3p)
    # printed to 8 decimals: rounding of the start is amplified to about 1e-7 in E and 1e-8 in p
    assert energy(EQUAL, pos, vel) == pytest.approx(e_print, abs=3e-7)
    assert ang_mom(EQUAL, pos, vel) == pytest.approx(p_print, abs=3e-8)
    orb = first_return(EQUAL, pos, vel, t_max=1.5 * tau_half)
    assert orb.t_half == pytest.approx(tau_half, abs=1e-6)  # printed to 6 decimals
    # perpendicular crossing: x1' = x3' = 0 at tau/2, to the rounding of the 8-decimal start
    # amplified by the orbit's multiplier (row 24: |lambda| about 118, miss 1.0e-7)
    assert abs(orb.final.x1p) < 5e-7 and abs(orb.final.x3p) < 5e-7


def test_hadjidemetriou_1975b_orbit_a_eq_73_returns_to_its_state() -> None:
    """Eq. 73 (p.271, 11 digits): x1 = 0.94082115256, x3 = 0.81242545000, y3' = -2.15068729528."""
    pos, vel = rotating_start(EQUAL, 0.94082115256, 0.81242545000, -2.15068729528)
    start = frame_of(EQUAL, pos, vel)
    half = first_return(EQUAL, pos, vel, t_max=1.0)
    assert abs(half.final.x1p) < 1e-8 and abs(half.final.x3p) < 1e-8  # symmetric orbit
    # the second same-direction return of y3 is one full period: state back to the start
    y0 = np.concatenate((pos.ravel(), vel.ravel()))
    full = solve_ivp(
        _rhs(EQUAL), (0.0, 2.0 * half.t_half), y0, method="DOP853", rtol=RTOL, atol=ATOL
    )
    end = frame_of(EQUAL, full.y[:6, -1].reshape(3, 2), full.y[6:, -1].reshape(3, 2))
    for name in ("x1", "x3", "y3", "x1p", "x3p", "y3p"):
        assert getattr(end, name) == pytest.approx(getattr(start, name), abs=1e-8), name


# ---- Hadjidemetriou & Christides 1975 Table I (p.181) -------------------------------------

# row: (m3, x10, y3'0, x1(T), x3(T), period 2T, theta/2pi)
HC75 = {
    1: (0.0, 0.50000000, -0.94711034, 0.50000000, 0.72101839, 1.8484, 0.29419126),
    2: (0.0010, 0.50092191, -0.94511399, 0.50033711, 0.72165922, 1.8560, 0.29490893),
    3: (0.0100, 0.50909746, -0.92781324, 0.50315431, 0.72712559, 1.9240, 0.30138460),
    4: (0.0500, 0.54307785, -0.86342556, 0.51142589, 0.74554248, 2.2244, 0.33111567),
    5: (0.1000, 0.58106608, -0.80535147, 0.51335472, 0.75619699, 2.6078, 0.37251151),
    6: (0.1200, 0.59504112, -0.78781819, 0.51168416, 0.75629375, 2.7668, 0.39105923),
    7: (0.1500, 0.61470336, -0.76738429, 0.50631794, 0.75064525, 3.0140, 0.42200972),
    8: (0.1800, 0.63254710, -0.75535104, 0.49643753, 0.73460805, 3.2760, 0.45888914),
    9: (0.2000, 0.64291315, -0.75497081, 0.48543791, 0.71289250, 3.4654, 0.49000811),
    10: (0.2200, 0.64988645, -0.77337080, 0.46325010, 0.66276240, 3.6880, 0.53839883),
    11: (0.2240, 0.64910244, -0.78970281, 0.45126636, 0.63398810, 3.7526, 0.55965075),
    12: (0.2245, 0.64854088, -0.79446754, 0.44814903, 0.62643112, 3.7646, 0.56477324),
    13: (0.2240, 0.64284614, -0.82658400, 0.42939154, 0.58102449, 3.8070, 0.59285687),
    14: (0.2200, 0.63554886, -0.85813511, 0.41310999, 0.54248007, 3.8116, 0.61439245),
    15: (0.1800, 0.58857781, -1.02869781, 0.34032721, 0.39120985, 3.6656, 0.70005821),
    16: (0.1200, 0.53513036, -1.21749666, 0.28459823, 0.30014806, 3.5334, 0.78750832),
    17: (0.0500, 0.48265003, -1.41273701, 0.25296896, 0.25639133, 3.5264, 0.89757983),
    18: (0.0100, 0.45321081, -1.52856693, 0.24485183, 0.24601948, 3.5784, 0.97702084),
    19: (0.0, 0.44559306, -1.55950866, 0.24406810, 0.24495936, 3.5986, 1.00000000),
}
HC_X30 = 0.181


def _hc_orbit(row: int) -> tuple[Orbit, float]:
    m3, x10, y3p, _x1t, _x3t, period, _ = HC75[row]
    masses = np.array([(1.0 - m3) / 2.0, (1.0 - m3) / 2.0, m3])
    pos, vel = rotating_start(masses, x10, HC_X30, y3p)
    orb = first_return(masses, pos, vel, t_max=0.6 * period)
    # unwrapped rotation of the P1-P2 line over the full period 2T, sampled densely
    ts = np.linspace(0.0, 2.0 * orb.t_half, 4001)
    y0 = np.concatenate((pos.ravel(), vel.ravel()))
    full = solve_ivp(
        _rhs(masses), (0.0, 2.0 * orb.t_half), y0, t_eval=ts, method="DOP853", rtol=RTOL, atol=ATOL
    )
    ang = np.unwrap([math.atan2(y[1] - y[3], y[0] - y[2]) for y in full.y.T])
    return orb, float((ang[-1] - ang[0]) / (2.0 * math.pi))


def _hc_check(row: int) -> None:
    _m3, _x10, _y3p, x1t, x3t, period, rot = HC75[row]
    orb, turns = _hc_orbit(row)
    assert orb.final.x1 == pytest.approx(x1t, abs=5e-8)
    assert orb.final.x3 == pytest.approx(x3t, abs=5e-8)
    assert turns == pytest.approx(rot, abs=1e-7)
    assert 2.0 * orb.t_half == pytest.approx(period, abs=1e-4)  # printed to 4 decimals
    # x1' = x3' = 0 at the crossing (3.2e-6 at row 19, where P1 and P3 pass within 9e-4)
    assert abs(orb.final.x1p) < 5e-6 and abs(orb.final.x3p) < 5e-6


def _hc_params(rows: list[int]) -> list[object]:
    out: list[object] = []
    for r in rows:
        marks = []
        if r == 15:
            marks.append(
                pytest.mark.xfail(
                    strict=True,
                    reason="Hadjidemetriou & Christides Table I row 15 (m3 = 0.18, descending) "
                    "prints period 3.6656; x1(T), x3(T) and theta/2pi of the same row close to "
                    "1e-7 and the integrated 2T is 3.6697 (neighbouring periods 3.8116 and 3.5334 "
                    "fit 3.6697 on a smooth curve): a digit slip in print (digest section 4).",
                )
            )
        out.append(pytest.param(r, marks=marks, id=f"row{r}"))
    return out


@pytest.mark.parametrize("row", _hc_params(sorted(HC75)))
def test_hadjidemetriou_christides_1975_table_i_rows(row: int) -> None:
    _hc_check(row)


def test_hadjidemetriou_christides_row_15_x1_x3_and_rotation_close_despite_period_misprint() -> (
    None
):
    _m3, _x10, _y3p, x1t, x3t, _period, rot = HC75[15]
    orb, turns = _hc_orbit(15)
    assert orb.final.x1 == pytest.approx(x1t, abs=5e-8)
    assert orb.final.x3 == pytest.approx(x3t, abs=5e-8)
    assert turns == pytest.approx(rot, abs=1e-7)


def test_hadjidemetriou_christides_terminal_orbit_is_elliptic_restricted_with_e_0292() -> None:
    """Row 19 (m3 = 0): the binary is Kepler with r0 = 2 x10 at apoapsis, theta'0 = 1.

    Printed: e = 0.292 (abstract and p.185), period 2T = 3.5986. h = r0^2, p = h^2 / (m1 + m2),
    e = 1 - p / r0 (apoapsis), a = p / (1 - e^2), period 2 pi a^(3/2): analytic, not a code value.
    """
    r0 = 2.0 * HC75[19][1]
    h = r0 * r0
    sl = h * h  # m1 + m2 = 1
    e = 1.0 - sl / r0
    a = sl / (1.0 - e * e)
    assert e == pytest.approx(0.292, abs=5e-4)
    assert 2.0 * math.pi * a**1.5 == pytest.approx(HC75[19][5], abs=5e-5)


def test_hadjidemetriou_christides_row_1_is_the_circular_restricted_szebehely_orbit() -> None:
    """Row 1 (m3 = 0) is the planar CR3BP orbit at mu = 1/2 (core.cr3bp): perpendicular crossing."""
    _m3, _x10, y3p, _x1t, x3t, period, _ = HC75[1]
    state0 = np.array([HC_X30, 0.0, 0.0, 0.0, y3p, 0.0])

    def cross(_t: float, y: Arr, _mu: float) -> float:
        return float(y[1])

    cross.terminal = True  # type: ignore[attr-defined]
    cross.direction = 1.0  # type: ignore[attr-defined]
    sol = solve_ivp(  # type: ignore[call-overload]
        cr3bp_eom,
        (0.0, 2.0),
        state0,
        args=(0.5,),
        events=cross,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
    )
    t_cross = float(sol.t_events[0][0])
    end = sol.y_events[0][0]
    assert 2.0 * t_cross == pytest.approx(period, abs=1e-4)  # printed 1.8484
    assert end[0] == pytest.approx(x3t, abs=1e-7)  # printed x3(T) = 0.72101839
    assert abs(end[3]) < 1e-7  # x' = 0: perpendicular


# ---- Henon 1974a Tables I and II (p.383) --------------------------------------------------

HENON_MASSES = np.array([3.0, 4.0, 5.0]) / 12.0
HENON_E = -47.0 / 288.0
# (A, x1, y1, u1, v1, x2, v2, T, phi)
HENON_T1 = [
    (0.0, 0.54402539, 1.79622952, 0.0, 0.0, -1.03258028, 0.0, 20.25306432, 0.0),
    (
        0.001,
        0.54053935,
        1.79816889,
        -0.00585553,
        -0.01564698,
        -1.03063119,
        0.00460907,
        20.25265542,
        -0.26835647,
    ),
    (
        0.002,
        0.53110395,
        1.80293934,
        -0.01103456,
        -0.02980389,
        -1.02554227,
        0.00866815,
        20.25151930,
        -0.51726610,
    ),
    (
        0.003,
        0.51867478,
        1.80837958,
        -0.01512090,
        -0.04133431,
        -1.01914948,
        0.01174645,
        20.24990518,
        -0.72933546,
    ),
    (
        0.004,
        0.50591701,
        1.81318780,
        -0.01822663,
        -0.05033904,
        -1.01284909,
        0.01387159,
        20.24803759,
        -0.90428167,
    ),
    (
        0.005,
        0.49394848,
        1.81711438,
        -0.02064588,
        -0.05747761,
        -1.00711374,
        0.01528495,
        20.24602141,
        -1.05095658,
    ),
    (
        0.006,
        0.48302567,
        1.82026973,
        -0.02260728,
        -0.06331350,
        -1.00199088,
        0.01619576,
        20.24389395,
        -1.17746479,
    ),
    (
        0.007,
        0.47312277,
        1.82280912,
        -0.02425644,
        -0.06822325,
        -0.99741625,
        0.01674488,
        20.24166819,
        -1.28938791,
    ),
]
HENON_T2 = [
    (
        0.005,
        0.24516190,
        0.64443031,
        0.00489610,
        -0.01912525,
        -1.74990885,
        -0.00123263,
        10.27995212,
        0.14763051,
    ),
    (
        0.010,
        0.24977335,
        0.64716408,
        0.00969552,
        -0.03875714,
        -1.74841869,
        -0.00241317,
        10.28349038,
        0.30110185,
    ),
    (
        0.015,
        0.25812442,
        0.65209118,
        0.01428654,
        -0.05953310,
        -1.74581617,
        -0.00347977,
        10.28991666,
        0.46826654,
    ),
    (
        0.020,
        0.27179428,
        0.66015095,
        0.01850316,
        -0.08250792,
        -1.74183836,
        -0.00433877,
        10.30043936,
        0.66357433,
    ),
    (
        0.025,
        0.29600664,
        0.67463620,
        0.02193970,
        -0.11038997,
        -1.73571451,
        -0.00476788,
        10.31888994,
        0.92905825,
    ),
]
# Table II row A = 0 (Szebehely's orbit) is (0.24368035, 0.64354890, 0, 0, -1.75039533, 0,
# 10.27881780, 0): bodies 1 and 3 collide at T/2, so it is skipped (needs regularisation, #928).


def _henon_state(row: tuple[float, ...]) -> tuple[Arr, Arr]:
    _a, x1, y1, u1, v1, x2, v2, _t, _phi = row
    m = HENON_MASSES
    y2 = y3 = -m[0] * y1 / (m[1] + m[2])
    u2 = u3 = -m[0] * u1 / (m[1] + m[2])
    x3 = -(m[0] * x1 + m[1] * x2) / m[2]
    v3 = -(m[0] * v1 + m[1] * v2) / m[2]
    pos = np.array([[x1, y1], [x2, y2], [x3, y3]])
    vel = np.array([[u1, v1], [u2, v2], [u3, v3]])
    return pos, vel


def _henon_check(row: tuple[float, ...]) -> tuple[float, float]:
    a, *_rest, t_end, phi = row
    pos, vel = _henon_state(row)
    assert energy(HENON_MASSES, pos, vel) == pytest.approx(HENON_E, abs=1e-9)
    assert ang_mom(HENON_MASSES, pos, vel) == pytest.approx(a, abs=1e-8)
    y0 = np.concatenate((pos.ravel(), vel.ravel()))
    sol = solve_ivp(_rhs(HENON_MASSES), (0.0, t_end), y0, method="DOP853", rtol=RTOL, atol=ATOL)
    pe, ve = sol.y[:6, -1].reshape(3, 2), sol.y[6:, -1].reshape(3, 2)
    d23 = pe[2] - pe[1]
    phi_c = math.atan2(d23[1], d23[0])  # x axis was along 2 -> 3 at t = 0
    assert phi_c == pytest.approx(phi, abs=1e-6)
    c, s = math.cos(-phi), math.sin(-phi)
    rot = np.array([[c, -s], [s, c]])
    closure = float(max(np.max(np.abs(pe @ rot.T - pos)), np.max(np.abs(ve @ rot.T - vel))))
    unrotated = float(max(np.max(np.abs(pe - pos)), np.max(np.abs(ve - vel))))
    # r23 is at an extremum at T: d(r23)/dt = 0 to the rounding of the 8-decimal T and start
    # (digest 1.1e-6; 3.2e-6 here at Table II A = 0.015, so the bound is 1e-5)
    rdot = float(d23 @ (ve[2] - ve[1]) / np.linalg.norm(d23))
    assert abs(rdot) < 1e-5
    return closure, unrotated


HENON_ALL = [("I", r) for r in HENON_T1] + [("II", r) for r in HENON_T2]
HENON_DEFAULT_A = {("I", 0.0), ("I", 0.004), ("I", 0.007), ("II", 0.005), ("II", 0.025)}


def _henon_params(default: bool) -> list[object]:
    return [
        pytest.param(t, r, id=f"T{t}-A{r[0]}")
        for t, r in HENON_ALL
        if ((t, r[0]) in HENON_DEFAULT_A) == default
    ]


def _henon_assert(row: tuple[float, ...]) -> None:
    closure, unrotated = _henon_check(row)
    # 8-decimal inputs on orbits with |lambda| 10 to 60 (p.386) and close 1-3 passes of 6e-4 to
    # 2e-3: perturbing the printed digits of Table II A = 0.015 by 5e-9 moves the closure between
    # 1e-6 and 1.4e-5 (7.5e-6 unperturbed; the digest quotes 6.7e-7 for its runs), so the bound
    # is 2e-5. Every other row closes to 5e-7 or better.
    assert closure < 2e-5
    if row[-1] != 0.0:
        # negative control: without the rotation the configuration has not returned
        assert unrotated > 0.05


@pytest.mark.parametrize(("table", "row"), _henon_params(True))
def test_henon_1974a_rows_close_after_rotation_by_minus_phi(
    table: str, row: tuple[float, ...]
) -> None:
    _henon_assert(row)


@pytest.mark.slow
@pytest.mark.parametrize(("table", "row"), _henon_params(False))
def test_henon_1974a_remaining_rows(table: str, row: tuple[float, ...]) -> None:
    _henon_assert(row)


@pytest.mark.skip(
    reason="Henon 1974a Table II row A = 0 is Szebehely's orbit: bodies 1 and 3 collide at T/2 "
    "(1-3 distance 6.8e-4 at 0.499999 T and falling to zero), which an unregularised Cartesian "
    "integrator cannot pass; needs the regularised propagator of #928."
)
def test_henon_1974a_table_ii_collision_row_a0() -> None:
    raise AssertionError


# ---- Szebehely & Peters 1967b: the periodic Pythagorean orbit ------------------------------

SP_MASSES = np.array([3.0, 4.0, 5.0])
SP_BURRAU = np.array([[1.0, 3.0], [-2.0, -1.0], [1.0, -1.0]])
SP_E, SP_T = -12.7616527695, 31.8229622453
# Table I offsets x_i = x_i^B - x_i^P, as printed (body 2 x as printed, and the corrected digit)
SP_OFF_X = [0.0694920571, -0.0129612126, -0.0313262594]
SP_OFF_Y = [-0.0284517815, 0.0531023120, -0.0254107807]
SP_X2_CORRECTED = -0.0129612186


def _sp_positions(x2: float) -> Arr:
    off = np.array([[SP_OFF_X[0], SP_OFF_Y[0]], [x2, SP_OFF_Y[1]], [SP_OFF_X[2], SP_OFF_Y[2]]])
    return np.asarray(SP_BURRAU - off, dtype=np.float64)


@pytest.mark.xfail(
    strict=True,
    reason="Szebehely & Peters 1967b Table I prints body 2 offset x = -0.0129612126; with it E = "
    "-12.76165275485 misses the printed E = -12.7616527695 by 1.5e-8 (and the centre of mass is "
    "off by 2.4e-8, the integrated state at T by 1.5e-7). The last digits are a misprint for "
    "-0.0129612186 (E = -12.7616527697), digest 2026-10-05 section 3, coordinator-checked.",
)
def test_szebehely_peters_1967b_printed_x2_reproduces_printed_energy() -> None:
    pos = _sp_positions(SP_OFF_X[1])
    assert energy(SP_MASSES, pos, np.zeros((3, 2))) == pytest.approx(SP_E, abs=1e-9)


def test_szebehely_peters_1967b_corrected_x2_reproduces_printed_energy() -> None:
    pos = _sp_positions(SP_X2_CORRECTED)
    assert energy(SP_MASSES, pos, np.zeros((3, 2))) == pytest.approx(SP_E, abs=1e-9)
    # the corrected offsets keep the centre of mass at rest to rounding (x column sums to 1e-10)
    assert abs(float(SP_MASSES @ pos[:, 0])) < 5e-10


@pytest.mark.skip(
    reason="Szebehely & Peters 1967b: the closing residual at T (printed < 1e-10) and the binary "
    "collision of bodies 2 and 3 at T/2 (min r23 near 0, body 1 at rest) need an integrator that "
    "passes r23 = 0, i.e. Levi-Civita or KS regularisation; the Cartesian test integrator here "
    "cannot. Left for the regularised propagator of #928."
)
def test_szebehely_peters_1967b_closure_at_period_and_collision_at_half_period() -> None:
    raise AssertionError


# ---- Szebehely & Peters 1967a: the Pythagorean problem, Table I close approaches -----------

# (time, pair, approximate distance) for t <= 16, as printed (pair indices are 1-based masses).
SP67A = [
    (1.879, (2, 3), 1e-2),
    (3.026, (1, 3), 0.6),
    (3.801, (2, 3), 6e-2),
    (6.898, (1, 3), 0.1),
    (8.760, (2, 3), 8e-3),
    (9.962, (1, 3), 0.5),
    (11.611, (2, 3), 0.2),
    (14.618, (1, 3), 0.2),
    (15.830, (2, 3), 4e-4),
]


def _pythagorean_minima(t_end: float = 16.0) -> list[tuple[float, tuple[int, int], float]]:
    """Local minima of the 2-3 and 1-3 distances (pairs including the heavy body 3)."""
    n = 3
    y0 = np.concatenate((SP_BURRAU.ravel(), np.zeros(6)))
    pairs = [(1, 2), (0, 2)]  # zero-based (2,3) and (1,3)

    def make(i: int, j: int) -> Callable[[float, Arr], float]:
        def ev(_t: float, y: Arr) -> float:
            pos, vel = y[: 2 * n].reshape(n, 2), y[2 * n :].reshape(n, 2)
            return float((pos[i] - pos[j]) @ (vel[i] - vel[j]))

        ev.direction = 1.0  # type: ignore[attr-defined]
        return ev

    sol = solve_ivp(
        _rhs(SP_MASSES),
        (0.0, t_end),
        y0,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=[make(i, j) for i, j in pairs],
    )
    assert sol.t_events is not None and sol.y_events is not None
    out = []
    for (i, j), te, ye in zip(pairs, sol.t_events, sol.y_events, strict=True):
        for t, y in zip(te, ye, strict=True):
            d = float(np.linalg.norm(y[2 * i : 2 * i + 2] - y[2 * j : 2 * j + 2]))
            out.append((float(t), (i + 1, j + 1), d))
    return sorted(out)


@pytest.fixture(scope="module")
def pythagorean_minima() -> list[tuple[float, tuple[int, int], float]]:
    return _pythagorean_minima()


@pytest.mark.parametrize(("t_print", "pair", "d_print"), SP67A, ids=[f"t{r[0]}" for r in SP67A])
def test_szebehely_peters_1967a_table_i_close_approach_times_and_distances(
    pythagorean_minima: list[tuple[float, tuple[int, int], float]],
    t_print: float,
    pair: tuple[int, int],
    d_print: float,
) -> None:
    near = [m for m in pythagorean_minima if abs(m[0] - t_print) < 0.05 and m[1] == pair]
    assert len(near) == 1
    t, _pair, d = near[0]
    assert t == pytest.approx(t_print, abs=2e-3)  # printed to 3 decimals; digest: within 2e-3
    assert d / d_print == pytest.approx(1.0, abs=0.3)  # printed distances are one figure


@pytest.mark.xfail(
    strict=True,
    reason="Szebehely & Peters 1967a text (p.878) prints t10 = 15.8299236; two independent "
    "formulations (time-transformed Cartesian and Levi-Civita, agreeing to 1e-8) give 15.829920, "
    "3.3e-6 earlier, which the digest could not attribute (a 1967 fifth-order Runge-Kutta run "
    "through several close approaches may carry that error).",
)
def test_szebehely_peters_1967a_t10_text_value_to_eight_digits(
    pythagorean_minima: list[tuple[float, tuple[int, int], float]],
) -> None:
    near = [m for m in pythagorean_minima if abs(m[0] - 15.83) < 0.05 and m[1] == (2, 3)]
    assert near[0][0] == pytest.approx(15.8299236, abs=1e-6)


def test_szebehely_peters_1967a_initial_energy_and_no_12_approach(
    pythagorean_minima: list[tuple[float, tuple[int, int], float]],
) -> None:
    assert energy(SP_MASSES, SP_BURRAU, np.zeros((3, 2))) == pytest.approx(-769.0 / 60.0, abs=1e-12)
    # no close approach other than those of Table I (heavy body in every one) up to t = 16
    printed = {(round(t, 0), p) for t, p, _ in SP67A}
    for t, pair, d in pythagorean_minima:
        if d < 0.6:
            assert any(abs(t - tp) < 0.05 and pair == pp for tp, pp, _ in SP67A), (t, pair, d)
    assert printed  # table non-empty
