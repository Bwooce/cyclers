"""#928 stage 2: the moon-centred KS propagator against exact and published controls.

Formulation: ``core/cr3bp_ks.py`` (Stiefel & Scheifele 1971 total-energy form (9,53)). Controls,
each with its source:

* pointwise: the KS right-hand side rebuilt into a physical acceleration equals
  ``core.cr3bp.cr3bp_eom``; ``V`` equals its direct form; ``h = C/2 - (1 - mu)^2/2 - (1 - mu)``.
* P = 0 Kepler motion against ``core.kepler.propagate`` (universal variables; exact for every
  conic), elliptic, hyperbolic, near-parabolic, and the e = 1 - 1e-8 ellipse through its
  pericentre; hyperbolic passes at q = 1e-2 to 1e-8 with the error flat in q.
* radial fall from rest (S&S section 5, eq. (5,30); digest section 14 item 3): collision at
  t = pi r0^(3/2) / (2 sqrt(2) K), then the reflected branch.
* Llibre 1982 (Celest. Mech. 26), mu = 0 rotating Kepler problem, Fig. 1a and 1c printed circular
  radii (roots of 4r^3 - C^2 r^2 + 2Cr - 1) and the zero-angular-momentum collision orbit with
  apocentre 2/C (digest ``2026-10-05-digest-llibre-1982-restricted-problem-small-mu.md`` sections
  3 and 8.1).
* Rodriguez del Rio 2021 thesis (digest
  ``2026-10-05-digest-rodriguez-del-rio-2021-thesis-part-a.md``): exactly four 1-EC orbits at
  mu = 0.1, C = 5 (pp. 42-43, T5); the printed Levi-Civita equations 2.26 integrated
  independently in this file agree with the KS run. The four angles are measured here, not
  printed (the thesis prints none), and are not asserted as goldens.
* agreement with ``cr3bp_eom`` away from the Moon; passes at q = 1e-2, 1e-3, 1e-4 against a tight
  Cartesian DOP853 run; Jacobi constant and reversibility flat in q down to 1e-6.
* Burdet 1968 (ZAMP 19:345) L4 benchmark, p.354-358: a time-dependent perturbation with an exact
  solution (the second body's Kepler orbit rotated by 60 degrees); printed periods 6.252003 and
  5.89817010367 (digest ``2026-10-04-digest-burdet-1968-theory-kepler-motion-perturbed-two-body.md``
  section 3). Tests the swappable force term and the integrated energy h.

Measured on 2026-10-05 (rtol 1e-13, atol 1e-15) and recorded at each test; bounds are set from
the measurement with margin.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from cyclerfinder.core import kepler
from cyclerfinder.core.cr3bp import cr3bp_eom, jacobi_constant
from cyclerfinder.core.cr3bp_ks import (
    KSModel,
    MoonCentredCR3BP,
    integrate_ks,
    ks_jacobian,
    ks_monitors,
    ks_rhs,
    ks_state_from_physical,
    physical_from_ks_state,
    propagate_ks,
)
from cyclerfinder.core.ks import ks_matrix, ks_position

FloatArray = NDArray[np.float64]

MU_EM = 0.0121505856  # Earth-Moon test value (the S&S digest's)
RTOL, ATOL = 1e-13, 1e-15


def _l3(u: FloatArray) -> FloatArray:
    return ks_matrix(u)[:3]


# --------------------------------------------------------------------------------------------
# pointwise


def test_potential_direct_form_and_value_at_centre() -> None:
    """V (cancellation-free form) equals its direct form; V(0) = 0 and grad V(0) = 0."""
    m = MoonCentredCR3BP(MU_EM)
    a = 1.0 - MU_EM
    rng = np.random.default_rng(928)
    for _ in range(5):
        x = rng.normal(size=3) * 0.3
        r1 = float(np.linalg.norm(x + np.array([1.0, 0.0, 0.0])))
        direct = -0.5 * ((x[0] + a) ** 2 + x[1] ** 2) - a / r1 + 0.5 * a * a + a
        assert m.potential(0.0, x) == pytest.approx(direct, abs=2e-15)
        g_direct = -np.array([x[0] + a, x[1], 0.0]) + a * (x + np.array([1.0, 0.0, 0.0])) / r1**3
        assert np.allclose(m.potential_grad(0.0, x), g_direct, rtol=0, atol=2e-15)
    assert m.potential(0.0, np.zeros(3)) == 0.0
    assert np.array_equal(m.potential_grad(0.0, np.zeros(3)), np.zeros(3))
    # near the centre V is O(r^2) and computed without cancellation: tidal limit
    x = np.array([3e-7, -2e-7, 1e-7])
    tidal = -0.5 * (x[0] ** 2 + x[1] ** 2) - a * (1.5 * x[0] ** 2 - 0.5 * float(x @ x))
    assert m.potential(0.0, x) == pytest.approx(tidal, rel=1e-6)


def test_pointwise_acceleration_matches_cr3bp_eom() -> None:
    """KS right-hand side rebuilt into xddot equals cr3bp_eom; monitors vanish; h = C/2 - c."""
    m = MoonCentredCR3BP(MU_EM)
    rng = np.random.default_rng(1971)
    for k in range(6):
        st = np.concatenate([rng.normal(size=3) * 0.3 + [1 - MU_EM, 0, 0], rng.normal(size=3)])
        y = ks_state_from_physical(m, st, fibre_angle=0.5 * k)
        f = ks_rhs(0.0, y, m)
        u, w, wp = y[:4], y[4:8], f[4:8]
        r = float(u @ u)
        xp = 2.0 * _l3(u) @ w
        xpp = 2.0 * _l3(w) @ w + 2.0 * _l3(u) @ wp
        acc = (xpp / r - xp * (2.0 * float(u @ w)) / r**2) / r
        ref = cr3bp_eom(0.0, st, MU_EM)[3:]
        assert np.allclose(acc, ref, rtol=0, atol=4e-15 * max(1.0, float(np.abs(ref).max())))
        assert np.allclose(physical_from_ks_state(m, y), st, rtol=0, atol=2e-15)
        k_res, bil = ks_monitors(m, y)
        assert abs(k_res) < 1e-14 and abs(bil) < 1e-15
        c = jacobi_constant(st, MU_EM)
        assert y[8] == pytest.approx(0.5 * c - m.energy_constant(), abs=4e-15)
        assert f[8] == 0.0  # autonomous: h' = 0 exactly
        assert f[9] == pytest.approx(r, rel=1e-15)


def test_analytic_jacobian_matches_finite_differences() -> None:
    """ks_jacobian against central differences of ks_rhs (step 1e-6): 1e-8 relative."""
    rng = np.random.default_rng(10)
    for model in (MoonCentredCR3BP(MU_EM), MoonCentredCR3BP(0.3), KSModel(0.7)):
        for _ in range(3):
            st = np.concatenate([rng.normal(size=3) * 0.2 + model.centre, rng.normal(size=3)])
            y = ks_state_from_physical(model, st, fibre_angle=float(rng.normal()))
            jac = ks_jacobian(y, model)
            fd = np.zeros((10, 10))
            for j in range(10):
                e = np.zeros(10)
                e[j] = 1e-6
                fd[:, j] = (ks_rhs(0.0, y + e, model) - ks_rhs(0.0, y - e, model)) / 2e-6
            assert np.abs(jac - fd).max() < 1e-8 * max(1.0, float(np.abs(jac).max()))


# --------------------------------------------------------------------------------------------
# P = 0 controls


@pytest.mark.parametrize(
    ("k2", "r0", "v0", "dt"),
    [
        (1.0, [1.0, 0.0, 0.0], [0.0, 1.2, 0.1], 3.0),  # ellipse e ~ 0.45
        (1.0, [1.0, 0.2, 0.0], [0.3, 1.5, 0.2], 2.0),  # hyperbola
        (0.0121505856, [0.3, -0.1, 0.05], [-0.2, 0.25, 0.1], 4.0),  # hyperbola, lunar k2
        (1.0, [1.0, 0.0, 0.0], [0.0, math.sqrt(2.0) * (1 + 1e-9), 0.0], 5.0),  # near-parabolic
        (1.0, [1.0, 0.0, 0.0], [0.0, math.sqrt(2.0) * (1 - 1e-9), 0.0], 5.0),  # near-parabolic
    ],
)
def test_kepler_against_universal_variables(
    k2: float, r0: list[float], v0: list[float], dt: float
) -> None:
    """P = 0: KS equals core.kepler.propagate. Measured 3e-16 to 1.1e-14 relative."""
    r0a, v0a = np.array(r0), np.array(v0)
    arc = propagate_ks(KSModel(k2), np.concatenate([r0a, v0a]), dt, rtol=RTOL, atol=ATOL)
    rr, vv = kepler.propagate(r0a, v0a, dt, mu=k2)
    assert np.abs(arc.state[:3] - rr).max() < 5e-14 * np.linalg.norm(rr)
    assert np.abs(arc.state[3:] - vv).max() < 5e-14 * np.linalg.norm(vv)
    assert abs(arc.energy_drift) < 1e-14 and arc.bilinear_max < 1e-15


def test_near_collision_ellipse_through_pericentre() -> None:
    """e = 1 - 1e-8, a = 1: through the pericentre (q = 1e-8) and out; against Kepler.

    Measured 2.2e-13 in position and 3.9e-13 relative in velocity at the end (the pericentre
    speed is 1.4e4); closest approach 1.000000005e-8.
    """
    k2, a, e = 1.0, 1.0, 1.0 - 1e-8
    # start at apocentre, run half a period plus 0.1 (through the pericentre)
    ra = a * (1 + e)
    va = math.sqrt(k2 * (1 - e) / ra)
    st = np.array([-ra, 0.0, 0.0, 0.0, -va, 0.0])
    dt = math.pi * a**1.5 / math.sqrt(k2) + 0.1
    arc = propagate_ks(KSModel(k2), st, dt, rtol=RTOL, atol=ATOL)
    rr, vv = kepler.propagate(st[:3], st[3:], dt, mu=k2)
    assert arc.r_min == pytest.approx(a * (1 - e), rel=1e-6)
    assert np.abs(arc.state[:3] - rr).max() < 1e-12
    assert np.abs(arc.state[3:] - vv).max() < 1e-12 * np.linalg.norm(vv)


@pytest.mark.parametrize("q", [1e-2, 1e-4, 1e-6, 1e-8])
def test_hyperbolic_pass_error_flat_in_pericentre(q: float) -> None:
    """Lunar k2, v_inf = 0.3: from r ~ 0.15 through pericentre q and out, against Kepler from
    the same start. Measured below 1.3e-14 in position and 2.4e-14 relative in velocity for every
    q from 1e-2 to 1e-8 (300 to 312 function evaluations, about 25 steps, for all q)."""
    k2 = MU_EM
    vp = math.sqrt(0.09 + 2 * k2 / q)
    rp = np.array([q, 0.0, 0.0])
    vpv = np.array([0.0, vp * math.cos(0.3), vp * math.sin(0.3)])
    r0, v0 = kepler.propagate(rp, vpv, -0.5, mu=k2)
    arc = propagate_ks(KSModel(k2), np.concatenate([r0, v0]), 1.0, rtol=RTOL, atol=ATOL)
    rr, vv = kepler.propagate(r0, v0, 1.0, mu=k2)
    assert np.abs(arc.state[:3] - rr).max() < 2e-14
    assert np.abs(arc.state[3:] - vv).max() < 5e-14 * np.linalg.norm(vv)
    assert arc.r_min == pytest.approx(q, rel=1e-6)
    assert arc.nfev < 500


def test_radial_fall_from_rest() -> None:
    """S&S (5,30): from rest at r0 = 1, K = 1, collision at t = pi/(2 sqrt 2) = 1.11072073453959,
    with u passing smoothly through 0 (measured: r = 0 at the event u1 = 0, t to 1e-14); then the
    reflected branch returns to r0 at t = pi/sqrt(2) (measured 8e-14)."""
    m = KSModel(1.0)
    y0 = ks_state_from_physical(m, np.array([1.0, 0, 0, 0, 0, 0]))
    assert y0[8] == 1.0  # h = K^2/r0

    def collide(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[0])

    collide.terminal = True  # type: ignore[attr-defined]
    sol = integrate_ks(m, y0, (0.0, 10.0), rtol=RTOL, atol=ATOL, events=collide)
    assert sol.y_events is not None
    yc = sol.y_events[0][0]
    assert float(yc[:4] @ yc[:4]) < 1e-26
    assert yc[9] == pytest.approx(math.pi / (2 * math.sqrt(2)), abs=1e-13)
    sol2 = integrate_ks(m, y0, (0.0, 2.0 * float(sol.t_events[0][0])), rtol=RTOL, atol=ATOL)  # type: ignore[index]
    yb = sol2.y[:, -1]
    assert np.allclose(ks_position(yb[:4]), [1.0, 0.0, 0.0], atol=1e-12)
    assert yb[9] == pytest.approx(math.pi / math.sqrt(2), abs=1e-12)


# --------------------------------------------------------------------------------------------
# mu = 0 rotating-frame controls (Llibre 1982)

_LLIBRE_RADII = [
    # (C, printed radius, sense: -1 retrograde r_1, +1 direct r_2, r_3), Fig. 1a and 1c
    (3.25, 0.2367865, -1),
    (3.25, 0.5783759, 1),
    (3.25, 1.825462, 1),
    (3.1, 0.2445555, -1),
    (3.1, 0.7022517, 1),
    (3.1, 1.4556927, 1),
]


@pytest.mark.parametrize(("c", "radius", "sense"), _LLIBRE_RADII)
def test_llibre_circular_orbits(c: float, radius: float, sense: int) -> None:
    """A circular Kepler orbit of the printed radius has the printed Jacobi constant (to the
    7-digit rounding of the radius: measured 1e-8 to 5e-7) and stays circular in the rotating
    frame, turning at n - 1 (measured: radius to 6e-15, angle to 2.7e-12 after t = 3)."""
    m = MoonCentredCR3BP(1.0)  # mu_core = 1: the mu = 0 problem about the origin
    n = sense * radius**-1.5
    st = np.array([radius, 0.0, 0.0, 0.0, (n - 1.0) * radius, 0.0])
    assert jacobi_constant(st, 1.0) == pytest.approx(c, abs=1e-6)
    wrong = np.array([radius, 0.0, 0.0, 0.0, (-n - 1.0) * radius, 0.0])
    assert abs(jacobi_constant(wrong, 1.0) - c) > 0.1  # the sense is a real check
    arc = propagate_ks(m, st, 3.0, rtol=RTOL, atol=ATOL)
    assert float(np.linalg.norm(arc.state[:3])) == pytest.approx(radius, abs=1e-13)
    ang = math.atan2(arc.state[1], arc.state[0])
    assert abs(math.remainder(ang - (n - 1.0) * 3.0, 2 * math.pi)) < 1e-11


def test_llibre_radial_orbit_returns_to_two_over_c() -> None:
    """Zero sidereal angular momentum ejection at C = 3.25: apocentre 2/C, rotating-frame
    velocity -e_z x x there, re-collision after the Kepler period 2 pi C^(-3/2) along the
    ejection direction turned by -t (the frame turns by +t). Measured: apocentre 2.4e-15,
    times 3.6e-15, direction 1.2e-14 (the wrong Coriolis sign would give -2.1 rad)."""
    m = MoonCentredCR3BP(1.0)
    c, theta = 3.25, 0.4
    y0 = np.zeros(10)
    y0[4], y0[5] = math.sqrt(0.5) * math.cos(theta), math.sqrt(0.5) * math.sin(theta)
    y0[8] = 0.5 * c  # h = C/2 (the constant c is 0 at mu_core = 1)
    assert abs(ks_monitors(m, y0)[0]) < 1e-15

    def extremum(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[:4] @ y[4:8])

    extremum.terminal = True  # type: ignore[attr-defined]
    extremum.direction = -1.0  # type: ignore[attr-defined]
    sol = integrate_ks(m, y0, (0.0, 100.0), rtol=RTOL, atol=ATOL, events=extremum)
    assert sol.y_events is not None and sol.t_events is not None
    ya = sol.y_events[0][0]
    assert float(ya[:4] @ ya[:4]) == pytest.approx(2.0 / c, abs=1e-13)
    st = physical_from_ks_state(m, ya)
    assert np.allclose(st[3:], [st[1], -st[0], 0.0], atol=1e-13)
    assert ya[9] == pytest.approx(math.pi * c**-1.5, abs=1e-13)
    extremum.direction = 1.0  # type: ignore[attr-defined]
    sol2 = integrate_ks(
        m, ya, (float(sol.t_events[0][0]), 200.0), rtol=RTOL, atol=ATOL, events=extremum
    )
    assert sol2.y_events is not None
    yc = sol2.y_events[0][0]
    assert float(yc[:4] @ yc[:4]) < 1e-26
    t_c = float(yc[9])
    assert t_c == pytest.approx(2.0 * math.pi * c**-1.5, abs=1e-13)
    xd = ks_position(yc[4:8])  # direction of x at the collision is that of L(w) w
    ang = math.atan2(xd[1], xd[0])
    assert abs(math.remainder(ang - (2 * theta - t_c), 2 * math.pi)) < 1e-12
    assert abs(math.remainder(ang - (2 * theta + t_c), 2 * math.pi)) > 1.0


# --------------------------------------------------------------------------------------------
# Rodriguez del Rio 2021: 1-EC orbits at mu = 0.1, C = 5

MU_T, C_T = 0.1, 5.0  # thesis convention: mu_core = 1 - mu_T, C_core = C_T - mu_T (1 - mu_T)


def _ec_model() -> tuple[MoonCentredCR3BP, float]:
    m = MoonCentredCR3BP(1.0 - MU_T)
    return m, 0.5 * (C_T - MU_T * (1.0 - MU_T)) - m.energy_constant()


def _first_minimum(theta: float, rtol: float) -> FloatArray:
    """Eject at Levi-Civita angle theta (thesis eq. 3.9), return the state at the first
    minimum of r (thesis section Sigma_m = {u.u' = 0, increasing})."""
    m, h = _ec_model()
    y0 = np.zeros(10)
    w = math.sqrt(0.5 * m.k2)  # thesis |u'_T| = 2 sqrt(2(1 - mu_T)); ours = that / 4 (dt = r ds)
    y0[4], y0[5], y0[8] = w * math.cos(theta), w * math.sin(theta), h

    def ext(s: float, y: FloatArray, model: KSModel) -> float:
        return float(y[:4] @ y[4:8])

    ext.terminal = True  # type: ignore[attr-defined]
    ext.direction = -1.0  # type: ignore[attr-defined]
    s1 = integrate_ks(m, y0, (0.0, 50.0), rtol=rtol, atol=ATOL, events=ext)
    assert s1.y_events is not None and s1.t_events is not None
    ext.direction = 1.0  # type: ignore[attr-defined]
    s2 = integrate_ks(
        m, s1.y_events[0][0], (float(s1.t_events[0][0]), 100.0), rtol=rtol, atol=ATOL, events=ext
    )
    assert s2.y_events is not None
    return np.asarray(s2.y_events[0][0], dtype=np.float64)


def _m1(theta: float, rtol: float = RTOL) -> float:
    y = _first_minimum(theta, rtol)
    return float(y[0] * y[5] - y[1] * y[4])  # thesis eq. 3.12, M_LC = u v' - v u'


def test_rodriguez_del_rio_four_one_ec_orbits() -> None:
    """Exactly four zeros of M_1 on [0, pi) (thesis T5, pp. 42-43; M_1 has period pi because the
    Levi-Civita plane double-covers), each a true ejection-collision orbit: |u| at the minimum is
    below 1e-12 (measured 3e-17 to 1.8e-16, i.e. r ~ 1e-32). Measured angles theta/pi =
    0.0446757, 0.3194493, 0.5455145, 0.7680409 (not printed in the thesis; the digest agent's
    independent Levi-Civita run found the same to six digits)."""
    n = 24
    grid = [math.pi * k / n for k in range(n)]
    vals = [_m1(t, 1e-12) for t in grid]
    roots = []
    for k in range(n):
        if vals[k] * vals[(k + 1) % n] < 0:
            a, b = grid[k], grid[k] + math.pi / n
            roots.append(brentq(_m1, a, b, xtol=1e-15))
    assert len(roots) == 4
    for theta in roots:
        y = _first_minimum(theta, RTOL)
        assert math.sqrt(float(y[:4] @ y[:4])) < 1e-12
        assert 0.5 < y[9] < 0.6  # one excursion: physical time 0.55 to 0.57 measured


def test_ks_matches_printed_levi_civita_equations() -> None:
    """Rodriguez del Rio eq. 2.26 (p.22) integrated as printed (a = 4, dt = 4 r ds_T, thesis
    frame = core frame with mu_core = 1 - mu_T) equals the KS run at s = 4 s_T: measured 1.1e-15."""

    def lc_rhs(s: float, z: FloatArray) -> FloatArray:
        u, v, up, vp = z[:4]
        mu, c = MU_T, C_T
        rho = u * u + v * v
        r2 = math.sqrt((1 + u * u - v * v) ** 2 + 4 * u * u * v * v)
        fu = (
            4 * mu * u
            + 16 * mu * u**3
            + 12 * rho**2 * u
            + 8 * mu * u / r2
            - 8 * mu * u * rho * (rho + 1) / r2**3
            - 4 * c * u
        )
        fv = (
            4 * mu * v
            - 16 * mu * v**3
            + 12 * rho**2 * v
            + 8 * mu * v / r2
            - 8 * mu * v * rho * (rho - 1) / r2**3
            - 4 * c * v
        )
        return np.array([up, vp, fu + 8 * rho * vp, fv - 8 * rho * up, 4 * rho])

    theta = 0.9  # any ejection angle
    lc0 = 2 * math.sqrt(2 * (1 - MU_T))  # thesis eq. 3.9
    z0 = np.array([0.0, 0.0, lc0 * math.cos(theta), lc0 * math.sin(theta), 0.0])
    s_t = 0.15
    lc = solve_ivp(lc_rhs, (0.0, s_t), z0, method="DOP853", rtol=RTOL, atol=ATOL).y[:, -1]
    m, h = _ec_model()
    y0 = np.zeros(10)
    y0[4:6], y0[8] = z0[2:4] / 4.0, h
    ks = integrate_ks(m, y0, (0.0, 4.0 * s_t), rtol=RTOL, atol=ATOL).y[:, -1]
    assert np.allclose(ks[:2], lc[:2], rtol=0, atol=1e-13)
    assert np.allclose(4.0 * ks[4:6], lc[2:4], rtol=0, atol=1e-13)
    assert ks[9] == pytest.approx(lc[4], abs=1e-13)
    assert ks[2] == 0.0 and ks[3] == 0.0 and ks[6] == 0.0 and ks[7] == 0.0
    # thesis Jacobi integral (2.28): u'^2 + v'^2 = 8 (u^2 + v^2) U, U from p.22 with C_T
    u, v = lc[0], lc[1]
    rho = u * u + v * v
    r2 = math.sqrt((1 + u * u - v * v) ** 2 + 4 * u * u * v * v)
    big_u = (
        0.5 * ((1 - MU_T) * rho**2 + MU_T * ((1 + u * u - v * v) ** 2 + 4 * u * u * v * v))
        + (1 - MU_T) / rho
        + MU_T / r2
        - 0.5 * C_T
    )
    assert lc[2] ** 2 + lc[3] ** 2 == pytest.approx(8 * rho * big_u, rel=1e-12)


# --------------------------------------------------------------------------------------------
# CR3BP against Cartesian


def test_agrees_with_cr3bp_eom_away_from_the_moon() -> None:
    """3D arc staying beyond 0.2 of the Moon: measured 2.1e-14 (state), 3.8e-14 (energy)."""
    m = MoonCentredCR3BP(MU_EM)
    st = np.array([0.8, 0.1, 0.05, 0.1, -0.2, 0.05])
    arc = propagate_ks(m, st, 1.0, rtol=RTOL, atol=ATOL)
    ref = solve_ivp(cr3bp_eom, (0, 1), st, method="DOP853", rtol=RTOL, atol=ATOL, args=(MU_EM,)).y[
        :, -1
    ]
    assert np.abs(arc.state - ref).max() < 2e-13
    assert arc.r_min > 0.19
    assert abs(arc.energy_drift) < 2e-13


def _pass_start(q: float) -> FloatArray:
    """A 3D state 0.06-0.07 from the Moon whose pass has pericentre about q (built by running
    the KS propagator back 0.08 from a pericentre state; the start is then just a state)."""
    m = MoonCentredCR3BP(MU_EM)
    vp = math.sqrt(0.05**2 + 2 * MU_EM / q)
    c7, s7, c4, s4 = math.cos(0.7), math.sin(0.7), math.cos(0.4), math.sin(0.4)
    peri = np.array([1 - MU_EM + q * c7, q * s7, 0.0, -vp * s7 * c4, vp * c7 * c4, vp * s4])
    return propagate_ks(m, peri, -0.08, rtol=RTOL, atol=ATOL).state


@pytest.mark.parametrize(("q", "bound"), [(1e-2, 2e-13), (1e-3, 5e-12), (1e-4, 1e-9)])
def test_pass_matches_tight_cartesian_reference(q: float, bound: float) -> None:
    """Through the pass and out (t = 0.16) against cr3bp_eom DOP853 at rtol 2.3e-14 (the
    tightest DOP853 accepts). Measured |KS - Cartesian| in two runs whose start states differ
    only in the last bits: 1.7e-14 and 4.8e-14, 9.2e-14 and 7.0e-13, 3.7e-11 and 1.2e-10; the
    Cartesian run's own Jacobi error moves with it (1.1e-14 and 1.7e-14, 2.1e-14 and 1.1e-12,
    5.8e-11 and 1.9e-10), so from q = 1e-3 down the reference is the limit, and the bounds are
    set by it."""
    m = MoonCentredCR3BP(MU_EM)
    s0 = _pass_start(q)
    arc = propagate_ks(m, s0, 0.16, rtol=RTOL, atol=ATOL)
    ref = solve_ivp(
        cr3bp_eom, (0, 0.16), s0, method="DOP853", rtol=2.3e-14, atol=2.3e-17, args=(MU_EM,)
    ).y[:, -1]
    assert arc.r_min == pytest.approx(q, rel=1e-6)
    assert np.abs(arc.state - ref).max() < bound
    c_ref = abs(jacobi_constant(ref, MU_EM) - jacobi_constant(s0, MU_EM))
    assert np.abs(arc.state - ref).max() < max(1e-13, 10 * c_ref)


@pytest.mark.parametrize("q", [1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
def test_pass_invariants_flat_in_pericentre(q: float) -> None:
    """KS Jacobi error, reversibility and tolerance self-consistency do not grow as q falls
    (measured for q = 1e-2..1e-6: |dC| 2e-14..7e-14, forward-back 1.3e-14..4e-14, rtol 1e-13
    against 1e-11 about 3e-12, 300-324 function evaluations). The Cartesian run's |dC| at rtol
    1e-13 grows from 7e-14 to 4.5e-6 over the same range (not asserted here; see the #929
    comparison)."""
    m = MoonCentredCR3BP(MU_EM)
    s0 = _pass_start(q)
    c0 = jacobi_constant(s0, MU_EM)
    arc = propagate_ks(m, s0, 0.16, rtol=RTOL, atol=ATOL)
    assert arc.r_min == pytest.approx(q, rel=1e-6)
    assert abs(jacobi_constant(arc.state, MU_EM) - c0) < 3e-13
    back = propagate_ks(m, arc.state, -0.16, rtol=RTOL, atol=ATOL)
    assert np.abs(back.state - s0).max() < 3e-13
    loose = propagate_ks(m, s0, 0.16, rtol=1e-11, atol=1e-13)
    assert np.abs(loose.state - arc.state).max() < 2e-11
    assert arc.hamiltonian_max < 1e-13 and arc.bilinear_max < 1e-15
    assert arc.nfev < 600


def test_planar_data_stay_planar_and_fibre_angle_is_invisible() -> None:
    m = MoonCentredCR3BP(MU_EM)
    st = np.array([1 - MU_EM + 0.003, 0.002, 0.0, 0.3, -2.6, 0.0])
    arc = propagate_ks(m, st, 0.1, rtol=RTOL, atol=ATOL)
    assert arc.y[2] == 0.0 and arc.y[3] == 0.0 and arc.y[6] == 0.0 and arc.y[7] == 0.0
    assert arc.state[2] == 0.0 and arc.state[5] == 0.0
    turned = propagate_ks(m, st, 0.1, rtol=RTOL, atol=ATOL, fibre_angle=1.1)
    assert np.abs(turned.state - arc.state).max() < 1e-12


# --------------------------------------------------------------------------------------------
# Burdet 1968 L4 benchmark (time-dependent perturbation, swapped force term)

_M2, _MU2 = 0.01, 1.01  # mass ratio m2/m1 and G (m1 + m2), units G m1 = 1 (p.354, p.358)
_R60 = np.array([[0.5, -math.sqrt(3) / 2, 0.0], [math.sqrt(3) / 2, 0.5, 0.0], [0.0, 0.0, 1.0]])


class _L4Model(KSModel):
    """Particle about m1 (K^2 = G m1 = 1) perturbed by m2 on a Kepler orbit, Burdet eq. (300):
    P = -m2 ((x - x2)/rho^3 + x2/r2^3) (direct plus indirect term)."""

    def __init__(self, r2: FloatArray, v2: FloatArray) -> None:
        super().__init__(1.0)
        self.autonomous = False
        self.r2, self.v2 = r2, v2

    def force(self, t: float, x: FloatArray) -> FloatArray:
        x2 = self.r2 if t == 0.0 else kepler.propagate(self.r2, self.v2, t, mu=_MU2)[0]
        d = x - x2
        return np.asarray(
            -_M2 * (d / np.linalg.norm(d) ** 3 + x2 / np.linalg.norm(x2) ** 3), dtype=np.float64
        )


@pytest.mark.parametrize(
    ("a", "e", "printed_period", "revs", "bound"),
    [
        (1.0, 0.0, "6.252003", 2, 1e-12),  # measured 1.2e-13 (position), 1.0e-13 (velocity)
        (101 / 105, 59 / 101, "5.89817010367", 1, 1e-12),  # measured 3.0e-13, 1.0e-12 rel
        (1.0, 0.95, "6.252003", 1, 3e-11),  # measured 1.4e-12, 1.4e-11 relative velocity
    ],
)
def test_burdet_l4_benchmark(
    a: float, e: float, printed_period: str, revs: int, bound: float
) -> None:
    """The L4 particle follows m2's orbit turned by 60 degrees (exact solution; initial state at
    m2's pericentre, derived, not printed). The period 2 pi a^(3/2)/sqrt(1.01) of the exact
    fractions matches the printed one to the printed digits (5.4e-8 at 7 digits; 7.3e-10 for the
    ellipse, whose printed digits carry the rounding of e noted in the digest)."""
    q = a * (1 - e)
    r2 = np.array([q, 0.0, 0.0])
    v2 = np.array([0.0, math.sqrt(_MU2 * (1 + e) / q), 0.0])
    period = 2 * math.pi * a**1.5 / math.sqrt(_MU2)
    digits = len(printed_period.split(".")[1])
    assert abs(period - float(printed_period)) < 10.0 ** (-digits) + 1e-9
    st = np.concatenate([_R60 @ r2, _R60 @ v2])
    arc = propagate_ks(_L4Model(r2, v2), st, revs * period, rtol=RTOL, atol=ATOL)
    xe, ve = kepler.propagate(r2, v2, revs * period, mu=_MU2)
    assert np.abs(arc.state[:3] - _R60 @ xe).max() < bound
    assert np.abs(arc.state[3:] - _R60 @ ve).max() < bound * np.linalg.norm(ve)
    assert abs(arc.energy_drift) < 1e-12  # integrated h against algebraic h
    assert arc.hamiltonian_max < 1e-13
