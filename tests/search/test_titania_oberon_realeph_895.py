"""#895: the Uranus-centred force model, its integrator, and its two positive controls.

Pre-registration: ``docs/notes/2026-10-04-895-titania-oberon-realeph.md`` section 1.2.

* P2 (no kernel needed, always runs): the planar `#890` four-body model rebuilt in this module's
  Uranus-centred non-rotating frame, through the same force function and integrator as the
  real-ephemeris arcs, reproduces the stored `#890` periodic orbit (half-cycle symmetry, the two
  flyby altitudes, one-cycle return), and a 17 km change of the model (Oberon about Uranus
  instead of the barycentre) destroys it. This is also the regression test on the `#890` orbit
  from code that does not share that build's right-hand side.
* P1 (needs the URA111 kernel, skips without it): each moon propagated as a test body in the
  force model stays within the pre-registered bounds of URA111, and the control fails without
  the zonal field.
* P3: gradient against finite differences, zonal closed forms, integrator against scipy, STM
  against finite differences, header constants against the kernel, interpolation against SPICE.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.search.titania_oberon_realeph_895 as m

REFINE = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "found"
    / "890_titania_oberon_candidate"
    / "refine.json"
)
KERNEL = pytest.mark.skipif(not m.kernels_present(), reason="URA111 SPICE kernel not installed")


@pytest.fixture(scope="module")
def c890() -> m.Model890:
    return m.registry_890()


@pytest.fixture(scope="module")
def s890() -> np.ndarray:
    return np.asarray(json.loads(REFINE.read_text())["refined_state"], dtype=np.float64)


# ------------------------------------------------------------------------------------------- #
# P3: physics pieces
# ------------------------------------------------------------------------------------------- #


def test_model890_constants_match_the_890_module(c890: m.Model890) -> None:
    from cyclerfinder.search.two_moon_periodic_890 import TwoMoonModel

    tm = TwoMoonModel()
    assert c890.gm_planet == pytest.approx(tm.gm_planet, rel=1e-15)
    assert c890.n_t == pytest.approx(tm.n_base, rel=1e-15)
    assert c890.mu == pytest.approx(tm.system.mu, rel=1e-14)
    assert c890.omega_o == pytest.approx(tm.system.omega_gan, rel=1e-13)
    assert c890.cycle_s == pytest.approx(5.0 * tm.forcing_period * tm.time_unit_s, rel=1e-13)


def _zonal_model(j2: float, j4: float, pole: np.ndarray) -> m.ForceModel:
    fm = m.model_890(m.registry_890())
    fm.gm[:] = 0.0
    fm.central = m.central_vector(m.GM_URANUS, j2, j4, pole)
    return fm


def test_zonal_closed_forms_on_axis_and_in_equator() -> None:
    pole = m.pole_vector(30.0, 50.0)
    fm2 = _zonal_model(m.J2_URA111, 0.0, pole)
    fm4 = _zonal_model(0.0, m.J4_URA111, pole)
    fm0 = _zonal_model(0.0, 0.0, pole)
    r = 100000.0
    gm, j2, j4, rr = m.GM_URANUS, m.J2_URA111, m.J4_URA111, m.R_REF_KM
    eq = np.cross(pole, [1.0, 0.0, 0.0])
    eq /= np.linalg.norm(eq)
    cases = ((r * pole, 1.0, 1.0, 3.0, 10.0), (r * eq, -0.5, 3.0 / 8.0, 0.0, 0.0))
    for x, p2, p4, dp2, dp4 in cases:
        a2 = fm2.accel(0.0, x) - fm0.accel(0.0, x)
        a4 = fm4.accel(0.0, x) - fm0.accel(0.0, x)
        u = float(np.dot(x, pole)) / r
        # a_n = GM J_n R^n / r^(n+2) [((n+1) P_n + u P_n') rhat - P_n' pole]
        e2 = gm * j2 * rr**2 / r**4 * ((3 * p2 + u * dp2) * x / r - dp2 * pole)
        e4 = gm * j4 * rr**4 / r**6 * ((5 * p4 + u * dp4) * x / r - dp4 * pole)
        np.testing.assert_allclose(a2, e2, rtol=1e-9, atol=1e-22)
        np.testing.assert_allclose(a4, e4, rtol=1e-9, atol=1e-24)
    # textbook J2 values: on the axis +3 GM J2 R^2/r^4 outward, in the equator -1.5 inward
    a_ax = fm2.accel(0.0, r * pole) - fm0.accel(0.0, r * pole)
    assert float(np.dot(a_ax, pole)) == pytest.approx(3 * gm * j2 * rr**2 / r**4, rel=1e-12)
    a_eq = fm2.accel(0.0, r * eq) - fm0.accel(0.0, r * eq)
    assert float(np.dot(a_eq, eq)) == pytest.approx(-1.5 * gm * j2 * rr**2 / r**4, rel=1e-12)


def test_gravity_gradient_matches_finite_differences(c890: m.Model890) -> None:
    fm = m.model_890(c890)
    fm.central = m.central_vector(c890.gm_planet, m.J2_URA111, m.J4_URA111, m.POLE_IAU)
    rng = np.random.default_rng(895)
    pts = [np.array([c890.a_t + 2500.0, 300.0, -150.0]), np.array([40000.0, -9000.0, 20000.0])]
    pts += [rng.normal(size=3) * 3e5 for _ in range(4)]
    for t, x in zip((0.0, 1e5, 3e5, 7e5, 1e6, 2e6), pts, strict=True):
        _, g = fm.accel_grad(t, x)
        gf = np.zeros((3, 3))
        h = 1e-4 * max(1.0, float(np.linalg.norm(x)) * 1e-5)
        for j in range(3):
            e = np.zeros(3)
            e[j] = h
            gf[:, j] = (fm.accel(t, x + e) - fm.accel(t, x - e)) / (2 * h)
        assert np.max(np.abs(g - gf)) <= 1e-6 * np.max(np.abs(g))


def test_integrator_matches_scipy_dop853_and_stm_matches_differences(
    c890: m.Model890, s890: np.ndarray
) -> None:
    fm = m.model_890(c890)
    fm.central = m.central_vector(c890.gm_planet, m.J2_URA111, 0.0, m.POLE_IAU)
    y0 = m.state_890_to_inertial(c890, s890)
    y0[2], y0[5] = 300.0, 0.002  # leave the plane so that every block of the STM is exercised
    t1 = 5.0 * m.DAY_S  # a node interval that contains the Titania flyby at t = 0 .. 0.3 d
    p = m.propagate(fm, -0.5 * m.DAY_S, t1, y0, stm=True, rtol=1e-13)
    sol = m.propagate_scipy(fm, -0.5 * m.DAY_S, t1, y0, rtol=1e-13)
    assert np.linalg.norm(p.state[:3] - sol.y[:3, -1]) < 1e-4  # km
    assert np.linalg.norm(p.state[3:] - sol.y[3:, -1]) < 1e-10  # km/s
    assert p.stm is not None
    fd = np.zeros((6, 6))
    steps = [1e-2] * 3 + [1e-8] * 3
    for j in range(6):
        e = np.zeros(6)
        e[j] = steps[j]
        yp = m.propagate(fm, -0.5 * m.DAY_S, t1, y0 + e, rtol=1e-13).state
        ym = m.propagate(fm, -0.5 * m.DAY_S, t1, y0 - e, rtol=1e-13).state
        fd[:, j] = (yp - ym) / (2 * steps[j])
    rel = np.abs(p.stm - fd) / (np.abs(fd).max(axis=0, keepdims=True))
    assert rel.max() < 1e-5


# ------------------------------------------------------------------------------------------- #
# P2: the #890 orbit in this frame (no kernel)
# ------------------------------------------------------------------------------------------- #


def _flyby_alt(fm: m.ForceModel, prop: m.Propagation, body: int, t0: float, t1: float) -> float:
    assert prop.rec_t is not None and prop.rec_y is not None
    tr = m.hermite_traj(prop.rec_t, prop.rec_y, fm)

    def bstate(t: float) -> tuple[np.ndarray, np.ndarray]:
        p, v = fm.body_states(t)
        return p[body], v[body]

    enc = m.find_minima(tr, bstate, t0, t1, dt=600.0, body=body, max_dist=20000.0)
    assert len(enc) == 1
    return enc[0].dist_km - m.BODY_RADIUS_KM[body]


def test_p2_890_orbit_reproduced_in_the_uranus_centred_frame(
    c890: m.Model890, s890: np.ndarray
) -> None:
    fm = m.model_890(c890)
    y0 = m.state_890_to_inertial(c890, s890)
    t_cyc = c890.cycle_s
    # (i) half-cycle symmetry conditions (the Oberon crossing of the x-axis, perpendicular)
    half = m.propagate(fm, 0.0, t_cyc / 2, y0, rtol=1e-13, atol_r=1e-10, atol_v=1e-16)
    rot = m.inertial_to_rot_890(c890, t_cyc / 2, half.state)
    assert abs(rot[1]) * c890.a_t <= 0.010  # km
    assert abs(rot[2]) * c890.a_t * c890.n_t <= 1e-6  # km/s
    # (ii) the two flyby altitudes
    full = m.propagate(fm, 0.0, t_cyc, y0, record=True, rtol=1e-13, atol_r=1e-10, atol_v=1e-16)
    d = m.DAY_S
    assert _flyby_alt(fm, full, m.I_OBERON, 55 * d, 68 * d) == pytest.approx(1364.2, abs=1.0)
    assert _flyby_alt(fm, full, m.I_TITANIA, t_cyc - 6 * d, t_cyc) == pytest.approx(1976.9, abs=1.0)
    # (iii) one-cycle return in the rotating frame
    end = m.inertial_to_rot_890(c890, t_cyc, full.state)
    assert np.linalg.norm(end[:2] - s890[:2]) * c890.a_t <= 1.0
    assert np.linalg.norm(end[2:] - s890[2:]) * c890.a_t * c890.n_t <= 1e-5


def test_p2_negative_control_oberon_about_uranus_destroys_the_orbit(
    c890: m.Model890, s890: np.ndarray
) -> None:
    good = m.model_890(c890)
    bad = m.model_890(c890, oberon_about_barycentre=False)
    y0 = m.state_890_to_inertial(c890, s890)
    t_cyc = c890.cycle_s
    ref = m.propagate(good, 0.0, 0.97 * t_cyc, y0, rtol=1e-13)
    alt = m.propagate(bad, 0.0, 0.97 * t_cyc, y0, rtol=1e-13)
    assert np.linalg.norm(alt.state[:3] - ref.state[:3]) > 100.0
    ret = m.propagate(bad, 0.0, t_cyc, y0, rtol=1e-13, max_steps=20000, allow_incomplete=True)
    if ret.complete:
        end = m.inertial_to_rot_890(c890, t_cyc, ret.state)
        assert np.linalg.norm(end[:2] - s890[:2]) * c890.a_t > 100.0


# ------------------------------------------------------------------------------------------- #
# Kernel-dependent: header constants, interpolation, P1
# ------------------------------------------------------------------------------------------- #


@KERNEL
def test_constants_are_the_ura111_header_values() -> None:
    txt = m.kernel_header_text()

    def num(pattern: str) -> float:
        hit = re.search(pattern, txt)
        assert hit is not None, pattern
        return float(hit.group(1).replace("D", "E"))

    for name, val in (
        ("Ariel", m.GM_ARIEL),
        ("Umbriel", m.GM_UMBRIEL),
        ("Titania", m.GM_TITANIA),
        ("Oberon", m.GM_OBERON),
        ("Miranda", m.GM_MIRANDA),
        ("Uranus", m.GM_URANUS),
    ):
        assert num(rf"{name}\s+7\d\d\s+([0-9.E+-]+)") == val
    assert num(r"J702\s+([0-9.E+-]+)") == m.J2_URA111
    assert num(r"J704\s+([0-9.E+-]+)") == m.J4_URA111
    assert num(r"RADIUS\s+([0-9.E+-]+)") == m.R_REF_KM
    assert num(r"GM10\s+([0-9.E+-]+)") == m.GM_SUN
    assert num(r"ZACPL7\s+([0-9.E+-]+)") == m.POLE_URA111_RA_DEC_DEG[0]
    assert num(r"ZDEPL7\s+([0-9.E+-]+)") == m.POLE_URA111_RA_DEC_DEG[1]


@pytest.fixture(scope="module")
def table_2030() -> m.EphemerisTable:
    t_ref = m.et_of("2030-01-12 00:00:00 TDB")
    return m.build_table(t_ref, -1 * m.DAY_S, 32 * m.DAY_S)


@KERNEL
def test_ephemeris_table_matches_spice_at_mid_knots(table_2030: m.EphemerisTable) -> None:
    fm = m.full_model(table_2030)
    bound_m = (5.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    for i in range(0, 4000, 37):
        t = table_2030.t0 + (i + 0.5) * table_2030.h
        pos, _ = fm.body_states(t)
        for k, naif in enumerate(m.NAIF_IDS):
            ref = m.spice_state(naif, table_2030.t_ref_et + t)
            assert np.linalg.norm(pos[k] - ref[:3]) * 1e3 < bound_m[k]


@KERNEL
def test_p1_moons_reproduce_ura111_and_the_control_discriminates(
    table_2030: m.EphemerisTable,
) -> None:
    t0 = 0.0
    for k in m.MOONS:
        (err,) = m.moon_positive_control(table_2030, k, t0, (30.0,))
        assert err <= 10.0, (m.BODY_NAMES[k], err)
    for k in (m.I_MIRANDA, m.I_ARIEL):
        for kw in ({"j2": 0.0, "j4": 0.0}, {"j2": m.J2_FORMER_REPO, "j4": 0.0}):
            (err,) = m.moon_positive_control(table_2030, k, t0, (30.0,), **kw)
            assert err > 10.0, (m.BODY_NAMES[k], kw, err)
    assert math.isfinite(err)
