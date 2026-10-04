"""Identity tests for the #890 two-moon CCR4BP tools (frames, two-body limit, symmetry)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.search import two_moon_periodic_890 as tm
from cyclerfinder.verify.turn_gate_closures import symmetric_closure

TOF = 61.58023222936292  # 2.5 Titania-Oberon synodic periods (#888 candidates.jsonl)


@pytest.fixture(scope="module")
def geom() -> tm.ClosureGeometry:
    c = symmetric_closure(
        "Uranus",
        "Titania",
        "Oberon",
        tof_days=TOF,
        n_rev=(5, 5),
        rel_offset_deg=0.0,
        branches=("high", "high", "high"),
    )
    assert c is not None
    return tm.closure_geometry(c)


def test_frame_round_trip() -> None:
    m = tm.TwoMoonModel(lam=1.0)
    rng = np.random.default_rng(1)
    for _ in range(5):
        tau = float(rng.uniform(-20, 20))
        s = rng.normal(size=4)
        r, v = tm.inertial_from_rot(m, tau, s)
        assert np.allclose(tm.rot_from_inertial(m, tau, r, v), s, rtol=0, atol=1e-13)


def test_moons_fixed_and_circular_in_inertial_frame() -> None:
    """Titania sits on its circle at n_T tau; Oberon at angle (1 + omega) tau, circular speed."""
    m = tm.TwoMoonModel(lam=1.0)
    tau = 7.3
    r, _v = m.moon_inertial("Titania", tau)
    assert math.isclose(math.atan2(r[1], r[0]), tau - 2 * math.pi, abs_tol=1e-12)
    ro, vo = m.moon_inertial("Oberon", tau)
    rate = 1.0 + m.system.omega_gan
    assert math.isclose(float(np.linalg.norm(vo)) / float(np.linalg.norm(ro)), rate * m.n_base)
    assert abs(float(ro @ vo)) < 1e-9 * float(np.linalg.norm(ro) * np.linalg.norm(vo))


def test_two_body_limit_reproduces_kepler_over_one_leg(geom: tm.ClosureGeometry) -> None:
    """lam = 0: the CCR4BP flow, through the frame conversions, is the Kepler leg (1e-9)."""
    m0 = tm.TwoMoonModel(lam=0.0)
    k1 = tm.kepler_propagate(geom.r_dep, geom.v_dep, 86400.0, geom.gm_sys)  # off Titania's centre
    t1 = m0.tau(1.0)
    s = tm.rot_from_inertial(m0, t1, k1[:2], k1[2:])
    t_end = m0.tau(geom.tof_days)
    a = tm.propagate(m0, s, t1, t_end)
    r, v = tm.inertial_from_rot(m0, t_end, a.state)
    k = tm.kepler_propagate(k1[:2], k1[2:], (geom.tof_days - 1.0) * 86400.0, geom.gm_sys)
    assert np.max(np.abs(r - k[:2])) / m0.length_km < 1e-9
    assert np.max(np.abs(v - k[2:])) / m0.vel_unit_kms < 1e-9
    # and the leg ends on Oberon (the patched-conic closure geometry)
    ro, _ = m0.moon_inertial("Oberon", t_end)
    assert float(np.linalg.norm(r - ro)) < 1e-3


def test_forcing_period_matches_patched_conic_at_lam_zero() -> None:
    m0 = tm.TwoMoonModel(lam=0.0)
    assert math.isclose(m0.days(2.5 * m0.forcing_period), TOF, rel_tol=1e-12)


def test_perturber_on_axis_at_half_period() -> None:
    m = tm.TwoMoonModel(lam=1.0)
    th = m.system.omega_gan * 2.5 * m.forcing_period
    assert math.isclose(math.cos(th), -1.0, abs_tol=1e-12)
    assert abs(math.sin(th)) < 1e-12


def test_mirror_time_reversal_symmetry() -> None:
    """theta0 = 0: M s(-tau) is a solution; theta0 = 0.3 breaks it (control)."""
    s0 = np.array([1.1, 0.03, 0.02, -0.4])
    tau = 3.7
    for theta0, holds in ((0.0, True), (0.3, False)):
        m = tm.TwoMoonModel(lam=1.0, theta0=theta0)
        f = tm.propagate(m, s0, 0.0, tau).state
        b = tm.propagate(m, tm.mirror(f), -tau, 0.0).state
        err = float(np.max(np.abs(b - tm.mirror(s0))))
        assert (err < 1e-10) is holds


def test_flyby_hyperbola_turns_v_in_onto_v_out() -> None:
    """Two-body propagation through the periapsis reproduces the requested asymptotes."""
    gm = 226.9
    v_in = np.array([0.1, 0.25])
    ang = math.radians(50.0)
    c, s = math.cos(ang), math.sin(ang)
    v_out = np.array([c * v_in[0] - s * v_in[1], s * v_in[0] + c * v_in[1]])
    h = tm.flyby_hyperbola(v_in, v_out, gm)
    a_in, a_out = tm.osculating_asymptotes(gm, h.rp_vec, h.vp_vec)
    assert np.allclose(a_in, v_in, atol=1e-12)
    assert np.allclose(a_out, v_out, atol=1e-12)
    osc = tm.osculating_flyby(gm, h.rp_vec, h.vp_vec)
    assert math.isclose(osc["turn_deg"], 50.0, rel_tol=1e-12)


def test_shooter_residual_zero_on_exact_trajectory() -> None:
    """Nodes sampled from one integration of a symmetric start give zero interior residual."""
    m = tm.TwoMoonModel(lam=1.0)
    sh = tm.SymmetricShooter(m, 2.0, 3)
    s0 = np.array([1.2, 0.0, 0.0, -0.35])
    mids = [tm.propagate(m, s0, 0.0, float(t)).state for t in sh.taus[1:-1]]
    s_end = tm.propagate(m, s0, 0.0, sh.T).state
    z = tm.SymmetricShooter.pack(s0, mids, s_end)
    res, jac = sh.evaluate(z)
    assert jac is not None
    assert np.max(np.abs(res[:-4])) < 1e-11
    # the last block measures only the dropped (y, vx) of the end state, propagated back
    assert np.max(np.abs(res[-4:])) > 0.0
    # Jacobian against finite differences on one column
    eps = 1e-7
    dz = np.zeros_like(z)
    dz[1] = eps
    rp, _ = sh.evaluate(z + dz, with_jac=False)
    rm, _ = sh.evaluate(z - dz, with_jac=False)
    assert np.allclose((rp - rm) / (2 * eps), jac[:, 1], atol=1e-6)


def test_rescale_flyby_nodes_scales_periapsis_offset(geom: tm.ClosureGeometry) -> None:
    """Identity at equal lam; periapsis offset proportional to lam otherwise."""
    m1 = tm.TwoMoonModel(lam=0.1)
    m2 = tm.TwoMoonModel(lam=0.2)
    _sh, z, _ = tm.symmetric_guess(m1, geom, 3)
    same = tm.rescale_flyby_nodes(z, m1, m1, 3)
    assert np.allclose(same, z, atol=1e-13)
    z2 = tm.rescale_flyby_nodes(z, m1, m2, 3)
    d1 = float(
        np.linalg.norm(tm.relative_inertial(m1, "Titania", 0.0, np.array([z[0], 0, 0, z[1]]))[0])
    )
    d2 = float(
        np.linalg.norm(tm.relative_inertial(m2, "Titania", 0.0, np.array([z2[0], 0, 0, z2[1]]))[0])
    )
    assert math.isclose(d2 / d1, 2.0, rel_tol=1e-9)
