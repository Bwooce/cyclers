"""CCR4BP core EOM / STM / propagator tests (#689, Stage B sub-task 1 of #686).

The Jupiter-Europa-Ganymede Concentric Circular Restricted 4-Body Problem
(:mod:`cyclerfinder.core.ccr4bp`), the structural analogue of
:mod:`cyclerfinder.core.bcr4bp`. Gates:

  1. mu_Ganymede -> 0 reduction (STRUCTURAL, the EASY limit): at mu_gan = 0 the
     CCR4BP EOM / STM / propagator match the Jupiter-Europa CR3BP
     (:mod:`cyclerfinder.core.cr3bp`) at the Europa mass ratio, to
     floating-point precision -- exactly mirroring bcr4bp's mu_S -> 0 -> CR3BP
     test. Trivial because Ganymede's direct+indirect terms vanish pointwise.

  2. mu_Europa -> 0 reduction (STRUCTURAL, the SUBTLE limit): at mu = 0 the
     CCR4BP EOM reduces to the Jupiter-Ganymede CR3BP -- but ONLY after a
     frame/units reinterpretation (rotate to the Ganymede-fixed frame, shift the
     origin to the Jupiter-Ganymede barycentre, rescale length by 1/a_gan, and
     rescale time so the frame rate is 1). See :func:`_europa_frame_to_ganymede`.
     The transformation is validated independently by two kinematic sub-checks
     (Jupiter and Ganymede map to the CR3BP fixed points at rest) before the
     dynamical check, so a pass cannot be a coincidence.

  3. STM finite-difference consistency (full model): the variational 6x6 STM
     agrees with a central-difference Jacobian over a moderate horizon. This is
     the independent correctness check for the STM in the full CCR4BP (the
     hard mu -> 0 reduction is validated at the EOM level only; conjugating the
     STM through the time-rescaling transformation is out of scope here).

  4. Indirect-term cancellation at the Jupiter-Europa barycentre (r = 0): the
     net synodic-frame Ganymede acceleration is zero, mirroring bcr4bp Gate 2.

  5. Sourced-constant / commensurability cross-check against #688's own
     ``genome/composed_moon_map.py`` (same JPL SSD registry): radii, mass
     ratios, and the ~2:1 Laplace period ratio must agree.

All tests assert STRUCTURE (reduction limits, FD<->STM agreement) or sourced
parameter VALUES; no value computed by our own code sits on the EXPECTED side of
a structural assertion (per ``feedback_golden_tests_sourced_only``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.composed_moon_map import moon_config

# Sample states away from the primaries (no close-approach singularities).
_SAMPLE_STATE_PLANAR = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)
_SAMPLE_STATE_3D = np.array([0.6, 0.2, 0.04, 0.05, 0.35, -0.02], dtype=np.float64)


# ---------------------------------------------------------------------------
# Gate 1: mu_Ganymede -> 0 (the easy, pointwise reduction to Jupiter-Europa CR3BP).
# ---------------------------------------------------------------------------


def _system_mugan0() -> ccr4bp.CCR4BPSystem:
    d = ccr4bp.jupiter_europa_ganymede_default()
    return ccr4bp.CCR4BPSystem(mu=d.mu, mu_gan=0.0, a_gan=d.a_gan, omega_gan=d.omega_gan)


def test_mugan0_limit_eom_planar() -> None:
    """At mu_gan = 0 the CCR4BP planar EOM matches Jupiter-Europa CR3BP exactly."""
    sys_c = _system_mugan0()
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys_c)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys_c.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14), (
        f"CCR4BP@mu_gan=0 deviates from CR3BP: max|delta|={np.max(np.abs(rhs_c - rhs_3)):.3e}"
    )


def test_mugan0_limit_eom_3d() -> None:
    """At mu_gan = 0 the CCR4BP 3D EOM matches CR3BP exactly (time is irrelevant)."""
    sys_c = _system_mugan0()
    rhs_c = ccr4bp.ccr4bp_eom(2.3, _SAMPLE_STATE_3D, sys_c)  # nonzero time: Ganymede phase ignored
    rhs_3 = cr3bp.cr3bp_eom(2.3, _SAMPLE_STATE_3D, sys_c.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_mugan0_limit_stm() -> None:
    """At mu_gan = 0 the CCR4BP STM EOM matches the CR3BP STM EOM exactly."""
    sys_c = _system_mugan0()
    y42 = np.concatenate([_SAMPLE_STATE_3D, np.eye(6).reshape(36)])
    rhs_c = ccr4bp.ccr4bp_stm_eom(0.5, y42, sys_c)
    rhs_3 = cr3bp.cr3bp_stm_eom(0.5, y42, sys_c.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-13), (
        f"CCR4BP STM@mu_gan=0 deviates from CR3BP STM: "
        f"max|delta|={np.max(np.abs(rhs_c - rhs_3)):.3e}"
    )


def test_mugan0_limit_propagator() -> None:
    """At mu_gan = 0 the CCR4BP propagator agrees with CR3BP over a finite arc."""
    sys_c = _system_mugan0()
    sys_3 = cr3bp.CR3BPSystem(
        mu=sys_c.mu, primary="jupiter", secondary="europa", l_km=671100.0, t_s=1.0
    )
    arc_c = ccr4bp.propagate_ccr4bp(sys_c, _SAMPLE_STATE_3D, 1.0)
    arc_3 = cr3bp.propagate(sys_3, _SAMPLE_STATE_3D, 1.0)
    assert np.allclose(arc_c.state_f, arc_3.state_f, rtol=1e-11, atol=1e-11), (
        f"CCR4BP@mu_gan=0 propagation deviates: max|delta|="
        f"{np.max(np.abs(arc_c.state_f - arc_3.state_f)):.3e}"
    )


# ---------------------------------------------------------------------------
# Gate 2: mu_Europa -> 0 (the subtle reduction to Jupiter-Ganymede CR3BP).
# ---------------------------------------------------------------------------


def _europa_frame_to_ganymede(
    t: float,
    pos: np.ndarray,
    vel: np.ndarray,
    accel: np.ndarray,
    system: ccr4bp.CCR4BPSystem,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map an Europa-synodic-frame (pos, vel, accel) to Jupiter-Ganymede CR3BP units.

    Steps (see the module docstring): (1) rotate by -theta_gan(t) so Ganymede
    lands on +x (a steadily-rotating frame at rate omega_gan, so velocity and
    acceleration pick up the standard rotating-frame terms); (2) shift the origin
    from Jupiter to the Jupiter-Ganymede barycentre (constant offset mu3*a_gan);
    (3) rescale length by 1/a_gan; (4) rescale time by k = 1 + omega_gan so the
    reduced frame rate is 1. In-plane rotation leaves z scaling only.

    Returns ``(X, V, A)`` in Jupiter-Ganymede CR3BP nondimensional coordinates.
    """
    omega = system.omega_gan
    k = 1.0 + omega
    a_gan = system.a_gan
    mu3 = system.mu3_reduction
    theta = system.theta_gan0 + omega * t
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, s], [-s, c]], dtype=np.float64)  # R(-theta)
    jmat = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=np.float64)

    p2 = pos[:2]
    v2 = vel[:2]
    a2 = accel[:2]
    p_rot = rot @ p2
    v_rot = rot @ v2
    a_rot = rot @ a2

    # Rotating-frame (rate omega) velocity / acceleration transport terms.
    v_prime = v_rot - omega * (jmat @ p_rot)
    a_prime = a_rot - 2.0 * omega * (jmat @ v_rot) - omega * omega * p_rot

    # Origin shift to the J-G barycentre (constant; does not affect v/a), then
    # length scale 1/a_gan and time scale k.
    x_xy = (p_rot - np.array([mu3 * a_gan, 0.0])) / a_gan
    v_xy = v_prime / (a_gan * k)
    a_xy = a_prime / (a_gan * k * k)

    x_out = np.array([x_xy[0], x_xy[1], pos[2] / a_gan], dtype=np.float64)
    v_out = np.array([v_xy[0], v_xy[1], vel[2] / (a_gan * k)], dtype=np.float64)
    a_out = np.array([a_xy[0], a_xy[1], accel[2] / (a_gan * k * k)], dtype=np.float64)
    return x_out, v_out, a_out


def _reduced_system_mu0() -> ccr4bp.CCR4BPSystem:
    """A self-consistent mu = 0 CCR4BP: Europa massless, omega_gan set to the
    two-body value for the reduced Jupiter-Ganymede pair so the reduction is
    EXACT (the reduced frame rate matches the J-G CR3BP unit rate)."""
    d = ccr4bp.jupiter_europa_ganymede_default()
    omega = ccr4bp.two_body_synodic_rate(0.0, d.mu_gan, d.a_gan)
    return ccr4bp.CCR4BPSystem(mu=0.0, mu_gan=d.mu_gan, a_gan=d.a_gan, omega_gan=omega)


def test_mu0_transform_maps_primaries_to_cr3bp_fixed_points() -> None:
    """Kinematic sub-check: Jupiter and Ganymede map to the J-G CR3BP fixed
    points (-mu3, 0) and (1-mu3, 0) AT REST -- validates the transformation
    independently of any dynamics."""
    sys_c = _reduced_system_mu0()
    mu3 = sys_c.mu3_reduction
    for t in [0.0, 1.1, 3.7]:
        # Jupiter: at rest at the Europa-frame origin.
        xj, vj, _ = _europa_frame_to_ganymede(t, np.zeros(3), np.zeros(3), np.zeros(3), sys_c)
        assert np.allclose(xj, [-mu3, 0.0, 0.0], atol=1e-13)
        assert np.allclose(vj, 0.0, atol=1e-13)
        # Ganymede: on its synodic circle, moving with it.
        gx, gy, _ = ccr4bp._ganymede_position(t, sys_c)
        theta = sys_c.theta_gan0 + sys_c.omega_gan * t
        gpos = np.array([gx, gy, 0.0])
        gvel = sys_c.a_gan * sys_c.omega_gan * np.array([-math.sin(theta), math.cos(theta), 0.0])
        xg, vg, _ = _europa_frame_to_ganymede(t, gpos, gvel, np.zeros(3), sys_c)
        assert np.allclose(xg, [1.0 - mu3, 0.0, 0.0], atol=1e-13)
        assert np.allclose(vg, 0.0, atol=1e-12)


def test_mu0_limit_reduces_to_jupiter_ganymede_cr3bp() -> None:
    """At mu = 0, the CCR4BP EOM equals the Jupiter-Ganymede CR3BP EOM after the
    documented frame/units transformation, to ~1e-11 (exact up to solver-free
    floating-point)."""
    sys_c = _reduced_system_mu0()
    mu3 = sys_c.mu3_reduction
    # Several times (nonzero Ganymede phases) and both a planar and a 3D state.
    for state in (_SAMPLE_STATE_PLANAR, _SAMPLE_STATE_3D):
        for t in [0.0, 0.8, 2.9]:
            rhs = ccr4bp.ccr4bp_eom(t, state, sys_c)
            pos, vel, accel = state[:3], state[3:6], rhs[3:6]
            x_g, v_g, a_g = _europa_frame_to_ganymede(t, pos, vel, accel, sys_c)
            state_g = np.concatenate([x_g, v_g])
            cr = cr3bp.cr3bp_eom(0.0, state_g, mu3)
            # velocity part is the identity check; acceleration part is the physics.
            assert np.allclose(cr[:3], v_g, atol=1e-12)
            assert np.allclose(cr[3:6], a_g, rtol=1e-9, atol=1e-11), (
                f"mu=0 reduction failed at t={t}, state={state}: "
                f"max|delta|={np.max(np.abs(cr[3:6] - a_g)):.3e}"
            )


# ---------------------------------------------------------------------------
# Gate 3: STM finite-difference consistency (full model).
# ---------------------------------------------------------------------------


def test_stm_finite_difference_consistency() -> None:
    """Variational STM matches a central-finite-difference Jacobian (full CCR4BP)."""
    sys_c = ccr4bp.jupiter_europa_ganymede_default()
    state0 = _SAMPLE_STATE_3D.copy()
    t_horizon = 0.3

    arc = ccr4bp.propagate_ccr4bp(sys_c, state0, t_horizon, with_stm=True)
    assert arc.stm is not None
    stm = arc.stm

    eps = 1e-6
    fd_jac = np.zeros((6, 6), dtype=np.float64)
    for i in range(6):
        sp = state0.copy()
        sp[i] += eps
        sm = state0.copy()
        sm[i] -= eps
        arc_p = ccr4bp.propagate_ccr4bp(sys_c, sp, t_horizon)
        arc_m = ccr4bp.propagate_ccr4bp(sys_c, sm, t_horizon)
        fd_jac[:, i] = (arc_p.state_f - arc_m.state_f) / (2.0 * eps)

    scale = max(1.0, float(np.max(np.abs(stm))))
    rel_delta = float(np.max(np.abs(fd_jac - stm))) / scale
    assert rel_delta < 1e-6, (
        f"STM vs central-FD relative disagreement: {rel_delta:.3e} "
        f"(abs max|delta|={np.max(np.abs(fd_jac - stm)):.3e}, stm_scale={scale:.3e})"
    )


# ---------------------------------------------------------------------------
# Gate 4: Indirect-term cancellation at the Jupiter-Europa barycentre.
# ---------------------------------------------------------------------------


def test_indirect_term_cancels_at_je_barycenter() -> None:
    """At r = 0 the direct + indirect Ganymede contributions cancel exactly."""
    sys_c = ccr4bp.jupiter_europa_ganymede_default()
    for frac in [0.0, 0.3, 0.7]:
        t = frac * sys_c.ganymede_synodic_period
        ax, ay, az = ccr4bp._ganymede_acceleration(0.0, 0.0, 0.0, t, sys_c)
        assert abs(ax) < 1e-12 and abs(ay) < 1e-12 and abs(az) < 1e-12, (
            f"Indirect cancellation failed at t={t}: ({ax:.3e}, {ay:.3e}, {az:.3e})"
        )


# ---------------------------------------------------------------------------
# Gate 5: Sourced-constant / commensurability cross-check vs #688 composed map.
# ---------------------------------------------------------------------------


def test_constants_cross_check_against_composed_moon_map() -> None:
    """The CCR4BP default parameters trace to the SAME JPL SSD registry values
    #688's composed_moon_map uses -- radii, mass ratios, and ~2:1 period ratio."""
    sys_c = ccr4bp.jupiter_europa_ganymede_default()
    europa = moon_config("Europa")
    ganymede = moon_config("Ganymede")

    # Europa base mu computed identically to composed_moon_map (exact match).
    assert sys_c.mu == europa.mu
    # a_gan is the pure radius ratio -- independent of any mean-motion convention.
    assert sys_c.a_gan == pytest.approx(ganymede.sma_km / europa.sma_km, rel=0.0, abs=1e-15)
    # Reduced Ganymede mass ratio matches composed's Ganymede mu to ~1e-4: both
    # are 7.804e-5; they differ by rel ~2.5e-5 = GM_Europa/GM_Jupiter because the
    # CCR4BP mass unit is G(m1+m2) (carries GM_Europa in the denominator) whereas
    # composed's Ganymede mu uses GM_G/(GM_J+GM_G). Physically immaterial.
    assert sys_c.mu3_reduction == pytest.approx(ganymede.mu, rel=1e-4)

    # ~2:1 Laplace commensurability: the physical period ratio T_G/T_E ~ 2.014.
    # composed_moon_map uses the geometric mean motion; the CCR4BP two-body rate
    # differs by ~1.3e-5, so agree to 1e-3 (the physics, not the tiny formula gap).
    ratio_composed = europa.mean_motion_rad_s / ganymede.mean_motion_rad_s
    period_ratio_ccr4bp = (1.0 + sys_c.omega_gan) ** (-1.0)  # n2/n3 = T_G/T_E
    assert ratio_composed == pytest.approx(2.0144, abs=5e-4)
    assert period_ratio_ccr4bp == pytest.approx(ratio_composed, abs=1e-3)
    assert 2.0 < period_ratio_ccr4bp < 2.03  # near the 2:1 Laplace resonance


def test_ganymede_synodic_period_and_position() -> None:
    """Ganymede lies on a planar circle of radius a_gan and returns after one
    synodic period."""
    sys_c = ccr4bp.jupiter_europa_ganymede_default()
    for t in np.linspace(0.0, sys_c.ganymede_synodic_period, 11):
        gx, gy, gz = ccr4bp._ganymede_position(float(t), sys_c)
        assert gz == 0.0
        assert abs(math.hypot(gx, gy) - sys_c.a_gan) < 1e-12
    p0 = ccr4bp._ganymede_position(0.0, sys_c)
    p1 = ccr4bp._ganymede_position(sys_c.ganymede_synodic_period, sys_c)
    assert np.allclose(p0, p1, atol=1e-10)


def test_ganymede_commensurate_period() -> None:
    """n=1 commensurate period equals the synodic period; nonpositive n rejected."""
    sys_c = ccr4bp.jupiter_europa_ganymede_default()
    t1 = ccr4bp.ganymede_commensurate_period(sys_c.omega_gan, n=1)
    assert t1 == pytest.approx(sys_c.ganymede_synodic_period, rel=0.0, abs=1e-13)
    with pytest.raises(ValueError):
        ccr4bp.ganymede_commensurate_period(sys_c.omega_gan, n=0)
    with pytest.raises(ValueError):
        ccr4bp.ganymede_commensurate_period(sys_c.omega_gan, n=-1)
