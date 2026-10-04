"""#891: the Sun's sense of motion in the bicircular model.

In the Earth-Moon rotating frame the Sun moves CLOCKWISE: its inertial mean motion
``n_S = 1 - omega_S`` is prograde and slower than the frame's own rate of 1, so its angle in
the frame is ``theta_sun0 - omega_S * t``. Andreu (1998), "The Quasi-bicircular Problem",
section 1.3: the angle is "clockwise measured ... in synodical coordinates, the Sun rotates in
reverse sense".

From 2026-06-16 to 2026-10-04 ``core/bcr4bp.py`` advanced the Sun counter-clockwise, which is a
Sun with an inertial period of 14.19 days. Nothing caught it because every test compared the
model with itself. The tests here compare it with an independent derivation: the same physics
written in a NON-rotating frame, where the sense is not a matter of convention.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp

FloatArray = NDArray[np.float64]

# Sidereal month and year, days (IAU nominal values, to the precision needed here).
_SIDEREAL_MONTH_D = 27.321661
_SIDEREAL_YEAR_D = 365.256363


def _rot_to_inertial(state: FloatArray, t: float) -> FloatArray:
    c, s = math.cos(t), math.sin(t)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    r = rot @ state[:3]
    v = rot @ state[3:] + np.cross([0.0, 0.0, 1.0], r)
    return np.concatenate([r, v])


def _inertial_to_rot(state: FloatArray, t: float) -> FloatArray:
    c, s = math.cos(t), math.sin(t)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    r = rot.T @ state[:3]
    v = rot.T @ (state[3:] - np.cross([0.0, 0.0, 1.0], state[:3]))
    return np.concatenate([r, v])


def _inertial_bicircular_rhs(
    t: float, state: FloatArray, system: bcr4bp.BCR4BPSystem
) -> FloatArray:
    """Bicircular model in a non-rotating frame centred on the Earth-Moon barycentre.

    Earth and Moon move counter-clockwise on circles at rate 1. The Sun moves counter-clockwise
    at its inertial rate ``1 - omega_S`` (prograde, one revolution per year). The Sun acts on
    the particle directly and through the barycentre's own acceleration (indirect term).
    """
    mu, mu_s, a_s = system.mu, system.mu_sun, system.a_sun_nondim
    n_sun = 1.0 - system.omega_sun_nondim
    r = state[:3]
    c, s = math.cos(t), math.sin(t)
    r_e = np.array([-mu * c, -mu * s, 0.0])
    r_m = np.array([(1.0 - mu) * c, (1.0 - mu) * s, 0.0])
    th = system.theta_sun0 + n_sun * t
    r_s = a_s * np.array([math.cos(th), math.sin(th), 0.0])
    acc = (
        -(1.0 - mu) * (r - r_e) / np.linalg.norm(r - r_e) ** 3
        - mu * (r - r_m) / np.linalg.norm(r - r_m) ** 3
        - mu_s * (r - r_s) / np.linalg.norm(r - r_s) ** 3
        - mu_s * r_s / a_s**3
    )
    return np.concatenate([state[3:], acc])


def _propagate_inertial(x0: FloatArray, t_end: float, system: bcr4bp.BCR4BPSystem) -> FloatArray:
    sol = solve_ivp(
        _inertial_bicircular_rhs,
        (0.0, t_end),
        _rot_to_inertial(x0, 0.0),
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    return _inertial_to_rot(sol.y[:, -1], t_end)


def _propagate_rotating(x0: FloatArray, t_end: float, system: bcr4bp.BCR4BPSystem) -> FloatArray:
    sol = solve_ivp(
        bcr4bp.bcr4bp_eom,
        (0.0, t_end),
        x0,
        args=(system,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    out: FloatArray = sol.y[:, -1]
    return out


_X0 = np.array([0.5, 0.1, 0.02, 0.05, 0.6, 0.01])
_T_END = 6.0


def test_inertial_frame_control_with_the_sun_off() -> None:
    """The two frames agree when there is no Sun, so the comparison itself is sound."""
    system = dataclasses.replace(bcr4bp.andreu_default(), mu_sun=0.0)
    diff = np.linalg.norm(
        _propagate_rotating(_X0, _T_END, system) - _propagate_inertial(_X0, _T_END, system)
    )
    assert diff < 1e-10


def test_rotating_model_matches_the_inertial_bicircular_model() -> None:
    """With the Sun on, the rotating-frame model reproduces the non-rotating derivation.

    With the Sun advancing the wrong way the difference after 6 time units is 0.116, larger
    than the 0.017 obtained by leaving the Sun out altogether.
    """
    system = bcr4bp.andreu_default()
    truth = _propagate_inertial(_X0, _T_END, system)
    diff = np.linalg.norm(_propagate_rotating(_X0, _T_END, system) - truth)
    assert diff < 1e-10
    no_sun = dataclasses.replace(system, mu_sun=0.0)
    assert np.linalg.norm(_propagate_rotating(_X0, _T_END, no_sun) - truth) > 1e-3


def test_rotating_model_matches_inertial_at_a_nonzero_sun_phase() -> None:
    system = dataclasses.replace(bcr4bp.andreu_default(), theta_sun0=1.1)
    diff = np.linalg.norm(
        _propagate_rotating(_X0, _T_END, system) - _propagate_inertial(_X0, _T_END, system)
    )
    assert diff < 1e-10


def test_sun_regresses_in_the_rotating_frame_and_takes_a_year_inertially() -> None:
    system = bcr4bp.andreu_default()
    dt = 1e-3
    x0, y0, _ = bcr4bp._sun_position(0.0, system)
    x1, y1, _ = bcr4bp._sun_position(dt, system)
    rate_rot = (math.atan2(y1, x1) - math.atan2(y0, x0)) / dt
    assert rate_rot < 0.0  # clockwise in the rotating frame
    rate_inertial = 1.0 + rate_rot  # the frame itself turns at rate 1
    period_days = _SIDEREAL_MONTH_D / rate_inertial
    assert math.isclose(period_days, _SIDEREAL_YEAR_D, rel_tol=1e-3)


def test_theta_sun0_is_the_suns_angle_at_time_zero() -> None:
    system = dataclasses.replace(bcr4bp.andreu_default(), theta_sun0=0.7)
    x, y, z = bcr4bp._sun_position(0.0, system)
    assert z == 0.0
    assert math.isclose(math.atan2(y, x), 0.7, abs_tol=1e-14)
    assert math.isclose(math.hypot(x, y), system.a_sun_nondim, rel_tol=1e-14)
