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
from scipy.optimize import brentq

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


# ---------------------------------------------------------------------------
# Published positive control: the L1 replacement orbit of the bicircular problem
# ---------------------------------------------------------------------------

# Jorba, Jorba-Cusco & Rosales (2020), "The vicinity of the Earth-Moon L1 point in the bicircular
# problem", Celestial Mechanics and Dynamical Astronomy 132:11. Table 1 (constants) and sections
# 3.1-3.2: the Sun is at (a_S cos theta, -a_S sin theta), theta = omega_S t (clockwise); the
# monodromy matrix of the periodic orbit that replaces L1 "has an hyperbolic eigenvalue close to
# 4.287 x 10^8"; the normalized logarithms of its elliptic eigenvalues are omega_1 =
# 2.32981963603288 and omega_2 = 2.26695149158478 (each defined up to +-(omega + k omega_S)); its
# (x, y) projection "revolves L1 twice in T_S units of time".
_JJR_2020_MU = 0.012150581623433623
_JJR_2020_M_SUN = 328900.54999999906
_JJR_2020_OMEGA_SUN = 0.925195985518289646
_JJR_2020_A_SUN = 388.81114302335106
_JJR_2020_OMEGA_1 = 2.32981963603288
_JJR_2020_OMEGA_2 = 2.26695149158478


def _l1_replacement_monodromy(n_seg: int = 24) -> tuple[list[FloatArray], FloatArray, float]:
    """The period-T_S orbit that replaces L1, by multiple shooting continued in the Sun's mass."""
    mu = _JJR_2020_MU
    period = 2.0 * math.pi / _JJR_2020_OMEGA_SUN
    dt = period / n_seg

    def d_omega(x: float) -> float:
        return (
            x
            - (1.0 - mu) * (x + mu) / abs(x + mu) ** 3
            - mu * (x - 1.0 + mu) / abs(x - 1.0 + mu) ** 3
        )

    x_l1 = float(brentq(d_omega, 0.5, 0.95, xtol=1e-15))
    nodes = [np.array([x_l1, 0.0, 0.0, 0.0, 0.0, 0.0]) for _ in range(n_seg)]
    stms: list[FloatArray] = []
    for eps in np.concatenate([[0.0], np.geomspace(1e-4, 1.0, 25)]):
        # The paper's frame has the Earth at (mu, 0); this module's is rotated by pi, which puts
        # the Sun at angle pi at t = 0 and leaves its sense unchanged.
        system = bcr4bp.BCR4BPSystem(
            mu=mu,
            mu_sun=float(eps) * _JJR_2020_M_SUN,
            a_sun_nondim=_JJR_2020_A_SUN,
            omega_sun_nondim=_JJR_2020_OMEGA_SUN,
            theta_sun0=math.pi,
        )
        for _ in range(40):
            residual = np.zeros(6 * n_seg)
            jac = np.zeros((6 * n_seg, 6 * n_seg))
            stms = []
            for i in range(n_seg):
                sol = solve_ivp(
                    bcr4bp.bcr4bp_stm_eom,
                    (i * dt, (i + 1) * dt),
                    np.concatenate([nodes[i], np.eye(6).ravel()]),
                    args=(system,),
                    method="DOP853",
                    rtol=1e-13,
                    atol=1e-13,
                )
                end = sol.y[:, -1]
                stm = end[6:].reshape(6, 6)
                stms.append(stm)
                j = (i + 1) % n_seg
                residual[6 * i : 6 * i + 6] = end[:6] - nodes[j]
                jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = stm
                jac[6 * i : 6 * i + 6, 6 * j : 6 * j + 6] -= np.eye(6)
            if float(np.linalg.norm(residual)) < 1e-11:
                break
            delta = np.linalg.solve(jac, -residual)
            nodes = [nodes[i] + delta[6 * i : 6 * i + 6] for i in range(n_seg)]
        else:
            raise AssertionError(f"L1 replacement did not converge at eps = {eps}")
    monodromy = np.eye(6)
    for stm in stms:
        monodromy = stm @ monodromy
    return nodes, monodromy, x_l1


def test_l1_replacement_orbit_matches_jorba_2020() -> None:
    """Published positive control for the bicircular model (#891).

    Measured with the corrected sense: unstable multiplier 4.287389e8, frequencies
    2.3298196303 and 2.266951491584771, two revolutions about L1 per period. With the Sun
    advancing counter-clockwise, as before 2026-10-04, the same computation gives 4.310e8,
    2.33046 and 2.26711, so this test discriminates the sense at the 1e-4 level.
    """
    nodes, monodromy, x_l1 = _l1_replacement_monodromy()
    period = 2.0 * math.pi / _JJR_2020_OMEGA_SUN
    eigenvalues = np.linalg.eigvals(monodromy)

    unstable = float(np.max(np.abs(eigenvalues)))
    assert abs(unstable / 1e8 - 4.287) < 5e-4  # printed as "close to 4.287 x 10^8"

    elliptic = [v for v in eigenvalues if abs(abs(v) - 1.0) < 1e-3 and v.imag > 0.0]
    assert len(elliptic) == 2
    printed = (_JJR_2020_OMEGA_1, _JJR_2020_OMEGA_2)
    found = []
    for value in elliptic:
        base = math.atan2(value.imag, value.real) / period
        candidates = [
            abs(sign * base + k * _JJR_2020_OMEGA_SUN) for sign in (1.0, -1.0) for k in range(-4, 5)
        ]
        found.append(min(candidates, key=lambda c: min(abs(c - p) for p in printed)))
    omega_1, omega_2 = sorted(found, reverse=True)
    # The in-plane frequency shares a block with the 4e8 multiplier and is resolved to about 1e-8;
    # the out-of-plane one is decoupled and agrees to the printed digits.
    assert math.isclose(omega_1, _JJR_2020_OMEGA_1, rel_tol=1e-7)
    assert math.isclose(omega_2, _JJR_2020_OMEGA_2, rel_tol=1e-12)

    angles = np.unwrap([math.atan2(s[1], s[0] - x_l1) for s in [*nodes, nodes[0]]])
    assert round(abs(angles[-1] - angles[0]) / (2.0 * math.pi)) == 2
