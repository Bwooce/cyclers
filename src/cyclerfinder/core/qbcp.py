"""Quasi-Bicircular Restricted 4-Body Problem (QBCP) Sun-Earth-Moon synodic dynamics.

This implements Andreu's (1998) Quasi-Bicircular Problem (QBCP) using the 8 Fourier
coefficient tables from Gimeno-Jorba (2018) Table 4.

To maintain coordinate consistency with the repository's BCR4BP and CR3BP conventions
(where Earth is at -mu and Moon is at 1-mu), the coordinates of the model are reflected
(x -> -x, y -> -y) relative to the paper's default convention. This maps the Sun's
rotation to the standard counter-clockwise direction and aligns the primary positions
perfectly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

# Fourier coefficients for alpha_1 to alpha_8 from Gimeno-Jorba (2018) Table 4.
# Each entry is a list starting at k=0.
# We apply the coordinate reflections to the coefficients directly:
# - alpha_1, alpha_2, alpha_3, alpha_6: unchanged (even/odd symmetry is preserved).
# - alpha_4, alpha_5, alpha_7, alpha_8: reflected (multiplied by -1.0).

_COEFFS_ALPHA1 = [
    1.001841608924835e00,
    5.767517726198399e-04,
    1.438777025507630e-02,
    -2.630362974972015e-06,
    1.176278356118933e-04,
    -38.068581391005552e-08,
    9.843249766501285e-07,
    -1.172054394418197e-09,
    8.311905970879588e-09,
    -1.408584238695393e-11,
    7.050713786466840e-11,
    -1.494259634910463e-13,
    5.982418979451232e-13,
]

_COEFFS_ALPHA2 = [
    0.0,
    -2.644376028499938e-04,
    -1.328686903400173e-02,
    9.386093208089751e-06,
    -1.218509057517414e-04,
    1.522127598557008e-07,
    -1.072102664277996e-06,
    1.889371261374048e-09,
    -9.324985038927486e-09,
    2.114490981280258e-11,
    -8.071111743144353e-11,
    2.218118050420168e-13,
    -7.036155161882012e-13,
]

_COEFFS_ALPHA3 = [
    9.999999999999983e-01,
    5.634125997553694e-04,
    1.889687440172882e-02,
    -9.911758802567132e-06,
    1.568708136031134e-04,
    -1.707762576173484e-07,
    1.319613679707437e-06,
    -2.136550041985646e-09,
    1.117168916673893e-08,
    -2.387253631031108e-11,
    9.490879622095902e-11,
    -2.462732581558427e-13,
    8.101067708009743e-13,
]

# Reflected (multiplied by -1.0)
_COEFFS_ALPHA4 = [
    9.755242327484885e-04,
    -2.154764362707107e00,
    -3.657484468968697e-04,
    -3.295673376166588e-03,
    -3.301031400812427e-07,
    -1.278840687376320e-05,
    2.623797952127926e-09,
    -6.533805514561511e-08,
    3.891720707783511e-11,
    -3.812275838944432e-10,
    3.907906049834876e-13,
    -2.407471187576443e-12,
]

# Reflected (multiplied by -1.0)
_COEFFS_ALPHA5 = [
    0.0,
    2.192570751040067e00,
    3.337210485472868e-04,
    3.295001430200974e-03,
    3.100635053052634e-07,
    1.277777336854128e-05,
    -2.652806405498111e-09,
    6.528479245085066e-08,
    -3.891720707783511e-11,
    3.812275838944432e-10,
    -3.907906049834876e-13,
    2.407471187576443e-12,
]

_COEFFS_ALPHA6 = [
    1.000907457708158e00,
    2.870921750053134e-04,
    7.187177998612875e-03,
    -2.351183147213254e-06,
    4.585758971122060e-05,
    -3.848683620107037e-08,
    3.270677504935666e-07,
    -4.406966481041876e-10,
    2.452600662570259e-09,
    -4.542938800673444e-12,
    1.892348855112616e-11,
    -4.178420101480123e-14,
    1.480048946961583e-13,
]

# Reflected (multiplied by -1.0)
_COEFFS_ALPHA7 = [
    6.314069568006227e-02,
    -3.885638623098048e02,
    -1.736910203345558e-01,
    -3.382908071669699e00,
    -1.574837565380491e-04,
    -0.02936360489004438,
    1.224434550116014e-05,
    -2.538935434262443e-04,
    2.278929040007574e-07,
    -2.190432706181655e-06,
    3.033311961234353e-09,
    -1.886971545290216e-08,
    -3.432375106898453e-11,
    1.611513703999101e-10,
]

# Reflected (multiplied by -1.0)
_COEFFS_ALPHA8 = [
    0.0,
    389.7437256237654,
    0.1734279166322518,
    3.38569648664212,
    0.0001555886632413398,
    0.02937582671967532,
    -1.225851213107933e-05,
    2.539596887692642e-04,
    -2.280029220202363e-07,
    2.19083462442904e-06,
    -3.036109035120856e-09,
    -1.887457647579322e-08,
    3.432375106898453e-11,
    -1.631723641506449e-10,
]

# Constants from Gimeno-Jorba (2018) Table 3
_QBCP_MU_EM: float = 0.012150581623433623
_QBCP_MU_S: float = 328900.54999999906
_QBCP_A_S: float = 388.81114302335106
_QBCP_OMEGA_S: float = 0.92519598551829646


@dataclass(frozen=True)
class QBCPSystem:
    """Quasi-Bicircular Restricted 4-Body Problem parameters.

    Attributes
    ----------
    mu :
        Earth-Moon mass ratio.
    mu_sun :
        Sun mass in units where Earth + Moon = 1.
    a_sun_nondim :
        Sun distance parameter in EM-distance units.
    omega_sun_nondim :
        Sun angular frequency in the synodic frame.
    theta_sun0 :
        Sun initial phase angle at t=0.
    """

    mu: float
    mu_sun: float
    a_sun_nondim: float
    omega_sun_nondim: float
    theta_sun0: float = 0.0

    @property
    def sun_period_tu(self) -> float:
        """Synodic-frame Sun period."""
        return 2.0 * math.pi / self.omega_sun_nondim


def qbcp_default() -> QBCPSystem:
    """QBCP parameter set matching Gimeno-Jorba (2018) Table 3."""
    return QBCPSystem(
        mu=_QBCP_MU_EM,
        mu_sun=_QBCP_MU_S,
        a_sun_nondim=_QBCP_A_S,
        omega_sun_nondim=_QBCP_OMEGA_S,
        theta_sun0=0.0,
    )


def _evaluate_alpha(theta: float, coeffs: list[float], is_even: bool) -> float:
    """Evaluate Fourier series for a given alpha coefficient at phase theta."""
    val = coeffs[0]
    if is_even:
        for k in range(1, len(coeffs)):
            val += coeffs[k] * math.cos(k * theta)
    else:
        for k in range(1, len(coeffs)):
            val += coeffs[k] * math.sin(k * theta)
    return val


def evaluate_alphas(t: float, system: QBCPSystem) -> NDArray[np.float64]:
    """Evaluate all 8 alpha(t) functions at time t.

    Returns array of size 9 (1-indexed, element 0 is unused).
    """
    theta = system.theta_sun0 + system.omega_sun_nondim * t
    alphas = np.zeros(9, dtype=np.float64)
    alphas[1] = _evaluate_alpha(theta, _COEFFS_ALPHA1, is_even=True)
    alphas[2] = _evaluate_alpha(theta, _COEFFS_ALPHA2, is_even=True)
    alphas[3] = _evaluate_alpha(theta, _COEFFS_ALPHA3, is_even=False)
    alphas[4] = _evaluate_alpha(theta, _COEFFS_ALPHA4, is_even=True)
    alphas[5] = _evaluate_alpha(theta, _COEFFS_ALPHA5, is_even=False)
    alphas[6] = _evaluate_alpha(theta, _COEFFS_ALPHA6, is_even=True)
    alphas[7] = _evaluate_alpha(theta, _COEFFS_ALPHA7, is_even=True)
    alphas[8] = _evaluate_alpha(theta, _COEFFS_ALPHA8, is_even=False)
    return alphas


def state_pv_to_pm(
    state_pv: NDArray[np.float64], t: float, system: QBCPSystem
) -> NDArray[np.float64]:
    """Convert state from position-velocity (PV) to position-momentum (PM) representation."""
    x, y, z, vx, vy, vz = state_pv
    alphas = evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]
    px = (vx - a2 * x - a3 * y) / a1
    py = (vy - a2 * y + a3 * x) / a1
    pz = (vz - a2 * z) / a1
    return np.array([x, y, z, px, py, pz], dtype=np.float64)


def state_pm_to_pv(
    state_pm: NDArray[np.float64], t: float, system: QBCPSystem
) -> NDArray[np.float64]:
    """Convert state from position-momentum (PM) to position-velocity (PV) representation."""
    x, y, z, px, py, pz = state_pm
    alphas = evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]
    vx = a1 * px + a2 * x + a3 * y
    vy = a1 * py + a2 * y - a3 * x
    vz = a1 * pz + a2 * z
    return np.array([x, y, z, vx, vy, vz], dtype=np.float64)


def transformation_jacobian(t: float, system: QBCPSystem) -> NDArray[np.float64]:
    """Get the Jacobian matrix M(t) mapping d(S_PV)/d(S_PM)."""
    alphas = evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]
    jac_m = np.zeros((6, 6), dtype=np.float64)
    # Positions are unchanged
    jac_m[0, 0] = 1.0
    jac_m[1, 1] = 1.0
    jac_m[2, 2] = 1.0
    # Velocities derivatives w.r.t (x, y, z, px, py, pz)
    jac_m[3, 0] = a2
    jac_m[3, 1] = a3
    jac_m[3, 3] = a1

    jac_m[4, 0] = -a3
    jac_m[4, 1] = a2
    jac_m[4, 4] = a1

    jac_m[5, 2] = a2
    jac_m[5, 5] = a1
    return jac_m


def transformation_jacobian_inverse(t: float, system: QBCPSystem) -> NDArray[np.float64]:
    """Get the inverse Jacobian matrix M_inv(t) mapping d(S_PM)/d(S_PV)."""
    alphas = evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]
    jac_minv = np.zeros((6, 6), dtype=np.float64)
    jac_minv[0, 0] = 1.0
    jac_minv[1, 1] = 1.0
    jac_minv[2, 2] = 1.0

    jac_minv[3, 0] = -a2 / a1
    jac_minv[3, 1] = -a3 / a1
    jac_minv[3, 3] = 1.0 / a1

    jac_minv[4, 0] = a3 / a1
    jac_minv[4, 1] = -a2 / a1
    jac_minv[4, 4] = 1.0 / a1

    jac_minv[5, 2] = -a2 / a1
    jac_minv[5, 5] = 1.0 / a1
    return jac_minv


def qbcp_eom(t: float, state_pm: NDArray[np.float64], system: QBCPSystem) -> NDArray[np.float64]:
    """QBCP equations of motion in canonical variables (PM representation)."""
    x, y, z, px, py, pz = state_pm
    alphas = evaluate_alphas(t, system)
    a1, a2, a3, a4, a5, a6, xs, ys = alphas[1:9]
    mu = system.mu

    # Position derivatives
    dx = a1 * px + a2 * x + a3 * y
    dy = a1 * py + a2 * y - a3 * x
    dz = a1 * pz + a2 * z

    # Distances
    rpe2 = (x + mu) ** 2 + y * y + z * z
    rpm2 = (x - 1.0 + mu) ** 2 + y * y + z * z
    rps2 = (x - xs) ** 2 + (y - ys) ** 2 + z * z

    rpe3 = rpe2 * math.sqrt(rpe2)
    rpm3 = rpm2 * math.sqrt(rpm2)
    rps3 = rps2 * math.sqrt(rps2)

    # Potential term derivatives w.r.t (x, y, z)
    pot_x = (
        (1.0 - mu) * (x + mu) / rpe3 + mu * (x - 1.0 + mu) / rpm3 + system.mu_sun * (x - xs) / rps3
    )
    pot_y = (1.0 - mu) * y / rpe3 + mu * y / rpm3 + system.mu_sun * (y - ys) / rps3
    pot_z = (1.0 - mu) * z / rpe3 + mu * z / rpm3 + system.mu_sun * z / rps3

    # Momenta derivatives
    dpx = -a2 * px + a3 * py - a4 - a6 * pot_x
    dpy = -a2 * py - a3 * px - a5 - a6 * pot_y
    dpz = -a2 * pz - a6 * pot_z

    return np.array([dx, dy, dz, dpx, dpy, dpz], dtype=np.float64)


def qbcp_potential_second_derivatives(
    x: float, y: float, z: float, t: float, system: QBCPSystem
) -> tuple[float, float, float, float, float, float]:
    """Get the potential second derivatives (uxx, uyy, uzz, uxy, uxz, uyz)."""
    alphas = evaluate_alphas(t, system)
    a6, xs, ys = alphas[6], alphas[7], alphas[8]
    mu = system.mu

    rpe2 = (x + mu) ** 2 + y * y + z * z
    rpm2 = (x - 1.0 + mu) ** 2 + y * y + z * z
    rps2 = (x - xs) ** 2 + (y - ys) ** 2 + z * z

    rpe3 = rpe2 * math.sqrt(rpe2)
    rpm3 = rpm2 * math.sqrt(rpm2)
    rps3 = rps2 * math.sqrt(rps2)

    rpe5 = rpe3 * rpe2
    rpm5 = rpm3 * rpm2
    rps5 = rps3 * rps2

    om1 = 1.0 - mu
    # Newtonian potential term second derivatives
    uxx = (
        -om1 * (1.0 / rpe3 - 3.0 * (x + mu) ** 2 / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * (x - 1.0 + mu) ** 2 / rpm5)
        - system.mu_sun * (1.0 / rps3 - 3.0 * (x - xs) ** 2 / rps5)
    )

    uyy = (
        -om1 * (1.0 / rpe3 - 3.0 * y * y / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * y * y / rpm5)
        - system.mu_sun * (1.0 / rps3 - 3.0 * (y - ys) ** 2 / rps5)
    )

    uzz = (
        -om1 * (1.0 / rpe3 - 3.0 * z * z / rpe5)
        - mu * (1.0 / rpm3 - 3.0 * z * z / rpm5)
        - system.mu_sun * (1.0 / rps3 - 3.0 * z * z / rps5)
    )

    uxy = (
        3.0 * om1 * (x + mu) * y / rpe5
        + 3.0 * mu * (x - 1.0 + mu) * y / rpm5
        + 3.0 * system.mu_sun * (x - xs) * (y - ys) / rps5
    )

    uxz = (
        3.0 * om1 * (x + mu) * z / rpe5
        + 3.0 * mu * (x - 1.0 + mu) * z / rpm5
        + 3.0 * system.mu_sun * (x - xs) * z / rps5
    )

    uyz = (
        3.0 * om1 * y * z / rpe5
        + 3.0 * mu * y * z / rpm5
        + 3.0 * system.mu_sun * (y - ys) * z / rps5
    )

    return (
        a6 * uxx,
        a6 * uyy,
        a6 * uzz,
        a6 * uxy,
        a6 * uxz,
        a6 * uyz,
    )


def qbcp_stm_eom(
    t: float, state_and_stm: NDArray[np.float64], system: QBCPSystem
) -> NDArray[np.float64]:
    """QBCP variational equations: state (6) + flattened 6x6 STM (36) in PM."""
    state_pm = state_and_stm[:6]
    stm_flat = state_and_stm[6:]
    x, y, z = state_pm[:3]

    # Derivative of state
    dstate = qbcp_eom(t, state_pm, system)

    # Coefficients for Jacobian matrix
    alphas = evaluate_alphas(t, system)
    a1, a2, a3 = alphas[1], alphas[2], alphas[3]

    # Potential derivatives
    uxx, uyy, uzz, uxy, uxz, uyz = qbcp_potential_second_derivatives(x, y, z, t, system)

    # Jacobian matrix jac_a
    jac_a = np.zeros((6, 6), dtype=np.float64)
    jac_a[0, 0] = a2
    jac_a[0, 1] = a3
    jac_a[0, 3] = a1

    jac_a[1, 0] = -a3
    jac_a[1, 1] = a2
    jac_a[1, 4] = a1

    jac_a[2, 2] = a2
    jac_a[2, 5] = a1

    jac_a[3, 0] = uxx
    jac_a[3, 1] = uxy
    jac_a[3, 2] = uxz
    jac_a[3, 3] = -a2
    jac_a[3, 4] = a3

    jac_a[4, 0] = uxy
    jac_a[4, 1] = uyy
    jac_a[4, 2] = uyz
    jac_a[4, 3] = -a3
    jac_a[4, 4] = -a2

    jac_a[5, 0] = uxz
    jac_a[5, 1] = uyz
    jac_a[5, 2] = uzz
    jac_a[5, 5] = -a2

    # STM derivative
    phi = stm_flat.reshape((6, 6))
    dphi = jac_a @ phi

    return np.concatenate([dstate, dphi.flatten()])


def propagate_qbcp_pv(
    state_pv0: NDArray[np.float64],
    t_span: tuple[float, float],
    system: QBCPSystem,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    with_stm: bool = False,
    collision_check: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Propagate QBCP in PV representation.

    If with_stm=True, also propagates and maps the 6x6 STM.

    ``collision_check`` (default ``False``) terminates integration early on a
    close approach to either primary (r < ~7,700 km of the first, r < ~3,800
    km of the second) -- an opt-in guard for exploratory manifold sweeps.
    Default OFF: the #538 cislunar cycler correction this model exists for
    deliberately targets close Earth/Moon approaches (that is the point of a
    cislunar transport orbit), and an unconditional guard at this radius
    silently truncated an analogous BCR4BP capture-state propagation
    (`core/bcr4bp.py`, found 2026-07-09) -- callers that want the guard must
    opt in explicitly.

    Returns:
        times: Array of times.
        states_pv: Array of states of shape (N, 6) or (N, 42) if with_stm=True.
    """
    t0, _tf = t_span
    state_pm0 = state_pv_to_pm(state_pv0, t0, system)

    def _collision_event(t: float, y: NDArray[np.float64]) -> float:
        mu = system.mu
        x, yc, z = y[0], y[1], y[2]
        r1_sq = (x + mu) ** 2 + yc**2 + z**2
        r2_sq = (x - 1.0 + mu) ** 2 + yc**2 + z**2
        if r1_sq < 0.0004 or r2_sq < 0.0001:
            return 0.0
        return 1.0

    _collision_event.terminal = True  # type: ignore[attr-defined]
    events = _collision_event if collision_check else None

    if with_stm:
        # Initialize identity STM
        stm0 = np.eye(6, dtype=np.float64).flatten()
        y0 = np.concatenate([state_pm0, stm0])

        def fun(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
            return qbcp_stm_eom(t, y, system)
    else:
        y0 = state_pm0

        def fun(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
            return qbcp_eom(t, y, system)

    sol = solve_ivp(fun, t_span, y0, method="DOP853", rtol=rtol, atol=atol, events=events)

    if not sol.success:
        raise RuntimeError(f"QBCP propagation failed: {sol.message}")

    times = sol.t
    states_pm = sol.y.T
    n_steps = len(times)

    if with_stm:
        states_pv = np.zeros((n_steps, 42), dtype=np.float64)
        for i in range(n_steps):
            ti = times[i]
            spm = states_pm[i, :6]
            stm_pm = states_pm[i, 6:].reshape((6, 6))

            # Map state
            spv = state_pm_to_pv(spm, ti, system)
            states_pv[i, :6] = spv

            # Map STM: Phi_PV = M(ti) * Phi_PM * M_inv(t0)
            m_ti = transformation_jacobian(ti, system)
            minv_t0 = transformation_jacobian_inverse(t0, system)
            stm_pv = m_ti @ stm_pm @ minv_t0
            states_pv[i, 6:] = stm_pv.flatten()
    else:
        states_pv = np.zeros((n_steps, 6), dtype=np.float64)
        for i in range(n_steps):
            states_pv[i] = state_pm_to_pv(states_pm[i], times[i], system)

    return times, states_pv
