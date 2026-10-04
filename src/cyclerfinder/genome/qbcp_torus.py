"""QBCP invariant 2-tori corrector and coordinate transform bridge.

Parameters and algorithms are designed for the non-autonomous time-periodic QBCP
model. Per the design guidelines, all exclamation marks are avoided in comments
and docstrings.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.genome.qp_tori import evaluate_invariant_circle
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit


@dataclass(frozen=True)
class QBCPTorus:
    """Quasi-periodic invariant 2-torus in the QBCP rotating frame."""

    system: qbcp.QBCPSystem
    omega_long: float
    omega_trans: float
    rho: float
    t_strob: float
    fourier_coeffs: NDArray[np.complex128]
    n_modes: int
    n_samples: int
    invariance_residual: float
    converged: bool
    n_iter: int
    notes: str = ""
    extras: dict[str, float] = field(default_factory=dict)


def se_to_em_transform(
    state_se: NDArray[np.float64], t: float, qbcp_sys: qbcp.QBCPSystem, mu_se: float
) -> NDArray[np.float64]:
    """Map a Sun-Earth CR3BP state to the QBCP Earth-Moon frame at time t (a SEED, not exact).

    The Sun-Earth secondary (mass ``mu_se = 1 / (mu_sun + 1)``, the Earth+Moon mass point)
    maps to the origin of the QBCP frame (the Earth-Moon barycentre) and the Sun-Earth
    primary to the QBCP's own Sun position ``(alpha_7, alpha_8)(t)``, which starts near
    angle pi in this module's frame and regresses. With ``D`` and ``theta_S`` the distance
    and angle of that Sun position, ``R`` the rotation by ``theta_S - pi`` and ``p`` the
    secondary-centred Sun-Earth position:

        r = D R p,
        v = dD/dt R p + D dtheta_S/dt J R p + D R (n_S v_SE),    J = [[0, -1], [1, 0]],

    the derivatives of ``(alpha_7, alpha_8)`` being taken by central difference. The Sun
    and the secondary map exactly onto the QBCP Sun (position and velocity) and the
    origin. The time scale is not exact: Sun-Earth time is taken to advance at the mean
    inertial rate ``n_S = 1 - omega_S``, while in the coherent model the Sun's motion is
    not uniform, and the Moon has no Sun-Earth counterpart. There is no exact identity
    for the flow; the measured bound (2026-10-04, ``tests/genome/test_qbcp_torus.py``) is
    that with the Moon's mass set to zero the QBCP flow and the transformed Sun-Earth
    flow agree to 7.5e-7 in position after 6 time units, a few Earth-Moon distances
    from the Earth (7.2e-4 at the physical Moon mass). Treat it as a seed generator.
    (#891/#892: until 2026-10-04 this placed the Sun on a prograde circle starting at
    angle 0, used the factor ``1 + omega_S`` and put the secondary at the Earth.)
    """
    alphas = qbcp.evaluate_alphas(t, qbcp_sys)
    xs, ys = float(alphas[7]), float(alphas[8])
    h = 1e-5
    alphas_p = qbcp.evaluate_alphas(t + h, qbcp_sys)
    alphas_m = qbcp.evaluate_alphas(t - h, qbcp_sys)
    dxs = float(alphas_p[7] - alphas_m[7]) / (2.0 * h)
    dys = float(alphas_p[8] - alphas_m[8]) / (2.0 * h)

    d_val = math.hypot(xs, ys)
    dot_d = (xs * dxs + ys * dys) / d_val
    dot_theta = (xs * dys - ys * dxs) / (d_val * d_val)

    alpha = math.atan2(ys, xs) - math.pi
    cos_a = math.cos(alpha)
    sin_a = math.sin(alpha)
    rot_mat = np.array([[cos_a, -sin_a, 0.0], [sin_a, cos_a, 0.0], [0.0, 0.0, 1.0]])
    j_mat = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    pos_rel = rot_mat @ (state_se[:3] - np.array([1.0 - mu_se, 0.0, 0.0]))
    n_s = 1.0 - qbcp_sys.omega_sun_nondim
    pos_em = d_val * pos_rel
    vel_em = dot_d * pos_rel + d_val * dot_theta * (j_mat @ pos_rel)
    vel_em += d_val * n_s * (rot_mat @ state_se[3:])

    return np.concatenate([pos_em, vel_em])


def se_lyapunov_to_qbcp_torus_seed(
    orbit_se: SymmetricOrbit,
    qbcp_sys: qbcp.QBCPSystem,
    mu_se: float,
    n_samples: int = 5,
) -> tuple[NDArray[np.float64], int, float]:
    """Seed a QBCP invariant circle from a Sun-Earth L2 Lyapunov orbit (unvalidated seed).

    As :func:`cyclerfinder.genome.bcr4bp_torus.se_lyapunov_to_bcr4bp_torus_seed`: every
    sample is transformed at the section time t = 0 and ``rho = 2 pi T_s / T_em`` with
    ``T_em = period / (1 - omega_S)``.
    """
    t_em = orbit_se.period / (1.0 - qbcp_sys.omega_sun_nondim)

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit_se.period),
        np.array([orbit_se.x0, 0.0, 0.0, 0.0, orbit_se.ydot0, 0.0]),
        args=(mu_se,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, orbit_se.period, n_samples, endpoint=False),
    )

    u_samples = np.zeros((n_samples, 6))
    for j in range(n_samples):
        u_samples[j] = se_to_em_transform(sol.y[:, j], 0.0, qbcp_sys, mu_se)

    coeffs = np.fft.fft(u_samples, axis=0) / n_samples
    t_s = 2.0 * math.pi / qbcp_sys.omega_sun_nondim
    rho = (2.0 * math.pi * t_s / t_em) % (2.0 * math.pi)
    if rho > math.pi:
        rho -= 2.0 * math.pi

    n_modes = (n_samples - 1) // 2
    n_unk = 6 + 12 * n_modes
    x0 = np.zeros(n_unk + 1)
    x0[0:6] = np.real(coeffs[0, :])
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        x0[i0 : i0 + 6] = np.real(coeffs[n, :])
        x0[i0 + 6 : i0 + 12] = np.imag(coeffs[n, :])
    x0[-1] = rho

    # Phase-pin gauge: pin ``Im(c_1[phase_pin_idx]) = 0`` to kill the rotation
    # invariance of the invariant-circle parameterization. The gauge derivative
    # with respect to a circle rotation phi is ``Re(c_1[phase_pin_idx])``, so the
    # pin coordinate MUST be one with a large real part or the gauge is singular
    # (the corrector then cannot fix the rotational phase and stalls). Pick the
    # coordinate with the largest real part, matching qp_tori.correct_qp_torus.
    # Selecting argmax|Im| instead would pin a near-purely-imaginary coordinate
    # (Re ~ 0), a singular gauge -- the #544 Earth-Moon L2 non-convergence bug.
    phase_pin_idx = int(np.argmax(np.abs(np.real(coeffs[1, :]))))
    amplitude_pin = float(np.linalg.norm(coeffs[1, :]))

    return x0, phase_pin_idx, amplitude_pin


def qbcp_torus_residual(
    x_unk: NDArray[np.float64],
    system: qbcp.QBCPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    amplitude_pin: float,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> NDArray[np.float64]:
    """Calculate the GMOS residual for QBCP torus."""
    coeffs_u, rho = x_unk[:-1], x_unk[-1]
    t_s = 2.0 * math.pi / system.omega_sun_nondim

    n_total = 2 * n_modes + 1
    coeffs = np.zeros((n_total, 6), dtype=np.complex128)
    coeffs[0, :] = coeffs_u[0:6].astype(np.complex128)
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        coeffs[n, :] = coeffs_u[i0 : i0 + 6] + 1j * coeffs_u[i0 + 6 : i0 + 12]
        coeffs[n_total - n, :] = np.conj(coeffs[n, :])

    thetas = 2 * math.pi * np.arange(n_samples) / n_samples
    u_s = evaluate_invariant_circle(coeffs, thetas)

    phi_samples = np.zeros_like(u_s)
    for j in range(n_samples):
        try:
            _, states_pv = qbcp.propagate_qbcp_pv(
                u_s[j], (0.0, t_s), system, with_stm=False, rtol=rtol, atol=atol
            )
            phi_samples[j] = states_pv[-1]
        except RuntimeError:
            return np.full(8 + 24 * n_modes, 1e5)

    phi_fft = np.fft.fft(phi_samples, axis=0) / n_samples
    expected = np.zeros((n_samples, 6), dtype=np.complex128)
    expected[0, :] = coeffs[0, :]
    for n in range(1, n_modes + 1):
        expected[n, :] = coeffs[n, :] * np.exp(1j * n * rho)
    for n in range(1, n_modes + 1):
        expected[n_samples - n, :] = coeffs[n_total - n, :] * np.exp(-1j * n * rho)

    f_res = phi_fft - expected
    if n_samples > n_total:
        for n in range(n_modes + 1, n_samples - n_modes):
            f_res[n, :] = 0.0

    parts = []
    parts.append(np.real(f_res[0, :]))
    for n in range(1, n_modes + 1):
        parts.append(np.real(f_res[n, :]))
        parts.append(np.imag(f_res[n, :]))
        parts.append(np.real(f_res[n_total - n, :]))
        parts.append(np.imag(f_res[n_total - n, :]))
    parts.append(np.array([float(np.imag(coeffs[1, phase_pin_idx]))]))
    amp = float(np.linalg.norm(coeffs[1, :]))
    parts.append(np.array([amp - amplitude_pin]))
    return np.concatenate(parts)


def correct_qbcp_torus(
    system: qbcp.QBCPSystem,
    x0: NDArray[np.float64],
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    amplitude_pin: float,
    *,
    tol: float = 1e-8,
    max_iter: int = 30,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> QBCPTorus:
    """Newton-correct a QBCP torus at a fixed stroboscopic period T_s."""
    res_gmos = least_squares(
        qbcp_torus_residual,
        x0,
        args=(system, n_modes, n_samples, phase_pin_idx, amplitude_pin),
        kwargs={"rtol": rtol, "atol": atol},
        method="lm",
        xtol=tol * 1e-2,
        ftol=tol * 1e-2,
        max_nfev=max_iter * (len(x0) + 1),
    )

    converged = bool(np.linalg.norm(res_gmos.fun) < tol)
    x_final = res_gmos.x
    coeffs_u, rho_final = x_final[:-1], x_final[-1]

    n_total = 2 * n_modes + 1
    coeffs = np.zeros((n_total, 6), dtype=np.complex128)
    coeffs[0, :] = coeffs_u[0:6].astype(np.complex128)
    for n in range(1, n_modes + 1):
        i0 = 6 + (n - 1) * 12
        coeffs[n, :] = coeffs_u[i0 : i0 + 6] + 1j * coeffs_u[i0 + 6 : i0 + 12]
        coeffs[n_total - n, :] = np.conj(coeffs[n, :])

    t_s = 2.0 * math.pi / system.omega_sun_nondim
    omega_long = 2.0 * math.pi / t_s
    omega_trans = rho_final * omega_long / (2.0 * math.pi)

    return QBCPTorus(
        system=system,
        omega_long=omega_long,
        omega_trans=omega_trans,
        rho=rho_final,
        t_strob=t_s,
        fourier_coeffs=coeffs,
        n_modes=n_modes,
        n_samples=n_samples,
        invariance_residual=float(np.linalg.norm(res_gmos.fun)),
        converged=converged,
        n_iter=int(res_gmos.nfev),
    )


def evaluate_qbcp_torus(
    torus: QBCPTorus,
    theta_long: float,
    theta_trans: float,
) -> NDArray[np.float64]:
    """Return the state on the QBCP torus at (theta_long, theta_trans)."""
    theta_long_red = float(theta_long) % (2.0 * math.pi)
    theta_trans_red = float(theta_trans) % (2.0 * math.pi)
    dt = theta_long_red / torus.omega_long
    theta_seed = (theta_trans_red - torus.rho * theta_long_red / (2.0 * math.pi)) % (2.0 * math.pi)
    u0 = evaluate_invariant_circle(torus.fourier_coeffs, theta_seed)
    if dt == 0.0:
        return u0
    _, states_pv = qbcp.propagate_qbcp_pv(u0, (0.0, dt), torus.system, rtol=1e-8, atol=1e-8)
    return np.asarray(states_pv[-1], dtype=np.float64)
