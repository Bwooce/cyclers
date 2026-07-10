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
    """Transform Sun-Earth state to Earth-Moon rotating frame at time t."""
    omega_s = qbcp_sys.omega_sun_nondim
    mu_em = qbcp_sys.mu
    a_s = qbcp_sys.a_sun_nondim

    # Sun position in EM rotating frame
    theta_s = qbcp_sys.theta_sun0 + omega_s * t
    sx = a_s * math.cos(theta_s)
    sy = a_s * math.sin(theta_s)

    # Earth position is fixed at (-mu_em, 0, 0)
    ex = -mu_em
    ey = 0.0

    # Vector from Earth to Sun
    dx = sx - ex
    dy = sy - ey
    d_val = math.sqrt(dx * dx + dy * dy)

    # Rotation angle from SE to EM frame
    theta_rel = math.atan2(dy, dx)
    alpha = theta_rel - math.pi

    cos_a = math.cos(alpha)
    sin_a = math.sin(alpha)
    rot_mat = np.array([[cos_a, -sin_a, 0.0], [sin_a, cos_a, 0.0], [0.0, 0.0, 1.0]])

    dot_dx = -a_s * omega_s * math.sin(theta_s)
    dot_dy = a_s * omega_s * math.cos(theta_s)
    dot_d = (dx * dot_dx + dy * dot_dy) / d_val

    dot_alpha = (dx * dot_dy - dy * dot_dx) / (d_val * d_val)
    dot_rot = dot_alpha * np.array([[-sin_a, -cos_a, 0.0], [cos_a, -sin_a, 0.0], [0.0, 0.0, 0.0]])

    pos_se = state_se[:3]
    vel_se = state_se[3:]
    pos_rel_se = pos_se - np.array([1.0 - mu_se, 0.0, 0.0])

    pos_rel_em = d_val * (rot_mat @ pos_rel_se)
    pos_em = pos_rel_em + np.array([-mu_em, 0.0, 0.0])

    gamma = 1.0 + omega_s
    vel_rel_em = (
        dot_d * (rot_mat @ pos_rel_se)
        + d_val * (dot_rot @ pos_rel_se)
        + d_val * (rot_mat @ (gamma * vel_se))
    )

    return np.concatenate([pos_em, vel_rel_em])


def se_lyapunov_to_qbcp_torus_seed(
    orbit_se: SymmetricOrbit,
    qbcp_sys: qbcp.QBCPSystem,
    mu_se: float,
    n_samples: int = 5,
) -> tuple[NDArray[np.float64], int, float]:
    """Sample a Sun-Earth L2 Lyapunov orbit and transform it to EM frame."""
    t_em = orbit_se.period / (1.0 + qbcp_sys.omega_sun_nondim)

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
        state_se = sol.y[:, j]
        t_val = j * t_em / n_samples
        # Use mu = 0.0 for the starting system representation
        sys_mu0 = qbcp.QBCPSystem(
            mu=0.0,
            mu_sun=qbcp_sys.mu_sun,
            a_sun_nondim=qbcp_sys.a_sun_nondim,
            omega_sun_nondim=qbcp_sys.omega_sun_nondim,
            theta_sun0=qbcp_sys.theta_sun0,
        )
        u_samples[j] = se_to_em_transform(state_se, t_val, sys_mu0, mu_se)

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
