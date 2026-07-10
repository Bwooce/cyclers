"""BCR4BP invariant 2-tori corrector and coordinate transform bridge.

Parameters and algorithms are designed for the non-autonomous time-periodic BCR4BP
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

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import evaluate_invariant_circle
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit


@dataclass(frozen=True)
class BCR4BPTorus:
    """Quasi-periodic invariant 2-torus in the BCR4BP rotating frame."""

    system: bcr4bp.BCR4BPSystem
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
    state_se: NDArray[np.float64], t: float, bcr_sys: bcr4bp.BCR4BPSystem, mu_SE: float
) -> NDArray[np.float64]:
    """Transform Sun-Earth state to Earth-Moon rotating frame at time t."""
    omega_S = bcr_sys.omega_sun_nondim
    mu_EM = bcr_sys.mu
    a_S = bcr_sys.a_sun_nondim

    # Sun position in EM rotating frame
    theta_S = bcr_sys.theta_sun0 + omega_S * t
    sx = a_S * math.cos(theta_S)
    sy = a_S * math.sin(theta_S)

    # Earth position is fixed at (-mu_EM, 0, 0)
    ex = -mu_EM
    ey = 0.0

    # Vector from Earth to Sun
    dx = sx - ex
    dy = sy - ey
    D = math.sqrt(dx * dx + dy * dy)

    # Rotation angle from SE to EM frame
    theta_rel = math.atan2(dy, dx)
    alpha = theta_rel - math.pi

    cos_a = math.cos(alpha)
    sin_a = math.sin(alpha)
    R = np.array([[cos_a, -sin_a, 0.0], [sin_a, cos_a, 0.0], [0.0, 0.0, 1.0]])

    dot_dx = -a_S * omega_S * math.sin(theta_S)
    dot_dy = a_S * omega_S * math.cos(theta_S)
    dot_D = (dx * dot_dx + dy * dot_dy) / D

    dot_alpha = (dx * dot_dy - dy * dot_dx) / (D * D)
    dot_R = dot_alpha * np.array([[-sin_a, -cos_a, 0.0], [cos_a, -sin_a, 0.0], [0.0, 0.0, 0.0]])

    pos_se = state_se[:3]
    vel_se = state_se[3:]
    pos_rel_se = pos_se - np.array([1.0 - mu_SE, 0.0, 0.0])

    pos_rel_em = D * (R @ pos_rel_se)
    pos_em = pos_rel_em + np.array([-mu_EM, 0.0, 0.0])

    gamma = 1.0 + omega_S
    vel_rel_em = dot_D * (R @ pos_rel_se) + D * (dot_R @ pos_rel_se) + D * (R @ (gamma * vel_se))

    return np.concatenate([pos_em, vel_rel_em])


def se_lyapunov_to_bcr4bp_torus_seed(
    orbit_se: SymmetricOrbit,
    bcr_sys: bcr4bp.BCR4BPSystem,
    mu_SE: float,
    n_samples: int = 5,
) -> tuple[NDArray[np.float64], int, float]:
    """Sample a Sun-Earth L2 Lyapunov orbit and transform it to EM frame."""
    T_em = orbit_se.period / (1.0 + bcr_sys.omega_sun_nondim)

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit_se.period),
        np.array([orbit_se.x0, 0.0, 0.0, 0.0, orbit_se.ydot0, 0.0]),
        args=(mu_SE,),
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, orbit_se.period, n_samples, endpoint=False),
    )

    u_samples = np.zeros((n_samples, 6))
    for j in range(n_samples):
        state_se = sol.y[:, j]
        t_em = j * T_em / n_samples
        # Use mu = 0.0 for the starting system representation
        sys_mu0 = bcr4bp.BCR4BPSystem(
            mu=0.0,
            mu_sun=bcr_sys.mu_sun,
            a_sun_nondim=bcr_sys.a_sun_nondim,
            omega_sun_nondim=bcr_sys.omega_sun_nondim,
        )
        u_samples[j] = se_to_em_transform(state_se, t_em, sys_mu0, mu_SE)

    coeffs = np.fft.fft(u_samples, axis=0) / n_samples
    T_s = 2.0 * math.pi / bcr_sys.omega_sun_nondim
    rho = (2.0 * math.pi * T_s / T_em) % (2.0 * math.pi)
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


def bcr4bp_torus_residual(
    x_unk: NDArray[np.float64],
    system: bcr4bp.BCR4BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    amplitude_pin: float,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> NDArray[np.float64]:
    """Calculate the GMOS residual for BCR4BP torus."""
    coeffs_u, rho = x_unk[:-1], x_unk[-1]
    T_s = 2.0 * math.pi / system.omega_sun_nondim

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
        arc = bcr4bp.propagate_bcr4bp(
            system, u_s[j], T_s, with_stm=False, t0=0.0, rtol=rtol, atol=atol
        )
        phi_samples[j] = arc.state_f

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


def correct_bcr4bp_torus(
    system: bcr4bp.BCR4BPSystem,
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
) -> BCR4BPTorus:
    """Newton-correct a BCR4BP torus at a fixed stroboscopic period T_s."""
    res_gmos = least_squares(
        bcr4bp_torus_residual,
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

    T_s = 2.0 * math.pi / system.omega_sun_nondim
    omega_long = 2.0 * math.pi / T_s
    omega_trans = rho_final * omega_long / (2.0 * math.pi)

    return BCR4BPTorus(
        system=system,
        omega_long=omega_long,
        omega_trans=omega_trans,
        rho=rho_final,
        t_strob=T_s,
        fourier_coeffs=coeffs,
        n_modes=n_modes,
        n_samples=n_samples,
        invariance_residual=float(np.linalg.norm(res_gmos.fun)),
        converged=converged,
        n_iter=int(res_gmos.nfev),
    )


def evaluate_bcr4bp_torus(
    torus: BCR4BPTorus,
    theta_long: float,
    theta_trans: float,
) -> NDArray[np.float64]:
    """Return the state on the BCR4BP torus at (theta_long, theta_trans)."""
    theta_long_red = float(theta_long) % (2.0 * math.pi)
    theta_trans_red = float(theta_trans) % (2.0 * math.pi)
    dt = theta_long_red / torus.omega_long
    theta_seed = (theta_trans_red - torus.rho * theta_long_red / (2.0 * math.pi)) % (2.0 * math.pi)
    u0 = evaluate_invariant_circle(torus.fourier_coeffs, theta_seed)
    if dt == 0.0:
        return u0
    arc = bcr4bp.propagate_bcr4bp(torus.system, u0, dt, with_stm=False, t0=0.0)
    return arc.state_f
