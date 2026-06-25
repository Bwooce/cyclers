"""Pseudo-arclength continuation of GMOS quasi-periodic 2-tori (#333 / #290 Phase 2).

Mirrors the proven ER3BP arclength walker (genome/er3bp_continuation.py):
augmented unknowns z = [pack_unknowns(modes, rho, t_strob), C_J], SVD null
tangent prediction, augmented-Jacobian Newton correction onto the GMOS residual
plus an energy-tie row plus the arclength constraint. Replaces the fold-blind
amplitude stub (genome/qp_tori_continuation.py), which is superseded not deleted.

Report-only: NO catalogue writeback, NO novelty claims. The method is
Olikara-Scheeres 2010 / Olikara 2016 GMOS continuation; the tori returned are
OUR computation.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.genome.qp_tori import (
    _enforce_reality,
    _gmos_residual,
    _pack_unknowns,
    _unpack_unknowns,
)


def _pack_augmented(
    coeffs: NDArray[np.complex128], rho: float, t_strob: float, cj: float
) -> NDArray[np.float64]:
    x = _pack_unknowns(coeffs, rho, t_strob)
    return np.concatenate([x, [float(cj)]])


def _unpack_augmented(
    z: NDArray[np.float64], n_modes: int
) -> tuple[NDArray[np.complex128], float, float, float]:
    coeffs, rho, t_strob = _unpack_unknowns(z[:-1], n_modes)
    return coeffs, rho, t_strob, float(z[-1])


def _residual_size(n_modes: int) -> int:
    """GMOS rows (6 + 4*6*N) + phase pin (1) + energy tie (1)."""
    return 6 + 4 * 6 * n_modes + 2


def _augmented_residual(
    z: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
) -> NDArray[np.float64]:
    coeffs, rho, t_strob, cj = _unpack_augmented(z, n_modes)
    if t_strob <= 0 or not math.isfinite(t_strob) or not math.isfinite(rho):
        return np.full(_residual_size(n_modes), 1e10)
    try:
        f_res = _gmos_residual(system, coeffs, rho, t_strob, n_samples=n_samples)
    except (RuntimeError, ValueError):
        return np.full(_residual_size(n_modes), 1e10)
    parts: list[NDArray[np.float64]] = [np.real(f_res[0, :])]
    n_total_sig = f_res.shape[0]
    for n in range(1, n_modes + 1):
        parts.append(np.real(f_res[n, :]))
        parts.append(np.imag(f_res[n, :]))
        parts.append(np.real(f_res[n_total_sig - n, :]))
        parts.append(np.imag(f_res[n_total_sig - n, :]))
    parts.append(np.array([float(np.imag(coeffs[1, phase_pin_idx]))]))
    cj_state = jacobi_constant(np.real(coeffs[0, :]), system.mu)
    parts.append(np.array([cj_state - cj]))
    return np.concatenate(parts)


def _jacobi_state_grad(state6: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Analytic gradient of jacobi_constant w.r.t. the 6 state components.

    C = (x^2+y^2) + 2(1-mu)/r1 + 2 mu/r2 - (vx^2+vy^2+vz^2),
    r1=||(x+mu,y,z)||, r2=||(x-1+mu,y,z)||.
    """
    x, y, zc, vx, vy, vz = (float(v) for v in state6)
    dx1, dx2 = x + mu, x - 1.0 + mu
    r1 = math.sqrt(dx1 * dx1 + y * y + zc * zc)
    r2 = math.sqrt(dx2 * dx2 + y * y + zc * zc)
    a1 = -2.0 * (1.0 - mu) / r1**3
    a2 = -2.0 * mu / r2**3
    gx = 2.0 * x + a1 * dx1 + a2 * dx2
    gy = 2.0 * y + a1 * y + a2 * y
    gz = a1 * zc + a2 * zc
    return np.array([gx, gy, gz, -2.0 * vx, -2.0 * vy, -2.0 * vz], dtype=np.float64)


def _gmos_residual_and_jac(
    z: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    *,
    analytic: bool = True,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    r0 = _augmented_residual(z, system, n_modes, n_samples, phase_pin_idx)
    m, n = r0.size, z.size
    jac = np.zeros((m, n))
    # One-sided FD column-by-column. The "analytic" path reuses the same FD
    # mechanism for the GMOS block (per-sample propagation already dominates
    # cost); the energy + phase rows below are overwritten closed-form so the
    # analytic Jacobian is exact on those structurally-known rows and FD-parity
    # holds on the GMOS block. (A full STM-variational GMOS block is a future
    # speedup; this keeps Phase-2 correctness-first while removing the
    # least_squares diff_step coupling that masked folds in Phase 1.)
    step = 1e-7
    for j in range(n):
        dz = z.copy()
        h = step * max(1.0, abs(z[j]))
        dz[j] += h
        rj = _augmented_residual(dz, system, n_modes, n_samples, phase_pin_idx)
        jac[:, j] = (rj - r0) / h
    if analytic:
        # Energy row (row -1): d/d(c_0 real components) = grad jacobi; d/dC_J = -1.
        coeffs, _, _, _ = _unpack_augmented(z, n_modes)
        c0 = np.real(coeffs[0, :]).astype(np.float64)
        g = _jacobi_state_grad(c0, system.mu)
        jac[-1, :] = 0.0
        jac[-1, 0:6] = g
        jac[-1, -1] = -1.0
        # Phase row (row -2): Im(c_1[phase_pin_idx]) -> 1.0 on that imag unknown.
        # c_1 imag block starts at offset 6 + 0*12 + 6 = 12 (n=1 mode).
        jac[-2, :] = 0.0
        jac[-2, 12 + phase_pin_idx] = 1.0
    return r0, jac


# ---------------------------------------------------------------------------
# Task 2: SVD null tangent + augmented arclength corrector.
# ---------------------------------------------------------------------------


def _arclength_tangent(
    jac: NDArray[np.float64], prev: NDArray[np.float64] | None
) -> NDArray[np.float64] | None:
    """Unit null tangent = last right-singular vector of the residual Jacobian.

    Mirrors er3bp_continuation._arclength_tangent. ``jac`` is the residual
    Jacobian of shape ``(len(z)-1, len(z))`` -- rank ``len(z)-1`` with a
    one-dimensional null space. Oriented by ``+dot`` with ``prev``.
    """
    try:
        _u, _s, vt = np.linalg.svd(jac, full_matrices=True)
    except np.linalg.LinAlgError:
        return None
    if vt.shape[0] != jac.shape[1]:
        return None
    tau = np.asarray(vt[-1], dtype=np.float64)
    norm = float(np.linalg.norm(tau))
    if not np.all(np.isfinite(tau)) or norm < 1e-12:
        return None
    tau = tau / norm
    if prev is not None and float(np.dot(tau, prev)) < 0.0:
        tau = -tau
    return tau


def _apply_step_caps(
    dz: NDArray[np.float64], n_modes: int, mode_cap: float, rho_cap: float, cj_cap: float
) -> NDArray[np.float64]:
    out = dz.copy()
    out[: 6 + 12 * n_modes] = np.clip(out[: 6 + 12 * n_modes], -mode_cap, mode_cap)
    out[-3] = float(np.clip(out[-3], -rho_cap, rho_cap))  # rho
    # t_strob (out[-2]) capped at mode_cap-scale; cj (out[-1]) at cj_cap
    out[-2] = float(np.clip(out[-2], -mode_cap, mode_cap))
    out[-1] = float(np.clip(out[-1], -cj_cap, cj_cap))
    return out


def _correct_arclength_torus(
    z_pred: NDArray[np.float64],
    tau: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    *,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    tol: float,
    max_iter: int = 60,
    mode_cap: float = 0.1,
    rho_cap: float = 0.05,
    cj_cap: float = 1e-2,
) -> NDArray[np.float64] | None:
    """Newton onto {R(z)=0, tau.(z - z_pred)=0}. Returns converged z or None."""
    z = z_pred.copy()
    for _ in range(max_iter):
        r0, grad = _gmos_residual_and_jac(z, system, n_modes, n_samples, phase_pin_idx)
        arc = float(np.dot(tau, z - z_pred))
        if float(np.linalg.norm(r0)) < tol and abs(arc) < 1e-10:
            return z
        jmat = np.vstack([grad, tau.reshape(1, -1)])
        rhs = -np.concatenate([r0, np.array([arc])])
        try:
            dz = np.linalg.solve(jmat, rhs)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jmat, rhs, rcond=None)
        dz = _apply_step_caps(np.asarray(dz, dtype=np.float64), n_modes, mode_cap, rho_cap, cj_cap)
        z = z + dz
        coeffs, rho, t_strob, cj = _unpack_augmented(z, n_modes)
        coeffs = _enforce_reality(coeffs)
        z = _pack_augmented(coeffs, rho, t_strob, cj)
    return None
