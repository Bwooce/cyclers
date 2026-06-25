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
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    _enforce_reality,
    _gmos_residual,
    _pack_unknowns,
    _unpack_unknowns,
    evaluate_invariant_circle,
    is_practically_irrational,
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


def _z_scale(z: NDArray[np.float64], n_modes: int) -> NDArray[np.float64]:
    """Per-component characteristic scale of the augmented unknown vector.

    The augmented z mixes O(amplitude ~ 5e-4) Fourier-mode unknowns with
    O(1-6) rho / t_strob / C_J unknowns. A unit tangent in this raw, badly
    scaled space (scaled by ``ds``) perturbs the tiny modes by a HUGE relative
    amount, blowing up the truncation tail (|c_2|/|c_1| -> 0.3 in one ds=5e-3
    step, verified). Continuation is run in the NORMALISED coordinate
    ``w = z / scale`` so a single ``ds`` is a uniform fractional step across
    every variable; the mode block then moves by ``ds * scale_mode`` rather than
    ``ds``, keeping the torus shape stable along the family. The mode scale uses
    the constant-mode amplitude ``|c_1|`` as the characteristic mode magnitude;
    rho / t_strob / C_J use their own magnitudes (floored at 1).
    """
    coeffs, rho, t_strob, cj = _unpack_augmented(z, n_modes)
    mode_scale = max(float(np.linalg.norm(coeffs[1, :])), 1e-6)
    scale = np.empty_like(z)
    scale[: 6 + 12 * n_modes] = mode_scale
    scale[-3] = max(abs(rho), 1.0)
    scale[-2] = max(abs(t_strob), 1.0)
    scale[-1] = max(abs(cj), 1.0)
    return scale


def _energy_row_scale(c0: NDArray[np.float64], mu: float) -> float:
    """Row-scaling factor for the energy-tie constraint.

    The Jacobi gradient has near-Moon ``1/r2**3`` terms that make the energy
    row's gradient O(1e3-1e4) while every GMOS / phase row is O(1). Left
    unscaled, that single row drives the augmented Jacobian's condition number
    to ~1e8 (verified on the #290 smoke torus) and the Newton step explodes.
    Dividing the energy residual + its Jacobian row by ``||grad jacobi||``
    brings the row to O(1) and restores a well-conditioned solve. The reported
    ``C_J`` coordinate (read straight from ``z[-1]``) is unaffected by this
    purely numerical rescaling of the constraint equation.
    """
    g = _jacobi_state_grad(c0, mu)
    return max(float(np.linalg.norm(g)), 1.0)


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
    c0 = np.real(coeffs[0, :]).astype(np.float64)
    cj_state = jacobi_constant(c0, system.mu)
    scale = _energy_row_scale(c0, system.mu)
    parts.append(np.array([(cj_state - cj) / scale]))
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
        # Row-scaled by ||grad jacobi|| (see _energy_row_scale) so the residual
        # and Jacobian energy row stay consistent and the augmented solve is
        # well-conditioned.
        coeffs, _, _, _ = _unpack_augmented(z, n_modes)
        c0 = np.real(coeffs[0, :]).astype(np.float64)
        g = _jacobi_state_grad(c0, system.mu)
        scale = _energy_row_scale(c0, system.mu)
        jac[-1, :] = 0.0
        jac[-1, 0:6] = g / scale
        jac[-1, -1] = -1.0 / scale
        # Phase row (row -2): Im(c_1[phase_pin_idx]) -> 1.0 on that imag unknown.
        # c_1 imag block starts at offset 6 + 0*12 + 6 = 12 (n=1 mode).
        jac[-2, :] = 0.0
        jac[-2, 12 + phase_pin_idx] = 1.0
    return r0, jac


# ---------------------------------------------------------------------------
# Task 2: SVD null tangent + augmented arclength corrector.
# ---------------------------------------------------------------------------


def _arclength_tangent(
    jac: NDArray[np.float64],
    prev: NDArray[np.float64] | None,
    scale: NDArray[np.float64] | None = None,
) -> NDArray[np.float64] | None:
    """Unit null tangent = last right-singular vector of the residual Jacobian.

    Mirrors er3bp_continuation._arclength_tangent. ``jac`` is the residual
    Jacobian of shape ``(len(z)-1, len(z))`` -- rank ``len(z)-1`` with a
    one-dimensional null space. Oriented by ``+dot`` with ``prev``.

    When ``scale`` is given the tangent is the null vector of the SCALED
    Jacobian ``jac @ diag(scale)`` (the derivative w.r.t. the normalised
    coordinate ``w = z / scale``), returned unit-norm in w-space. ``prev`` is
    then also a w-space tangent. With ``scale=None`` the behaviour is the
    original raw-space tangent (kept for the Task-1/2 unit tests).
    """
    work = jac if scale is None else jac * scale[np.newaxis, :]
    try:
        _u, _s, vt = np.linalg.svd(work, full_matrices=True)
    except np.linalg.LinAlgError:
        return None
    if vt.shape[0] != work.shape[1]:
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
    scale: NDArray[np.float64] | None = None,
    max_iter: int = 60,
    mode_cap: float = 0.1,
    rho_cap: float = 0.05,
    cj_cap: float = 1e-2,
    residual_floor: float = 1e-5,
) -> NDArray[np.float64] | None:
    """Newton onto {R(z)=0, arclength=0}. Returns converged z or None.

    With ``scale`` given, ``tau`` is a w-space (normalised) tangent and the
    arclength row is ``tau . ((z - z_pred) / scale) = 0`` -- continuation runs
    in ``w = z / scale`` so the step is uniform across the badly-scaled mode vs
    energy/period unknowns (see _z_scale). With ``scale=None`` the constraint is
    the raw ``tau . (z - z_pred) = 0`` (Task-2 unit-test path).

    The GMOS invariance residual bottoms out at the FD-Jacobian / DOP853 noise
    floor (~4e-7 on the #290 smoke torus, the documented Phase-1 limit), so the
    convergence gate is ``max(tol, residual_floor)`` -- mirroring the legacy
    amplitude stub's ``max(tol, 1e-5)`` -- rather than the (unreachable) raw
    ``tol``. A tighter ``tol`` than the floor simply pins the gate at the floor.
    """
    gate = max(tol, residual_floor)
    inv_scale = None if scale is None else 1.0 / scale

    def arclength(zz: NDArray[np.float64]) -> float:
        if inv_scale is None:
            return float(np.dot(tau, zz - z_pred))
        return float(np.dot(tau, (zz - z_pred) * inv_scale))

    def arclength_grad() -> NDArray[np.float64]:
        return tau if inv_scale is None else tau * inv_scale

    def full_residual(zz: NDArray[np.float64]) -> NDArray[np.float64]:
        rr = _augmented_residual(zz, system, n_modes, n_samples, phase_pin_idx)
        return np.concatenate([rr, np.array([arclength(zz)])])

    z = z_pred.copy()
    f = full_residual(z)
    lam = 1e-3
    for _ in range(max_iter):
        if float(np.linalg.norm(f[:-1])) < gate and abs(f[-1]) < 1e-8:
            return z
        _, grad = _gmos_residual_and_jac(z, system, n_modes, n_samples, phase_pin_idx)
        jmat = np.vstack([grad, arclength_grad().reshape(1, -1)])
        # Levenberg-Marquardt damped Gauss-Newton: the augmented GMOS Jacobian
        # is intrinsically ill-conditioned (cond ~5e8 on the #290 smoke torus --
        # a soft mode direction), so a plain Newton solve amplifies noise and
        # diverges. LM damping (matching the trust-region behaviour of the
        # Phase-1 scipy.least_squares corrector) keeps the step stable; lambda
        # adapts down on accepted steps and up on rejected ones.
        jtj = jmat.T @ jmat
        jtf = jmat.T @ f
        diag = np.diag(jtj).copy()
        diag[diag <= 0] = 1.0
        accepted = False
        for _inner in range(12):
            try:
                dz = np.linalg.solve(jtj + lam * np.diag(diag), -jtf)
            except np.linalg.LinAlgError:
                dz, *_ = np.linalg.lstsq(jtj + lam * np.diag(diag), -jtf, rcond=None)
            dz = _apply_step_caps(
                np.asarray(dz, dtype=np.float64), n_modes, mode_cap, rho_cap, cj_cap
            )
            z_try = z + dz
            coeffs, rho, t_strob, cj = _unpack_augmented(z_try, n_modes)
            coeffs = _enforce_reality(coeffs)
            z_try = _pack_augmented(coeffs, rho, t_strob, cj)
            f_try = full_residual(z_try)
            if float(np.linalg.norm(f_try)) < float(np.linalg.norm(f)):
                z, f = z_try, f_try
                lam = max(lam * 0.5, 1e-12)
                accepted = True
                break
            lam = min(lam * 4.0, 1e8)
        if not accepted:
            return None
    return None


# ---------------------------------------------------------------------------
# Task 3: family-member records + single forward arclength step.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResonanceFlag:
    p: int
    q: int
    distance: float


@dataclass(frozen=True)
class ResonanceCrossing:
    member_index: int
    p: int
    q: int
    freq_ratio: float


@dataclass(frozen=True)
class QPTorusFold:
    member_index: int
    param_at_fold: float
    tangent_param_component: float


@dataclass(frozen=True)
class QPTorusFamilyMember:
    torus: QPTorus
    jacobi: float
    arclength_s: float
    tangent: NDArray[np.float64]
    rho: float
    freq_ratio: float
    is_practically_irrational: bool
    near_resonance: ResonanceFlag | None
    fold_index: int | None
    residual_norm: float
    extras: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class QPFamily:
    members: list[QPTorusFamilyMember]
    folds: list[QPTorusFold]
    resonance_crossings: list[ResonanceCrossing]
    terminated_reason: str
    seed_torus_id: str


def _nearest_rational(ratio: float, max_denominator: int) -> tuple[int, int, float]:
    """Closest p/q (|q| <= max_denominator) to ratio; returns (p, q, distance)."""
    from fractions import Fraction

    fr = Fraction(ratio).limit_denominator(max_denominator)
    return fr.numerator, fr.denominator, abs(ratio - float(fr))


def _independent_residual(
    coeffs: NDArray[np.complex128],
    rho: float,
    t_strob: float,
    system: cr3bp.CR3BPSystem,
    *,
    n_off_grid: int = 16,
) -> float:
    rng = np.random.default_rng(seed=0xC0FFEE)
    grid_thetas = 2 * math.pi * np.arange(2 * coeffs.shape[0]) / (2 * coeffs.shape[0])
    max_err = 0.0
    for _ in range(n_off_grid):
        theta = rng.uniform(0.0, 2 * math.pi)
        while np.any(np.abs(theta - grid_thetas) < 1e-6):
            theta = rng.uniform(0.0, 2 * math.pi)
        u0 = evaluate_invariant_circle(coeffs, theta)
        try:
            arc = cr3bp.propagate(system, u0, t_strob, with_stm=False)
        except RuntimeError:
            return float("inf")
        u_target = evaluate_invariant_circle(coeffs, theta + rho)
        max_err = max(max_err, float(np.linalg.norm(arc.state_f - u_target)))
    return max_err


def _member_from_z(
    z: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    *,
    tau: NDArray[np.float64],
    arclength_s: float,
    fold_index: int | None,
    resonance_max_denominator: int = 12,
    resonance_tol: float = 1e-4,
) -> QPTorusFamilyMember:
    coeffs, rho, t_strob, cj = _unpack_augmented(z, n_modes)
    coeffs = _enforce_reality(coeffs)
    omega_long = 2 * math.pi / t_strob
    omega_trans = rho / t_strob
    freq_ratio = omega_trans / omega_long if omega_long != 0.0 else float("nan")
    r = _augmented_residual(z, system, n_modes, n_samples, phase_pin_idx)
    residual_norm = float(np.linalg.norm(r))
    indep = _independent_residual(coeffs, rho, t_strob, system)
    irrational = is_practically_irrational(
        freq_ratio, max_denominator=resonance_max_denominator, tol=resonance_tol
    )
    p, q, dist = _nearest_rational(freq_ratio, resonance_max_denominator)
    near = ResonanceFlag(p, q, dist) if dist < resonance_tol else None
    torus = QPTorus(
        system=system,
        omega_long=omega_long,
        omega_trans=omega_trans,
        rho=rho,
        t_strob=t_strob,
        fourier_coeffs=coeffs,
        n_modes=n_modes,
        n_samples=n_samples,
        invariance_residual=residual_norm,
        independent_closure_residual=indep,
        converged=(residual_norm < 1e-5 and indep < 1e-4),
        n_iter=0,
        notes="333_arclength_member",
    )
    return QPTorusFamilyMember(
        torus=torus,
        jacobi=float(cj),
        arclength_s=float(arclength_s),
        tangent=tau.copy(),
        rho=float(rho),
        freq_ratio=float(freq_ratio),
        is_practically_irrational=bool(irrational),
        near_resonance=near,
        fold_index=fold_index,
        residual_norm=residual_norm,
        extras={"independent_residual": indep},
    )


# ---------------------------------------------------------------------------
# Task 4: both-directions arclength driver.
# ---------------------------------------------------------------------------


def _seed_augmented_z(
    seed_torus: QPTorus, system: cr3bp.CR3BPSystem
) -> tuple[NDArray[np.float64], int]:
    cj = jacobi_constant(np.real(seed_torus.fourier_coeffs[0, :]), system.mu)
    phase_pin_idx = int(np.argmax(np.abs(np.real(seed_torus.fourier_coeffs[1, :]))))
    z = _pack_augmented(seed_torus.fourier_coeffs, seed_torus.rho, seed_torus.t_strob, cj)
    return z, phase_pin_idx


def _tail_energy_ratio(coeffs: NDArray[np.complex128], n_modes: int) -> float:
    c1 = float(np.linalg.norm(coeffs[1, :]))
    c_n = float(np.linalg.norm(coeffs[n_modes, :]))
    return c_n / c1 if c1 > 0 else float("inf")


def _walk_direction(
    z0: NDArray[np.float64],
    tau0: NDArray[np.float64],
    sign: float,
    system: cr3bp.CR3BPSystem,
    n_modes: int,
    n_samples: int,
    phase_pin_idx: int,
    *,
    ds: float,
    max_steps: int,
    corrector_tol: float,
    fold_detection: bool,
    fold_eps: float,
    resonance_max_denominator: int,
    resonance_tol: float,
    mode_truncation_guard: float,
    on_step: Callable[[QPTorusFamilyMember], None] | None,
    folds: list[QPTorusFold],
    crossings: list[ResonanceCrossing],
    member_offset: int,
) -> tuple[list[QPTorusFamilyMember], str]:
    members: list[QPTorusFamilyMember] = []
    z_cur = z0.copy()
    tau = sign * tau0
    s_acc = 0.0
    reason = "max_steps"
    for _step in range(max_steps):
        # Normalised-coordinate continuation: scale is the per-component
        # characteristic magnitude of z_cur, tau is the w-space tangent, so the
        # physical predictor step is ds * scale * tau (uniform fractional step).
        scale = _z_scale(z_cur, n_modes)
        z_pred = z_cur + ds * scale * tau
        z_next = _correct_arclength_torus(
            z_pred,
            tau,
            system,
            n_modes=n_modes,
            n_samples=n_samples,
            phase_pin_idx=phase_pin_idx,
            tol=corrector_tol,
            scale=scale,
        )
        if z_next is None:
            z_next = _correct_arclength_torus(
                z_cur + 0.5 * ds * scale * tau,
                tau,
                system,
                n_modes=n_modes,
                n_samples=n_samples,
                phase_pin_idx=phase_pin_idx,
                tol=corrector_tol,
                scale=scale,
            )
            if z_next is None:
                reason = "corrector_fail"
                break
        _, jac_next = _gmos_residual_and_jac(z_next, system, n_modes, n_samples, phase_pin_idx)
        scale_next = _z_scale(z_next, n_modes)
        tau_next = _arclength_tangent(jac_next, tau, scale_next)
        if tau_next is None:
            reason = "corrector_fail"
            break
        fold_index = None
        # A genuine fold is a turning point in the C_J continuation coordinate:
        # the tangent's C_J component (tau[-1]) flips sign with non-trivial
        # magnitude. The smoke family barely moves in energy (|tau[-1]| ~ 3e-3,
        # noise level), so a raw sign-flip there fires spuriously every step.
        # Require the magnitude on at least one side to exceed fold_eps so only
        # a real energy turning point (where |tau[-1]| sweeps up before crossing
        # zero) is recorded -- the analytic-Jacobian de-noising of Phase 2.
        if (
            fold_detection
            and tau[-1] * tau_next[-1] < 0.0
            and max(abs(float(tau[-1])), abs(float(tau_next[-1]))) > fold_eps
        ):
            fold_index = member_offset + len(members)
            folds.append(QPTorusFold(fold_index, float(z_next[-1]), float(tau_next[-1])))
        s_acc += ds
        member = _member_from_z(
            z_next,
            system,
            n_modes,
            n_samples,
            phase_pin_idx,
            tau=tau_next,
            arclength_s=s_acc,
            fold_index=fold_index,
            resonance_max_denominator=resonance_max_denominator,
            resonance_tol=resonance_tol,
        )
        coeffs, _, _, _ = _unpack_augmented(z_next, n_modes)
        if _tail_energy_ratio(coeffs, n_modes) > mode_truncation_guard:
            members.append(member)
            if on_step is not None:
                on_step(member)
            reason = "mode_truncation_breach"
            break
        if member.near_resonance is not None and not member.is_practically_irrational:
            crossings.append(
                ResonanceCrossing(
                    member_offset + len(members),
                    member.near_resonance.p,
                    member.near_resonance.q,
                    member.freq_ratio,
                )
            )
            reason = "resonance_lock"
            members.append(member)
            if on_step is not None:
                on_step(member)
            break
        members.append(member)
        if on_step is not None:
            on_step(member)
        z_cur, tau = z_next, tau_next
    return members, reason


def continue_qp_family_arclength(
    seed_torus: QPTorus,
    *,
    param: Literal["jacobi", "rho"] = "jacobi",
    ds: float = 5e-3,
    max_steps: int = 200,
    direction: Literal["both", "fwd", "rev"] = "both",
    corrector_tol: float = 1e-8,
    phase_pin_idx: int | None = None,
    fold_detection: bool = True,
    fold_eps: float = 1e-2,
    resonance_max_denominator: int = 12,
    resonance_tol: float = 1e-4,
    mode_truncation_guard: float = 0.1,
    on_step: Callable[[QPTorusFamilyMember], None] | None = None,
) -> QPFamily:
    """Pseudo-arclength continuation of a converged QP 2-torus into a family.

    Predictor z_pred = z_cur + ds*tau (SVD null tangent); corrector
    _correct_arclength_torus. Walks BOTH directions by default. param="rho"
    monitors/drives the rotation number instead of energy (the resonance
    monitor fires off freq_ratio either way); "jacobi" (default) crosses Arnold
    tongues transversally. Report-only -- NO catalogue writeback.

    ``resonance_tol`` defaults to 1e-4 (not 1e-3): a thin Neimark-Sacker torus
    born off a k:1 bracket sits at ``freq_ratio = 1/k + O(amplitude^2)``. For
    the #290 smoke seed (k=4, amp=5e-4) the drift from 1/4 is ~3e-4, so a 1e-3
    band would falsely classify the genuine seed AND every family member as
    phase-locked. 1e-4 distinguishes that real torus drift from a true rational
    lock (which sits at ~1e-15 from p/q), matching the Phase-1 smoke test's
    "drift > 1e-6 => genuine torus, not the bifurcation periodic orbit" gate.
    """
    if param not in ("jacobi", "rho"):
        raise ValueError(f"param must be 'jacobi' or 'rho', got {param!r}")
    system = seed_torus.system
    n_modes = seed_torus.n_modes
    n_samples = seed_torus.n_samples
    z0, auto_pin = _seed_augmented_z(seed_torus, system)
    pin = auto_pin if phase_pin_idx is None else phase_pin_idx
    _, jac0 = _gmos_residual_and_jac(z0, system, n_modes, n_samples, pin)
    scale0 = _z_scale(z0, n_modes)
    tau0 = _arclength_tangent(jac0, None, scale0)
    if tau0 is None:
        return QPFamily([], [], [], "corrector_fail", seed_torus.notes)
    folds: list[QPTorusFold] = []
    crossings: list[ResonanceCrossing] = []
    seed_member = _member_from_z(
        z0,
        system,
        n_modes,
        n_samples,
        pin,
        tau=tau0,
        arclength_s=0.0,
        fold_index=None,
        resonance_max_denominator=resonance_max_denominator,
        resonance_tol=resonance_tol,
    )
    rev: list[QPTorusFamilyMember] = []
    fwd: list[QPTorusFamilyMember] = []
    reason = "max_steps"
    kw: dict[str, object] = dict(
        ds=ds,
        max_steps=max_steps,
        corrector_tol=corrector_tol,
        fold_detection=fold_detection,
        fold_eps=fold_eps,
        resonance_max_denominator=resonance_max_denominator,
        resonance_tol=resonance_tol,
        mode_truncation_guard=mode_truncation_guard,
        on_step=on_step,
    )
    if direction in ("both", "rev"):
        rev, r1 = _walk_direction(
            z0,
            tau0,
            -1.0,
            system,
            n_modes,
            n_samples,
            pin,
            folds=folds,
            crossings=crossings,
            member_offset=0,
            **kw,  # type: ignore[arg-type]
        )
        reason = r1
    if direction in ("both", "fwd"):
        fwd, r2 = _walk_direction(
            z0,
            tau0,
            +1.0,
            system,
            n_modes,
            n_samples,
            pin,
            folds=folds,
            crossings=crossings,
            member_offset=len(rev) + 1,
            **kw,  # type: ignore[arg-type]
        )
        reason = r2
    members = [*reversed(rev), seed_member, *fwd]
    return QPFamily(members, folds, crossings, reason, seed_torus.notes)
