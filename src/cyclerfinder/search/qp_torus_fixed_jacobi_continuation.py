"""Free-rho, fixed-Jacobi pseudo-arclength torus continuation + halo-family
rotation-number profiling (#615).

Why this module exists
----------------------
The `#548`->`#553`->`#555`->`#556`->`#612`->`#613` chain tried to reproduce Owen
& Baresi (2024)'s EM L1<->L2 quasi-halo heteroclinic pair: L1 latitudinal
frequency (rotation number) 0.2739 and L2 0.02163, both at Jacobi constant
C=3.15, mu=0.012153643. Every prior task explored the rotation number either
(a) at the near-bifurcation LINEAR estimate (holding rho fixed and varying
amplitude -- `#612`'s pin), or (b) versus C at the near-bifurcation estimate
(`#613`). A direct re-read of the paper surfaced two untested avenues, which
this module implements and `#615` used to characterise:

**Hypothesis 1 (this module's headline capability).** The paper's own Section
2.3 states that at a FIXED Jacobi constant, quasi-periodic tori "exist in
1-parameter families defined by their unique ratio of fundamental frequencies"
-- i.e. there is NOT one natural rotation number at C=3.15, there is a
continuum. `#612`'s corrector PINS the rotation number rho as a hard gauge row
(needed to stop the least-squares escaping to a degenerate omega2->0 branch),
so it can only walk amplitude at a caller-fixed rho. To follow the genuine
free-rho branch at fixed C, :func:`continue_qp_torus_fixed_jacobi` runs a
**pseudo-arclength continuation in the full pseudospectral solution space**
``(coeffs, omega1, omega2)``: rho is NOT pinned; instead the Jacobi constant is
held fixed by a constraint row and the family is walked by an arclength row. The
arclength constraint (not a rho pin) is what prevents the omega2->0 collapse --
that degenerate branch has a different local tangent, so stepping along the
genuine tangent stays on the physical family. The invariance-PDE residual and
its analytic Jacobian are reused UNCHANGED from
:mod:`cyclerfinder.search.variational_qp_torus`; only the *continuation
strategy* (which extra rows are appended and what is stepped) is new.

**Hypothesis 2.** The paper's illustrative Fig. 5 shows an L2 southern NRHO, so
the actual §4.1.1 orbits might be large-amplitude near-rectilinear halos rather
than the small-amplitude near-bifurcation halos every prior task assumed.
:func:`halo_family_rotation_profile` marches the symmetric-halo family in x0
into the large-amplitude / NRHO regime and records the Neimark-Sacker
rotation number at each member, so the orbit-type question can be checked
directly.

`#615`'s findings (both DECISIVE NEGATIVES; see ``data/OUTSTANDING.md``):
* H1: at C=3.15 the free-rho L1 quasi-halo family's rotation number is confined
  to |rho| in ~[0.061, 0.075]; as amplitude grows |rho| *decreases* (away from
  0.2739) to an amplitude FOLD near transverse-amp ~0.092, then turns back. The
  paper's 1-parameter family IS real (rho varies, it is not the flat pin
  `#612` imposed) but its range does not contain 0.2739 at C=3.15.
* H2: the halo family's Jacobi constant is MAXIMISED at its planar-Lyapunov
  bifurcation (L1 C~3.1745, L2 C~3.152), reached at ~zero amplitude; the
  large-amplitude / NRHO members sit at C well below 3.15 (L1 NRHO branch
  C~2.98-3.02). So no NRHO exists at C=3.15 -- O&B's C=3.15 quasi-halos cannot
  be NRHOs -- and along the family the rotation number never realises the
  (0.2739, 0.02163) pair at any single C.

Discipline
----------
No catalogue writeback: this is a capability/characterisation module, not a
discovery result. The returned tori are OUR computation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.nrho_continuation import SymmetricNRHO, correct_symmetric_nrho
from cyclerfinder.search.variational_qp_torus import (
    _N_STATE,
    _basis_matrices,
    _jacobian,
    _k2_first_harmonic_cols,
    _pack,
    _residual,
    _transverse_amplitude,
    _unpack,
    project_gmos_torus_to_2d,
)

# ---------------------------------------------------------------------------
# Jacobi constant of a 2D-Fourier torus (codebase convention C = 2*Omega - v^2).
# ---------------------------------------------------------------------------


def torus_jacobi_and_gradient(
    coeffs: NDArray[np.float64], n1: int, n2: int, mu: float
) -> tuple[float, NDArray[np.float64]]:
    """Jacobi constant at the reference torus angle ``(theta1, theta2) = (0, 0)``
    and its gradient w.r.t. the packed free vector ``z = (coeffs.ravel, w1, w2)``.

    Uses the codebase convention ``C = 2*Omega - v^2`` (twice the paper's
    ``C = Omega - v^2/2``), matching :func:`cyclerfinder.core.cr3bp.jacobi_constant`.
    At ``theta = (0, 0)`` only the constant/cosine coefficient block contributes
    (all sine basis functions vanish), which makes the gradient sparse and cheap.
    The ``[w1, w2]`` gradient entries are zero (C is velocity/position only).
    """
    a1, a2 = 2 * n1 + 1, 2 * n2 + 1
    mask1 = np.zeros(a1)
    mask1[: n1 + 1] = 1.0  # constant + cos block evaluate to 1 at theta1=0
    mask2 = np.zeros(a2)
    mask2[: n2 + 1] = 1.0
    outer = np.outer(mask1, mask2)  # (a1, a2)
    u = np.array([float(np.sum(coeffs[c] * outer)) for c in range(_N_STATE)])
    x, y, z, vx, vy, vz = u
    om = 1.0 - mu
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    d_omega_dx = x - om * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3
    d_omega_dy = y - om * y / r1**3 - mu * y / r2**3
    d_omega_dz = -om * z / r1**3 - mu * z / r2**3
    omega = 0.5 * (x * x + y * y) + om / r1 + mu / r2
    c_val = 2.0 * omega - (vx * vx + vy * vy + vz * vz)
    dc_dstate = 2.0 * np.array([d_omega_dx, d_omega_dy, d_omega_dz, -vx, -vy, -vz])
    n_coef = _N_STATE * a1 * a2
    grad = np.zeros(n_coef + 2)
    for c in range(_N_STATE):
        grad[c * a1 * a2 : (c + 1) * a1 * a2] = (dc_dstate[c] * outer).reshape(-1)
    return c_val, grad


# ---------------------------------------------------------------------------
# Fixed-Jacobi, free-rho pseudo-arclength continuation.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FixedJacobiTorusStep:
    """One accepted member of a fixed-Jacobi free-rho torus continuation.

    ``rotation_number = omega2 / omega1`` is FREE (not pinned) -- it is the
    quantity `#615` tracks. ``transverse_amplitude`` is the arclength-family
    coordinate (the k2=1 coefficient L2 norm). ``jacobi`` is held at the
    caller's target to ``|C - C_target| < 1e-4`` by the continuation's Jacobi
    constraint row. ``residual_rms`` is the RMS of the invariance-PDE residual
    (the physical closure quantity; the constraint rows are separate).
    """

    coeffs: NDArray[np.float64]
    omega1: float
    omega2: float
    rotation_number: float
    transverse_amplitude: float
    jacobi: float
    residual_rms: float
    n1: int
    n2: int


def _first_tangent(
    constraint_jac: NDArray[np.float64], amp_grad: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Unit null-vector of the constraint Jacobian (the family tangent),
    oriented toward increasing transverse amplitude."""
    _, _, vt = np.linalg.svd(constraint_jac, full_matrices=True)
    tangent = vt[-1]
    tangent = tangent / np.linalg.norm(tangent)
    if float(np.dot(tangent, amp_grad)) < 0.0:
        tangent = -tangent
    return tangent


def continue_qp_torus_fixed_jacobi(
    system: cr3bp.CR3BPSystem,
    seed_gmos_torus: QPTorus,
    jacobi_target: float,
    *,
    n1: int = 10,
    n2: int = 4,
    ds: float = 0.01,
    n_steps: int = 30,
    constraint_weight: float = 1.0e3,
    max_nfev: int = 120,
    pde_tol: float = 5e-5,
    jacobi_tol: float = 1e-4,
    rho_collapse_tol: float = 5e-3,
    max_shrink: int = 12,
) -> list[FixedJacobiTorusStep]:
    """Pseudo-arclength continuation of the quasi-halo torus family at FIXED
    Jacobi constant, with the rotation number FREE (Hypothesis 1).

    Bootstraps from a small-amplitude GMOS torus (projected onto the 2D Fourier
    basis, exactly as :func:`variational_qp_torus.discover_qp_torus`), then walks
    the family. The corrected system at each step is

        [ invariance-PDE residual        (reused UNCHANGED, no rho pin) ]
        [ longitudinal phase gauge                                     ]
        [ transverse phase gauge                                       ]
        [ Jacobi constraint   C(u) - jacobi_target                     ]
        [ pseudo-arclength    tangent . (z - z_prev) - ds              ]

    The 4 constraint rows are weighted by ``constraint_weight`` so they are not
    outvoted by the ~1500 PDE rows in the least-squares objective (found the
    hard way: an unweighted single Jacobi row lets the solver trade energy for
    PDE fit and drift off-shell / collapse to omega2->0). Steps that fail to
    converge, drift off-energy, or collapse rho toward zero are rejected and
    ``ds`` is halved; the continuation stops when ``ds`` underflows or the
    shrink budget is exhausted (typically the family's amplitude fold).

    Returns the list of accepted :class:`FixedJacobiTorusStep` (increasing
    arclength; it will turn an amplitude fold and begin retracing, which is the
    correct pseudo-arclength behaviour at the family's maximum amplitude).
    """
    mu = system.mu
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    coeffs0 = project_gmos_torus_to_2d(seed_gmos_torus, n1, n2, m1, m2)
    z = _pack(coeffs0, seed_gmos_torus.omega_long, seed_gmos_torus.omega_trans)

    t1 = 2 * np.pi * np.arange(m1) / m1
    t2 = 2 * np.pi * np.arange(m2) / m2
    p1, p1d = _basis_matrices(n1, t1)
    p2, p2d = _basis_matrices(n2, t2)
    cos1_col, sin1_col = _k2_first_harmonic_cols(n2)
    trans_by_coord = np.sqrt(
        np.sum(coeffs0[:, :, cos1_col] ** 2 + coeffs0[:, :, sin1_col] ** 2, axis=1)
    )
    phase2_coord = int(np.argmax(trans_by_coord))
    phase1_coord = 0
    sin1_col_t1 = n1 + 1
    a1, a2 = 2 * n1 + 1, 2 * n2 + 1
    n_pde = _N_STATE * m1 * m2

    def _coef_index(c: int, a: int, b: int) -> int:
        return c * a1 * a2 + a * a2 + b

    ph1_idx = _coef_index(phase1_coord, sin1_col_t1, 0)
    ph2_idx = _coef_index(phase2_coord, 0, sin1_col)

    # Reuse the PDE block of the existing residual/Jacobian with ALL gauge rows
    # zeroed (gauge_weight=0, rho_weight=0): we append our own constraint rows.
    base_args = (mu, n1, n2, p1, p1d, p2, p2d, phase1_coord, phase2_coord, 0.0, 0.0, 0.0, 0.0)

    def _amp_grad(zv: NDArray[np.float64]) -> NDArray[np.float64]:
        coeffs, _, _ = _unpack(zv, n1, n2)
        grad = np.zeros_like(zv)
        amp = _transverse_amplitude(coeffs, n2)
        if amp <= 0.0:
            return grad
        for c in range(_N_STATE):
            for a in range(a1):
                grad[_coef_index(c, a, cos1_col)] = coeffs[c, a, cos1_col] / amp
                grad[_coef_index(c, a, sin1_col)] = coeffs[c, a, sin1_col] / amp
        return grad

    def _residual_fn(
        zv: NDArray[np.float64],
        tangent: NDArray[np.float64],
        z_prev: NDArray[np.float64],
        step: float,
    ) -> NDArray[np.float64]:
        pde = _residual(zv, *base_args)[:n_pde]
        coeffs, _, _ = _unpack(zv, n1, n2)
        c_val, _ = torus_jacobi_and_gradient(coeffs, n1, n2, mu)
        cw = constraint_weight
        rows = np.array(
            [
                cw * coeffs[phase1_coord, sin1_col_t1, 0],
                cw * coeffs[phase2_coord, 0, sin1_col],
                cw * (c_val - jacobi_target),
                cw * (float(np.dot(zv - z_prev, tangent)) - step),
            ]
        )
        return np.concatenate([pde, rows])

    def _jac_fn(
        zv: NDArray[np.float64],
        tangent: NDArray[np.float64],
        z_prev: NDArray[np.float64],
        step: float,
    ) -> NDArray[np.float64]:
        jpde = _jacobian(zv, *base_args)[:n_pde]
        coeffs, _, _ = _unpack(zv, n1, n2)
        _, gc = torus_jacobi_and_gradient(coeffs, n1, n2, mu)
        cw = constraint_weight
        jac = np.zeros((n_pde + 4, zv.size))
        jac[:n_pde] = jpde
        jac[n_pde + 0, ph1_idx] = cw
        jac[n_pde + 1, ph2_idx] = cw
        jac[n_pde + 2, :] = cw * gc
        jac[n_pde + 3, :] = cw * tangent
        return jac

    # First tangent: null-space of the constraint Jacobian (PDE + phase + Jacobi).
    j0 = _jac_fn(z, np.zeros_like(z), z, 0.0)[:-1]
    tangent = _first_tangent(j0, _amp_grad(z))

    steps: list[FixedJacobiTorusStep] = []
    z_prev = z.copy()
    step_ds = ds
    accepted = 0
    shrinks = 0
    while accepted < n_steps:
        z_pred = z_prev + step_ds * tangent
        sol = least_squares(
            _residual_fn,
            z_pred,
            jac=_jac_fn,
            args=(tangent, z_prev, step_ds),
            method="trf",
            x_scale="jac",
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
            max_nfev=max_nfev,
        )
        z_new = sol.x
        coeffs, w1, w2 = _unpack(z_new, n1, n2)
        pde_rms = float(np.sqrt(np.sum(sol.fun[:n_pde] ** 2) / n_pde))
        c_val, _ = torus_jacobi_and_gradient(coeffs, n1, n2, mu)
        rot = w2 / w1 if w1 != 0.0 else float("nan")
        good = (
            pde_rms < pde_tol
            and abs(c_val - jacobi_target) < jacobi_tol
            and abs(rot) > rho_collapse_tol
            and math.isfinite(rot)
        )
        if not good:
            shrinks += 1
            step_ds *= 0.5
            if step_ds < 1e-4 or shrinks > max_shrink:
                break
            continue
        shrinks = 0
        steps.append(
            FixedJacobiTorusStep(
                coeffs=coeffs,
                omega1=w1,
                omega2=w2,
                rotation_number=rot,
                transverse_amplitude=_transverse_amplitude(coeffs, n2),
                jacobi=c_val,
                residual_rms=pde_rms,
                n1=n1,
                n2=n2,
            )
        )
        new_tan = z_new - z_prev
        nrm = float(np.linalg.norm(new_tan))
        if nrm < 1e-14:
            break
        new_tan = new_tan / nrm
        if float(np.dot(new_tan, tangent)) < 0.0:
            new_tan = -new_tan
        tangent = new_tan
        z_prev = z_new
        accepted += 1
    return steps


# ---------------------------------------------------------------------------
# Halo-family rotation-number profile into the large-amplitude / NRHO regime.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HaloFamilyRotationSample:
    """One member of a symmetric-halo family sweep: its IC, Jacobi constant,
    Neimark-Sacker center-pair rotation number (``|arg(lambda)|/2pi``, the
    linear/infinitesimal-amplitude estimate `#612`'s positive control matched to
    an actual GMOS torus), and Floquet spectral radius."""

    x0: float
    z0: float
    ydot0: float
    T_TU: float
    jacobi: float
    rotation_number: float | None
    spectral_radius: float


def _center_rotation(
    system: cr3bp.CR3BPSystem, x0: float, z0: float, ydot0: float, t_period: float
) -> tuple[float | None, float]:
    state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0])
    eigs = floquet_multipliers(monodromy(system, state0, t_period))
    spectral_radius = max(abs(e) for e in eigs)
    cands = [complex(e) for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    if not cands:
        return None, spectral_radius
    return abs(math.atan2(cands[0].imag, cands[0].real)) / (2 * math.pi), spectral_radius


def halo_family_rotation_profile(
    system: cr3bp.CR3BPSystem,
    seed: SymmetricNRHO,
    *,
    dx0: float,
    n_steps: int,
    z0_floor: float = 3e-3,
    max_shrink: int = 3,
) -> list[HaloFamilyRotationSample]:
    """March the symmetric-halo family in x0 (secant prediction) recording each
    member's Jacobi constant and Neimark-Sacker rotation number (Hypothesis 2).

    ``dx0`` sign selects the continuation direction; positive/negative reach the
    two family branches (toward the planar bifurcation vs. toward large
    amplitude / NRHOs). Stops cleanly when the corrector fails ``max_shrink``
    times in a row (halving ``dx0`` between retries) or the family returns to the
    planar orbit (``|z0| < z0_floor``). Each retained member carries a computed
    monodromy, so the rotation number is never optional.
    """
    samples: list[HaloFamilyRotationSample] = []
    prev: SymmetricNRHO | None = None
    x0, z0, ydot0, t_period = seed.x0, seed.z0, seed.ydot0, seed.T_TU
    step = dx0
    fails = 0
    for _ in range(n_steps):
        member = correct_symmetric_nrho(system, float(x0), z0, ydot0, t_period, with_monodromy=True)
        if not member.converged or member.monodromy is None or abs(member.z0) < z0_floor:
            fails += 1
            if fails >= max_shrink or prev is None:
                break
            step *= 0.5
            x0 = prev.x0 + step
            z0, ydot0, t_period = prev.z0, prev.ydot0, prev.T_TU
            continue
        fails = 0
        rot, spectral_radius = _center_rotation(
            system, member.x0, member.z0, member.ydot0, member.T_TU
        )
        samples.append(
            HaloFamilyRotationSample(
                x0=member.x0,
                z0=member.z0,
                ydot0=member.ydot0,
                T_TU=member.T_TU,
                jacobi=member.jacobi,
                rotation_number=rot,
                spectral_radius=spectral_radius,
            )
        )
        if prev is not None:
            denom = member.x0 - prev.x0
            if abs(denom) > 1e-14:
                z0 = member.z0 + (member.z0 - prev.z0) / denom * step
                ydot0 = member.ydot0 + (member.ydot0 - prev.ydot0) / denom * step
                t_period = member.T_TU + (member.T_TU - prev.T_TU) / denom * step
            else:
                z0, ydot0, t_period = member.z0, member.ydot0, member.T_TU
        else:
            z0, ydot0, t_period = member.z0, member.ydot0, member.T_TU
        prev = member
        x0 = member.x0 + step
    return samples


__all__ = [
    "FixedJacobiTorusStep",
    "HaloFamilyRotationSample",
    "continue_qp_torus_fixed_jacobi",
    "halo_family_rotation_profile",
    "torus_jacobi_and_gradient",
]
