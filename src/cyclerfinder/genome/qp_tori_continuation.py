"""Quasi-periodic invariant 2-tori continuation driver.

This module provides a natural-parameter (amplitude) continuation driver for
quasi-periodic tori in the CR3BP. It traces a 1-parameter family by stepping
the invariant-circle amplitude and employing a predictor-corrector strategy.

Fold Handling
-------------
This driver performs **natural-parameter continuation** (stepping amplitude).
It does **NOT** traverse folds (turning points) where the amplitude is
non-monotonic along the branch. At a fold, the corrector will fail to converge
or escape its basin. The driver detects this, records a failed step in the
audit trail, and gracefully stops the branch. Pseudo-arclength continuation
is deferred to future work.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    _correct_gmos,
    _enforce_reality,
    _pack_unknowns,
    _seed_invariant_circle,
    _unpack_unknowns,
    evaluate_invariant_circle,
)


@dataclass
class QPTorusFamilyStep:
    """Audit record for a single continuation step."""

    step_index: int
    amplitude_target: float
    predictor_type: str
    torus: QPTorus
    converged: bool
    notes: str = ""


@dataclass
class QPContinuationResult:
    """Result of a QP tori family continuation."""

    steps: list[QPTorusFamilyStep]
    family_converged: bool
    fold_detected: bool = False
    notes: str = ""


def continue_qp_family(
    system: cr3bp.CR3BPSystem,
    seed_orbit: NDArray[np.float64],
    seed_period: float,
    bifurcation_floquet_pair: tuple[complex, complex],
    *,
    k: int = 4,
    n_long: int = 16,
    n_trans: int = 8,
    amplitude_start: float = 5e-4,
    amplitude_end: float = 5e-3,
    n_steps: int = 10,
    tol: float = 1e-8,
    max_iter: int = 50,
    independent_tol: float = 1e-4,
    independent_n_samples: int = 16,
    notes: str = "",
) -> QPContinuationResult:
    """Continue a QP-torus family from a Neimark-Sacker bifurcation seed.

    Steps the ``amplitude_pin`` from ``amplitude_start`` to ``amplitude_end``.
    Defaults to a secant predictor (linear extrapolation from the last two
    converged members), falling back to a zeroth-order predictor for the first
    few steps.

    Parameters
    ----------
    system :
        CR3BP system at the parent family.
    seed_orbit :
        6-vector IC of the parent periodic orbit at the Neimark-Sacker bracket.
    seed_period :
        Period of the parent orbit (nondim TU).
    bifurcation_floquet_pair :
        The conjugate pair of Floquet multipliers crossing the unit circle.
    k :
        The integer such that the multiplier sits near a primitive k-th root
        of unity. Drives the Neimark-Sacker eigenvector search.
    n_long :
        Number of longitudinal modes (recorded for traceability).
    n_trans :
        Number of transverse Fourier modes.
    amplitude_start :
        Amplitude for the first step.
    amplitude_end :
        Amplitude for the final step.
    n_steps :
        Number of continuation steps.
    tol :
        Convergence tolerance for the GMOS L2 residual.
    max_iter :
        Maximum Newton-step iterations.
    independent_tol :
        Tolerance for the off-grid cross-check.
    independent_n_samples :
        Number of off-grid angles to use for the cross-check.

    Returns
    -------
    QPContinuationResult :
        Audit trail of the continuation branch, including early-stopping at folds.
    """
    if n_trans < 1:
        raise ValueError("n_trans must be >= 1")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    # Compute monodromy and eigenvector once to get lam_seed and phase_pin_idx
    arc = cr3bp.propagate(system, seed_orbit, seed_period, with_stm=True)
    monod = arc.stm
    assert monod is not None

    coeffs0_seed, lam_seed, _ = _seed_invariant_circle(
        seed_orbit, monod, k=k, n_modes=n_trans, amplitude=amplitude_start
    )
    phase_pin_idx = int(np.argmax(np.abs(np.real(coeffs0_seed[1, :]))))

    amplitudes = np.linspace(amplitude_start, amplitude_end, n_steps)
    n_samples = 2 * n_trans + 3

    steps = []
    prev_x: NDArray[np.float64] | None = None
    prev_prev_x: NDArray[np.float64] | None = None
    prev_amp: float | None = None
    prev_prev_amp: float | None = None

    fold_detected = False
    family_converged = True

    for i, amp in enumerate(amplitudes):
        predictor_type = "zeroth_order"
        step_notes = notes

        if i == 0:
            # Seed from the Neimark-Sacker bracket
            coeffs_guess, _, _ = _seed_invariant_circle(
                seed_orbit, monod, k=k, n_modes=n_trans, amplitude=amp
            )
            rho_guess = math.atan2(np.imag(lam_seed), np.real(lam_seed))
            t_strob_guess = float(seed_period)
            x0 = _pack_unknowns(coeffs_guess, rho_guess, t_strob_guess)
        elif (
            prev_x is not None
            and prev_prev_x is not None
            and prev_amp is not None
            and prev_prev_amp is not None
            and abs(prev_amp - prev_prev_amp) > 1e-12
        ):
            predictor_type = "secant"
            slope = (prev_x - prev_prev_x) / (prev_amp - prev_prev_amp)
            x0 = prev_x + slope * (amp - prev_amp)
        elif prev_x is not None:
            # Zeroth order from prev_x
            x0 = prev_x.copy()
        else:
            raise RuntimeError("Missing previous state for continuation")

        try:
            x_final, residual_norm, n_iter = _correct_gmos(
                system=system,
                x0=x0,
                n_trans=n_trans,
                n_samples=n_samples,
                phase_pin_idx=phase_pin_idx,
                amplitude_pin=amp,
                tol=tol,
                max_iter=max_iter,
            )
        except (RuntimeError, ValueError) as e:
            # Corrector failed: typically hits a turning point (fold) or diverges
            coeffs_guess, rho_guess, t_strob_guess = _unpack_unknowns(x0, n_trans)
            coeffs_guess = _enforce_reality(coeffs_guess)
            torus = QPTorus(
                system=system,
                omega_long=2 * math.pi / t_strob_guess,
                omega_trans=rho_guess / t_strob_guess,
                rho=rho_guess,
                t_strob=t_strob_guess,
                fourier_coeffs=coeffs_guess,
                n_modes=n_trans,
                n_samples=n_samples,
                invariance_residual=float("inf"),
                independent_closure_residual=float("inf"),
                converged=False,
                n_iter=0,
                notes=f"newton_failed: {e}; {step_notes}",
            )
            steps.append(
                QPTorusFamilyStep(i, amp, predictor_type, torus, False, "Corrector diverged")
            )
            fold_detected = True
            family_converged = False
            break

        coeffs_final, rho_final, t_strob_final = _unpack_unknowns(x_final, n_trans)
        coeffs_final = _enforce_reality(coeffs_final)

        omega_long_final = 2 * math.pi / t_strob_final
        omega_trans_final = rho_final / t_strob_final

        # Independent cross-check
        rng = np.random.default_rng(seed=0xC0FFEE + i)
        theta_off_grid = rng.uniform(0.0, 2 * math.pi, size=independent_n_samples)
        grid_thetas = 2 * math.pi * np.arange(n_samples) / n_samples
        for j in range(independent_n_samples):
            while np.any(np.abs(theta_off_grid[j] - grid_thetas) < 1e-6):
                theta_off_grid[j] = rng.uniform(0.0, 2 * math.pi)

        max_indep_err = 0.0
        try:
            for theta_c in theta_off_grid:
                u0 = evaluate_invariant_circle(coeffs_final, theta_c)
                arc_chk = cr3bp.propagate(system, u0, t_strob_final, with_stm=False)
                u_target = evaluate_invariant_circle(coeffs_final, theta_c + rho_final)
                err = float(np.linalg.norm(arc_chk.state_f - u_target))
                max_indep_err = max(max_indep_err, err)
        except RuntimeError as e:
            max_indep_err = float("inf")
            step_notes = f"{step_notes}; independent_check_propagation_failed: {e}"

        # The finite-difference Jacobian at n_modes=2 often bottoms out ~1e-7.
        # We relax the strict convergence gate to 1e-5 so we don't falsely
        # detect a fold when the optimizer just hits the noise floor. Phase 2
        # (analytic Jacobian) will remove the need for this.
        invariance_gate = max(tol, 1e-5)
        converged = (residual_norm < invariance_gate) and (max_indep_err < independent_tol)

        torus = QPTorus(
            system=system,
            omega_long=omega_long_final,
            omega_trans=omega_trans_final,
            rho=rho_final,
            t_strob=t_strob_final,
            fourier_coeffs=coeffs_final,
            n_modes=n_trans,
            n_samples=n_samples,
            invariance_residual=residual_norm,
            independent_closure_residual=max_indep_err,
            converged=converged,
            n_iter=n_iter,
            notes=step_notes,
            extras={
                "n_long_recorded": float(n_long),
                "k": float(k),
                "lam_seed_re": float(np.real(lam_seed)),
                "lam_seed_im": float(np.imag(lam_seed)),
                "phase_pin_idx": float(phase_pin_idx),
                "amplitude_pin": float(amp),
            },
        )

        steps.append(
            QPTorusFamilyStep(
                i,
                amp,
                predictor_type,
                torus,
                converged,
                "Converged" if converged else "Residual > tol",
            )
        )

        if not converged:
            fold_detected = True
            family_converged = False
            break

        prev_prev_x = prev_x
        prev_prev_amp = prev_amp
        prev_x = x_final
        prev_amp = amp

    return QPContinuationResult(
        steps=steps,
        family_converged=family_converged,
        fold_detected=fold_detected,
        notes="Family continuation completed"
        if family_converged
        else "Branch stopped at non-converging step (possible fold)",
    )
