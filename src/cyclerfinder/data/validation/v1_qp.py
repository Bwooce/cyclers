"""V1 same-model QP-torus gauntlet (#319 Phase 1 Part A).

Spec reinterpretation for QP-tori
---------------------------------
The strict-periodic V1 in :mod:`v1_3d` asserts ``||X(T) - X(0)|| < 1 m/s`` as
the same-model closure gate. A quasi-periodic 2-torus does NOT close in that
sense: it is invariant under the stroboscopic flow modulo a rotation by the
rotation number ``rho`` on the invariant circle. Strict ``X(T) - X(0) = 0``
holds nowhere on the torus (except the degenerate phase-locked limit).

For QP-tori the same-model gate is therefore **Fourier-mode invariance**::

    ||T(theta + omega * Delta_t) - phi(T(theta), Delta_t)||  <  V1_qp_floor

where ``T(theta)`` is the parameterized torus and ``phi`` the natural flow.

For the invariant-circle representation in :mod:`cyclerfinder.genome.qp_tori`
this reduces to two numbers already emitted by the corrector:

  * ``invariance_residual`` -- L2 norm of the GMOS residual in Fourier-mode
    space after Newton convergence; the "is this a closed torus in the
    truncation we chose?" gate.
  * ``independent_closure_residual`` -- maximum over OFF-GRID sample angles
    of ``||phi_{t_strob}(u(theta)) - u(theta + rho)||``; the independent
    cross-check that guards against the corrector matching aliased Fourier
    modes rather than a genuine torus.

V1_qp wraps both into a single verdict, plus runs ANOTHER round of off-grid
samples with a different RNG seed (independent of the corrector's own check)
so the audit trail records a fresh confirmation rather than re-quoting the
corrector's number.

Floor rationale (sourced)
-------------------------
Olikara & Howell 2014 ("Quasi-periodic invariant tori of an inertially
oriented spacecraft", Howell 2014 / Olikara dissertation 2016) report GMOS
invariance residuals in the 1e-4 to 1e-6 band depending on truncation order
``N`` and rotation-number value (small denominators in the continued fraction
of ``rho/(2*pi)`` push the residual up). Olikara 2016 §3.3 documents the
truncation-error scaling as ``O(|c_{N+1}|)``, the magnitude of the first
truncated Fourier mode.

We adopt **1e-5 as the V1_qp Fourier-norm floor** -- one order more permissive
than periodic V1's 1e-6, matching the worst case Olikara reports for the
moderate truncation orders the project uses (``N=2..8``).

The independent cross-check floor 1e-4 nondim mirrors the QPTorus default
``independent_tol`` (the corrector's own off-grid check; Olikara 2016 §4
recommends this).

What V1_qp is NOT
-----------------
* It does NOT re-solve the GMOS system from scratch -- it accepts a converged
  :class:`QPTorus` and re-verifies the residual at the stored frequencies.
  Re-running Newton from a different seed is a TIER-PROMOTION question
  (V2 / V3), not a V1 question.
* It does NOT assert km/s units. QP-tori, unlike strictly periodic cyclers,
  do not have a meaningful "closure delta-v" -- the residual lives in
  Fourier-mode coefficient space, which has units of state-vector norm
  rather than velocity. The km equivalent of the off-grid residual is
  REPORTED for human triage but NOT gated against the spec km/s floor.

Discipline
----------
* READ-ONLY on :mod:`cyclerfinder.genome.qp_tori` (wrap, don't modify).
* The floors are sourced module constants; an override is for caller
  reflexivity only (e.g. tightening to 1e-6 if a higher truncation order
  ``N`` justifies it).
* Independent cross-check uses a DIFFERENT RNG seed from the corrector's
  internal check -- the corrector seeds at ``0xC0FFEE``; we seed at
  ``0xDECAFBAD`` so the two checks are statistically independent.
* NO catalogue writeback. V1_qp pass alone does NOT admit a torus to the
  catalogue as a ``quasi_cycler`` -- V2 + #328 lit-check follow.

References
----------
* Olikara, Z. (2016). "Computation of Quasi-Periodic Tori and Heteroclinic
  Connections in Astrodynamics." PhD dissertation, Purdue University.
* Howell, K. C., & Howell, A. R. (2014). "Survey of Quasi-Periodic Motion
  in Cislunar Space for Transfer Design." AAS Astrodynamics Specialist.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, evaluate_invariant_circle

V1_QP_FOURIER_FLOOR: Final[float] = 1.0e-5
"""V1 same-model Fourier-mode invariance floor for QP-tori (nondim).

One order more permissive than periodic V1 (1e-6); per Olikara-Howell 2014
the residual saturates at ``O(|c_{N+1}|)`` from Fourier truncation. The
sourced practice for ``N=2..8`` lands in the 1e-4 to 1e-6 band; 1e-5 sits
in the middle and matches the floor adopted by Olikara 2016 §3.3 examples."""

V1_QP_INDEPENDENT_FLOOR: Final[float] = 1.0e-4
"""V1 independent off-grid cross-check floor (nondim).

Matches the default ``independent_tol`` in :func:`correct_qp_torus`. The
off-grid check is the corrector-independent topology gate -- a corrector
that converged in Fourier-mode space might still have aliasing artefacts;
the off-grid check refuses those."""

V1_QP_INDEPENDENT_N_SAMPLES: Final[int] = 100
"""Default number of off-grid angles used by the independent cross-check.

Olikara 2016 §4 recommends >=16 samples for the off-grid check; 100 gives
~6x oversampling, well below the cost of the corrector itself."""


@dataclass(frozen=True)
class V1VerdictQP:
    """Frozen V1 verdict for a QP-torus candidate.

    Attributes
    ----------
    candidate_id :
        Identifier of the candidate (carried for the audit trail).
    invariance_residual_fourier_norm :
        L2 norm of the GMOS residual in Fourier-mode space, RE-EVALUATED at
        the torus's stored frequencies. Equals
        :attr:`QPTorus.invariance_residual` for a freshly converged torus;
        will differ if the torus was loaded from disk after a bit-flip or
        was corrupted between Newton convergence and V1 evaluation.
    independent_invariance_residual_nondim :
        Maximum off-grid invariance error: for ``N`` random angles NOT on
        the corrector's FFT grid, the maximum of ``||phi_{t_strob}(u(theta))
        - u(theta + rho)||`` (nondim). Uses an INDEPENDENT RNG seed
        (``0xDECAFBAD``) from the corrector's own off-grid check
        (``0xC0FFEE``) -- so the V1 verdict carries a statistically
        independent confirmation rather than re-quoting the corrector.
    independent_residual_km :
        ``independent_invariance_residual_nondim`` converted to km via the
        system's length unit. Reported for human triage; NOT gated.
    v1_floor_invariance :
        The Fourier-norm floor this verdict was held against. Default
        :data:`V1_QP_FOURIER_FLOOR`.
    independent_floor_nondim :
        The off-grid nondim floor this verdict was held against. Default
        :data:`V1_QP_INDEPENDENT_FLOOR`.
    converged_corrector :
        ``True`` iff ``invariance_residual_fourier_norm <= v1_floor_invariance``.
    converged_independent :
        ``True`` iff ``independent_invariance_residual_nondim <=
        independent_floor_nondim``.
    passes_v1_qp :
        ``converged_corrector AND converged_independent``. The headline
        boolean.
    n_off_grid_samples :
        Number of off-grid angles used. Recorded for audit.
    n_modes :
        Fourier truncation order ``N`` carried from the torus. The floor
        rationale depends on this (Olikara-Howell 2014); a small ``N``
        justifies the relaxed floor.
    notes :
        Free-form audit string.
    """

    candidate_id: str
    invariance_residual_fourier_norm: float
    independent_invariance_residual_nondim: float
    independent_residual_km: float
    v1_floor_invariance: float
    independent_floor_nondim: float
    converged_corrector: bool
    converged_independent: bool
    passes_v1_qp: bool
    n_off_grid_samples: int
    n_modes: int
    notes: str = ""


def _off_grid_thetas(
    n_samples: int,
    grid_thetas: np.ndarray,
    *,
    rng_seed: int,
    avoid_tol: float = 1e-6,
) -> np.ndarray:
    """Draw ``n_samples`` random angles in ``[0, 2*pi)`` that are NOT within
    ``avoid_tol`` of any value in ``grid_thetas``.

    Distinct RNG seed from the corrector's own off-grid draw so the V1
    independent check is statistically independent of what the corrector
    saw at convergence time.
    """
    rng = np.random.default_rng(seed=rng_seed)
    out = rng.uniform(0.0, 2 * math.pi, size=n_samples)
    for i in range(n_samples):
        # Retry up to a generous cap; for n_samples=100 and a 7-point grid
        # the probability of collision is ~7e-5 per draw, so retries are
        # rare but possible.
        retries = 0
        while np.any(np.abs(out[i] - grid_thetas) < avoid_tol):
            out[i] = rng.uniform(0.0, 2 * math.pi)
            retries += 1
            if retries > 100:
                raise RuntimeError(
                    f"off-grid sampler stuck; cannot avoid grid_thetas with "
                    f"n_grid={len(grid_thetas)} after 100 retries"
                )
    return out


def _evaluate_invariance_residual_at_grid(
    torus: QPTorus,
) -> float:
    """Re-evaluate the GMOS invariance residual in Fourier-mode space at the
    torus's stored frequencies.

    For a freshly converged torus this returns ``torus.invariance_residual``
    (modulo integrator non-determinism). For a corrupted torus the two will
    differ -- the V1 check catches a torus that has been mutated since
    convergence.

    The residual is computed identically to :func:`_gmos_residual` in
    :mod:`qp_tori`, but inlined here because that helper is private and
    we MUST NOT modify the genome module per the discipline.
    """
    coeffs = torus.fourier_coeffs
    rho = torus.rho
    t_strob = torus.t_strob
    n_samples = torus.n_samples
    n_total = coeffs.shape[0]
    n_modes = (n_total - 1) // 2
    if n_samples < n_total:
        # Defensive: shouldn't happen for a torus that came from
        # correct_qp_torus, but guard against tampered tori.
        raise ValueError(f"torus has n_samples={n_samples} < 2 N + 1 = {n_total}")
    thetas = 2 * math.pi * np.arange(n_samples) / n_samples
    u_samples = evaluate_invariant_circle(coeffs, thetas)
    phi_samples = np.zeros_like(u_samples)
    for j in range(n_samples):
        arc = cr3bp.propagate(torus.system, u_samples[j], t_strob, with_stm=False)
        phi_samples[j] = arc.state_f
    phi_fft = np.fft.fft(phi_samples, axis=0) / n_samples
    expected = np.zeros((n_samples, 6), dtype=np.complex128)
    expected[0, :] = coeffs[0, :]
    for n in range(1, n_modes + 1):
        expected[n, :] = coeffs[n, :] * np.exp(1j * n * rho)
    for n in range(1, n_modes + 1):
        expected[n_samples - n, :] = coeffs[n_total - n, :] * np.exp(-1j * n * rho)
    residual = phi_fft - expected
    # Mask out the high-frequency tail (modes |n| > n_modes); they are
    # truncation indicators, not solver state. Matches qp_tori._gmos_residual.
    if n_samples > n_total:
        for n in range(n_modes + 1, n_samples - n_modes):
            residual[n, :] = 0.0
    return float(np.linalg.norm(residual))


def _independent_off_grid_residual(
    torus: QPTorus,
    *,
    n_samples: int,
    rng_seed: int,
) -> float:
    """Independent off-grid invariance check:
    Pick ``n_samples`` random angles NOT on the FFT grid, propagate each by
    ``t_strob``, and verify the result lies on the resampled invariant
    circle to within the Olikara-recommended off-grid floor.

    Returns the MAX (not mean) over the samples -- the V1 gate is L_infinity
    in topology checks, matching the QPTorus corrector's own discipline.
    """
    n_grid = torus.n_samples
    grid_thetas = 2 * math.pi * np.arange(n_grid) / n_grid
    off_thetas = _off_grid_thetas(n_samples, grid_thetas, rng_seed=rng_seed)
    max_err = 0.0
    for theta in off_thetas:
        u0 = evaluate_invariant_circle(torus.fourier_coeffs, theta)
        try:
            arc = cr3bp.propagate(torus.system, u0, torus.t_strob, with_stm=False)
        except RuntimeError:
            return float("inf")
        u_target = evaluate_invariant_circle(torus.fourier_coeffs, theta + torus.rho)
        err = float(np.linalg.norm(arc.state_f - u_target))
        if err > max_err:
            max_err = err
    return max_err


def run_v1_qp(
    candidate_id: str,
    torus: QPTorus,
    *,
    invariance_floor: float = V1_QP_FOURIER_FLOOR,
    independent_tol_nondim: float = V1_QP_INDEPENDENT_FLOOR,
    n_off_grid_samples: int = V1_QP_INDEPENDENT_N_SAMPLES,
    rng_seed: int = 0xDECAFBAD,
    notes: str = "",
) -> V1VerdictQP:
    """Run the V1 same-model QP-torus gauntlet.

    Pipeline:
      1. Re-evaluate the torus's GMOS invariance residual at the stored
         frequencies (no Newton iteration; pure residual check).
      2. Run an INDEPENDENT off-grid sample-point invariance test with a
         distinct RNG seed from the corrector's own check.
      3. Gate the Fourier-norm residual against ``invariance_floor``
         (default 1e-5; Olikara-Howell 2014).
      4. Gate the off-grid residual against ``independent_tol_nondim``
         (default 1e-4; matches QPTorus default).
      5. PASS iff both gates hold.

    Parameters
    ----------
    candidate_id :
        Identifier carried into the verdict for the audit trail.
    torus :
        A :class:`QPTorus` from :mod:`cyclerfinder.genome.qp_tori`. Must
        carry a converged state (a torus with ``invariance_residual=inf``
        will trivially fail, which is the correct outcome).
    invariance_floor :
        Fourier-norm floor. Default :data:`V1_QP_FOURIER_FLOOR` (1e-5).
    independent_tol_nondim :
        Off-grid nondim floor. Default :data:`V1_QP_INDEPENDENT_FLOOR`
        (1e-4).
    n_off_grid_samples :
        Number of off-grid angles. Default
        :data:`V1_QP_INDEPENDENT_N_SAMPLES` (100).
    rng_seed :
        Seed for the off-grid sampler. Default ``0xDECAFBAD`` (deliberately
        distinct from the corrector's ``0xC0FFEE`` so the V1 check is
        independent of the corrector's own off-grid check).
    notes :
        Free-form audit note.

    Returns
    -------
    V1VerdictQP
        The frozen verdict. ``passes_v1_qp`` is the headline boolean.

    Notes
    -----
    A V1_qp PASS does NOT admit to the catalogue as a ``quasi_cycler``
    (v4.7 ``orbit_class``). The orbit-closure discipline mandates V2 +
    #328 lit-check follow.
    """
    if not isinstance(torus, QPTorus):
        raise TypeError(f"torus must be a QPTorus instance; got {type(torus).__name__}")
    if invariance_floor <= 0.0:
        raise ValueError(f"invariance_floor must be > 0; got {invariance_floor}")
    if independent_tol_nondim <= 0.0:
        raise ValueError(f"independent_tol_nondim must be > 0; got {independent_tol_nondim}")
    if n_off_grid_samples < 1:
        raise ValueError(f"n_off_grid_samples must be >= 1; got {n_off_grid_samples}")
    if torus.system.l_km <= 0.0 or torus.system.t_s <= 0.0:
        raise ValueError(
            f"invalid CR3BP system for V1-qp km conversion: "
            f"l_km={torus.system.l_km} t_s={torus.system.t_s}"
        )

    invariance_norm = _evaluate_invariance_residual_at_grid(torus)
    indep_nondim = _independent_off_grid_residual(
        torus,
        n_samples=n_off_grid_samples,
        rng_seed=rng_seed,
    )
    indep_km = indep_nondim * float(torus.system.l_km)

    converged_corrector = bool(invariance_norm <= invariance_floor)
    converged_independent = bool(indep_nondim <= independent_tol_nondim)
    passes = bool(converged_corrector and converged_independent)

    return V1VerdictQP(
        candidate_id=candidate_id,
        invariance_residual_fourier_norm=float(invariance_norm),
        independent_invariance_residual_nondim=float(indep_nondim),
        independent_residual_km=float(indep_km),
        v1_floor_invariance=float(invariance_floor),
        independent_floor_nondim=float(independent_tol_nondim),
        converged_corrector=converged_corrector,
        converged_independent=converged_independent,
        passes_v1_qp=passes,
        n_off_grid_samples=int(n_off_grid_samples),
        n_modes=int(torus.n_modes),
        notes=notes,
    )


__all__ = [
    "V1_QP_FOURIER_FLOOR",
    "V1_QP_INDEPENDENT_FLOOR",
    "V1_QP_INDEPENDENT_N_SAMPLES",
    "V1VerdictQP",
    "run_v1_qp",
]
