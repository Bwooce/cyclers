"""V1 same-model BCR4BP gauntlet tests (#305 Part A).

Validation against a SOURCED BCR4BP family member: the Andreu / Rosales-Jorba
2023 Table 4 POL1 (L1 dynamical substitute), closed at the Sun-commensurate
period ``T = 2*pi / omega_sun`` (n=1) with the canonical-momentum -> velocity
conversion ``vy = py - x`` (the known #292 gotcha). The published POL1 IC is
the SEED (QBCP), not the EXPECTED side of any equality — per
``feedback_golden_tests_sourced_only`` the test asserts CLOSURE under an
independent integrator, not specific IC numbers.

Gates:
  1. POL1 (commensurate, converged) PASSES V1.
  2. A deliberately perturbed IC (off the periodic basin) FAILS the Radau
     same-model floor.

These run in the DEFAULT (non-slow) suite — the POL1 closure + the Radau
re-propagation complete in a few seconds — so the V1 evidence is verified by
CI, never @pytest.mark.slow-skipped.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v1_bcr4bp import (
    V1_BCR4BP_FLOOR_KMS,
    V1_BCR4BP_FLOOR_NONDIM,
    V1VerdictBCR4BP,
    run_v1_bcr4bp,
)
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    BCR4BPPeriodicOrbit,
    correct_bcr4bp_periodic,
)

# Rosales-Jorba (2023) Table 4 POL1 IC (canonical-momentum (x, py)); convert to
# velocity via vy = py - x (CR3BP-limit Hamiltonian convention). SEED only.
_POL1_X = -0.8369141677649317
_POL1_PY = -0.8391311559808445
_POL1_VY = _POL1_PY - _POL1_X
_POL1_SEED = np.array([_POL1_X, 0.0, 0.0, 0.0, _POL1_VY, 0.0], dtype=np.float64)


def _close_pol1() -> BCR4BPPeriodicOrbit:
    """Close POL1 to a nearby BCR4BP periodic orbit (T fixed, Sun-commensurate)."""
    sys_bcr = bcr4bp.andreu_default()
    period_fixed = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)
    return correct_bcr4bp_periodic(
        sys_bcr,
        _POL1_SEED,
        period_fixed,
        sun_commensurate_n=1,
        free_vars=(0, IDX_YDOT),  # IDX_X, IDX_YDOT; T FIXED at commensurate value
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-10,
        independent_tol=1e-6,
        state_step_cap=0.2,
        require_monotone_decrease=False,
        max_iter=80,
        notes="POL1 -> BCR4BP (T fixed at Sun-commensurate n=1) for V1 evidence",
    )


def test_v1_bcr4bp_pol1_passes() -> None:
    """The sourced Andreu POL1 substitute PASSES V1 under a fresh Radau check."""
    orbit = _close_pol1()
    assert orbit.converged, (
        f"POL1 failed to converge (corrector={orbit.corrector_residual:.3e}, "
        f"independent={orbit.independent_closure_residual:.3e}); cannot test V1"
    )
    verdict = run_v1_bcr4bp("andreu-pol1-bcr4bp", orbit)

    assert isinstance(verdict, V1VerdictBCR4BP)
    assert verdict.converged_corrector, (
        f"corrector residual {verdict.corrector_residual_nondim:.3e} above floor"
    )
    assert verdict.converged_independent, (
        f"fresh Radau closure {verdict.independent_closure_nondim:.3e} > "
        f"{verdict.v1_floor_nondim:.3e} nondim"
    )
    assert verdict.passes_v1_bcr4bp, (
        f"POL1 should pass V1: independent_kms={verdict.independent_closure_kms:.3e} "
        f"vs floor {verdict.v1_floor_kms:.3e}"
    )
    # The fresh Radau km/s closure must clear the spec §14 1 m/s bar.
    assert verdict.independent_closure_kms <= V1_BCR4BP_FLOOR_KMS
    # n=1 commensurate: V1 records the drift but does not gate on it.
    assert verdict.sun_commensurate_n == 1
    assert verdict.sun_phase_drift < 1e-10, (
        f"T held fixed at commensurate value but drift={verdict.sun_phase_drift:.3e}"
    )


def test_v1_bcr4bp_perturbed_ic_fails() -> None:
    """A deliberately perturbed (non-periodic) IC FAILS the same-model Radau floor.

    Take the converged POL1 orbit and shove the IC off the periodic basin by a
    finite displacement. The fresh Radau re-propagation no longer closes, so V1
    must FAIL — guarding against a verdict that passes any nearby state.
    """
    orbit = _close_pol1()
    assert orbit.converged

    perturbed_ic = np.asarray(orbit.state_initial, dtype=np.float64).copy()
    perturbed_ic[0] += 1.0e-3  # off-family x displacement (~384 km)
    bad_orbit = BCR4BPPeriodicOrbit(
        state_initial=perturbed_ic,
        period_nondim=orbit.period_nondim,
        sun_commensurate_n=orbit.sun_commensurate_n,
        sun_phase_drift=orbit.sun_phase_drift,
        converged=True,  # pretend the corrector flagged it converged
        corrector_residual=0.0,  # and that its residual was clean
        independent_closure_residual=orbit.independent_closure_residual,
        n_iter=orbit.n_iter,
        system=orbit.system,
        free_vars=orbit.free_vars,
        residual_indices=orbit.residual_indices,
        is_half_period_residual=orbit.is_half_period_residual,
        notes="perturbed off-family IC",
    )
    verdict = run_v1_bcr4bp("perturbed-off-family", bad_orbit)

    assert not verdict.passes_v1_bcr4bp, (
        f"A perturbed off-family IC must FAIL V1; got independent_nondim="
        f"{verdict.independent_closure_nondim:.3e}"
    )
    assert not verdict.converged_independent, (
        "fresh Radau re-propagation should NOT close on a perturbed IC"
    )
    assert verdict.independent_closure_nondim > V1_BCR4BP_FLOOR_NONDIM


def test_v1_bcr4bp_floor_constants_match_spec() -> None:
    """The V1-BCR4BP floors reuse the sourced spec / periodic-V1 values."""
    assert V1_BCR4BP_FLOOR_KMS == 1.0e-3  # spec §14 V1 (1 m/s)
    assert V1_BCR4BP_FLOOR_NONDIM == 1.0e-6  # periodic V1 nondim band
