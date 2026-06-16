"""BCR4BP periodic-orbit corrector tests (#292 Phase 1 Part B).

Five gates:

  1. **CR3BP-limit closure**: at mu_sun = 0 the BCR4BP corrector reproduces
     a sourced CR3BP halo IC (Howell 1984 / NASA TN D-1949 style; used in
     ``tests/search/test_cr3bp_periodic.py``). The structural correctness
     anchor: if mu_sun = 0 closure fails, the BCR4BP -> CR3BP reduction is
     broken regardless of what the EOM tests say at the algebraic level.

  2. **Andreu POL1 seed closes to a NEARBY BCR4BP periodic orbit**: starting
     from the Rosales-Jorba 2023 Table 4 POL1 IC (interpreted as canonical
     momentum ``py``; converted to velocity by ``vy = py - x`` in the
     CR3BP-limit convention), at the Sun-commensurate period
     ``T = 2*pi / omega_sun`` (n=1), with ``T`` FIXED and ``(x, vy)`` free,
     the corrector converges to a BCR4BP L1 dynamical substitute. The
     published POL1 IC is for the QBCP -- the implemented model is the
     standard BCR4BP -- so the converged IC is NEAR but not at POL1. Per
     ``feedback_golden_tests_sourced_only``: the published value is the
     SEED, not the EXPECTED side of an equality assertion. The test
     asserts CLOSURE under both the corrector and an independent Radau
     re-propagation.

  3. **Weak-perturbation halo (mu_sun_eps): at mu_sun = 0.01 * mu_sun_andreu
     the BCR4BP corrector with a CR3BP halo seed still closes. Demonstrates
     the corrector works WITH the Sun term enabled, complementing the
     mu_sun=0 anchor.

  4. **Bogus seed REJECTED**: a deliberately bad IC (state near a primary,
     way off any periodic basin) does NOT silently "converge" -- the
     corrector flags it with ``converged = False``.

  5. **STM Jacobian sanity (corrector-internal)**: the analytic Jacobian
     built from the BCR4BP STM matches the central-FD Jacobian on the
     same residual function at 1e-5 relative tolerance.

Per ``feedback_golden_tests_sourced_only``: no value computed by our own code
is on the EXPECTED side of any equality assertion. All convergence checks are
on RESIDUAL norms or STRUCTURAL closure (state at T = state at 0), where the
EXPECTED side is zero (the mathematical definition of "periodic").

Per ``feedback_orbit_closure_discipline``: every converged orbit also passes
an independent (Radau) full-period re-propagation closure check before being
flagged ``converged = True`` by the corrector itself.
"""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    correct_bcr4bp_periodic,
)

# ---------------------------------------------------------------------------
# Sourced ICs from the Andreu / Rosales-Jorba digest.
# ---------------------------------------------------------------------------

# Rosales-Jorba (2023) Table 4 -- POL1 (L1 dynamical substitute). The published
# IC is in canonical-momentum coordinates ``(x, py)``; the corresponding
# velocity in the CR3BP-limit Hamiltonian convention ``vy = py - x`` (with
# alpha_1 = alpha_3 = 1, alpha_2 = 0). The Phase 1 BCR4BP module does NOT
# carry the QBCP alpha coefficients -- the conversion is therefore exact only
# in the unperturbed limit, with corrections of order alpha_2 ~ O(eps) at the
# implemented model's parameter values. Used as a SEED for the corrector, NOT
# as a closure-equality golden.
_POL1_X = -0.8369141677649317
_POL1_PY = -0.8391311559808445
_POL1_VY = _POL1_PY - _POL1_X  # ~ -0.0022, the velocity in the CR3BP-limit convention

_POL1_SEED = np.array([_POL1_X, 0.0, 0.0, 0.0, _POL1_VY, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Gate 1: CR3BP-limit closure -- structural correctness anchor.
# ---------------------------------------------------------------------------


def test_cr3bp_limit_closure_recovers_cr3bp_halo() -> None:
    """At mu_sun = 0 the BCR4BP corrector recovers a known CR3BP halo orbit.

    Seed: Earth-Moon L1 southern halo IC ``(0.824024728136525, 0,
    -0.054501847320725, 0, 0.164671964079122, 0)`` with T = 2.7549 (Howell
    1984 / NASA TN D-1949 family; used as the sourced seed in
    ``tests/search/test_cr3bp_periodic.py::test_cr3bp_periodic_halo_l1_southern``).

    At mu_sun = 0 the BCR4BP propagator IS the CR3BP propagator (verified
    structurally in ``tests/core/test_bcr4bp.py``). This test verifies the
    Newton corrector + STM-Jacobian assembly also functions in that limit.
    The free vars ``(x, ydot, T)`` and half-period residual ``(y, xdot, zdot)``
    match the classical symmetric-halo perpendicular-crossing setup.
    """
    sys_bcr = bcr4bp.BCR4BPSystem(
        mu=0.012150581600000,
        mu_sun=0.0,  # CR3BP limit
        a_sun_nondim=388.0,
        omega_sun_nondim=0.9252,  # bookkeeping only
    )
    state_seed = np.array(
        [0.824024728136525, 0.0, -0.054501847320725, 0.0, 0.164671964079122, 0.0],
        dtype=np.float64,
    )
    period_guess = 2.7549
    result = correct_bcr4bp_periodic(
        sys_bcr,
        state_seed,
        period_guess,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        independent_tol=1e-6,
        tol=1e-10,
    )
    assert result.converged, (
        f"BCR4BP@mu_sun=0 corrector failed on a sourced CR3BP halo seed. "
        f"corrector_residual={result.corrector_residual:.3e}, "
        f"independent_closure={result.independent_closure_residual:.3e}, "
        f"n_iter={result.n_iter}"
    )
    # Cross-check: at mu_sun=0 the converged orbit must also close under the
    # CR3BP propagator (the two integrators are identical in this limit).
    sys_cr3 = cr3bp.CR3BPSystem(
        mu=0.012150581600000,
        primary="earth",
        secondary="moon",
        l_km=384400.0,
        t_s=375190.0,
    )
    arc = cr3bp.propagate(sys_cr3, result.state_initial, result.period_nondim)
    cr3bp_closure = float(np.linalg.norm(arc.state_f - result.state_initial))
    assert cr3bp_closure < 1e-6, (
        f"BCR4BP@mu_sun=0 converged orbit fails to close under CR3BP: residual={cr3bp_closure:.3e}"
    )
    # The corrected x0 should not have moved far from the seed (Howell halo
    # is locally 1-parameter; the symmetric corrector slides x0 along the
    # family but not by much for a few-iteration Newton).
    assert abs(result.state_initial[IDX_X] - state_seed[IDX_X]) < 0.01


# ---------------------------------------------------------------------------
# Gate 2: Andreu POL1 seed closes to a NEARBY BCR4BP periodic orbit.
# ---------------------------------------------------------------------------


def test_andreu_pol1_seed_closes_in_bcr4bp() -> None:
    """The Rosales-Jorba POL1 IC (QBCP) seeds a closing BCR4BP periodic orbit.

    Setup:
      * Seed: ``_POL1_SEED`` (x = -0.8369..., vy from py via the CR3BP-limit
        ``vy = py - x``; see _POL1_VY).
      * Period: FIXED at ``T = 2*pi / omega_sun`` (n=1 Sun synodic period).
      * Free vars: ``(x, vy)``. The published POL1 has y = z = vx = vz = 0
        at t = 0 (perpendicular crossing on the x-axis); the symmetric
        Newton corrector preserves that structure.
      * Residual: ``(y, vx, vz)`` at T/2 (perpendicular re-crossing closure).

    The published POL1 IC is for the QBCP -- not the standard BCR4BP this
    module implements. The model gap means the converged IC is NEAR (within
    a few percent in nondim units) but NOT AT the published POL1 numbers.
    Per ``feedback_golden_tests_sourced_only``: the published value is the
    SEED; the test asserts CLOSURE (corrector residual + independent Radau
    closure both below thresholds) and TOPOLOGY (the BCR4BP analogue stays
    in the L1 region near x = -0.84), not specific IC values.
    """
    sys_bcr = bcr4bp.andreu_default()
    period_fixed = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)
    result = correct_bcr4bp_periodic(
        sys_bcr,
        _POL1_SEED,
        period_fixed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT),  # T FIXED at the commensurate value
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-10,
        independent_tol=1e-6,
        state_step_cap=0.2,
        require_monotone_decrease=False,
        max_iter=80,
        notes="Andreu POL1 -> BCR4BP closure (T fixed at Sun-commensurate n=1)",
    )
    assert result.converged, (
        f"POL1 seed failed to converge to a BCR4BP periodic orbit. "
        f"corrector_residual={result.corrector_residual:.3e}, "
        f"independent_closure={result.independent_closure_residual:.3e}, "
        f"sun_phase_drift={result.sun_phase_drift:.3e}, n_iter={result.n_iter}, "
        f"state={result.state_initial}, T={result.period_nondim}"
    )
    # Independent (Radau) closure tight gate.
    assert result.independent_closure_residual < 1e-6, (
        f"POL1 -> BCR4BP independent closure={result.independent_closure_residual:.3e} > 1e-6"
    )
    # Strict-periodic BCR4BP requires the period to commensurate with the Sun
    # synodic frequency; we held T fixed at the n=1 value, so the drift must
    # be zero (within floating-point).
    assert result.sun_phase_drift < 1e-12, (
        f"T held fixed at the commensurate value but sun_phase_drift="
        f"{result.sun_phase_drift:.3e} is nonzero -- bookkeeping bug?"
    )
    # Topology: the converged orbit must stay in the L1 region (x near -0.84,
    # NOT have wandered out of the libration neighbourhood).
    x_converged = float(result.state_initial[IDX_X])
    assert -0.90 < x_converged < -0.80, (
        f"Converged x0={x_converged} outside the L1-substitute neighbourhood "
        f"[-0.90, -0.80] -- POL1 seed walked off-family."
    )


# ---------------------------------------------------------------------------
# Gate 3: Weak-perturbation halo: BCR4BP at small mu_sun closes a halo.
# ---------------------------------------------------------------------------


def test_weak_sun_perturbation_halo_closes() -> None:
    """At small but nonzero mu_sun the corrector still closes a CR3BP halo seed.

    Complements Gate 1 (mu_sun = 0 exactly) by exercising the Sun-term code
    path. mu_sun = 1% of the Andreu value is small enough that the BCR4BP
    halo lies close to the CR3BP halo (a regular perturbation problem) but
    nonzero so the EOM, STM, and indirect-term code paths are all active.
    """
    sys_andreu = bcr4bp.andreu_default()
    # 1% of the Andreu Sun mass -- weak enough to keep the CR3BP halo seed
    # in the convergence basin, strong enough to exercise the Sun code path.
    sys_weak = bcr4bp.BCR4BPSystem(
        mu=sys_andreu.mu,
        mu_sun=0.01 * sys_andreu.mu_sun,
        a_sun_nondim=sys_andreu.a_sun_nondim,
        omega_sun_nondim=sys_andreu.omega_sun_nondim,
        theta_sun0=0.0,
    )
    state_seed = np.array(
        [0.824024728136525, 0.0, -0.054501847320725, 0.0, 0.164671964079122, 0.0],
        dtype=np.float64,
    )
    period_guess = 2.7549
    result = correct_bcr4bp_periodic(
        sys_weak,
        state_seed,
        period_guess,
        sun_commensurate_n=1,  # bookkeeping (period isn't strictly commensurate here)
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-9,
        independent_tol=1e-3,  # the period isn't Sun-commensurate; full-T closure is loose
        state_step_cap=0.05,
        max_iter=60,
        notes="Weak-Sun halo (1% mu_sun_andreu)",
    )
    # The corrector residual gate is tight (half-period symmetric closure
    # closes well even at moderate mu_sun); the independent-closure gate is
    # loose because the half-period symmetric closure does NOT imply
    # full-period closure when T is not commensurate with 2*pi/omega_sun.
    assert result.corrector_residual < 1e-7, (
        f"Weak-Sun halo corrector residual={result.corrector_residual:.3e} > 1e-7"
    )
    # Topology: still a halo.
    assert -0.06 < float(result.state_initial[2]) < -0.03


# ---------------------------------------------------------------------------
# Gate 4: Bogus seed REJECTED.
# ---------------------------------------------------------------------------


def test_bogus_seed_rejected_not_silent_success() -> None:
    """A deliberately bad IC fails to converge -- the corrector flags it.

    Seed: a state very close to the Moon (well inside its Hill sphere) with
    a hugely off-family velocity. There's no obvious local periodic orbit at
    the chosen period, and the propagator may fail near the Moon. The
    corrector must NEVER return ``converged = True`` for such a seed.
    """
    sys_bcr = bcr4bp.andreu_default()
    bogus_seed = np.array(
        [1.0 - 0.012150581600000 + 0.001, 0.0, 0.0, 0.0, 5.0, 0.0], dtype=np.float64
    )
    period_seed = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)

    result = correct_bcr4bp_periodic(
        sys_bcr,
        bogus_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-10,
        independent_tol=1e-6,
        max_iter=20,
    )
    assert not result.converged, (
        f"Bogus seed silently 'converged' -- this is a FALSE POSITIVE. "
        f"corrector_residual={result.corrector_residual:.3e}, "
        f"independent_closure={result.independent_closure_residual:.3e}. "
        "The corrector must NEVER return converged=True from a seed with "
        "no nearby periodic basin."
    )


# ---------------------------------------------------------------------------
# Gate 5: STM Jacobian sanity (analytic vs central-difference).
# ---------------------------------------------------------------------------


def test_stm_jacobian_matches_finite_difference() -> None:
    """The corrector's analytic Jacobian matches central-FD at 1e-5.

    Builds the Jacobian as the corrector does it -- ``J[i,j] = STM[i,j] -
    delta_{i,j}`` for state columns; ``J[i,T-col] = f_end[i] * t_scale`` for
    the period column -- then compares to a central-difference Jacobian on
    the same residual function. The model parameters are the Andreu default
    BCR4BP, state away from primaries, T = 1 TU (well below one Sun period).
    """
    sys_bcr = bcr4bp.andreu_default()
    state0 = np.array([0.85, 0.0, 0.0, 0.0, 0.5, 0.0], dtype=np.float64)
    period = 1.0
    free_vars = (IDX_X, IDX_YDOT, IDX_T)
    residual_indices = (IDX_X, IDX_Y, IDX_XDOT, IDX_YDOT)

    arc = bcr4bp.propagate_bcr4bp(sys_bcr, state0, period, with_stm=True)
    assert arc.stm is not None
    sf, stm = arc.state_f, arc.stm

    # Analytic Jacobian (same construction the corrector uses).
    n_res = len(residual_indices)
    n_free = len(free_vars)
    jac = np.zeros((n_res, n_free), dtype=np.float64)
    f_end = bcr4bp.bcr4bp_eom(period, sf, sys_bcr)
    for col, unknown in enumerate(free_vars):
        if unknown == IDX_T:
            for row, ridx in enumerate(residual_indices):
                jac[row, col] = float(f_end[ridx])
        else:
            for row, ridx in enumerate(residual_indices):
                val = float(stm[ridx, unknown])
                if ridx == unknown:
                    val -= 1.0
                jac[row, col] = val

    # Central-FD Jacobian.
    fd_jac = np.zeros((n_res, n_free), dtype=np.float64)
    eps = 1e-5
    for col, unknown in enumerate(free_vars):
        if unknown == IDX_T:
            arc_p = bcr4bp.propagate_bcr4bp(sys_bcr, state0, period + eps)
            arc_m = bcr4bp.propagate_bcr4bp(sys_bcr, state0, period - eps)
            d_p = arc_p.state_f - state0
            d_m = arc_m.state_f - state0
        else:
            sp = state0.copy()
            sp[unknown] += eps
            sm = state0.copy()
            sm[unknown] -= eps
            arc_p = bcr4bp.propagate_bcr4bp(sys_bcr, sp, period)
            arc_m = bcr4bp.propagate_bcr4bp(sys_bcr, sm, period)
            d_p = arc_p.state_f - sp
            d_m = arc_m.state_f - sm
        fd_jac[:, col] = (d_p[list(residual_indices)] - d_m[list(residual_indices)]) / (2.0 * eps)

    scale = max(1.0, float(np.max(np.abs(jac))))
    rel_delta = float(np.max(np.abs(jac - fd_jac))) / scale
    assert rel_delta < 1e-5, (
        f"Jacobian (analytic vs central-FD) rel_delta={rel_delta:.3e} > 1e-5 "
        f"(abs max|delta|={np.max(np.abs(jac - fd_jac)):.3e}, scale={scale:.3e})"
    )


# ---------------------------------------------------------------------------
# Bonus: Sun-phase drift bookkeeping sanity.
# ---------------------------------------------------------------------------


def test_sun_phase_drift_at_commensurate_period() -> None:
    """A period set exactly to the Sun-commensurate value gives zero drift."""
    sys_bcr = bcr4bp.andreu_default()
    period_n3 = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=3)
    expected_drift = abs(sys_bcr.omega_sun_nondim * period_n3 - 2.0 * math.pi * 3)
    assert expected_drift < 1e-13
