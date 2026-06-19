"""Tests for the BVP integral-constraint corrector core (#380 Steps 1-3).

Sourced-golden discipline (per ``feedback_golden_tests_sourced_only``)
----------------------------------------------------------------------
EXPECTED sides trace to a PUBLISHED source or are ANALYTICAL:

  * The CR3BP periodic-orbit IC + period come from the catalogue row
    ``braik-ross-c11a-cycler-2026`` (the same sourced golden used in
    ``tests/search/test_cr3bp_3d_corrector.py``).
  * The time-integral constraint's target is the SOURCED period T (not a value
    our code computed): q = integral 1 dt = T, EXPECTED = T.
  * The augmented-Jacobian-vs-finite-difference test is a self-consistency
    check (no external source needed) — same pattern as existing gradient tests.

NB: the Jacobi-drift constraint (and its CR3BP/BCR4BP tests) was removed — a
conserved quantity (dC/dt ~ 0) integrated as an augmented quadrature at tight
tol is numerically pathological (noisy ~0 => solve_ivp step collapse). See
docs/notes/2026-06-19-380-bvp-integral-corrector-blueprint.md.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bvp_integral import (
    IDX_X,
    IDX_YDOT,
    AugmentedArc,
    IntegralConstraint,
    correct_with_integral_constraints,
    propagate_augmented_cr3bp,
    time_integral_constraint,
)

# ---------------------------------------------------------------------------
# Sourced golden: Braik-Ross C11a planar Earth-Moon (1,1) cycler.
# data/catalogue.yaml row braik-ross-c11a-cycler-2026 (also used in
# tests/search/test_cr3bp_3d_corrector.py).
# ---------------------------------------------------------------------------
C11A_X0 = -0.8116406668238195
C11A_YDOT0 = -0.11859055759763637
C11A_PERIOD_TU = 9.69107744379376
EM_MU = 1.2150584270572e-2  # Braik-Ross 2026 Table 1
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _cr3bp_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=EM_MU, primary="earth", secondary="moon", l_km=EM_L_KM, t_s=EM_T_S)


def _c11a_state() -> np.ndarray:
    return np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Augmented Jacobian (dq_dx0 + an STM column) matches central FD.
# Self-consistency check (no external source). Uses a well-behaved
# state-dependent constraint h = x (analytic gradient [1,0,0,0,0,0]) so dq_dx0
# is non-trivial. NB: the original Jacobi-drift constraint was removed — a
# conserved quantity (dC/dt ~ 0) is numerically pathological to integrate as an
# augmented quadrature at tight tol (noisy ~0 => solve_ivp step collapse / hang).
# The actionable machinery is fully validated here without it.
# ---------------------------------------------------------------------------


def test_augmented_jacobian_matches_fd() -> None:
    """``dq_dx0`` and the STM agree with central finite differences over X0."""
    system = _cr3bp_system()
    state0 = _c11a_state()
    t = 2.0  # short arc; FD on a full sensitive period is noisier than needed.

    # State-dependent constraint h = x with analytic gradient -> non-trivial dq_dx0.
    def _hx(t_: float, x6: np.ndarray, params: object) -> float:
        return float(x6[0])

    def _hx_grad(t_: float, x6: np.ndarray, params: object) -> np.ndarray:
        g = np.zeros(6, dtype=np.float64)
        g[0] = 1.0
        return g

    constraint = IntegralConstraint(
        label="x_integral", integrand=_hx, integrand_grad=_hx_grad, target=0.0
    )
    base = propagate_augmented_cr3bp(system, state0, t, [constraint], rtol=1e-12, atol=1e-12)

    eps = 1e-6
    fd_dq = np.zeros(6, dtype=np.float64)
    fd_stm_col0 = np.zeros(6, dtype=np.float64)
    for k in range(6):
        sp = state0.copy()
        sm = state0.copy()
        h = eps * max(1.0, abs(float(state0[k])))
        sp[k] += h
        sm[k] -= h
        arc_p = propagate_augmented_cr3bp(system, sp, t, [constraint], rtol=1e-12, atol=1e-12)
        arc_m = propagate_augmented_cr3bp(system, sm, t, [constraint], rtol=1e-12, atol=1e-12)
        fd_dq[k] = (float(arc_p.q_values[0]) - float(arc_m.q_values[0])) / (2.0 * h)
        # STM column 0: d(state_f)/d(X0_k) for the first state component.
        fd_stm_col0[k] = (arc_p.state_f[0] - arc_m.state_f[0]) / (2.0 * h)

    analytic_dq = base.dq_dx0[0]
    denom = np.maximum(np.abs(fd_dq), 1e-8)
    rel = np.abs(analytic_dq - fd_dq) / denom
    assert np.max(rel) < 1e-3, (
        f"dq_dx0 vs FD mismatch: analytic={analytic_dq}, fd={fd_dq}, rel={rel}"
    )

    # STM row 0 (d state_f[0] / d X0[k]) vs FD.
    stm_row0 = base.stm[0, :]
    denom_s = np.maximum(np.abs(fd_stm_col0), 1e-6)
    rel_s = np.abs(stm_row0 - fd_stm_col0) / denom_s
    assert np.max(rel_s) < 1e-6, (
        f"STM row0 vs FD mismatch: analytic={stm_row0}, fd={fd_stm_col0}, rel={rel_s}"
    )


# ---------------------------------------------------------------------------
# Test 3: corrector converges with a time-integral (period) constraint.
# Period target = SOURCED catalogue period (not our-computed).
# ---------------------------------------------------------------------------


def test_corrector_time_integral() -> None:
    """The corrector closes the point residual AND pins the period to the
    sourced catalogue period via a time-integral constraint (q = integral 1 dt = T)."""
    system = _cr3bp_system()
    state0 = _c11a_state()
    constraint = time_integral_constraint(C11A_PERIOD_TU)

    constraints = [constraint]

    def propagate_fn(s: np.ndarray, t: float) -> AugmentedArc:
        return propagate_augmented_cr3bp(system, s, t, constraints, rtol=1e-12, atol=1e-12)

    def eom(t: float, x: np.ndarray) -> np.ndarray:
        return cr3bp.cr3bp_eom(t, x, system.mu)

    result = correct_with_integral_constraints(
        propagate_fn,
        state0,
        C11A_PERIOD_TU,
        # Symmetric perpendicular-crossing closure at T/2 + period free, plus
        # the time-integral pinning T to the sourced value.
        free_vars=(IDX_X, IDX_YDOT, 6),
        point_residual_indices=(1, 3, 5),  # y, xdot, zdot
        is_half_period_residual=True,
        integral_constraints=constraints,
        independent_eom=eom,
        tol=1e-9,
        max_iter=60,
    )
    assert result.converged, (
        f"corrector must close with the time-integral constraint; "
        f"point={result.point_residual:.3e}, integral={result.integral_residual:.3e}, "
        f"independent={result.independent_closure_residual:.3e}, n_iter={result.n_iter}"
    )
    # The integral residual = |T_final - T_target| ~ 0.
    assert result.integral_residual < 1e-7
    # Period pinned to the sourced golden.
    assert abs(result.period - C11A_PERIOD_TU) < 1e-6
    # q = integral 1 dt = T.
    assert abs(float(result.q_values[0]) - result.period) < 1e-7
    # Point closure preserved (independent Radau).
    assert result.independent_closure_residual < 1e-6


# ---------------------------------------------------------------------------
# Test 5: empty-constraints fast path matches the plain STM propagation.
# ---------------------------------------------------------------------------


def test_empty_constraints_matches_plain_stm() -> None:
    """With no constraints the augmented propagator must equal the existing
    plain STM path exactly (identical result / cost)."""
    system = _cr3bp_system()
    state0 = _c11a_state()
    t = 3.0
    arc_aug = propagate_augmented_cr3bp(system, state0, t, [], rtol=1e-12, atol=1e-12)
    arc_plain = cr3bp.propagate(system, state0, t, with_stm=True, rtol=1e-12, atol=1e-12)
    assert arc_plain.stm is not None
    assert np.allclose(arc_aug.state_f, arc_plain.state_f, rtol=0, atol=0)
    assert np.allclose(arc_aug.stm, arc_plain.stm, rtol=0, atol=0)
    assert arc_aug.q_values.shape == (0,)
    assert arc_aug.dq_dx0.shape == (0, 6)


# ---------------------------------------------------------------------------
# Test 6: deferred Sun-commensurate constraint is a clear stub.
# ---------------------------------------------------------------------------


def test_sun_commensurate_constraint_is_deferred_stub() -> None:
    from cyclerfinder.search.bvp_integral import sun_commensurate_period_constraint

    with pytest.raises(NotImplementedError, match="Step 6"):
        sun_commensurate_period_constraint()
