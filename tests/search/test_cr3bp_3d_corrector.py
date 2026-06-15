"""Tests for the full-3D broken-plane CR3BP corrector (#291 Phase 1).

Sourced golden (per ``feedback_golden_tests_sourced_only``)
-----------------------------------------------------------
The PLANAR seed comes from ``data/catalogue.yaml`` row
``braik-ross-c11a-cycler-2026``:

    state_nd: [-0.8116406668238195, 0.0, 0.0, 0.0, -0.11859055759763637, 0.0]
    period_nd: 9.69107744379376
    mass_ratio: 1.2150584270572e-2  # Braik-Ross 2026 Table 1
    jacobi_constant: 3.1294

This planar state IS the golden -- it traces to a published Braik-Ross 2026
member, reproduced in-repo at +0.0011% (#249). For the 3D extension we CANNOT
write a "expected 3D state == X" assertion against a number our own code
produced (that would be circular). Instead the 3D tests assert on
**topology + closure**:

  * ``|z0| > eps`` confirms the orbit is genuinely 3D (not a planar collapse).
  * Independent re-propagation under Radau closes the IC to <1e-9 nondim
    (~384 m in Earth-Moon units).
  * Period is finite and within an order of magnitude of the planar seed
    (the spike showed the 3D member at T = 44.37 d vs planar 42.14 d).
  * Jacobi constant is consistent (computed at the corrected IC).

Tests
-----
  1. Planar reproduction (z0_guess = 0): recovers the catalogue Braik-Ross
     C11a state, ``degenerate_planar=True``, closure < 1e-6.
  2. Non-trivial 3D closure: ``z0_guess = 0.05`` finds a 3D orbit with
     ``|z0| > 0.05``, independent closure < 1e-8.
  3. Small-z0 collapse: ``z0_guess = 1e-3`` correctly collapses to the planar
     manifold (the corrector returns ``degenerate_planar=True``); this is the
     mathematically correct behaviour, not a bug.
  4. STM coupling sanity: verify that for a 3D IC the STM's z-coupling
     entries (``uzz``, ``uxz``, ``uyz``) are exercised in the variational
     integration -- i.e. perturbing z0 propagates to non-zero
     ``state_f - state0`` even when ``y0 = ydot0 = 0`` initially.
  5. Bogus seed: a totally implausible seed (far from any periodic orbit)
     does NOT converge silently -- the independent closure check rejects it.
  6. Full-asymmetric mode: the default 7-unknown asymmetric API converges on
     the planar seed and returns a state matching the planar member.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    FREE_VARS_SYMMETRIC_TULIP,
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_Z,
    IDX_ZDOT,
    RESIDUAL_FULL_STATE_AT_T,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)

# ---------------------------------------------------------------------------
# Sourced golden: Braik-Ross C11a planar Earth-Moon (1,1) cycler.
# data/catalogue.yaml row braik-ross-c11a-cycler-2026 (planar state IS the
# sourced golden; the 3D extension's correctness is verified by topology +
# closure, NOT by a hand-computed 3D number).
# ---------------------------------------------------------------------------
C11A_X0 = -0.8116406668238195
C11A_YDOT0 = -0.11859055759763637
C11A_PERIOD_TU = 9.69107744379376
C11A_JACOBI = 3.1294
EM_MU = 1.2150584270572e-2  # Braik-Ross 2026 Table 1
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _make_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP at the catalogued (sourced) Braik-Ross mu."""
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


# ---------------------------------------------------------------------------
# Test 1: planar reproduction (z0_guess = 0).
# ---------------------------------------------------------------------------


def test_planar_seed_recovers_braik_ross_c11a() -> None:
    """Seed with the planar Braik-Ross C11a IC; symmetric corrector mode
    should recover the same orbit with ``degenerate_planar=True`` and
    a tight independent closure.

    Sourced golden: the planar state itself is the catalogue row, so this is
    an honest reproduce-before-trust gate. The corrector recovers the same
    IC (within numerical noise) and confirms the planar manifold is invariant.
    """
    system = _make_system()
    state0_planar = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    result = correct_general_periodic_3d(
        system,
        state0_planar,
        C11A_PERIOD_TU,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    assert isinstance(result, Periodic3DOrbit)
    assert result.converged, (
        f"planar seed must reproduce; corrector={result.corrector_residual:.3e}, "
        f"independent={result.independent_closure_residual:.3e}"
    )
    # State recovered to numerical noise.
    assert abs(float(result.state0[IDX_X]) - C11A_X0) < 1e-9
    assert abs(float(result.state0[IDX_YDOT]) - C11A_YDOT0) < 1e-9
    assert abs(result.T_TU - C11A_PERIOD_TU) < 1e-7
    # Planar manifold invariant: z0 stays ~0.
    assert result.degenerate_planar
    assert abs(float(result.state0[IDX_Z])) < 1e-9
    assert abs(float(result.state0[IDX_ZDOT])) < 1e-9
    # Jacobi held.
    assert abs(result.jacobi - C11A_JACOBI) < 1e-9
    # Independent closure tight (planar member is moderately unstable).
    assert result.independent_closure_residual < 1e-6


# ---------------------------------------------------------------------------
# Test 2: non-trivial 3D closure with z0_guess = 0.05 (the spike's payload).
# ---------------------------------------------------------------------------


def test_z0_guess_p05_finds_nontrivial_3d_orbit() -> None:
    """z0_guess = 0.05 from the planar C11a seed lands a 3D orbit.

    Topology assertion: the corrected orbit has |z0| >> 0.05 (the spike found
    z0 ~ -0.241). Closure assertion: independent Radau re-propagation is
    tight. This does NOT assert "z0 == -0.241" because that number was
    computed by our own code; instead it asserts the orbit is genuinely 3D
    (planar manifold escape) and physically closed.
    """
    system = _make_system()
    state0_guess = np.array([C11A_X0, 0.0, 0.05, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    result = correct_general_periodic_3d(
        system,
        state0_guess,
        C11A_PERIOD_TU,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
        max_iter=80,
    )
    assert result.converged, (
        f"3D corrector must close; corrector={result.corrector_residual:.3e}, "
        f"independent={result.independent_closure_residual:.3e}, n_iter={result.n_iter}"
    )
    # Genuinely 3D, not a planar collapse.
    assert not result.degenerate_planar
    assert abs(float(result.state0[IDX_Z])) > 0.05, (
        f"3D orbit must escape the planar manifold; z0={result.state0[IDX_Z]:.3e}"
    )
    # Tight independent closure.
    assert result.independent_closure_residual < 1e-8
    # Period plausibility (the spike found T = 44.37 d; planar was 42.14 d --
    # same order of magnitude, ~10% offset is family-internal).
    t_days = result.T_TU * (EM_T_S / 86400.0)
    assert 30.0 < t_days < 60.0, f"3D period outside plausible range: {t_days} d"
    # Jacobi diagnostic: the 3D member sits at a DIFFERENT Jacobi than the
    # planar seed (the corrector doesn't pin C in this mode). The spike found
    # C ~ 3.027 at the 3D member; we just assert finite + reasonable.
    assert 2.5 < result.jacobi < 3.5


# ---------------------------------------------------------------------------
# Test 3: small-z0 collapse to planar manifold.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("z0_guess", [1e-4, 1e-3, 5e-3])
def test_small_z0_collapses_to_planar(z0_guess: float) -> None:
    """Small z0_guess collapses to the planar member (manifold invariance).

    This is mathematically correct -- the planar manifold (z = zdot = 0) is
    dynamically invariant under the CR3BP, so any IC inside the planar
    manifold's basin of attraction for the corrector lands back on it. The
    test asserts this happens cleanly and the ``degenerate_planar`` flag is
    set so callers searching for "genuinely 3D" orbits can screen.

    Honest limitation: the spike (data/spike_287.jsonl) shows this collapse at
    z0_guess up to ~0.01 in the Braik-Ross C11a basin. Larger guesses
    (z0_guess >= 0.05) escape to genuinely 3D members.
    """
    system = _make_system()
    state0_guess = np.array([C11A_X0, 0.0, z0_guess, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    result = correct_general_periodic_3d(
        system,
        state0_guess,
        C11A_PERIOD_TU,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    assert result.converged
    # Collapsed to the planar manifold.
    assert result.degenerate_planar, (
        f"small z0_guess={z0_guess:.1e} should collapse to planar; "
        f"got z0={result.state0[IDX_Z]:.3e}"
    )


# ---------------------------------------------------------------------------
# Test 4: STM z-coupling sanity check.
# ---------------------------------------------------------------------------


def test_stm_z_coupling_is_exercised() -> None:
    """Sanity check: for a 3D IC the STM properly propagates z perturbations.

    Confirms that ``uzz``, ``uxz``, ``uyz`` in :func:`cr3bp.cr3bp_stm_eom`
    are integrated -- a perturbation in z0 must induce a non-zero change in
    BOTH the position AND velocity components at a finite time.

    This is NOT a corrector test, but a sanity check on the underlying
    machinery (cr3bp.cr3bp_stm_eom). If this fails, the corrector's Jacobian
    is structurally wrong.
    """
    system = _make_system()
    # Use the 3D spike-converged member as a stand-in for a 3D IC; we don't
    # need it to be exactly periodic, just genuinely 3D.
    state0 = np.array([C11A_X0, 0.0, -0.241, 0.0, -0.106, 0.0], dtype=np.float64)
    arc = cr3bp.propagate(system, state0, 1.0, with_stm=True, rtol=1e-12, atol=1e-12)
    assert arc.stm is not None
    stm = arc.stm
    # z perturbation must propagate to x (in-plane coupling via uxz):
    # stm[IDX_X, IDX_Z] = dx(t)/dz(0). For a planar IC this would be zero by
    # symmetry; for a 3D IC it must be non-zero.
    assert abs(float(stm[IDX_X, IDX_Z])) > 1e-3, (
        f"STM z->x coupling not exercised; stm[x,z]={stm[IDX_X, IDX_Z]:.3e}"
    )
    # z perturbation propagates to ydot (uyz coupling):
    assert abs(float(stm[IDX_YDOT, IDX_Z])) > 1e-3
    # Symmetric: x perturbation propagates to z (uxz coupling):
    assert abs(float(stm[IDX_Z, IDX_X])) > 1e-3


# ---------------------------------------------------------------------------
# Test 5: bogus seed does not converge silently.
# ---------------------------------------------------------------------------


def test_bogus_seed_does_not_converge_silently() -> None:
    """A wildly off-family seed must NOT report converged.

    The independent closure check is the gate: even if the corrector's own
    residual happens to land at a near-zero L2 by luck (it shouldn't), the
    re-propagation under Radau will catch the discrepancy.
    """
    system = _make_system()
    # Seed far from any periodic orbit: random-ish position, no velocity.
    # This should fail either (a) the corrector residual or (b) the
    # independent closure check.
    state0_bogus = np.array([0.5, 0.3, 0.4, 0.7, 0.2, 0.1], dtype=np.float64)
    result = correct_general_periodic_3d(
        system,
        state0_bogus,
        period_guess=3.0,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        residual_indices=RESIDUAL_FULL_STATE_AT_T,
        is_half_period_residual=False,
        max_iter=30,
        tol=1e-11,
        independent_tol=1e-6,
    )
    # The bogus seed must NOT report converged. Either the corrector didn't
    # land at all, or it landed somewhere whose Radau re-propagation doesn't
    # close.
    assert not result.converged, (
        f"bogus seed must not silently converge; "
        f"corrector={result.corrector_residual:.3e}, "
        f"independent={result.independent_closure_residual:.3e}"
    )


# ---------------------------------------------------------------------------
# Test 6: full-asymmetric mode reproduces the planar Braik-Ross member.
# ---------------------------------------------------------------------------


def test_full_asymmetric_mode_reproduces_planar_member() -> None:
    """The 7-unknown asymmetric API converges on the planar Braik-Ross IC.

    Mode: free vars = all 6 state components + T, residual = full 6D closure
    at T. The planar IC ``(x0, 0, 0, 0, ydot0, 0)`` is a periodic orbit of
    the full CR3BP so the corrector must hold the IC + period to numerical
    noise.
    """
    system = _make_system()
    state0_planar = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    result = correct_general_periodic_3d(
        system,
        state0_planar,
        C11A_PERIOD_TU,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        residual_indices=RESIDUAL_FULL_STATE_AT_T,
        is_half_period_residual=False,
    )
    assert result.converged
    # Recovered the planar member to numerical noise. The min-norm step on
    # the 6x7 system may permit small drift along the local family tangent,
    # so we use a generous tolerance.
    assert np.linalg.norm(result.state0 - state0_planar) < 1e-6
    assert abs(result.T_TU - C11A_PERIOD_TU) < 1e-6
    assert result.degenerate_planar
    assert abs(result.jacobi - C11A_JACOBI) < 1e-9
    assert result.independent_closure_residual < 1e-6


# ---------------------------------------------------------------------------
# Test 7: API sanity (input validation).
# ---------------------------------------------------------------------------


def test_input_validation() -> None:
    """The corrector rejects malformed inputs cleanly."""
    system = _make_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    # Wrong state shape.
    with pytest.raises(ValueError, match="shape"):
        correct_general_periodic_3d(system, np.zeros(5), C11A_PERIOD_TU)
    # Empty free_vars.
    with pytest.raises(ValueError, match="free_vars"):
        correct_general_periodic_3d(system, state0, C11A_PERIOD_TU, free_vars=())
    # Empty residual_indices.
    with pytest.raises(ValueError, match="residual_indices"):
        correct_general_periodic_3d(system, state0, C11A_PERIOD_TU, residual_indices=())
    # Out-of-range index.
    with pytest.raises(ValueError, match="free_vars must be in"):
        correct_general_periodic_3d(system, state0, C11A_PERIOD_TU, free_vars=(0, 7))
    with pytest.raises(ValueError, match="residual_indices must be in"):
        correct_general_periodic_3d(system, state0, C11A_PERIOD_TU, residual_indices=(0, 6))
    # Non-positive period.
    with pytest.raises(ValueError, match="period_guess"):
        correct_general_periodic_3d(system, state0, -1.0)


# ---------------------------------------------------------------------------
# Constants re-export sanity (index aliases unchanged).
# ---------------------------------------------------------------------------


def test_index_aliases() -> None:
    """The IDX_* constants match standard 6D state ordering."""
    assert (IDX_X, IDX_Y, IDX_Z, IDX_XDOT, IDX_YDOT, IDX_ZDOT, IDX_T) == (0, 1, 2, 3, 4, 5, 6)
