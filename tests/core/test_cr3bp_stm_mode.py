"""Tests for the ``stm_mode`` parameter on :func:`cyclerfinder.core.cr3bp.propagate`.

#372 P372.5 — exercises the Pellegrini-Russell 2016 fixed-path STM mitigation
introduced in commit ``d901c90``. Reference:

  Pellegrini, E. & Russell, R.P. (2016), "On the Computation and Accuracy of
  Trajectory State Transition Matrices", *Journal of Guidance, Control, and
  Dynamics*, DOI 10.2514/1.G001920.

Coverage:

  * ``stm_mode='variable'`` is the default and matches the existing
    variable-step variational behaviour exactly (backward-compat).
  * ``stm_mode='fixed_path'`` produces a finite, non-NaN STM on a typical
    Earth-Moon orbit.
  * For a *well-conditioned* orbit (small Arenstorf at lambda_max ~ 500),
    both modes agree at the multiplier level to high precision (~ 1e-9
    relative).
  * For a *highly sensitive* orbit (the #347 Phase 1 (3,2) cluster point at
    lambda_max ~ 3e4) the disagreement metric is finite and bounded — the
    comparator works but we do not require a specific magnitude (the
    cross-check artifact :file:`data/372_stm_mode_crosscheck.jsonl` records
    the actual values).
  * Invalid ``stm_mode`` strings raise ``ValueError``.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp


@pytest.fixture
def arenstorf_system() -> cr3bp.CR3BPSystem:
    """The classical Arenstorf test mu (mu=0.012277471)."""
    return cr3bp.CR3BPSystem(mu=0.012277471, primary="test", secondary="test", l_km=1.0, t_s=1.0)


@pytest.fixture
def arenstorf_state() -> np.ndarray:
    """Well-conditioned IC used by other CR3BP tests."""
    return np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0])


@pytest.fixture
def em_system() -> cr3bp.CR3BPSystem:
    """Physical Earth-Moon CR3BP mu per Ross-Roberts-Tsoukkas 2025 (AAS 25-621, p.3)."""
    return cr3bp.CR3BPSystem(
        mu=1.2150584270572e-2,
        primary="Earth",
        secondary="Moon",
        l_km=1.0,
        t_s=1.0,
    )


@pytest.fixture
def c32_anchor_state() -> np.ndarray:
    """#347 Phase 1 (3,2) Earth-Moon C32 post-saddle-center anchor IC.

    Source: ``data/floquet_phase1_reproduction.jsonl`` first record
    (``parent_state0``). C = 3.14180. lambda_max ~ 3.0e4.
    """
    return np.array(
        [
            -0.28434935291806585,
            0.0,
            0.0,
            0.0,
            -2.0533998110844403,
            0.0,
        ]
    )


C32_ANCHOR_PERIOD = 17.4642064064264  # TU, from the same record.


def test_stm_mode_default_is_variable_and_matches_existing(
    arenstorf_system: cr3bp.CR3BPSystem, arenstorf_state: np.ndarray
) -> None:
    """The default ``stm_mode='variable'`` keyword must be backward-compatible.

    Calling ``propagate(with_stm=True)`` with no ``stm_mode`` keyword and with
    ``stm_mode='variable'`` must produce IDENTICAL STMs (same code path).
    """
    arc_default = cr3bp.propagate(arenstorf_system, arenstorf_state, 5.0, with_stm=True)
    arc_explicit = cr3bp.propagate(
        arenstorf_system, arenstorf_state, 5.0, with_stm=True, stm_mode="variable"
    )
    assert arc_default.stm is not None
    assert arc_explicit.stm is not None
    np.testing.assert_array_equal(arc_default.state_f, arc_explicit.state_f)
    np.testing.assert_array_equal(arc_default.stm, arc_explicit.stm)


def test_stm_mode_fixed_path_runs_clean(
    arenstorf_system: cr3bp.CR3BPSystem, arenstorf_state: np.ndarray
) -> None:
    """Fixed-path mode must produce a finite, real, non-NaN STM of shape (6,6)."""
    arc = cr3bp.propagate(
        arenstorf_system, arenstorf_state, 5.0, with_stm=True, stm_mode="fixed_path"
    )
    assert arc.stm is not None
    assert arc.stm.shape == (6, 6)
    assert np.all(np.isfinite(arc.stm))
    assert arc.stm.dtype == np.float64
    # The state must also propagate to the right time.
    assert arc.t == 5.0
    assert arc.state_f.shape == (6,)
    assert np.all(np.isfinite(arc.state_f))


def test_stm_mode_agreement_on_well_conditioned_orbit(
    arenstorf_system: cr3bp.CR3BPSystem, arenstorf_state: np.ndarray
) -> None:
    """On a well-conditioned orbit, variable and fixed_path agree at high precision.

    The Arenstorf test orbit (lambda_max ~ 537) has none of the conditioning
    issues that motivate fixed_path. Pellegrini-Russell 2016 §III.A.2 (p. 6)
    predicts the two methods should agree to integrator precision in this
    regime. We check the dominant Floquet multiplier (the integration-error-
    sensitive quantity) to ~1e-9 relative agreement.
    """
    arc_var = cr3bp.propagate(
        arenstorf_system, arenstorf_state, 5.0, with_stm=True, stm_mode="variable"
    )
    arc_fix = cr3bp.propagate(
        arenstorf_system, arenstorf_state, 5.0, with_stm=True, stm_mode="fixed_path"
    )
    assert arc_var.stm is not None
    assert arc_fix.stm is not None
    eigs_var = np.sort(np.abs(np.linalg.eigvals(arc_var.stm)))[::-1]
    eigs_fix = np.sort(np.abs(np.linalg.eigvals(arc_fix.stm)))[::-1]
    rel_disagree_lam_max = abs(eigs_fix[0] - eigs_var[0]) / eigs_var[0]
    # Well-conditioned: empirical agreement is ~1.5e-13 (measured 2026-06-17).
    # Floor at 1e-9 for portability across BLAS implementations.
    assert rel_disagree_lam_max < 1e-9, (
        f"well-conditioned-orbit disagreement {rel_disagree_lam_max:.3e} exceeds "
        "1e-9; either the variable-step or fixed-path code path has regressed"
    )


def test_stm_mode_disagreement_on_sensitive_orbit_is_bounded(
    em_system: cr3bp.CR3BPSystem, c32_anchor_state: np.ndarray
) -> None:
    """On the (3,2) cluster point, the comparator works; magnitude bounded.

    Per the #372 P372.3 cross-check (``data/372_stm_mode_crosscheck.jsonl``),
    the (3,2) C32 anchor has lambda_max ~ 3e4 and shows ~0.07% relative
    disagreement at lambda_max + ~0.35% trivial-pair smear under the
    variable-step path (the Pellegrini-Russell 2016 eq. 17 contamination on
    display).

    This test verifies the COMPARATOR works: both modes return finite STMs,
    eigenvalues are computable, and the disagreement is in a bounded, sane
    range. It does NOT pin a specific magnitude (which is integrator-version
    dependent) -- the JSONL artifact carries the recorded values for
    reproducibility.
    """
    arc_var = cr3bp.propagate(
        em_system,
        c32_anchor_state,
        C32_ANCHOR_PERIOD,
        with_stm=True,
        stm_mode="variable",
    )
    arc_fix = cr3bp.propagate(
        em_system,
        c32_anchor_state,
        C32_ANCHOR_PERIOD,
        with_stm=True,
        stm_mode="fixed_path",
    )
    assert arc_var.stm is not None
    assert arc_fix.stm is not None
    assert np.all(np.isfinite(arc_var.stm))
    assert np.all(np.isfinite(arc_fix.stm))
    eigs_var = np.sort(np.abs(np.linalg.eigvals(arc_var.stm)))[::-1]
    eigs_fix = np.sort(np.abs(np.linalg.eigvals(arc_fix.stm)))[::-1]
    # Both lambda_max values must be at the right order of magnitude (~3e4 per the
    # P372.3 cross-check); a regression in either path would yield wildly different
    # multiplier sets.
    assert 1e3 < eigs_var[0] < 1e6, f"variable-step lambda_max out of range: {eigs_var[0]}"
    assert 1e3 < eigs_fix[0] < 1e6, f"fixed-path lambda_max out of range: {eigs_fix[0]}"
    # The disagreement comparator: a finite, real, non-negative number.
    rel_disagree = abs(eigs_fix[0] - eigs_var[0]) / eigs_var[0]
    assert np.isfinite(rel_disagree)
    assert rel_disagree >= 0.0
    # Recorded magnitude per data/372_stm_mode_crosscheck.jsonl is 7.2e-4 (0.07%);
    # accept anything < 1% (the plan's acceptance gate). If this assertion ever
    # fails, the substrate has drifted and #372 P372.3 needs to be re-run.
    assert rel_disagree < 1e-2, (
        f"(3,2) lambda_max disagreement {rel_disagree:.3e} exceeds 1%; "
        "STM contamination has grown beyond the #372 P372.3 acceptance gate"
    )


def test_stm_mode_invalid_raises_value_error(
    arenstorf_system: cr3bp.CR3BPSystem, arenstorf_state: np.ndarray
) -> None:
    """Unsupported ``stm_mode`` strings must surface a ValueError, not a silent fallthrough."""
    with pytest.raises(ValueError, match="stm_mode"):
        cr3bp.propagate(
            arenstorf_system,
            arenstorf_state,
            1.0,
            with_stm=True,
            stm_mode="bicomplex",  # type: ignore[arg-type]
        )


def test_stm_mode_irrelevant_when_with_stm_false(
    arenstorf_system: cr3bp.CR3BPSystem, arenstorf_state: np.ndarray
) -> None:
    """When ``with_stm=False``, the ``stm_mode`` keyword is ignored (state-only).

    The fixed-path argument is meaningless without STM integration; the
    function should still return a clean state-only arc.
    """
    arc = cr3bp.propagate(
        arenstorf_system, arenstorf_state, 1.0, with_stm=False, stm_mode="fixed_path"
    )
    assert arc.stm is None
    assert arc.state_f.shape == (6,)
    assert np.all(np.isfinite(arc.state_f))
