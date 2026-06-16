"""BCR4BP core EOM / STM / propagator tests (#292 Phase 1).

Four gates:

  1. CR3BP-limit (STRUCTURAL): at mu_S = 0 the BCR4BP EOM matches
     :func:`cyclerfinder.core.cr3bp.cr3bp_eom` to floating-point precision
     at sampled states. Same for the STM EOM.

  2. Indirect-term cancellation at the EM barycenter: at r = 0 the net Sun
     acceleration is zero (direct + indirect cancel by construction). This
     is the physical meaning of the indirect term.

  3. STM consistency: variational integration of the STM (6x6) agrees with
     finite-difference Jacobian at 1e-5 over a moderate horizon.

  4. Andreu / Rosales-Jorba parameter constants: the four published
     constants are wired through :func:`andreu_default` exactly to the
     digest values.

All tests assert STRUCTURE (CR3BP limit, FD<->STM agreement) or sourced
parameter VALUES (digest table). No value computed by our own code is on the
EXPECTED side of any assertion (per ``feedback_golden_tests_sourced_only``).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp

# ---------------------------------------------------------------------------
# Sample states used across tests. Picked to be away from primaries (no
# close-approach singularities) and out of the planar manifold (so z-coupling
# is exercised).
# ---------------------------------------------------------------------------

_SAMPLE_STATE_PLANAR = np.array([0.5, 0.1, 0.0, 0.0, 0.3, 0.0], dtype=np.float64)
_SAMPLE_STATE_3D = np.array([0.5, 0.1, 0.05, 0.02, 0.3, -0.01], dtype=np.float64)
_SAMPLE_STATE_ANDREU_POL1 = np.array(
    [-0.8369141677649317, 0.0, 0.0, 0.0, -0.8391311559808445, 0.0], dtype=np.float64
)


# ---------------------------------------------------------------------------
# Gate 1: CR3BP-limit (mu_S = 0) -- structural correctness anchor.
# ---------------------------------------------------------------------------


def test_cr3bp_limit_eom_planar() -> None:
    """At mu_S = 0 the BCR4BP planar EOM matches CR3BP to floating-point precision."""
    mu = 0.012150581600000
    sys_bcr = bcr4bp.BCR4BPSystem(mu=mu, mu_sun=0.0, a_sun_nondim=388.0, omega_sun_nondim=0.9252)
    rhs_bcr = bcr4bp.bcr4bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys_bcr)
    rhs_cr3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE_PLANAR, mu)
    assert np.allclose(rhs_bcr, rhs_cr3, rtol=0.0, atol=1e-14), (
        f"BCR4BP@mu_S=0 deviates from CR3BP: max|delta|={np.max(np.abs(rhs_bcr - rhs_cr3)):.3e}"
    )


def test_cr3bp_limit_eom_3d() -> None:
    """At mu_S = 0 the BCR4BP 3D EOM matches CR3BP to floating-point precision."""
    mu = 0.012150581600000
    sys_bcr = bcr4bp.BCR4BPSystem(mu=mu, mu_sun=0.0, a_sun_nondim=388.0, omega_sun_nondim=0.9252)
    # Time non-trivial too: the Sun-phase argument should be ignored at mu_S=0.
    rhs_bcr = bcr4bp.bcr4bp_eom(1.7, _SAMPLE_STATE_3D, sys_bcr)
    rhs_cr3 = cr3bp.cr3bp_eom(1.7, _SAMPLE_STATE_3D, mu)
    assert np.allclose(rhs_bcr, rhs_cr3, rtol=0.0, atol=1e-14)


def test_cr3bp_limit_stm() -> None:
    """At mu_S = 0 the BCR4BP STM EOM matches CR3BP STM EOM to floating-point precision."""
    mu = 0.012150581600000
    sys_bcr = bcr4bp.BCR4BPSystem(mu=mu, mu_sun=0.0, a_sun_nondim=388.0, omega_sun_nondim=0.9252)
    y42 = np.concatenate([_SAMPLE_STATE_3D, np.eye(6).reshape(36)])
    rhs_bcr = bcr4bp.bcr4bp_stm_eom(0.5, y42, sys_bcr)
    rhs_cr3 = cr3bp.cr3bp_stm_eom(0.5, y42, mu)
    assert np.allclose(rhs_bcr, rhs_cr3, rtol=0.0, atol=1e-13), (
        f"BCR4BP STM@mu_S=0 deviates from CR3BP STM: "
        f"max|delta|={np.max(np.abs(rhs_bcr - rhs_cr3)):.3e}"
    )


def test_cr3bp_limit_propagator() -> None:
    """At mu_S = 0 the BCR4BP propagator agrees with CR3BP over a finite arc."""
    mu = 0.012150581600000
    sys_bcr = bcr4bp.BCR4BPSystem(mu=mu, mu_sun=0.0, a_sun_nondim=388.0, omega_sun_nondim=0.9252)
    sys_cr3 = cr3bp.CR3BPSystem(
        mu=mu, primary="earth", secondary="moon", l_km=384400.0, t_s=375190.0
    )
    arc_bcr = bcr4bp.propagate_bcr4bp(sys_bcr, _SAMPLE_STATE_3D, 1.0)
    arc_cr3 = cr3bp.propagate(sys_cr3, _SAMPLE_STATE_3D, 1.0)
    assert np.allclose(arc_bcr.state_f, arc_cr3.state_f, rtol=1e-12, atol=1e-12), (
        f"BCR4BP@mu_S=0 propagation deviates from CR3BP after t=1: "
        f"max|delta|={np.max(np.abs(arc_bcr.state_f - arc_cr3.state_f)):.3e}"
    )


# ---------------------------------------------------------------------------
# Gate 2: Indirect-term cancellation at the EM barycenter.
# ---------------------------------------------------------------------------


def test_indirect_term_cancels_at_em_barycenter() -> None:
    """At r = 0 the direct + indirect Sun contributions cancel exactly.

    Physical meaning: the synodic frame is centred on the EM barycenter, so
    the *synodic-frame* Sun acceleration on a particle co-located with the
    barycenter must be zero (that particle moves with the barycenter, which
    is in free fall toward the Sun). The indirect term is engineered to
    enforce this; if it does not cancel at r = 0 the term is wrong.
    """
    sys_bcr = bcr4bp.andreu_default()
    # Evaluate at a few different Sun phases.
    for t in [0.0, 0.3 * sys_bcr.sun_period_tu, 0.7 * sys_bcr.sun_period_tu]:
        ax, ay, az = bcr4bp._sun_acceleration(0.0, 0.0, 0.0, t, sys_bcr)
        assert abs(ax) < 1e-10 and abs(ay) < 1e-10 and abs(az) < 1e-10, (
            f"Indirect cancellation failed at t={t}: (ax,ay,az)=({ax:.3e}, {ay:.3e}, {az:.3e})"
        )


# ---------------------------------------------------------------------------
# Gate 3: STM finite-difference consistency.
# ---------------------------------------------------------------------------


def test_stm_finite_difference_consistency() -> None:
    """Variational STM matches the central-finite-difference Jacobian.

    Propagate the state for a moderate horizon under BCR4BP with STM. Then
    propagate 12 perturbed states (eps in each component, both signs) and
    compare central-difference columns
    ``(X(T;X0+eps*e_i) - X(T;X0-eps*e_i))/(2*eps)`` to STM[:,i].

    Central FD floor scales as eps^2; with eps = 1e-5 the floor is ~1e-10
    AND the truncation error is ~ |D^3 X / D X^3| * eps^2 ~ O(1) * 1e-10.
    We gate at 1e-5 (relative tolerance for the largest entry, which can be
    ~10 over half a TU under Sun-perturbed dynamics).
    """
    sys_bcr = bcr4bp.andreu_default()
    state0 = _SAMPLE_STATE_3D.copy()
    t_horizon = 0.3  # short enough to keep nonlinearity controllable

    arc = bcr4bp.propagate_bcr4bp(sys_bcr, state0, t_horizon, with_stm=True)
    assert arc.stm is not None
    stm = arc.stm

    eps = 1e-5
    fd_jac = np.zeros((6, 6), dtype=np.float64)
    for i in range(6):
        state_p = state0.copy()
        state_p[i] += eps
        state_m = state0.copy()
        state_m[i] -= eps
        arc_p = bcr4bp.propagate_bcr4bp(sys_bcr, state_p, t_horizon)
        arc_m = bcr4bp.propagate_bcr4bp(sys_bcr, state_m, t_horizon)
        fd_jac[:, i] = (arc_p.state_f - arc_m.state_f) / (2.0 * eps)

    # Compare relative to the scale of the STM itself: STM entries can reach
    # O(10) for divergent directions even over t=0.3, so a relative tol is
    # the meaningful check.
    scale = max(1.0, float(np.max(np.abs(stm))))
    rel_delta = float(np.max(np.abs(fd_jac - stm))) / scale
    assert rel_delta < 1e-5, (
        f"STM vs central-FD relative disagreement: rel_delta={rel_delta:.3e} > 1e-5 "
        f"(absolute max|delta|={np.max(np.abs(fd_jac - stm)):.3e}, "
        f"stm_scale={scale:.3e})"
    )


# ---------------------------------------------------------------------------
# Gate 4: Sourced parameter constants from the Andreu / Rosales-Jorba digest.
# ---------------------------------------------------------------------------


def test_andreu_default_parameter_constants() -> None:
    """The four BCR4BP constants in :func:`andreu_default` trace exactly to the digest.

    EXPECTED side values are from docs/notes/2026-06-14-andreu-quasi-bicircular-digest.md
    "Sourced parameters (candidate goldens)" table -- the Rosales/Jorba (2023)
    Table 3 column. Per ``feedback_golden_tests_sourced_only`` the EXPECTED
    side is published, not code-computed.
    """
    sys_bcr = bcr4bp.andreu_default()
    # Rosales/Jorba 2023 Table 3 (digest):
    assert sys_bcr.mu == 0.012150581600000
    assert sys_bcr.mu_sun == 328900.5423094043
    assert sys_bcr.a_sun_nondim == 388.8111430233511
    assert sys_bcr.omega_sun_nondim == 0.925195985520347
    # Cross-check vs the digest's Gimeno 2018 column at *printed precision*:
    # mu agrees to 5e-12, mu_S to ~1e-2 (different rounding), a_S to 5e-13,
    # omega_S to 3e-12. These are within-table consistency checks.
    assert abs(sys_bcr.mu - 0.012150581623433623) < 3e-11
    assert abs(sys_bcr.a_sun_nondim - 388.81114302335106) < 1e-12
    assert abs(sys_bcr.omega_sun_nondim - 0.92519598551829646) < 1e-11
    # Sun synodic period in EM TU: digest says ~6.7912.
    assert 6.79 < sys_bcr.sun_period_tu < 6.80


def test_sun_commensurate_period_n1() -> None:
    """n=1 commensurate period equals 2*pi/omega_S (one Sun-synodic revolution)."""
    sys_bcr = bcr4bp.andreu_default()
    t_n1 = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)
    assert abs(t_n1 - sys_bcr.sun_period_tu) < 1e-14
    assert t_n1 == pytest.approx(2.0 * np.pi / sys_bcr.omega_sun_nondim, rel=0.0, abs=1e-14)


def test_sun_commensurate_period_rejects_nonpositive() -> None:
    sys_bcr = bcr4bp.andreu_default()
    with pytest.raises(ValueError):
        bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=0)
    with pytest.raises(ValueError):
        bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=-1)


# ---------------------------------------------------------------------------
# Bonus: Sun position is at distance a_S and on a planar circle.
# ---------------------------------------------------------------------------


def test_sun_position_is_planar_circle() -> None:
    """Sun position lies on a planar circle of radius a_sun in the synodic frame."""
    sys_bcr = bcr4bp.andreu_default()
    for t in np.linspace(0.0, sys_bcr.sun_period_tu, 13):
        sx, sy, sz = bcr4bp._sun_position(float(t), sys_bcr)
        radius = float(np.sqrt(sx * sx + sy * sy))
        assert sz == 0.0
        assert abs(radius - sys_bcr.a_sun_nondim) < 1e-12


def test_sun_position_period() -> None:
    """Sun returns to its initial position after exactly one synodic period."""
    sys_bcr = bcr4bp.andreu_default()
    s_init = bcr4bp._sun_position(0.0, sys_bcr)
    s_final = bcr4bp._sun_position(sys_bcr.sun_period_tu, sys_bcr)
    assert np.allclose(s_init, s_final, atol=1e-10)
