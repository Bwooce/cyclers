"""CR3BP analytic linear-seed generator (#435 Task 1).

Validates that the collinear-Lyapunov linear seed, after the fixed-Jacobi
symmetric corrector, is a genuine planar periodic orbit (it closes on
re-propagation over one full period). Golden = closure + the linear-regime
period band; the collinear linear-Lyapunov construction is textbook
(Szebehely / Koon-Lo-Marsden-Ross), so we pin behaviour, not an unsourced
exact initial condition.

The Richardson (1980) halo-branch tests below (#580) are different: Table I
(p. 253 of the source paper) gives fully-sourced, published numeric golden
values for every Appendix I coefficient at Sun-Earth mu=3.04036e-6, L1/L2/L3
-- so those ARE exact-value goldens, traced directly to a fresh read of the
primary-source PDF (``cyclers_pdf/papers/richardson-1980-...pdf``), not to
this codebase's own computed output (see
``docs/notes/2026-07-12-digest-richardson-1980-collinear-halo-analytic.md``
section 4). NOTE: that digest's Table I transcription has several confirmed
transcription errors (a lambda_L3/k_L1 cell swap, a dropped minus sign on
c3_L2, a wrong exponent on s2, a digit transposition on l1_L3, a dropped
1e-1 scale on a24, and a copy-paste of c2_L2 into d32_L2) caught while
building this golden test by re-reading the source PDF's embedded text
layer directly (``pdftotext -layout``) and cross-checking against an
independent third-party implementation (jacobwilliams/Fortran-
Astrodynamics-Toolkit's ``halo_orbit_module.f90``, also citing Richardson
1980). The GOLDEN_TABLE_I values below are the corrected, PDF-verified
figures -- see the #580 OUTSTANDING resolution note for the full diff.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import cr3bp_system, jacobi_constant, propagate
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    correct_general_periodic_3d,
)
from cyclerfinder.search.cr3bp_seed_generator import (
    _collinear_c2,
    dro_seed,
    lyapunov_seed,
    lyapunov_seed_3d,
    richardson_halo_coefficients,
    richardson_halo_ic,
    richardson_halo_seed,
)
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x

# Sun-Earth system constants used by Richardson (1980) Table I, p. 253.
SUN_EARTH_MU = 3.04036e-6

# Sourced golden values, Table I p. 253 (see module docstring for the
# corrected-transcription provenance note). rtol=1e-4 comfortably covers the
# 6-significant-figure published precision (this implementation's actual
# max relative deviation across all 3 points x 28 fields is ~4e-6).
GOLDEN_TABLE_I: dict[str, dict[str, float]] = {
    "L1": dict(
        gamma_l=1.00109e-2,
        lam=2.08645,
        k=3.22927,
        delta=2.92214e-1,
        c2=4.06107,
        c3=3.02001,
        c4=3.03054,
        s1=-8.24661e-1,
        s2=1.21099e-1,
        l1=-1.59656e1,
        l2=1.74090,
        a21=2.09270,
        a22=2.48298e-1,
        b21=-4.92446e-1,
        b22=6.07465e-2,
        d21=-3.46865e-1,
    ),
    "L2": dict(
        gamma_l=1.00782e-2,
        lam=2.05701,
        k=3.18723,
        delta=2.90785e-1,
        c2=3.94052,
        c3=-2.97984,
        c4=2.97026,
        s1=-7.44452e-1,
        s2=1.25047e-1,
        l1=-1.48288e1,
        l2=1.67369,
        a21=-2.05304,
        a22=-2.51646e-1,
        b21=4.91357e-1,
        b22=-6.27190e-2,
        d21=3.52118e-1,
    ),
    "L3": dict(
        gamma_l=9.99998e-1,
        lam=1.00000,
        k=2.00000,
        delta=2.66029e-6,
        c2=1.00000,
        c3=1.00000,
        c4=1.00000,
        s1=-1.59141e-6,
        s2=6.29433e-6,
        l1=-1.57717e-5,
        l2=1.40258e-5,
        a21=5.00000e-1,
        a22=2.50000e-1,
        b21=-2.50000e-1,
        b22=2.49998e-1,
        d21=-4.99999e-1,
    ),
}


def test_lyapunov_seed_closes_on_earth_moon_l1() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-3)
    assert state0.shape == (6,)
    assert 2.0 < period < 3.5  # EM L1 linear-Lyapunov period band (TU)
    # Closure: re-propagate state0 one full period and compare to state0.
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-7)


def test_lyapunov_seed_jacobi_consistent() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-3)
    # Symmetric form: y=z=vx=vz=0 at the corrected start.
    assert abs(state0[1]) < 1e-12
    assert abs(state0[2]) < 1e-12
    assert abs(state0[3]) < 1e-12
    assert abs(state0[5]) < 1e-12
    # Jacobi conserved across the full period.
    c0 = jacobi_constant(state0, sysem.mu)
    arc = propagate(sysem, state0, period)
    assert abs(jacobi_constant(arc.state_f, sysem.mu) - c0) < 1e-9


def test_lyapunov_seed_arbitrary_mu_sun_earth_l1() -> None:
    # Arbitrary (small) mu path: Sun-Earth L1. At this very small mu (~3e-6) the
    # corrector converges on its |xdot(t_half)| residual (tol 1e-8) but the
    # full-period state-closure floors near ~4e-7 (the half-crossing residual
    # does not directly bound full-period closure); 1e-6 is still a tight,
    # meaningful closure and confirms the arbitrary-mu seed path.
    sysem = cr3bp_system("Sun", "Earth")
    state0, period = lyapunov_seed(sysem, point="L1", amplitude=1e-4)
    assert state0.shape == (6,)
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-6)


def test_lyapunov_seed_bad_point_raises() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError):
        lyapunov_seed(sysem, point="L9", amplitude=1e-3)


def test_lyapunov_seed_3d_earth_moon_l1_closes_nonplanar() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = lyapunov_seed_3d(sysem, point="L1", amplitude_z=0.02)
    assert state0.shape == (6,)
    assert abs(state0[2]) > 0.01  # genuinely out-of-plane, not a planar collapse
    assert period > 0
    # Independent closure: propagate one full period, assert it returns near itself.
    arc = propagate(sysem, state0, period)
    assert np.linalg.norm(arc.state_f - state0) < 1e-6


def test_dro_seed_closes_on_earth_moon() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = dro_seed(sysem, amplitude=5e-2)
    assert state0.shape == (6,)
    assert period > 0.0
    arc = propagate(sysem, state0, period)
    assert np.allclose(arc.state_f, state0, atol=1e-7)


# --- Richardson (1980) analytic halo construction (#580) -------------------


@pytest.mark.parametrize("point", ["L1", "L2", "L3"])
def test_richardson_halo_coefficients_sun_earth_table1_golden(point: str) -> None:
    """Every Appendix I coefficient reproduces Table I (p. 253) at Sun-Earth mu.

    Sourced, non-circular golden: GOLDEN_TABLE_I traces to a fresh read of
    the published paper's PDF (not to this module's own output; see module
    docstring for the corrected-transcription provenance, including 6
    confirmed digest-transcription errors caught while building this test).
    """
    c = richardson_halo_coefficients(SUN_EARTH_MU, point)
    for field, expected in GOLDEN_TABLE_I[point].items():
        actual = getattr(c, field)
        assert actual == pytest.approx(expected, rel=1e-4), (
            f"{point}.{field}: computed {actual!r} vs Table I {expected!r}"
        )


def test_richardson_halo_coefficients_bad_point_raises() -> None:
    with pytest.raises(ValueError):
        richardson_halo_coefficients(SUN_EARTH_MU, "L9")


def test_richardson_halo_ic_bad_branch_raises() -> None:
    sysem = cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError):
        richardson_halo_ic(sysem.mu, "L1", 0.08, branch="III")


def test_richardson_halo_ic_phi_zero_is_exact_perpendicular_crossing() -> None:
    """At phi=0 every sin(n*phi) term in Eq. 20 vanishes identically, so
    y0 = xdot0 = zdot0 = 0 EXACTLY -- the analytic construction's own reason
    for defaulting to phi=0 when feeding the tulip/NRHO-style corrector
    (FREE_VARS_SYMMETRIC_TULIP assumes this exact perpendicular-crossing
    form)."""
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = richardson_halo_ic(sysem.mu, "L1", 0.08, branch="I", phi=0.0)
    assert state0.shape == (6,)
    assert state0[1] == 0.0
    assert state0[3] == 0.0
    assert state0[5] == 0.0
    assert period > 0.0


@pytest.mark.parametrize("point,branch", [("L1", "I"), ("L1", "II"), ("L2", "I"), ("L2", "II")])
def test_richardson_halo_seed_closes_nonplanar_earth_moon(point: str, branch: str) -> None:
    """The Richardson-seeded corrector converges to a genuinely non-planar
    (|z0| significantly nonzero) periodic orbit, and the IC closes on
    independent re-propagation over one full period."""
    sysem = cr3bp_system("Earth", "Moon")
    state0, period = richardson_halo_seed(sysem, point=point, amplitude_z=0.08, branch=branch)
    assert state0.shape == (6,)
    assert abs(state0[2]) > 0.05  # genuinely out-of-plane halo, not a planar collapse
    assert period > 0
    arc = propagate(sysem, state0, period)
    assert np.linalg.norm(arc.state_f - state0) < 1e-6


def test_richardson_halo_seed_small_amplitude_raises_not_silently_planar() -> None:
    """Below Richardson's own Ax_minimum bifurcation (Eq. 22), NO genuine
    halo exists -- the true periodic member at that fixed x0 really is
    planar. richardson_halo_seed must raise rather than silently return the
    (degenerate) planar solution as if it were a halo."""
    sysem = cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError):
        richardson_halo_seed(sysem, point="L1", amplitude_z=0.02, branch="I")


def test_richardson_halo_seed_beats_linear_seed_at_large_amplitude() -> None:
    """Fable-mandated seed comparison (#580): don't assume the payoff, measure
    it. At amplitude_z=0.08 -- in the regime the existing
    :func:`lyapunov_seed_3d` ladder's own comment documents as fragile
    ("non-monotonic in the vertical amplitude") -- feed BOTH seeds' analytic
    ICs directly to the SAME corrector call (no ladder fallback for either,
    so the comparison is apples-to-apples at one fixed amplitude), and
    compare convergence to a genuinely non-planar orbit.

    Measured result (honest, not forced positive): at this amplitude the
    linear vertical-Lyapunov seed's single-shot corrector call does NOT
    reach a genuinely non-planar solution (it collapses to the planar
    member), while the Richardson third-order seed converges cleanly in a
    handful of iterations. This is a genuine, reproducible advantage in the
    documented-fragile large-amplitude regime -- not true at every
    amplitude (both seeds are comparable or both fail at very small/very
    large amplitude_z; see the #580 OUTSTANDING resolution note for the
    fuller amplitude scan), but real in the regime this test targets.
    """
    sysem = cr3bp_system("Earth", "Moon")
    amplitude_z = 0.08

    # Linear vertical-Lyapunov seed, IDENTICAL construction to
    # lyapunov_seed_3d's internal seed (reusing its own helpers), fed
    # directly to the corrector at this ONE amplitude -- no ladder retry.
    x_l = lagrange_collinear_x(sysem.mu, "L1")
    c2 = _collinear_c2(x_l, sysem.mu)
    omega_z = math.sqrt(c2)
    linear_state0 = np.array([x_l, 0.0, amplitude_z, 0.0, 1e-3, 0.0], dtype=np.float64)
    linear_orbit = correct_general_periodic_3d(
        sysem,
        state0_guess=linear_state0,
        period_guess=2.0 * math.pi / omega_z,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    linear_ok = linear_orbit.converged and not linear_orbit.degenerate_planar

    richardson_state0, richardson_period = richardson_halo_ic(
        sysem.mu, "L1", amplitude_z, branch="I", phi=0.0
    )
    richardson_orbit = correct_general_periodic_3d(
        sysem,
        state0_guess=richardson_state0,
        period_guess=richardson_period,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    richardson_ok = richardson_orbit.converged and not richardson_orbit.degenerate_planar

    assert not linear_ok, "expected the linear seed to NOT reach a non-planar solution here"
    assert richardson_ok, "expected the Richardson seed to converge to a non-planar halo here"
    assert richardson_orbit.n_iter < 15
    assert abs(richardson_orbit.state0[2]) > 0.05
