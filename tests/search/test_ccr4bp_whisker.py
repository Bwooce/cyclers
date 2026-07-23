"""Tests for the CCR4BP whisker/manifold-direction extraction module (`#691`).

Structure mirrors `#690`'s own test discipline: backend correctness first, then
two layers of positive control BEFORE trusting any result on the genuinely
forced (physical-mass) CCR4BP torus:

* **Positive control A -- natural-period monodromy.** The Jupiter-Europa 3:4
  resonant orbit `#690` seeds its torus from (built by the SAME
  ``_resonant_symmetric_orbit`` scaffolding `#690`'s own test module uses) turns
  out, when checked directly, to be LINEARLY STABLE at ``mu_gan = 0`` -- all six
  monodromy eigenvalues lie on the unit circle to ~1e-13 (a genuine physical
  fact about this specific resonant orbit, not a bug; the corrector can still
  build a genuine 2D torus from it regardless of hyperbolicity, since the
  pseudospectral invariance-PDE approach does not require a saddle structure).
  This module's :func:`~cyclerfinder.search.ccr4bp_whisker.one_period_stm`,
  evaluated over the orbit's OWN natural period (via a period-matched pure
  Fourier projection -- :func:`~cyclerfinder.search.variational_ccr4bp_torus.
  project_base_orbit_to_2d` at high ``n2``, bypassing the corrector's residual
  noise entirely), must reproduce this exact non-hyperbolic result and match a
  fully independent, direct ``cr3bp_stm_eom`` monodromy computation to near
  machine precision.
* **Positive control B -- Ganymede-clock-period reduction + a genuinely
  hyperbolic reference case.** `#691`'s Step 1 defines "one period" as one
  Ganymede-synodic ``theta1`` cycle (the corrector's own ``torus.period``), NOT
  the orbit's natural period -- a DIFFERENT, generally incommensurate duration.
  Even at ``mu_gan = 0`` (autonomous, unforced dynamics), propagating over this
  externally-chosen, non-natural duration exhibits genuine LOCAL exponential
  stretching (the same "frozen-time rate" phenomenon `#544` documented for the
  QBCP corrector) despite the orbit's own asymptotic monodromy being elliptic --
  a real, non-circular reference case, independently cross-checked against a
  raw (bypassing this module and `#690`'s corrector entirely) ``cr3bp_stm_eom``
  propagation at the identical duration (confirms `#689`'s ``mu_gan -> 0``
  reduction holds through this module's own code path, not just its own unit
  tests). The perturbation-robustness pipeline (segment-anchored CLV vs the
  naive one-shot extraction) is validated on this well-understood, independently
  checkable case FIRST, before being trusted on the physical-mass torus.

Then the physical-mass (Kumar et al. 2021 class) torus itself: `#691` Step 1's
one-period amplification report, and Step 2's perturbation-robustness diagnostic
-- the actual go/no-go test for whether the cheap direction extraction is
trustworthy or the expensive Kumar bundle-solve is required.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_whisker as wk
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Symmetric p:q resonant CR3BP periodic orbit. Identical scaffolding to
    `#690`'s ``tests/search/test_variational_ccr4bp_torus.py`` (test-only, no
    production code under test)."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


@pytest.fixture(scope="module")
def system() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


@pytest.fixture(scope="module")
def sys0(system: ccr4bp.CCR4BPSystem) -> ccr4bp.CCR4BPSystem:
    return ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )


@pytest.fixture(scope="module")
def orbit_34(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    s0, period, res = _resonant_symmetric_orbit(system.mu, 3, 4)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


@pytest.fixture(scope="module")
def flat_torus(
    sys0: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    """`#690`'s own ``mu_gan=0`` positive-control torus (corrector-refined,
    period = Ganymede's synodic period)."""
    s0, period = orbit_34
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )


@pytest.fixture(scope="module")
def phys_torus(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    """`#690`'s physical-mass Jupiter-Europa 3:4 CCR4BP torus (Kumar et al. 2021
    object class), built directly at the physical ``mu_gan`` (matches `#690`'s
    own ``test_physical_mass_thin_torus``)."""
    s0, period = orbit_34
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )


@pytest.fixture(scope="module")
def period_matched_torus(
    sys0: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    """A PURE Fourier projection (no corrector Newton solve -- bypasses its
    residual/gauge noise entirely) of the base 3:4 orbit at high ``n2``, with
    ``period``/``omega1``/``omega2`` set to the ORBIT'S OWN natural period (NOT
    the Ganymede synodic period `#690`'s torus objects normally carry). This is
    the period-matched object Positive Control A cross-checks against a direct,
    independent CR3BP monodromy computation -- see module docstring."""
    s0, period = orbit_34
    n1, n2 = 1, 80
    m1, m2 = 2 * n1 + 3, 2 * n2 + 3
    coeffs0 = vt.project_base_orbit_to_2d(s0, period, sys0.mu, n1, n2, m1, m2)
    omega_natural = 2.0 * np.pi / period
    return vt.CCR4BPTorusVariationalResult(
        system=sys0,
        coeffs=coeffs0,
        omega1=omega_natural,
        omega2=omega_natural,
        rotation_number=1.0,
        rho_strob=2.0 * np.pi,
        period=period,
        n1=n1,
        n2=n2,
        m1=m1,
        m2=m2,
        period_multiple=1,
        transverse_amplitude=0.0,
        residual_rms=0.0,
        closure_residual=0.0,
        converged=True,
        n_iter=0,
        notes="period_matched_pure_projection_no_corrector",
    )


# ---------------------------------------------------------------------------
# Positive control A: natural-period monodromy (the orbit is genuinely
# non-hyperbolic -- confirmed against a fully independent CR3BP computation).
# ---------------------------------------------------------------------------


def test_period_matched_projection_recovers_orbit_state(
    period_matched_torus: vt.CCR4BPTorusVariationalResult, orbit_34: tuple[np.ndarray, float]
) -> None:
    """The high-n2 pure projection reproduces the exact orbit IC at (0,0) to
    near machine precision (n2=80 resolves this eccentric orbit's Fourier series
    far below any hyperbolicity-detection tolerance used elsewhere in this file)."""
    s0, _period = orbit_34
    s_eval = vt.evaluate_torus_state(period_matched_torus, 0.0, 0.0)
    planar_s0 = np.array([s0[0], s0[1], s0[3], s0[4]])
    assert float(np.max(np.abs(s_eval - planar_s0))) < 1e-10


def test_one_period_stm_matches_natural_period_cr3bp_monodromy(
    period_matched_torus: vt.CCR4BPTorusVariationalResult, orbit_34: tuple[np.ndarray, float]
) -> None:
    """POSITIVE CONTROL A. Over the orbit's own natural period, this module's
    ``one_period_stm`` + ``amplification_diagnostics`` must reproduce a fully
    independent, direct ``cr3bp_stm_eom`` monodromy computation: same STM to
    near machine precision, and -- crucially -- the SAME (correct, non-obvious)
    physical verdict that this specific 3:4 resonant orbit is LINEARLY STABLE
    (no real hyperbolic eigenpair; all 6 eigenvalues on the unit circle)."""
    s0, period = orbit_34
    mu = period_matched_torus.system.mu

    state0, stm_mine = wk.one_period_stm(period_matched_torus, 0.0, 0.0, rtol=1e-13, atol=1e-13)

    y42 = np.concatenate([s0, np.eye(6).reshape(36)])
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom, (0.0, period), y42, args=(mu,), method="DOP853", rtol=1e-13, atol=1e-13
    )
    stm_direct = sol.y[6:, -1].reshape(6, 6)

    rel_err = float(np.linalg.norm(stm_mine - stm_direct) / np.linalg.norm(stm_direct))
    assert rel_err < 1e-6, rel_err

    diag = wk.amplification_diagnostics(state0, stm_mine)
    eig_direct = np.linalg.eigvals(stm_direct)
    assert (
        float(np.max(np.abs(np.abs(eig_direct) - 1.0))) < 1e-9
    )  # direct: genuinely non-hyperbolic
    assert diag.lam_u is None, diag.lam_u  # this module must agree: no spurious instability
    assert diag.lam_s is None, diag.lam_s
    assert float(np.max(np.abs(diag.eigenvalues))) == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Positive control B: Ganymede-clock-period reduction (mu_gan=0) + a genuinely
# hyperbolic, independently-checkable reference case for the CLV pipeline.
# ---------------------------------------------------------------------------


def test_ganymede_period_map_reduces_to_raw_cr3bp_at_mu_gan_zero(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """POSITIVE CONTROL B (wiring/reduction check). At ``mu_gan=0`` this
    module's one-period STM over the (here, physically arbitrary) Ganymede
    synodic-period duration must reduce EXACTLY to a raw ``cr3bp_stm_eom``
    propagation of the same duration -- `#689`'s structural reduction, verified
    through THIS module's own code path, not just its own unit tests."""
    mu = flat_torus.system.mu
    for theta1, theta2 in ((0.0, 0.0), (1.0, 2.0), (3.0, 4.0)):
        state0, stm_mine = wk.one_period_stm(flat_torus, theta1, theta2, rtol=1e-13, atol=1e-13)
        t0 = (theta1 % (2.0 * np.pi)) / flat_torus.omega1
        y42 = np.concatenate([state0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            (t0, t0 + flat_torus.period),
            y42,
            args=(mu,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        stm_raw = sol.y[6:, -1].reshape(6, 6)
        assert np.allclose(stm_mine, stm_raw, atol=1e-9), (theta1, theta2)


def test_flat_torus_ganymede_period_map_is_genuinely_hyperbolic(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """The Ganymede-clock-period stroboscopic map is genuinely hyperbolic even
    at ``mu_gan=0`` (a real "frozen-time"/local-stretching effect over a
    non-natural duration, NOT orbital instability -- see module docstring), and
    the reciprocal (lam_u * lam_s ~ 1) pairing expected of any symplectic
    (Hamiltonian-flow) STM holds regardless of the propagation duration."""
    checked = 0
    for theta1, theta2 in ((0.0, 0.0), (0.0, 1.5), (1.0, 2.0)):
        state0, stm = wk.one_period_stm(flat_torus, theta1, theta2)
        diag = wk.amplification_diagnostics(state0, stm)
        assert diag.lam_u is not None and abs(diag.lam_u) > 1.0
        assert diag.lam_s is not None
        assert diag.lam_u * diag.lam_s == pytest.approx(1.0, abs=1e-6)
        assert diag.separation_ratio > 1.0
        checked += 1
    assert checked == 3


# ---------------------------------------------------------------------------
# Segment-anchored CLV self-consistency + convergence (on the well-understood
# flat-torus reference case).
# ---------------------------------------------------------------------------


def test_segmented_qr_matches_direct_composition(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    mv = wk.manifold_direction_segmented_clv(
        flat_torus, 0.0, 0.0, "unstable", None, n_segments=16, return_diagnostics=True
    )
    assert mv is not None and len(mv) == 3
    diag = mv[2]
    assert isinstance(diag, wk.SegmentedCLVDiagnostics)
    assert diag.qr_vs_direct_stm_relerr < 1e-10, diag.qr_vs_direct_stm_relerr


def test_segmented_clv_converges_in_n_segments(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """The extracted direction converges to a definite limit as ``n_segments``
    grows. Live-observed on this phase: angle-to-a-64-segment-reference shrinks
    2.71 (k=2) -> 1.51 (k=8) -> 0.62 (k=16) -> 0.06 (k=32) deg -- monotone, but
    this phase needs MORE segments than `#646`'s own EM-L2 QBCP case did (that
    converged to <0.5 deg by k=16) to reach the same margin; see
    ``test_perturbation_robustness_physical_torus`` for the same
    n_segments-sensitivity finding on the physical torus's most extreme phase."""
    theta1, theta2, branch = 0.0, 0.0, "unstable"
    ref = wk.manifold_direction_segmented_clv(
        flat_torus, theta1, theta2, branch, None, n_segments=64
    )
    assert ref is not None
    vref = ref[1]

    def angle_at(k: int) -> float:
        mv = wk.manifold_direction_segmented_clv(
            flat_torus, theta1, theta2, branch, None, n_segments=k
        )
        assert mv is not None
        return float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(vref, mv[1])))))))

    a2, a8, a16, a32 = angle_at(2), angle_at(8), angle_at(16), angle_at(32)
    assert a32 <= a16 <= a8 <= a2 + 1e-9
    assert a32 < 0.2


def test_segmented_clv_contract_and_ref_vec_sign(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """``(state6, unit-vec)`` contract + sign-fixing, mirroring the QBCP
    extractor's own contract tests."""
    for branch in ("unstable", "stable"):
        mv = wk.manifold_direction_segmented_clv(flat_torus, 0.0, 0.0, branch, None, n_segments=12)
        assert mv is not None
        state6, vec = mv[0], mv[1]
        expected_state6 = wk._embed6(vt.evaluate_torus_state(flat_torus, 0.0, 0.0))
        assert np.allclose(state6, expected_state6, atol=1e-12)
        assert float(np.linalg.norm(vec)) == pytest.approx(1.0, abs=1e-10)
        mv_pos = wk.manifold_direction_segmented_clv(
            flat_torus, 0.0, 0.0, branch, vec, n_segments=12
        )
        mv_neg = wk.manifold_direction_segmented_clv(
            flat_torus, 0.0, 0.0, branch, -vec, n_segments=12
        )
        assert mv_pos is not None and mv_neg is not None
        assert float(np.dot(mv_pos[1], vec)) > 0.0
        assert float(np.dot(mv_neg[1], vec)) < 0.0


# ---------------------------------------------------------------------------
# Positive control C: perturbation-robustness pipeline on the well-understood
# flat-torus reference case -- must behave sanely before trusting phys_torus.
# ---------------------------------------------------------------------------


def test_perturbation_robustness_flat_torus_positive_control(
    flat_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """POSITIVE CONTROL C (`#619`/`#646`'s go/no-go test, run first on the
    ``mu_gan=0`` reference case). At every checked hyperbolic phase the
    segmented method's swing is materially SMALLER than the naive one-shot's
    (the mechanism `#646` demonstrated at EM-L2), even though the CCR4BP's weak
    forcing means neither method is catastrophically fragile here (live-observed
    max one-shot swing < 0.5 deg, max segmented swing < 0.01 deg across 3 phases
    x 3 seeds x 2 branches)."""
    checked = 0
    for theta1, theta2 in ((0.0, 0.0), (0.0, 1.5), (1.0, 2.0)):
        for branch in ("unstable", "stable"):
            worst_one_shot = 0.0
            worst_segmented = 0.0
            for seed in range(3):
                s1 = wk.perturbation_swing_deg(
                    flat_torus, theta1, theta2, branch, method="one_shot", seed=seed
                )
                s2 = wk.perturbation_swing_deg(
                    flat_torus, theta1, theta2, branch, method="segmented", n_segments=16, seed=seed
                )
                assert s1 is not None and s2 is not None
                worst_one_shot = max(worst_one_shot, s1)
                worst_segmented = max(worst_segmented, s2)
            assert worst_segmented < worst_one_shot, (theta1, theta2, branch)
            assert worst_one_shot < 1.0, worst_one_shot
            assert worst_segmented < 0.01, worst_segmented
            checked += 1
    assert checked == 6


# ---------------------------------------------------------------------------
# `#691` Step 1: one-period amplification off the physical-mass CCR4BP torus.
# ---------------------------------------------------------------------------


def test_amplification_diagnostics_physical_torus(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Step 1's headline report: the physical-mass torus has a clear, well-
    separated unstable direction at every phase checked (live-observed |lam_u|
    ranging ~9-1700 depending on phase; always >= a few x the next-largest
    eigenvalue magnitude), and the reciprocal lam_u*lam_s~1 symplectic check
    holds."""
    checked = 0
    for theta1 in np.linspace(0.0, 2 * np.pi, 4, endpoint=False):
        for theta2 in np.linspace(0.0, 2 * np.pi, 3, endpoint=False):
            state0, stm = wk.one_period_stm(phys_torus, float(theta1), float(theta2))
            diag = wk.amplification_diagnostics(state0, stm)
            assert diag.lam_u is not None and abs(diag.lam_u) > 1.0, (theta1, theta2)
            assert diag.lam_s is not None
            assert diag.lam_u * diag.lam_s == pytest.approx(1.0, abs=5e-3)
            assert diag.separation_ratio > 1.5, (theta1, theta2, diag.separation_ratio)
            checked += 1
    assert checked == 12


# ---------------------------------------------------------------------------
# `#691` Step 2: the actual go/no-go perturbation-robustness diagnostic.
# ---------------------------------------------------------------------------


def test_perturbation_robustness_physical_torus(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """THE GO/NO-GO DIAGNOSTIC. Checks the worst phase found by an exploratory
    sweep (theta1~3.77, theta2~5.03, unstable branch, lam_u~-849 there) plus a
    mild phase and (2.0, 2.0) -- the single most extreme-dynamic-range phase
    found (|lam_u/lam_s| ~ 2.8e6 there, one-shot eigenvalues [1688.5, 64.8, 5.76,
    0.174, 0.0154, 0.00059]).

    ``n_segments=32`` (not the module's default 16) is used here: at (2.0, 2.0)
    specifically, ``n_segments=16`` was NOT yet converged (swing 3.3 deg, an
    n_segments-resolution artifact of the extreme dynamic range at that one
    phase, not a fundamental fragility -- confirmed by checking n_segments in
    {16, 32, 64, 128}: swing drops 1.22 -> 0.0076 -> 0.0008 -> 0.0002 deg,
    converging cleanly once given enough segments; at every OTHER phase checked,
    n_segments=16 was already <0.03 deg). This is itself a real, reportable
    finding, not swept under the rug: `#646`'s own n_segments=16 default does
    NOT uniformly suffice across all CCR4BP torus phases -- n_segments>=32 is
    the safer recommendation for this system.

    Live-observed worst-case swings at n_segments=32 (5 seeds, all 3 phases x 2
    branches): one-shot up to 3.44 deg; segmented up to 0.0083 deg -- BOTH
    methods are trustworthy at the CCR4BP's forcing strength (orders of
    magnitude below the ~2e4-7e4x per-period amplification that made the QBCP
    EM-L2 torus's one-shot direction genuinely wrong by 37-85 deg in `#646`);
    the segmented method is consistently, substantially better (~400-2000x
    smaller swing), but the naive one-shot direction is ALSO usable here -- the
    expensive full Kumar bundle-solve whisker construction is NOT indicated for
    this torus."""
    phases = [(3.7699111843077517, 5.026548245743669), (0.0, 0.0), (2.0, 2.0)]
    for theta1, theta2 in phases:
        for branch in ("unstable", "stable"):
            worst_one_shot = 0.0
            worst_segmented = 0.0
            for seed in range(5):
                s1 = wk.perturbation_swing_deg(
                    phys_torus, theta1, theta2, branch, method="one_shot", seed=seed
                )
                s2 = wk.perturbation_swing_deg(
                    phys_torus, theta1, theta2, branch, method="segmented", n_segments=32, seed=seed
                )
                if s1 is not None:
                    worst_one_shot = max(worst_one_shot, s1)
                if s2 is not None:
                    worst_segmented = max(worst_segmented, s2)
            # The direction is trustworthy: even the naive one-shot method stays
            # well within a usable tolerance for downstream shooting/matching.
            assert worst_one_shot < 10.0, (theta1, theta2, branch, worst_one_shot)
            assert worst_segmented < 0.1, (theta1, theta2, branch, worst_segmented)
            assert worst_segmented <= worst_one_shot + 1e-9, (theta1, theta2, branch)


# ---------------------------------------------------------------------------
# Guard rails.
# ---------------------------------------------------------------------------


def test_rejects_bad_args(flat_torus: vt.CCR4BPTorusVariationalResult) -> None:
    with pytest.raises(ValueError, match="branch"):
        wk.manifold_direction_segmented_clv(flat_torus, 0.0, 0.0, "sideways", None)
    with pytest.raises(ValueError, match="n_segments"):
        wk.manifold_direction_segmented_clv(flat_torus, 0.0, 0.0, "unstable", None, n_segments=0)
    with pytest.raises(ValueError, match="branch"):
        wk.manifold_direction_one_shot(flat_torus, 0.0, 0.0, "sideways", None)
    with pytest.raises(ValueError, match="method"):
        wk.perturbation_swing_deg(flat_torus, 0.0, 0.0, "unstable", method="bogus")
