"""Tests for `#738`'s Radau independent-integrator ghost-guard cross-check of a converged
N=5 CRNBP quasi-periodic torus (:mod:`cyclerfinder.search.crnbp_torus_ghost_guard`).

Structure, mirroring `#701`'s ``test_ccr4bp_heteroclinic_search.py`` convention: cheap
mechanism/guard-rail tests first (a tiny, DELIBERATELY UNCONVERGED n1=1/n2=2 toy torus --
its physical quality is irrelevant to what these tests check, only that
``correct_crnbp_torus_pseudospectral`` populates a real ``closure_residual`` fast), then the
module's own real evidence-grade positive control: the ACTUAL physical-mass n1=2/n2=20
Jupiter-Europa 3:4 N=5 torus (`#720`/`#723`/`#724`'s own object, the one written to
``europa-3-4-crnbp-torus-jupiter-2026``), reconstructed here via the identical pipeline
``tests/search/test_variational_crnbp_torus.py::test_continuation_reaches_physical_mu_io_at_n1_2``
uses. That reconstruction is NOT marked ``slow`` (this project's own
``feedback_delegation_fresh_agent_not_fork`` / `#701` precedent: a V-gauntlet evidence test
must not be hidden behind a marker CI skips by default) -- it costs several minutes, matching
the sibling file's own already-accepted cost for the same pipeline.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.variational_ccr4bp_torus as vt
import cyclerfinder.search.variational_crnbp_torus as vc
from cyclerfinder.genome.composed_moon_map import resonance_semimajor
from cyclerfinder.search.crnbp_torus_ghost_guard import radau_ghost_guard

# ---------------------------------------------------------------------------
# Cheap fixture: a tiny, deliberately UNCONVERGED n1=1/n2=2 toy torus. Only used
# to exercise the ghost-guard's own comparison/anchoring MECHANICS fast (< 1s) --
# its dynamical quality is irrelevant to what these tests check.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def toy_torus() -> vc.CRNBPTorusVariationalResult:
    mu = 0.008
    io = crnbp.CRNBPPerturber(mu=0.0, a=0.6, omega=1.0, theta0=np.pi)
    gan = crnbp.CRNBPPerturber(mu=0.0005, a=1.5, omega=-0.5, theta0=0.0)
    system = crnbp.CRNBPSystem(mu=mu, perturbers=(io, gan))
    n1, n2 = 1, 2
    rng = np.random.default_rng(0)
    coeffs0 = np.zeros((4, 2 * n1 + 1, 2 * n2 + 1))
    coeffs0[0, 0, 0] = 1.0 - mu + 0.1
    coeffs0[3, 0, 0] = 1.0
    coeffs0 += rng.normal(scale=1e-3, size=coeffs0.shape)
    return vc.correct_crnbp_torus_pseudospectral(
        system,
        coeffs0,
        omega2_0=0.6,
        n1=n1,
        n2=n2,
        amplitude_anchor=0.05,
        tol=1e-6,
        closure_tol=1e-4,
        max_nfev=100,
        gauge_weight=1.0,
        rho_weight=1.0,
        tr_solver="exact",
    )


def test_toy_torus_has_a_populated_closure_residual(
    toy_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    assert toy_torus.closure_residual is not None
    assert np.isfinite(toy_torus.closure_residual)
    assert toy_torus.closure_residual > 0.0


def test_dop853_side_reproduces_stored_closure_residual_exactly(
    toy_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """Same-anchor sanity gate: :func:`radau_ghost_guard`'s own DOP853 side must reproduce
    ``result.closure_residual`` to (near) machine precision -- it draws the SAME
    (theta1, theta2) anchors, at the SAME dt/seed, that
    ``variational_crnbp_torus._independent_closure`` used to compute that stored value in the
    first place. This is the module's own guard against the `#702` anchoring-mismatch failure
    class: if this test ever fails, the sampling here has silently diverged from the
    corrector's own anchors and nothing else in the report can be trusted (see
    ``CRNBPTorusGhostGuardReport.reproduces_reported_closure``)."""
    report = radau_ghost_guard(toy_torus)
    assert report.dop853_closure_residual == pytest.approx(toy_torus.closure_residual, rel=1e-9)
    assert report.reproduces_reported_closure


def test_radau_agrees_with_dop853_to_a_comfortable_margin_on_the_toy_torus(
    toy_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """Radau-vs-DOP853 agreement is an INTEGRATOR property of the short forward flow, not a
    property of how well-converged the torus is -- even this deliberately unconverged toy
    torus's forward flow should see both integrators agree to a tight margin (well under the
    `#701`/`#708` Umbriel-Titania row's own 1.0 km ``ghost_guard`` gate)."""
    report = radau_ghost_guard(toy_torus)
    assert report.radau_consistent
    assert report.integrator_delta_max_km < 1e-3, report.integrator_delta_max_km
    assert report.genuine


def test_mismatched_seed_anchor_is_detected_not_silently_trusted(
    toy_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """Regression for the `#702`-style anchoring discipline this module's docstring commits
    to, adapted to this object class: cross-checking against a DIFFERENT random sample set
    than the one ``result.closure_residual`` was computed from is not a legitimate
    apples-to-apples comparison -- it's two unrelated measurements. ``radau_ghost_guard`` must
    DETECT this (via ``reproduces_reported_closure``) rather than silently report a possibly
    spurious pass or fail on the Radau side, exactly as `#702`'s fixed ``ghost_guard`` refuses
    to compare two different manifold lobes without noticing."""
    mismatched = radau_ghost_guard(toy_torus, seed=12345)
    assert not mismatched.reproduces_reported_closure
    assert not mismatched.genuine
    assert "anchors" in mismatched.notes or "diverged" in mismatched.notes


def test_failed_reproduction_gate_short_circuits_the_genuine_verdict(
    toy_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """Even if the mismatched-anchor Radau side happened to numerically agree with the
    mismatched-anchor DOP853 side (both computed self-consistently against the WRONG target),
    ``genuine`` must still be False -- ``reproduces_reported_closure`` is a hard precondition,
    not merely one signal among several, per this module's own docstring."""
    mismatched = radau_ghost_guard(toy_torus, seed=999)
    if mismatched.radau_consistent:
        # Even where the wrong-anchor comparison happens to be internally consistent between
        # the two integrators, the overall verdict must not be genuine.
        assert not mismatched.genuine
    assert not mismatched.reproduces_reported_closure


# ---------------------------------------------------------------------------
# The real evidence-grade positive control: the ACTUAL physical-mass n1=2/n2=20 N=5 torus
# (`#720`/`#723`/`#724`'s object, the row this task's cross-check is FOR). Mirrors
# tests/search/test_variational_crnbp_torus.py's own construction exactly.
# ---------------------------------------------------------------------------


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical scaffolding to `#690`'s/`#720`'s own test modules."""
    from scipy.integrate import solve_ivp

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
def physical_n5_torus() -> vc.CRNBPTorusVariationalResult:
    """Reconstructs `#724`'s own delivered object bit-for-bit: Europa 3:4 symmetric resonant
    CR3BP orbit -> `#690` CCR4BP torus (n1=2/n2=20) -> mu_io=0 CRNBP seed -> 8-step mu_Io
    continuation to the physical mass (4.7043400305e-05), with `#724`'s own solver settings
    (``gauge_weight=30.0``, ``rho_weight=100.0`` throughout, including the continuation --
    the earlier, unpublished `#724` draft used the ``correct_crnbp_torus_pseudospectral``
    defaults of 1.0/1.0 for the continuation step by mistake and silently reproduced a
    DIFFERENT, if still-converged, torus; this fixture pins the correct weights explicitly so
    that regression cannot recur here)."""
    sys4 = ccr4bp.jupiter_europa_ganymede_default()
    target = crnbp.jupiter_europa_io_ganymede_default()
    io, _gan = target.perturbers

    s0, period, res = _resonant_symmetric_orbit(sys4.mu, 3, 4)
    assert res < 1e-10

    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys4,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert len(steps) == 8, "continuation diverged before reaching the physical mu_io"
    return steps[-1]


# #763: this test and the one below hit CI's global 600s pytest-timeout cap
# (run 30360261381) despite measuring only ~37-41s LOCALLY -- both are
# _LEVEL_EVIDENCE V1 evidence for europa-3-4-crnbp-torus-jupiter-2026 (this
# file is cited directly: "physical-torus positive control" / the ghost-
# guard gate itself). The ``physical_n5_torus`` module fixture rebuilds
# independently per xdist worker (each test can land on a different
# worker), so BOTH tests bear the full construction cost under -n auto, not
# just the first. A per-test @pytest.mark.timeout(1800) override (matching
# tests/search/test_variational_crnbp_torus.py's #763 fix, established
# precedent: tests/search/test_png_lane_recovery.py's
# @pytest.mark.timeout(1200)) gives ~44x margin over the local measurement
# without touching the pinned solver settings this fixture exists to
# regression-guard.
@pytest.mark.timeout(1800)
def test_physical_torus_reproduces_the_catalogued_provenance_numbers(
    physical_n5_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """Pins reproduction of the numbers already committed to ``data/catalogue.yaml``'s
    ``europa-3-4-crnbp-torus-jupiter-2026`` row (``crnbp_provenance.torus``) -- if this ever
    drifts by more than solver-path floating-point noise, the row's own provenance has gone
    stale. ``rel=1e-5`` (not exact-equality): a fresh ``least_squares(tr_solver="exact")``
    process, even with identical inputs, was observed to differ from the committed value by
    ~2.4e-8 relative (thread-scheduling-order floating-point noise in the dense linear
    algebra, not a genuine divergence -- both this run and the original are FAR inside every
    convergence gate)."""
    assert physical_n5_torus.residual_rms == pytest.approx(1.2343143649e-04, rel=1e-5)
    assert physical_n5_torus.closure_residual == pytest.approx(2.3212508597912776e-3, rel=1e-5)
    assert physical_n5_torus.rotation_number == pytest.approx(0.4964682691957192, rel=1e-5)


@pytest.mark.timeout(1800)  # #763: see the note above
def test_radau_ghost_guard_passes_the_physical_n5_torus_with_a_comfortable_margin(
    physical_n5_torus: vc.CRNBPTorusVariationalResult,
) -> None:
    """The `#738` headline result: the actual catalogued torus's short-time flow-consistency
    claim survives an independent Radau cross-check with the SAME anchors, DIFFERENT
    integrator, to a comfortable (not razor-thin) margin -- live-observed
    ``integrator_delta_max_km`` ~7e-8 km, six orders of magnitude below the `#701`/`#708`
    Umbriel-Titania row's own 1.0 km V1 gate. Loose bounds asserted (pinning behaviour, not
    an exact float, per this project's own convention)."""
    report = radau_ghost_guard(physical_n5_torus)
    assert report.reproduces_reported_closure, report.notes
    assert report.dop853_closure_residual == pytest.approx(
        physical_n5_torus.closure_residual, rel=1e-9
    )
    assert report.radau_consistent, report.notes
    assert report.integrator_delta_max_km < 1e-4, report.integrator_delta_max_km
    assert report.radau_closure_residual == pytest.approx(
        physical_n5_torus.closure_residual, rel=1e-6
    )
    assert report.genuine, report.notes
