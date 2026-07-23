"""Tests for the CCR4BP manifold-tube intersection search + ghost-guard (`#694`).

Structure: coarse-search mechanics and guard rails first (small, cheap
grids), then the module's own positive-control finding -- a genuine,
ghost-guard-PASSED near-homoclinic connection on the Jupiter-Europa 3:4
CCR4BP torus (the same physical-mass torus `#690`/`#691` already validated),
refined from a specific, reproducible coarse-search seed. This is a
DETERMINISTIC computation (no randomness anywhere in the pipeline), so the
live-observed numbers below are pinned as regression bounds, not just
narrative -- see the `#694` commit message / final report for the full
numeric detail (this is the module-level distillation of that finding).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs
import cyclerfinder.search.ccr4bp_manifold_globalize as mg
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to `#690`/`#691`'s own test modules."""
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
def orbit_34(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    s0, period, res = _resonant_symmetric_orbit(system.mu, 3, 4)
    assert res < 1e-10
    return s0, period


@pytest.fixture(scope="module")
def phys_torus(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
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
def small_tubes(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> tuple[mg.ManifoldTube, mg.ManifoldTube]:
    """Small, cheap tube pair for coarse-search mechanics tests."""
    tube_u = mg.globalize_manifold_tube(
        phys_torus, "unstable", n_theta2=10, t_max_periods=1.0, n_time=20, n_segments_dir=16
    )
    tube_s = mg.globalize_manifold_tube(
        phys_torus, "stable", n_theta2=10, t_max_periods=1.0, n_time=20, n_segments_dir=16
    )
    return tube_u, tube_s


# ---------------------------------------------------------------------------
# coarse_candidates: input validation + the trivial-near-departure guard.
# ---------------------------------------------------------------------------


def test_coarse_candidates_rejects_wrong_branches(
    small_tubes: tuple[mg.ManifoldTube, mg.ManifoldTube],
) -> None:
    tube_u, tube_s = small_tubes
    with pytest.raises(ValueError, match="unstable"):
        hs.coarse_candidates(tube_s, tube_s)
    with pytest.raises(ValueError, match="stable"):
        hs.coarse_candidates(tube_u, tube_u)


def test_coarse_candidates_excludes_trivial_near_departure_match(
    phys_torus: vt.CCR4BPTorusVariationalResult,
    small_tubes: tuple[mg.ManifoldTube, mg.ManifoldTube],
) -> None:
    """Without the ``t_min_frac`` exclusion, the tightest "candidate" a naive
    nearest-neighbour search finds is always the trivial ``t_u = t_s = 0``
    match (both branches depart from nearly the same torus point) -- see the
    module docstring's account of this being FOUND during `#694`'s own
    development. Confirms it stays excluded for every returned candidate."""
    tube_u, tube_s = small_tubes
    t_min_frac = 0.15
    candidates = hs.coarse_candidates(tube_u, tube_s, n_candidates=8, t_min_frac=t_min_frac)
    assert len(candidates) > 0
    period = phys_torus.period
    for c in candidates:
        assert c.t_u >= t_min_frac * period - 1e-9
        assert c.t_s >= t_min_frac * period - 1e-9


def test_coarse_candidates_sorted_ascending_gap(
    small_tubes: tuple[mg.ManifoldTube, mg.ManifoldTube],
) -> None:
    tube_u, tube_s = small_tubes
    candidates = hs.coarse_candidates(tube_u, tube_s, n_candidates=6)
    gaps = [c.gap_planar for c in candidates]
    assert gaps == sorted(gaps)


def test_coarse_candidates_theta2_deduplication(
    small_tubes: tuple[mg.ManifoldTube, mg.ManifoldTube],
) -> None:
    tube_u, tube_s = small_tubes
    candidates = hs.coarse_candidates(tube_u, tube_s, n_candidates=5, min_theta2_separation=0.3)
    thetas = [c.theta2_u for c in candidates]
    for i in range(len(thetas)):
        for j in range(i + 1, len(thetas)):
            assert abs(thetas[i] - thetas[j]) >= 0.3 - 1e-9


# ---------------------------------------------------------------------------
# quasi_jacobi_gap.
# ---------------------------------------------------------------------------


def test_quasi_jacobi_gap_zero_for_identical_states(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    state = np.array([0.7, 0.6, 0.0, -0.3, 0.2, 0.0])
    gap = hs.quasi_jacobi_gap(state, state, phys_torus.system.mu)
    assert gap == pytest.approx(0.0, abs=1e-14)


# ---------------------------------------------------------------------------
# ghost_guard: the trivial near-departure "connection" must be REJECTED.
# ---------------------------------------------------------------------------


def test_ghost_guard_rejects_trivial_near_departure_pseudo_connection(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """A "refined connection" sitting right at the torus itself (t_u = t_s ~
    0, both branches essentially still AT their departure point) must fail
    the off-torus sanity check -- it is the trivial self-match, not a genuine
    excursion-and-return."""
    su = mg.manifold_state_at(phys_torus, "unstable", 0.0, 1.0, 1e-6, n_segments_dir=16)
    ss = mg.manifold_state_at(phys_torus, "stable", 0.0, 1.0, 1e-6, n_segments_dir=16)
    assert su is not None and ss is not None
    fake = hs.RefinedConnection(
        theta2_u=1.0,
        t_u=1e-6,
        theta2_s=1.0,
        t_s=1e-6,
        state_u=su,
        state_s=ss,
        pos_gap_km=float(np.linalg.norm(su[:3] - ss[:3])) * hs._L_KM,
        vel_gap_km_s=0.0,
        residual_norm=0.0,
        converged=True,
        n_iter=0,
    )
    guard = hs.ghost_guard(phys_torus, phys_torus, fake)
    assert not guard.genuine
    assert guard.off_torus_km < 1000.0
    assert "torus" in guard.notes


# ---------------------------------------------------------------------------
# The module's own positive control: a real, ghost-guard-PASSED near-
# homoclinic connection on the physical-mass Jupiter-Europa 3:4 torus.
# ---------------------------------------------------------------------------


def test_refine_candidate_reaches_tight_homoclinic_near_connection(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """DETERMINISTIC finding (reproduced across repeated calls during `#694`'s
    development -- see the module docstring / commit message): starting from
    the coarse-search seed ``(theta2_u=3.613, t_u=5.243, theta2_s=0.157,
    t_s=6.082)`` -- one of the closer candidates
    :func:`~cyclerfinder.search.ccr4bp_heteroclinic_search.coarse_candidates`
    finds on a ``n_theta2=40`` unstable/stable tube pair off the SAME torus
    -- continuous refinement converges to a position gap of ~37.8 km and a
    velocity gap of ~0.38 m/s: far tighter than any prior connection-search
    floor in this project's manifold lineage (`#619`'s ~166,000 km,
    `#620`'s ~351,000 km, `#626`'s ~1.3-1.9e6 km on the SE<->EM QBCP arc).
    Loose bounds (not exact-equality) are asserted -- pinning behaviour, not
    an exact float."""
    candidate = hs.ManifoldCandidate(
        theta2_u=3.613, t_u=5.243, theta2_s=0.157, t_s=6.082, gap_planar=0.00093
    )
    refined = hs.refine_candidate(phys_torus, phys_torus, candidate)
    assert refined is not None
    assert refined.pos_gap_km < 100.0, refined.pos_gap_km
    assert refined.vel_gap_km_s < 0.01, refined.vel_gap_km_s
    assert refined.converged


def test_ghost_guard_passes_the_positive_control_connection(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """The SAME refined connection above, run through the mandatory ghost
    guard: Radau-vs-DOP853 agreement to well under a metre, the candidate
    genuinely far from the torus (live-observed ~2000 km -- the EXACT
    divergence-from-the-unperturbed-torus-evolution metric, not a coarse-grid
    proxy; well above the 1000 km trivial-near-departure gate), and the
    guard's own ``genuine`` verdict True -- this connection SURVIVES the
    exact scrutiny that killed `#620`/`#626`'s "it closed!" ghost minima."""
    candidate = hs.ManifoldCandidate(
        theta2_u=3.613, t_u=5.243, theta2_s=0.157, t_s=6.082, gap_planar=0.00093
    )
    refined = hs.refine_candidate(phys_torus, phys_torus, candidate)
    assert refined is not None
    guard = hs.ghost_guard(phys_torus, phys_torus, refined)
    assert guard.radau_consistent, guard
    assert guard.integrator_delta_km < 1e-3, guard.integrator_delta_km
    assert guard.off_torus_km > 1_000.0, guard.off_torus_km
    assert abs(guard.quasi_jacobi_gap) < 1e-2, guard.quasi_jacobi_gap
    assert guard.genuine, guard.notes


def test_refine_candidate_seed_failure_returns_none(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """A seed so far off (both branches' extraction fails or the states are
    wildly separated) that the initial residual saturates the ``1e3`` sentinel
    must return ``None``, not a garbage "refined" result."""
    bad = hs.ManifoldCandidate(theta2_u=0.0, t_u=-5.0, theta2_s=0.0, t_s=-5.0, gap_planar=0.0)
    refined = hs.refine_candidate(phys_torus, phys_torus, bad)
    assert refined is None
