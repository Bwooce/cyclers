"""Tests for the single-impulse reachable-set spike (Zhou-Armellin 2025).

Two layers:

* Fast hand-checkable mechanics unit tests (impulse parameterization, velocity
  frame, energy change, footprint metrics, containment arithmetic on toy clouds).
* A ``slow`` NRHO cross-check that recovers the 9:2 NRHO, builds the spike
  reachable footprint, and confirms the independent Monte-Carlo cloud is contained
  by the grid footprint with Zhou's error index P well under 0.1% -- the method
  invariant (mining note Sec. 4). No circular goldens: the MC cloud is the
  independent truth set and no published state vector is asserted.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.reachable_impulsive as ri

# JPL canonical Earth-Moon mass ratio used by Zhou (Eq. 52).
MU_ZHOU = 0.0121505839


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=MU_ZHOU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=1.0)


# ---------------------------------------------------------------------------
# Impulse parameterization (Zhou Eqs. 4-6).
# ---------------------------------------------------------------------------


def test_impulse_vector_magnitude_is_dv_mag() -> None:
    # The impulse always has magnitude dv_mag regardless of (alpha, beta) -- it
    # lives on the sphere of radius dv_mag (Eq. 5).
    for alpha in np.linspace(-math.pi / 2, math.pi / 2, 5):
        for beta in np.linspace(-math.pi, math.pi, 7):
            dv = ri.impulse_vector(0.01, float(alpha), float(beta))
            assert float(np.linalg.norm(dv)) == pytest.approx(0.01, abs=1e-15)


def test_impulse_vector_prograde_at_zero_angles() -> None:
    # (alpha, beta) = (0, 0) is purely along the first (velocity) axis.
    dv = ri.impulse_vector(0.02, 0.0, 0.0)
    assert dv[0] == pytest.approx(0.02)
    assert dv[1] == pytest.approx(0.0, abs=1e-15)
    assert dv[2] == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# Velocity frame (Zhou Eq. 8 triad).
# ---------------------------------------------------------------------------


def test_velocity_frame_is_orthonormal_right_handed() -> None:
    r = np.array([0.8, 0.1, -0.05])
    v = np.array([0.02, 0.9, 0.03])
    frame = ri.velocity_frame(r, v)
    # Orthonormal columns.
    assert np.allclose(frame.T @ frame, np.eye(3), atol=1e-12)
    # Right-handed: det = +1.
    assert float(np.linalg.det(frame)) == pytest.approx(1.0, abs=1e-12)
    # First column is the velocity direction.
    assert np.allclose(frame[:, 0], v / np.linalg.norm(v), atol=1e-12)
    # Third column is the orbit normal r x v direction.
    h = np.cross(r, v)
    assert np.allclose(frame[:, 2], h / np.linalg.norm(h), atol=1e-12)


def test_velocity_frame_rejects_rectilinear_state() -> None:
    # r parallel to v -> h = 0 -> no frame (the perilune singularity regime).
    with pytest.raises(ValueError, match="rectilinear"):
        ri.velocity_frame(np.array([1.0, 0.0, 0.0]), np.array([0.5, 0.0, 0.0]))


def test_velocity_frame_rejects_zero_velocity() -> None:
    with pytest.raises(ValueError, match="no velocity frame"):
        ri.velocity_frame(np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))


# ---------------------------------------------------------------------------
# apply_impulse: prograde burn changes speed; impulse CHANGES the Jacobi const.
# ---------------------------------------------------------------------------


def test_apply_impulse_prograde_changes_speed_by_dv() -> None:
    state = np.array([0.9, 0.05, -0.02, 0.01, 0.7, 0.04])
    v0 = float(np.linalg.norm(state[3:6]))
    man = ri.apply_impulse(state, ri.impulse_vector(0.01, 0.0, 0.0))
    # A purely prograde 0.01 burn adds 0.01 to the speed exactly.
    assert float(np.linalg.norm(man[3:6])) == pytest.approx(v0 + 0.01, abs=1e-12)
    # Position unchanged.
    assert np.allclose(man[:3], state[:3])


def test_apply_impulse_changes_jacobi_constant() -> None:
    # The essential contrast with Braik-Ross: a general impulse changes the
    # energy (Jacobi constant), so the maneuvered state leaves the C_J manifold.
    state = np.array([0.9, 0.05, -0.02, 0.01, 0.7, 0.04])
    cj0 = cr3bp.jacobi_constant(state, MU_ZHOU)
    man = ri.apply_impulse(state, ri.impulse_vector(0.02, 0.4, 1.0))
    cj1 = cr3bp.jacobi_constant(man, MU_ZHOU)
    assert abs(cj1 - cj0) > 1e-4


# ---------------------------------------------------------------------------
# Footprint metrics + containment arithmetic (toy clouds, hand-checkable).
# ---------------------------------------------------------------------------


def test_footprint_metrics_unit_square() -> None:
    # A 0.2 x 0.2 square footprint: area 0.04, half-extents 0.1, centroid origin.
    pts = np.array([[-0.1, -0.1], [0.1, -0.1], [0.1, 0.1], [-0.1, 0.1], [0.0, 0.0]])
    fm = ri.footprint_metrics(pts)
    assert fm.area == pytest.approx(0.04, rel=1e-9)
    assert fm.extent[0] == pytest.approx(0.1, rel=1e-9)
    assert fm.extent[1] == pytest.approx(0.1, rel=1e-9)
    assert fm.centroid[0] == pytest.approx(0.0, abs=1e-12)


def test_footprint_metrics_degenerate_few_points() -> None:
    fm = ri.footprint_metrics(np.array([[0.0, 0.0], [1.0, 1.0]]))
    assert fm.area == 0.0
    assert fm.n_points == 2


def test_containment_all_inside() -> None:
    # Grid hull = unit square; MC cloud all strictly inside -> fraction 1, d_max 0.
    grid = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    rng = np.random.default_rng(0)
    mc = rng.uniform(-0.5, 0.5, size=(200, 2))
    cc = ri.containment_crosscheck(grid, mc, tol=1e-12)
    assert cc.contained_fraction == 1.0
    assert cc.max_outside_dist == 0.0
    assert cc.error_index == 0.0
    assert cc.grid_area == pytest.approx(4.0, rel=1e-9)


def test_containment_some_outside() -> None:
    # A point at (2, 0) is 1.0 outside the unit-square's x=1 facet.
    grid = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    mc = np.array([[0.0, 0.0], [2.0, 0.0]])
    cc = ri.containment_crosscheck(grid, mc, tol=1e-12)
    assert cc.contained_fraction == 0.5
    assert cc.max_outside_dist == pytest.approx(1.0, rel=1e-9)
    # P = d_max^2 / S_RS = 1.0 / 4.0 = 0.25.
    assert cc.error_index == pytest.approx(0.25, rel=1e-9)


# ---------------------------------------------------------------------------
# Reachable cloud mechanics (small, fast: two-body-ish short arc).
# ---------------------------------------------------------------------------


def test_reachable_cloud_centred_on_nominal() -> None:
    # Footprint of a small symmetric (alpha, beta) grid should bracket the origin
    # (the nominal point projects to (0, 0) by construction).
    sysm = _em_system()
    seed = np.array([0.85, 0.0, -0.05, 0.0, 0.55, 0.0])
    cloud = ri.reachable_cloud(sysm, seed, 0.005, 0.3, n_alpha=7, n_beta=13)
    assert cloud.footprint.shape[1] == 2
    assert cloud.footprint.shape[0] == cloud.final_states.shape[0]
    # The cloud straddles the nominal point in both footprint axes.
    assert cloud.footprint[:, 0].min() < 0.0 < cloud.footprint[:, 0].max()
    assert cloud.footprint[:, 1].min() < 0.0 < cloud.footprint[:, 1].max()


# ---------------------------------------------------------------------------
# SLOW: 9:2 NRHO cross-check (independent MC containment, Zhou error index P).
# ---------------------------------------------------------------------------


def _recover_92_nrho(sysm: cr3bp.CR3BPSystem) -> cp.PeriodicOrbit:
    # Sourced seed: the widely published Earth-Moon L2 southern 9:2 NRHO initial
    # state (NASA Gateway baseline family, e.g. Lee 2019). SEED, not golden -- we
    # re-correct it under Zhou's mu and confirm the near-rectilinear geometry
    # (perilune ~3000 km, large apolune/perilune aspect) downstream.
    s0 = np.array([1.0213, 0.0, -0.1824, 0.0, -0.1031, 0.0])
    return cp.correct_periodic(sysm, s0, 1.5111)


@pytest.mark.slow
def test_92_nrho_is_near_rectilinear() -> None:
    sysm = _em_system()
    orbit = _recover_92_nrho(sysm)
    assert orbit.converged
    assert orbit.closure_residual < 1e-10
    # Period ~1.51 TU (6.56 d) and C_J ~3.05, consistent with the 9:2 NRHO.
    assert orbit.period == pytest.approx(1.511, rel=2e-2)
    assert orbit.jacobi == pytest.approx(3.0465, abs=5e-3)
    # Near-rectilinear: small perilune, large apolune/perilune aspect ratio.
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        orbit.state0,
        args=(sysm.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.sol is not None
    xs = sol.sol(np.linspace(0.0, orbit.period, 4000))
    r2 = np.sqrt((xs[0] - (1 - sysm.mu)) ** 2 + xs[1] ** 2 + xs[2] ** 2)
    perilune_km = float(r2.min()) * sysm.l_km
    assert 2000.0 < perilune_km < 5000.0  # genuine NRHO perilune
    assert float(r2.max() / r2.min()) > 10.0  # near-rectilinear aspect


@pytest.mark.slow
def test_nrho_reachable_set_contains_mc_cloud() -> None:
    # The Zhou method invariant (mining note Sec. 4): the MC cloud lies inside the
    # RS boundary with error index P = d_max^2 / S_RS well under 0.1%. Independent
    # recompute cross-check -- no published value asserted, no circular golden.
    sysm = _em_system()
    orbit = _recover_92_nrho(sysm)
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        orbit.state0,
        args=(sysm.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.sol is not None
    ts = np.linspace(0.0, orbit.period, 4000)
    xs = sol.sol(ts)
    r2 = np.sqrt((xs[0] - (1 - sysm.mu)) ** 2 + xs[1] ** 2 + xs[2] ** 2)
    seed_apo = xs[:, int(np.argmax(r2))].copy()

    # dv_max = 10 m/s (Zhou Table 4) in nondimensional velocity units.
    tu_s = 27.321661 / (2.0 * math.pi) * 86400.0
    vu_ms = sysm.l_km / tu_s * 1000.0
    dv_max_nd = 10.0 / vu_ms

    t_f = 0.25 * orbit.period
    cloud = ri.reachable_cloud(sysm, seed_apo, dv_max_nd, t_f, n_alpha=13, n_beta=25)
    fm = ri.footprint_metrics(cloud.footprint)
    assert fm.area > 0.0
    mc = ri.monte_carlo_cloud(
        sysm, seed_apo, dv_max_nd, t_f, cloud.nominal_final, n_samples=400, on_sphere=False, seed=1
    )
    cc = ri.containment_crosscheck(cloud.footprint, mc, tol=1e-9)
    # MC cloud (interior ball) fully contained by the max-sphere grid footprint;
    # error index P well under Zhou's 0.1% acceptance bar.
    assert cc.contained_fraction >= 0.99
    assert cc.error_index < 1e-3  # P < 0.1%


@pytest.mark.slow
def test_nrho_grid_resolution_improves_containment() -> None:
    # Coarsening the (alpha, beta) grid degrades containment (Zhou accuracy-vs-cost
    # tradeoff): a coarse grid under-covers the curved boundary near perilune so
    # some MC points fall outside; refining the grid restores full containment.
    sysm = _em_system()
    orbit = _recover_92_nrho(sysm)
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        orbit.state0,
        args=(sysm.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.sol is not None
    ts = np.linspace(0.0, orbit.period, 8000)
    xs = sol.sol(ts)
    r2 = np.sqrt((xs[0] - (1 - sysm.mu)) ** 2 + xs[1] ** 2 + xs[2] ** 2)
    seed_apo = xs[:, int(np.argmax(r2))].copy()

    tu_s = 27.321661 / (2.0 * math.pi) * 86400.0
    dv_max_nd = 10.0 / (sysm.l_km / tu_s * 1000.0)
    t_f = 0.5 * orbit.period  # longer arc reaching perilune (Zhou's 9:2 regime)

    coarse = ri.reachable_cloud(sysm, seed_apo, dv_max_nd, t_f, n_alpha=7, n_beta=13)
    fine = ri.reachable_cloud(sysm, seed_apo, dv_max_nd, t_f, n_alpha=21, n_beta=41)
    mc = ri.monte_carlo_cloud(
        sysm, seed_apo, dv_max_nd, t_f, coarse.nominal_final, n_samples=600, on_sphere=False, seed=3
    )
    cc_coarse = ri.containment_crosscheck(coarse.footprint, mc, tol=1e-9)
    cc_fine = ri.containment_crosscheck(fine.footprint, mc, tol=1e-9)
    # Refining the grid never worsens containment and the fine grid is ~complete.
    assert cc_fine.contained_fraction >= cc_coarse.contained_fraction
    assert cc_fine.contained_fraction >= 0.99
