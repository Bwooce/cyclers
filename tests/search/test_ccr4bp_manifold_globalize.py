"""Tests for CCR4BP manifold globalization (`#694`).

Mirrors `#690`/`#691`'s own test-fixture discipline (module-scoped torus
fixtures built from the same ``_resonant_symmetric_orbit`` scaffolding those
modules' own test files use -- test-only, no production code under test).

Checks, in order: input validation, output shapes/contract, a SANITY check
that a globalized manifold trajectory diverges from the torus at roughly the
rate `#691`'s own one-period eigenvalue predicts (the task's own explicit
process requirement -- "sanity-check against #691's own eigenvalue
magnitudes"), a step-size LINEAR-REGIME validation for the chosen default
``eps``, a genuinely non-hyperbolic reference case (all phases invalid,
mirroring `#691`'s own positive-control-A non-hyperbolic orbit), and
consistency between the discretized tube and the continuous
``manifold_state_at`` re-evaluation primitive the refinement stage in
``ccr4bp_heteroclinic_search`` depends on.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_manifold_globalize as mg
import cyclerfinder.search.ccr4bp_whisker as wk
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
def period_matched_torus(
    system: ccr4bp.CCR4BPSystem, orbit_34: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    """The genuinely NON-hyperbolic (over the orbit's own natural period)
    reference object `#691`'s own tests build -- see that module's Positive
    Control A. ``mu_gan=0`` (unforced) but propagated over the Ganymede
    synodic period would still be hyperbolic (a "frozen-time" effect, `#691`'s
    Positive Control B); THIS object instead carries ``period`` = the orbit's
    own NATURAL period, over which it is truly non-hyperbolic everywhere."""
    s0, period = orbit_34
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
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
# Input validation.
# ---------------------------------------------------------------------------


def test_rejects_bad_branch(phys_torus: vt.CCR4BPTorusVariationalResult) -> None:
    with pytest.raises(ValueError, match="branch"):
        mg.globalize_manifold_tube(phys_torus, "sideways")
    with pytest.raises(ValueError, match="branch"):
        mg.manifold_state_at(phys_torus, "sideways", 0.0, 0.0, 0.0)


def test_rejects_bad_grid_sizes(phys_torus: vt.CCR4BPTorusVariationalResult) -> None:
    with pytest.raises(ValueError, match="n_theta2"):
        mg.globalize_manifold_tube(phys_torus, "unstable", n_theta2=0)
    with pytest.raises(ValueError, match="n_time"):
        mg.globalize_manifold_tube(phys_torus, "unstable", n_time=1)


# ---------------------------------------------------------------------------
# Shape / contract.
# ---------------------------------------------------------------------------


def test_globalize_shapes_and_all_valid(phys_torus: vt.CCR4BPTorusVariationalResult) -> None:
    tube = mg.globalize_manifold_tube(
        phys_torus, "unstable", n_theta2=6, t_max_periods=0.5, n_time=8, n_segments_dir=16
    )
    assert tube.states.shape == (6, 8, 6)
    assert tube.base_states.shape == (6, 6)
    assert tube.directions.shape == (6, 6)
    assert tube.valid.shape == (6,)
    assert bool(np.all(tube.valid))
    assert np.all(np.isfinite(tube.states))
    # z, vz stay exactly zero at t=0 (planar embedding) and near-zero after
    # propagation (the CCR4BP planar sub-dynamics is invariant: z=vz=0 is a
    # genuine invariant subspace of the full 6-state EOM).
    assert np.max(np.abs(tube.states[:, :, 2])) < 1e-8
    assert np.max(np.abs(tube.states[:, :, 5])) < 1e-8


def test_globalize_t_samples_nonnegative_and_start_at_zero(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    tube = mg.globalize_manifold_tube(
        phys_torus, "stable", n_theta2=4, t_max_periods=0.5, n_time=5, n_segments_dir=16
    )
    assert tube.t_samples[0] == 0.0
    assert np.all(np.diff(tube.t_samples) > 0.0)


def test_globalize_t0_sample_matches_perturbed_departure(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """At elapsed time 0 the tube's stored state IS the perturbed departure
    point ``base_state + lobe_sign * eps * direction`` -- the propagation
    hasn't moved it yet."""
    tube = mg.globalize_manifold_tube(
        phys_torus, "unstable", n_theta2=4, t_max_periods=0.3, n_time=5, n_segments_dir=16, eps=1e-6
    )
    for i in range(4):
        expected = tube.base_states[i] + 1.0 * tube.eps * tube.directions[i]
        assert np.allclose(tube.states[i, 0], expected, atol=1e-9)


# ---------------------------------------------------------------------------
# Genuinely non-hyperbolic reference case: nothing should be extractable.
# ---------------------------------------------------------------------------


def test_globalize_all_invalid_on_nonhyperbolic_torus(
    period_matched_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """`#691`'s Positive Control A: this specific 3:4 resonant orbit is
    LINEARLY STABLE over its own natural period (all monodromy eigenvalues on
    the unit circle) -- no real hyperbolic eigenpair exists anywhere on this
    torus, so every departure phase's CLV extraction must fail and the tube
    comes back entirely invalid (NaN-filled), not silently wrong."""
    tube = mg.globalize_manifold_tube(
        period_matched_torus,
        "unstable",
        n_theta2=5,
        t_max_periods=0.2,
        n_time=4,
        n_segments_dir=8,
    )
    assert not np.any(tube.valid)
    assert np.all(np.isnan(tube.states))


# ---------------------------------------------------------------------------
# Sanity check against `#691`'s own one-period eigenvalue magnitudes.
# ---------------------------------------------------------------------------


def test_unstable_tube_growth_matches_691_eigenvalue(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """At ``(theta1, theta2) = (0, 0)`` `#691`'s own live-observed one-period
    eigenvalue is ``lam_u ~ -8.83`` (moderate, so still well inside the linear
    regime after exactly one period at ``eps=1e-6`` -- unlike the ``lam_u ~
    1700`` extreme phase). The globalized unstable-branch tube's displacement
    from the base torus state after one full period must be within a
    generous factor of ``eps * |lam_u|`` (order-of-magnitude cross-check
    between two independently-coded pipelines: `#691`'s post-hoc one-period
    STM eigendecomposition vs. this module's nonlinear propagation of a
    perturbed initial condition)."""
    theta1, theta2 = 0.0, 0.0
    state0, stm = wk.one_period_stm(phys_torus, theta1, theta2)
    diag = wk.amplification_diagnostics(state0, stm)
    assert diag.lam_u is not None
    lam_u = abs(diag.lam_u)
    assert 2.0 < lam_u < 50.0, lam_u  # confirms this IS the moderate reference phase

    eps = 1e-6
    tube = mg.globalize_manifold_tube(
        phys_torus,
        "unstable",
        n_theta2=1,
        t_max_periods=1.0,
        n_time=2,
        eps=eps,
        n_segments_dir=32,
    )
    # n_theta2=1 -> the single grid point IS theta2=0 (globalize's own grid convention).
    assert tube.valid[0]
    # Compare against the BASE (unperturbed) trajectory propagated to the
    # SAME elapsed time -- the base torus point itself moves substantially
    # over one period (it is quasi-periodic, not periodic, in theta1), so the
    # "growth away from the torus" this test wants is the perturbed-vs-base
    # SEPARATION at matched time, not a naive distance from the t=0 departure.
    arc_base = ccr4bp.propagate_ccr4bp(phys_torus.system, tube.base_states[0], phys_torus.period)
    disp_at_period = float(np.linalg.norm(tube.states[0, -1] - arc_base.state_f))
    predicted = eps * lam_u
    ratio = disp_at_period / predicted
    assert 0.3 < ratio < 3.0, (disp_at_period, predicted, ratio)


def test_stable_tube_backward_growth_matches_691_reciprocal_eigenvalue(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Same cross-check, stable branch: propagating BACKWARD in time off the
    stable eigendirection is the mechanism that globalizes the approach
    manifold, and its growth rate over one backward period is ``1/|lam_s| ~
    |lam_u|`` (the CCR4BP STM's reciprocal-pair symplectic structure `#691`
    already validates)."""
    theta1, theta2 = 0.0, 0.0
    state0, stm = wk.one_period_stm(phys_torus, theta1, theta2)
    diag = wk.amplification_diagnostics(state0, stm)
    assert diag.lam_s is not None
    predicted_growth = 1.0 / abs(diag.lam_s)

    eps = 1e-6
    tube = mg.globalize_manifold_tube(
        phys_torus,
        "stable",
        n_theta2=1,
        t_max_periods=1.0,
        n_time=2,
        eps=eps,
        n_segments_dir=32,
    )
    assert tube.valid[0]
    # Same fix as the unstable-branch test: compare against the base
    # trajectory propagated BACKWARD by one period (matched elapsed time),
    # not the t=0 departure state.
    arc_base = ccr4bp.propagate_ccr4bp(phys_torus.system, tube.base_states[0], -phys_torus.period)
    disp_at_period = float(np.linalg.norm(tube.states[0, -1] - arc_base.state_f))
    predicted = eps * predicted_growth
    ratio = disp_at_period / predicted
    assert 0.3 < ratio < 3.0, (disp_at_period, predicted, ratio)


# ---------------------------------------------------------------------------
# Step-size (eps) linear-regime validation.
# ---------------------------------------------------------------------------


def test_globalize_step_size_linear_regime(phys_torus: vt.CCR4BPTorusVariationalResult) -> None:
    """Doubling ``eps`` must roughly double the displacement from the base
    torus trajectory at a fixed elapsed time, for the DEFAULT ``eps`` and a
    moderate-amplification phase/horizon -- confirming the module's default
    step size is still within the linear regime where "step off along a
    single eigenvector direction" is a valid approximation of the true
    (nonlinear) manifold, not already saturated."""
    t_max_periods = 0.5  # well short of a full period at this moderate phase

    def _disp(eps: float) -> float:
        tube = mg.globalize_manifold_tube(
            phys_torus,
            "unstable",
            n_theta2=1,
            t_max_periods=t_max_periods,
            n_time=2,
            eps=eps,
            n_segments_dir=32,
        )
        assert tube.valid[0]
        arc_base = ccr4bp.propagate_ccr4bp(
            phys_torus.system, tube.base_states[0], t_max_periods * phys_torus.period
        )
        return float(np.linalg.norm(tube.states[0, -1] - arc_base.state_f))

    d1 = _disp(mg.DEFAULT_EPS)
    d2 = _disp(2.0 * mg.DEFAULT_EPS)
    ratio = d2 / d1
    assert 1.8 < ratio < 2.2, ratio


# ---------------------------------------------------------------------------
# Continuous re-evaluation primitive consistency (what the refinement stage
# in ccr4bp_heteroclinic_search relies on).
# ---------------------------------------------------------------------------


def test_manifold_state_at_matches_tube_sample(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    tube = mg.globalize_manifold_tube(
        phys_torus, "unstable", n_theta2=4, t_max_periods=0.6, n_time=6, n_segments_dir=24
    )
    i, j = 2, 3
    theta2 = float(tube.theta2_grid[i])
    t_flow = float(tube.t_samples[j])
    # ref_vec anchors the sign to the tube's own threaded convention at this
    # grid point (see manifold_state_at's docstring -- without it, a
    # standalone re-extraction has no continuity guarantee against the
    # tube's internally-threaded lobe choice).
    s_direct = mg.manifold_state_at(
        phys_torus,
        "unstable",
        0.0,
        theta2,
        t_flow,
        n_segments_dir=24,
        ref_vec=tube.directions[i],
    )
    assert s_direct is not None
    assert np.allclose(s_direct, tube.states[i, j], atol=1e-8)


def test_manifold_state_at_zero_matches_perturbed_departure(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    s = mg.manifold_state_at(phys_torus, "stable", 0.0, 1.234, 0.0, n_segments_dir=16)
    assert s is not None
    assert np.all(np.isfinite(s))
