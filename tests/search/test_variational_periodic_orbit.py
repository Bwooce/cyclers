"""Seedless spectral periodic-orbit discovery tests (#606).

Positive control (MANDATORY, load-bearing): reproduces the Earth-Moon L1
planar Lyapunov orbit -- an already-validated family in this codebase
(``cyclerfinder.search.cr3bp_seed_generator.lyapunov_seed``, itself sourced
from Koon, Lo, Marsden & Ross 2011 Ch. 2.5-2.7's collinear linearization) --
via :func:`cyclerfinder.search.variational_periodic_orbit.discover_periodic_orbit`
started COLD: a location offset from ``x_L1``, a period guess ~19% off the
true value, and random small higher-harmonic Fourier coefficients -- NOT
seeded from the shooting corrector's own converged state or eigenvector.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_periodic import correct_periodic, crosscheck_periodic
from cyclerfinder.search.cr3bp_seed_generator import lyapunov_seed
from cyclerfinder.search.variational_periodic_orbit import discover_periodic_orbit

# --- Ground truth reference, computed ONCE via independent, already-validated
# machinery (KLMR2011 linear-theory seed + this codebase's own STM corrector +
# an independent-integrator crosscheck) and pinned as a golden. Recomputed
# 2026-07-16: state0, period, jacobi below are exact printed values from that
# reference build (closure_residual 1.36e-14, crosscheck dj 1.78e-15).
_REF_STATE0 = np.array(
    [8.33992655e-01, 5.97485404e-09, 1.50765113e-39, 4.19326632e-09, 2.50150753e-02, 2.88911379e-31]
)
_REF_PERIOD = 2.693446543350158
_REF_JACOBI = 3.187810771959082

# The amplitude anchor for the new method's phase/scale gauge, extracted via an
# INDEPENDENT Fourier fit of the reference orbit's own x(t) (not a value the new
# solver computed): a0=0.83705516, a1 (cos, k=1) = -0.0030062993582869327. This
# only pins WHICH family member (amplitude) to look for -- not its location,
# period, or shape -- exactly analogous to how ``lyapunov_seed`` itself is
# parameterized by a caller-chosen amplitude ``Ax``.
_ANCHOR_X1 = -0.0030062993582869327


def _build_reference() -> None:
    """Sanity-check the pinned golden still reproduces from source (not just hardcoded)."""
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    seed_state0, seed_period = lyapunov_seed(sysm, point="L1", amplitude=3e-3)
    ref = correct_periodic(sysm, seed_state0, seed_period, tol=1e-13, max_iter=50)
    assert ref.converged
    ok, _dj = crosscheck_periodic(sysm, ref)
    assert ok
    assert np.linalg.norm(ref.state0 - _REF_STATE0) < 1e-6
    assert ref.period == pytest.approx(_REF_PERIOD, abs=1e-6)
    assert ref.jacobi == pytest.approx(_REF_JACOBI, abs=1e-6)


def test_reference_orbit_still_reproduces_from_source() -> None:
    _build_reference()


def test_positive_control_cold_start_reproduces_em_l1_lyapunov() -> None:
    """Seedless spectral method, cold-started, matches the validated EM L1
    Lyapunov orbit's state/period/Jacobi constant to a tight tolerance.

    "Cold": ``center_guess=(0.80, ...)`` is offset from the true x0=0.834 (and
    from x_L1=0.8369); ``period_guess=3.2`` is ~19% off the true period
    2.6934; the higher Fourier harmonics start at small random values (not
    the true shape). None of this is derived from ``_REF_STATE0`` -- only the
    scalar amplitude anchor (see module docstring) is shared, exactly the
    "how big an orbit do you want" input any periodic-orbit continuation
    needs to select a family member.
    """
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    rng = np.random.default_rng(123)
    res = discover_periodic_orbit(
        sysm,
        n_harmonics=8,
        anchor_x1=_ANCHOR_X1,
        center_guess=(0.80, 0.0, 0.0),
        period_guess=3.2,
        coefficient_noise=0.01,
        n_restarts=8,
        rng=rng,
        tol=1e-9,
    )
    assert res.converged
    assert res.residual_rms < 1e-9
    # Independent check: propagating state0 through the TRUE nonlinear EOM
    # (not the truncated Fourier series) for the discovered period closes
    # tightly -- this is not circular with residual_rms.
    assert res.closure_residual < 1e-8

    assert np.linalg.norm(res.state0 - _REF_STATE0) < 1e-6
    assert res.period == pytest.approx(_REF_PERIOD, abs=1e-9)
    assert res.jacobi == pytest.approx(_REF_JACOBI, abs=1e-9)

    # Second, fully independent confirmation: a different integrator (Radau)
    # over the discovered period, from the discovered state.
    sysm_mu = sysm.mu
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, res.period),
        res.state0,
        args=(sysm_mu,),
        method="Radau",
        rtol=1e-12,
        atol=1e-12,
    )
    closure_radau = float(np.linalg.norm(sol.y[:, -1] - res.state0))
    assert closure_radau < 1e-6


def test_cold_multistart_basin_is_broad_not_seed_specific() -> None:
    """The converged answer is reached from SEVERAL genuinely different cold
    starts (broad box for center/period), not just one lucky guess -- the
    core claim that this method's basin is not restricted to states already
    near the answer, unlike the existing shooting correctors.
    """
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    n_ok = 0
    n_try = 4
    for seed in range(n_try):
        rng = np.random.default_rng(1000 + seed)
        center_x = float(rng.uniform(0.55, 0.95))
        period_guess = float(rng.uniform(1.0, 6.0))
        res = discover_periodic_orbit(
            sysm,
            n_harmonics=8,
            anchor_x1=_ANCHOR_X1,
            center_guess=(center_x, 0.0, 0.0),
            period_guess=period_guess,
            coefficient_noise=0.01,
            n_restarts=3,
            rng=rng,
            tol=1e-8,
            # Bound wasted iteration on doomed random attempts -- a genuine
            # convergence in this problem takes O(10) nfev; capping keeps
            # this test fast without weakening the claim under test.
            max_nfev=2000,
        )
        if res.converged and np.linalg.norm(res.state0 - _REF_STATE0) < 1e-4:
            n_ok += 1
    # Not every random box position converges (this is a local method, not a
    # global one); a majority landing on the SAME known orbit from broadly
    # scattered cold starts is the claim under test.
    assert n_ok >= n_try // 2


def test_missing_amplitude_anchor_raises() -> None:
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError, match="amplitude anchor"):
        discover_periodic_orbit(sysm, n_harmonics=4)


def test_planar_result_has_zero_z_component() -> None:
    """With only anchor_x1 set (no z anchor), the discovered orbit is planar."""
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    rng = np.random.default_rng(5)
    res = discover_periodic_orbit(
        sysm,
        n_harmonics=8,
        anchor_x1=_ANCHOR_X1,
        center_guess=(0.80, 0.0, 0.0),
        period_guess=3.2,
        n_restarts=4,
        rng=rng,
        tol=1e-8,
    )
    assert abs(res.state0[2]) < 1e-8  # z0
    assert abs(res.state0[5]) < 1e-8  # zdot0


def test_pilot_warm_started_continuation_crosses_the_wall() -> None:
    """Bounded regression of ``scripts/run_606_variational_pilot.py``'s headline
    result: warm-started (from ONE existing-tool success point, Fourier-fit
    into this module's own coefficients) continuation using ONLY
    :func:`discover_periodic_orbit`'s own solves crosses several steps into
    the #556 near-bifurcation wall region with near-machine-precision
    closure. Bounded to 5 steps (vs. the pilot script's 20) to keep this
    test fast; the pilot script itself is the full-range record.
    """
    from cyclerfinder.search.cr3bp_seed_generator import richardson_halo_seed

    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    mu = sysm.mu
    n_harm = 16

    def fourier_fit_all(
        state0: NDArray[np.float64], period: float, n: int = n_harm, n_samples: int = 512
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period),
            state0,
            args=(mu,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
            dense_output=True,
        )
        assert sol.success
        # dense_output=True above (+ the success check just above) guarantees
        # a dense interpolant -- see the same pattern's comment in
        # src/cyclerfinder/search/perimoon_passage.py.
        assert sol.sol is not None
        tgrid = np.linspace(0.0, period, n_samples, endpoint=False)
        traj = sol.sol(tgrid)
        theta = 2.0 * np.pi * tgrid / period
        k = np.arange(1, n + 1)
        ang = np.outer(theta, k)
        cosm, sinm = np.cos(ang), np.sin(ang)
        dc = traj[:3].mean(axis=1)
        cosc = np.zeros((3, n))
        sinc = np.zeros((3, n))
        for c in range(3):
            f = traj[c] - dc[c]
            cosc[c] = 2.0 / n_samples * (cosm.T @ f)
            sinc[c] = 2.0 / n_samples * (sinm.T @ f)
        return dc, cosc, sinc

    def pack_z1_anchor(
        dc: NDArray[np.float64],
        cosc: NDArray[np.float64],
        sinc: NDArray[np.float64],
        period: float,
    ) -> NDArray[np.float64]:
        parts = [dc[0], dc[1], dc[2], cosc[0, 0]]
        parts.extend(cosc[0, 1:])
        parts.extend(sinc[0, 1:])
        parts.append(sinc[1, 0])
        parts.append(cosc[1, 0])
        parts.extend(cosc[1, 1:])
        parts.extend(sinc[1, 1:])
        parts.append(sinc[2, 0])
        parts.extend(cosc[2, 1:])
        parts.extend(sinc[2, 1:])
        parts.append(np.log(period))
        return np.array(parts, dtype=np.float64)

    ref_state0, ref_period = richardson_halo_seed(sysm, point="L1", amplitude_z=-0.038, branch="I")
    dc, cosc, sinc = fourier_fit_all(ref_state0, ref_period)
    az1 = float(cosc[2, 0])
    warm = pack_z1_anchor(dc, cosc, sinc, ref_period)
    center = (float(ref_state0[0]), float(ref_state0[1]), float(ref_state0[2]))
    period_guess = ref_period

    step = 0.0025
    for i in range(5):
        az1 += step
        res = discover_periodic_orbit(
            sysm,
            n_harmonics=n_harm,
            anchor_z1=az1,
            center_guess=center,
            period_guess=period_guess,
            n_restarts=1,
            tol=1e-9,
            max_nfev=10000,
            warm_start=warm,
            rng=np.random.default_rng(i),
        )
        assert res.converged
        assert res.closure_residual < 1e-6
        warm = res.raw_coeffs
        center = (float(res.state0[0]), float(res.state0[1]), float(res.state0[2]))
        period_guess = res.period

    # After 5 steps of 0.0025 the Jacobi constant should have moved measurably
    # TOWARD the confirmed bifurcation (3.1745, #555) from the start point
    # (3.145298), while staying strictly inside the documented wall region.
    assert res.jacobi > 3.145298 + 0.005
    assert res.jacobi < 3.1745


def test_pilot_l1_halo_wall_region_richardson_seed_fails() -> None:
    """Documents the #556 wall this pilot targets: the existing
    Richardson-analytic-seed + shooting corrector
    (:func:`cyclerfinder.search.cr3bp_seed_generator.richardson_halo_seed`)
    fails to converge to a genuine non-planar L1 halo across essentially the
    entire near-bifurcation amplitude range (independently reproduced here,
    not merely asserted from `data/OUTSTANDING.md`'s #548/#555 narrative).
    """
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    from cyclerfinder.search.cr3bp_seed_generator import richardson_halo_seed

    n_fail = 0
    n_total = 0
    for az in np.arange(0.015, 0.036, 0.004):
        n_total += 1
        try:
            richardson_halo_seed(sysm, point="L1", amplitude_z=-float(az), branch="I")
        except ValueError:
            n_fail += 1
    # The wall is essentially total in this range (observed 5/5 in scratchpad
    # exploration at finer spacing); require a strong majority so this test
    # is robust to the exact grid without being vacuous.
    assert n_fail >= n_total - 1
