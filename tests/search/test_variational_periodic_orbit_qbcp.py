"""Seedless spectral (harmonic-balance) QBCP periodic-orbit discovery tests (#611).

Positive control (MANDATORY, load-bearing): cold-starts
:func:`cyclerfinder.search.variational_periodic_orbit_qbcp.discover_qbcp_periodic_orbit`
-- no ``warm_start``, no CR3BP-L1-to-BCR4BP-to-QBCP continuation bootstrap,
just the module's own default rough center guess near EM-L1 and the known
fixed period ``T_s`` -- targeting the published POL1 "dynamical substitute"
(Rosales & Jorba 2023, Table 4) in the EM-L1/L2 region #538/#544 root-caused
as "violently unstable" (every shooting-based QBCP corrector in this
codebase needs an expensive multi-stage bootstrap chain, or fails outright).

Numbers below are NOT copied from the module's docstring: they were
independently reproduced live in this session (2026-07-16) by running
:func:`discover_qbcp_periodic_orbit` directly and cross-checking against
the project's own independent 12-segment multiple-shooting corrector
(``scripts/analyze_593_qbcp_l1_substitute_reconciliation.py``).

Two bugs were found and fixed in the previous agent's implementation while
producing this test file (both in ``variational_periodic_orbit_qbcp.py``):

1. The shipped defaults (``n_harmonics=8``, ``n_restarts=8``,
   ``coefficient_noise=0.02``, ``tol=1e-6``) did NOT reliably converge at
   all (measured: 8/8 random restarts landed on the exact same
   ``residual_rms=9.4e-5`` floor, > tol, after 236s) -- ``n_harmonics=8`` is
   simply too coarse a Fourier truncation for this violently unstable
   region: it can satisfy the collocation-point residual while still being
   an O(1)-per-period-amplified, non-periodic loop (``closure_residual``
   as large as ``0.63``, i.e. NOT a real periodic orbit despite a passing
   residual). The default was raised to ``n_harmonics=32``, which
   empirically drives ``closure_residual`` down to ``~1e-6``, matching the
   independent multi-shooting corrector to near machine precision (see
   below). ``test_low_harmonics_converges_residual_but_not_closure`` below
   pins this exact failure mode as a regression.
2. ``scipy.optimize.least_squares(method="lm")`` (MINPACK's ``lmdif``) does
   NOT respect ``max_nfev`` as a real wall-clock bound on this problem:
   measured directly, requesting ``max_nfev=100`` cost 38,710 actual
   residual evaluations (~390x over). This is what caused the previous
   agent's session to appear to hang for hours: an unlucky random restart
   at ``n_harmonics=8`` (or worse, a would-be higher ``n_harmonics``) can
   run for a very long time under ``"lm"`` with no way to bound it via
   ``max_nfev``. Switched to ``method="trf"`` (Trust Region Reflective,
   pure Python), which still overshoots the requested ``max_nfev``
   somewhat (measured ~4x-23x) but stays bounded to tens of seconds per
   attempt rather than open-ended, and reduced the default ``max_nfev``
   from 30000 to 1500 accordingly.

Even after both fixes, this remains a STOCHASTIC method: most random cold
starts converge in ~9s, but an occasional unlucky seed can take several
minutes before landing on the (same) answer -- mirroring the CR3BP sibling
module's own documented "not every cold start converges" property. The
primary positive-control test below pins one specific, verified-fast seed
(0) to keep the test suite fast; it is not the ONLY seed that converges.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.search.variational_periodic_orbit_qbcp import (
    _eval_series,
    _n_free,
    _reconstruct_state0,
    _unpack,
    discover_qbcp_periodic_orbit,
)

# Published POL1 "dynamical substitute" golden (Rosales & Jorba 2023, Table 4),
# in this repository's reflected (x -> -x) frame -- identical constant to
# tests/core/test_qbcp.py's _POL1_REFLECTED (see that file for provenance).
_POL1 = np.array([0.8369141677649317, 0.0, 0.0, 0.0, 0.8391311559808445, 0.0])

# Independent cross-check golden: the project's OWN from-scratch 12-segment
# multiple-shooting corrector (scripts/analyze_593_qbcp_l1_substitute_reconciliation.py
# -- CR3BP-L1 (own quintic solve) -> BCR4BP via sequential mu_sun continuation
# -> QBCP handoff with real time-varying alphas, analytic STM, plain Newton),
# converging to periodicity resnorm 1.1292257521847358e-14. This is the FULL
# 6-state (the script itself only prints x0/py0); recomputed directly from
# its own build_l1_substitute(use_buggy_qbcp=False) internals in this session
# (2026-07-16), independently of and prior to running the harmonic-balance
# method below.
_MULTISHOOT_STATE0 = np.array(
    [
        0.8358664822188155,
        -0.01016652946503905,
        0.0,
        0.01022379318549257,
        0.8211145677631191,
        0.0,
    ]
)


def test_reconstruct_state0_matches_series_at_theta_zero() -> None:
    """``_reconstruct_state0`` must equal each component's Fourier series
    evaluated at theta=0 (cos(0)=1, sin(0)=0 -- so state0 = dc + sum(cos))."""
    n_harm = 3
    z = np.zeros(_n_free(n_harm))
    rng = np.random.default_rng(7)
    z[:] = rng.normal(scale=0.1, size=z.size)
    dc, cosc, sinc = _unpack(z, n_harm)
    state0 = _reconstruct_state0(z, n_harm)
    for v in range(6):
        theta0 = np.array([0.0])
        f0, _fp0 = _eval_series(theta0, dc[v], cosc[v], sinc[v])
        assert state0[v] == pytest.approx(float(f0[0]), abs=1e-14)


def test_eval_series_derivative_matches_known_sinusoid() -> None:
    """A single-harmonic series f(theta) = dc + a*cos(theta) + b*sin(theta)
    has the known analytic derivative fp(theta) = -a*sin(theta) + b*cos(theta);
    check ``_eval_series`` reproduces both at several points, not just theta=0."""
    dc, a, b = 0.5, 0.3, -0.2
    ccos = np.array([a])
    csin = np.array([b])
    theta = np.linspace(0.0, 2.0 * np.pi, 11, endpoint=False)
    f, fp = _eval_series(theta, dc, ccos, csin)
    expected_f = dc + a * np.cos(theta) + b * np.sin(theta)
    expected_fp = -a * np.sin(theta) + b * np.cos(theta)
    assert np.allclose(f, expected_f, atol=1e-14)
    assert np.allclose(fp, expected_fp, atol=1e-14)


def test_unpack_round_trips_raw_coefficients() -> None:
    """``_unpack`` must place each free-variable-vector entry into the
    documented (dc, cos, sin) row/column layout, not silently reorder it."""
    n_harm = 2
    n_free = _n_free(n_harm)
    z = np.arange(n_free, dtype=np.float64)
    dc, cosc, sinc = _unpack(z, n_harm)
    # Row order is (x, y, z, px, py, pz); each row consumes 1 + 2*n_harm
    # consecutive entries (dc, then cos[1..n], then sin[1..n]).
    block = 1 + 2 * n_harm
    for v in range(6):
        base = v * block
        assert dc[v] == base
        assert list(cosc[v]) == [base + 1, base + 2]
        assert list(sinc[v]) == [base + 3, base + 4]


def test_n_harmonics_and_period_multiple_validation() -> None:
    system = qbcp.qbcp_default()
    with pytest.raises(ValueError, match="n_harmonics"):
        discover_qbcp_periodic_orbit(system, n_harmonics=0)
    with pytest.raises(ValueError, match="period_multiple"):
        discover_qbcp_periodic_orbit(system, period_multiple=0)


def test_warm_start_shape_mismatch_raises() -> None:
    system = qbcp.qbcp_default()
    bad_warm_start = np.zeros(3)
    with pytest.raises(ValueError, match="warm_start"):
        discover_qbcp_periodic_orbit(system, n_harmonics=4, warm_start=bad_warm_start)


def test_low_harmonics_converges_residual_but_not_closure() -> None:
    """Documents WHY the default ``n_harmonics`` is 32, not the CR3BP
    sibling's 8 (module docstring's cited numbers, reproduced live here):
    at ``n_harmonics=8`` the harmonic-balance residual can plateau just
    above ``tol`` while ``closure_residual`` -- an independent check via
    real nonlinear propagation -- is O(1), i.e. NOT a periodic orbit at
    all. This is the harmonic-balance signature of #544's violent
    instability: too few Fourier degrees of freedom let the truncated
    series zero the residual AT collocation points while diverging wildly
    BETWEEN them once the true unstable flow amplifies the gap.
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(
        system,
        n_harmonics=8,
        n_restarts=1,
        coefficient_noise=0.0,
        rng=np.random.default_rng(0),
        tol=1e-6,
    )
    assert not res.converged
    assert res.residual_rms == pytest.approx(9.396109459538388e-05, rel=1e-3)
    assert res.closure_residual > 0.1  # NOT a genuine periodic orbit


def test_positive_control_cold_start_reproduces_qbcp_l1_substitute() -> None:
    """Seedless spectral method, cold-started with the module's own default
    center guess (no warm_start, no continuation bootstrap), converges to
    the EM-L1 QBCP periodic orbit anchoring the POL1 dynamical substitute --
    and does so to essentially machine precision agreement with this
    project's independently-built 12-segment multiple-shooting corrector.

    "Cold": only ``rng=np.random.default_rng(0)`` is supplied; the caller
    provides no state derived from ``_MULTISHOOT_STATE0`` or ``_POL1`` at
    all, only the module's own default ``center_guess=(0.85, ..., 0.8, ...)``
    (itself ~0.041 away from POL1's (x, py), not already-converged).
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(system, rng=np.random.default_rng(0))

    assert res.converged
    assert res.residual_rms < 1e-6
    # Independent check: propagating state0_pm through the TRUE nonlinear
    # QBCP EOM (not the truncated Fourier series) for the discovered fixed
    # period closes tightly -- not circular with residual_rms.
    assert res.closure_residual < 1e-4

    # Cross-check against the independent 12-segment multi-shooting
    # corrector: essentially machine-precision agreement (~5e-15 observed),
    # not merely "in the same neighborhood".
    diff_multishoot = np.linalg.norm(res.state0_pm - _MULTISHOOT_STATE0)
    assert diff_multishoot < 1e-8

    # Regression-pinned distances to the PUBLISHED POL1 golden. These are
    # NOT ~0 -- #544's own investigation attributes the gap to a
    # Gimeno-2018-vs-Rosales-2023 model-instance Fourier-refit difference,
    # not a corrector defect (see module docstring); both this method and
    # the independent multi-shooting corrector land at the SAME nonzero
    # distance, which is the actual claim under test here.
    xy_dist = math.hypot(res.state0_pm[0] - _POL1[0], res.state0_pm[4] - _POL1[4])
    full_dist = float(np.linalg.norm(res.state0_pm - _POL1))
    assert xy_dist == pytest.approx(0.018047024, abs=1e-6)
    assert full_dist == pytest.approx(0.023099337, abs=1e-6)

    # Second, fully independent confirmation: a different integrator
    # (Radau, not the module's own DOP853 closure check) over the
    # discovered period, from the discovered state.
    sol = solve_ivp(
        qbcp.qbcp_eom,
        (0.0, res.period),
        res.state0_pm,
        args=(system,),
        method="Radau",
        rtol=1e-12,
        atol=1e-12,
    )
    closure_radau = float(np.linalg.norm(sol.y[:, -1] - res.state0_pm))
    assert closure_radau < 1e-4


def test_planar_symmetry_components_are_near_zero() -> None:
    """The published POL1 substitute has y=z=px=pz=0 (planar, symmetric);
    this method's converged state has small but genuinely NONZERO y/px
    (~1e-2, matching the independent multi-shooting result exactly -- see
    the positive control above), while z/pz are near machine-zero. Pin
    that z/pz specifically stay at the expected near-zero floor (the
    planar collinear-point family), distinguishing "genuinely planar" from
    "small but real y/px offset from the published golden".
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(system, rng=np.random.default_rng(0))
    assert abs(res.state0_pm[2]) < 1e-10  # z0
    assert abs(res.state0_pm[5]) < 1e-10  # pz0
    assert abs(res.state0_pm[1]) > 1e-4  # y0 -- genuinely nonzero, not a rounding artifact
    assert abs(res.state0_pm[3]) > 1e-4  # px0 -- likewise
