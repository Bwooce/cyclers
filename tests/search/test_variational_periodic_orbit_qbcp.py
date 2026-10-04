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

# #892 (2026-10-04): this file used to cross-check against a 12-segment multiple-shooting
# state computed in the defective model (``_MULTISHOOT_STATE0``, 1.8e-2 from POL1, with
# y = -0.0102 and px = +0.0102). In the corrected model the method reproduces the published
# point itself, so the positive control now compares with ``_POL1`` directly.


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
    sibling's 8: at a low harmonic count the harmonic-balance residual is
    driven below ``tol`` (``converged=True``) while ``closure_residual`` --
    an independent check via real nonlinear propagation -- is O(1), i.e.
    NOT a periodic orbit at all. This is the harmonic-balance signature of
    #544's violent instability: too few Fourier degrees of freedom let the
    truncated series zero the residual AT collocation points while
    diverging wildly BETWEEN them once the true unstable flow amplifies the
    gap.

    #892 (2026-10-04): the test used ``n_harmonics=8`` and pinned a residual
    plateau of 9.396e-5 above ``tol``; both belonged to the defective model.
    In the corrected model (measured): ``n_harmonics`` 4 gives residual
    1.428e-7 with closure 0.26, 6 gives 1.0e-9 with closure 6.3e-2, 8 gives
    7.9e-12 with closure 5.9e-4. 4 is used so the closure gate still has an
    O(1) failure to catch.
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(
        system,
        n_harmonics=4,
        n_restarts=1,
        coefficient_noise=0.0,
        rng=np.random.default_rng(0),
        tol=1e-6,
    )
    assert res.converged  # the residual criterion alone is satisfied ...
    assert res.residual_rms == pytest.approx(1.427801e-07, rel=1e-3)
    assert res.closure_residual > 0.1  # ... but it is NOT a genuine periodic orbit


def test_positive_control_cold_start_reproduces_qbcp_l1_substitute() -> None:
    """Seedless spectral method, cold-started with the module's own default
    center guess (no warm_start, no continuation bootstrap), converges to
    the EM-L1 QBCP periodic orbit anchoring the POL1 dynamical substitute --
    and lands on the published POL1 point (#892: before the model was
    corrected it agreed only with a multiple-shooting state computed in the
    same defective model, 1.8e-2 from POL1).

    "Cold": only ``rng=np.random.default_rng(0)`` is supplied; the caller
    provides no state derived from ``_POL1`` at all, only the module's own
    default ``center_guess=(0.85, ..., 0.8, ...)`` (itself ~0.041 away from
    POL1's (x, py), not already-converged).
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(system, rng=np.random.default_rng(0))

    assert res.converged
    assert res.residual_rms < 1e-6
    # Independent check: propagating state0_pm through the TRUE nonlinear
    # QBCP EOM (not the truncated Fourier series) for the discovered fixed
    # period closes tightly -- not circular with residual_rms.
    assert res.closure_residual < 1e-4

    # Published positive control (#892): the published POL1 point (Rosales &
    # Jorba 2023 Table 4), reproduced to 1.7e-8 (measured 2026-10-04); the
    # core module's own multiple-shooting check in tests/core/test_qbcp.py
    # holds it to 1e-7 in x and py, the same bound used here.
    assert np.max(np.abs(res.state0_pm - _POL1)) < 1e-7

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
    """The published POL1 substitute has y=z=px=pz=0 (planar, and symmetric
    under the model's reversing symmetry (x, -y, z, -px, py, -pz, -t)).

    #892 (2026-10-04): this test used to assert that y and px were NONZERO
    (about 1e-2). That asymmetry was a symptom of the defective model, whose
    reversing symmetry failed; in the corrected model the symmetry holds
    (tests/core/test_qbcp.py) and the converged state has y, px of order
    1e-17 (measured), so they are now required to vanish.
    """
    system = qbcp.qbcp_default()
    res = discover_qbcp_periodic_orbit(system, rng=np.random.default_rng(0))
    assert abs(res.state0_pm[2]) < 1e-10  # z0
    assert abs(res.state0_pm[5]) < 1e-10  # pz0
    assert abs(res.state0_pm[1]) < 1e-10  # y0
    assert abs(res.state0_pm[3]) < 1e-10  # px0
