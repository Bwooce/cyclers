"""Tests for the free-Jacobi continuation residual/Jacobian (#523 rework).

POSITIVE CONTROL: the analytic dR/dC column must match a finite-difference
estimate at the SAME sourced positive-control seed
(``run_523_earth_coorbital_search.py``'s 2006-RH120-derived orbit) already
used elsewhere in this codebase -- an internal-consistency check, not a
value this module invented and then asserted against itself; the finite
difference is computed independently of the analytic formula being tested.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic
from cyclerfinder.search.cr3bp_general_periodic_free_c import (
    make_residual_and_jacobian_fns,
    residual_and_jacobian_free_c,
)
from cyclerfinder.search.pseudo_arclength import ContinuationStopReason, continue_curve

# Same #523 positive-control seed (2006 RH120-derived), certified once here
# via the EXISTING corrector to get a genuine on-family point to test at.
SEED_C = 2.9998797409719242


@pytest.fixture(scope="module")
def certified_seed() -> tuple[cr3bp.CR3BPSystem, float, float, float]:
    """Certify the #523 positive-control seed ONCE and share it across
    tests -- certification itself costs real wall-clock time (a long-horizon
    STM propagation per Newton iteration), so re-certifying per test would
    triple this file's runtime for no additional coverage.
    """
    system = cr3bp.cr3bp_system("Sun", "Earth")
    # A nearby coarse guess known (this session) to certify at SEED_C.
    orbit = correct_general_periodic(
        system,
        0.93000,
        0.03600,
        SEED_C,
        period_guess=100.0,
        half_crossings=2,
        ydot0_sign=1.0,
        tol=1e-11,
        t_hi_frac=1.15,
        max_iter=60,
    )
    assert orbit.converged and orbit.residual < 1e-9
    return system, orbit.x0, orbit.xdot0, orbit.period


def test_analytic_dc_column_matches_finite_difference(
    certified_seed: tuple[cr3bp.CR3BPSystem, float, float, float],
) -> None:
    system, x0, xdot0, period = certified_seed
    mu = float(system.mu)
    t_hi = 1.15 * period

    out = residual_and_jacobian_free_c(x0, xdot0, SEED_C, mu, 1.0, 2, t_hi)
    assert out is not None
    _r0, jac = out

    h = 1e-6
    out_plus = residual_and_jacobian_free_c(x0, xdot0, SEED_C + h, mu, 1.0, 2, t_hi)
    out_minus = residual_and_jacobian_free_c(x0, xdot0, SEED_C - h, mu, 1.0, 2, t_hi)
    assert out_plus is not None
    assert out_minus is not None
    fd_dr_dc = (out_plus[0] - out_minus[0]) / (2.0 * h)

    assert np.allclose(jac[:, 2], fd_dr_dc, atol=1e-4), (
        f"analytic dR/dC {jac[:, 2]} vs finite-difference {fd_dr_dc}"
    )


def test_analytic_dx0_dxdot0_columns_match_finite_difference(
    certified_seed: tuple[cr3bp.CR3BPSystem, float, float, float],
) -> None:
    """Sanity: the (x0, xdot0) columns should reproduce the SAME analytic
    result ``_residual_jac_stm`` already produces and this project already
    trusts (cross-check against finite differences here, independently).
    """
    system, x0, xdot0, period = certified_seed
    mu = float(system.mu)
    t_hi = 1.15 * period

    out = residual_and_jacobian_free_c(x0, xdot0, SEED_C, mu, 1.0, 2, t_hi)
    assert out is not None
    _r0, jac = out

    h = 1e-6
    for col, (dx0, dxdot0) in enumerate([(h, 0.0), (0.0, h)]):
        out_plus = residual_and_jacobian_free_c(x0 + dx0, xdot0 + dxdot0, SEED_C, mu, 1.0, 2, t_hi)
        out_minus = residual_and_jacobian_free_c(x0 - dx0, xdot0 - dxdot0, SEED_C, mu, 1.0, 2, t_hi)
        assert out_plus is not None
        assert out_minus is not None
        fd_col = (out_plus[0] - out_minus[0]) / (2.0 * h)
        assert np.allclose(jac[:, col], fd_col, atol=1e-4), (
            f"column {col}: analytic {jac[:, col]} vs finite-difference {fd_col}"
        )


def test_continuation_walks_the_family_with_shrinking_residual(
    certified_seed: tuple[cr3bp.CR3BPSystem, float, float, float],
) -> None:
    """Mechanical smoke test: continuation from the certified seed produces
    a sequence of well-converged points with C actually advancing.
    """
    system, x0, xdot0, period = certified_seed
    mu = float(system.mu)

    residual_fn, jacobian_fn = make_residual_and_jacobian_fns(
        mu, 1.0, 2, 1.15, period, rtol=1e-12, atol=1e-12
    )
    z0 = np.array([x0, xdot0, SEED_C])
    curve = continue_curve(
        residual_fn,
        z0,
        jacobian_fn=jacobian_fn,
        step_size=0.0005,
        max_steps=6,
        tol=1e-9,
        max_iter=30,
        step_caps=[0.01, 0.01, 0.001],
    )
    assert curve.stop_reason == ContinuationStopReason.MAX_STEPS
    assert len(curve.points) == 7
    for p in curve.points:
        assert p.residual_norm < 1e-8
    c_values = [p.z[2] for p in curve.points]
    assert abs(c_values[-1] - c_values[0]) > 1e-5, "C did not actually advance along the family"
