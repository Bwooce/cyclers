"""Tests for :mod:`cyclerfinder.core.lambert`.

Includes the M1 gate: ``lambert_crosscheck(...)["max_diff_mps"] < 1e-3`` on
three distinct legs (Aldrin medium, short Earth-to-Earth, long Earth-to-Mars).

Plan: ``docs/phases/m1-core-mechanics/plan.md`` §4.4, §4.1.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.lambert import (
    LambertGeometryError,
    LambertSolution,
    lambert,
    lambert_crosscheck,
)

from .conftest import Leg

# ---------------------------------------------------------------------------
# Standalone behaviour
# ---------------------------------------------------------------------------


def test_lambert_returns_list_singleton(leg_aldrin: Leg) -> None:
    """Single-rev request returns a length-1 list with the expected metadata."""
    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof)
    assert isinstance(sols, list)
    assert len(sols) == 1
    sol = sols[0]
    assert sol.n_revs == 0
    assert sol.branch == "single"
    assert sol.v1.shape == (3,)
    assert sol.v2.shape == (3,)
    assert sol.v1.dtype == np.float64
    assert sol.v2.dtype == np.float64


def test_lambert_max_revs_too_short_returns_single_rev(leg_aldrin: Leg) -> None:
    """A 146 d Earth->Mars arc is below t_min(1), so no full revolution fits.

    ``max_revs=2`` is honoured but every revolution n>=1 is infeasible at this
    short tof, so only the single-revolution solution is returned.
    """
    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof, max_revs=2)
    assert len(sols) == 1
    assert sols[0].n_revs == 0


def test_lambert_retrograde(leg_aldrin: Leg) -> None:
    """Retrograde request matches lamberthub with the same flag set."""
    from lamberthub import izzo2015  # type: ignore[import-untyped]

    sols = lambert(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof, prograde=False)
    sol = sols[0]
    v1_ref, v2_ref = izzo2015(
        MU_SUN_KM3_S2,
        np.asarray(leg_aldrin.r1, dtype=np.float64),
        np.asarray(leg_aldrin.r2, dtype=np.float64),
        leg_aldrin.tof,
        M=0,
        prograde=False,
    )
    diff = max(
        float(np.linalg.norm(sol.v1 - v1_ref)),
        float(np.linalg.norm(sol.v2 - v2_ref)),
    )
    assert diff * 1000.0 < 1.0e-3


def test_lambert_zero_tof_raises(leg_aldrin: Leg) -> None:
    """Non-positive ``tof`` is a :class:`ValueError`."""
    with pytest.raises(ValueError):
        lambert(leg_aldrin.r1, leg_aldrin.r2, 0.0)
    with pytest.raises(ValueError):
        lambert(leg_aldrin.r1, leg_aldrin.r2, -100.0)


def test_lambert_180_deg_raises() -> None:
    """A pure 180-degree transfer raises :class:`LambertGeometryError`."""
    r1 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([-1.5 * AU_KM, 0.0, 0.0], dtype=np.float64)
    with pytest.raises(LambertGeometryError):
        lambert(r1, r2, 200.0 * SECONDS_PER_DAY)


def test_lambert_zero_magnitude_raises() -> None:
    """Zero-magnitude endpoints are a :class:`ValueError`."""
    r1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    with pytest.raises(ValueError):
        lambert(r1, r2, 100.0 * SECONDS_PER_DAY)


def test_lambert_solution_dataclass_frozen() -> None:
    """:class:`LambertSolution` is frozen — direct assignment is rejected."""
    sol = LambertSolution(
        n_revs=0,
        branch="single",
        v1=np.zeros(3, dtype=np.float64),
        v2=np.zeros(3, dtype=np.float64),
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        sol.n_revs = 1  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Gate: lamberthub cross-check on three legs
# ---------------------------------------------------------------------------


def test_aldrin_leg_cross_check(leg_aldrin: Leg) -> None:
    """Aldrin E->M ~146 d: agreement with izzo+gooding < 1e-3 m/s."""
    res = lambert_crosscheck(leg_aldrin.r1, leg_aldrin.r2, leg_aldrin.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]


def test_short_arc_cross_check(leg_short: Leg) -> None:
    """Earth->Earth short arc 50 d: agreement < 1e-3 m/s."""
    res = lambert_crosscheck(leg_short.r1, leg_short.r2, leg_short.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]


def test_long_arc_cross_check(leg_long: Leg) -> None:
    """Earth->Mars long arc 500 d: agreement < 1e-3 m/s."""
    res = lambert_crosscheck(leg_long.r1, leg_long.r2, leg_long.tof)
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]


# ---------------------------------------------------------------------------
# Bracket-finder robustness (task #56)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("tof_days", [5.0, 10.0])
def test_short_tof_high_energy_brackets_and_matches_lamberthub(tof_days: float) -> None:
    """A very-short, high-energy hyperbolic transfer is bracketed and solved.

    The valid universal-variable domain for a near-1 AU -> 1.52 AU, 0.8 rad
    transfer has its ``y(z) >= 0`` floor close to ``z = 0``; the feasible root
    for a 5-10 d time-of-flight sits in the narrow window between that floor and
    ``z = 0``. The prior fixed-start (``z_lo = -50``) widen walk halved toward
    ``z = 0`` from the invalid hyperbolic side and oscillated until the
    ``_BRACKET_MAX_WIDEN_ITERS`` cap, raising ``LambertConvergenceError`` on a
    transfer that physically has a solution. The floor-anchored bracket finder
    locates that window directly.

    EXPECTED values are sourced from ``lamberthub.izzo2015`` (an independent,
    published Lambert implementation) -- golden cross-check, not self-computed.
    """
    from lamberthub import izzo2015

    r2_n = 1.52 * AU_KM
    dnu = 0.8
    r1 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([r2_n * np.cos(dnu), r2_n * np.sin(dnu), 0.0], dtype=np.float64)
    tof = tof_days * SECONDS_PER_DAY

    sols = lambert(r1, r2, tof)
    assert len(sols) == 1
    sol = sols[0]

    v1_ref, v2_ref = izzo2015(MU_SUN_KM3_S2, r1, r2, tof, M=0, prograde=True)
    diff = max(
        float(np.linalg.norm(sol.v1 - v1_ref)),
        float(np.linalg.norm(sol.v2 - v2_ref)),
    )
    assert diff * 1000.0 < 1.0e-3, diff


def test_deep_floor_geometry_brackets_within_a_few_iters() -> None:
    """A deep-negative-floor geometry brackets quickly via floor bisection.

    For an Earth -> Jupiter-distance (1 -> 5.2 AU, ~2.5 rad) transfer the
    ``y(z) >= 0`` floor lies near ``z = -18``; a fixed-step linear widen walk
    from ``z_lo = -50`` would need many doublings/halvings to land inside the
    valid window. The bracket finder is instrumented to report its widen-loop
    iteration count via the private ``_bracket_diagnostics`` hook; assert it is
    well under the historical cap.
    """
    from cyclerfinder.core.lambert import _bracket_diagnostics

    r2_n = 5.2 * AU_KM
    dnu = 2.5
    r1 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    r2 = np.array([r2_n * np.cos(dnu), r2_n * np.sin(dnu), 0.0], dtype=np.float64)
    tof = 300.0 * SECONDS_PER_DAY

    sols = lambert(r1, r2, tof)
    assert len(sols) == 1

    diag = _bracket_diagnostics(r1, r2, tof)
    # Floor bisection converges in O(log2(range/tol)) ~ 60 steps worst case but
    # never spins at the cap; the prior linear walk could exhaust 100 here.
    assert diag["widen_iters"] < 80, diag


def test_bracket_diagnostics_degenerate_geometry_raises() -> None:
    """``_bracket_diagnostics`` guards singular geometry like ``lambert()`` does.

    At ``cos_dnu = 1`` (collinear, same direction) the ``a_coef`` expression
    divides by ``(1 - cos_dnu) = 0``; before the #200 fix the helper produced a
    NaN ``a_coef`` instead of raising :class:`LambertGeometryError`.
    """
    from cyclerfinder.core.lambert import _bracket_diagnostics

    r1 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    tof = 200.0 * SECONDS_PER_DAY
    # dnu = 0 (r2 parallel to r1, cos_dnu = 1).
    with pytest.raises(LambertGeometryError):
        _bracket_diagnostics(r1, 2.0 * r1, tof)
    # dnu = pi (anti-parallel) is singular too, same as lambert().
    with pytest.raises(LambertGeometryError):
        _bracket_diagnostics(r1, -1.5 * r1, tof)


# ---------------------------------------------------------------------------
# _dt_dz algebra fix (task #205 defect A)
# ---------------------------------------------------------------------------


def test_long_way_newton_converges_and_bvp_validates() -> None:
    """Pinning case for the ``_dt_dz`` term2 algebra error (task #205 defect A).

    Boundary states are built closed-form (independent of the Lambert solver)
    on a heliocentric ellipse at ``nu1`` and ``nu1 + dnu`` with ``dnu > pi``
    (long-way). The pre-fix ``_dt_dz`` carried a spurious ``sqrt(C)`` in
    term2 (``3*S*sqrt(y)/C^1.5`` instead of ``3*S*sqrt(y)/C``), making Newton
    oscillate at a contraction ratio ~0.873/step until the 60-iteration cap:
    this exact case raised :class:`LambertConvergenceError` (true root
    ``z ~ 10.349184``) even though the transfer is feasible. Post-fix it must
    converge and satisfy the boundary-value problem: propagating ``(r1, v1)``
    by ``tof`` with the in-house Kepler propagator reproduces ``r2`` to
    ``<= 1e-9`` relative.
    """
    from cyclerfinder.core.kepler import propagate

    from .conftest import coe3d_to_rv

    a_km = 1.991442e8
    e = 0.210298
    raan = 2.072523
    inc = 0.085380
    argp = 2.874090
    nu1 = 4.729555
    dnu = 3.624510  # > pi: long-way transfer
    tof = 2.1689e7  # s

    r1, _v1_gen = coe3d_to_rv(a_km, e, raan, inc, argp, nu1)
    r2, _v2_gen = coe3d_to_rv(a_km, e, raan, inc, argp, nu1 + dnu)

    # Pre-fix: raises LambertConvergenceError. Post-fix: converges.
    sols = lambert(r1, r2, tof)
    assert len(sols) == 1
    sol = sols[0]

    r2_prop, _v2_prop = propagate(r1, sol.v1, tof)
    rel = float(np.linalg.norm(r2_prop - r2) / np.linalg.norm(r2))
    assert rel <= 1.0e-9, rel


def test_dt_dz_matches_finite_difference() -> None:
    """Property test: analytic ``_dt_dz`` agrees with a central difference.

    Seeded random sweep over ``(z, A)`` (both ``A`` signs, elliptic and
    hyperbolic ``z``); the analytic derivative must match the central FD of
    ``_t_of_z`` to high relative accuracy everywhere ``y(z) > 0``. The pre-fix
    term2 algebra error produced FD ratios anywhere in ``-6.2 .. +4.7`` over
    this domain; the corrected form sits at FD truncation level (~1e-7).
    """
    from cyclerfinder.core.lambert import _dt_dz, _t_of_z

    rng = np.random.default_rng(20260611)
    n_checked = 0
    for _trial in range(400):
        r1_n = float(rng.uniform(0.5, 4.0)) * 1.5e8
        r2_n = float(rng.uniform(0.5, 4.0)) * 1.5e8
        dnu = float(rng.uniform(0.1, 2.0 * np.pi - 0.1))
        if abs(dnu - np.pi) < 0.05:
            continue
        cos_dnu = float(np.cos(dnu))
        a_coef = float(np.sin(dnu)) * float(np.sqrt(r1_n * r2_n / (1.0 - cos_dnu)))
        z = float(rng.uniform(-20.0, 35.0))
        if abs(z) < 1.0e-3:
            continue  # the near-zero Maclaurin window is exercised separately
        h = max(1.0e-7 * abs(z), 1.0e-8)
        try:
            _t, y = _t_of_z(z, a_coef, r1_n, r2_n, MU_SUN_KM3_S2)
            t_plus, y_plus = _t_of_z(z + h, a_coef, r1_n, r2_n, MU_SUN_KM3_S2)
            t_minus, y_minus = _t_of_z(z - h, a_coef, r1_n, r2_n, MU_SUN_KM3_S2)
        except ValueError:
            continue  # y < 0 somewhere in the stencil: outside the valid domain
        if y <= 0.0 or y_plus <= 0.0 or y_minus <= 0.0:
            continue
        fd = (t_plus - t_minus) / (2.0 * h)
        if fd == 0.0:
            continue
        analytic = _dt_dz(z, y, a_coef, MU_SUN_KM3_S2)
        assert abs(analytic / fd - 1.0) < 1.0e-4, (z, a_coef, analytic, fd)
        n_checked += 1
    assert n_checked > 100  # the sweep must actually exercise the domain
