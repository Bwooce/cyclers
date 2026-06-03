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
