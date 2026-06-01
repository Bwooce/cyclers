"""Multi-revolution Lambert: unit tests + lamberthub crosscheck gate.

The correctness gate compares the in-house multi-rev solver against
lamberthub.izzo2015 / gooding1990 with matching M (revs) and low_path
(branch) -- the same sourced-crosscheck discipline as the single-rev gate
in test_lambert.py. EXPECTED values trace to an independent third-party
solver, never to a value our own solver computed.

Plan: ``docs/superpowers/plans/2026-06-02-ml-multirev-lambert.md``.
"""

from __future__ import annotations

from math import acos, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import _min_time_of_revolution

from .conftest import Leg


@pytest.fixture(scope="session")
def leg_multirev() -> Leg:
    """Earth -> Mars, 780 d -- long enough that revolution n=1 is feasible.

    At 780 d the n=1 time-of-flight minimum (~448 d) is below the tof, so both
    the low and high branches exist (confirmed independently with lamberthub).
    The 500 d ``leg_long`` is too short for any full revolution.
    """
    eph = Ephemeris(model="circular")
    r1, _ = eph.state("E", 0.0)
    r2, _ = eph.state("M", 780.0 * SECONDS_PER_DAY)
    return Leg(name="multirev-EM-780d", r1=r1, r2=r2, tof=780.0 * SECONDS_PER_DAY)


def _a_coef_for(leg: Leg) -> tuple[float, float, float]:
    """Return ``(a_coef, r1_n, r2_n)`` for a leg (Vallado's A and the radii)."""
    r1_n = float(np.linalg.norm(leg.r1))
    r2_n = float(np.linalg.norm(leg.r2))
    cos_dnu = float(np.dot(leg.r1, leg.r2) / (r1_n * r2_n))
    cos_dnu = max(min(cos_dnu, 1.0), -1.0)
    dnu = acos(cos_dnu)
    cross_z = float(leg.r1[0] * leg.r2[1] - leg.r1[1] * leg.r2[0])
    if cross_z < 0.0:
        dnu = 2.0 * np.pi - dnu
    sin_dnu = float(np.sin(dnu))
    a_coef = sin_dnu * sqrt(r1_n * r2_n / (1.0 - cos_dnu))
    return a_coef, r1_n, r2_n


def test_min_time_of_revolution_is_a_lower_bound(leg_long: Leg) -> None:
    """t_min(n=1) is finite, positive, inside the n=1 domain.

    The endpoints z=(2*pi*n)^2 and z=(2*pi*(n+1))^2 are time singularities
    (t -> +inf); the interior minimum z_min must sit strictly between them.
    """
    a_coef, r1_n, r2_n = _a_coef_for(leg_long)
    z_min, t_min = _min_time_of_revolution(1, a_coef, r1_n, r2_n, MU_SUN_KM3_S2)
    assert (2.0 * np.pi) ** 2 < z_min < (4.0 * np.pi) ** 2
    assert t_min > 0.0
    assert np.isfinite(t_min)


def test_lambert_max_revs_returns_n0_plus_branches(leg_multirev: Leg) -> None:
    """A 780-day Earth->Mars leg admits n=0 and n=1 low/high.

    Asserts STRUCTURE only -- counts, n_revs labels, branch labels, shapes.
    Numerical correctness is the crosscheck gate (test below), whose
    EXPECTED side is lamberthub, not our own solver.
    """
    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_multirev.r1, leg_multirev.r2, leg_multirev.tof, max_revs=1)
    assert sols[0].n_revs == 0 and sols[0].branch == "single"
    n1 = [s for s in sols if s.n_revs == 1]
    branches = sorted(s.branch for s in n1)
    assert branches == ["high", "low"]  # both branches feasible for 780 d
    for s in sols:
        assert s.v1.shape == (3,) and s.v2.shape == (3,)
        assert s.v1.dtype == np.float64


def test_lambert_infeasible_revolution_is_skipped(leg_short: Leg) -> None:
    """A 50-day Earth->Earth arc is far too short for even n=1: only single-rev.

    t_min(1) for a 50-day short arc vastly exceeds 50 days, so revolution 1
    is infeasible and contributes zero solutions (skipped, not an error).
    """
    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_short.r1, leg_short.r2, leg_short.tof, max_revs=3)
    assert all(s.n_revs == 0 for s in sols)
    assert len(sols) == 1


@pytest.mark.parametrize("branch,low_path", [("low", True), ("high", False)])
def test_multirev_crosscheck_against_lamberthub(
    leg_multirev: Leg, branch: str, low_path: bool
) -> None:
    """In-house n=1 branch agrees with lamberthub izzo+gooding < 1e-3 m/s.

    EXPECTED values come from lamberthub (independent third-party solvers),
    not from our own solver -- satisfies the sourced-golden discipline.
    If this fails only by swapping branch<->low_path, the low/high labels in
    lambert() are inverted: swap them there, not here.
    """
    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    from cyclerfinder.core.lambert import lambert

    sols = lambert(leg_multirev.r1, leg_multirev.r2, leg_multirev.tof, max_revs=1)
    mine = next(s for s in sols if s.n_revs == 1 and s.branch == branch)

    r1 = np.asarray(leg_multirev.r1, dtype=np.float64)
    r2 = np.asarray(leg_multirev.r2, dtype=np.float64)
    v1_izzo, v2_izzo = izzo2015(
        MU_SUN_KM3_S2, r1, r2, leg_multirev.tof, M=1, prograde=True, low_path=low_path
    )
    v1_g, v2_g = gooding1990(
        MU_SUN_KM3_S2, r1, r2, leg_multirev.tof, M=1, prograde=True, low_path=low_path
    )
    worst_mps = 1000.0 * max(
        float(np.linalg.norm(mine.v1 - v1_izzo)),
        float(np.linalg.norm(mine.v2 - v2_izzo)),
        float(np.linalg.norm(mine.v1 - v1_g)),
        float(np.linalg.norm(mine.v2 - v2_g)),
    )
    assert worst_mps < 1.0e-3, (branch, worst_mps)


def test_lambert_crosscheck_multirev(leg_multirev: Leg) -> None:
    """lambert_crosscheck(..., n_revs=1, branch='low') agrees < 1e-3 m/s."""
    from cyclerfinder.core.lambert import lambert_crosscheck

    res = lambert_crosscheck(
        leg_multirev.r1, leg_multirev.r2, leg_multirev.tof, n_revs=1, branch="low"
    )
    assert res["max_diff_mps"] < 1.0e-3, res["max_diff_mps"]
