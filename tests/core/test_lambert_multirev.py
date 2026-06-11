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


# ---------------------------------------------------------------------------
# Illinois stall against the rev-boundary singularity (task #205 defect B)
# ---------------------------------------------------------------------------


def test_high_branch_near_rev_boundary_is_returned() -> None:
    """Pinning case for the silently-dropped n=1 high branch (task #205 B).

    Boundary states and the time of flight are built closed-form on a
    generating ellipse: ``r1, r2`` at ``nu1`` and ``nu1 + dnu``, ``tof`` the
    Kepler sweep time from ``nu1`` to ``nu1 + dnu`` plus one full period
    (~7.5 yr) — so the generating orbit itself *is* an exact ``n=1`` solution,
    with universal variable ``z ~ 138.7573``, close to the ``z = (4*pi)^2``
    revolution boundary where ``t(z)`` blows up (~1e29 s at the bracket
    endpoint). Pre-fix, Illinois iterated the raw residual ``t - tof``: false
    position against the ~1e29 s endpoint stalls (the Illinois halving needs
    ~70 iterations to deflate it, past the 60-iteration cap), the iterate sat
    AT the root with a 3.66e-4 s residual (1.5e-12 relative) but missed the
    1e-12 threshold, and ``lambert()`` caught the raise and silently omitted
    the feasible branch. Post-fix (log-compressed residual) the high branch
    must be returned, BVP-validate to ``<= 1e-9`` relative, and match the
    generating orbit's velocity (the closed-form expected value, independent
    of the Lambert solver).
    """
    from math import atan2, cos, pi, sin, sqrt

    from cyclerfinder.core.kepler import propagate
    from cyclerfinder.core.lambert import lambert

    from .conftest import coe3d_to_rv

    a_km = 373755302.85
    e = 0.279159
    raan = 2.810967
    inc = 0.405913
    argp = 2.335005
    nu1 = 5.955822
    dnu = 5.314586  # > pi: long-way transfer
    nrevs = 1

    def _time_since_periapsis(nu: float, ecc: float, n_mean: float) -> float:
        ecc_anom = 2.0 * atan2(sqrt(1.0 - ecc) * sin(nu / 2.0), sqrt(1.0 + ecc) * cos(nu / 2.0))
        mean_anom = ecc_anom - ecc * sin(ecc_anom)
        return mean_anom / n_mean

    n_mean = sqrt(MU_SUN_KM3_S2 / a_km**3)
    period = 2.0 * pi / n_mean
    nu2 = (nu1 + dnu) % (2.0 * pi)
    sweep = _time_since_periapsis(nu2, e, n_mean) - _time_since_periapsis(nu1, e, n_mean)
    if sweep < 0.0:
        sweep += period
    tof = sweep + nrevs * period  # ~7.5 yr

    r1, v1_gen = coe3d_to_rv(a_km, e, raan, inc, argp, nu1)
    r2, _v2_gen = coe3d_to_rv(a_km, e, raan, inc, argp, nu1 + dnu)

    sols = lambert(r1, r2, tof, max_revs=nrevs)
    high = [s for s in sols if s.n_revs == 1 and s.branch == "high"]
    assert high, [(s.n_revs, s.branch) for s in sols]  # pre-fix: branch silently dropped
    sol = high[0]

    # BVP validity against the in-house propagator.
    r2_prop, _v2_prop = propagate(r1, sol.v1, tof)
    rel = float(np.linalg.norm(r2_prop - r2) / np.linalg.norm(r2))
    assert rel <= 1.0e-9, rel

    # The generating orbit is the expected solution (closed-form source).
    assert 1000.0 * float(np.linalg.norm(sol.v1 - v1_gen)) < 1.0e-3
