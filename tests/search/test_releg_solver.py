"""Tests for the releg leg-solver seam (#449, DSM branch).

The releg genome replaces the per-leg ballistic ``lambert(...)`` call of the
moon-tour closure with a swappable leg solver. The contract is the ``Releg``
protocol returning a :class:`RelegResult` (vinf_out / vinf_in / dv_kms /
feasible). ``BallisticReleg`` is the zero-ΔV backend that must reproduce today's
ballistic Lambert leg bit-for-bit (the regression lock); ``DsmReleg`` is the
powered backend whose delivered ΔV is golden-anchored to the Campagnola-Russell
VILM leveraging floor (``search.vilm``), which is already golden-validated.

Source discipline (``feedback_golden_tests_sourced_only``): the powered ΔV
golden's EXPECTED side is the published VILM floor (``vilm.vilm_dv_min``,
reproducing Endgame Part-1 Tables 1/2), never a number the releg itself
computed. The ballistic regression test's EXPECTED side is the in-repo
``core.lambert`` path, the object the seam swaps OUT.
"""

from __future__ import annotations

import math

import numpy as np

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day, _moon_state
from cyclerfinder.search.releg_solver import BallisticReleg, RelegResult

DAY_S = 86400.0


def _jovian_leg(
    moon_a: str,
    moon_b: str,
    *,
    tof_days: float,
    theta_a: float = 0.0,
    theta_b: float = 1.1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Build a fixed planet-frame two-moon leg (positions/velocities, tof, mu).

    Circular-coplanar moon states from the registry (the same ``_moon_state``
    the discovery/validation seams use), evaluated at a fixed phasing so the
    leg geometry is deterministic.
    """
    mu = PRIMARIES["Jupiter"]
    sma_a = SATELLITES[moon_a].sma_km
    sma_b = SATELLITES[moon_b].sma_km
    n_a = _mean_motion_rad_day(mu, sma_a)
    n_b = _mean_motion_rad_day(mu, sma_b)
    r_a, v_a = _moon_state(theta_a, n_a, 0.0, sma_a, mu)
    # Arrival moon advanced to the leg's arrival epoch.
    r_b, v_b = _moon_state(theta_b, n_b, tof_days, sma_b, mu)
    return r_a, v_a, r_b, v_b, tof_days * DAY_S, mu


def _ballistic_vinf_reference(
    r_a: np.ndarray,
    v_a: np.ndarray,
    r_b: np.ndarray,
    v_b: np.ndarray,
    tof_s: float,
    mu: float,
    n_rev: int,
) -> tuple[float, float]:
    """The exact ``_close_one_phasing`` lowest-energy-branch V_inf reference.

    Lifted verbatim from ``discovery_campaign._close_one_phasing`` lines
    494-501: solve Lambert at the requested n_rev, keep that rev, pick the
    branch minimising departure V_inf, return ``(vinf_out, vinf_in)``.
    """
    sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
    wanted = [s for s in sols if s.n_revs == n_rev]
    best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
    vinf_out = float(np.linalg.norm(best.v1 - v_a))
    vinf_in = float(np.linalg.norm(best.v2 - v_b))
    return vinf_out, vinf_in


# ---------------------------------------------------------------------------
# Task 1 — BallisticReleg reproduces the ballistic Lambert leg bit-for-bit
# ---------------------------------------------------------------------------


def test_ballistic_releg_matches_lambert_path() -> None:
    """``BallisticReleg`` reproduces the lowest-energy Lambert branch exactly."""
    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Io", "Europa", tof_days=3.5)
    vinf_out_ref, vinf_in_ref = _ballistic_vinf_reference(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)

    result = BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)

    assert isinstance(result, RelegResult)
    assert result.feasible is True
    assert result.dv_kms == 0.0
    assert result.vinf_out == vinf_out_ref
    assert result.vinf_in == vinf_in_ref


def test_ballistic_releg_infeasible_when_no_branch() -> None:
    """No Lambert solution at the requested n_rev → not feasible (no crash)."""
    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Io", "Europa", tof_days=3.5)
    # A high multi-rev count that the short tof cannot admit.
    result = BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=8)
    assert result.feasible is False
    assert math.isinf(result.vinf_out) or result.vinf_out == 0.0
