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
import pytest

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day, _moon_state
from cyclerfinder.search.releg_solver import (
    BallisticReleg,
    DsmReleg,
    LowThrustReleg,
    RelegResult,
)
from cyclerfinder.search.vilm import vilm_dv_min

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


# ---------------------------------------------------------------------------
# Task 2 — DsmReleg: V_inf retarget capability + VILM-floor golden
# ---------------------------------------------------------------------------


def test_dsm_releg_retargets_vinf() -> None:
    """The powered DSM leg drives arrival V_inf to a requested value, dv > 0.

    On a Jovian Io->Europa leg whose ballistic arrival V_inf differs from the
    requested target, ``DsmReleg`` retargets to the target (within tolerance)
    and spends a strictly-positive ΔV to do so.
    """
    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Io", "Europa", tof_days=5.0)
    ballistic = BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)
    # A target arrival V_inf well away from the ballistic value.
    target = 5.0
    assert abs(ballistic.vinf_in - target) > 0.5  # the leg actually needs retargeting

    result = DsmReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0, vinf_target_in=target)
    assert result.feasible is True
    assert result.vinf_in == pytest.approx(target, abs=0.15)
    assert result.dv_kms > 0.0


def test_dsm_releg_min_dv_recovers_ballistic_close() -> None:
    """With no target, the cheapest powered close is ~ballistic (dv ≈ 0).

    A powered leg strictly subsumes the ballistic leg: where the ballistic
    Lambert already closes, the minimum-ΔV DSM delivers a near-zero impulse.
    """
    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Io", "Europa", tof_days=5.0)
    result = DsmReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)
    assert result.feasible is True
    assert result.dv_kms == pytest.approx(0.0, abs=1.0e-2)


def test_dsm_releg_dv_geq_vilm_floor() -> None:
    """GOLDEN (sourced): a Ganymede->Europa powered V_inf-lowering releg's
    delivered ΔV is bounded below by the published VILM leveraging floor.

    EXPECTED side: ``vilm.vilm_dv_min("Ganymede","Europa")`` — the analytic
    theoretical-minimum VILM transfer ΔV that reproduces Campagnola-Russell
    Endgame Part-1 Table 1 (Ganymede-Europa ΔV_min = 1.71 km/s, mining note
    ``2026-06-05-endgame-tisserand-mining.md`` A2). A single impulsive in-leg
    DSM that performs the inter-moon V_inf change CANNOT beat the
    theoretical-minimum multi-VILM ΔV, so the delivered ΔV must be ≥ that floor
    (the floor is the sourced lower bound; the DSM is the over-cost realisation).
    The golden is non-circular: the EXPECTED value is the published table value,
    never a number the releg itself produced.
    """
    floor = vilm_dv_min("Ganymede", "Europa")
    # Sanity: the floor reproduces the published 1.71 km/s Table 1 value.
    assert floor == pytest.approx(1.71, abs=0.17)

    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Ganymede", "Europa", tof_days=6.0)
    ballistic = BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)
    # A real V_inf-lowering endgame-like step: retarget arrival V_inf well down.
    target = 0.3 * ballistic.vinf_in
    result = DsmReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0, vinf_target_in=target)
    assert result.feasible is True
    assert result.vinf_in == pytest.approx(target, abs=0.2)
    # The delivered ΔV is at or above the sourced VILM floor (an in-leg DSM is
    # the over-cost realisation of the same inter-moon V_inf change).
    assert result.dv_kms >= floor - 1.0e-3
    # ...and inside the powered band's sanity ceiling (3.5 km/s/cycle x reach).
    assert result.dv_kms < 25.0


# ---------------------------------------------------------------------------
# Task 7 — LowThrustReleg: Sims-Flanagan backend, regression + bracket golden
# ---------------------------------------------------------------------------


def test_lowthrust_releg_zero_retarget_matches_ballistic() -> None:
    """REGRESSION (coplanar limit): a zero-retarget SF leg ≈ ballistic (dv ≈ 0).

    With neither departure nor arrival V_inf retargeted, the SF leg's boundary
    states ARE the ballistic Lambert endpoints, so its all-zero thrust schedule
    closes at ΔV ≈ 0 — the SF leg's zero-thrust limit is the ballistic leg (the
    same regression anchor :class:`BallisticReleg` provides). The reported V_inf
    chain matches the ballistic leg's lowest-energy branch.
    """
    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Io", "Europa", tof_days=5.0)
    ballistic = BallisticReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)

    result = LowThrustReleg().solve(r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0)
    assert result.feasible is True
    # The SF leg's zero-thrust limit closes at ~zero ΔV (the ballistic leg).
    assert result.dv_kms == pytest.approx(0.0, abs=1.0e-2)
    # V_inf chain matches the ballistic lowest-energy branch (same boundary states).
    assert result.vinf_out == pytest.approx(ballistic.vinf_out, abs=1.0e-3)
    assert result.vinf_in == pytest.approx(ballistic.vinf_in, abs=1.0e-3)


def test_lowthrust_releg_brackets_dsm() -> None:
    """BRACKET GOLDEN (sourced floor): the SF low-thrust ΔV brackets DSM/VILM-floor.

    On the Ganymede->Europa leg, retarget both departure and arrival V_inf to a
    common ``T`` (the moon-tour driver's continuity-by-construction contract).
    The SF low-thrust delivered ΔV must:

    * be ≥ the published VILM leveraging floor ``vilm.vilm_dv_min`` (a powered leg
      cannot beat the theoretical-minimum VILM ΔV for the transfer — the SOURCED,
      non-circular lower bound that reproduces Campagnola-Russell Part-1 Table 1
      Ganymede-Europa ΔV_min = 1.71 km/s);
    * be ≤ the single-impulse :class:`DsmReleg` cost (distributing the ΔV across a
      thrust train is strictly more efficient than one impulse — the physics the
      whole #449 low-thrust bet rests on).

    The SF leg model has no clean state-level literature anchor (the Vasile-
    Campagnola 2009 DFET transcription "DOES NOT MAP" to our SF leg, digest
    ``2026-06-07-vasile-campagnola-dfet-method-mining.md`` §2.6), so this is a
    bracket, not an equality — and the bracket's lower edge is the SAME sourced
    VILM floor the tight DSM golden uses, never a number SF itself computed.
    """
    floor = vilm_dv_min("Ganymede", "Europa")
    assert floor == pytest.approx(1.71, abs=0.17)  # reproduces the published value

    r_a, v_a, r_b, v_b, tof_s, mu = _jovian_leg("Ganymede", "Europa", tof_days=6.0)
    target = 4.0
    dsm = DsmReleg().solve(
        r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0, vinf_target_in=target, vinf_depart_mag=target
    )
    assert dsm.feasible is True

    result = LowThrustReleg().solve(
        r_a, v_a, r_b, v_b, tof_s, mu, n_rev=0, vinf_target_in=target, vinf_depart_mag=target
    )
    assert result.feasible is True
    # Lower edge: the sourced VILM floor (cannot beat the theoretical minimum).
    assert result.dv_kms >= floor - 1.0e-3
    # Upper edge: at or below the single-impulse DSM cost (distributed is cheaper).
    # A 10% slack absorbs the different-transcription / convergence-tolerance gap.
    assert result.dv_kms <= dsm.dv_kms * 1.1
    # The retargeted V_inf chain is pinned to the common flyby target by construction.
    assert result.vinf_out == pytest.approx(target, abs=1.0e-6)
    assert result.vinf_in == pytest.approx(target, abs=1.0e-6)
