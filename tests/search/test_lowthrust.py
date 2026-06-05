"""Tests for the Sims-Flanagan low-thrust optimiser integration (Phase 3).

Phase 3 of ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md``:
wire the leg model into the M5 DE+SLSQP pattern. The two-phase solve mirrors
Yam §1 — Phase 1 minimises total ΔV (at constant mass), Phase 2 re-optimises to
maximise final mass via the rocket equation.

Golden discipline: no literature anchor exists for the leg model (see the Phase
1/2 test headers and the plan's Phase 5 section). The gates here are physics
invariants only:

* the optimiser drives the match-point defect below tolerance on a
  self-consistent synthetic leg (plan §5 invariant 5);
* a powered solve never beats the zero-thrust closure that already exists when
  the leg is ballistically closed (it should find ~zero ΔV);
* Phase 2 final mass >= Phase 1 final mass (re-optimisation for mass cannot do
  worse) and the converged ΔV never exceeds the per-segment thrust bound.

The convergence tests run the full DE+SLSQP stack and are marked ``slow``; a
couple of fast structural tests stay in the default suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.sims_flanagan import (
    SimsFlanaganLeg,
    final_mass,
    leg_feasible,
    match_point_defect,
    segment_dv_bounds,
)
from cyclerfinder.search.lowthrust import (
    LowThrustLegResult,
    solve_leg_max_mass,
    solve_leg_min_dv,
)


def _ballistic_leg(n_segments: int = 6, tof_days: float = 200.0) -> SimsFlanaganLeg:
    """A leg whose endpoints already close ballistically (Kepler arc)."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = tof_days * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    return SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n_segments,
        m0_kg=1000.0,
        isp_s=3000.0,
        tmax_kn=0.5,
    )


def _offset_leg(n_segments: int = 6, tof_days: float = 200.0) -> SimsFlanaganLeg:
    """A leg whose end velocity is perturbed so a small burn is needed to close.

    The end state is the Kepler-propagated start state with a small velocity
    nudge, so the leg does NOT close ballistically but is reachable with a
    modest thrust budget. This is the self-consistent synthetic leg the
    optimiser must close (plan §5 invariant 5).
    """
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = tof_days * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    vf = vf + np.array([0.0, 0.05, 0.0], dtype=np.float64)  # 50 m/s nudge
    return SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n_segments,
        m0_kg=1000.0,
        isp_s=3000.0,
        tmax_kn=0.5,
    )


# ---------------------------------------------------------------------------
# Fast structural tests
# ---------------------------------------------------------------------------


def test_result_shape() -> None:
    """A solve returns a schedule of the right shape and a finite ΔV."""
    leg = _ballistic_leg(n_segments=4)
    result = solve_leg_min_dv(leg, seed=0, use_de=False, n_starts=1)
    assert isinstance(result, LowThrustLegResult)
    assert result.dvs.shape == (4, 3)
    assert np.isfinite(result.total_dv_kms)
    assert result.total_dv_kms >= 0.0


def test_zero_thrust_leg_rejects_solve() -> None:
    """A coast-only leg (tmax=0) admits only the zero schedule."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = 200.0 * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=4,
        m0_kg=1000.0,
        isp_s=3000.0,
        tmax_kn=0.0,
    )
    result = solve_leg_min_dv(leg, seed=0, use_de=False, n_starts=1)
    # No thrust capability => every segment bound is 0 => only zero schedule.
    assert float(np.max(np.linalg.norm(result.dvs, axis=1))) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Convergence invariants (slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_ballistic_leg_finds_near_zero_dv() -> None:
    """A ballistically-closed leg's min-ΔV solution is ~zero ΔV.

    The leg already closes with no thrust; the optimiser must not invent a
    powered solution that beats the trivial zero-ΔV closure (plan §5: a powered
    solve must never beat the theoretical bound).
    """
    leg = _ballistic_leg(n_segments=6)
    result = solve_leg_min_dv(leg, seed=0)
    assert result.converged
    # Feasible at the solver's tolerance (the residual floor is the Kepler /
    # half-segment-coast numerical scale, ~1e-2 km here, not a real burn).
    assert leg_feasible(leg, result.dvs, pos_tol_km=1.0, vel_tol_kms=1.0e-3)
    # Total ΔV should be tiny — the leg needs no real burn to close.
    assert result.total_dv_kms < 1.0e-2


@pytest.mark.slow
def test_offset_leg_defect_driven_to_tolerance() -> None:
    """The optimiser drives the match-point defect below tolerance.

    Plan §5 invariant 5: ``norm(S_mf - S_mb)`` below tol on a self-consistent
    synthetic leg.
    """
    leg = _offset_leg(n_segments=8)
    result = solve_leg_min_dv(leg, seed=0)
    assert result.converged
    assert leg_feasible(leg, result.dvs)
    defect = match_point_defect(leg, result.dvs, final_mass(leg, result.dvs))
    assert float(np.linalg.norm(defect[0:3])) < 1.0  # km
    assert float(np.linalg.norm(defect[3:6])) < 1.0e-3  # km/s


@pytest.mark.slow
def test_converged_schedule_respects_thrust_bound() -> None:
    """Every per-segment ΔV in the solution respects its capability bound."""
    leg = _offset_leg(n_segments=8)
    result = solve_leg_min_dv(leg, seed=0)
    assert result.converged
    bounds = segment_dv_bounds(leg, result.dvs)
    mags = np.linalg.norm(result.dvs, axis=1)
    # Allow a small numerical slack on the boundary.
    assert np.all(mags <= bounds + 1.0e-9)


@pytest.mark.slow
def test_max_mass_phase_not_worse_than_min_dv() -> None:
    """Phase 2 (max final mass) does not reduce final mass vs Phase 1.

    Re-optimising for final mass starting from the min-ΔV solution can only
    keep or improve the final mass while staying feasible.
    """
    leg = _offset_leg(n_segments=8)
    phase1 = solve_leg_min_dv(leg, seed=0)
    assert phase1.converged
    phase2 = solve_leg_max_mass(leg, phase1, seed=0)
    assert phase2.converged
    assert leg_feasible(leg, phase2.dvs)
    m1 = final_mass(leg, phase1.dvs)
    m2 = final_mass(leg, phase2.dvs)
    assert m2 >= m1 - 1.0e-6  # never meaningfully worse
