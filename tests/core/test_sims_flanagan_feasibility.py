"""Physics-invariant tests for the Sims-Flanagan feasibility / defect layer.

Phase 2 of ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md``.

These tests rest on physics invariants and *structural* (paper-dimension)
checks only — there is no usable literature anchor for the leg model yet (see
the Phase 1 test header and the plan's Phase 5 section). The load-bearing
invariants here are:

* a ballistically-closed leg (zero ΔV) is ``leg_feasible``;
* a leg with a real burn on only the forward half is **not** feasible;
* the assembled chain defect vector reduces to per-leg
  :func:`match_point_defect` blocks;
* the flyby-bend chain constraint reuses
  :func:`cyclerfinder.core.flyby.is_ballistic_feasible` exactly;
* the NLP dimension bookkeeping reproduces Yam's ``(8 + 3N)·M`` /
  ``neq·M`` (recorded in ``docs/v2-future-references.md`` §1).
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import (
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import is_ballistic_feasible
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.sims_flanagan import (
    FlybyBoundary,
    SimsFlanaganError,
    SimsFlanaganLeg,
    chain_defect,
    flyby_bend_slacks,
    leg_feasible,
    nlp_dimensions,
)


def _closed_leg(n_segments: int = 10, tof_days: float = 200.0) -> SimsFlanaganLeg:
    """A ballistically-closed leg: end state = Kepler-propagated start state."""
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
        m0_kg=500.0,
        isp_s=2500.0,
        tmax_kn=0.1,
    )


def _zero_schedule(n: int) -> np.ndarray:
    return np.zeros((n, 3), dtype=np.float64)


# ---------------------------------------------------------------------------
# leg_feasible
# ---------------------------------------------------------------------------


def test_zero_thrust_closed_leg_is_feasible() -> None:
    leg = _closed_leg()
    dvs = _zero_schedule(leg.n_segments)
    assert leg_feasible(leg, dvs) is True


def test_forward_burn_breaks_feasibility() -> None:
    leg = _closed_leg()
    dvs = _zero_schedule(leg.n_segments)
    dvs[1] = np.array([0.02, 0.0, 0.0])  # burn on the forward half only
    assert leg_feasible(leg, dvs) is False


def test_leg_feasible_respects_tolerance() -> None:
    """A loose enough tolerance accepts a leg a tight one rejects."""
    leg = _closed_leg()
    dvs = _zero_schedule(leg.n_segments)
    dvs[1] = np.array([1.0e-7, 0.0, 0.0])  # a tiny burn
    assert leg_feasible(leg, dvs, pos_tol_km=1.0e-9, vel_tol_kms=1.0e-12) is False
    assert leg_feasible(leg, dvs, pos_tol_km=1.0e6, vel_tol_kms=1.0e6) is True


def test_leg_feasible_rejects_mass_mismatch() -> None:
    """Supplying an inconsistent end mass opens the mass defect."""
    leg = _closed_leg()
    dvs = _zero_schedule(leg.n_segments)
    # Consistent mf is m0 (no burns); pass a wrong mf and demand a tight mass tol.
    assert leg_feasible(leg, dvs, mf_kg=leg.m0_kg * 0.5, mass_tol_kg=1.0e-6) is False
    assert leg_feasible(leg, dvs, mf_kg=leg.m0_kg, mass_tol_kg=1.0e-9) is True


# ---------------------------------------------------------------------------
# chain_defect
# ---------------------------------------------------------------------------


def test_chain_defect_concatenates_leg_defects() -> None:
    """The chain defect over M legs is the stack of each leg's 7-vector."""
    leg_a = _closed_leg(n_segments=8, tof_days=180.0)
    leg_b = _closed_leg(n_segments=6, tof_days=220.0)
    schedules = [_zero_schedule(8), _zero_schedule(6)]
    defect = chain_defect([leg_a, leg_b], schedules)
    assert defect.shape == (14,)  # 7 per leg, 2 legs
    # All-zero schedules on closed legs => defect ~ 0.
    assert float(np.linalg.norm(defect)) < 1.0e-2


def test_chain_defect_rejects_length_mismatch() -> None:
    leg = _closed_leg(n_segments=8)
    with pytest.raises(SimsFlanaganError, match="legs"):
        chain_defect([leg], [_zero_schedule(8), _zero_schedule(8)])


# ---------------------------------------------------------------------------
# flyby_bend_slacks (reuses flyby.is_ballistic_feasible / max_bend)
# ---------------------------------------------------------------------------


def test_flyby_bend_slack_zero_bend_is_feasible() -> None:
    """An identical in/out V_inf needs no bend => positive (feasible) slack."""
    vinf = np.array([3.0, 0.0, 0.0], dtype=np.float64)
    boundary = FlybyBoundary(
        body="E",
        vinf_in=vinf,
        vinf_out=vinf,
        mu_planet=PLANETS["E"].mu_km3_s2,
        rp_min=SAFE_PERIHELION_KM["E"],
    )
    slacks = flyby_bend_slacks([boundary])
    assert slacks.shape == (1,)
    assert slacks[0] >= 0.0
    # Cross-check: the same predicate flyby.py exposes.
    assert is_ballistic_feasible(vinf, vinf, PLANETS["E"].mu_km3_s2, SAFE_PERIHELION_KM["E"])


def test_flyby_bend_slack_excess_bend_is_infeasible() -> None:
    """A near-reversal at high V_inf cannot be bent => negative slack."""
    vin = np.array([20.0, 0.0, 0.0], dtype=np.float64)
    vout = np.array([-20.0, 0.5, 0.0], dtype=np.float64)  # ~180 deg, same speed-ish
    boundary = FlybyBoundary(
        body="E",
        vinf_in=vin,
        vinf_out=vout,
        mu_planet=PLANETS["E"].mu_km3_s2,
        rp_min=SAFE_PERIHELION_KM["E"],
    )
    slacks = flyby_bend_slacks([boundary])
    assert slacks[0] < 0.0
    assert not is_ballistic_feasible(vin, vout, PLANETS["E"].mu_km3_s2, SAFE_PERIHELION_KM["E"])


def test_flyby_bend_slack_sign_matches_predicate() -> None:
    """slack >= 0 iff is_ballistic_feasible, across a sweep of bend angles."""
    speed = 5.0
    mu = PLANETS["E"].mu_km3_s2
    rp = SAFE_PERIHELION_KM["E"]
    for deg in (0.0, 10.0, 30.0, 60.0, 90.0, 150.0):
        ang = np.radians(deg)
        vin = np.array([speed, 0.0, 0.0], dtype=np.float64)
        vout = np.array([speed * np.cos(ang), speed * np.sin(ang), 0.0], dtype=np.float64)
        boundary = FlybyBoundary(body="E", vinf_in=vin, vinf_out=vout, mu_planet=mu, rp_min=rp)
        slack = float(flyby_bend_slacks([boundary])[0])
        feasible = is_ballistic_feasible(vin, vout, mu, rp)
        assert (slack >= 0.0) == feasible, f"mismatch at {deg} deg"


# ---------------------------------------------------------------------------
# NLP dimension bookkeeping (Yam: (8 + 3N)*M variables, neq*M constraints)
# ---------------------------------------------------------------------------


def test_nlp_dimensions_match_yam_formula() -> None:
    n_segments, n_legs = 20, 1  # Yam E-E-J: NLP dim 75, 35 constraints
    dims = nlp_dimensions(n_segments=n_segments, n_legs=n_legs)
    # Variables: (8 + 3N)*M ; with N segments and M legs.
    assert dims.n_variables == (8 + 3 * n_segments) * n_legs
    # Constraints: neq*M with neq = 7 (3D position + velocity + mass).
    assert dims.n_constraints == 7 * n_legs
    assert dims.neq == 7


def test_nlp_dimensions_two_legs() -> None:
    dims = nlp_dimensions(n_segments=10, n_legs=2)
    assert dims.n_variables == (8 + 30) * 2
    assert dims.n_constraints == 14


def test_nlp_dimensions_rejects_nonpositive() -> None:
    with pytest.raises(SimsFlanaganError, match="n_segments"):
        nlp_dimensions(n_segments=0, n_legs=1)
    with pytest.raises(SimsFlanaganError, match="n_legs"):
        nlp_dimensions(n_segments=10, n_legs=0)
