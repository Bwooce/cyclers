"""Tests for the circular -> ephemeris continuation driver (task #158).

Two layers:

* FAST mechanics gates (circular / ramped-element only, no astropy):
  - λ=0 reproduces ``Ephemeris('circular')`` bit-identically;
  - a single e-ramp step moves the planet state continuously (no jump);
  - the i-ramp tilts the plane while preserving the orbit radius;
  - ladder bookkeeping (schedule shape, skipped recording, keep-best).

* SLOW GOLDEN gates (``@pytest.mark.slow``, DE440 via astropy) — the #158
  acceptance tests, sourced from Russell 2004 dissertation Table 5.5 (p.178):
  - Aldrin 6.399G1, launch Aug 2003: total maintenance Δv = 0 m/s / 7 cycles;
  - 4.991gG2 (S1L1), launch Jun 2025: 0 m/s.
  GOLDEN DISCIPLINE: the EXPECTED side is the sourced "ballistic (0 m/s)" claim;
  our evidence is the continued solution's final-step residual being below the
  documented closure tolerance (``TOL_KMS``), i.e. "≈0 within tol" — never a
  fabricated digit. The emerged V_inf (compared to the INDEPENDENTLY sourced
  Russell anchors) is corroborating evidence, never imposed.

The "tolerance" for the ballistic claim is the corrector's km/s closure floor
(``free_return_correct`` default ``tol_kms=0.1``): Russell's "0 m/s" means below
his SNOPT/post-processing residual floor (deep-dive §7: exact SNOPT tolerances
are not printed), so we state OUR floor explicitly and report the achieved
residual as evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime, timedelta

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import continuation as cont

# Closure tolerance (km/s) for the "ballistic ≈ 0 within tol" golden gate. This
# is OUR documented floor (the free-return corrector default); Russell's "0 m/s"
# is below his unprinted SNOPT/post-processing floor.
TOL_KMS = 0.1

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def _epoch_sec(dt: datetime) -> float:
    return (dt - _J2000).total_seconds()


# ---------------------------------------------------------------------------
# FAST mechanics gates
# ---------------------------------------------------------------------------


def test_lambda_zero_reproduces_circular_bit_identically() -> None:
    """λ_e = λ_i = λ_p = 0 backend == ``Ephemeris('circular')``, byte for byte."""
    circ = Ephemeris("circular")
    ramp0 = cont.ramped_ephemeris(0.0, 0.0, 0.0)
    for body in ("E", "M"):
        for t in (0.0, 1.0e7, -3.3e7, 2.5e8):
            rc, vc = circ.state(body, t)
            rr, vr = ramp0.state(body, t)
            assert np.array_equal(rc, rr), (body, t)
            assert np.array_equal(vc, vr), (body, t)


def test_ramped_ephemeris_default_is_circular() -> None:
    """The public ``ramped_ephemeris`` defaults (lam_p=0) to the circular model."""
    circ = Ephemeris("circular")
    ramp = cont.ramped_ephemeris(0.0, 0.0)
    rc, vc = circ.state("M", 1.234e8)
    rr, vr = ramp.state("M", 1.234e8)
    assert np.array_equal(rc, rr)
    assert np.array_equal(vc, vr)


def test_e_ramp_step_moves_state_continuously() -> None:
    """A single small e-ramp step perturbs Mars's position smoothly (no jump).

    Mechanics (not golden): successive equal λ_e increments produce bounded,
    monotone-in-magnitude position deltas — there is no discontinuity between the
    circular model and the first eccentric step (the continuation premise).
    """
    t = 1.0e7
    prev = None
    deltas = []
    for lam in np.linspace(0.0, 1.0, 11):
        r, _ = cont.ramped_ephemeris(float(lam), 0.0, 0.0).state("M", t)
        # e-ramp at λ_i=λ_p=0 stays in the ecliptic plane.
        assert r[2] == 0.0
        if prev is not None:
            deltas.append(float(np.linalg.norm(r - prev)))
        prev = r
    # No single step is more than ~3x any other (continuity, no jump).
    assert max(deltas) < 3.0 * min(deltas)
    # And the cumulative move is non-trivial (e actually got introduced).
    assert sum(deltas) > 1.0e6  # km


def test_i_ramp_tilts_plane_preserving_radius() -> None:
    """The i-ramp lifts Mars off the ecliptic continuously, |r| invariant."""
    t = 5.0e7
    r0, _ = cont.ramped_ephemeris(1.0, 0.0, 1.0).state("M", t)
    assert abs(r0[2]) < 1.0  # in-plane at λ_i=0 (km)
    rmag0 = float(np.linalg.norm(r0))
    prev_z = 0.0
    for lam in np.linspace(0.1, 1.0, 10):
        r, _ = cont.ramped_ephemeris(1.0, float(lam), 1.0).state("M", t)
        # |r| preserved (the conjugated tilt is a pure rotation).
        assert abs(float(np.linalg.norm(r)) - rmag0) < 1.0  # km
        # z grows monotonically off the ecliptic.
        assert abs(r[2]) >= abs(prev_z)
        prev_z = r[2]
    assert abs(prev_z) > 1.0e6  # km, a real tilt at λ_i=1


def test_phase_ramp_aligns_with_de440_anchor_in_plane() -> None:
    """λ_p=1 sets Mars's J2000 longitude to the frozen Table-5.4 mean longitude.

    Mechanics gate (the phase-anchor is what lets the final DE440 step be small):
    at t=0 the ramped (λ_p=1) model's Mars heliocentric longitude equals
    ``Omega + argp + nu`` from Table 5.4 (~355.4°), not the circular 0°.
    """
    r, _ = cont.ramped_ephemeris(0.0, 0.0, 1.0).state("M", 0.0)
    lon_deg = float(np.degrees(np.arctan2(r[1], r[0]))) % 360.0
    assert abs(((lon_deg - 355.43 + 180.0) % 360.0) - 180.0) < 1.0


def test_ramp_schedule_shape() -> None:
    """The schedule is p-ramp(nstep) + e-ramp(nstep) + i-ramp(nstep)."""
    sched = cont._ramp_schedule(3)
    phases = [p for p, *_ in sched]
    assert phases == ["p-ramp"] * 3 + ["e-ramp"] * 3 + ["i-ramp"] * 3
    # Each phase ends at its target λ; the others are at their held value.
    assert sched[2] == ("p-ramp", 0.0, 0.0, 1.0)
    assert sched[5] == ("e-ramp", 1.0, 0.0, 1.0)
    assert sched[8] == ("i-ramp", 1.0, 1.0, 1.0)


def test_ladder_records_skipped_and_keeps_best() -> None:
    """Driver runs the requested ladder, records the skipped rungs, and keeps the
    lowest-residual completed final. Uses a ramped (1,1) target (fast, no DE440)."""
    # A circular-coplanar seed that the free-return genome closes (S1L1 ellipse).
    a, e = 1.30, 0.257
    period_sec = 4.27 * 365.25 * 86400.0
    result = cont.continuation_correct(
        0.0,
        a,
        e,
        period_sec,
        ladder=(1, 3),
        final_ephemeris=cont.ramped_ephemeris(1.0, 1.0, 1.0),
    )
    assert result.ladder == (1, 3)
    # 243 plus the un-run middle rungs are recorded, never silently dropped.
    assert set(result.skipped) == {9, 27, 81, 243}
    assert len(result.rungs) == 2
    # best_final is the minimum-residual completed rung.
    completed = [r for r in result.rungs if r.completed]
    assert completed
    best = min(completed, key=lambda r: r.final.max_residual_kms)
    assert result.winning_nstep == best.nstep
    assert result.best_final.max_residual_kms == best.final.max_residual_kms


def test_step_audit_trail_is_complete() -> None:
    """Every rung carries one ContinuationStep per schedule entry plus the
    ephemeris step, each with its λ trio and residual (the audit trail)."""
    result = cont.continuation_correct(
        0.0,
        1.30,
        0.257,
        4.27 * 365.25 * 86400.0,
        ladder=(1,),
        final_ephemeris=cont.ramped_ephemeris(1.0, 1.0, 1.0),
    )
    rung = result.rungs[0]
    assert rung.completed
    # nstep=1 -> 1 p + 1 e + 1 i + 1 ephemeris = 4 steps.
    assert [s.phase for s in rung.steps] == ["p-ramp", "e-ramp", "i-ramp", "ephemeris"]
    for s in rung.steps:
        assert np.isfinite(s.max_residual_kms)
        assert 0.0 <= s.lam_e <= 1.0 and 0.0 <= s.lam_i <= 1.0 and 0.0 <= s.lam_p <= 1.0


# ---------------------------------------------------------------------------
# SLOW GOLDEN gates (DE440 via astropy) — Russell 2004 Table 5.5, p.178
# ---------------------------------------------------------------------------


# Sourced ellipse + launch window + V_inf anchors per golden. The ellipse (a,e)
# is the SOURCED constraint; the launch epoch is the sourced window; the V_inf
# anchors are INDEPENDENTLY sourced (Russell Table 3.4 / 4.9) corroborating
# evidence. Period is the heliocentric-ellipse period (Kepler from a), the
# free-return genome's natural period.
@dataclass(frozen=True)
class _Golden:
    rid: str
    a_au: float
    e: float
    period_yr: float
    launch: datetime
    seed: datetime
    vinf_e: float
    vinf_m: float


_ALDRIN = _Golden(
    rid="6.399G1 (Aldrin)",
    a_au=1.60,  # Rogers 2012 Table 1 (catalogue aldrin-classic-em-k1)
    e=0.393,
    period_yr=2.0239,  # 2π√(a³/μ), the ellipse period
    launch=datetime(2003, 8, 6, tzinfo=UTC),  # Russell p.176 / Fig.5.7
    seed=datetime(2003, 7, 11, tzinfo=UTC),  # DE440 basin near the sourced window
    vinf_e=6.5,  # Russell Table 3.4 cycler 1.0.1.-1
    vinf_m=9.7,
)
_S1L1 = _Golden(
    rid="4.991gG2 (S1L1)",
    a_au=1.30,  # Rogers 2012 Table 1 (catalogue russell-ch4-4.991gG2 sister)
    e=0.257,
    period_yr=4.27,
    launch=datetime(2025, 6, 15, tzinfo=UTC),  # Russell p.176 (Jun 2025)
    seed=datetime(2026, 12, 10, tzinfo=UTC),  # DE440 basin in the repeat window
    vinf_e=4.99,  # Russell Table 4.9
    vinf_m=5.10,
)


@pytest.mark.slow
@pytest.mark.parametrize("g", [_ALDRIN, _S1L1], ids=["aldrin-6.399G1", "s1l1-4.991gG2"])
def test_golden_continuation_is_ballistic_within_tol(g: _Golden) -> None:
    """#158 GOLDEN: the continued solution at the true-ephemeris (DE440) step is
    ballistic within the documented closure tolerance for the sourced launch
    window — the evidence for Russell's "total maintenance Δv = 0 m/s / 7 cycles".

    EXPECTED (sourced, Russell Table 5.5): 0 m/s.  EVIDENCE: ``best_final``
    residual < TOL_KMS (≈0 within tol). If this fails it is a reportable
    scientific result (recorded in the results note with the break point), NOT a
    gate to soften.
    """
    period_sec = g.period_yr * 365.25 * 86400.0
    result = cont.continuation_correct(
        _epoch_sec(g.seed),
        g.a_au,
        g.e,
        period_sec,
        ladder=(1, 3, 9, 27, 81),  # 243 skipped (wall-time; recorded in .skipped)
        final_ephemeris=Ephemeris("astropy"),
    )
    bf = result.best_final

    # The headline gate: ballistic ≈ 0 within the documented tolerance.
    assert bf.converged, (
        f"{g.rid}: not ballistic within tol; best-final residual "
        f"{bf.max_residual_kms:.4f} km/s > {TOL_KMS} km/s at nstep={result.winning_nstep}"
    )
    assert bf.max_residual_kms < TOL_KMS

    # The continued solution lands in the sourced launch window (within ~6 months).
    launched = _J2000 + timedelta(seconds=bf.t0_sec)
    assert abs((launched - g.launch).days) < 200

    # Corroborating EVIDENCE (derived, never imposed): emerged V_inf is in the
    # neighbourhood of the independently sourced anchors. Coplanar single-ellipse
    # free-return vs the real 3-D ephemeris: a ~1.5 km/s tolerance acknowledges
    # the genome/model gap while still pinning the right family.
    assert abs(bf.vinf_kms["E"] - g.vinf_e) <= 1.5
    assert abs(bf.vinf_kms["M"] - g.vinf_m) <= 1.5

    # 243 was skipped for wall-time and must be recorded (no silent cap).
    assert 243 in result.skipped
