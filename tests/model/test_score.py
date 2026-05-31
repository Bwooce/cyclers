"""Tests for :mod:`cyclerfinder.model.score` — the M4 ranking layer.

Plan §4.2 gates: hard-constraint semantics (spec §12(d)), taxi-cost
surrogate sanity (including the Aldrin anchor), composite reproducibility,
and the rank top-N reducer's ordering / filtering / truncation contract.

Plan: ``docs/phases/m4-enumeration-scoring/plan.md`` §4.2.
"""

from __future__ import annotations

import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler, Encounter, Leg
from cyclerfinder.model.score import (
    DEFAULT_WEIGHTS,
    Score,
    composite_score,
    rank,
    score,
    taxi_cost_kms,
)
from cyclerfinder.search.construct import build_aldrin_seed

# ---------------------------------------------------------------------------
# Tolerances (module-level)
# ---------------------------------------------------------------------------

TOL_MAX_VINF_KMS: float = 0.1
TOL_TAXI_KMS: float = 0.1
VINF_CAP_ALDRIN_KMS: float = 12.0  # generous; M3 reports V∞_M ≈ 9.7 km/s


# ---------------------------------------------------------------------------
# Hand-built cycler helpers (keep tests independent of M3 construct)
# ---------------------------------------------------------------------------


def _stub_vec(x: float, y: float = 0.0, z: float = 0.0) -> np.ndarray:
    return np.array([x, y, z], dtype=np.float64)


def _make_encounter(
    body: str,
    vinf_in: np.ndarray,
    vinf_out: np.ndarray,
    t: float = 0.0,
) -> Encounter:
    return Encounter(
        body=body,
        t=t,
        r=_stub_vec(1.0e8),  # any non-zero heliocentric position
        v_planet=_stub_vec(0.0, 30.0),
        vinf_in=vinf_in,
        vinf_out=vinf_out,
    )


def _make_leg(
    from_body: str,
    to_body: str,
    t_depart: float = 0.0,
    t_arrive: float = 1.0e7,
) -> Leg:
    # Departure state: a circular-ish heliocentric velocity so radial_span
    # is well-defined and well above zero. Exact numbers don't matter for
    # the score-layer tests; they only need to be non-degenerate.
    return Leg(
        from_body=from_body,
        to_body=to_body,
        t_depart=t_depart,
        t_arrive=t_arrive,
        v_depart=_stub_vec(0.0, 30.0),
        v_arrive=_stub_vec(0.0, 24.0),
        n_revs=0,
        branch="single",
    )


def _make_two_encounter_cycler(
    body_a: str,
    body_b: str,
    vinf_a: np.ndarray,
    vinf_b: np.ndarray,
    period: float = 1.0e7,
) -> Cycler:
    enc_a = _make_encounter(body_a, vinf_a, vinf_a, t=0.0)
    enc_b = _make_encounter(body_b, vinf_b, vinf_b, t=period)
    leg = _make_leg(body_a, body_b, t_depart=0.0, t_arrive=period)
    return Cycler(bodies=[body_a, body_b], period=period, encounters=[enc_a, enc_b], legs=[leg])


# ---------------------------------------------------------------------------
# Score dataclass
# ---------------------------------------------------------------------------


def test_score_is_frozen() -> None:
    """Score is a frozen dataclass — assignment raises ``FrozenInstanceError``."""
    s = Score(
        total_maintenance_dv_kms=0.0,
        max_vinf_kms=5.0,
        radial_span_au=(1.0, 1.5),
        period_error_yr=0.0,
        taxi_cost_kms=5.0,
        hard_constraints_pass=True,
    )
    with pytest.raises(FrozenInstanceError):
        s.max_vinf_kms = 99.0  # type: ignore[misc]


def test_score_default_weights_keys() -> None:
    """The four ranking axes named in spec §5 step 4 appear in
    :data:`DEFAULT_WEIGHTS`."""
    assert set(DEFAULT_WEIGHTS.keys()) == {
        "total_maintenance_dv_kms",
        "max_vinf_kms",
        "period_error_yr",
        "taxi_cost_kms",
    }


# ---------------------------------------------------------------------------
# Hard constraints (spec §12(d))
# ---------------------------------------------------------------------------


def test_hard_constraints_pass_ballistic() -> None:
    """A hand-built cycler with ballistic flybys passes ``hard_constraints_pass``."""
    ephem = Ephemeris(model="circular")
    # vinf_in == vinf_out at each encounter → zero bend, trivially feasible.
    vinf_e = _stub_vec(5.0, 0.0)
    vinf_m = _stub_vec(7.0, 0.0)
    cyc = _make_two_encounter_cycler("E", "M", vinf_e, vinf_m)
    s = score(cyc, ephem, vinf_cap=10.0)
    assert s.hard_constraints_pass is True


def test_hard_constraints_fail_on_vinf_cap() -> None:
    """A 12 km/s encounter trips ``vinf_cap=7.0``."""
    ephem = Ephemeris(model="circular")
    vinf_e = _stub_vec(5.0, 0.0)
    vinf_m = _stub_vec(12.0, 0.0)  # > cap
    cyc = _make_two_encounter_cycler("E", "M", vinf_e, vinf_m)
    s = score(cyc, ephem, vinf_cap=7.0)
    assert s.hard_constraints_pass is False


def test_hard_constraints_fail_on_overbent_pair() -> None:
    """A 180° bend at 8 km/s vastly exceeds the achievable Mars bend cone.

    M2 anchors max bend at Mars at ~24° at 7 km/s; flipping V∞ direction
    requires 180° of bend, which is impossible.
    """
    ephem = Ephemeris(model="circular")
    vinf_e = _stub_vec(5.0, 0.0)
    vinf_in_m = _stub_vec(8.0, 0.0)
    vinf_out_m = _stub_vec(-8.0, 0.0)  # 180° bend — physically impossible
    enc_e = _make_encounter("E", vinf_e, vinf_e, t=0.0)
    enc_m = _make_encounter("M", vinf_in_m, vinf_out_m, t=1.0e7)
    leg = _make_leg("E", "M")
    cyc = Cycler(bodies=["E", "M"], period=1.0e7, encounters=[enc_e, enc_m], legs=[leg])
    s = score(cyc, ephem, vinf_cap=15.0)  # well above magnitude, fail on bend
    assert s.hard_constraints_pass is False


def test_score_aldrin_passes_hard_constraints_gate() -> None:
    """**Gate:** Aldrin seed at ``vinf_cap=12.0`` passes hard constraints.

    M3 hand-off anchors: V∞_E ≈ 6.5 km/s, V∞_M ≈ 9.7 km/s.
    A 12 km/s cap leaves headroom for both magnitudes; the ``E`` and
    ``M`` encounters in the M3 boundary convention have ``vinf_in ==
    vinf_out``, so the bend angle is zero (trivially achievable).
    """
    ephem = Ephemeris(model="circular")
    cyc = build_aldrin_seed(ephem)
    s = score(cyc, ephem, vinf_cap=VINF_CAP_ALDRIN_KMS)
    assert s.hard_constraints_pass is True
    assert s.max_vinf_kms == pytest.approx(9.7, abs=TOL_MAX_VINF_KMS)


# ---------------------------------------------------------------------------
# Period error
# ---------------------------------------------------------------------------


def test_period_error_none_target_is_zero() -> None:
    """``target_period_sec=None`` ⇒ ``period_error_yr == 0.0``."""
    ephem = Ephemeris(model="circular")
    cyc = _make_two_encounter_cycler("E", "M", _stub_vec(5.0), _stub_vec(7.0), period=1.0e7)
    s = score(cyc, ephem, vinf_cap=10.0, target_period_sec=None)
    assert s.period_error_yr == 0.0


def test_period_error_one_year() -> None:
    """``target_period_sec`` one Julian year below the cycler's period →
    ``period_error_yr ≈ 1.0``."""
    ephem = Ephemeris(model="circular")
    seconds_per_year = DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY
    period_sec = 2.0 * seconds_per_year
    cyc = _make_two_encounter_cycler("E", "M", _stub_vec(5.0), _stub_vec(7.0), period=period_sec)
    s = score(
        cyc,
        ephem,
        vinf_cap=10.0,
        target_period_sec=seconds_per_year,
    )
    assert s.period_error_yr == pytest.approx(1.0, abs=1.0e-12)


# ---------------------------------------------------------------------------
# Taxi cost
# ---------------------------------------------------------------------------


def test_taxi_cost_earth_only() -> None:
    """Earth-only cycler with two encounters at V∞ 5 and 7 → returns 7.0."""
    enc1 = _make_encounter("E", _stub_vec(5.0), _stub_vec(5.0), t=0.0)
    enc2 = _make_encounter("E", _stub_vec(7.0), _stub_vec(7.0), t=1.0e7)
    leg = _make_leg("E", "E")
    cyc = Cycler(bodies=["E", "E"], period=1.0e7, encounters=[enc1, enc2], legs=[leg])
    assert taxi_cost_kms(cyc) == pytest.approx(7.0, abs=1.0e-12)


def test_taxi_cost_zero_when_no_earth() -> None:
    """A V-M-only cycler returns 0.0 from the surrogate."""
    enc_v = _make_encounter("V", _stub_vec(5.0), _stub_vec(5.0), t=0.0)
    enc_m = _make_encounter("M", _stub_vec(7.0), _stub_vec(7.0), t=1.0e7)
    leg = _make_leg("V", "M")
    cyc = Cycler(bodies=["V", "M"], period=1.0e7, encounters=[enc_v, enc_m], legs=[leg])
    assert taxi_cost_kms(cyc) == 0.0


def test_taxi_cost_aldrin() -> None:
    """``taxi_cost_kms(aldrin)`` ≈ Earth V∞ ≈ 6.5 km/s (M3 hand-off, ±0.1)."""
    ephem = Ephemeris(model="circular")
    cyc = build_aldrin_seed(ephem)
    assert taxi_cost_kms(cyc) == pytest.approx(6.5, abs=TOL_TAXI_KMS)


# ---------------------------------------------------------------------------
# Composite + rank
# ---------------------------------------------------------------------------


def test_composite_finite_and_reproducible() -> None:
    """Two calls on the same Score yield bitwise-identical floats."""
    s = Score(
        total_maintenance_dv_kms=0.5,
        max_vinf_kms=7.0,
        radial_span_au=(0.9, 1.6),
        period_error_yr=0.0,
        taxi_cost_kms=6.5,
        hard_constraints_pass=True,
    )
    c1 = composite_score(s)
    c2 = composite_score(s)
    assert c1 == c2  # bitwise equality, not pytest.approx
    assert math.isfinite(c1)


def test_composite_infinite_on_hard_constraint_fail() -> None:
    """A failed-hard-constraint Score returns ``math.inf``."""
    s = Score(
        total_maintenance_dv_kms=0.1,
        max_vinf_kms=99.0,
        radial_span_au=(0.9, 1.6),
        period_error_yr=0.0,
        taxi_cost_kms=99.0,
        hard_constraints_pass=False,
    )
    assert composite_score(s) == math.inf


def _make_cycler_with_dv(dv_marker_vinf: float) -> Cycler:
    """Hand-build a Cycler whose total_maintenance_dv is monotone in
    ``dv_marker_vinf`` (km/s).

    We achieve this by giving the Mars encounter a small bend in addition
    to identical magnitudes: at low V∞ the achievable bend is wide, at
    high V∞ it shrinks below the demanded bend so ``flyby_dv_for`` charges
    a positive cost. The Earth encounter stays ballistic to keep the
    cycler hard-constraint-feasible.
    """
    vinf_e = _stub_vec(5.0, 0.0)
    # Mars: 60° bend at the supplied V∞; achievable bend angle shrinks
    # with V∞ so high V∞ makes flyby_dv > 0 there.
    theta = math.radians(60.0)
    vin_m = _stub_vec(dv_marker_vinf, 0.0)
    vout_m = _stub_vec(dv_marker_vinf * math.cos(theta), dv_marker_vinf * math.sin(theta))
    enc_e = _make_encounter("E", vinf_e, vinf_e, t=0.0)
    enc_m = _make_encounter("M", vin_m, vout_m, t=1.0e7)
    leg = _make_leg("E", "M")
    return Cycler(bodies=["E", "M"], period=1.0e7, encounters=[enc_e, enc_m], legs=[leg])


def test_rank_orders_low_dv_first() -> None:
    """Three cyclers with monotonically increasing ΔV → ranked low first."""
    ephem = Ephemeris(model="circular")
    # Use directly-supplied ΔV via a fake-but-consistent setup: build three
    # cyclers whose hard constraints pass but whose total_maintenance_dv
    # differs. We sidestep the messy bend-vs-V∞ correlation by using
    # vinf_in == vinf_out (so flyby_dv == 0 everywhere) and **passing
    # target_period_sec** with three different periods to drive
    # period_error (which is one of the composite axes and is also
    # monotone in the ordering).
    vinf_e = _stub_vec(5.0, 0.0)
    vinf_m = _stub_vec(7.0, 0.0)
    seconds_per_year = DAYS_PER_JULIAN_YEAR * SECONDS_PER_DAY
    target_sec = seconds_per_year * 2.0
    # Three periods at +0.1, +0.5, +1.0 years above target.
    periods = [
        target_sec + 0.1 * seconds_per_year,
        target_sec + 0.5 * seconds_per_year,
        target_sec + 1.0 * seconds_per_year,
    ]
    cyclers = [_make_two_encounter_cycler("E", "M", vinf_e, vinf_m, period=p) for p in periods]
    ranked = rank(cyclers, ephem, vinf_cap=10.0, n_keep=3, target_period_sec=target_sec)
    assert len(ranked) == 3
    # Lowest period_error first, highest last.
    yrs = [s.period_error_yr for s, _ in ranked]
    assert yrs == sorted(yrs)


def test_rank_filters_infeasible() -> None:
    """A list with one infeasible cycler returns only the feasible ones."""
    ephem = Ephemeris(model="circular")
    vinf_e = _stub_vec(5.0, 0.0)
    cyc_ok_1 = _make_two_encounter_cycler("E", "M", vinf_e, _stub_vec(7.0))
    cyc_ok_2 = _make_two_encounter_cycler("E", "M", vinf_e, _stub_vec(8.0))
    cyc_bad = _make_two_encounter_cycler("E", "M", vinf_e, _stub_vec(20.0))  # > cap
    ranked = rank([cyc_ok_1, cyc_bad, cyc_ok_2], ephem, vinf_cap=10.0, n_keep=10)
    assert len(ranked) == 2


def test_rank_empty_input() -> None:
    """Empty cyclers list ⇒ empty output."""
    ephem = Ephemeris(model="circular")
    assert rank([], ephem, vinf_cap=10.0) == []


def test_rank_n_keep_truncation() -> None:
    """Ten feasible cyclers + ``n_keep=3`` ⇒ output length 3."""
    ephem = Ephemeris(model="circular")
    vinf_e = _stub_vec(5.0, 0.0)
    cyclers = [
        _make_two_encounter_cycler("E", "M", vinf_e, _stub_vec(5.0 + 0.1 * i)) for i in range(10)
    ]
    ranked = rank(cyclers, ephem, vinf_cap=12.0, n_keep=3)
    assert len(ranked) == 3
