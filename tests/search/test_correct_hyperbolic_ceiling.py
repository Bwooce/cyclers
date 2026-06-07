"""Solver-layer elliptic-periodicity ceiling flag (task #127, deliverable 1).

A synthetic over-ceiling case must raise ``hyperbolic_impossible`` on the result
object — a loud flag, never a filter. Legitimate high-energy magnitudes (Russell
reaches 20.3 km/s at Earth) must NOT trip it.
"""

from __future__ import annotations

import math

from cyclerfinder.core.constants import VINF_CEILING_KMS, vinf_ceiling_kms
from cyclerfinder.search.correct import BallisticClosureResult, _hyperbolic_impossible


def test_earth_ceiling_matches_derivation() -> None:
    # v_esc_sun(1 AU) + v_Earth ~= 42.12 + 29.79 = 71.91 km/s (brief anchor).
    assert math.isclose(vinf_ceiling_kms("E"), 71.9, abs_tol=0.1)
    assert math.isclose(VINF_CEILING_KMS["E"], 71.9, abs_tol=0.1)


def test_ceiling_ordering_inner_higher_than_outer() -> None:
    # Closer to the Sun -> faster escape + orbital speed -> higher ceiling.
    assert VINF_CEILING_KMS["Me"] > VINF_CEILING_KMS["V"] > VINF_CEILING_KMS["E"]
    assert VINF_CEILING_KMS["E"] > VINF_CEILING_KMS["M"] > VINF_CEILING_KMS["J"]


def test_high_energy_russell_magnitude_does_not_trip_flag() -> None:
    # Russell-Ocampo legitimately reaches 20.3 km/s at Earth -- far below 71.9.
    seq = ("E", "M", "E")
    vinf = (20.3, 14.4, 19.0)
    assert _hyperbolic_impossible(seq, vinf) is False


def test_synthetic_over_ceiling_trips_flag() -> None:
    # 55.32 km/s at Earth is below the ceiling (it is the maintenance-dv class,
    # not a vinf-ceiling breach); a genuinely over-ceiling Earth encounter is
    # > 71.9. Use 80 km/s at Earth -- impossible for a periodic heliocentric leg.
    seq = ("E", "M", "E")
    vinf = (80.0, 12.0, 14.0)
    assert _hyperbolic_impossible(seq, vinf) is True


def test_per_body_ceiling_is_what_trips() -> None:
    # 60 km/s is fine at Earth (< 71.9) but impossible at Mars (< 58.25).
    assert _hyperbolic_impossible(("E",), (60.0,)) is False
    assert _hyperbolic_impossible(("M",), (60.0,)) is True


def test_unknown_body_never_treated_as_breach() -> None:
    # A body absent from the ceiling table cannot be assessed -> never a breach.
    assert _hyperbolic_impossible(("X",), (999.0,)) is False


def test_result_dataclass_carries_flag_default_false() -> None:
    r = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(100.0,),
        max_residual_kms=0.01,
        vinf_per_encounter_kms=(5.0,),
        converged=True,
        bend_feasible=True,
    )
    assert r.hyperbolic_impossible is False
