"""Unit tests for the paper's zero-radius-SOI flyby-maneuver ΔV (Hernandez 2017 Eqs 3-5).

:func:`cyclerfinder.nbody.jovian.flyby_maneuver_dv` charges ΔV only for the V∞-magnitude
mismatch (the bend is ballistic). The hand-check below is a *forward construction*: pick
``r_p`` and the two V∞ magnitudes, compute the bend ``δ`` Eq 4 *produces* at that ``r_p``,
feed vectors with that bend back in, and assert the solver recovers ``r_p``, the altitude,
and the Eq-5 ΔV. The expected numbers are independent hand arithmetic, not code output.
"""

from __future__ import annotations

import math

import numpy as np

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.nbody.jovian import flyby_maneuver_dv


def _bend_for(r_p: float, vi: float, vo: float, mu: float) -> float:
    """Eq-4 bend the flyby delivers at periapsis radius ``r_p`` (rad)."""
    return math.asin(mu / (mu + r_p * vi * vi)) + math.asin(mu / (mu + r_p * vo * vo))


def test_hand_checked_unequal_magnitude_case() -> None:
    """Ganymede, r_p=3000 km, |V∞in|=8, |V∞out|=7 km/s → ΔV≈945.9 m/s, alt≈368.8 km.

    Hand arithmetic (μ=9887.834, R=2631.2):
      δ = asin(μ/(μ+3000·64)) + asin(μ/(μ+3000·49)) = 0.112064 rad
      v_p,in  = sqrt(64 + 2μ/3000) = 8.40190 km/s
      v_p,out = sqrt(49 + 2μ/3000) = 7.45600 km/s
      ΔV = 945.90 m/s ; altitude = 3000 - 2631.2 = 368.8 km
    """
    mu = SATELLITES["Ganymede"].mu_km3_s2
    r_p, vi, vo = 3000.0, 8.0, 7.0
    delta = _bend_for(r_p, vi, vo, mu)
    vinf_in = np.array([vi, 0.0, 0.0])
    vinf_out = vo * np.array([math.cos(delta), math.sin(delta), 0.0])

    dv_ms, alt_km, feasible = flyby_maneuver_dv(vinf_in, vinf_out, "Ganymede")

    assert feasible
    assert abs(alt_km - 368.8) < 1e-6  # 3000 - 2631.2
    assert abs(dv_ms - 945.90) < 0.5


def test_equal_magnitude_is_ballistic() -> None:
    """Equal V∞ magnitudes → ΔV exactly 0 at a feasible r_p (ballistic-cycler property).

    Forward: r_p=4000 km at Ganymede gives a bend of 0.094282 rad; with both magnitudes
    7.07 km/s, v_p,in == v_p,out, so ΔV is identically zero. alt = 4000 - 2631.2 = 1368.8.
    """
    mu = SATELLITES["Ganymede"].mu_km3_s2
    r_p, v = 4000.0, 7.07
    delta = _bend_for(r_p, v, v, mu)
    vinf_in = np.array([v, 0.0, 0.0])
    vinf_out = v * np.array([math.cos(delta), math.sin(delta), 0.0])

    dv_ms, alt_km, feasible = flyby_maneuver_dv(vinf_in, vinf_out, "Ganymede")

    assert feasible
    assert dv_ms < 1e-9
    assert abs(alt_km - 1368.8) < 1e-3


def test_tiny_bend_is_infeasible_above_window() -> None:
    """A near-zero bend pushes r_p above the altitude window → infeasible (no real flyby)."""
    v = 7.0
    delta = 1.0e-5
    vinf_in = np.array([v, 0.0, 0.0])
    vinf_out = v * np.array([math.cos(delta), math.sin(delta), 0.0])
    _dv, alt_km, feasible = flyby_maneuver_dv(vinf_in, vinf_out, "Ganymede")
    assert not feasible
    assert alt_km > 70000.0


def test_excessive_bend_is_infeasible_below_window() -> None:
    """A bend larger than the deepest pass can supply → r_p below the window → infeasible."""
    v = 7.0
    delta = 3.0  # ~172°, near the π ceiling
    vinf_in = np.array([v, 0.0, 0.0])
    vinf_out = v * np.array([math.cos(delta), math.sin(delta), 0.0])
    _dv, alt_km, feasible = flyby_maneuver_dv(vinf_in, vinf_out, "Io")
    assert not feasible
    assert alt_km < 25.0
