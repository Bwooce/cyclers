# tests/search/test_leveraging_leg.py
"""Phase-full VILM single-leg evaluator (plan 2026-06-09, Component 1)."""

from __future__ import annotations

import math

import pytest

from cyclerfinder.search import leveraging_leg as ll


def test_resonant_sma_canonical() -> None:
    assert ll.resonant_sma(1, 1) == pytest.approx(1.0)
    assert ll.resonant_sma(2, 1) == pytest.approx(2.0 ** (2.0 / 3.0))


def test_tisserand_vinf_roundtrip() -> None:
    assert ll.tisserand_vinf(a=1.0, e=0.0) == pytest.approx(0.0, abs=1e-12)
    # a=1.6, e=0.4: periapsis 0.96 <= 1 <= apoapsis 2.24, so it crosses the moon.
    a, e = 1.6, 0.4
    v = ll.tisserand_vinf(a=a, e=e)
    assert ll.eccentricity_from_vinf(a=a, vinf=v) == pytest.approx(e, abs=1e-9)


def test_eccentricity_infeasible_returns_nan() -> None:
    assert math.isnan(ll.eccentricity_from_vinf(a=1.0, vinf=5.0))


def test_endgame_leg_lowers_vinf_at_europa_with_leverage() -> None:
    # Interior endgame VILM at Europa lowering V∞ 2.5 -> 2.2 km/s on a 2:1 orbit.
    res = ll.evaluate_leveraging_leg(
        moon="Europa",
        n_moon_revs=2,
        m_sc_revs=1,
        vinf_in_kms=2.5,
        vinf_out_target_kms=2.2,
        exterior=False,
    )
    assert res.converged
    assert res.vinf_out_kms == pytest.approx(2.2, abs=0.05)
    # Leverage realised: the near-root burn is small (~0.3 km/s), NOT the
    # flip-orbit far root (~4.7 km/s). Guards against the root-selection bug.
    assert 0.0 < res.dv_dsm_kms < 1.0
    # The realised ΔV must be >= the Γ-quadrature analytic floor (cross-check).
    assert res.gamma_floor_ok


def test_leg_noncrossing_preburn_orbit_does_not_converge() -> None:
    # A 3:1 resonant orbit (a=2.08) at a low V∞ has periapsis > 1: it never
    # reaches the moon, so there is no encounter and no leg -> not converged.
    # (No fabricated leg; the crossing constraint is honoured.)
    res = ll.evaluate_leveraging_leg(
        moon="Europa",
        n_moon_revs=3,
        m_sc_revs=1,
        vinf_in_kms=1.0,
        vinf_out_target_kms=0.7,
        exterior=True,
    )
    assert not res.converged
