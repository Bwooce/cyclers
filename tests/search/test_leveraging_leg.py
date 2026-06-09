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
