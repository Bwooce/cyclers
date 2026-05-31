"""Tests for :mod:`cyclerfinder.search.resonance`.

Covers three of the four M2 spec §9 gate anchors:

* Earth-Mars synodic period = 2.135 yr (tol 0.001 yr).
* Earth-Venus synodic period = 1.599 yr (tol 0.001 yr).
* Venus-Earth-Mars beat tuple ``(3, 4)`` (or ordering-equivalent) at
  ~6.406 yr (tol 0.01 yr).
"""

from __future__ import annotations

from itertools import pairwise

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.search.resonance import (
    beat_period_days,
    k_synodic_periods_days,
    multi_body_beat_days,
    synodic_period_days,
    synodic_period_years,
)

# ---------------------------------------------------------------------------
# Gate anchors (spec §9)
# ---------------------------------------------------------------------------


def test_em_synodic_2135yr() -> None:
    """Earth-Mars synodic = 2.135 yr (tol 0.001 yr)."""
    yr = synodic_period_years("E", "M")
    assert yr == pytest.approx(2.135, abs=0.001), f"E-M synodic = {yr:.5f} yr"


def test_ev_synodic_1599yr() -> None:
    """Earth-Venus synodic = 1.599 yr (tol 0.001 yr)."""
    yr = synodic_period_years("E", "V")
    assert yr == pytest.approx(1.599, abs=0.001), f"E-V synodic = {yr:.5f} yr"


def test_vem_beat_yields_3_4() -> None:
    """VEM beat search returns the ``3*E-M ~ 4*E-V`` tuple.

    The returned tuple ordering follows ``[b for b in bodies if b != ref]``
    with ``ref = bodies[1] = E``, giving non-reference ``[V, M]`` and the
    tuple ``(k_V, k_M) = (4, 3)``. The reverse ordering ``(3, 4)`` would
    apply if the convention were ``[M, V]``; the test accepts either, the
    period assertion below pins the physics.
    """
    tuples = multi_body_beat_days(["V", "E", "M"], k_max=6)
    assert tuples, "no VEM beat tuples within tolerance"
    assert (4, 3) in tuples or (3, 4) in tuples


def test_vem_beat_period_6_406yr() -> None:
    """The top VEM beat tuple corresponds to ~6.406 yr (tol 0.01 yr)."""
    tuples = multi_body_beat_days(["V", "E", "M"], k_max=6)
    assert tuples, "no VEM beat tuples within tolerance"
    top = tuples[0]
    period_yr = beat_period_days(["V", "E", "M"], top) / DAYS_PER_JULIAN_YEAR
    assert period_yr == pytest.approx(6.406, abs=0.01), (
        f"top VEM beat period = {period_yr:.4f} yr (tuple = {top})"
    )


# ---------------------------------------------------------------------------
# API contract
# ---------------------------------------------------------------------------


def test_synodic_symmetric() -> None:
    """``T_syn(a, b) == T_syn(b, a)`` for all pairs."""
    pairs = [("E", "M"), ("E", "V"), ("V", "M")]
    for a, b in pairs:
        assert synodic_period_days(a, b) == synodic_period_days(b, a)


def test_synodic_self_raises() -> None:
    """``T_syn(body, body)`` is undefined and must raise."""
    with pytest.raises(ValueError, match="distinct bodies"):
        synodic_period_days("E", "E")


def test_synodic_unknown_body_raises() -> None:
    """Unknown one-letter body code propagates as ``KeyError``."""
    with pytest.raises(KeyError):
        synodic_period_days("E", "X")


def test_k_synodic_monotone() -> None:
    """``k_synodic_periods_days("E","M", 5)`` strictly increasing."""
    periods = k_synodic_periods_days("E", "M", 5)
    assert len(periods) == 5
    assert all(b > a for a, b in pairwise(periods))


def test_k_synodic_k_max_zero_raises() -> None:
    """k_max < 1 is nonsensical."""
    with pytest.raises(ValueError, match="k_max"):
        k_synodic_periods_days("E", "M", 0)


def test_multi_body_beat_k_max_zero_raises() -> None:
    """k_max < 1 is nonsensical for the beat search too."""
    with pytest.raises(ValueError, match="k_max"):
        multi_body_beat_days(["V", "E", "M"], k_max=0)


def test_multi_body_beat_single_body_raises() -> None:
    """Need at least 2 bodies for a beat."""
    with pytest.raises(ValueError, match=">= 2 bodies"):
        multi_body_beat_days(["E"], k_max=6)


def test_beat_period_days_wrong_tuple_length_raises() -> None:
    """k_tuple must match the count of non-reference bodies."""
    with pytest.raises(ValueError, match="k_tuple length"):
        beat_period_days(["V", "E", "M"], (1,))  # need 2 entries
