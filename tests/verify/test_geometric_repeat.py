"""Geometric repeat L at the owner's 0.02-rev tolerance (#1007; owner decision 2026-10-11).

Expected values are those the #1061 sweep computed independently (by hand, not with this module):
docs/notes/2026-10-11-1061-v2-negatives-sweep.md secs. 5.2-5.3 (gc k = 3: 4; Earth-Venus k = 2: 5;
Earth-Venus k = 3: 49, the 0.0201-rev edge; Earth-Mars k = 1 / k = 3: 37 / 32, i.e. 74 / 64
cycles for a 2L chain).
"""

from __future__ import annotations

import pytest

from cyclerfinder.verify.geometric_repeat import (
    GEOMETRIC_REPEAT_TOL_REV,
    REAL_PERIODS_D,
    geometric_repeat,
    repeat_for_pair,
    synodic_period_d,
)


def test_tolerance_is_the_owner_value() -> None:
    assert GEOMETRIC_REPEAT_TOL_REV == 0.02


@pytest.mark.parametrize(
    ("k", "a", "b", "expected"),
    [
        (3, "Ganymede", "Callisto", 4),  # gc-1, gc-2
        (2, "E", "V", 5),  # ev-A, ev-C, H&M orbit 1
        (3, "E", "V", 49),  # ev-B: 0.0201 rev at L = 5, just outside the tolerance
        (1, "E", "M", 37),
        (3, "E", "M", 32),
    ],
)
def test_repeat_for_pair(k: int, a: str, b: str, expected: int) -> None:
    assert repeat_for_pair(k, a, b) == expected


def test_ev_b_edge_is_just_outside() -> None:
    p = 3 * synodic_period_d("E", "V")
    for body in ("E", "V"):
        x = 5 * p / REAL_PERIODS_D[body]
        assert 0.02 < abs(x - round(x)) < 0.021, (body, x)
    assert geometric_repeat(p, ("E", "V"), tol_rev=0.021) == 5
