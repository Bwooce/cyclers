"""Geometric repeat L of a real-ephemeris chain (spec §14 V2-ballistic same-phase drift, #1007).

L is the smallest integer such that, for every encounter body, ``L P / T_body`` lies within
:data:`GEOMETRIC_REPEAT_TOL_REV` of an integer, with ``P`` the cycle period (``k`` real synodic
periods of the pair) and ``T`` the body's real mean sidereal period. The same-phase drift test
compares cycle ``i + L`` with cycle ``i`` and needs a chain of at least ``2 L`` cycles.

The tolerance is the OWNER's decision of 2026-10-11 ("Set 0.02 rev now"), for ALL systems. The
periods are those of the #1061 sweep (docs/notes/2026-10-11-1061-v2-negatives-sweep.md sec. 5.2).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Final

GEOMETRIC_REPEAT_TOL_REV: Final[float] = 0.02
"""Owner decision 2026-10-11: a body is back in phase when ``L P / T`` is within 0.02 rev of an
integer."""

#: Real mean sidereal periods (days), as in the #1061 sweep sec. 5.2.
REAL_PERIODS_D: Final[dict[str, float]] = {
    "E": 365.256363,
    "V": 224.7008,
    "M": 686.979586,
    "Io": 1.769138,
    "Europa": 3.551181,
    "Ganymede": 7.154553,
    "Callisto": 16.689018,
}


def synodic_period_d(a: str, b: str) -> float:
    """Real synodic period (days) of the pair ``a``-``b``."""
    return 1.0 / abs(1.0 / REAL_PERIODS_D[a] - 1.0 / REAL_PERIODS_D[b])


def geometric_repeat(
    period_d: float,
    bodies: Sequence[str],
    *,
    tol_rev: float = GEOMETRIC_REPEAT_TOL_REV,
    max_l: int = 1000,
) -> int:
    """The smallest ``L`` with ``L period_d / T_b`` within ``tol_rev`` of an integer for every body.

    Raises ``ValueError`` if none exists up to ``max_l``.
    """
    for big_l in range(1, max_l + 1):
        if all(
            abs(big_l * period_d / REAL_PERIODS_D[b] - round(big_l * period_d / REAL_PERIODS_D[b]))
            <= tol_rev
            for b in bodies
        ):
            return big_l
    raise ValueError(f"no geometric repeat up to L = {max_l} at tol {tol_rev} rev")


def repeat_for_pair(k: int, a: str, b: str, *, tol_rev: float = GEOMETRIC_REPEAT_TOL_REV) -> int:
    """``L`` for a cycle of ``k`` real synodic periods of the pair ``a``-``b``."""
    return geometric_repeat(k * synodic_period_d(a, b), (a, b), tol_rev=tol_rev)
