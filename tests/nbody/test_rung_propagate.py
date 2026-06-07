"""N-body Phase B: full-period propagation + terminal closure error (plan Phase B).

NON-GOLDEN: asserts a finite, recorded closure error for a held SILVER candidate;
the candidate V_inf is OUR computation, never an EXPECTED.
"""

from __future__ import annotations

import math

import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.data.review_queue import ReviewQueueEntry  # noqa: E402
from cyclerfinder.nbody.rung import RungArc, propagate_one_period  # noqa: E402


@pytest.mark.slow
def test_propagate_one_period_reports_finite_closure(
    silver_fixture: ReviewQueueEntry,
) -> None:
    arc = propagate_one_period(silver_fixture, Ephemeris("astropy"), accuracy=1e-10)
    assert isinstance(arc, RungArc)
    assert math.isfinite(arc.terminal_closure_km)
    assert arc.terminal_closure_km >= 0.0
