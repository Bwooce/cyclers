"""M-ED Phase 2: arc ToF seed in days (plan Phase 2; spec §16.7.7)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.search.descriptor import arc_tof_seed_days


def test_generic_tof_years_to_days() -> None:
    # g(1.4612,...) -> 1.4612 yr (spec §16.7.7; NOT the 526.02 deg psi param).
    assert arc_tof_seed_days("generic", tof_years=1.4612, resonance=None) == pytest.approx(
        1.4612 * DAYS_PER_JULIAN_YEAR
    )


def test_full_rev_tof_from_resonance() -> None:
    # F(3:2,...) -> resonant E-E interval ~ N=2 Earth years (seed only).
    assert arc_tof_seed_days("full-rev", tof_years=None, resonance="3:2") == pytest.approx(
        2.0 * DAYS_PER_JULIAN_YEAR
    )
