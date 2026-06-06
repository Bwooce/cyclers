"""M-ED Phase 1: ballistic_correct closes the S1L1 two-arc chain on DE440.

NON-GOLDEN non-regression fixture (spec §5, project memory): the asserted V_inf
floor (Mars ~6.4 km/s) is OUR prior computation, pinned from the live prototype,
NOT a published anchor. This guards the SOLVER against regression; it is NOT a
rediscovery of any sourced number. See project memory
project_s1l1_realeph_closure_blocker.md for why this family floors at ~6.4.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct


@pytest.mark.slow
def test_ballistic_correct_closes_s1l1_two_arc() -> None:
    ephem = Ephemeris("astropy")
    period_days = (1.4612 + 2.8096) * 365.25
    t0_seed = (
        np.datetime64("2030-03-22T00:00:00") - np.datetime64("2000-01-01T12:00:00")
    ) / np.timedelta64(1, "s")
    r = ballistic_correct(
        sequence=("E", "M", "E", "E"),
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
        t0_seed_sec=float(t0_seed),
        tof_seed_days=(154.0, 379.0),  # the two free legs; E-E is slack
        period_sec=period_days * 86400.0,
        ephem=ephem,
        vinf_cap=9.0,
        slack_leg=2,
    )
    assert r.converged  # ballistic closure reached (residual < 0.1 km/s)
    # The closed S1L1 family floors Mars V_inf ~6.4 km/s (OUR computation, not
    # a sourced anchor). Assert the regime, not a sourced value.
    vinf_mars = r.vinf_per_encounter_kms[1]  # encounter index 1 = Mars
    assert vinf_mars > 5.5
