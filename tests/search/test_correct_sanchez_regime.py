"""M-ED Phase 5 Task 5.2: Sanchez-regime near-ballistic closure gate.

Promotes the prototype's demonstrated DE440 2030-2034 near-ballistic E-M-E-E
closure (scripts/correct_s1l1_twoarc.py; spec §5) to a regression: the N-arc
ballistic corrector closes >=1 chain in that regime with V_inf <= cap and all
flybys bend-feasible.

NON-GOLDEN (project memory feedback_golden_tests_sourced_only): the asserted V_inf
is OUR computation, NOT a sourced anchor. The gate is *a closed, feasible,
sub-cap chain exists in the near-ballistic regime* — regime, cap, feasibility
only. The closed family floors well above the (unverified) Sanchez/S1L1 anchors
(project memory project_s1l1_realeph_closure_blocker.md), which is why no sourced
V_inf is asserted here.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
_PERIOD_DAYS = (1.4612 + 2.8096) * 365.25  # 2-synodic, sourced Russell arcs


@pytest.mark.slow
def test_corrector_closes_feasible_chain_in_sanchez_regime() -> None:
    ephem = Ephemeris("astropy")
    # 2030-2034 near-ballistic launch window (prototype scan regime, spec §5).
    t0 = ((datetime(2030, 3, 22, tzinfo=UTC) + timedelta(days=-20)) - _J2000).total_seconds()
    result = ballistic_correct(
        sequence=("E", "M", "E", "E"),
        per_leg_revs=(0, 0, 1),
        per_leg_branch=("single", "single", "low"),
        t0_seed_sec=t0,
        tof_seed_days=(154.0, 379.0),  # E-E slack leg eliminated
        period_sec=_PERIOD_DAYS * 86400.0,
        ephem=ephem,
        vinf_cap=14.0,
        slack_leg=2,
    )
    # Regime + cap + feasibility only — NO sourced V_inf assertion.
    assert result.converged  # ballistic closure reached (residual < 0.1 km/s)
    assert result.bend_feasible  # every flyby's required turn fits its V_inf limit
    assert result.vinf_cap_ok  # max V_inf <= cap
    assert result.constraints_satisfied
