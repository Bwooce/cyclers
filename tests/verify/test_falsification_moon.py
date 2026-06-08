"""Tier-1 Phase 7 Task 7.2 (slow): Axis-D falsification — a Jovian moon chain
solved with the WRONG central mu (mu_Sun instead of mu_Jupiter) on Jupiter-centred
states must be refuted, not spuriously accepted.

This is the "deliberately bogus" Axis-D guard (gauntlet.py): a physically
inconsistent model must fail the verdict. Solving Lambert about the Sun on
Jupiter-centred (km-scale) moon states is the wrong centre; the resulting chain
must NOT report a bend-feasible closed cycler. NON-GOLDEN (a falsification
predicate, not a sourced value)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct


@pytest.mark.slow
def test_wrong_central_mu_does_not_spuriously_close() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    r = ballistic_correct(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
        t0_seed_sec=0.6 * 86400.0,
        tof_seed_days=(4.0, 3.4),
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        vinf_cap=20.0,
        slack_leg=2,
        ephem=ephem,
        mu_central=MU_SUN_KM3_S2,  # WRONG centre for a Jovicentric chain
    )
    # Lambert about the Sun on Jupiter-centred states is physically inconsistent;
    # it must not report a bend-feasible closed chain (Axis-D refutation).
    assert not (r.converged and r.bend_feasible)
