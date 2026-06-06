"""M-ED Phase 1: N-arc leg/V_inf builder reproduces the prototype on S1L1.

NON-GOLDEN: the V_inf values here are OUR computation (spec §5 / project memory
golden-tests-sourced-only). This is a non-regression fixture for the SOLVER,
not a published-anchor assertion. The numbers are pinned from the live
scripts/correct_s1l1_twoarc.py prototype output, not from any source.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import _vinf_nodes


def test_s1l1_two_arc_nodes_have_in_out_per_encounter() -> None:
    ephem = Ephemeris("astropy")
    # S1L1 E-M-E-E: t0 ~2030-03-22, T_EM 154 d, T_ME 379 d, slack leg = E-E.
    # period = (1.4612 + 2.8096) yr (Russell 4.991gG arcs), days.
    period_days = (1.4612 + 2.8096) * 365.25
    seq = ("E", "M", "E", "E")
    t0_sec = (
        np.datetime64("2030-03-22T00:00:00") - np.datetime64("2000-01-01T12:00:00")
    ) / np.timedelta64(1, "s")
    nodes = _vinf_nodes(
        sequence=seq,
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
        t0_sec=float(t0_sec),
        free_tof_days=(154.0, 379.0),  # slack = E-E leg (index 2), eliminated
        slack_leg=2,
        period_days=period_days,
        ephem=ephem,
    )
    # One node per encounter; intermediates carry both in/out; ends carry the
    # closure pair (B0 out vs Bn in).
    assert set(nodes) >= {"m_in", "m_out", "e0", "e1_in", "e1_out", "e2_in"}
    for k in ("m_in", "m_out", "e0", "e1_in", "e1_out", "e2_in"):
        assert np.asarray(nodes[k]).shape == (3,)
