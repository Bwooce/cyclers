"""S1L1 real-ephemeris rediscovery — the idealised model could NOT host the
5.65/3.05 km/s anchors (see scripts/characterise_s1l1.py); this gate tests
whether the real ephemeris (Mars eccentricity) can, like the Aldrin cycler.
Anchors are sourced (spec §9); epoch and leg ToFs are computed."""

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.4  # sourced anchors; real-eph band


def _vinf_by_body(result: object) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:  # type: ignore[attr-defined]
        m = max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


@pytest.mark.slow
@pytest.mark.xfail(
    strict=False,
    reason=(
        "TOPOLOGY MISMATCH resolved by STAGE 1 multi-rev plumbing; flip to pass "
        "once the topology is confirmed (2026-06-03). STAGE 1 threads "
        "per_leg_revs / per_leg_branch end-to-end (optimise_cell_ephemeris -> "
        "optimise_maintenance_dv -> _build_chain) and adds a multi-rev ToF floor, "
        "so a 4-encounter E-M-E-E cell with a multi-rev Earth-to-Earth resonant "
        "interval is now reachable. This test still uses the 3-encounter E-M-E "
        "cell with a direct Mars->Earth return leg, which the multi-seed resolver "
        "confirms cannot host 3.05 km/s at Mars at any epoch/ToF: the S/L labels "
        "are Earth-to-Earth resonant intervals (see [[s1l1-nomenclature]]), not "
        "an E-M-E return leg. Re-modelling as outbound E->M plus the S1/L1 "
        "resonant intervals and pinning the winning per_leg_revs (see "
        "test_s1l1_idealised_multirev.py and scripts/characterise_s1l1.py) flips "
        "this gate. The SOURCED 5.65/3.05 anchors remain the only assertion "
        "targets."
    ),
)
def test_s1l1_real_ephemeris_rediscovers_anchors() -> None:
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    # Asymmetric family-appropriate seed: the sourced ~154 d E->M outbound +
    # the remainder of the 2-synodic period, so the phase-match resolves the
    # S1L1 launch epoch instead of a symmetric degenerate basin.
    t_syn_em_days = 779.9
    period_days = 2 * t_syn_em_days
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=8.15,
        priority_date_iso="2002-08-05",
        vinf_targets_kms={"E": VINF_E, "M": VINF_M},
        n_starts=5,
        seed=0,
        tof_seed_days=[154.0, period_days - 154.0],
    )
    assert result.constraints_satisfied
    v = _vinf_by_body(result)
    assert abs(v["E"] - VINF_E) < TOL
    assert abs(v["M"] - VINF_M) < TOL
