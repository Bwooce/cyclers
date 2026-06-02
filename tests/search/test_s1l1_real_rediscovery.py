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
        "aspirational: real-ephemeris S1L1 closure is the open goal. UPDATE "
        "2026-06-02: an asymmetric family-appropriate ToF seed (tof_seed_days "
        "= [154 d outbound, ~1406 d return]) was tried and does NOT close it "
        "either — the phase-match resolves no launch window (sentinel path), "
        "because a direct ~1406 d Mars->Earth leg cannot match the 3.05 km/s "
        "Mars anchor at any epoch. Root cause is therefore TOPOLOGY, not "
        "seeding: S1L1 is not a 3-encounter E-M-E cycler with a direct return "
        "leg; the S/L labels are Earth-to-Earth resonant intervals (see "
        "[[s1l1-nomenclature]]). Flips to passing once the entry is re-modelled "
        "as outbound E->M plus the S1/L1 Earth-to-Earth resonant intervals and "
        "optimised with that cell."
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
