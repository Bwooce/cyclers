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
        "aspirational: real-ephemeris S1L1 closure is the open goal this plan "
        "targets; flips to passing once the optimiser reaches the 5.65/3.05 "
        "basin. Currently the equispaced cold-start seed resolves a launch "
        "epoch from a symmetric leg signature, which lands in a degenerate "
        "high-V_inf basin (achieved V_inf_E~21.7, V_inf_M~16.8 km/s vs the "
        "sourced 5.65/3.05 anchors) rather than the S1L1 family; reaching the "
        "basin needs a family-appropriate (asymmetric) ToF seed for the "
        "phase-match epoch resolution (see plan 'Open risk', line 287)."
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
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=8.15,
        priority_date_iso="2002-08-05",
        vinf_targets_kms={"E": VINF_E, "M": VINF_M},
        n_starts=5,
        seed=0,
    )
    assert result.constraints_satisfied
    v = _vinf_by_body(result)
    assert abs(v["E"] - VINF_E) < TOL
    assert abs(v["M"] - VINF_M) < TOL
