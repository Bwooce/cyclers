"""#942/#943 writeback rows: the spec sec. 14 V1 evidence, recomputed.

For each written-back two-working-body cycler, from its representative ideal-model zero
(the cell gauntlet files):
- CONSISTENCY gate: every Lambert leg re-solved with lamberthub izzo2015 + gooding1990 agrees with
  the in-house solver to < 1e-3 m/s (spec sec. 14 V1);
- re-propagation gate: every leg (Lambert and fixed) re-flown with scipy DOP853 from the
  solved/chosen departure state meets the arrival body to < 1 km, and the turn gate passes on the
  integrated vectors (scripts/gauntlet_942.py cross_check).

Expected thresholds are the spec's (1e-3 m/s; the 1-km encounter tolerance of the results note sec.
6.1), not values computed by this code.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.lambert import lambert_crosscheck
from cyclerfinder.search.two_working_body import LambertLeg, cycle_flybys, eval_lambert_legs

REPO = Path(__file__).resolve().parents[2]
DAY = 86400.0


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ROWS = {
    "ganymede-callisto-two-working-body-cycler-gc1-2026": (
        "gc",
        "data/943_cell_gc_gauntlet.json",
        {"Ganymede": 2.397},
    ),
    "ganymede-callisto-two-working-body-cycler-gc2-2026": (
        "gc",
        "data/943_cell_gc_gauntlet.json",
        {"Ganymede": 3.617},
    ),
    "earth-venus-venus-hosted-cycler-evc-2026": (
        "ev",
        "data/942_cell_ev_gauntlet.json",
        {"E": 9.075, "V": 13.166},
    ),
    "earth-venus-two-working-body-cycler-eva-2026": (
        "ev",
        "data/942_cell_ev_gauntlet.json",
        {"E": 4.893, "V": 10.364},
    ),
}


@pytest.mark.parametrize("row_id", sorted(ROWS))
def test_writeback_row_v1_evidence(row_id: str) -> None:
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    gauntlet = _load("gauntlet_942", REPO / "scripts" / "gauntlet_942.py")
    cell, path, pick = ROWS[row_id]
    g = json.loads((REPO / path).read_text())
    rows = [
        c for c in g["candidates"] if all(abs(c["vinf_kms"][q] - v) < 0.01 for q, v in pick.items())
    ]
    assert len(rows) == 1
    c = rows[0]
    system, a, b = enum.cell_system(cell)
    _, cyc = enum.parse_cycle_key(c["key"], system, a, b)
    x = np.asarray(c["x_days"]) * DAY
    legs = eval_lambert_legs(system, cyc, x)
    assert legs is not None
    for ev, li in zip(legs, cyc.lambert_index, strict=True):
        leg = cyc.legs[li]
        assert isinstance(leg, LambertLeg)
        r1, _ = system.state(leg.frm, ev.t_dep)
        r2, _ = system.state(leg.to, ev.t_arr)
        res = lambert_crosscheck(
            r1, r2, ev.t_arr - ev.t_dep, mu=system.mu, n_revs=leg.nrev, branch=leg.branch
        )
        assert float(res["max_diff_mps"]) < 1e-3
    fl = cycle_flybys(system, cyc, x)
    assert fl is not None
    xc = gauntlet.cross_check(system, cyc, x, fl)
    assert xc["max_arrival_miss_km"] < 1.0
    assert xc["gate_status_integrated"] == "pass"
