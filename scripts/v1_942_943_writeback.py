"""#942/#943 writeback V1 evidence: lamberthub (Izzo 2015, Gooding 1990) agreement on every Lambert
leg of a candidate's representative ideal-model zero, plus the independent DOP853 re-fly of every
leg (Lambert and fixed) from the gauntlet's cross_check, for spec sec. 14 V1.

Usage: uv run python scripts/v1_942_943_writeback.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.lambert import lambert_crosscheck
from cyclerfinder.search.two_working_body import LambertLeg, eval_lambert_legs

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")

#: name: (cell, gauntlet file, V_inf that identifies the row)
CANDIDATES = {
    "gc-1": ("gc", "data/943_cell_gc_gauntlet.json", {"Ganymede": 2.397}),
    "gc-2": ("gc", "data/943_cell_gc_gauntlet.json", {"Ganymede": 3.617}),
    "ev-A": ("ev", "data/942_cell_ev_gauntlet.json", {"E": 4.893, "V": 10.364}),
    "ev-C": ("ev", "data/942_cell_ev_gauntlet.json", {"E": 9.075, "V": 13.166}),
}


def main() -> None:
    out = {}
    for name, (cell, path, vt) in CANDIDATES.items():
        g = json.loads((REPO / path).read_text())
        rows = [
            c
            for c in g["candidates"]
            if all(abs(c["vinf_kms"][q] - v) < 0.01 for q, v in vt.items())
        ]
        assert len(rows) == 1, (name, len(rows))
        c = rows[0]
        system, a, b = ENUM.cell_system(cell)
        _, cyc = ENUM.parse_cycle_key(c["key"], system, a, b)
        x = np.asarray(c["x_days"]) * DAY
        legs = eval_lambert_legs(system, cyc, x)
        assert legs is not None
        diffs = []
        for ev, li in zip(legs, cyc.lambert_index, strict=True):
            leg = cyc.legs[li]
            assert isinstance(leg, LambertLeg)
            r1, _ = system.state(leg.frm, ev.t_dep)
            r2, _ = system.state(leg.to, ev.t_arr)
            res = lambert_crosscheck(
                r1, r2, ev.t_arr - ev.t_dep, mu=system.mu, n_revs=leg.nrev, branch=leg.branch
            )
            diffs.append(float(res["max_diff_mps"]))
        out[name] = {
            "key": c["key"],
            "x_days": c["x_days"],
            "lamberthub_max_diff_mps_per_leg": diffs,
            "dop853_refly_max_miss_km": c["cross_check"]["max_arrival_miss_km"],
            "dop853_max_vinf_vector_error_kms": c["cross_check"]["max_vinf_vector_error_kms"],
            "gate_on_integrated_vectors": c["cross_check"]["gate_status_integrated"],
        }
        print(name, out[name])
    (REPO / "data" / "942_943_writeback_v1.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
