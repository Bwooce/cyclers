"""#942/#943: re-report r_min/r_max of the reported candidates after the leg_extent fix.

Before the fix the extent sampled the Lambert legs only (results note sec. 6.30). Each candidate is
taken from its cell's gauntlet file (matched by key and V_inf), reassessed with the current code,
and written with the old and new values to data/942_943_extent_rereport.json.

Usage: uv run python scripts/extent_942_943_rereport.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search.two_working_body_enum import Zero, assess

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
AU = 149_597_870.7

#: name: (cell, gauntlet file, key, V_inf that identifies the row)
CANDIDATES: dict[str, tuple[str, str, str, dict[str, float]]] = {
    "gc-1": (
        "gc",
        "data/943_cell_gc_gauntlet.json",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|RCallisto/1:1|LCallisto>Ganymede/0s",
        {"Ganymede": 2.397},
    ),
    "gc-2": (
        "gc",
        "data/943_cell_gc_gauntlet.json",
        "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s",
        {"Ganymede": 3.617},
    ),
    "ev-A": (
        "ev",
        "data/942_cell_ev_gauntlet.json",
        "k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s",
        {"E": 4.89, "V": 10.36},
    ),
    "ev-B": (
        "ev",
        "data/942_cell_ev_gauntlet.json",
        "k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s",
        {"E": 8.01, "V": 10.93},
    ),
    "ev-C": (
        "ev",
        "data/942_cell_ev_gauntlet.json",
        "k2|LE>V/0s|LV>V/1h|LV>E/0s",
        {"E": 9.07, "V": 13.17},
    ),
}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    out = {}
    for name, (cell, path, key, vt) in CANDIDATES.items():
        g = json.loads((REPO / path).read_text())
        system, a, b = enum.cell_system(cell)
        rows = [
            c
            for c in g["candidates"]
            if key in c["keys"] and all(abs(c["vinf_kms"][q] - v) < 0.01 for q, v in vt.items())
        ]
        assert len(rows) == 1, (name, len(rows))
        c = rows[0]
        _, cyc = enum.parse_cycle_key(c["key"], system, a, b)
        ass = assess(system, Zero(cyc, np.asarray(c["x_days"]) * DAY, 0.0))
        out[name] = {
            "key": c["key"],
            "status": ass.status,
            "vinf_kms": ass.vinf_kms,
            "old_r_min_km": c["r_min_km"],
            "old_r_max_km": c["r_max_km"],
            "new_r_min_km": ass.r_min_km,
            "new_r_max_km": ass.r_max_km,
        }
        u, unit = (AU, "AU") if cell == "ev" else (1.0, "km")
        print(
            f"{name}: {ass.status} r_min {c['r_min_km'] / u:.5g} -> {ass.r_min_km / u:.5g}, "
            f"r_max {c['r_max_km'] / u:.5g} -> {ass.r_max_km / u:.5g} {unit}"
        )
    (REPO / "data" / "942_943_extent_rereport.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
