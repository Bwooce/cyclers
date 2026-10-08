"""#1025: the literature gate (v3-A1) on the 42 clean #973 members.

Signatures are built exactly as scripts/litcheck_942_943_scope.py builds them for gc-1 and ev-A
(its ``sig_of``: time-ordered encounter sequence, period k, V_inf per encounter, topology
{"repeated-moon"}, working_bodies from the demanded turns at 0.05 deg, return types from the cycle
key and the leg flight times). The offline gate runs once per member.

In the same run, as checks of the gate state (expected values from data/942_943_litcheck_scope.json,
the v3-A1 record): the GanCal#5 and Hollister 1H controls (published) and the five #942/#943
candidates (gc-1, gc-2, ev-C inconclusive; ev-A, ev-B not-found).

Pre-registration: docs/notes/2026-10-08-1025-literature-step.md sec. 2.

Usage: uv run python scripts/litcheck_1025_973_members.py [--labels-only] --out FILE
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SCOPE = _load("litcheck_942_943_scope", REPO / "scripts" / "litcheck_942_943_scope.py")

#: The clean members (#973 note secs. 10-14): gauntlet file and candidate indices.
MEMBERS = {
    "gc4": ("data/973_gc/k4_gauntlet.json", [1, 2, 6, 8, 10, 19, 25, 28, 38, 41]),
    "gc5": ("data/973_gc/k5_gauntlet.json", [0, 3, 4, 5, 6, 7]),
    "gc6": ("data/973_gc/k6_gauntlet.json", [0, 1, 2, 3]),
    "ev4": ("data/973_ev/k4_gauntlet.json", [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
    "ev5": ("data/973_ev/k5_gauntlet.json", list(range(11))),
}

#: Gate-state checks: (name, gauntlet, V_inf pick, key or None, expected status).
CHECKS = [
    ("GanCal#5", "data/943_cell_gc_gauntlet.json", {"Ganymede": 3.238}, None, "published"),
    (
        "Hollister 1H",
        "data/942_cell_ev_gauntlet.json",
        {"E": 2.994, "V": 3.19},
        "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s",
        "published",
    ),
    ("gc-1", "data/943_cell_gc_gauntlet.json", {"Ganymede": 2.397}, None, "inconclusive"),
    ("gc-2", "data/943_cell_gc_gauntlet.json", {"Ganymede": 3.617}, None, "inconclusive"),
    ("ev-A", "data/942_cell_ev_gauntlet.json", {"E": 4.893, "V": 10.364}, None, "not-found"),
    ("ev-B", "data/942_cell_ev_gauntlet.json", {"E": 8.012, "V": 10.932}, None, "not-found"),
    ("ev-C", "data/942_cell_ev_gauntlet.json", {"E": 9.075, "V": 13.166}, None, "inconclusive"),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels-only", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    out: dict[str, list[dict[str, Any]]] = {"checks": [], "members": []}
    if not args.labels_only:
        for name, path, pick, key, expect in CHECKS:
            rec = SCOPE.run(name, SCOPE.sig_of(*SCOPE.from_gauntlet(path, pick, key)))
            rec["expected"] = expect
            rec["as_expected"] = rec["status"] == expect
            out["checks"].append(rec)
    screen4 = json.loads((REPO / "data/973_gc/k4_screen.json").read_text())
    for tag, (path, idx) in MEMBERS.items():
        g = json.loads((REPO / path).read_text())
        for i in idx:
            c = g["candidates"][i]
            # gc k = 4 predates the in-gauntlet screen; its 2.6/2.7 statuses are in k4_screen.json
            status = c.get("screen_status") or screen4["candidates"][i]["screen_status"]
            assert status == "pass", (tag, i, status)
            sig = SCOPE.sig_of(g["cell"], c["key"], c["x_days"])
            name = f"{tag}-{i}"
            if args.labels_only:
                rec = {
                    "name": name,
                    "sequence": list(sig.sequence),
                    "working_bodies": sig.working_bodies,
                    "return_types": sorted(sig.return_types or ()),
                }
                print(rec)
            else:
                rec = SCOPE.run(name, sig)
            rec |= {"key": c["key"], "period_k": sig.period_k}
            out["members"].append(rec)
    args.out.write_text(json.dumps(out, indent=1, default=str))
    if not args.labels_only:
        bad = [c["name"] for c in out["checks"] if not c["as_expected"]]
        print("CHECKS NOT AS EXPECTED:", bad or "none")


if __name__ == "__main__":
    main()
