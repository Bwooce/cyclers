"""#942 em-cell control: blind recall of published Russell-Ocampo Earth-Mars cyclers in cell em.

The em enumeration (k <= 3, <= 2 returns per block, 0-rev transfers) held no published row, so its
results had no in-run recall control (results note sec. 6.36). This run solves the published rows'
own structures with the em cell's model (Earth and Mars circular coplanar, Mars 1.875 yr, both
massive) and the production solver and assessment.

Rows (Russell 2004 Table 3.4 and Tables 3.5; McConaghy, Russell & Longuski 2005 Table 2 labels):
- 2.5.1.+0 = 2g(1-11/14, 11/14 rev, U) f(1:1) h(0.5, 0, U, +-15.081 deg) f(1:1): Earth block of
  three returns (full-rev 1:1, the half-year tilted-circle half-rev, full-rev 1:1), then the
  generic Earth-Earth return through Mars (11/14 rev: both pieces 0-rev). Expected: V_inf E 7.8 /
  M 9.9, E->M 94 d, turns 54, 54, 54 deg, turn ratio (max/required) 1.12.
- 2.3.1.+1 (Byrnes et al. "case 3") = 2g(2-11/14, 1-11/14 rev, U) f(1:1) h(0.5, 0, U, +-10.388):
  Earth block of two returns, then the generic return through Mars with one complete revolution
  (split 0+1 or 1+0). Expected: V_inf E 5.4 / M 5.3, E->M 143 d, turns 93, 93 deg, turn ratio 0.92.

Every cyclic order of the Earth block and both Lambert branches are tried, k = 2, production seeds.

Usage: uv run python scripts/recall_942_ro.py --out FILE
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import time
from pathlib import Path
from typing import Any

from cyclerfinder.search.two_working_body import Cycle, HalfRevLeg, LambertLeg, Leg, ResonantLeg
from cyclerfinder.search.two_working_body_enum import assess, flyby_table, solve_structure

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")

TRANSFERS = [(0, "single"), (1, "low"), (1, "high")]


def _orders(block: tuple[Leg, ...]) -> list[tuple[Leg, ...]]:
    return sorted(set(itertools.permutations(block)), key=str)


def candidate_cycles(period_s: float) -> dict[str, list[Cycle]]:
    f, h = ResonantLeg("E", 1, 1), HalfRevLeg("E", 1, 0, False, 0)
    out: dict[str, list[Cycle]] = {"2.5.1.+0": [], "2.3.1.+1": []}
    for blk in _orders((f, h, f)):
        out["2.5.1.+0"].append(
            Cycle(
                (*blk, LambertLeg("E", "M", 0, "single"), LambertLeg("M", "E", 0, "single")),
                period_s,
            )
        )
    for blk in _orders((f, h)):
        for (n1, b1), (n2, b2) in itertools.product(TRANSFERS, TRANSFERS):
            if n1 + n2 != 1:
                continue
            legs = (*blk, LambertLeg("E", "M", n1, b1), LambertLeg("M", "E", n2, b2))
            out["2.3.1.+1"].append(Cycle(legs, period_s))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    system, a, b = ENUM.cell_system("em")
    syn = system.synodic_s(a, b)
    recs = []
    t0 = time.time()
    for name, cycles in candidate_cycles(2 * syn).items():
        for cyc in cycles:
            key = ENUM.cycle_key(cyc, 2)
            zeros = solve_structure(
                system, cyc, phase_period_s=syn, n_phase=36, n_split=12, n_refine=40
            )
            for z in zeros:
                ass = assess(system, z)
                recs.append(
                    {
                        "row": name,
                        "key": key,
                        "x_days": [float(v) / DAY for v in z.x],
                        "status": ass.status,
                        "vinf_kms": ass.vinf_kms,
                        "worst_ratio": ass.report.gate.worst_ratio if ass.report else None,
                        "max_encounter_miss_km": ass.max_encounter_miss_km,
                        "flybys": flyby_table(system, ass),
                    }
                )
            print(f"{time.time() - t0:.0f}s {name} {key}: {len(zeros)} zeros", flush=True)
    args.out.write_text(json.dumps(recs, indent=1, default=float))
    for name, (ve, vm) in {"2.5.1.+0": (7.8, 9.9), "2.3.1.+1": (5.4, 5.3)}.items():
        hits = [
            r
            for r in recs
            if r["row"] == name
            and r["vinf_kms"]
            and abs(r["vinf_kms"].get("E", 0) - ve) <= 0.05
            and abs(r["vinf_kms"].get("M", 0) - vm) <= 0.05
        ]
        print(f"{name}: {len(hits)} zeros within 0.05 km/s of {ve}/{vm}")
        for r in hits:
            turns = [(f["body"], round(float(f["turn_deg"]), 1)) for f in r["flybys"]]
            print(f"  {r['key']} {r['status']} worst {r['worst_ratio']} turns {turns}")


if __name__ == "__main__":
    main()
