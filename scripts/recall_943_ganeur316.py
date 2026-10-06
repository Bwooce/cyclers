"""#943 control: blind recall of Russell & Strange GanEur#316 in cell ge (both moons massive).

GanEur#316 (AAS 07-118 Tables 3 and 5; R-S 2009 Tables 5 and 6) is
g(1.31322,472.76044,U) h(1.5,540.0,L,-3.98557) g(1.31322,472.76044,U) G(2.77217,1357.97970,U):
two 1-rev generic Ganymede returns around a 3-pi half-rev Ganymede return, then the
Ganymede -> Europa -> Ganymede transit. It is the published ballistic 10-cycle patched-conic
ephemeris cycler of Fig. 10(a), and the real-ephemeris control for the ge cell (results note
sec. 6.27).

Structures tried (all of them, blind): Ganymede block (LG>G/1 low|high, HGanymede/3,1,peri|apo,
LG>G/1 low|high), then LG>E/n1 and LE>G/n2 with n1, n2 in 0..2 (both branches when n > 0),
k = 7 (49.36 d). Production seeds (n_phase 36, n_refine 40) with n_split 6 (four Lambert legs).

Expected (Table 3): V_inf G/E 3.20/3.81 km/s, period 49.4 d, minimum altitude at Ganymede
1,447 km, distance to Jupiter 592,969-1,496,829 km, transits G->E 7.60 d and E->G 12.23 d,
Europa turn 0 (Europa is the massless target in R-S's ideal model).

Usage: uv run python scripts/recall_943_ganeur316.py --out FILE [--limit N]
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import time
from pathlib import Path
from typing import Any

from cyclerfinder.search.two_working_body import Cycle, HalfRevLeg, LambertLeg
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


def candidate_cycles(period_s: float) -> list[Cycle]:
    g, e = "Ganymede", "Europa"
    transfers = [(0, "single"), (1, "low"), (1, "high"), (2, "low"), (2, "high")]
    out = []
    for b1, peri, b2 in itertools.product(("low", "high"), (True, False), ("low", "high")):
        block = (LambertLeg(g, g, 1, b1), HalfRevLeg(g, 3, 1, peri, 0), LambertLeg(g, g, 1, b2))
        for (n1, br1), (n2, br2) in itertools.product(transfers, transfers):
            legs = (*block, LambertLeg(g, e, n1, br1), LambertLeg(e, g, n2, br2))
            out.append(Cycle(legs, period_s))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=0, help="first N structures only (timing)")
    args = ap.parse_args()
    system, a, b = ENUM.cell_system("ge")
    syn = system.synodic_s(a, b)
    cycles = candidate_cycles(7 * syn)
    if args.limit:
        cycles = cycles[: args.limit]
    t0 = time.time()
    recs = []
    for i, cyc in enumerate(cycles):
        key = ENUM.cycle_key(cyc, 7)
        ts = time.time()
        try:
            zeros = solve_structure(
                system, cyc, phase_period_s=syn, n_phase=36, n_split=6, n_refine=40
            )
        except Exception as exc:  # recorded, never counted as "absent"
            recs.append({"key": key, "error": repr(exc)})
            print(f"[{i + 1}/{len(cycles)}] {key}: ERROR {exc!r}", flush=True)
            continue
        for z in zeros:
            ass = assess(system, z)
            rec = {
                "key": key,
                "x_days": [float(v) / DAY for v in z.x],
                "residual_kms": z.residual_kms,
                "status": ass.status,
                "vinf_kms": ass.vinf_kms,
                "r_min_km": ass.r_min_km,
                "r_max_km": ass.r_max_km,
                "max_encounter_miss_km": ass.max_encounter_miss_km,
                "min_required_alt_km": ass.report.gate.min_required_alt_km if ass.report else None,
                "flybys": flyby_table(system, ass),
            }
            recs.append(rec)
        el = time.time() - t0
        print(
            f"{time.strftime('%H:%M:%S')} [{i + 1}/{len(cycles)}] {key}: {len(zeros)} zeros "
            f"({time.time() - ts:.1f}s) eta {el / (i + 1) * (len(cycles) - i - 1) / 60:.1f} min",
            flush=True,
        )
    args.out.write_text(json.dumps(recs, indent=1, default=float))
    hits = [
        r
        for r in recs
        if "vinf_kms" in r
        and abs(r["vinf_kms"].get("Ganymede", 0) - 3.20) <= 0.05
        and abs(r["vinf_kms"].get("Europa", 0) - 3.81) <= 0.05
    ]
    print(f"DONE {len(cycles)} structures, {len(recs)} records, {len(hits)} near 3.20/3.81")
    for r in hits:
        print(json.dumps({k: r[k] for k in ("key", "status", "vinf_kms", "x_days")}, default=float))


if __name__ == "__main__":
    main()
