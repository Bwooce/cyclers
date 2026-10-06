"""#942/#943: re-run the assessment of stored zeros after a correctness fix.

The zeros (date vectors) come from the Lambert-only date residual and are unchanged;
the assessment (flyby directions, half-rev arrivals, gate, encounter self-consistency,
extents) is recomputed with the current code. Writes ``<out>/zeros.jsonl`` and copies
``structures.jsonl``; the original run directory is not modified.

Usage: uv run python scripts/reassess_942_zeros.py --cell vm DIR [DIR ...] --out OUTDIR
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from cyclerfinder.search.two_working_body_enum import Zero, assess, flyby_table

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0


def _load(name: str, path: Path):  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--key-contains", default="", help="reassess only zeros whose key contains this text"
    )
    args = ap.parse_args()
    system, a, b = ENUM.cell_system(args.cell)
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "structures.jsonl").open("w") as fs:
        for d in args.dirs:
            fs.write((d / "structures.jsonl").read_text())
    recs = [json.loads(line) for d in args.dirs for line in (d / "zeros.jsonl").open()]
    recs = [r for r in recs if args.key_contains in r["key"]]
    t0 = time.time()
    changed = 0
    with (args.out / "zeros.jsonl").open("w") as fz:
        for i, r in enumerate(recs):
            _, cycle = ENUM.parse_cycle_key(r["key"], system, a, b)
            x = np.asarray(r["x_days"]) * DAY
            try:
                ass = assess(system, Zero(cycle, x, r["residual_kms"]))
            except Exception as exc:  # recorded, never counted as a fail
                fz.write(json.dumps(r | {"status": "error", "error": repr(exc)}) + "\n")
                print(f"assessment ERROR {r['key']}: {exc!r}", flush=True)
                continue
            rep = ass.report
            new = r | {
                "status": ass.status,
                "status_hm_floor": ass.report_hm_floor.status if ass.report_hm_floor else None,
                "worst_ratio": rep.gate.worst_ratio if rep else None,
                "max_turn_deg": rep.max_turn_deg if rep else None,
                "min_required_alt_km": rep.gate.min_required_alt_km if rep else None,
                "max_encounter_miss_km": ass.max_encounter_miss_km,
                "r_min_km": ass.r_min_km,
                "r_max_km": ass.r_max_km,
                "vinf_kms": ass.vinf_kms,
                "flybys": flyby_table(system, ass),
                "reassessed": True,
            }
            if new["status"] != r["status"] or (
                (r.get("max_encounter_miss_km", 0.0) >= 1.0)
                != (new["max_encounter_miss_km"] >= 1.0)
            ):
                changed += 1
                print(
                    f"CHANGED {r['key']} x={[round(v, 4) for v in r['x_days']]}: "
                    f"{r['status']} -> {new['status']} "
                    f"(worst {r.get('worst_ratio')} -> {new['worst_ratio']})",
                    flush=True,
                )
            fz.write(json.dumps(new, default=float) + "\n")
            if (i + 1) % 200 == 0:
                el = time.time() - t0
                print(
                    f"{i + 1}/{len(recs)} changed={changed} "
                    f"eta {el / (i + 1) * (len(recs) - i - 1):.0f}s",
                    flush=True,
                )
    print(f"DONE {len(recs)} zeros, verdict changed on {changed}", flush=True)


if __name__ == "__main__":
    main()
