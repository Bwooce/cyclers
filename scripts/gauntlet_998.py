"""#998 gauntlet for the Pluto small-moon cells: de-duplicate, re-fly (DOP853), filter, compare.

For each cell (ps pn pk ph) in data/998_pluto_smallmoons/<cell>/ it takes every gate-passing merged
group (scripts/analyse_942_enumeration.py rules), re-solves the extent of the best member, re-flies
every leg with DOP853 (scripts/gauntlet_942.py cross_check, imported with the literature module
stubbed: the literature step is deferred to #972 and nothing here calls it), and records:

  period, V_inf at Charon and at the target, Charon turn demanded / available and the ratio at the
  registry floor, required altitude, r_min / r_max of the whole cycle, the Pluto-surface check
  (r_min against Pluto's radius 1188.3 km), dwell inside Charon's sphere of influence, the SOI
  fraction of the re-fly miss, and literal-collision checks against the #320 Pluto rows and the
  catalogue's Pluto rows.

Usage: uv run python scripts/gauntlet_998.py --cell pn \
    --out data/998_pluto_smallmoons/pn_gauntlet.json
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.two_working_body import cycle_flybys, sphere_of_influence_km
from cyclerfinder.search.two_working_body_enum import leg_extent

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
CHARON_RADIUS_KM = 606.0  # registry
PLUTO_RADIUS_KM = 1188.3  # Nimmo et al. 2017 (mean radius); registry has no Pluto body radius


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_gauntlet_942() -> Any:
    stub = types.ModuleType("cyclerfinder.search.literature_check")
    stub.CandidateSignature = object  # type: ignore[attr-defined]
    stub.check_literature = None  # type: ignore[attr-defined]
    stub.offline_corpus_search = None  # type: ignore[attr-defined]
    real = sys.modules.get("cyclerfinder.search.literature_check")
    sys.modules["cyclerfinder.search.literature_check"] = stub
    try:
        return _load("gauntlet_942", REPO / "scripts" / "gauntlet_942.py")
    finally:
        if real is not None:
            sys.modules["cyclerfinder.search.literature_check"] = real
        else:
            sys.modules.pop("cyclerfinder.search.literature_check", None)


def scan320_rows() -> list[dict[str, Any]]:
    rows = []
    for line in (REPO / "data" / "scan_320_epoch_aware_pluto.jsonl").read_text().splitlines():
        d = json.loads(line)
        if not d.get("_meta") and "Charon" in d["sequence"]:
            rows.append(d)
    return rows


def catalogue_pluto_rows() -> list[str]:
    cat = yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())
    rows = cat if isinstance(cat, list) else cat.get("cyclers", [])
    return [
        str(r.get("id"))
        for r in rows
        if "Pluto" in (r.get("bodies") or []) or r.get("primary") == "Pluto"
    ]


def collide_320(moon: str, vinf_charon: float, vinf_moon: float, period_d: float) -> list[str]:
    """#320 Pluto rows with Charon and this moon, by tier. #320 rows are two-leg cycles that do
    not close exactly (residual 0.0026-0.044 km/s) and none passes its physical gate, so only a
    close match of both V_inf values and the total flight time is called LITERAL?: each of the
    Charon and moon V_inf within 0.02 km/s (the row's Charon value taken at its nearest Charon
    encounter, the moon value at its moon encounter) and the total flight time within 1 percent
    of this candidate's period. NEAR: Charon V_inf within 0.05 km/s and time within 2 percent."""
    out = []
    for r in scan320_rows():
        if moon not in r["sequence"]:
            continue
        tof = sum(r["tof_days"])
        vc = [
            v
            for body, v in zip(r["sequence"], r["vinf_per_encounter_kms"], strict=True)
            if body == "Charon"
        ]
        vm = [
            v
            for body, v in zip(r["sequence"], r["vinf_per_encounter_kms"], strict=True)
            if body == moon
        ]
        dc = min(abs(v - vinf_charon) for v in vc)
        dm = min(abs(v - vinf_moon) for v in vm)
        rel = abs(tof - period_d) / period_d
        label = f"{r['sequence']} n_rev {r['n_rev']} residual {r['residual_kms']:.4f}"
        if dc <= 0.02 and dm <= 0.02 and rel <= 0.01:
            out.append(f"#320 LITERAL? {label}")
        elif dc <= 0.05 and rel <= 0.02:
            out.append(f"#320 NEAR {label}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    drv = _load("run_998_enumerate", REPO / "scripts" / "run_998_enumerate.py")
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    ana = _load("analyse_942_enumeration", REPO / "scripts" / "analyse_942_enumeration.py")
    g942 = _load_gauntlet_942()
    moon = drv.CELLS[args.cell]
    system = drv.pluto_system(moon)
    a, b = "Charon", moon
    d = REPO / "data" / "998_pluto_smallmoons" / args.cell
    n_struct, recs, groups, passing = ana.collect([d])
    n_zero_err = sum(1 for r in recs if r.get("status") == "error")
    print(
        f"cell {args.cell}: structures {n_struct}, zeros {len(recs)}, physical {len(groups)}, "
        f"gate-passing {len(passing)}",
        flush=True,
    )
    cat_rows = catalogue_pluto_rows()
    results = []
    items = sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    if args.limit:
        items = items[: args.limit]
    for i, ((k, _seq), g) in enumerate(items):
        best = min(
            (m for m in g["members"] if ana.is_pass(m)),
            key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9,
        )
        kk, cycle = enum.parse_cycle_key(best["key"], system, a, b)
        x = np.asarray(best["x_days"]) * DAY
        fl = cycle_flybys(system, cycle, x)
        assert fl is not None
        r_min, r_max = leg_extent(system, cycle, x, flybys=fl)
        xc = g942.cross_check(system, cycle, x, fl)
        soi_charon = sphere_of_influence_km(system, "Charon")
        soi_moon = sphere_of_influence_km(system, moon)
        vc = best["vinf_kms"]["Charon"]
        vm = best["vinf_kms"][moon]
        period_d = kk * system.synodic_s(a, b) / DAY
        charon_fb = [f for f in best["flybys"] if f["body"] == "Charon"]
        rec = {
            "k": k,
            "key": best["key"],
            "keys": sorted({m["key"] for m in g["members"]}),
            "mirror_pair": g["mirror_pair"],
            "x_days": best["x_days"],
            "n_member_zeros": len(g["members"]),
            "period_days": period_d,
            "vinf_charon_kms": vc,
            "vinf_moon_kms": vm,
            "n_charon_flybys": len(charon_fb),
            "charon_flybys": charon_fb,
            "worst_ratio": best["worst_ratio"],
            "max_turn_deg": best["max_turn_deg"],
            "min_required_alt_km": best["min_required_alt_km"],
            "status_hm_floor": best["status_hm_floor"],
            "r_min_km": r_min,
            "r_max_km": r_max,
            "pluto_surface_clear": bool(r_min > PLUTO_RADIUS_KM),
            "max_required_rp_charon_km": max(f["required_alt_km"] for f in charon_fb)
            + CHARON_RADIUS_KM,
            "rp_within_charon_soi": bool(
                max(f["required_alt_km"] for f in charon_fb) + CHARON_RADIUS_KM <= soi_charon
            ),
            "r_min_over_pluto_radius": r_min / PLUTO_RADIUS_KM,
            "soi_charon_km": soi_charon,
            "soi_moon_km": soi_moon,
            "dwell_in_charon_soi_days": 2.0 * soi_charon / vc / DAY,
            "cross_check": xc,
            "soi_fraction": xc["max_arrival_miss_km"] / min(soi_charon, soi_moon),
            "collisions_320": collide_320(moon, vc, vm, period_d),
        }
        results.append(rec)
        if (i + 1) % 20 == 0 or i + 1 == len(items):
            print(f"[{i + 1}/{len(items)}]", flush=True)
    out = {
        "cell": args.cell,
        "moon": moon,
        "n_structures": n_struct,
        "n_zeros": len(recs),
        "n_physical": len(groups),
        "n_gate_passing": len(passing),
        "n_zero_assessment_errors": n_zero_err,
        "catalogue_pluto_rows": cat_rows,
        "candidates": results,
    }
    for r in results:
        r["strong"] = bool(r["pluto_surface_clear"] and r["rp_within_charon_soi"])
    clear = [r for r in results if r["pluto_surface_clear"]]
    args.out.write_text(json.dumps(out, indent=1, default=float))
    print(
        f"gate-passing {len(results)}; strong {sum(r['strong'] for r in results)}; "
        f"Pluto-surface clear {len(clear)}; "
        f"re-fly gate disagrees "
        f"{sum(1 for r in results if r['cross_check']['gate_status_integrated'] != 'pass')}",
        flush=True,
    )


if __name__ == "__main__":
    main()
