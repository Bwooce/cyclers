"""#1000 known-class gates (note sec. 0.5) for every cycler-class candidate of method (1) and of
the method (2) branches: the #997 families (symmetric orbits; C-interpolated crossing match as
in scripts/run_997_gate.py), the Restrepo-Russell 2018 Earth-Moon files (symmetric orbits; the
#997 sec. 1.8 rule), Franz-Russell 2022 (maximum Moon distance > 350,000 km), and the
asymmetric families' labels. Output: data/1000_complement/gate_summary.json.
"""

from __future__ import annotations

import glob
import itertools
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_997_rr_crossmatch as rrx
import run_1000_complement as base

REPO = Path(__file__).resolve().parent.parent
GRID = REPO / "data" / "1000_complement" / "shooting" / "grid.jsonl"
OUT = REPO / "data" / "1000_complement" / "gate_summary.json"


def chains997() -> dict[str, list[dict[str, Any]]]:
    out = {}
    for f in sorted(glob.glob(str(REPO / "data" / "997_lineage" / "F*_[pm].jsonl"))):
        out[Path(f).name] = [json.loads(line) for line in Path(f).read_text().splitlines()]
    return out


def match997(
    chains: dict[str, list[dict[str, Any]]], c: float, t: float, perp: list[Any]
) -> str | None:
    for name, ch in chains.items():
        for a, b in itertools.pairwise(ch):
            if (a["C"] - c) * (b["C"] - c) > 0 or a["C"] == b["C"]:
                continue
            w = (c - a["C"]) / (b["C"] - a["C"])
            x0 = a["x0"] + w * (b["x0"] - a["x0"])
            yd = a["ydot0"] + w * (b["ydot0"] - a["ydot0"])
            tt = a["T"] + w * (b["T"] - a["T"])
            if abs(tt - t) / t > 1e-3:
                continue
            if any(abs(p[0] - x0) < 1e-3 and abs(p[1] - yd) < 1e-3 for p in perp):
                return name
    return None


def main() -> None:
    preflight_search(
        task_no=1000,
        region_id="em-complement-known-class-gate",
        method=MethodCapability(
            genome="#1000 candidates vs #997 families, RR 2018, FR 2022",
            corrector="none (lookup)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    rr, names = rrx.load_rr(rrx.DEFAULT_DB)
    j_all = rr[:, 2] + 3.0
    chains = chains997()
    recs = [json.loads(line) for line in GRID.read_text().splitlines()]
    cands = [r for r in recs if r["status"] == "converged" and r["cycler_class_candidate"]]
    rows = []
    for r in cands:
        key = (
            f"{r['p']}_{r['q']}_{'pro' if r['sense'] > 0 else 'ret'}"
            f"_wE{round(r['wind_E'])}_wM{round(r['wind_M'])}"
        )
        rec: dict[str, Any] = {
            "family_key": key,
            "symmetric": r["symmetric"],
            "C": r["C"],
            "T": r["T"],
            "perigee_alt_km": r["perigee_alt_km"],
            "periselene_alt_km": r["periselene_alt_km"],
            "max_moon_km": r["max_moon_km"],
            "fr_excluded": bool(r["max_moon_km"] > 350000.0),
            "in_rr_J_range": bool(j_all.min() <= r["C"] <= j_all.max()),
        }
        if r["symmetric"]:
            s = np.array(r["s_perigee"])
            # perpendicular crossings from the classification run
            sol_perp = base.classify(s, r["T"])["perp_crossings"]
            perp = [(p[1], p[2]) for p in sol_perp]
            rec["match_997"] = match997(chains, r["C"], r["T"], perp)
            m = rrx.match(rr, names, r["C"], r["T"], perp)
            rec["rr_matches"] = m["n_match"]
        else:
            rec["rr_matches"] = None
            rec["match_997"] = None
        rows.append(rec)
    fams: dict[str, dict[str, Any]] = {}
    for x in rows:
        f = fams.setdefault(
            x["family_key"] + ("_sym" if x["symmetric"] else "_asym"),
            {
                "n": 0,
                "C": [],
                "match_997": set(),
                "rr_hits": 0,
                "in_rr_J": 0,
                "fr_excluded": 0,
                "pe": [],
                "ps": [],
            },
        )
        f["n"] += 1
        f["C"].append(x["C"])
        if x["match_997"]:
            f["match_997"].add(x["match_997"])
        f["rr_hits"] += 1 if (x["rr_matches"] or 0) > 0 else 0
        f["in_rr_J"] += 1 if x["in_rr_J_range"] else 0
        f["fr_excluded"] += 1 if x["fr_excluded"] else 0
        f["pe"].append(x["perigee_alt_km"])
        f["ps"].append(x["periselene_alt_km"])
    summary = {}
    for k, f in sorted(fams.items()):
        summary[k] = {
            "n_candidates": f["n"],
            "C_range": [min(f["C"]), max(f["C"])],
            "perigee_alt_km": [min(f["pe"]), max(f["pe"])],
            "periselene_alt_km": [min(f["ps"]), max(f["ps"])],
            "matched_997_families": sorted(f["match_997"]),
            "n_rr_matched": f["rr_hits"],
            "n_in_rr_J_range": f["in_rr_J"],
            "n_fr_excluded": f["fr_excluded"],
        }
        print(k, summary[k], flush=True)
    OUT.write_text(json.dumps({"families": summary, "candidates": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
