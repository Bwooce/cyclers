"""#1054: the family keys that the full F1 grid found and the #1000 grid did not (pre-registration:
docs/notes/2026-10-11-1054-seed-design-fix.md sec. 8), checked against the published record.

For each key, its cycler-class members in data/1054_seed_design/grid_f1.jsonl (and the pilot arm E):
- p:q and sense recomputed from N = w_E + T / (2 pi);
- (a) the #1047 literal table rule (scripts/run_1047_literal.py) against every published Earth-Moon
  p:q orbit held as numbers (data/1000_complement/literal_check.json);
- (b) the Restrepo-Russell 2018 database, the #997 match rule and tolerances
  (scripts/run_997_rr_crossmatch.py), for symmetric keys; an asymmetric key is outside RR by
  construction (RR holds axisymmetric orbits only), stated and not counted;
- (c) Franz-Russell 2022: maximum Moon distance against its cutoff of about 350,000 km;
- (d) the #997 families (data/997_lineage/F*_[pm].jsonl): same p:q, sense and windings, and C
  overlap;
- (e) the v3 literature gate (class level; it does not compare p:q).
Labels: known-reproduction (a literal member), known-class (a matched RR record or #997 family, or
class level only), unmatched. Never "new". Positive control first (a filter must pass a known
case): the #1000 grid's key "2:1 pro wE1 wM-1 sym" must give at least one RR match (#1000 note sec.
5.2: 1 of 4) and the #997 F3 family; otherwise stop.
Output: data/1054_seed_design/new_keys_check.json.
"""

from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    check_literature,
    offline_corpus_search,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_997_rr_crossmatch as rrx
import run_1047_literal as lit
import run_1054_seed_design as sd

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1054_seed_design" / "new_keys_check.json"
MU = 0.01215058439469525
FR_CUTOFF_KM = 350_000.0


CONTROL_KEY = "2:1 pro wE1 wM-1 sym"


def members(key: str) -> list[dict[str, Any]]:
    out = []
    files = (
        [sd.SRC / "grid.jsonl"]
        if key == CONTROL_KEY
        else [sd.OUT / "grid_f1.jsonl", sd.OUT / "arm_E.jsonl"]
    )
    for path in files:
        name = path.name
        for line in path.read_text().splitlines():
            r = json.loads(line)
            if r.get("cycler_class_candidate") and sd.family_key(r) == key:
                r["_file"] = name
                out.append(r)
    return out


def perp_crossings(s4: list[float], period: float) -> list[tuple[float, float]]:
    """Perpendicular x-axis crossings (x, ydot) over one period (|xdot| < 1e-6 at y = 0)."""

    def ev(t: float, y: np.ndarray, m: float) -> float:
        return float(y[1])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, period),
        np.array([s4[0], s4[1], 0.0, s4[2], s4[3], 0.0]),
        args=(MU,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=ev,
    )
    return [(float(y[0]), float(y[4])) for y in sol.y_events[0] if abs(y[3]) < 1e-6]


def lineage_families() -> dict[str, dict[str, Any]]:
    fams: dict[str, dict[str, Any]] = {}
    for f in sorted(glob.glob(str(REPO / "data" / "997_lineage" / "F*_[pm].jsonl"))):
        for line in Path(f).read_text().splitlines():
            m = json.loads(line)
            if "wind_E" not in m:
                continue
            n_in = m["wind_E"] + m["T"] / (2 * math.pi)
            key = (
                abs(round(n_in)),
                round(m["T"] / (2 * math.pi)),
                1 if n_in > 0 else -1,
                round(m["wind_E"]),
                round(m["wind_M"]),
            )
            fam = fams.setdefault(Path(f).name, {"keys": set(), "C": []})
            fam["keys"].add(key)
            fam["C"].append(m["C"])
    return fams


def main() -> None:
    preflight_search(
        task_no=1054,
        region_id="grid-f1-new-family-keys-literature",
        method=MethodCapability(
            genome="#1054 F1 grid new family keys vs published record",
            corrector="none (lookup)",
            capability_tags=frozenset({"literature", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    report = json.loads((sd.OUT / "grid_f1_report.json").read_text())
    keys = report["family_keys_not_in_old_grid"]
    rows = json.loads((REPO / "data" / "1000_complement" / "literal_check.json").read_text())
    rr, names = rrx.load_rr(rrx.DEFAULT_DB)
    j_all = rr[:, 2] + 3.0
    lineage = lineage_families()
    out: dict[str, Any] = {}
    for key in [CONTROL_KEY, *keys]:
        ms = members(key)
        info = []
        for r in ms:
            n_in = r["wind_E"] + r["T"] / (2 * math.pi)
            info.append(
                {
                    "file": r["_file"],
                    "C": r["C"],
                    "T": r["T"],
                    "pe": r["perigee_alt_km"],
                    "ps": r["periselene_alt_km"],
                    "pq": (abs(round(n_in)), round(r["T"] / (2 * math.pi))),
                    "sense": 1 if n_in > 0 else -1,
                    "max_moon_km": r["max_moon_km"],
                    "s_perigee": r["s_perigee"],
                    "symmetric": r["symmetric"],
                    "wind": (round(r["wind_E"]), round(r["wind_M"])),
                }
            )
        info.sort(key=lambda m: m["C"])
        fam = {
            "pq": info[0]["pq"],
            "sense": info[0]["sense"],
            "C": (info[0]["C"], info[-1]["C"]),
            "members": info,
        }
        # (a) literal table
        verdicts = [
            {"source": x["source"], "orbit": x["orbit"], "verdict": lit.verdict(x, fam)}
            for x in rows
        ]
        literal = [v for v in verdicts if v["verdict"].startswith("LITERAL")]
        same_pq = [v for v in verdicts if "same p:q" in v["verdict"]]
        # (b) RR
        sym = all(m["symmetric"] for m in info)
        rr_res = []
        if sym:
            for m in info:
                cr = perp_crossings(m["s_perigee"], m["T"])
                res = rrx.match(rr, names, m["C"], m["T"], cr)
                rr_res.append(
                    {
                        "C": m["C"],
                        "inside_J_range": bool(j_all.min() <= m["C"] <= j_all.max()),
                        "n_match": res["n_match"],
                        "best": res["best"][:1],
                    }
                )
        # (c) FR
        fr_excluded = all(m["max_moon_km"] > FR_CUTOFF_KM for m in info)
        # (d) #997 families
        k997 = (fam["pq"][0], fam["pq"][1], fam["sense"], info[0]["wind"][0], info[0]["wind"][1])
        f997 = [
            {"file": f, "C_range": (min(v["C"]), max(v["C"]))}
            for f, v in lineage.items()
            if k997 in v["keys"]
            and min(v["C"]) - 0.01 <= fam["C"][1]
            and fam["C"][0] <= max(v["C"]) + 0.01
        ]
        f997_any_c = [f for f, v in lineage.items() if k997 in v["keys"]]
        # (e) v3 gate
        ts = [m["T"] for m in info]
        sig = CandidateSignature(
            primary="Earth",
            sequence=("Moon",),
            period_k=fam["pq"][1],
            resonances=(f"{fam['pq'][0]}:{fam['pq'][1]}",),
            period_band_tu=(min(ts) * 0.999, max(ts) * 1.001),
            topology_label=frozenset({"resonant"}),
            working_bodies="one",
        )
        gate = check_literature(sig, search=offline_corpus_search)
        rr_hit = any(x["n_match"] > 0 for x in rr_res)
        if literal:
            label = "known-reproduction"
        elif rr_hit or f997:
            label = (
                "known-class (matched: "
                + ("RR " if rr_hit else "")
                + ("#997 " + ", ".join(x["file"] for x in f997) if f997 else "")
                + ")"
            )
        elif gate.status == "published":
            label = "known-class (class level only: v3 gate; no matching orbit held)"
        else:
            label = "unmatched"
        out[key] = {
            "n_members": len(info),
            "computed_pq": fam["pq"],
            "sense": "pro" if fam["sense"] > 0 else "ret",
            "C_range": fam["C"],
            "T_range": (min(ts), max(ts)),
            "perigee_alt_km": (min(m["pe"] for m in info), max(m["pe"] for m in info)),
            "periselene_alt_km": (min(m["ps"] for m in info), max(m["ps"] for m in info)),
            "symmetric": sym,
            "a_literal_members": literal,
            "a_same_pq_rows": same_pq,
            "b_rr": rr_res if sym else "asymmetric: outside RR by construction",
            "c_fr_excluded_by_construction": fr_excluded,
            "max_moon_km": (
                min(m["max_moon_km"] for m in info),
                max(m["max_moon_km"] for m in info),
            ),
            "d_997_same_key_C_overlap": f997,
            "d_997_same_key_any_C": f997_any_c,
            "e_gate": {"status": gate.status, "citation": gate.citation},
            "label": label,
        }
        print(f"{key}: {len(info)} members; {label}", flush=True)
        if key == CONTROL_KEY:
            ok = rr_hit and any(x["file"].startswith("F3") for x in f997)
            out[key]["control"] = True
            out[key]["control_pass"] = ok
            print(f"CONTROL {'PASS' if ok else 'FAIL: stop'}", flush=True)
            if not ok:
                break
    OUT.write_text(json.dumps(out, indent=1, default=list) + "\n")


if __name__ == "__main__":
    main()
