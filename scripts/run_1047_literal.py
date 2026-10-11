"""#1047: explicit row-by-row literal collision table for the #1000 symmetric cycler-class families
and the #1048 3:1 / 7:2 families, against every published Earth-Moon p:q orbit held as numbers
(the rows already computed by scripts/run_1000_literal.py in
data/1000_complement/literal_check.json: Casoliva 2010 Table 3, Vaquero 2013, Newton 1959,
Casoliva 2008 seeds, Liang 2017, Liang 2020, Hoelker-Winston 1968, Genova-Aldrin 2015).

Families: the cycler-class candidates of data/1000_complement/shooting/grid.jsonl (#1000) and
gaps.jsonl (#1048), grouped by (skeleton p:q, inertial sense, winding about the Earth and about the
Moon) and split by symmetry. For each family, its p:q and sense are recomputed from the converged
orbits as for the published rows: N = w_E + T / (2 pi), p:q = |N| : round(T / 2 pi), sense = sign N.

Match rule (as #1000 sec. 5.1b): same p:q and sense; the published C inside the family's computed C
range +- 0.01; at that C (linear interpolation between the family's members) |dT|/T < 1e-3 and the
perigee and periselene altitudes within 10%; Earth-Moon mu. A published row at another mu, or with
no state, is class-level only. Output: data/1047_literature/literal_table.json and .md.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "data" / "1000_complement"
OUT = REPO / "data" / "1047_literature"


def families() -> dict[str, dict[str, Any]]:
    fams: dict[str, dict[str, Any]] = {}
    for fname, origin in (("grid.jsonl", "#1000"), ("gaps.jsonl", "#1048")):
        for line in (SRC / "shooting" / fname).read_text().splitlines():
            r = json.loads(line)
            if r.get("status") != "converged" or not r.get("cycler_class_candidate"):
                continue
            sense = "pro" if r["sense"] > 0 else "ret"
            key = (
                f"{r['p']}:{r['q']} {sense} wE{round(r['wind_E'])} wM{round(r['wind_M'])} "
                f"{'sym' if r['symmetric'] else 'asym'}"
            )
            n_in = r["wind_E"] + r["T"] / (2 * math.pi)
            f = fams.setdefault(key, {"origin": origin, "members": []})
            f["members"].append(
                {
                    "C": r["C"],
                    "T": r["T"],
                    "pe": r["perigee_alt_km"],
                    "ps": r["periselene_alt_km"],
                    "pq": (abs(round(n_in)), round(r["T"] / (2 * math.pi))),
                    "sense": 1 if n_in > 0 else -1,
                }
            )
    for f in fams.values():
        ms = sorted(f["members"], key=lambda m: m["C"])
        f["members"] = ms
        f["C"] = (ms[0]["C"], ms[-1]["C"])
        f["pq"] = ms[0]["pq"]
        f["sense"] = ms[0]["sense"]
        f["pq_consistent"] = all(m["pq"] == f["pq"] and m["sense"] == f["sense"] for m in ms)
    return fams


def interp(ms: list[dict[str, Any]], c: float) -> dict[str, float] | None:
    for a, b in itertools.pairwise(ms):
        if a["C"] <= c <= b["C"] and b["C"] > a["C"]:
            w = (c - a["C"]) / (b["C"] - a["C"])
            return {k: a[k] + w * (b[k] - a[k]) for k in ("T", "pe", "ps")}
    near = min(ms, key=lambda m: abs(m["C"] - c))
    if abs(near["C"] - c) <= 0.01:
        return {k: near[k] for k in ("T", "pe", "ps")}
    return None


def verdict(row: dict[str, Any], f: dict[str, Any]) -> str:
    pq = tuple(row.get("pq") or ())
    if not pq:
        return row.get("why_no_pq", "no p:q")
    if pq != tuple(f["pq"]):
        return f"p:q {pq[0]}:{pq[1]} differs"
    if row.get("sense") is not None and row["sense"] != f["sense"]:
        return "same p:q, opposite sense"
    if row.get("mu_note"):
        return f"same p:q; class-level only ({row['mu_note']})"
    c = row.get("C")
    same = "same p:q and sense" if row.get("sense") is not None else "same p:q (sense not printed)"
    if c is None or not (f["C"][0] - 0.01 <= c <= f["C"][1] + 0.01):
        return f"{same}; C {c} outside the family's computed range"
    ip = interp(f["members"], c)
    if ip is None or row.get("T") is None:
        return f"{same}; C inside the range but no T to compare (class-level)"
    d_t = abs(ip["T"] - row["T"]) / row["T"]
    d_pe = abs(ip["pe"] - row["perigee_alt_km"]) / max(abs(row["perigee_alt_km"]), 1.0)
    d_ps = abs(ip["ps"] - row["periselene_alt_km"]) / max(abs(row["periselene_alt_km"]), 1.0)
    if d_t < 1e-3 and d_pe < 0.10 and d_ps < 0.10:
        return f"LITERAL MEMBER (dT/T {d_t:.1e}, perigee {d_pe:.1%}, periselene {d_ps:.1%})"
    return (
        f"{same}; C inside the range but dT/T {d_t:.1e}, perigee {d_pe:.0%}, periselene {d_ps:.0%}"
    )


def main() -> None:
    preflight_search(
        task_no=1047,
        region_id="literal-collision-table-symmetric-and-gap-families",
        method=MethodCapability(
            genome="published Earth-Moon p:q orbits vs #1000/#1048 families",
            corrector="none (table lookup on computed members)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    rows = json.loads((SRC / "literal_check.json").read_text())
    fams = families()
    OUT.mkdir(parents=True, exist_ok=True)
    table: dict[str, Any] = {}
    md = []
    for key, f in sorted(fams.items()):
        verdicts = [
            {"source": r["source"], "orbit": r["orbit"], "verdict": verdict(r, f)} for r in rows
        ]
        hits = [v for v in verdicts if v["verdict"].startswith("LITERAL")]
        near = [v for v in verdicts if "same p:q" in v["verdict"]]
        table[key] = {
            "origin": f["origin"],
            "n_candidates": len(f["members"]),
            "computed_pq": f["pq"],
            "computed_sense": f["sense"],
            "pq_consistent": f["pq_consistent"],
            "C_range": f["C"],
            "literal_members": hits,
            "same_pq_rows": near,
            "rows": verdicts,
        }
        md.append(
            f"| {key} ({f['origin']}) | {len(f['members'])} | {f['pq'][0]}:{f['pq'][1]} "
            f"{'pro' if f['sense'] > 0 else 'ret'} | {f['C'][0]:.3f}-{f['C'][1]:.3f} | "
            + ("; ".join(f"{v['source']} {v['orbit']}" for v in hits) or "none")
            + " | "
            + ("; ".join(f"{v['source']} {v['orbit']}: {v['verdict']}" for v in near) or "none")
            + " |"
        )
    (OUT / "literal_table.json").write_text(json.dumps(table, indent=1) + "\n")
    head = [
        "| family (origin) | candidates | computed p:q, sense | C range | literal members | "
        "rows with the same p:q (verdict) |",
        "|---|---|---|---|---|---|",
    ]
    (OUT / "literal_table.md").write_text("\n".join(head + md) + "\n")
    print("\n".join(head + md))


if __name__ == "__main__":
    main()
