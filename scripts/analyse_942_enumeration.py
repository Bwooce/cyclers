"""#942/#943: collect enumeration shards, merge physically identical zeros, report.

Two zeros are the SAME cycler when the cyclic sequence of their massive flybys
(body, V-infinity to 1 m/s, demanded turn to 0.1 deg) agrees up to rotation;
they are MIRROR twins when one sequence is the reverse of the other (time
reversal in the circular model). Split labels of one conic at a massless
target collapse here too, since only massive flybys enter the key.

Usage: uv run python scripts/analyse_942_enumeration.py DIR [DIR ...] [--json OUT]
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def flyby_seq(rec: dict) -> tuple[tuple[str, float, float], ...]:
    return tuple(
        (f["body"], round(float(f["vinf_kms"]), 3), round(float(f["turn_deg"]), 1))
        for f in rec["flybys"]
        if "ratio" in f
    )


def canonical(seq: tuple) -> tuple:
    if not seq:
        return seq
    rots = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(rots)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    recs = []
    n_struct = 0
    for d in args.dirs:
        sp = d / "structures.jsonl"
        if sp.exists():
            n_struct += sum(1 for _ in sp.open())
        zp = d / "zeros.jsonl"
        if zp.exists():
            recs += [json.loads(line) for line in zp.open()]
    status = Counter(r["status"] for r in recs)
    groups: dict[tuple, dict] = {}
    for r in recs:
        seq = flyby_seq(r)
        fwd, rev = canonical(seq), canonical(tuple(reversed(seq)))
        key = (r["k"], min(fwd, rev))
        g = groups.setdefault(key, {"members": [], "mirror_pair": fwd != rev})
        g["members"].append(r)
    passing = {k: g for k, g in groups.items() if any(m["status"] == "pass" for m in g["members"])}
    print(f"structures {n_struct}; zeros {len(recs)}; status {dict(status)}")
    print(f"physical cyclers (merged incl. mirrors) {len(groups)}; with a gate pass {len(passing)}")
    out = []
    for (k, seq), g in sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        best = min(
            g["members"], key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9
        )
        line = {
            "k": k,
            "flybys": seq,
            "n_member_zeros": len(g["members"]),
            "mirror_pair": g["mirror_pair"],
            "keys": sorted({m["key"] for m in g["members"]}),
            "worst_ratio": best["worst_ratio"],
            "max_turn_deg": best["max_turn_deg"],
            "min_required_alt_km": best["min_required_alt_km"],
            "status_hm_floor": best["status_hm_floor"],
            "r_min_km": best["r_min_km"],
            "r_max_km": best["r_max_km"],
            "vinf_kms": best["vinf_kms"],
            "x_days": best["x_days"],
            "max_encounter_miss_km": best["max_encounter_miss_km"],
        }
        out.append(line)
        print(json.dumps(line, default=float))
    if args.json:
        args.json.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
