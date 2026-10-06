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
    if "flybys" not in rec:
        return ()
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


def is_pass(m: dict) -> bool:
    """Pre-registered pass (results note sec. 6.1): gate "pass" (near-180 included),
    and the independent re-propagation miss under 1 km at every encounter."""
    miss = m.get("max_encounter_miss_km")
    return m["status"] == "pass" and miss is not None and miss < 1.0


def collect(dirs: list[Path]) -> tuple[int, list[dict], dict[tuple, dict], dict[tuple, dict]]:
    """``(n_structures, zeros, physical groups, gate-passing groups)``."""
    recs: list[dict] = []
    n_struct = 0
    for d in dirs:
        sp = d / "structures.jsonl"
        if sp.exists():
            n_struct += sum(1 for _ in sp.open())
        zp = d / "zeros.jsonl"
        if zp.exists():
            recs += [json.loads(line) for line in zp.open()]
    groups: dict[tuple, dict] = {}
    for r in recs:
        if r.get("status") == "error":
            continue
        seq = flyby_seq(r)
        fwd, rev = canonical(seq), canonical(tuple(reversed(seq)))
        key = (r["k"], min(fwd, rev))
        g = groups.setdefault(key, {"members": [], "mirror_pair": fwd != rev})
        g["members"].append(r)
    passing = {k: g for k, g in groups.items() if any(is_pass(m) for m in g["members"])}
    return n_struct, recs, groups, passing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    n_struct, recs, groups, passing = collect(args.dirs)
    status = Counter(r["status"] for r in recs)
    n_struct_err = sum(
        1
        for d in args.dirs
        if (d / "structures.jsonl").exists()
        for line in (d / "structures.jsonl").open()
        if json.loads(line).get("error")
    )
    print(f"structures {n_struct}; zeros {len(recs)}; status {dict(status)}")
    print(
        f"ERRORS: structures that errored {n_struct_err} (not searched, NOT empty); "
        f"zero assessments that errored {status.get('error', 0)}"
    )
    print(f"physical cyclers (merged incl. mirrors) {len(groups)}; with a gate pass {len(passing)}")
    out = []
    for (k, seq), g in sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        best = min(
            (m for m in g["members"] if is_pass(m)),
            key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9,
        )
        line = {
            "k": k,
            "flybys": seq,
            "n_member_zeros": len(g["members"]),
            "mirror_pair": g["mirror_pair"],
            "keys": sorted({m["key"] for m in g["members"]}),
            "key": best["key"],
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
