"""#1048 (c): what the impacting solves of the #1000 multiple-shooting grid say about the seeds.
For each seed of data/1000_complement/shooting/grid.jsonl (and gaps.jsonl), the two-body
skeleton's closest approach to the Moon over one period (rotating frame, Moon at (1 - mu, 0))
and its perigee radius; status rates are tabulated against those and against p:q and sense.
Output: data/1000_complement/impact_analysis.json.
"""

from __future__ import annotations

import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_1000_shooting as sh

REPO = Path(__file__).resolve().parent.parent
D = REPO / "data" / "1000_complement"
L_KM = sh.base.L_KM


def skeleton_min_moon_km(p: int, q: int, sense: int, rp: float, om: float) -> float:
    a = (q / p) ** (2.0 / 3.0) * sh.GM_E ** (1.0 / 3.0)
    e = 1.0 - rp / a
    n = math.sqrt(sh.GM_E / a**3)
    t = np.linspace(0.0, 2 * math.pi * q, 6000)
    m_anom = n * t
    ecc = m_anom.copy()
    for _ in range(30):
        ecc = ecc - (ecc - e * np.sin(ecc) - m_anom) / (1 - e * np.cos(ecc))
    nu = 2 * np.arctan2(math.sqrt(1 + e) * np.sin(ecc / 2), math.sqrt(1 - e) * np.cos(ecc / 2))
    r = a * (1 - e * np.cos(ecc))
    ang = om + sense * nu - t  # inertial angle minus frame rotation
    x = r * np.cos(ang) - sh.MU
    y = r * np.sin(ang)
    return float(np.min(np.hypot(x - 1 + sh.MU, y)) * L_KM)


def main() -> None:
    preflight_search(
        task_no=1048,
        region_id="1000-grid-impact-analysis",
        method=MethodCapability(
            genome="#1000 skeleton seeds, analysis only",
            corrector="none",
            capability_tags=frozenset({"cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    out = {}
    for name in ("grid", "gaps"):
        recs = [json.loads(x) for x in (D / "shooting" / f"{name}.jsonl").read_text().splitlines()]
        recs = [r for r in recs if r["status"] != "no lunar-orbit crossing"]
        bins = {
            "<5000": (0, 5000),
            "5000-20000": (5000, 20000),
            "20000-60000": (20000, 60000),
            ">60000": (60000, 1e12),
        }
        by_bin: dict[str, Counter[str]] = defaultdict(Counter)
        by_rp: dict[int, Counter[str]] = defaultdict(Counter)
        by_pq: dict[str, Counter[str]] = defaultdict(Counter)
        for r in recs:
            st = r["status"] if r["status"] in ("converged", "impact") else "other"
            dmin = skeleton_min_moon_km(
                r["p"], r["q"], r["sense"], sh.RP[r["ir"]], sh.OMEGA[r["iw"]]
            )
            for b, (lo, hi) in bins.items():
                if lo <= dmin < hi:
                    by_bin[b][st] += 1
            by_rp[r["ir"]][st] += 1
            by_pq[f"{r['p']}:{r['q']} {'pro' if r['sense'] > 0 else 'ret'}"][st] += 1

        def rates(c: Counter[str]) -> dict[str, float | int]:
            n = sum(c.values())
            return {"n": n, **{k: round(v / n, 3) for k, v in c.items()}}

        out[name] = {
            "by_skeleton_min_moon_km": {b: rates(c) for b, c in by_bin.items()},
            "by_rp_km": {f"{sh.RP[i] * L_KM:.0f}": rates(c) for i, c in sorted(by_rp.items())},
            "by_pq_sense": {k: rates(c) for k, c in sorted(by_pq.items())},
        }
        print(name, json.dumps(out[name]["by_skeleton_min_moon_km"]))
        print(name, json.dumps(out[name]["by_rp_km"]))
    (D / "impact_analysis.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
