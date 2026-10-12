"""#975 P2 v2 coincidence check (note sec. 11.5 item 3): both traces closed, and each member of one
lies on the other's family.

Metric 1 (every member): distance in p = (x1(0), y2(0), T/2pi) to the other family by orthogonal
correction (`distance_to_loop` seeded from the three nearest members of the other trace).
Metric 2 ("matched T"): for each member whose unit tangent has |T component| >= 0.1, (x1(0), y2(0))
corrected by Newton at the member's fixed T, from the nearest member of the other trace; the
distance in (x1(0), y2(0)). Doubling-back check: no member j within 1e-6 of member n - j. Member
closure: the half residuals (x2, y1 at T/2) of every stored member, recomputed (criterion 1e-9).

    uv run python scripts/check_975_p2_coincidence.py --a p2_up_loop_v2 --b p2_down_loop_v2
Output: data/975_varin/p2_v2_coincidence.json.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search import varin_closed_families as vcf

LOOPS = Path(__file__).resolve().parent.parent / "data" / "975_varin" / "loops"


def _load(name: str) -> tuple[float, list[dict[str, Any]], np.ndarray]:
    rs = [json.loads(x) for x in (LOOPS / f"{name}.jsonl").read_text().splitlines()]
    return float(rs[0]["mu"]), rs, np.array([r["p"] for r in rs if "p" in r])


def _matched_t(p: np.ndarray, other: np.ndarray, mu: float) -> float:
    """Newton in (x1, y2) at the fixed T of p, from the nearest member of ``other``."""
    j = int(np.argmin(np.linalg.norm(other - p, axis=1)))
    x = other[j][:2].copy()
    t_full = 2.0 * math.pi * float(p[2])
    for _ in range(20):
        r, jac = vcf.half_residuals_jac(float(x[0]), float(x[1]), t_full, mu)
        if float(np.max(np.abs(r))) < 1e-11:
            break
        x = x + np.linalg.solve(jac[:, :2], -r)
    else:
        return float("nan")
    return float(np.linalg.norm(x - p[:2]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    args = ap.parse_args()
    out: dict[str, Any] = {}
    for a, b in ((args.a, args.b), (args.b, args.a)):
        mu, ra, pa = _load(a)
        _, _, pb = _load(b)
        ends = [r for r in ra if "end" in r]
        n = len(pa) - 1
        doubling = [j for j in range(1, n // 2) if float(np.linalg.norm(pa[j] - pa[n - j])) < 1e-6]
        m1, m2, res, n_skip = [], [], [], 0
        for i, p in enumerate(pa):
            near = np.argsort(np.linalg.norm(pb - p, axis=1))[:3]
            m1.append(vcf.distance_to_loop(p, [pb[k] for k in sorted(near)], mu, max_refine=20))
            r, jac = vcf._loop_eval(p, mu)
            res.append(float(np.max(np.abs(r))))
            t = vcf._tangent(jac, None)
            if abs(float(t[2])) >= 0.1:
                m2.append(_matched_t(p, pb, mu))
            else:
                n_skip += 1
            if i % 25 == 0:
                print(f"{time.strftime('%H:%M:%S')} {a} on {b}: {i}/{len(pa)}", flush=True)
        rec = {
            "members": len(pa),
            "end": ends,
            "doubling_back_pairs": doubling,
            "metric1_max": float(np.max(m1)),
            "metric1_argmax": int(np.argmax(m1)),
            "metric2_max": float(np.nanmax(m2)),
            "metric2_nan": int(np.sum(np.isnan(m2))),
            "metric2_n": len(m2),
            "metric2_skipped_small_T_component": n_skip,
            "max_half_residual": float(np.max(res)),
        }
        out[f"{a} on {b}"] = rec
        print(f"{a} on {b}: {rec}", flush=True)
    closed = all(
        v["end"] and v["end"][-1]["end"] == "closed" and v["end"][-1]["closure"] < 1e-6
        for v in out.values()
    )
    within = all(
        v["metric1_max"] < 1e-5
        and v["metric2_max"] < 1e-5
        and v["metric2_nan"] == 0
        and v["max_half_residual"] < 1e-9
        for v in out.values()
    )
    no_doubling = not any(v["doubling_back_pairs"] for v in out.values())
    out["P2_pass"] = bool(closed and within and no_doubling)
    path = LOOPS.parent / "p2_v2_coincidence.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"P2 {'PASS' if out['P2_pass'] else 'FAIL'} -> {path}", flush=True)


if __name__ == "__main__":
    main()
