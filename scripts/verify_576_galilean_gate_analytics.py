"""#576 step 2 -- independent two-sided #324 physical-bend gate-feasibility
check for all 6 Galilean moon pairs.

Mirrors ``scripts/verify_571_gate_analytics.py``'s method exactly (same
``core/flyby.py::max_bend`` bisection for the 5-deg-bend V_inf ceiling, same
Hohmann-tangent minimum-achievable-V_inf floor), generalized from #571's
one-anchor-vs-N-small-moons shape (Titan fixed, checking each smaller moon)
to a genuine PAIRWISE, TWO-SIDED check over all C(4,2)=6 Galilean pairs: for
pair (A, B), BOTH of the following must hold for the pair to clear the gate
at either flyby node:

* A's own minimum-achievable V_inf (Hohmann-tangent, reaching B's orbital
  radius) <= A's own 5-deg-bend ceiling V_inf.
* B's own minimum-achievable V_inf (Hohmann-tangent, reaching A's orbital
  radius) <= B's own 5-deg-bend ceiling V_inf.

This is "two-sided" because a symmetric anchor-flyby-anchor closure
(``A -> B -> A``) needs BOTH the departure/return bend at A and the
mid-tour bend at B to be physically useful -- the #324 gate as embodied in
``search/physical_sanity.py::candidate_passes_physical_gate`` checks every
node in the sequence, not just one.

The #576 Fable plan review states (as Hohmann-tangent hand-computation, NOT
run through this script directly) that all 6 pairs clear with 2-5x margin,
even the widest (Io-Callisto: min-achievable V_inf 4.82/3.24 km/s vs
ceilings 8.25/7.77 km/s). This script INDEPENDENTLY recomputes those numbers
via code, per the "verify, don't just trust the transcription" mandate.

Run as::

    uv run python scripts/verify_576_galilean_gate_analytics.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.flyby import max_bend  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "verify_576_galilean_gate_analytics.jsonl"

MU_JUPITER = PRIMARIES["Jupiter"]
GALILEAN_MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")
TARGET_BEND_DEG = 5.0  # search/physical_sanity.py::DEFAULT_MIN_USEFUL_BEND_DEG


def bend_ceiling_vinf(body: str, target_bend_deg: float = TARGET_BEND_DEG) -> float:
    """Solve ``max_bend(mu, rp, vinf) == target_bend_deg`` for ``vinf`` (bisection),
    verbatim method from ``verify_571_gate_analytics.py``."""
    sat = SATELLITES[body]
    mu = sat.mu_km3_s2
    rp = sat.radius_eq_km + sat.safe_alt_km
    target_rad = math.radians(target_bend_deg)

    lo, hi = 1e-6, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if max_bend(mu, rp, mid) > target_rad:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def hohmann_vinf_at_r1(r1_km: float, r2_km: float, mu_primary: float) -> float:
    """Minimum-achievable V_inf at radius r1 for ANY conic connecting circular
    orbits at r1 and r2, verbatim from ``verify_571_gate_analytics.py``."""
    a_t = 0.5 * (r1_km + r2_km)
    v_circ = math.sqrt(mu_primary / r1_km)
    v_transfer = math.sqrt(mu_primary * (2.0 / r1_km - 1.0 / a_t))
    return abs(v_transfer - v_circ)


def bend_at_vinf(body: str, vinf: float) -> float:
    sat = SATELLITES[body]
    return math.degrees(max_bend(sat.mu_km3_s2, sat.radius_eq_km + sat.safe_alt_km, vinf))


def one_sided_check(body: str, other: str) -> dict[str, Any]:
    """``body``'s own min-achievable V_inf reaching ``other``'s radius, vs
    ``body``'s own 5-deg-bend ceiling."""
    sat_body = SATELLITES[body]
    sat_other = SATELLITES[other]
    ceiling = bend_ceiling_vinf(body)
    min_vinf = hohmann_vinf_at_r1(sat_body.sma_km, sat_other.sma_km, MU_JUPITER)
    bend_at_min = bend_at_vinf(body, min_vinf)
    feasible = min_vinf <= ceiling
    margin = ceiling / min_vinf if min_vinf > 0 else math.inf
    return {
        "body": body,
        "other": other,
        "ceiling_vinf_kms": ceiling,
        "min_achievable_vinf_kms": min_vinf,
        "bend_at_min_vinf_deg": bend_at_min,
        "feasible": feasible,
        "margin_ratio": margin,
    }


def main() -> int:
    print(f"Jupiter: GM={MU_JUPITER} km^3/s^2 (JPL SSD gm_de440 system GM)")
    print()
    for m in GALILEAN_MOONS:
        sat = SATELLITES[m]
        print(
            f"{m:10s} GM={sat.mu_km3_s2:10.4f} sma_km={sat.sma_km:10.1f} "
            f"safe_alt_km={sat.safe_alt_km:6.1f}"
        )
    print()

    pairs = list(itertools.combinations(GALILEAN_MOONS, 2))
    header = (
        f"{'pair':22s} {'A-ceiling':>10s} {'A-min-vinf':>11s} {'A-margin':>9s} "
        f"{'B-ceiling':>10s} {'B-min-vinf':>11s} {'B-margin':>9s} {'both-feasible':>13s}"
    )
    print(header)

    records: list[dict[str, Any]] = []
    n_both_feasible = 0
    for a, b in pairs:
        side_a = one_sided_check(a, b)
        side_b = one_sided_check(b, a)
        both_feasible = side_a["feasible"] and side_b["feasible"]
        n_both_feasible += int(both_feasible)
        rec = {
            "pair": [a, b],
            "side_a": side_a,
            "side_b": side_b,
            "both_sides_feasible": both_feasible,
        }
        records.append(rec)
        print(
            f"{a + '-' + b:22s} "
            f"{side_a['ceiling_vinf_kms']:10.3f} {side_a['min_achievable_vinf_kms']:11.3f} "
            f"{side_a['margin_ratio']:8.2f}x "
            f"{side_b['ceiling_vinf_kms']:10.3f} {side_b['min_achievable_vinf_kms']:11.3f} "
            f"{side_b['margin_ratio']:8.2f}x "
            f"{both_feasible!s:>13s}",
            flush=True,
        )

    print(
        f"\n[576-gate] DONE: {n_both_feasible}/{len(pairs)} pairs clear the two-sided "
        "#324 gate with positive margin on BOTH sides.",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#576 step 2 -- two-sided #324 gate-feasibility check, "
                    "all 6 Galilean pairs",
                    "mu_jupiter": MU_JUPITER,
                    "target_bend_deg": TARGET_BEND_DEG,
                    "n_pairs": len(pairs),
                    "n_both_sides_feasible": n_both_feasible,
                }
            )
            + "\n"
        )
        for rec in records:
            fh.write(json.dumps({"kind": "pair_result", **rec}) + "\n")
    print(f"[576-gate] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
