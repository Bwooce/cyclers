"""S1L1 closure attempt with the CORRECT E-M-E-E topology + sourced leg seed.

Topology (per the S1/L1 nomenclature): outbound E->M (~154 d), M->E completing
the short S1 Earth-to-Earth interval (~379 d), then the long L1 Earth-to-Earth
resonant loop (~1030 d, multi-rev, no Mars encounter). Sourced anchors: V_inf
5.65 (E) / 3.05 (M) km/s, 2-synodic period. The S1/L1 split is a COMPUTED
hypothesis seed, not a golden value.

NOT a test. Records whether this topology closes to the anchors.

# FINDING (2026-06-02):
#   NO HIT. The correct E-M-E-E topology + sourced [154, 379, 1030] d seed does
#   NOT close in the circular-coplanar model: every (L1 n_revs, branch) variant
#   lands at V_inf_E ~25-39 km/s (degenerate), residual >26 km/s. This confirms
#   the blocker is the MODEL, not the topology: S1L1's 154-day E->M leg is
#   near-hyperbolic in circular-coplanar (same root cause as the Aldrin cycler).
#   Closing S1L1 requires the real ephemeris (Mars eccentricity) AND a real-eph
#   optimiser that handles the L1 multi-rev leg (the maintenance engine behind
#   optimise_cell_ephemeris is currently single-rev only).
"""

from __future__ import annotations

import itertools

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_idealized
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.3
VINF_CAP = max(VINF_E, VINF_M) + 2.5

# Seed leg ToFs (days): E->M 154, M->E ~379, E->E (L1 loop) ~1030.
SEED_TOFS_D = [154.0, 379.0, 1030.0]


def vinf_by_body(result) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def main() -> None:
    eph = Ephemeris(model="circular")
    seq = ("E", "M", "E", "E")
    hits = []
    for l1_revs, l1_branch in itertools.product((1, 2, 3), ("low", "high")):
        cell = Cell(
            bodies=("E", "M"),
            sequence=seq,
            period_k=2,
            per_leg_revs=(0, 0, l1_revs),
            per_leg_branch=("single", "single", l1_branch),
        )
        # interior encounter epochs (seconds): cumulative ToFs for M and E2.
        t_m = SEED_TOFS_D[0]
        t_e2 = SEED_TOFS_D[0] + SEED_TOFS_D[1]
        warm = [(t_m * SECONDS_PER_DAY, t_e2 * SECONDS_PER_DAY)]
        try:
            res = optimise_cell_idealized(
                cell, eph, vinf_cap=VINF_CAP, n_starts=5, seed=0, use_de=True, warm_starts=warm
            )
        except Exception as exc:
            print(f"  {cell.id}: raised {type(exc).__name__}: {str(exc)[:60]}")
            continue
        v = vinf_by_body(res)
        ok = (
            res.constraints_satisfied
            and abs(v.get("E", 1e9) - VINF_E) < TOL
            and abs(v.get("M", 1e9) - VINF_M) < TOL
        )
        tag = "HIT" if ok else "   "
        print(
            f"{tag} {cell.id}: E={v.get('E', float('nan')):.3f} M={v.get('M', float('nan')):.3f} "
            f"resid={res.closure_residual_kms:.4f} feasible={res.constraints_satisfied}"
        )
        if ok:
            tofs = [(leg.t_arrive - leg.t_depart) / SECONDS_PER_DAY for leg in res.best_cycler.legs]
            hits.append((cell.id, l1_revs, l1_branch, tofs))
    print("\n=== HITS ===")
    for cid, revs, branch, tofs in hits:
        print(f"{cid} L1_revs={revs} branch={branch} leg_tofs_days={[round(t, 1) for t in tofs]}")
    if not hits:
        print("(none) — E-M-E-E does not close to 5.65/3.05 in circular-coplanar under these seeds")


if __name__ == "__main__":
    main()
