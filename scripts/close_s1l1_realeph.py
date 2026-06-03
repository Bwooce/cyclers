"""Constrained real-ephemeris closure of the S1L1 two-arc cycler.

Final step of the two-arc programme. The free epoch scan
(scripts/construct_s1l1_twoarc_realeph.py) showed real Mars eccentricity reaches
Mars V_inf ~2.83 (near the published 3.05) but does not close the full cycle (the
Earth side lands off-family). This drives the existing constrained optimiser
(optimise_cell_ephemeris) with the two-arc structure as the seed and the real
DE440 ephemeris, to see whether simultaneous flyby V_inf + turn-angle continuity
closes near the published 5.65 / 3.05.

Differs from tests/search/test_s1l1_idealised_multirev.py only in the two things
the real-ephemeris result justified: Ephemeris("astropy") instead of "circular",
and the 2030-03 low-Mars epoch instead of 2002. Sourced anchors 5.65/3.05
(spec §9); leg ToFs sourced/derived (Russell 4.991gG arcs 1.4612+2.8096 yr).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.4
TOF_SEED = [154.0, 379.0, 1026.0]  # E->M, M->E, E->E(multi-rev); arcs 533+1026 d


def _vinf_by_body(result: object) -> dict[str, float]:
    out: dict[str, float] = {}
    cyc = result.best_cycler  # type: ignore[attr-defined]
    if cyc is None:
        return out
    for enc in cyc.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def main() -> None:
    eph = Ephemeris("astropy")
    print("Constrained real-ephemeris S1L1 closure (E-M-E-E two-arc seed)")
    print(f"Targets E={VINF_E} M={VINF_M} km/s; tof seed {TOF_SEED}\n")
    # Epochs near the low-Mars-V_inf result (2030-03) plus the next synodic window.
    epochs = ["2030-02-20", "2030-03-22", "2030-04-21", "2032-01-11"]
    revs_options = [(0, 0, 2), (0, 0, 1)]
    best = None
    for revs in revs_options:
        branch = ("single", "single", "low" if revs[2] else "single")
        cell = Cell(
            bodies=("E", "M"),
            sequence=("E", "M", "E", "E"),
            period_k=2,
            per_leg_revs=revs,
            per_leg_branch=branch,
        )
        for epoch in epochs:
            try:
                result = optimise_cell_ephemeris(
                    cell,
                    eph,
                    vinf_cap=8.15,
                    priority_date_iso=epoch,
                    vinf_targets_kms={"E": VINF_E, "M": VINF_M},
                    n_starts=5,
                    seed=0,
                    tof_seed_days=TOF_SEED,
                )
            except Exception as exc:
                print(f"  revs={revs} {epoch}: ERROR {type(exc).__name__}: {exc}")
                continue
            v = _vinf_by_body(result)
            ve, vm = v.get("E", float("nan")), v.get("M", float("nan"))
            err = abs(ve - VINF_E) + abs(vm - VINF_M) if v else float("inf")
            ok = result.constraints_satisfied
            print(
                f"  revs={revs} {epoch}: feasible={ok}  "
                f"V_inf E={ve:.2f} M={vm:.2f}  "
                f"resid={result.closure_residual_kms:.3f}  |err|={err:.3f}"
            )
            if best is None or err < best[0]:
                best = (err, revs, epoch, ve, vm, ok)
    if best is not None:
        err, revs, epoch, ve, vm, ok = best
        print(
            f"\nBEST: revs={revs} epoch={epoch} feasible={ok} "
            f"V_inf E={ve:.2f} (t {VINF_E}) M={vm:.2f} (t {VINF_M}) |err|={err:.3f}"
        )
        within = abs(ve - VINF_E) < TOL and abs(vm - VINF_M) < TOL
        print(f"Within +/-{TOL} of BOTH published anchors: {within}")


if __name__ == "__main__":
    main()
