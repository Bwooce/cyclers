"""Diagnostic: run S1L1's correct E-M-E-E topology through the NEW real-eph
optimiser (multi-rev + multi-encounter) with the Stage-3 multi-seed epoch
resolver and sourced V_inf targets. Goal: see how close it lands to the
sourced 5.65/3.05 anchors and WHY it lands off-family if it does.

NOT a test. Pure diagnostic.
"""

from __future__ import annotations

import itertools

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M = 5.65, 3.05
T_SYN_EM_D = 779.9
SEED = [154.0, 379.0, T_SYN_EM_D * 2 - 154.0 - 379.0]  # E->M, M->E, L1 E->E loop


def vinf_by_body(result) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def main() -> None:
    eph = Ephemeris(model="astropy")
    print(f"seed leg ToFs (days): {[round(s, 1) for s in SEED]}  (sum={sum(SEED):.0f})")
    for l1_revs, l1_branch in itertools.product((1, 2), ("low", "high")):
        cell = Cell(
            bodies=("E", "M"),
            sequence=("E", "M", "E", "E"),
            period_k=2,
            per_leg_revs=(0, 0, l1_revs),
            per_leg_branch=("single", "single", l1_branch),
        )
        try:
            res = optimise_cell_ephemeris(
                cell,
                eph,
                vinf_cap=8.15,
                priority_date_iso="2002-08-05",
                vinf_targets_kms={"E": VINF_E, "M": VINF_M},
                tof_seed_days=SEED,
                n_starts=5,
                seed=0,
            )
        except Exception as exc:
            print(f"  {cell.id}: RAISED {type(exc).__name__}: {str(exc)[:80]}")
            continue
        v = vinf_by_body(res)
        ve, vm = v.get("E", float("nan")), v.get("M", float("nan"))
        tofs = [(leg.t_arrive - leg.t_depart) / 86400.0 for leg in res.best_cycler.legs]
        print(
            f"  {cell.id}: converged={res.converged} feasible={res.constraints_satisfied} "
            f"V_E={ve:.2f}(d{ve - VINF_E:+.2f}) V_M={vm:.2f}(d{vm - VINF_M:+.2f}) "
            f"resid={res.closure_residual_kms:.3f} tofs={[round(t) for t in tofs]}"
        )


if __name__ == "__main__":
    main()

# FINDING (2026-06-03, original):
#   All E-M-E-E configs returned the "no real window" SENTINEL (converged=False,
#   resid=inf, equispaced tofs [520,520,520], V_inf ~17). The Stage-3 multi-seed
#   epoch resolver returned t0=None for the 4-encounter E-M-E-E cell because
#   epoch resolution (find_candidate_windows / phase_match) Lambert-solved every
#   leg SINGLE-REV and could not phase-match S1L1's L1 multi-rev same-body
#   Earth->Earth leg.
#
# UPDATE (2026-06-03, after per-leg n_revs/branch threading):
#   PhaseSignature now carries optional leg_revs/leg_branches and the phase-match
#   Lambert (_vinf_at_lambert / _mismatch_at_date) selects the matching
#   (n_revs, branch) solution; optimise_cell_ephemeris threads cell.per_leg_revs
#   / per_leg_branch into it. The epoch resolver no longer ALWAYS returns the
#   inf sentinel: the r001/high config now RESOLVES an epoch (converged=True,
#   finite resid, non-equispaced tofs ~[118,849,1213], V_E~26.5 / V_M~23.6).
#   The (1,low) / (2,*) configs still return the sentinel, but NOT because of
#   the single-rev bug: at the seed ToF (~1027 d) the (1,"high") branch is the
#   only feasible multi-rev branch, and even when a branch resolves the V_inf
#   mismatch is dominated by the still-mis-phased SINGLE-rev E->M (154 d,
#   V_inf~17) and M->E (379 d, V_inf~10) FRONT legs -- i.e. the E-M-E-E leg
#   ToF / encounter geometry, not the L1 multi-rev plumbing. The maintenance
#   solve lands off-family, so the published 5.65/3.05 anchors are NOT yet
#   reproduced and the S1L1 rediscovery gate stays xfail. Remaining work is the
#   E-M-E-E leg-ToF / encounter phasing, separate from the multi-rev epoch
#   resolver fixed here. (Aldrin's 2 direct legs resolve unchanged.)
