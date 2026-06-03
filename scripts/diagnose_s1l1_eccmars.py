"""Track B diagnostic: eccentric-Mars (real-ephemeris) search for S1L1's SPEC
anchors V_inf_E = 5.65, V_inf_M = 3.05 km/s.

Rationale
---------
The circular-coplanar resonance construction reproduces the *coplanar* S1L1
anchors (4.90 / 4.98 km/s == Russell 2004's 4.99 / 5.10) but NOT the spec §9
higher-fidelity pair 5.65 / 3.05: the Mars 3.05 km/s specifically needs an
*eccentric* Mars orbit (Mars e=0.093), which the circular model deliberately
does not carry. This script drives the real-ephemeris optimiser
(``Ephemeris('astropy')`` -> DE440 states, eccentric Mars) seeded from the
resonance construction, and sweeps a coarse grid of launch epochs / Mars-
encounter phases to find the lowest ``|V_inf_M - 3.05|`` with V_inf_E near 5.65.

Provenance / golden discipline
------------------------------
The seed leg ToFs come from ``construct_resonant_cycler(a_au=1.30, e=0.257)``
(E->M ~152.6 d, M->E ~388.8 d), themselves COMPUTED from the sourced (a, e).
The 5.65 / 3.05 anchors are used only as phase-match TARGETS for epoch
resolution; nothing is asserted back against them and no tolerance is loosened.

Topology: E-M-E-E (per the corrected S1L1 nomenclature -- S/L are Earth-to-Earth
resonant intervals, not E<->M legs; the outbound E->M carries a ballistic Mars
flyby, then the short S1 + long L1 Earth-to-Earth intervals carry the vehicle
back). The L1 Earth->Earth leg is multi-rev; we sweep its (n_revs, branch).

NOT a test. Pure diagnostic.
"""

from __future__ import annotations

import itertools

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.resonance import synodic_period_days
from cyclerfinder.search.resonant_construct import construct_resonant_cycler
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M = 5.65, 3.05
SPEC_TOL = 0.4  # km/s; "reaches the spec anchors" means both within this.
VINF_CAP = 8.15


def vinf_by_body(result) -> dict[str, float]:
    """Max |V_inf| (in or out) per body across all encounters of the cycler."""
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def main() -> None:
    eph = Ephemeris(model="astropy")

    # --- Seed leg ToFs from the sourced (a, e) resonance construction. -------
    rc = construct_resonant_cycler(a_au=1.30, e=0.257, bodies=("E", "M"))
    em_seed = rc.leg_tofs_days["E->M"]  # ~152.6 d (the sourced outbound transit)
    me_seed = rc.leg_tofs_days["M->E"]  # ~388.8 d (short S1 Earth-to-Earth)
    t_syn = synodic_period_days("E", "M")
    period_2syn = 2.0 * t_syn
    l1_seed = period_2syn - em_seed - me_seed  # ~1018 d long L1 Earth-to-Earth
    print(
        f"resonance seed (days): E->M={em_seed:.1f}  M->E={me_seed:.1f}  "
        f"L1(E->E)={l1_seed:.1f}  (2-syn period={period_2syn:.1f})"
    )
    print(f"coplanar construction V_inf: E={rc.vinf_kms['E']:.2f}  M={rc.vinf_kms['M']:.2f}")
    print(f"targets: V_E={VINF_E}  V_M={VINF_M}  (spec-reach tol={SPEC_TOL} km/s)\n")

    # --- Coarse sweep grid. --------------------------------------------------
    # (1) Launch epochs: the priority date centres a +/-10yr window resolver
    #     internally, but seeding from successive E-M synodic priority dates lets
    #     the resolver settle into different Mars-encounter geometries (Mars'
    #     eccentric true anomaly at the flyby drives V_inf_M most strongly).
    # NB each optimise_cell_ephemeris call is ~50 s on the astropy backend (the
    # +/-10 yr window resolver dominates), so the grid is deliberately coarse:
    # 3 priority epochs x 3 Mars-phase offsets x 2 L1 configs = 18 runs.
    priority_dates = [
        "2000-05-01",
        "2002-08-05",  # the literature priority epoch
        "2004-08-01",
    ]
    # (2) Mars-encounter phase: perturb the E->M outbound ToF around the sourced
    #     ~152.6 d (this moves where on its eccentric orbit Mars is met -> the
    #     lever on V_inf_M = 3.05). M->E flexes to keep S1+L1 ~ 2-synodic.
    em_offsets = [-30.0, 0.0, 30.0]
    # (3) L1 Earth->Earth multi-rev leg branch/revs.
    l1_configs = [(1, "low"), (1, "high")]

    best = None  # (objective, ve, vm, label, tofs, converged, feasible)
    for pdate, em_off, (l1_revs, l1_branch) in itertools.product(
        priority_dates, em_offsets, l1_configs
    ):
        em = em_seed + em_off
        # Keep S1 (M->E) proportional so the loop still sums to ~2-synodic.
        me = me_seed
        l1 = period_2syn - em - me
        if l1 <= 0:
            continue
        seed = [em, me, l1]
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
                vinf_cap=VINF_CAP,
                priority_date_iso=pdate,
                vinf_targets_kms={"E": VINF_E, "M": VINF_M},
                tof_seed_days=seed,
                n_starts=5,
                seed=0,
            )
        except Exception as exc:
            print(
                f"  [{pdate} em{em_off:+.0f} r{l1_revs}{l1_branch[0]}] RAISED "
                f"{type(exc).__name__}: {str(exc)[:70]}"
            )
            continue
        v = vinf_by_body(res)
        ve = v.get("E", float("nan"))
        vm = v.get("M", float("nan"))
        # Rank primarily by closeness to the Mars anchor (3.05 needs eccentric
        # Mars -- it is the hard one), secondarily by the Earth anchor.
        objective = abs(vm - VINF_M) + 0.5 * abs(ve - VINF_E)
        label = f"{pdate} em{em_off:+.0f} L1=r{l1_revs}{l1_branch}"
        tofs = [(leg.t_arrive - leg.t_depart) / 86400.0 for leg in res.best_cycler.legs]
        print(
            f"  [{label}]: conv={res.converged} feas={res.constraints_satisfied} "
            f"V_E={ve:.2f}(d{ve - VINF_E:+.2f}) V_M={vm:.2f}(d{vm - VINF_M:+.2f}) "
            f"resid={res.closure_residual_kms:.3f} tofs={[round(t) for t in tofs]}"
        )
        if np.isfinite(objective) and (best is None or objective < best[0]):
            best = (objective, ve, vm, label, tofs, res.converged, res.constraints_satisfied)

    print("\n=== BEST (lowest |V_M - 3.05| + 0.5|V_E - 5.65|) ===")
    if best is None:
        print("(no finite-V_inf result across the sweep)")
        return
    _, ve, vm, label, tofs, conv, feas = best
    reached = abs(ve - VINF_E) <= SPEC_TOL and abs(vm - VINF_M) <= SPEC_TOL
    print(f"  {label}")
    print(f"  V_E={ve:.3f} (d{ve - VINF_E:+.3f})   V_M={vm:.3f} (d{vm - VINF_M:+.3f})")
    print(f"  converged={conv} feasible={feas} tofs={[round(t) for t in tofs]}")
    print(f"  REACHES spec anchors (both within {SPEC_TOL} km/s): {reached}")


if __name__ == "__main__":
    main()


# FINDING (2026-06-03):
#   NO HIT. The eccentric-Mars real-ephemeris optimiser (Ephemeris('astropy'),
#   DE440 states) seeded from the resonance construction does NOT reach the spec
#   anchors 5.65 / 3.05 km/s. Sweep: 3 priority epochs (2000-05/2002-08/2004-08)
#   x 3 Mars-phase E->M offsets (-30/0/+30 d) x 2 L1 (n_revs,branch) = 18 runs.
#
#   Best achieved (lowest |V_M-3.05| + 0.5|V_E-5.65|):
#       2004-08-01, em+0, L1=(1,"high"):  V_E=17.22 (d+11.57)  V_M=10.28 (d+7.23)
#   -> off by +11.6 / +7.2 km/s; NOT within the 0.4 km/s spec-reach tolerance.
#
#   Why it misses (consistent with scripts/diagnose_s1l1_realeph.py's finding):
#   * Only the L1=(1,"high") config ever RESOLVES a real launch window
#     (converged=True). Every (1,"low") config returns the "no real window"
#     SENTINEL (equispaced tofs [520,520,520], V_inf~17, resid=inf) -- at the
#     ~1018 d seed L1 ToF the "high" branch is the only feasible multi-rev arc.
#   * Where a window IS resolved, the maintenance solve lands in the spec §10
#     degenerate high-V_inf basin (V_E~17, V_M~10): the dominant error is the
#     still-mis-phased SINGLE-rev E->M / M->E FRONT legs, not the L1 multi-rev
#     plumbing. The Mars-phase (em_offset) sweep does NOT move the resolver off
#     that basin -- it snaps to the same degenerate window regardless of offset.
#   * Mars eccentricity (the reason 3.05 needs eccentric Mars) IS present in the
#     astropy states and does shift V_M slightly across epochs (10.45 @ 2002 vs
#     10.28 @ 2004), but the shift is ~0.2 km/s -- far short of the ~7 km/s gap.
#     The blocker is therefore the E-M-E-E leg-ToF / encounter phasing inside
#     the maintenance solve, NOT the absence of eccentric Mars.
#
#   Conclusion: Track B's eccentric-Mars search does NOT reproduce 5.65 / 3.05
#   under the current epoch resolver + maintenance engine. The S1L1 spec-anchor
#   rediscovery gate stays xfail; no tolerance was loosened. Reaching the anchors
#   needs the asymmetric front-leg phasing fixed (the resolver must escape the
#   degenerate symmetric/high-V_inf basin) -- a maintenance-engine change that is
#   the Apply phase's job, not this diagnostic's.
