"""Galileo VEEGA reproduction with Phase 3 capabilities (#300 / #289 Phase 3).

This script combines all three Phase 3 capabilities:

  1. Multi-shell BFS over the Tisserand-Poincaré graph
     (:func:`find_mga_chains` with ``multi_shell=True``) — admits V_inf
     shifts across hetero flybys within the Strange-Longuski 2002 §12
     pump envelope.
  2. Per-leg TOF optimisation (:func:`optimise_chain_tofs`) — Nelder-Mead
     over (launch_epoch, leg_tofs) with loss = closure_residual +
     2*flyby_continuity. The alpha=2 weight reflects that ballistic
     continuity is a harder constraint than V_inf-match.
  3. DSM extension (:class:`DSMSpec` field on
     :class:`EpochLockedTrajectory`) — for the Earth-2 → Jupiter leg,
     which the ballistic Tisserand predicate cannot link below V_inf=11
     km/s; a small hand-placed DSM on that leg redirects to Jupiter.

The scan emits ``data/scan_300_galileo_veega_multishell.jsonl`` with one
record per (candidate, closure) pair. SOLID Phase-3 verdict if the
closure on Galileo VEEGA reaches <= 5 km/s; if the residual sits in
[5, 10] km/s the gap is documented as the DSM-placement-automation
deferral (Phase 4).

Per ``feedback_long_agents_commit_incrementally`` and the spec
discipline list — this script writes JSONL, not the catalogue.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import (
    DSMSpec,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.tisserand_mga_window import (
    MGAChainCandidate,
    find_mga_chains,
    optimise_chain_tofs,
)

OUTPUT_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "scan_300_galileo_veega_multishell.jsonl"
)

# Galileo published anchors (Diehl-Belbruno-Roberts 1986; KNOWN_CORPUS rev 568d8a4).
GALILEO_EARTH_LIFTOFF_UTC = "1989-10-18T16:53:40"
# In the Phase-3 substrate, ``launch_epoch_utc`` is the FIRST encounter of the
# sequence. For (V,E,E,J) that's the Venus flyby (1990-02-10), NOT the Earth
# liftoff. The Galileo public trajectory data:
#   Launch (Earth) : 1989-10-18T16:53:40 UTC
#   Venus flyby    : 1990-02-10
#   Earth-1 flyby  : 1990-12-08
#   Earth-2 flyby  : 1992-12-08
#   Jupiter        : 1995-12-07
GALILEO_FIRST_ENCOUNTER_UTC = "1990-02-10T00:00:00"
GALILEO_PUBLISHED_VINF = (5.3, 8.93, 8.93, 5.6)  # km/s at V / E1 / E2 / J


def _candidate_to_json(cand: MGAChainCandidate) -> dict:
    return {
        "sequence": list(cand.sequence),
        "vinf_tuple_kms": list(cand.vinf_tuple_kms),
        "leg_tofs_days": list(cand.leg_tofs_days),
        "launch_epoch_utc": cand.launch_epoch_utc,
        "tisserand_parameter": cand.tisserand_parameter,
        "chain_score": cand.chain_score,
    }


def main() -> None:
    print(f"[{time.strftime('%H:%M:%S')}] Galileo VEEGA Phase 3 scan starting")
    eph = Ephemeris("astropy")

    # Step 1: enumerate (V,E,E,J) candidates with multi-shell BFS.
    t0 = time.time()
    cands: list[MGAChainCandidate] = []
    for cand in find_mga_chains(
        launch_window=("1990-01-15T00:00:00", "1990-03-15T00:00:00"),
        planet_set=("V", "E", "J"),
        max_legs=4,
        vinf_grid_kms=(3.0, 5.0, 7.0, 9.0, 11.0),
        tof_box_days_per_leg=(60.0, 1200.0),
        epoch_step_days=10.0,
        a_range_au=(0.3, 8.0),
        multi_shell=True,
        pump_envelope_factor=1.0,
        start_body_filter=("V",),
    ):
        if cand.sequence == ("V", "E", "E", "J"):
            cands.append(cand)
    t1 = time.time()
    print(
        f"[{time.strftime('%H:%M:%S')}] Multi-shell BFS surfaced "
        f"{len(cands)} (V,E,E,J) candidates in {t1 - t0:.2f}s"
    )

    # Step 2: pick the candidate closest to the Galileo published launch + V_inf
    # pattern (lowest weighted score).
    def _score_against_galileo(cand: MGAChainCandidate) -> float:
        # First three V_inf entries should match (5.3, 8.93, 8.93). The last
        # entry (Jupiter) is structurally pinned by the Tisserand predicate
        # at >= 11 km/s; we ignore it in the seed-selection score.
        vinf_diff = sum(
            abs(a - b)
            for a, b in zip(
                cand.vinf_tuple_kms[:3],
                GALILEO_PUBLISHED_VINF[:3],
                strict=False,
            )
        )
        # Launch epoch closeness (days from 1990-02-10 Venus flyby).
        from astropy.time import Time

        dt_days = abs(
            (Time(cand.launch_epoch_utc) - Time(GALILEO_FIRST_ENCOUNTER_UTC)).to("day").value
        )
        return float(vinf_diff + 0.1 * dt_days)

    if not cands:
        print(f"[{time.strftime('%H:%M:%S')}] No (V,E,E,J) chains surfaced — bailing")
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_PATH.write_text("")
        return

    # Sort seeds by Galileo-alignment score; try them in order until one
    # gives a finite (non-Lambert-failing) baseline closure.
    cands.sort(key=_score_against_galileo)
    best_seed = None
    base_closure = None
    for cand_try in cands[:25]:
        try:
            total_tof_days = sum(cand_try.leg_tofs_days)
            from cyclerfinder.genome.epoch_aware_genome import _add_days_to_iso_utc

            end_utc_try = _add_days_to_iso_utc(cand_try.launch_epoch_utc, total_tof_days)
            traj_try = EpochLockedTrajectory(
                sequence=cand_try.sequence,
                leg_tofs_days=cand_try.leg_tofs_days,
                vinf_kms_at_encounters=cand_try.vinf_tuple_kms,
                launch_epoch_utc=cand_try.launch_epoch_utc,
                orbit_class="mga_tour",
                n_returns=1,
                validity_window_start_utc=cand_try.launch_epoch_utc,
                validity_window_end_utc=end_utc_try,
                periapsis_altitudes_km=(None, 965.0, 303.0, None),
            )
            closure_try = close_epoch_locked(
                traj_try,
                eph,
                closure_tol_kms=1.0e6,
                flyby_continuity_tol_kms=1.0e6,
                independent_cross_check=False,
                independent_tol_kms=1.0e6,
            )
            best_seed = cand_try
            base_closure = closure_try
            break
        except Exception as exc:
            # The Phase-2 geometric proposal is sometimes too far from a
            # well-conditioned Lambert geometry; try the next seed.
            print(
                f"[{time.strftime('%H:%M:%S')}] Seed "
                f"vinf={tuple(round(v, 2) for v in cand_try.vinf_tuple_kms)} "
                f"launch={cand_try.launch_epoch_utc} failed: {type(exc).__name__}"
            )
            continue

    if best_seed is None or base_closure is None:
        print(f"[{time.strftime('%H:%M:%S')}] No seed produced a finite baseline closure — bailing")
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_PATH.write_text("")
        return

    seed_score = _score_against_galileo(best_seed)
    print(
        f"[{time.strftime('%H:%M:%S')}] Picked seed: "
        f"vinf={tuple(round(v, 2) for v in best_seed.vinf_tuple_kms)} "
        f"launch={best_seed.launch_epoch_utc} score={seed_score:.2f}"
    )

    print(
        f"[{time.strftime('%H:%M:%S')}] Baseline ballistic closure: "
        f"residual={base_closure.closure_residual_kms:.3f} km/s "
        f"flyby_dv={base_closure.flyby_continuity_max_dv_kms:.3f} km/s "
        f"per_encounter_vinf={tuple(round(v, 2) for v in base_closure.per_encounter_vinf_kms)}"
    )

    # Step 4: TOF + epoch optimisation around the seed.
    t3 = time.time()
    opt_result = optimise_chain_tofs(
        best_seed,
        eph,
        periapsis_altitudes_km=(None, 965.0, 303.0, None),
        max_iter=80,
        epoch_search_half_width_days=15.0,
        tof_search_relative_half_width=0.4,
        alpha_flyby_continuity=2.0,
        accept_loss_kms=None,  # don't reject; record whatever the optimiser finds
        independent_cross_check=False,
    )
    if opt_result is None:
        print(f"[{time.strftime('%H:%M:%S')}] TOF optimiser returned None")
        opt_cand = best_seed
        opt_closure = base_closure
        opt_loss = (
            base_closure.closure_residual_kms + 2.0 * base_closure.flyby_continuity_max_dv_kms
        )
    else:
        opt_cand, opt_closure, opt_loss = opt_result
    print(
        f"[{time.strftime('%H:%M:%S')}] TOF optimised closure: "
        f"residual={opt_closure.closure_residual_kms:.3f} km/s "
        f"flyby_dv={opt_closure.flyby_continuity_max_dv_kms:.3f} km/s "
        f"loss={opt_loss:.3f} ({time.time() - t3:.1f}s)"
    )

    # Step 5: add a small hand-placed DSM on the E2 → J leg (leg index 2).
    # The Galileo trajectory's actual Earth-2 → Jupiter leg is NOT ballistic;
    # the ballistic Tisserand predicate only links E↔J above V_inf=11 km/s.
    # A small DSM mid-leg can redirect the E2 departure (V_inf=9 km/s) onto a
    # Jupiter-intercept arc.
    t4 = time.time()
    end_utc_opt = _add_days_to_iso_utc(opt_cand.launch_epoch_utc, sum(opt_cand.leg_tofs_days))
    dsm_traj = EpochLockedTrajectory(
        sequence=opt_cand.sequence,
        leg_tofs_days=opt_cand.leg_tofs_days,
        vinf_kms_at_encounters=opt_cand.vinf_tuple_kms,
        launch_epoch_utc=opt_cand.launch_epoch_utc,
        orbit_class="mga_tour",
        n_returns=1,
        validity_window_start_utc=opt_cand.launch_epoch_utc,
        validity_window_end_utc=end_utc_opt,
        periapsis_altitudes_km=(None, 965.0, 303.0, None),
        # Hand-placed DSM on leg 2 (E2 -> J) at fraction 0.3 (closer to
        # E2 to retain the heliocentric leveraging Galileo actually uses).
        dsm_specs=(
            DSMSpec(
                leg_index=2,
                fraction_along_leg=0.3,
                # ~50 m/s along heliocentric +X. The Galileo "perihelion ΔV"
                # was reported at ~62 m/s in the public mission profile.
                delta_v_kms=(0.062, 0.0, 0.0),
            ),
        ),
        notes="Phase 3 Galileo VEEGA + hand-placed E2->J DSM",
    )
    dsm_closure = close_epoch_locked(
        dsm_traj,
        eph,
        closure_tol_kms=1.0e6,
        flyby_continuity_tol_kms=1.0e6,
        independent_cross_check=True,
        independent_tol_kms=1.0e6,
    )
    print(
        f"[{time.strftime('%H:%M:%S')}] With DSM (hand-placed): "
        f"residual={dsm_closure.closure_residual_kms:.3f} km/s "
        f"flyby_dv={dsm_closure.flyby_continuity_max_dv_kms:.3f} km/s "
        f"dsm_dv={dsm_closure.dsm_delta_v_kms_per_leg!r} "
        f"({time.time() - t4:.1f}s)"
    )

    # Verdict.
    best_residual = min(
        base_closure.closure_residual_kms,
        opt_closure.closure_residual_kms,
        dsm_closure.closure_residual_kms,
    )
    if best_residual <= 5.0:
        verdict = "SOLID Phase-3 verdict: closure <= 5 km/s gate."
    elif best_residual <= 10.0:
        verdict = (
            f"Closure {best_residual:.2f} km/s in [5, 10]; remaining gap is "
            "automated DSM placement (Phase 4+) and higher-fidelity outer-planet "
            "ephemeris."
        )
    else:
        verdict = (
            f"Closure {best_residual:.2f} km/s > 10 km/s; the (V,E,E,J) "
            "ballistic-Tisserand chain at the Jupiter pinning V_inf=11 km/s "
            "is structurally not Galileo's trajectory at the 5 km/s level."
        )
    print(f"[{time.strftime('%H:%M:%S')}] {verdict}")

    # Emit JSONL.
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as f:
        for stage_label, cand, closure in (
            ("baseline_ballistic", best_seed, base_closure),
            ("tof_optimised", opt_cand, opt_closure),
            ("with_handplaced_dsm", opt_cand, dsm_closure),
        ):
            record = {
                "stage": stage_label,
                "candidate": _candidate_to_json(cand),
                "closure_residual_kms": closure.closure_residual_kms,
                "flyby_continuity_max_dv_kms": closure.flyby_continuity_max_dv_kms,
                "independent_check_residual_kms": closure.independent_check_residual_kms,
                "per_encounter_vinf_kms": list(closure.per_encounter_vinf_kms),
                "expected_published_vinf_kms": list(GALILEO_PUBLISHED_VINF),
                "dsm_delta_v_kms_per_leg": list(closure.dsm_delta_v_kms_per_leg),
                "converged_at_published_grade": closure.closure_residual_kms <= 5.0,
            }
            f.write(json.dumps(record) + "\n")
    print(
        f"[{time.strftime('%H:%M:%S')}] Wrote {OUTPUT_PATH.relative_to(OUTPUT_PATH.parents[1])} "
        f"(3 rows: baseline / tof_optimised / with_handplaced_dsm)"
    )


if __name__ == "__main__":
    main()
