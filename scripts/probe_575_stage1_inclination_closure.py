"""#575 step 6 -- pipeline stage 1: #572-equivalent inclination closure for the 9
genuine Titan-Iapetus symmetric-closure seeds.

Per the #575 spec: "if [the direct symmetric-closure search] produces genuine
symmetric closures, THEN run them back through the SAME #572->#573->#574 pipeline
(inclination closure test, population widening, eccentricity kill-gate...)". This
script is that first pipeline stage, reusing
``scripts/probe_572_titan_iapetus_3d_closure.py``'s own ``evaluate_point``/
``sweep_node_alignment``/``_smoke_test_reduction`` functions VERBATIM (imported, not
reimplemented) -- #572 is ALREADY Saturn/Titan/Iapetus-specific (``PRIMARY``/
``ANCHOR``/``FLYBY`` hardcoded), so no genericization is needed here, only its
``CANDIDATES`` tuple is swapped for the 9 #575 symmetric-enumeration seeds (in place
of #572's original 2 free-search seeds).

C3 discipline (LOAD-BEARING, #575 Fable correction)
----------------------------------------------------
#572's ``sweep_node_alignment`` refines EACH basin with a bounded Nelder-Mead over
BOTH Omega (+-15 deg) AND ``tof_scale`` (+-0.1) -- #573/#574 Stage A's own refinement
stages do the same. A refined ``tof_scale`` can silently walk away from the seed's
exact commensurate value (``tof = n*T_syn/2``), reproducing #574 Stage B's exact
0/15 failure mode (candidates that "close" but are no longer periodic) after all the
compute is spent. So for EVERY refined basin this script reports:

1. **Commensurability-drift diagnostic**: ``|tof_refined_days - n*T_syn/2|`` (days),
   where ``tof_refined_days = refined_tof_scale * sqrt(P_Titan*P_Iapetus)`` -- NOT
   hard-fixed to the commensurate value (a true inclined solution's tof need not be
   EXACTLY commensurate -- hard-fixing would false-kill genuine candidates), just
   measured and reported.
2. **Re-run of the #575 C2 repeat check** (``run_v2_saturn_3d``, reused from
   ``scripts/probe_575_titan_iapetus_repeat_check.py``, at the refined
   ``(omega_deg, tof_scale, n_rev, rel_offset_deg)`` and the REAL inclination
   (15.5 deg, still e=0 circular at this stage) -- if the refined basin stops
   repeating, that is the drift-off-the-manifold signal, flagged immediately here
   rather than only surfacing at a later gauntlet stage.

Discipline: NO catalogue writeback, reuses #572's gate machinery verbatim.

Run as::

    uv run python scripts/probe_575_stage1_inclination_closure.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from probe_572_titan_iapetus_3d_closure import (  # noqa: E402
    INCLINATION_DEG,
    _smoke_test_reduction,
    sweep_node_alignment,
)
from probe_575_titan_iapetus_repeat_check import (  # noqa: E402
    N_CYCLES,
    _t_syn_and_sqrt_papb,
    repeat_check,
)
from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.genome.titan_iapetus_corrector import TitanIapetusClosureParams  # noqa: E402
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

DATA_DIR = ROOT / "data"
ENUM_575_PATH = DATA_DIR / "enumerate_575_titan_iapetus_symmetric_closures.jsonl"
OUT_PATH = DATA_DIR / "probe_575_stage1_inclination_closure.jsonl"


def load_575_seeds() -> list[dict[str, Any]]:
    seeds = []
    with ENUM_575_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "pass" and d["anchor"] == "Titan":
                seeds.append(d)
    return seeds


def to_572_candidate(seed: dict[str, Any], sqrt_papb: float) -> dict[str, Any]:
    tof_scale = seed["tof_days"] / sqrt_papb
    label = f"n{seed['n_commensurate_int']}_nrev{seed['n_rev']}_rel{seed['rel_offset_deg']:.0f}"
    return {
        "label": label,
        "rel_offset_deg": seed["rel_offset_deg"],
        "tof_scale": tof_scale,
        "n_rev": tuple(seed["n_rev"]),
        "coplanar_residual_kms": seed["residual_kms"],
        "coplanar_iapetus_vinf_kms": seed["vinf_per_encounter_kms"][1],
        "coplanar_bend_deg": seed["max_bend_deg_per_encounter"][1],
        "_n_commensurate_int": seed["n_commensurate_int"],
    }


def main() -> int:
    print("[575-stage1] smoke test: 3D state generator reduces to coplanar...", flush=True)
    if not _smoke_test_reduction():
        print("[575-stage1] ABORT: smoke test failed.", flush=True)
        return 1

    t_syn, sqrt_papb = _t_syn_and_sqrt_papb()
    seeds = load_575_seeds()
    print(
        f"[575-stage1] {len(seeds)} seeds, T_syn={t_syn:.4f}d, inclination={INCLINATION_DEG}deg",
        flush=True,
    )

    out_records: list[dict[str, Any]] = [
        {
            "_meta": True,
            "task": "#575 step 6 pipeline stage 1 -- inclination closure (#572-equivalent)",
            "n_seeds": len(seeds),
            "t_syn_days": t_syn,
            "inclination_deg": INCLINATION_DEG,
        }
    ]

    n_closure_found = 0
    n_repeat_confirmed_after_refinement = 0
    for seed in seeds:
        cand = to_572_candidate(seed, sqrt_papb)
        n = seed["n_commensurate_int"]
        seed_target_tof_days = n * t_syn / 2.0
        print(
            f"[575-stage1] --- {cand['label']} (seed tof={seed_target_tof_days:.4f}d) ---",
            flush=True,
        )
        sweep = sweep_node_alignment(cand, n_omega=3600)

        best = sweep["best_overall"]
        record: dict[str, Any] = {
            "kind": "stage1_result",
            "label": cand["label"],
            "n_commensurate_int": n,
            "seed_rel_offset_deg": cand["rel_offset_deg"],
            "seed_tof_scale": cand["tof_scale"],
            "seed_target_tof_days": seed_target_tof_days,
            "n_basins_found": sweep["n_basins_found"],
            "n_feasible_omega_points": sweep["n_feasible_omega_points"],
        }

        if best is None:
            record["closure_found"] = False
            record["reason"] = "no_feasible_3d_point_at_any_omega"
            print(f"[575-stage1]   NO 3D closure found for {cand['label']}", flush=True)
            out_records.append(record)
            continue

        # Also require the #324 physical gate + near-coplanar residual (matching
        # #572's own "closure" definition) -- take the BEST GATE-PASSING basin, not
        # just the lowest-residual one (mirrors #572's own basin_evals logic).
        seq = ("Titan", "Iapetus", "Titan")
        closing_basins = []
        for b in sweep["basins"]:
            residual_near = b["residual_kms"] < GATE_RESIDUAL_KMS
            gate_pass, _ = candidate_passes_physical_gate(
                seq, tuple(b["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
            )
            if residual_near and gate_pass:
                closing_basins.append(b)

        if not closing_basins:
            record["closure_found"] = False
            record["reason"] = "no_basin_clears_both_gates"
            print(f"[575-stage1]   0/{sweep['n_basins_found']} basins clear both gates", flush=True)
            out_records.append(record)
            continue

        n_closure_found += 1
        refined = min(closing_basins, key=lambda b: b["residual_kms"])
        refined_tof_days = refined["tof_scale"] * sqrt_papb
        drift_days = abs(refined_tof_days - seed_target_tof_days)

        # C3 (i): commensurability-drift diagnostic.
        record["closure_found"] = True
        record["refined_omega_deg"] = refined["omega_deg"]
        record["refined_tof_scale"] = refined["tof_scale"]
        record["refined_tof_days"] = refined_tof_days
        record["refined_residual_kms"] = refined["residual_kms"]
        record["refined_vinf_kms"] = refined["vinf_kms"]
        record["c3_commensurability_drift_days"] = drift_days
        record["c3_commensurability_drift_fraction_of_t_syn_half"] = drift_days / (t_syn / 2.0)

        # C3 (ii): re-run the C2 repeat check at the refined point, real inclination,
        # e=0 (this stage's own model).
        params = TitanIapetusClosureParams(
            omega_deg=refined["omega_deg"],
            tof_scale=refined["tof_scale"],
            n_rev=tuple(seed["n_rev"]),
            m0_titan_deg=0.0,
            m0_iapetus_deg=cand["rel_offset_deg"],
            e_titan=0.0,
            e_iapetus=0.0,
            inclination_deg=INCLINATION_DEG,
        )
        c2_result = repeat_check(f"{cand['label']}_stage1", params)
        record["c2_repeat_check"] = {
            "n_cycles_completed": c2_result["n_cycles_completed"],
            "max_closure_residual_kms": c2_result["max_closure_residual_kms"],
            "max_drift_kms": c2_result["max_drift_kms"],
            "repeats_to_machine_precision": c2_result["repeats_to_machine_precision"],
        }
        if (
            c2_result["n_cycles_completed"] == N_CYCLES
            and c2_result["max_closure_residual_kms"] < 1e-6
        ):
            n_repeat_confirmed_after_refinement += 1

        print(
            f"[575-stage1]   CLOSURE: omega={refined['omega_deg']:.3f} "
            f"tof_scale={refined['tof_scale']:.4f} residual={refined['residual_kms']:.3e} "
            f"drift_from_commensurate={drift_days:.4f}d "
            f"({record['c3_commensurability_drift_fraction_of_t_syn_half'] * 100:.2f}% of T_syn/2) "
            f"c2_repeats={record['c2_repeat_check']['repeats_to_machine_precision']} "
            f"(n_cycles_completed={c2_result['n_cycles_completed']}, "
            f"closure_res={c2_result['max_closure_residual_kms']:.3e})",
            flush=True,
        )
        out_records.append(record)

    print(
        f"[575-stage1] DONE: {n_closure_found}/{len(seeds)} seeds find an inclined-closure "
        f"basin; {n_repeat_confirmed_after_refinement}/{n_closure_found} of those still "
        f"repeat to machine precision after refinement",
        flush=True,
    )

    out_records[0]["n_closure_found"] = n_closure_found
    out_records[0]["n_repeat_confirmed_after_refinement"] = n_repeat_confirmed_after_refinement
    with OUT_PATH.open("w", encoding="utf-8") as fh:
        for rec in out_records:
            fh.write(json.dumps(rec, default=str) + "\n")
    print(f"[575-stage1] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
