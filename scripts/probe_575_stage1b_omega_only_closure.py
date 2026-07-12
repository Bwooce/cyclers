"""#575 step 6 -- pipeline stage 1b: Omega-only inclined closure at the EXACT
commensurate tof (zero C3 drift by construction).

Stage 1 (``probe_575_stage1_inclination_closure.py``) found that #572's own
Nelder-Mead refinement (which moves BOTH Omega and ``tof_scale`` within a bounded
box) walks every one of the 6 basins it finds off the commensurate manifold
(0.23%-11.6% of T_syn/2), and every one of those 6 then FAILS the #575 C2 repeat
check post-refinement -- the exact C3 failure mode the #575 Fable correction warned
about.

This is the more conservative, complementary check: sweep ONLY Omega (Iapetus's
node alignment) at the SEED's EXACT commensurate ``tof_scale`` -- zero tof freedom,
so C3 drift is trivially 0 by construction -- and ask whether any Omega value alone
produces a gate-passing 3D-inclined closure. Reuses #572's own ``evaluate_point``
verbatim (no new gate/closure logic); the only difference from stage 1 is that
``tof_scale`` is never perturbed.

Run as::

    uv run python scripts/probe_575_stage1b_omega_only_closure.py
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
    evaluate_point,
)
from probe_575_stage1_inclination_closure import load_575_seeds, to_572_candidate  # noqa: E402
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
OUT_PATH = DATA_DIR / "probe_575_stage1b_omega_only_closure.jsonl"
N_OMEGA = 3600


def main() -> int:
    _t_syn, sqrt_papb = _t_syn_and_sqrt_papb()
    seeds = load_575_seeds()
    print(f"[575-stage1b] {len(seeds)} seeds, Omega-only sweep at n_omega={N_OMEGA}", flush=True)

    seq = ("Titan", "Iapetus", "Titan")
    out_records: list[dict[str, Any]] = [
        {
            "_meta": True,
            "task": "#575 step 6 pipeline stage 1b -- Omega-only closure (zero tof drift)",
            "n_seeds": len(seeds),
            "n_omega": N_OMEGA,
        }
    ]

    n_closure_found = 0
    n_repeat_confirmed = 0
    for seed in seeds:
        cand = to_572_candidate(seed, sqrt_papb)
        rel = cand["rel_offset_deg"]
        tof_scale = cand["tof_scale"]
        n_rev = cand["n_rev"]

        best_pt: dict[str, Any] | None = None
        for i in range(N_OMEGA):
            omega_deg = 360.0 * i / N_OMEGA
            pt = evaluate_point(rel, tof_scale, n_rev, omega_deg, INCLINATION_DEG)
            if pt is None:
                continue
            residual_near = pt["residual_kms"] < GATE_RESIDUAL_KMS
            if not residual_near:
                continue
            gate_pass, _ = candidate_passes_physical_gate(
                seq, tuple(pt["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
            )
            if not gate_pass:
                continue
            if best_pt is None or pt["residual_kms"] < best_pt["residual_kms"]:
                best_pt = pt

        record: dict[str, Any] = {
            "kind": "stage1b_result",
            "label": cand["label"],
            "n_commensurate_int": seed["n_commensurate_int"],
            "tof_scale_fixed": tof_scale,
            "closure_found": best_pt is not None,
        }
        if best_pt is None:
            print(f"[575-stage1b]   {cand['label']}: NO Omega-only closure found", flush=True)
            out_records.append(record)
            continue

        n_closure_found += 1
        record["omega_deg"] = best_pt["omega_deg"]
        record["residual_kms"] = best_pt["residual_kms"]
        record["vinf_kms"] = best_pt["vinf_kms"]
        record["c3_commensurability_drift_days"] = 0.0  # tof never perturbed, by construction

        params = TitanIapetusClosureParams(
            omega_deg=best_pt["omega_deg"],
            tof_scale=tof_scale,
            n_rev=tuple(seed["n_rev"]),
            m0_titan_deg=0.0,
            m0_iapetus_deg=rel,
            e_titan=0.0,
            e_iapetus=0.0,
            inclination_deg=INCLINATION_DEG,
        )
        c2 = repeat_check(f"{cand['label']}_stage1b", params)
        record["c2_repeat_check"] = {
            "n_cycles_completed": c2["n_cycles_completed"],
            "max_closure_residual_kms": c2["max_closure_residual_kms"],
            "repeats_to_machine_precision": c2["repeats_to_machine_precision"],
        }
        if c2["n_cycles_completed"] == N_CYCLES and c2["max_closure_residual_kms"] < 1e-6:
            n_repeat_confirmed += 1
        print(
            f"[575-stage1b]   {cand['label']}: omega={best_pt['omega_deg']:.3f} "
            f"residual={best_pt['residual_kms']:.3e} "
            f"c2_repeats={record['c2_repeat_check']['repeats_to_machine_precision']}",
            flush=True,
        )
        out_records.append(record)

    print(
        f"[575-stage1b] DONE: {n_closure_found}/{len(seeds)} seeds find an Omega-only "
        f"(zero-drift) closure; {n_repeat_confirmed}/{max(n_closure_found, 1)} repeat to "
        f"machine precision",
        flush=True,
    )
    out_records[0]["n_closure_found"] = n_closure_found
    out_records[0]["n_repeat_confirmed"] = n_repeat_confirmed
    with OUT_PATH.open("w", encoding="utf-8") as fh:
        for rec in out_records:
            fh.write(json.dumps(rec, default=str) + "\n")
    print(f"[575-stage1b] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
