"""#320 follow-on — run the V3_qp gauntlet on the two QP-tori SILVERs.

The #320 first quasi_cycler sweep surfaced 2 SILVER QP-tori (Earth-Moon 3D
Neimark-Sacker brackets 2 & 10, k=4) that passed V1_qp + V2_qp but were blocked on
the V3_qp infra (#319 V3 was deferred). This driver builds V3_qp
(`data/validation/v3_qp.py`, REBOUND IAS15 independent-integrator invariance check),
regenerates each SILVER torus from its #299 bracket parent state, and runs V3_qp —
the next gate on the project's only genuinely-fresh + closing #320 candidates (a
distinct class from the V0-known repeated-moon cyclers).

Output: `data/silver_320_qp_v3_verdicts.jsonl`. A V3_qp PASS does NOT admit to the
catalogue (V4 + human review follow).
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v2_qp import run_v2_qp
from cyclerfinder.data.validation.v3_qp import run_v3_qp
from cyclerfinder.genome.qp_tori import correct_qp_torus

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "data" / "scan_320_qp_tori_3d_brackets.jsonl"
OUT = ROOT / "data" / "silver_320_qp_v3_verdicts.jsonl"


def main() -> int:
    rows = [json.loads(line) for line in SCAN.read_text().splitlines() if line.strip()]
    silvers = [r for r in rows if r.get("silver_status") == "SILVER"]
    print(
        f"[v3-qp] {len(silvers)} SILVER QP-tori from #320: {[r['bracket_idx'] for r in silvers]}",
        flush=True,
    )
    system = cr3bp.cr3bp_system("Earth", "Moon")
    verdicts = []
    for r in silvers:
        bidx = r["bracket_idx"]
        cid = f"qp-torus-320-bracket-{bidx}"
        parent = np.asarray(r["parent_state_nd"], dtype=np.float64)
        lam_a = complex(r["lam_a"]["re"], r["lam_a"]["im"])
        lam_b = complex(r["lam_b"]["re"], r["lam_b"]["im"])
        print(f"\n[{cid}] k={r['k']} T={r['T_a_TU']:.4f} — regenerating torus...", flush=True)
        t0 = time.monotonic()
        torus = correct_qp_torus(
            system,
            parent,
            float(r["T_a_TU"]),
            (lam_a, lam_b),
            k=int(r["k"]),
            n_long=16,
            n_trans=2,
            initial_torus_amplitude=5e-4,
            tol=1e-8,
            max_iter=40,
            independent_tol=1e-3,
            notes=f"v3_gauntlet_bracket_{bidx}",
        )
        print(
            f"  regenerated in {time.monotonic() - t0:.1f}s "
            f"(rho={torus.rho:.4f}, t_strob={torus.t_strob:.4f})",
            flush=True,
        )
        # Re-confirm V2 reproduces the stored SILVER, then run V3.
        v2 = run_v2_qp(cid, torus)
        t0 = time.monotonic()
        v3 = run_v3_qp(cid, torus)
        print(
            f"  V2_qp: passes={v2.passes_v2_qp} max_drift={v2.max_invariance_drift:.3e}", flush=True
        )
        print(
            f"  V3_qp ({time.monotonic() - t0:.0f}s) [{v3.integrator}]: passes={v3.passes_v3_qp}",
            flush=True,
        )
        print(
            f"    IAS15 invariance max_drift={v3.max_invariance_drift_ias15:.3e} "
            f"(floor {v3.drift_floor})",
            flush=True,
        )
        print(
            f"    IAS15-vs-DOP853 disagreement={v3.max_integrator_disagreement:.3e} "
            f"(floor {v3.agreement_floor})",
            flush=True,
        )
        rec = {
            "bracket_idx": bidx,
            "candidate_id": cid,
            "v2_passes_qp": v2.passes_v2_qp,
            "v2_max_drift_nondim": v2.max_invariance_drift,
            **asdict(v3),
        }
        verdicts.append(rec)
    OUT.write_text("\n".join(json.dumps(v) for v in verdicts) + "\n")
    n_pass = sum(1 for v in verdicts if v["passes_v3_qp"])
    print(f"\n[v3-qp] DONE: {n_pass}/{len(verdicts)} SILVER tori PASS V3_qp -> {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
