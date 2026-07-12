"""#574 Stage B -- V2->V3->V4->V4-strict gauntlet on the 15 Stage-A eccentric survivors.

Background
----------
#574 Stage A (continuation-in-eccentricity kill gate, `scripts/run_574_titan_iapetus_
eccentric_kill_gate.py`) found 15 deduped survivors among the 17 eccentricity-robust
#573 Titan-Iapetus 3D-closure branches (branch ids 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12,
13, 16, 19, 20 -- verified directly against `data/probe_574_titan_iapetus_eccentric_
kill_gate.jsonl` below, not hand-transcribed). This script runs the productized
V2->V3->V4->V4-strict gauntlet built in Stage B
(`cyclerfinder.data.validation.{v2_saturn_3d,v3_saturn_3d,v4_saturn,v4_saturn_strict}`,
on top of `cyclerfinder.genome.titan_iapetus_corrector`) on all 15.

Load-bearing finding this run surfaces (read before interpreting results)
----------------------------------------------------------------------------
Unlike the Uranian #558-#569 family -- whose `tof_scale`/`rel_offset_deg` were found by
#563's DEDICATED symmetric/commensurate-closure enumeration, specifically so a fixed-TOF
multi-cycle repeat reproduces the SAME encounter geometry -- the Titan-Iapetus #571-#574
closures were found by a free (Omega, tof_scale, rel_offset) search for a SINGLE
V_inf-continuity closure, with NO periodicity/commensurability constraint at all. V2 here
re-solves the SAME fixed-TOF, fixed-n_rev legs at k*cycle_period advances (mirroring
`v2_moontour`'s own "cycle ToFs held fixed" discipline verbatim) -- and empirically, for
this family, the Titan-Iapetus transfer ANGLE changes by tens of degrees cycle-to-cycle
(confirmed by direct instrumentation before this script was written: branch 1's leg-0
transfer angle is 117.6 deg at cycle 0, 78.0 deg at cycle 1, 38.3 deg at cycle 2 -- even
at e=0). For n_rev>=1 branches this can make the SAME-n_rev Lambert transfer physically
NOT EXIST past cycle 0 (a genuine multi-rev feasibility floor, not a solver bug); for
n_rev=0 branches a solution always exists but the per-cycle V_inf-continuity residual
typically blows out to several km/s, far past even the #566 quasi-bounded 0.5 km/s floor.
This is reported here as a real, structurally-grounded finding -- NOT a bug in this
script's cycle-repeat logic (which is a direct, verified generalization of the
established `v2_moontour` pattern) and NOT evidence the underlying eccentric closures
themselves are wrong (branch 1's stage-0 residual is machine-precision, per #574 Stage A).
It means: most of Stage A's 15 "closures" are single V_inf-continuity transfers, not
repeating cyclers, under this literal multi-cycle test.

Gating discipline (mirrors #566's `v2_admits_downstream`)
-------------------------------------------------------------
* V2 is run at `n_cycles` in {3, 5, 10} (the #566 grid) for every candidate.
* `v2_status` per candidate: PASS (all three n_cycles pass the strict floor),
  FAIL_QUASI_BOUNDED (never strict-passes, but every completed n_cycles run stays under
  the 0.5 km/s quasi-bounded closure-residual floor AND completes all requested cycles),
  or FAIL_UNBOUNDED (anything else, INCLUDING a Lambert branch simply ceasing to exist
  before completing `n_cycles`).
* V3/V4/V4-strict are only attempted for a candidate at `n_cycles=3` if V2 completed 3
  cycles at n_cycles=3 (a hard requirement -- `run_v3_saturn_3d` needs a 3-cycle
  `v2_verdict.per_cycle` to compare against; candidates that Lambert-fail before cycle 3
  cannot feed a downstream stage at all, structurally, not by choice).
* V4-strict runs at ONE representative launch epoch (2000-06-21T00:00:00, the same
  #338/#566/#559-established reference epoch reused for consistency, not because Saturn
  has any special significance to that date) -- NOT a full annual/daily sweep (out of
  this dispatch's scope per the #574 Stage B spec).

Discipline
----------
* READ-ONLY on `data/probe_574_titan_iapetus_eccentric_kill_gate.jsonl` and every
  `cyclerfinder.data.validation.*`/`cyclerfinder.genome.titan_iapetus_corrector` module.
* NO catalogue writeback.
* Framing (mandatory): ANY result here is quasi-cycler-class evidence about our own
  idealized + real-ephemeris-tested search space, same standing as #312's Uranian
  family -- NOT a ballistic-cycler finding and NOT a novelty claim.

Run as::

    uv run python scripts/run_574_stageB_saturn_gauntlet.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.satellites import PRIMARIES  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.data.validation.v2_saturn_3d import (  # noqa: E402
    V2Saturn3DVerdict,
    run_v2_saturn_3d,
)
from cyclerfinder.data.validation.v3_saturn_3d import (  # noqa: E402
    run_v3_saturn_3d,
)
from cyclerfinder.data.validation.v4_saturn import (  # noqa: E402
    run_v4_saturn,
)
from cyclerfinder.data.validation.v4_saturn_strict import (  # noqa: E402
    run_v4_saturn_strict,
    verdict_to_jsonable,
)
from cyclerfinder.genome.titan_iapetus_corrector import TitanIapetusClosureParams  # noqa: E402
from cyclerfinder.verify.spice_kernels import ensure_sat441_kernel  # noqa: E402

PROBE_574_PATH = ROOT / "data" / "probe_574_titan_iapetus_eccentric_kill_gate.jsonl"
OUT_JSONL = ROOT / "data" / "gauntlet_574_saturn_stageB.jsonl"

SURVIVOR_IDS: tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 19, 20)
N_CYCLES_GRID: tuple[int, ...] = (3, 5, 10)
QUASI_BOUNDED_CLOSURE_FLOOR_KMS = 0.5
"""Mirrors #566's own quasi_cycler_3cycle test threshold verbatim (the same test
that let the #312/#327 SILVER itself proceed past a strict V2 drift-floor FAIL)."""
LAUNCH_EPOCH_UTC = "2000-06-21T00:00:00"

_REGION_ID = "saturn-titan-iapetus-574-stageB-v2v3v4v4strict-gauntlet-2026-07-12"
_METHOD = MethodCapability(
    genome=(
        "the 15 #574 Stage-A eccentric-Keplerian-robust Titan-Iapetus closures, run "
        "through the newly-productized V2->V3->V4->V4-strict Saturn validation chain"
    ),
    corrector=(
        "cyclerfinder.data.validation.{v2_saturn_3d,v3_saturn_3d,v4_saturn,"
        "v4_saturn_strict} on top of cyclerfinder.genome.titan_iapetus_corrector -- "
        "new Stage-B modules, not a discovery search"
    ),
    capability_tags=frozenset(
        {"saturn", "titan", "iapetus", "v2", "v3", "v4", "v4-strict", "spice", "gauntlet"}
    ),
    git_sha="working-tree",
)


@dataclass(frozen=True)
class Candidate:
    branch_id: int
    n_rev: tuple[int, int]
    params: TitanIapetusClosureParams
    circular_iapetus_bend_deg: float
    final_residual_kms: float


def load_candidates() -> list[Candidate]:
    """Load the 15 Stage-A survivors' FINAL (real-eccentricity) closure params
    directly from the Stage-A jsonl -- never hand-transcribed."""
    recs = []
    with PROBE_574_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                recs.append(json.loads(line))
    by_id = {r["branch_id"]: r for r in recs if r.get("kind") == "branch_result"}

    # The jsonl's raw `survives` flag is TRUE for BOTH the main eccentricity-robust
    # population AND the floor-hugger control group (4/5 of which also happened to
    # survive per the Stage-A result -- see OUTSTANDING.md #574 Stage-A RESULT). The
    # 15-member Stage-A PASS population is specifically the DEDUPED survivors of the
    # 17 ECCENTRICITY-ROBUST branches (`eccentricity_robust=True`), not the raw
    # survives-anywhere set -- filter on both fields together, matching the
    # OUTSTANDING.md #574 text verbatim ("17 eccentricity-robust branches: 15 deduped
    # survivors").
    survives_check = [
        bid for bid in by_id if by_id[bid]["survives"] and by_id[bid]["eccentricity_robust"]
    ]
    if sorted(survives_check) != sorted(SURVIVOR_IDS):
        raise RuntimeError(
            f"SURVIVOR_IDS {sorted(SURVIVOR_IDS)} does not match the jsonl's actual "
            f"eccentricity_robust+survives=True branch ids {sorted(survives_check)} -- "
            "refusing to proceed on a stale/mistranscribed candidate list."
        )

    candidates = []
    for bid in SURVIVOR_IDS:
        b = by_id[bid]
        final = b["stages"][-1]
        params = TitanIapetusClosureParams(
            omega_deg=final["omega_deg"],
            tof_scale=final["tof_scale"],
            n_rev=tuple(b["n_rev"]),
            m0_titan_deg=final["m0_titan_deg"],
            m0_iapetus_deg=final["m0_iapetus_deg"],
            e_titan=final["e_titan"],
            e_iapetus=final["e_iapetus"],
        )
        candidates.append(
            Candidate(
                branch_id=bid,
                n_rev=tuple(b["n_rev"]),
                params=params,
                circular_iapetus_bend_deg=b["circular_iapetus_bend_deg"],
                final_residual_kms=final["residual_kms"],
            )
        )
    return candidates


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _v2_row(cand: Candidate, nc: int, v: V2Saturn3DVerdict) -> dict[str, Any]:
    return {
        "kind": "v2_verdict",
        "branch_id": cand.branch_id,
        "n_cycles": nc,
        "n_cycles_completed": v.n_cycles_completed,
        "max_drift_kms": v.max_drift_kms,
        "max_closure_residual_kms": v.max_closure_residual_kms,
        "passes_v2": v.passes_v2,
    }


def main() -> int:
    t_start = time.time()
    sha = _git_sha()
    print(f"[574B] Saturn Stage-B V2->V3->V4->V4-strict gauntlet -- sha={sha}", flush=True)

    preflight_search(
        task_no=574,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(SURVIVOR_IDS),
        override_reason=(
            "read-only validation-reporting gauntlet on the 15 already-enumerated #574 "
            "Stage-A survivor closures using the newly-productized Saturn V2->V3->V4->"
            "V4-strict chain -- not a discovery sweep; 15 candidates x <=3 downstream "
            "stages measured at well under 2s/candidate total (see docstring), so this "
            "is seconds-scale, mirroring the #566 precedent's own justification."
        ),
    )

    try:
        ensure_sat441_kernel()
    except RuntimeError as exc:
        print(f"[574B] FATAL: {exc}", file=sys.stderr)
        return 1

    mu = PRIMARIES["Saturn"]
    candidates = load_candidates()
    print(f"[574B] loaded {len(candidates)} Stage-A survivors: {SURVIVOR_IDS}", flush=True)

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "_meta": True,
            "task": "#574 Stage B -- Saturn Titan-Iapetus V2->V3->V4->V4-strict gauntlet",
            "git_sha": sha,
            "survivor_branch_ids": list(SURVIVOR_IDS),
            "n_cycles_grid": list(N_CYCLES_GRID),
            "quasi_bounded_closure_floor_kms": QUASI_BOUNDED_CLOSURE_FLOOR_KMS,
            "v4_strict_launch_epoch_utc": LAUNCH_EPOCH_UTC,
        }
    )

    summary: list[dict[str, Any]] = []

    out_fh = OUT_JSONL.open("w", encoding="utf-8")

    def _write(rec: dict[str, Any]) -> None:
        out_fh.write(json.dumps(rec, default=str) + "\n")
        out_fh.flush()

    for rec in rows:
        _write(rec)

    for cand in candidates:
        t_cand = time.time()
        print(f"\n[574B] === branch {cand.branch_id} (n_rev={cand.n_rev}) ===", flush=True)

        v2_by_nc: dict[int, V2Saturn3DVerdict] = {}
        for nc in N_CYCLES_GRID:
            v2 = run_v2_saturn_3d(f"branch{cand.branch_id}", cand.params, mu=mu, n_cycles=nc)
            v2_by_nc[nc] = v2
            _write(_v2_row(cand, nc, v2))
            print(
                f"[574B]   V2(nc={nc}): completed={v2.n_cycles_completed} "
                f"passes={v2.passes_v2} max_drift={v2.max_drift_kms:.3e}km "
                f"max_closure={v2.max_closure_residual_kms:.3e}km/s",
                flush=True,
            )

        v2_all_pass = all(v2_by_nc[nc].passes_v2 for nc in N_CYCLES_GRID)
        v2_all_quasi_bounded = all(
            (not v2_by_nc[nc].passes_v2)
            and v2_by_nc[nc].n_cycles_completed == nc
            and v2_by_nc[nc].max_closure_residual_kms < QUASI_BOUNDED_CLOSURE_FLOOR_KMS
            for nc in N_CYCLES_GRID
        )
        if v2_all_pass:
            v2_status = "PASS"
        elif v2_all_quasi_bounded:
            v2_status = "FAIL_QUASI_BOUNDED"
        else:
            v2_status = "FAIL_UNBOUNDED"

        # Downstream (V3/V4/V4-strict) needs a 3-cycle V2 verdict to compare against --
        # a HARD structural requirement (see module docstring), not a discretionary gate.
        v2_nc3 = v2_by_nc[3]
        can_run_downstream = v2_nc3.n_cycles_completed == 3

        cand_summary: dict[str, Any] = {
            "branch_id": cand.branch_id,
            "n_rev": list(cand.n_rev),
            "circular_iapetus_bend_deg": cand.circular_iapetus_bend_deg,
            "v2_status": v2_status,
            "v2_all_pass": v2_all_pass,
            "v2_all_quasi_bounded": v2_all_quasi_bounded,
            "v2_nc3_max_closure_residual_kms": v2_nc3.max_closure_residual_kms,
            "v2_nc3_max_drift_kms": v2_nc3.max_drift_kms,
            "downstream_attempted": can_run_downstream,
        }

        if not can_run_downstream:
            cand_summary["chain_verdict"] = "FAIL_AT_V2_LAMBERT_INFEASIBLE"
            cand_summary["elapsed_s"] = time.time() - t_cand
            summary.append(cand_summary)
            print(
                f"[574B]   V2 did not complete 3 cycles (Lambert branch ceased to exist) "
                f"-- downstream V3/V4/V4-strict cannot run for branch {cand.branch_id}.",
                flush=True,
            )
            continue

        v3 = run_v3_saturn_3d(
            f"branch{cand.branch_id}", cand.params, mu=mu, v2_verdict=v2_nc3, n_cycles=3
        )
        _write(
            {
                "kind": "v3_verdict",
                "branch_id": cand.branch_id,
                "passes_v3": v3.passes_v3,
                "drift_agreement_kms": v3.drift_agreement_kms,
                "integrator": v3.integrator,
            }
        )
        print(
            f"[574B]   V3: passes={v3.passes_v3} agreement={v3.drift_agreement_kms:.3e}km "
            f"integrator={v3.integrator}",
            flush=True,
        )

        v4 = run_v4_saturn(
            f"branch{cand.branch_id}", cand.params, mu_primary=mu, v3_verdict=v3, n_cycles=3
        )
        _write(
            {
                "kind": "v4_verdict",
                "branch_id": cand.branch_id,
                "passes_v4": v4.passes_v4,
                "drift_agreement_kms": v4.drift_agreement_kms,
                "bounded_drift_survives": v4.bounded_drift_survives,
            }
        )
        print(
            f"[574B]   V4: passes={v4.passes_v4} agreement={v4.drift_agreement_kms:.3e}km "
            f"bounded={v4.bounded_drift_survives}",
            flush=True,
        )

        v4s = run_v4_saturn_strict(
            f"branch{cand.branch_id}",
            cand.params,
            LAUNCH_EPOCH_UTC,
            mu_primary=mu,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=3,
        )
        row = verdict_to_jsonable(v4s)
        row["branch_id"] = cand.branch_id
        _write(row)
        print(
            f"[574B]   V4-strict: passes={v4s.passes_v4_strict} "
            f"bounded={v4s.bounded_drift_survives} "
            f"agreement_vs_v3={v4s.drift_agreement_kms_vs_v3:.3e}km "
            f"(e_titan={v4s.eccentricity_used_e_titan:.4f} "
            f"e_iap={v4s.eccentricity_used_e_iapetus:.4f} "
            f"i_iap={v4s.inclination_used_deg_iapetus:.2f}deg)",
            flush=True,
        )
        fail_modes = [c.failure_mode for c in v4s.per_cycle if c.failure_mode != "converged"]

        cand_summary["v3_passes"] = v3.passes_v3
        cand_summary["v3_drift_agreement_kms"] = v3.drift_agreement_kms
        cand_summary["v4_passes"] = v4.passes_v4
        cand_summary["v4_drift_agreement_kms"] = v4.drift_agreement_kms
        cand_summary["v4_bounded_drift_survives"] = v4.bounded_drift_survives
        cand_summary["v4_strict_passes"] = v4s.passes_v4_strict
        cand_summary["v4_strict_drift_agreement_vs_v3_kms"] = v4s.drift_agreement_kms_vs_v3
        cand_summary["v4_strict_bounded_drift_survives"] = v4s.bounded_drift_survives
        cand_summary["v4_strict_failure_modes"] = fail_modes
        cand_summary["v4_strict_eccentricity_used_e_titan"] = v4s.eccentricity_used_e_titan
        cand_summary["v4_strict_eccentricity_used_e_iapetus"] = v4s.eccentricity_used_e_iapetus
        cand_summary["v4_strict_inclination_used_deg_iapetus"] = v4s.inclination_used_deg_iapetus

        if (
            v2_status in ("PASS", "FAIL_QUASI_BOUNDED")
            and v3.passes_v3
            and v4.passes_v4
            and v4s.passes_v4_strict
        ):
            chain_verdict = "PASS_AS_QUASI_CYCLER" if v2_status == "FAIL_QUASI_BOUNDED" else "PASS"
        elif v2_status in ("PASS", "FAIL_QUASI_BOUNDED") and v3.passes_v3 and v4.passes_v4:
            chain_verdict = "FAIL_AT_V4_STRICT"
        elif v2_status in ("PASS", "FAIL_QUASI_BOUNDED") and v3.passes_v3:
            chain_verdict = "FAIL_AT_V4"
        elif v2_status in ("PASS", "FAIL_QUASI_BOUNDED"):
            chain_verdict = "FAIL_AT_V3"
        else:
            chain_verdict = "FAIL_AT_V2_UNBOUNDED_DOWNSTREAM_COMPUTED_ANYWAY"
        cand_summary["chain_verdict"] = chain_verdict
        cand_summary["elapsed_s"] = time.time() - t_cand
        summary.append(cand_summary)
        print(f"[574B]   -> chain_verdict={chain_verdict}", flush=True)

    _write(
        {
            "_meta": True,
            "kind": "headline",
            "candidates": summary,
            "writeback_to_catalogue": False,
            "elapsed_s": time.time() - t_start,
        }
    )
    out_fh.close()

    print("\n[574B] === SUMMARY ===", flush=True)
    for s in summary:
        print(f"[574B] branch {s['branch_id']:2d}: {s['chain_verdict']}", flush=True)
    print(f"[574B] wrote {OUT_JSONL}", flush=True)
    print(f"[574B] total elapsed {time.time() - t_start:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
