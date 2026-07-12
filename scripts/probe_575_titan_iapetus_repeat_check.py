"""#575 C2 -- two-sided multi-cycle repeat-instrumentation control.

Reuses the per-cycle repeat check ALREADY BUILT for #574 Stage B
(``cyclerfinder.data.validation.v2_saturn_3d.run_v2_saturn_3d``, imported and driven by
``scripts/run_574_stageB_saturn_gauntlet.py``) -- not reimplemented. That module re-solves
the SAME fixed-tof, fixed-n_rev Lambert legs at ``k * cycle_period`` time advances and
reports, per cycle, whether the Lambert closure still converges, the V_inf-continuity
residual, and the rendezvous drift vs cycle 0. #574 Stage B used this exact mechanism
(at real Titan/Iapetus eccentricity) to discover that #571's free-closure candidates
do NOT repeat cycle-to-cycle (branch 1: Lambert infeasible past cycle 0). Here it is run
at ``e_titan=e_iapetus=inclination_deg=0`` (:func:`kepler_state_3d` is pinned, by its own
C2 positive control in ``tests/genome/test_titan_iapetus_corrector.py``, to reduce EXACTLY
to the circular-coplanar model #575's symmetric enumeration was constructed in), so this
check operates in the SAME model the new candidates were found in.

Two sides, per #575's C2 spec:

(a) POSITIVE -- every #575 symmetric-enumeration survivor (from
    ``data/enumerate_575_titan_iapetus_symmetric_closures.jsonl``) must repeat to
    machine precision at every cycle (this is GUARANTEED by the commensurate-tof
    symmetric construction, not a new claim -- a failure here is a red flag the
    construction has a bug).
(b) NEGATIVE -- a known-bad #571 candidate (branch 1, the #574 Stage B write-up's own
    documented example: Lambert infeasible past cycle 0, even at e=0) must NOT repeat.
    Loaded directly from ``data/probe_574_titan_iapetus_eccentric_kill_gate.jsonl``'s
    own ``ecc_fraction=0.0`` stage (never hand-transcribed), proving this instrumentation
    actually discriminates real periodicity from a one-off closure.

Also runs the cheap cross-check #575's C2 spec calls for: scan #571's original 187
raw candidates (``data/scan_571_saturn_titan_iapetus.jsonl`` +
``data/scan_571_saturn_iapetus_titan.jsonl``) for any point within tolerance of
rel_offset in {0, 180} deg AND commensurate tof (``|tof - n*T_syn/2| < tolerance`` for
some integer n) -- reporting whether this intersects the new #575 enumeration.

Discipline: NO catalogue writeback, read-only on every #571/#574 data file.

Run as::

    uv run python scripts/probe_575_titan_iapetus_repeat_check.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.data.validation.v2_saturn_3d import run_v2_saturn_3d  # noqa: E402
from cyclerfinder.genome.titan_iapetus_corrector import TitanIapetusClosureParams  # noqa: E402
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402

DATA_DIR = ROOT / "data"
ENUM_575_PATH = DATA_DIR / "enumerate_575_titan_iapetus_symmetric_closures.jsonl"
PROBE_574_PATH = DATA_DIR / "probe_574_titan_iapetus_eccentric_kill_gate.jsonl"
SCAN_571_TITAN_ANCHORED = DATA_DIR / "scan_571_saturn_titan_iapetus.jsonl"
SCAN_571_IAPETUS_ANCHORED = DATA_DIR / "scan_571_saturn_iapetus_titan.jsonl"
OUT_PATH = DATA_DIR / "probe_575_titan_iapetus_repeat_check.jsonl"

N_CYCLES = 3


def _t_syn_and_sqrt_papb() -> tuple[float, float]:
    mu = PRIMARIES["Saturn"]
    titan = SATELLITES["Titan"]
    iapetus = SATELLITES["Iapetus"]
    n_t = _mean_motion_rad_day(mu, titan.sma_km)
    n_i = _mean_motion_rad_day(mu, iapetus.sma_km)
    p_t = 2.0 * math.pi / n_t
    p_i = 2.0 * math.pi / n_i
    t_syn = 1.0 / abs(1.0 / p_t - 1.0 / p_i)
    return t_syn, math.sqrt(p_t * p_i)


def load_575_titan_anchored_survivors() -> list[dict[str, Any]]:
    """The 9 Titan-anchored (anchor='Titan') #575 symmetric-enumeration survivors --
    directly reusable as ``TitanIapetusClosureParams`` since the corrector module's
    ``ANCHOR`` constant is fixed to Titan. (The 9 Iapetus-anchored survivors are the
    SAME physical closures under the anchor/flyby swap convention #563's own data
    already established -- residuals match exactly by ``n``; see the #575 dispatch
    report.)"""
    recs = []
    with ENUM_575_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "pass" and d["anchor"] == "Titan":
                recs.append(d)
    return recs


def params_from_575_record(rec: dict[str, Any], sqrt_papb: float) -> TitanIapetusClosureParams:
    return TitanIapetusClosureParams(
        omega_deg=0.0,
        tof_scale=rec["tof_days"] / sqrt_papb,
        n_rev=tuple(rec["n_rev"]),
        m0_titan_deg=0.0,
        m0_iapetus_deg=rec["rel_offset_deg"],
        e_titan=0.0,
        e_iapetus=0.0,
        inclination_deg=0.0,
    )


def load_branch1_negative_control() -> TitanIapetusClosureParams:
    """Branch 1's own e=0 (``ecc_fraction=0.0``) stage, loaded directly from
    #574 Stage-A's own jsonl -- never hand-transcribed. This is the #574 Stage B
    write-up's own documented negative example (Lambert infeasible past cycle 0,
    117.6/78.0/38.3 deg transfer-angle drift)."""
    with PROBE_574_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "branch_result" and d.get("branch_id") == 1:
                stage0 = d["stages"][0]
                assert stage0["ecc_fraction"] == 0.0
                return TitanIapetusClosureParams(
                    omega_deg=stage0["omega_deg"],
                    tof_scale=stage0["tof_scale"],
                    n_rev=tuple(d["n_rev"]),
                    m0_titan_deg=stage0["m0_titan_deg"],
                    m0_iapetus_deg=stage0["m0_iapetus_deg"],
                    e_titan=0.0,
                    e_iapetus=0.0,
                    inclination_deg=0.0,
                )
    raise RuntimeError("branch_id 1 not found in probe_574 data")


def repeat_check(candidate_id: str, params: TitanIapetusClosureParams) -> dict[str, Any]:
    mu = PRIMARIES["Saturn"]
    v2 = run_v2_saturn_3d(candidate_id, params, mu=mu, n_cycles=N_CYCLES)
    return {
        "candidate_id": candidate_id,
        "n_rev": list(params.n_rev),
        "tof_scale": params.tof_scale,
        "m0_titan_deg": params.m0_titan_deg,
        "m0_iapetus_deg": params.m0_iapetus_deg,
        "n_cycles_completed": v2.n_cycles_completed,
        "max_drift_kms": v2.max_drift_kms,
        "max_closure_residual_kms": v2.max_closure_residual_kms,
        "per_cycle": [
            {
                "cycle_index": c.cycle_index,
                "converged_legs": c.converged_legs,
                "n_legs": c.n_legs,
                "rendezvous_drift_kms": c.rendezvous_drift_kms,
                "closure_residual_kms": c.closure_residual_kms,
            }
            for c in v2.per_cycle
        ],
        # NOTE (found empirically this run, see #575 dispatch report): this check is
        # deliberately based on Lambert-completion + CLOSURE-residual repeat ONLY, not
        # on run_v2_saturn_3d's inertial-frame `max_drift_kms` -- that metric measures
        # ABSOLUTE Titan position offset in the Saturn-inertial frame, which is large
        # (100k-2.4M km) for EVERY candidate here (both the genuine #575 symmetric
        # constructions AND, cross-checked below, even #312's OWN catalogued SILVER
        # under the identical v2_moontour mechanism: max_drift_kms=515,499 km,
        # passes_v2=False -- #312 itself is only ever a FAIL_QUASI_BOUNDED under this
        # exact metric, never a strict PASS). It reflects T_syn not being commensurate
        # with the anchor moon's OWN individual orbital period (a separate axis from
        # the rel_offset/tof symmetric-closure condition #575 constructs against), not
        # a defect in the discovery-stage construction. The discovery-stage "does this
        # repeat" question (#575's C2) is about whether the SAME Lambert transfer
        # geometry (n_rev branch, V_inf magnitudes) keeps closing cycle-to-cycle --
        # exactly what `closure_residual_kms` measures and what killed #571 branch 1
        # (Lambert infeasible past cycle 0, not a drift-floor issue).
        "repeats_to_machine_precision": bool(
            v2.n_cycles_completed == N_CYCLES and v2.max_closure_residual_kms < 1e-6
        ),
    }


def cross_check_571_against_symmetric_condition(
    t_syn: float, tol_deg: float = 2.0, tol_tof_days: float = 0.5
) -> dict[str, Any]:
    """Scan #571's original 187 candidates (``n_all_gates_passed`` summed across the
    Titan-anchored (69) + Iapetus-anchored (118) directions, per
    ``data/scan_571_saturn_titan_pairs_index.jsonl`` and the "187 Titan-Iapetus
    candidates" figure quoted throughout ``data/OUTSTANDING.md`` #571/#552) for any
    point within tolerance of rel_offset in {0, 180} deg AND commensurate tof --
    report the intersection (empty or not) honestly, per #575's C2 spec."""
    hits: list[dict[str, Any]] = []
    n_scanned = 0
    for path, direction in (
        (SCAN_571_TITAN_ANCHORED, "Titan-anchored"),
        (SCAN_571_IAPETUS_ANCHORED, "Iapetus-anchored"),
    ):
        if not path.exists():
            continue
        with path.open(encoding="utf-8") as fh:
            for line in fh:
                d = json.loads(line)
                if d.get("kind") != "gate_result" or not d.get("all_gates_passed"):
                    continue
                rec = d.get("record", {})
                rel = rec.get("rel_offset_deg")
                tof = rec.get("tof_days")
                if rel is None or tof is None:
                    continue
                n_scanned += 1
                near_0 = min(rel % 360.0, 360.0 - (rel % 360.0)) <= tol_deg
                near_180 = abs((rel % 360.0) - 180.0) <= tol_deg
                if not (near_0 or near_180):
                    continue
                nearest_n = round(2.0 * tof / t_syn)
                if nearest_n < 1:
                    continue
                target_tof = nearest_n * t_syn / 2.0
                if abs(tof - target_tof) <= tol_tof_days:
                    hits.append(
                        {
                            "direction": direction,
                            "anchor": d.get("anchor"),
                            "flyby": d.get("flyby"),
                            "rel_offset_deg": rel,
                            "tof_days": tof,
                            "nearest_n": nearest_n,
                            "target_commensurate_tof_days": target_tof,
                            "n_rev": rec.get("n_rev"),
                        }
                    )
    return {
        "n_571_records_scanned": n_scanned,
        "tol_deg": tol_deg,
        "tol_tof_days": tol_tof_days,
        "n_hits_near_symmetric_condition": len(hits),
        "hits": hits,
    }


def main() -> int:
    t_syn, sqrt_papb = _t_syn_and_sqrt_papb()

    print("[575-C2] (a) POSITIVE: 9 Titan-anchored #575 symmetric survivors...", flush=True)
    survivors = load_575_titan_anchored_survivors()
    positive_results = []
    n_positive_repeat = 0
    for rec in survivors:
        params = params_from_575_record(rec, sqrt_papb)
        cid = f"n{rec['n_commensurate_int']}_nrev{rec['n_rev']}_rel{rec['rel_offset_deg']:.0f}"
        res = repeat_check(cid, params)
        positive_results.append(res)
        if res["repeats_to_machine_precision"]:
            n_positive_repeat += 1
        print(
            f"[575-C2]   {cid}: n_cycles_completed={res['n_cycles_completed']} "
            f"max_drift_kms={res['max_drift_kms']:.3e} "
            f"max_closure_residual_kms={res['max_closure_residual_kms']:.3e} "
            f"repeats={res['repeats_to_machine_precision']}",
            flush=True,
        )
    print(
        f"[575-C2] POSITIVE control: {n_positive_repeat}/{len(survivors)} repeat to "
        f"machine precision",
        flush=True,
    )

    print("[575-C2] (b) NEGATIVE: #571 branch 1 (e=0 stage)...", flush=True)
    branch1_params = load_branch1_negative_control()
    negative_result = repeat_check("branch1_e0", branch1_params)
    print(
        f"[575-C2]   branch1_e0: n_cycles_completed={negative_result['n_cycles_completed']} "
        f"max_drift_kms={negative_result['max_drift_kms']:.3e} "
        f"repeats={negative_result['repeats_to_machine_precision']} "
        f"(expected: False -- known one-off closure)",
        flush=True,
    )

    print("[575-C2] Cross-check: #571's 187 raw candidates vs symmetric condition...", flush=True)
    cross = cross_check_571_against_symmetric_condition(t_syn)
    print(
        f"[575-C2]   scanned={cross['n_571_records_scanned']} "
        f"hits_near_symmetric_condition={cross['n_hits_near_symmetric_condition']}",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#575 C2 two-sided repeat-instrumentation control",
                    "t_syn_days": t_syn,
                    "n_cycles": N_CYCLES,
                    "n_positive_survivors_tested": len(survivors),
                    "n_positive_repeat_to_machine_precision": n_positive_repeat,
                    "negative_control_repeats": negative_result["repeats_to_machine_precision"],
                }
            )
            + "\n"
        )
        for res in positive_results:
            fh.write(json.dumps({"kind": "positive_control", **res}) + "\n")
        fh.write(json.dumps({"kind": "negative_control", **negative_result}) + "\n")
        fh.write(json.dumps({"kind": "cross_check_571", **cross}) + "\n")
    print(f"[575-C2] written: {OUT_PATH}", flush=True)

    ok = n_positive_repeat == len(survivors) and not negative_result["repeats_to_machine_precision"]
    print(f"[575-C2] C2 DISCIPLINE {'PASSED' if ok else 'FAILED'}", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
