"""#576 step 4 -- multi-cycle repeat-instrumentation control for the 36 Galilean
symmetric-closure survivors found in
``data/enumerate_576_jupiter_galilean_symmetric_closures.jsonl`` (step 3).

Reuses #575's own C2 mechanism (``probe_575_titan_iapetus_repeat_check.py``)'s
underlying driver, not reimplemented: that script wraps the Saturn-specific
3D/eccentric ``run_v2_saturn_3d`` (built ON TOP OF the generic
``v2_moontour.run_v2_moontour`` "generalized ... to the 3D eccentric model").
The Galilean symmetric-closure construction here stays CIRCULAR-COPLANAR
(the same idealized model #563's own Uranian construction and #575's coplanar
Stage-0 construction used) -- there is no inclination/eccentricity correction
in play at this stage (that is out of scope per the #576 dispatch: no
real-ephemeris gauntlet run here). So the correct, already-generic tool is
``v2_moontour.run_v2_moontour`` ITSELF -- the same driver #574/#575's own
Saturn-specific wrapper is built on, and the SAME driver the Uranian #330/#558
gauntlet already used for this exact idealized model. No genericization
needed: ``run_v2_moontour``'s ``system`` argument auto-resolves the primary
from the moon names via ``core/satellites.py``.

Two sides, per #575's own C2 spec, mirrored here:

(a) POSITIVE -- every one of the 36 #576 symmetric-enumeration survivors must
    repeat to machine precision at every cycle (guaranteed BY CONSTRUCTION,
    not a new claim -- a failure here would flag a bug in the construction).
(b) NEGATIVE -- per the #576 dispatch instruction ("use a non-symmetric
    reject from your OWN new Jupiter enumeration output... pick any candidate
    that failed the symmetric gate"): constructed here directly from the SAME
    ``residual_at_point``/``gate_candidate`` machinery step 3 used, on the
    SAME pair/n_rev/commensurate-tof point as a genuine survivor, but at
    rel_offset=90 deg -- a value OUTSIDE the {0 deg, 180 deg} symmetric set
    the construction enumerates over, i.e. a point step 3's own enumeration
    loop never visits and which fails "the symmetric gate" (membership in
    the constructed set) by definition. This is sourced from THIS dispatch's
    own new Jupiter parameter sweep (same pair/n_rev/tof/primary), not from
    #501's real-ephemeris data (wrong model for this repeat check) or #571's
    unrelated Saturn data.

Discipline: NO catalogue writeback, reuses #563/#575's own gate/repeat
machinery verbatim, read-only on the step-3 enumeration output.

Run as::

    uv run python scripts/probe_576_galilean_repeat_check.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scan_558_uranus_all_pairs_offset_sweep import residual_at_point  # noqa: E402

from cyclerfinder.data.validation.v2_moontour import (  # noqa: E402
    V2MoontourVerdict,
    run_v2_moontour,
)

DATA_DIR = ROOT / "data"
ENUM_PATH = DATA_DIR / "enumerate_576_jupiter_galilean_symmetric_closures.jsonl"
OUT_PATH = DATA_DIR / "probe_576_galilean_repeat_check.jsonl"

N_CYCLES = 3
PRIMARY = "Jupiter"


def load_576_survivors() -> list[dict[str, Any]]:
    recs = []
    with ENUM_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "pass":
                recs.append(d)
    return recs


def _verdict_summary(v: V2MoontourVerdict) -> dict[str, Any]:
    return {
        "candidate_id": v.candidate_id,
        "n_cycles_requested": v.n_cycles_requested,
        "n_cycles_completed": v.n_cycles_completed,
        "max_drift_kms": v.max_drift_kms,
        "max_closure_residual_kms": v.max_closure_residual_kms,
        "passes_v2": v.passes_v2,
        "repeats_to_machine_precision": bool(
            v.n_cycles_completed == v.n_cycles_requested and v.max_closure_residual_kms < 1e-6
        ),
    }


def repeat_check_survivor(rec: dict[str, Any]) -> dict[str, Any]:
    seq = (rec["anchor"], rec["flyby"], rec["anchor"])
    tof = rec["tof_days"]
    cid = (
        f"{rec['anchor']}-{rec['flyby']}_n{rec['n_commensurate_int']}_"
        f"nrev{rec['n_rev']}_rel{rec['rel_offset_deg']:.0f}"
    )
    verdict = run_v2_moontour(
        cid,
        seq,
        tuple(rec["vinf_per_encounter_kms"]),
        (tof, tof),
        rec["rel_offset_deg"],
        None,  # primary auto-resolved from sequence (Jupiter)
        n_cycles=N_CYCLES,
        n_revs=tuple(rec["n_rev"]),
        phase0_deg=0.0,
        notes="#576 step 4 positive control (genuine symmetric-enumeration survivor)",
    )
    return _verdict_summary(verdict)


def build_negative_control(
    survivors: list[dict[str, Any]],
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    """The #576 dispatch's negative control: same pair/n_rev/commensurate-tof
    point as the richest genuine Ganymede-Callisto survivor, but rel_offset=90
    deg -- a value OUTSIDE the {0,180} symmetric set step 3's enumeration
    ever visits, sourced entirely from THIS dispatch's own Jupiter parameter
    sweep (never #501's real-ephemeris data or #571's Saturn data). Returns
    ``(neg_pt_or_None, negative_result)``."""
    base = max(
        (r for r in survivors if {r["anchor"], r["flyby"]} == {"Ganymede", "Callisto"}),
        key=lambda r: r["n_commensurate_int"],
    )
    off_symmetric_rel = 90.0
    tof_scale = base["tof_days"] / _sqrt_papb(base["anchor"], base["flyby"])
    neg_pt = residual_at_point(
        base["anchor"],
        base["flyby"],
        rel_offset_deg=off_symmetric_rel,
        tof_scale=tof_scale,
        n_rev=tuple(base["n_rev"]),
        primary=PRIMARY,
    )
    if neg_pt is None:
        negative_result: dict[str, Any] = {
            "candidate_id": f"{base['anchor']}-{base['flyby']}_rel90_NEGATIVE",
            "closure_at_cycle0": False,
            "reason": "Lambert infeasible at rel_offset=90deg for this n_rev branch",
            "repeats_to_machine_precision": False,
        }
        return None, negative_result

    cid = f"{base['anchor']}-{base['flyby']}_rel90_NEGATIVE"
    verdict = run_v2_moontour(
        cid,
        (base["anchor"], base["flyby"], base["anchor"]),
        _vinf_tuple(neg_pt),
        (neg_pt["tof_days"], neg_pt["tof_days"]),
        off_symmetric_rel,
        None,
        n_cycles=N_CYCLES,
        n_revs=tuple(base["n_rev"]),
        phase0_deg=0.0,
        notes="#576 step 4 negative control (non-symmetric rel_offset=90deg reject, "
        "from this dispatch's own Jupiter parameter sweep)",
    )
    negative_result = _verdict_summary(verdict)
    negative_result["cycle0_residual_kms"] = neg_pt["residual_kms"]
    return neg_pt, negative_result


def main() -> int:
    survivors = load_576_survivors()
    print(f"[576-repeat] (a) POSITIVE: {len(survivors)} #576 symmetric survivors...", flush=True)

    positive_results = []
    n_positive_repeat = 0
    for rec in survivors:
        res = repeat_check_survivor(rec)
        positive_results.append(res)
        if res["repeats_to_machine_precision"]:
            n_positive_repeat += 1
        print(
            f"[576-repeat]   {res['candidate_id']}: "
            f"n_cycles_completed={res['n_cycles_completed']} "
            f"max_closure_residual_kms={res['max_closure_residual_kms']:.3e} "
            f"repeats={res['repeats_to_machine_precision']}",
            flush=True,
        )
    print(
        f"[576-repeat] POSITIVE control: {n_positive_repeat}/{len(survivors)} repeat to "
        "machine precision",
        flush=True,
    )

    print(
        "[576-repeat] (b) NEGATIVE: same pair/n_rev/commensurate-tof point as a genuine "
        "survivor, but rel_offset=90 deg (outside the {0,180} symmetric set)...",
        flush=True,
    )
    neg_pt, negative_result = build_negative_control(survivors)
    if neg_pt is None:
        print(
            f"[576-repeat]   NEGATIVE control: {negative_result['reason']} "
            "-- trivially does not repeat/close.",
            flush=True,
        )
    else:
        print(
            f"[576-repeat]   NEGATIVE control {negative_result['candidate_id']}: "
            f"cycle0_residual={neg_pt['residual_kms']:.4f} km/s "
            f"n_cycles_completed={negative_result['n_cycles_completed']} "
            f"repeats={negative_result['repeats_to_machine_precision']} "
            "(expected: False -- non-symmetric point)",
            flush=True,
        )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#576 step 4 -- Galilean multi-cycle repeat-instrumentation control",
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
    print(f"[576-repeat] written: {OUT_PATH}", flush=True)

    ok = n_positive_repeat == len(survivors) and not negative_result["repeats_to_machine_precision"]
    print(f"[576-repeat] C2-STYLE DISCIPLINE {'PASSED' if ok else 'FAILED'}", flush=True)
    return 0 if ok else 1


def _sqrt_papb(anchor: str, flyby: str) -> float:
    import math

    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
    from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day

    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    return math.sqrt(p_a * p_b)


def _vinf_tuple(pt: dict[str, Any]) -> tuple[float, float, float]:
    """Per-encounter V_inf magnitude, matching enumerate_563's own
    ``encounter_vinfs_kms`` extraction (max of in/out asymptote per node)."""
    vin = pt["vinf_in"]
    vout = pt["vinf_out"]
    return tuple(max(abs(vin[k]), abs(vout[k])) for k in range(3))


if __name__ == "__main__":
    raise SystemExit(main())
