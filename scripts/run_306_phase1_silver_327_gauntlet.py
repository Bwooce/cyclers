"""#306 Phase 1 — run the new V1 + V2 3D gauntlets on the #327 SILVER.

Produces ``data/silver_327_v1_v2_verdicts.jsonl`` with two records:

  * V1 verdict — same-model closure of the Lambert moontour, reinterpreted
    against the spec §14 1 m/s independent-cross-check floor (the dop853
    re-propagation arrival residual from the SILVER's own verification run).
  * V2 verdict — long-span bounded-drift report.

The #327 SILVER is a multi-leg Lambert moontour (Uranus-Oberon-Umbriel,
patched-conic in CR3BP frame), NOT a single 6D periodic-orbit IC. The V1/V2
modules built in Phase 1 Part A/B operate on a CR3BP periodic IC + period.
This script does NOT shoehorn the moontour into the periodic-orbit
interpretation; instead it reports the SILVER's existing same-model numbers
under the V1/V2 *philosophy* (closure vs independent re-propagation, bounded
drift) and explicitly flags where the moontour case differs from the CR3BP
periodic-orbit case.

Discipline
----------
* NO catalogue writeback.
* The verdict explicitly notes the moontour-vs-periodic-IC distinction so
  no later reader confuses the SILVER's V1 PASS with a CR3BP periodic-orbit
  V1 PASS.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.data.validation.v1_3d import V1_FLOOR_KMS
from cyclerfinder.data.validation.v2_3d import V2_DRIFT_FLOOR_KMS, V2_N_CYCLES_MIN

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SILVER_JSONL = _REPO_ROOT / "data" / "silver_327_verified.jsonl"
_SOURCE_JSONL = _REPO_ROOT / "data" / "scan_312_uranus_oberon_umbriel.jsonl"
_OUTPUT_JSONL = _REPO_ROOT / "data" / "silver_327_v1_v2_verdicts.jsonl"


def _load_silver_meta_row() -> dict[str, Any]:
    """Read the SILVER's meta header row (first line)."""
    with _SILVER_JSONL.open() as f:
        first = json.loads(f.readline())
    if not first.get("_meta"):
        raise RuntimeError(f"first row of {_SILVER_JSONL} is not a _meta header")
    return first


def _find_silver_source_row() -> dict[str, Any]:
    """Read the SILVER's source row from the original scan JSONL."""
    target_id = "repeated-moon-uranus-00000041"
    with _SOURCE_JSONL.open() as f:
        for raw in f:
            row = json.loads(raw)
            if row.get("candidate_id") == target_id:
                return row
    raise RuntimeError(f"candidate {target_id!r} not found in {_SOURCE_JSONL}")


def _git_sha() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=_REPO_ROOT)
        return out.decode().strip()
    except Exception:
        return "unknown"


def _build_v1_verdict(meta: dict[str, Any], source: dict[str, Any]) -> dict[str, Any]:
    """Build the V1-style verdict for the moontour SILVER.

    Honest framing:
      * The Lambert-leg patched-conic closure (``residual_kms``) is the
        analogue of the corrector's residual but is NOT the V1 spec quantity
        (which is solver-disagreement, not closure).
      * The dop853 re-propagation arrival residuals (``cross_check`` block)
        are the spec V1 independent-cross-check quantities — they MUST pass
        the 1 m/s floor for a V1 PASS, and they do (2.3e-10 km/s for the
        worst leg, 7+ orders of magnitude below the floor).
    """
    cross = source.get("cross_check", {})
    per_leg = cross.get("per_leg", [])
    max_dr_arrival_km = float(cross.get("max_dr_arrival_km", float("nan")))
    max_dv_arrival_kms = float(cross.get("max_dv_arrival_km_s", float("nan")))
    lambert_residual_kms = float(source.get("residual_kms", float("nan")))

    # Convert max_dr_arrival to nondim using Uranus-Oberon scale (Oberon
    # SMA ~ 583,520 km — the larger of the two moons in the sequence; this
    # is purely a reporting unit, the V1 spec floor is in km/s).
    oberon_sma_km = 583520.0
    max_dr_arrival_nondim = max_dr_arrival_km / oberon_sma_km

    # V1 PASS predicate: dop853 cross-check arrival velocity residual must
    # be below the spec §14 1 m/s floor on every leg.
    passes_v1 = bool(max_dv_arrival_kms <= V1_FLOOR_KMS) and bool(cross.get("all_passed", False))

    return {
        "kind": "v1_verdict_3d_moontour",
        "task": "#306 Phase 1 Part D — V1 verdict on #327 SILVER",
        "candidate_id": meta["candidate_id"],
        "primary": meta["primary"],
        "sequence": meta["sequence"],
        "n_rev": meta["n_rev"],
        # Inputs (sourced from the SILVER's own verification record).
        "lambert_closure_kms": lambert_residual_kms,
        "lambert_closure_note": (
            "Lambert leg-end-mismatch (the moontour's natural closure metric); "
            "not the V1 spec quantity, which is independent-cross-check arrival "
            "residual"
        ),
        "independent_max_dr_arrival_km": max_dr_arrival_km,
        "independent_max_dr_arrival_nondim": max_dr_arrival_nondim,
        "independent_max_dv_arrival_kms": max_dv_arrival_kms,
        "independent_cross_check_per_leg": per_leg,
        "v1_floor_kms": V1_FLOOR_KMS,
        "v1_floor_source": "spec §14 V1 — agreement < 1e-3 m/s",
        "passes_v1": passes_v1,
        "headline_kms": max_dv_arrival_kms,
        # Caveats — read these before quoting the PASS.
        "model_distinction_note": (
            "The SILVER is a multi-leg Lambert moontour (patched-conic, 2 legs, "
            "Uranus-moon CR3BP frame), NOT a 6D periodic CR3BP orbit. The V1_3D "
            "module shipped in Phase 1 Part A operates on (state0, period) for a "
            "periodic orbit. The SILVER's V1 verdict here uses the moontour's own "
            "independent cross-check residual (DOP853 re-propagation of each leg "
            "from the Lambert-derived velocities) against the same spec V1 km/s "
            "floor. A V1 PASS here means 'independent re-propagation arrived at "
            "the next encounter to < 1 m/s', NOT 'the orbit re-closes as a "
            "periodic IC under the 3D corrector'."
        ),
        "writeback_to_catalogue": False,
    }


def _build_v2_verdict(meta: dict[str, Any], source: dict[str, Any]) -> dict[str, Any]:
    """Build the V2-style verdict for the moontour SILVER.

    Honest framing: the SILVER is a single-cycle moontour solve. The V2
    long-span bounded-drift gate (spec §14 V2-ballistic: >=3 continuous laps)
    requires propagating ``n_cycles`` consecutive periods of a PERIODIC
    orbit. A single-cycle moontour has no such ``n_cycles`` semantics —
    we cannot run V2 in the same way.

    Three honest options:
      (a) Run V2 by *replicating* the moontour 3 times (3 consecutive
          Umbriel-Oberon-Umbriel sequences). Requires the Lambert solver to
          re-converge each cycle's legs, with phase advanced by the moons'
          orbital periods — non-trivial.
      (b) Treat the moontour as the FIRST CYCLE of a hypothetical
          periodic-orbit IC (its return state must equal its initial state).
          But the SILVER's stored ``residual_kms = 0.025`` is precisely the
          extent to which the SILVER FAILS this — return state differs from
          initial state by 25 m/s after one cycle.
      (c) Report V2 as NOT_RUN with a documented reason and flag this as
          Phase 2's V3 / multi-cycle moontour work.

    We take (c) and report the closure-residual (cycle-1 'drift') as a
    DIAGNOSTIC: extrapolating linearly to 3 cycles gives ~3 * 25 m/s in
    velocity space, which the moontour Lambert framework reports as the
    'returns' fidelity. That extrapolation is HEURISTIC; a faithful V2
    requires Phase 2's moontour V2 implementation (re-solving Lambert legs
    each cycle, allowing phase to evolve, gating on bounded position
    drift). Phase 1 explicitly excludes that.
    """
    lambert_residual_kms = float(source.get("residual_kms", float("nan")))
    tof_days = source.get("tof_days", [])
    one_cycle_days = sum(float(t) for t in tof_days) if tof_days else None

    return {
        "kind": "v2_verdict_3d_moontour",
        "task": "#306 Phase 1 Part D — V2 verdict on #327 SILVER",
        "candidate_id": meta["candidate_id"],
        "primary": meta["primary"],
        "sequence": meta["sequence"],
        "n_rev": meta["n_rev"],
        "one_cycle_tof_days": one_cycle_days,
        "lambert_one_cycle_velocity_residual_kms": lambert_residual_kms,
        "v2_n_cycles_min": V2_N_CYCLES_MIN,
        "v2_drift_floor_kms": V2_DRIFT_FLOOR_KMS,
        "v2_drift_floor_source": (
            "propagate.DRIFT_TOLERANCE_KM — V2-ballistic same-model 50,000 km"
        ),
        "v2_status": "NOT_RUN_PHASE_1",
        "v2_status_reason": (
            "The Phase 1 V2_3D module operates on a CR3BP periodic-orbit IC "
            "(6D state + period). The #327 SILVER is a single-cycle Lambert "
            "moontour with no closed periodic IC — its stored 25 m/s velocity "
            "residual at the end of cycle 1 is precisely the extent to which "
            "the moontour is NOT periodic. A faithful V2 for moontours requires "
            "re-solving the Lambert legs over 3 consecutive cycles with phase "
            "advanced by moon orbital periods; that is Phase 2 work, explicitly "
            "out of scope for Phase 1."
        ),
        "diagnostic_extrapolated_drift_kms_per_cycle": lambert_residual_kms,
        "diagnostic_note": (
            "lambert_one_cycle_velocity_residual_kms is the moontour's natural "
            "'V2 cycle-1 drift' analogue in km/s velocity space. Heuristic "
            "extrapolation to 3 cycles (linear in the absence of maintenance) "
            "would be ~3 * 25 m/s in velocity space. This is NOT a position "
            "drift in km and is NOT the V2 spec quantity — Phase 2 must "
            "implement the moontour-specific V2."
        ),
        "passes_v2": False,
        "writeback_to_catalogue": False,
    }


def main() -> None:
    meta = _load_silver_meta_row()
    source = _find_silver_source_row()
    v1 = _build_v1_verdict(meta, source)
    v2 = _build_v2_verdict(meta, source)

    header = {
        "_meta": True,
        "task": "#306 Phase 1 Part D — V1+V2 gauntlet verdicts on #327 SILVER",
        "phase": "phase1_v1_v2_3d_gauntlet",
        "candidate_id": meta["candidate_id"],
        "primary": meta["primary"],
        "sequence": meta["sequence"],
        "n_rev": meta["n_rev"],
        "tof_days": source.get("tof_days"),
        "source_silver_jsonl": str(_SILVER_JSONL.relative_to(_REPO_ROOT)),
        "source_scan_jsonl": str(_SOURCE_JSONL.relative_to(_REPO_ROOT)),
        "git_sha": _git_sha(),
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "discipline": (
            "V1 PASS does NOT admit to catalogue.yaml. Catalogue admission "
            "still blocked by #306 Phase 2-5 (V3+V4+V5) and the #328 lit-check "
            "verdict. Moontour-vs-CR3BP-periodic-IC distinction documented "
            "explicitly in each verdict's notes."
        ),
    }

    with _OUTPUT_JSONL.open("w") as f:
        f.write(json.dumps(header) + "\n")
        f.write(json.dumps(v1) + "\n")
        f.write(json.dumps(v2) + "\n")

    v1_head = v1["headline_kms"]
    v2_resid = v2["lambert_one_cycle_velocity_residual_kms"]
    print(f"V1 verdict: passes_v1={v1['passes_v1']} (headline {v1_head:.3e} km/s)")
    print(f"V2 verdict: status={v2['v2_status']} (one-cycle residual {v2_resid:.3e} km/s)")
    print(f"wrote {_OUTPUT_JSONL.relative_to(_REPO_ROOT)} (3 records)")


if __name__ == "__main__":
    main()
