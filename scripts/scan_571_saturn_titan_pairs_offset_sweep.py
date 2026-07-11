"""#571 -- Titan-centered Saturnian moon-pair sweep via the #558-corrected genome.

Applies the SAME discovery genome that found #312's entire 30-member Uranian
family (``scripts/scan_558_uranus_all_pairs_offset_sweep.py``'s relative-offset
x global-phase x tof_scale x n_rev sweep -- NOT holding relative offset fixed
by convention, the exact bug that nearly hid #312 itself) to a Saturnian pair
it has never touched: Titan-Rhea and Titan-Iapetus.

This is a genuine re-application, not a rebuild: #558's own ``sweep_pair`` /
``residual_at_point`` / ``gate_candidate`` functions are imported and reused
VERBATIM (the residual formula and the #324 physical-bend + DOP853 gate are
untouched); the only change made to that module (see its own #571 diff) is
adding an explicit ``primary: str = "Uranus"`` keyword so the SAME functions
can be pointed at ``PRIMARIES["Saturn"]`` / Saturn's moons instead of
re-deriving the Kepler-state / Lambert / residual logic here.

SCOPE (per the 2026-07-12 Fable plan-review correction to #571 in
``data/OUTSTANDING.md`` -- read that entry, this docstring only summarizes):

* Titan-Mimas / Titan-Enceladus / Titan-Tethys / Titan-Dione are NOT swept
  here -- Fable's independent max-bend-ceiling-vs-min-achievable-V_inf
  analysis (re-verified independently in this task, see
  ``scripts/verify_571_gate_analytics.py``... actually inline in the #571
  OUTSTANDING.md RESULT note) showed all four are analytically empty under
  this project's own two-sided #324 gate at EVERY grid point -- stamped
  directly into ``data/empty_regions.jsonl`` instead of swept.
* Titan-Hyperion is excluded (GM=0.37, mass-deficient, matches the same
  profile -- also confirmed analytically infeasible, see the same
  cross-check: ceiling 0.237 km/s vs min-achievable 0.249 km/s, 4.53 deg
  bend at minimum V_inf, just under the 5 deg floor).
* Titan-Rhea is V0-known (#489's existing adjudication) -- any basin found
  here is a CENSUS entry, not a novelty claim.
* Titan-Iapetus IS novelty-eligible if a real hit survives, BUT Iapetus's
  real ~15.5 deg inclination to Titan's plane means a coplanar closure here
  is NOT on its own a physical statement -- any sub-gate Titan-Iapetus basin
  MUST be flagged "coplanar-idealized only, needs 3D/inclined treatment"
  before being treated as a genuine candidate, exactly as any #558-style
  planar hit would NOT be.

Positive control (INTERNAL cross-lineage smoke test, NOT a literature
golden -- Russell-Strange 2009 published Titan-Enceladus data only, never
Titan-Rhea; per [[feedback_golden_tests_sourced_only]]): this project's own
#320 epoch-aware quasi_cycler scan already found a Titan-Rhea-Titan basin
(``data/scan_320_epoch_aware_saturn.jsonl`` idx 289: residual 0.031620 km/s,
V_inf 1.6769-1.7528 km/s, bends 49.37/6.83/50.49 deg, at rel_offset_deg=285,
phase0_deg=90, tof_scale=2.0, n_rev=(1,1)). This script's genuinely
DIFFERENT method (circular-coplanar-Kepler + patched-conic-Lambert, vs
#320's own epoch-aware real-ephemeris-adjacent construction) is checked for
whether it recovers a MATCHING/CONSISTENT basin near this same signature --
an independent cross-lineage consistency check, not a reproduction of the
identical method.

Discipline: NO catalogue writeback, NO V1-V4-strict gauntlet, NO Titan-Mimas/
Enceladus/Tethys/Dione/Hyperion sweep. Reuses #558's own gate machinery
verbatim.

Run as::

    uv run python scripts/scan_571_saturn_titan_pairs_offset_sweep.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

# Reuse #558's own gate machinery verbatim (residual formula + #324 physical
# bend gate + DOP853 cross-check) -- no new gate logic, only the `primary=`
# genericization added directly to that module for this task.
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    N_REV_MAX,
    _git_sha,
    gate_candidate,
    residual_at_point,
    sweep_pair,
)

PRIMARY = "Saturn"
PAIR_DIRECTIONS: tuple[tuple[str, str], ...] = (
    ("Titan", "Rhea"),
    ("Rhea", "Titan"),
    ("Titan", "Iapetus"),
    ("Iapetus", "Titan"),
)

DATA_DIR = ROOT / "data"

# #320's own Titan-Rhea-Titan basin (data/scan_320_epoch_aware_saturn.jsonl,
# idx 289) -- the INTERNAL cross-lineage positive control, not a literature
# golden. See module docstring.
INTERNAL_ANCHOR = {
    "sequence": ["Titan", "Rhea", "Titan"],
    "rel_offset_deg": 285.0,
    "phase0_deg": 90.0,
    "tof_scale": 2.0,
    "n_rev": (1, 1),
    "residual_kms": 0.03161954212289819,
    "vinf_per_encounter_kms": [1.752838353627724, 1.6769016336877292, 1.721218811504826],
    "max_bend_deg_per_enc": [49.37110969924298, 6.833550698112894, 50.493642373507285],
    "source": "data/scan_320_epoch_aware_saturn.jsonl idx 289 (#320 epoch-aware sweep)",
}


def main() -> int:
    sha = _git_sha()
    print(f"[571] Saturn Titan-Rhea/Titan-Iapetus offset sweep -- sha={sha}", flush=True)

    # --- Internal cross-lineage positive control -------------------------
    print(
        "[571] INTERNAL cross-lineage check (NOT a literature golden -- see module "
        "docstring): does this DIFFERENT method recover a basin near #320's own "
        "Titan-Rhea-Titan signature (rel_off=285, tof_scale=2.0, n_rev=(1,1), "
        "residual=0.031620 km/s)?",
        flush=True,
    )
    pt = residual_at_point(
        "Titan",
        "Rhea",
        rel_offset_deg=INTERNAL_ANCHOR["rel_offset_deg"],
        tof_scale=INTERNAL_ANCHOR["tof_scale"],
        n_rev=INTERNAL_ANCHOR["n_rev"],
        phase0_deg=INTERNAL_ANCHOR["phase0_deg"],
        primary=PRIMARY,
    )
    if pt is not None:
        print(
            f"[571]   exact-point (this method) residual = {pt['residual_kms']:.6f} km/s "
            f"vs #320's own 0.031620 km/s at the same (rel_off, tof_scale, n_rev)",
            flush=True,
        )
    else:
        print("[571]   exact-point evaluation INFEASIBLE (no Lambert solution)", flush=True)

    pc_out = DATA_DIR / "scan_571_internal_positive_control.jsonl"
    with pc_out.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#571 internal cross-lineage positive control -- NOT a "
                    "literature golden (see module docstring)",
                    "git_sha": sha,
                    "internal_anchor": INTERNAL_ANCHOR,
                }
            )
            + "\n"
        )
        if pt is not None:
            fh.write(json.dumps({"kind": "exact_point_this_method", **pt}) + "\n")
    print(f"[571] internal positive control written to {pc_out}", flush=True)

    # --- Titan-Rhea / Titan-Iapetus sweep (both directions each) ---------
    index_rows: list[dict[str, Any]] = []
    candidates_needing_adjudication: list[dict[str, Any]] = []
    tof_scales = tuple(round(0.5 + 0.05 * k, 2) for k in range(51))  # matches #558 'dense'

    for anchor, flyby in PAIR_DIRECTIONS:
        slug = f"{anchor.lower()}_{flyby.lower()}"
        out_path = DATA_DIR / f"scan_571_saturn_{slug}.jsonl"
        t0 = time.time()
        result = sweep_pair(
            anchor,
            flyby,
            n_phase=1,  # phase0 is provably redundant for fixed rel_offset/tof_scale (#558 proof)
            n_offset=360,
            tof_scales=tof_scales,
            n_rev_max=N_REV_MAX,
            keep_top=200,
            primary=PRIMARY,
        )
        elapsed = time.time() - t0
        best = result["best_overall_residual_kms"]
        print(
            f"[571] {anchor}-{flyby}-{anchor}: best={best} n_records={result['n_records_total']} "
            f"n_gate_passing={result['n_gate_passing_records']}  ({elapsed:.1f}s)",
            flush=True,
        )

        gated: list[dict[str, Any]] = []
        for rec in result["gate_passing_records"]:
            g = gate_candidate(anchor, flyby, rec, primary=PRIMARY)
            gated.append(g)
            if g["all_gates_passed"]:
                if "Iapetus" in (anchor, flyby):
                    g["iapetus_inclination_caveat"] = (
                        "COPLANAR-IDEALIZED ONLY. Iapetus's real inclination (~15.5 deg to "
                        "Titan's orbital plane) gives an out-of-plane relative velocity "
                        "(~0.85 km/s) comparable to the entire coplanar V_inf budget at this "
                        "pair (~0.93 km/s minimum-achievable) -- a coplanar closure at the "
                        "0.05 km/s gate floor is NOT a physical statement for Titan-Iapetus on "
                        "its own. This candidate needs 3D/inclined treatment (the #552-scoped "
                        "genome extension, not yet built) before being treated as a genuine "
                        "candidate; it is NOT equivalent in confidence to a #558-style planar "
                        "Uranian hit and must not be adjudicated as if it were."
                    )
                candidates_needing_adjudication.append(g)

        with out_path.open("w", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "_meta": True,
                        "task": "#571 Saturn Titan-Rhea/Titan-Iapetus offset sweep",
                        "anchor": anchor,
                        "flyby": flyby,
                        "sequence": [anchor, flyby, anchor],
                        "primary": PRIMARY,
                        "n_phase": result["n_phase"],
                        "n_offset": result["n_offset"],
                        "tof_scales": result["tof_scales"],
                        "n_rev_max": result["n_rev_max"],
                        "elapsed_s": elapsed,
                        "git_sha": sha,
                    }
                )
                + "\n"
            )
            fh.write(json.dumps({"kind": "best_overall", **result["best_overall_record"]}) + "\n")
            for rec in result["top_records"]:
                fh.write(json.dumps({"kind": "top_record", **rec}) + "\n")
            for g in gated:
                fh.write(json.dumps({"kind": "gate_result", **g}) + "\n")

        index_rows.append(
            {
                "anchor": anchor,
                "flyby": flyby,
                "path": str(out_path.relative_to(ROOT)),
                "best_overall_residual_kms": best,
                "n_gate_passing_records": result["n_gate_passing_records"],
                "n_all_gates_passed": sum(1 for g in gated if g["all_gates_passed"]),
                "elapsed_s": elapsed,
            }
        )

    index_out = DATA_DIR / "scan_571_saturn_titan_pairs_index.jsonl"
    with index_out.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#571 Saturn Titan-Rhea/Titan-Iapetus offset sweep -- index",
                    "git_sha": sha,
                    "n_pair_directions": len(PAIR_DIRECTIONS),
                    "gate_residual_kms": GATE_RESIDUAL_KMS,
                    "excluded_pairs_stamped_analytically_empty": [
                        "Titan-Mimas",
                        "Titan-Enceladus",
                        "Titan-Tethys",
                        "Titan-Dione",
                    ],
                    "excluded_titan_hyperion_reason": "GM=0.37, mass-deficient, confirmed "
                    "analytically infeasible (ceiling 0.237 km/s < min-achievable would need "
                    "to be, actual bend at min-achievable V_inf 4.53 deg < 5 deg floor)",
                }
            )
            + "\n"
        )
        for row in index_rows:
            fh.write(json.dumps(row) + "\n")
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "candidates_needing_adjudication",
                    "count": len(candidates_needing_adjudication),
                    "candidates": candidates_needing_adjudication,
                }
            )
            + "\n"
        )
    print(f"[571] DONE -- index: {index_out}", flush=True)
    print(
        f"[571] candidates passing ALL gates (residual+physical+DOP853): "
        f"{len(candidates_needing_adjudication)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
