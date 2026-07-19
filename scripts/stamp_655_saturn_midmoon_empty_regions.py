"""#655 -- stamp the 3 clean, unambiguous 0/N Saturn mid-moon-pair symmetric-closure
negatives (Dione-Rhea, Tethys-Dione, Enceladus-Tethys) to ``data/empty_regions.jsonl``.

Rhea-Titan is DELIBERATELY EXCLUDED here: the base #563-method enumeration found 3
gate-passing coplanar symmetric closures for that pair
(``data/enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl``). Per the #655
dispatch's explicit discipline, any pair with base-enumeration gate-passing closures
is NOT stamped/closed by this dispatch -- it is documented in full (including the
#655 inclination-extension + multi-cycle repeat-check result) and held for the
coordinating session's own decision on any literature-check/adjudication follow-up.
This script only touches the 3 pairs with an unambiguous 0/192, 0/128, 0/64 result
at the base gate -- no adjudication judgment call is needed for these three.

Run as::

    uv run python scripts/stamp_655_saturn_midmoon_empty_regions.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.empty_regions import (  # noqa: E402
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402

DATA_DIR = ROOT / "data"


def _git_sha() -> str:
    out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT)
    return out.strip()


def _load_meta_and_directions(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    meta = None
    dirs = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("_meta"):
                meta = d
            elif d.get("kind") == "direction_summary":
                dirs.append(d)
    assert meta is not None
    return meta, dirs


PAIRS: tuple[tuple[str, str, str], ...] = (
    (
        "dione-rhea",
        "Dione",
        "Rhea",
    ),
    (
        "tethys-dione",
        "Tethys",
        "Dione",
    ),
    (
        "enceladus-tethys",
        "Enceladus",
        "Tethys",
    ),
)


def main() -> int:
    sha = _git_sha()
    method = MethodCapability(
        genome=(
            "#563 direct symmetric-closure enumeration (rel_offset in {0,180deg}, "
            "exactly-commensurate tof = n*T_syn/2, matched leg-revolution counts -- the "
            "classical symmetric perpendicular-crossing construction), genericized to "
            "primary=Saturn (scripts/enumerate_563_symmetric_closures.py --primary Saturn "
            "--moons <pair>), same tool #575/#576/#599 already genericized and "
            "bit-for-bit validated against the 30-closure Uranian golden."
        ),
        corrector=(
            "residual_at_point/gate_candidate (scan_558_uranus_all_pairs_offset_sweep.py, "
            "reused verbatim) -- circular-coplanar Kepler Lambert-leg closure + #324 "
            "two-sided physical max-bend gate + DOP853 independent cross-check."
        ),
        capability_tags=frozenset(
            {
                "ballistic",
                "patched-conic",
                "symmetric-closure-construction",
                "circular",
                "coplanar",
                "single-arc",
            }
        ),
        git_sha=sha,
    )

    for region_slug, anchor, flyby in PAIRS:
        fname = f"enumerate_655_saturn_{anchor.lower()}_{flyby.lower()}_symmetric_closures.jsonl"
        path = DATA_DIR / fname
        meta, dirs = _load_meta_and_directions(path)
        n_max = dirs[0]["n_max"]
        t_syn = dirs[0]["t_syn_days"]
        n_infeasible_total = sum(d["n_infeasible"] for d in dirs)
        n_subgate_total = sum(d["n_subgate_residual_only"] for d in dirs)

        report = EmptyRegionReport(
            region_id=f"saturn-{region_slug}-symmetric-closure-empty-655",
            family=(
                f"{anchor}-{flyby}-{anchor} symmetric-closure quasi_cyclers (#563-method direct "
                "construction) -- the finite, exhaustive symmetric-closure enumeration finds "
                "ZERO candidates that clear the two-sided #324 physical bend gate at BOTH "
                "encounters simultaneously, despite a healthy sub-gate residual-closure population"
            ),
            centre="Saturn",
            topologies=(
                {"sequence": [anchor, flyby, anchor], "period_k": 1, "n_rev_range": [0, 3]},
            ),
            method_capability=method,
            search_extent={
                "points_total": meta["total_evaluated"],
                "n_epochs": 0,
                "span_days": 0.0,
                "ephem_model": "circular-coplanar Kepler (idealized construction)",
                "center": "Saturn",
                "directions": meta["directions"],
                "n_rev_combos": meta["n_rev_combos"],
                "rel_offsets_deg": meta["rel_offsets_deg"],
                "tof_scale_max_bound": meta["tof_scale_max_bound"],
                "n_max": n_max,
                "t_syn_days": t_syn,
                "n_lambert_infeasible": n_infeasible_total,
                "n_subgate_residual_only": n_subgate_total,
            },
            prune_gates=(
                "#563 symmetric-closure construction: rel_offset in {0,180deg} AND "
                "exactly-commensurate tof = n*T_syn/2 AND matched leg-revolution counts "
                "(perpendicular-crossing mirror-symmetry condition)",
                "residual_kms < GATE_RESIDUAL_KMS=0.05 km/s (Lambert V_inf-continuity closure)",
                "two-sided #324 physical max-bend gate (search/physical_sanity.py::"
                "candidate_passes_physical_gate), every encounter clears "
                "DEFAULT_MIN_USEFUL_BEND_DEG=5.0 deg",
                "DOP853 independent-integrator cross-check",
            ),
            result={
                "n_evaluated": meta["total_evaluated"],
                "n_lambert_infeasible": n_infeasible_total,
                "n_subgate_residual_only": n_subgate_total,
                "n_all_gates_passed": meta["total_all_gates_passed"],
                "failure_mode": (
                    "residual gate passes for a meaningful subset "
                    f"({n_subgate_total}/{meta['total_evaluated']}), but every one of those "
                    "candidates has at least one of its two encounters below the 5deg bend "
                    "floor -- the small tof_scale_max=3.0-bounded n_max (1-3 for these "
                    "close-in, short-period pairs) leaves too few commensurate points for "
                    "any of them to land on a V_inf combination that bends >=5deg at BOTH "
                    "encounters simultaneously (unlike Rhea-Titan's wider period-ratio "
                    "spacing, where 3/512 candidates DO)."
                ),
            },
            verdict="EMPTY",
            interpretation=(
                f"Zero genuine {anchor}-{flyby} symmetric-closure quasi_cyclers exist under the "
                "#563 direct-construction method within the tof_scale_max=3.0 bound (the same "
                "bound #558's own Uranian production sweep and #575's Titan-Iapetus sweep used). "
                "This is a real, gate-level negative, not a coverage gap: the finite enumeration "
                "is EXHAUSTIVE within the bound (every commensurate n/n_rev/rel_offset point is "
                "directly constructed and evaluated, nothing is grid-sampled or interpolated), so "
                "there is no basin this run could have missed inside that bound. Root cause "
                f"(#654's own live feasibility check correctly identified WIDE two-sided bend "
                "windows exist at the Hohmann-FLOOR V_inf for this pair, but the actual finite "
                "symmetric-closure family does not land its V_inf combinations anywhere near "
                "that floor -- the discrete n/n_rev grid only realizes a sparse set of specific "
                "V_inf pairs, and for this short-period, small-n_max pair none of them "
                "simultaneously bend >=5deg at both encounters). This is a materially different, "
                "and in one sense EARLIER, failure mode than #575's Titan-Iapetus result (which "
                "found a genuine coplanar family that then failed only at the inclination-"
                "extension stage): here the idealized COPLANAR model itself is empty at the "
                "physical-bend gate, before any inclination correction is even introduced."
            ),
            source_anchors=(
                "scripts/enumerate_563_symmetric_closures.py (#563 symmetric construction, "
                "genericized to primary=Saturn, C1-golden-validated bit-for-bit vs the "
                "30-closure Uranian table); core/satellites.py Dione/Tethys/Enceladus/Rhea "
                "GM/sma (JPL SSD gm_de440/phys_par+mean elements). "
                f"Data: {path.relative_to(ROOT)}. "
                "Companion pair Rhea-Titan (NOT stamped here, see #655's own OUTSTANDING.md "
                "bullet): 3/512 candidates pass this same base gate but 0/3 survive the "
                "#655 inclination-extension + multi-cycle repeat check "
                "(data/probe_655_rhea_titan_3d_closure.jsonl, "
                "data/probe_655_rhea_titan_repeat_check.jsonl)."
            ),
            run={
                "date": "2026-07-19",
                "cores": 1,
                "git_sha": sha,
                "task": 655,
                "note": (
                    f"{anchor}-{flyby} clean 0/{meta['total_evaluated']} negative -- no "
                    "adjudication judgment call needed (unlike the companion Rhea-Titan pair)."
                ),
            },
        )
        append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
        print(f"[655] stamped: {report.region_id}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
