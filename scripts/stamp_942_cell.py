"""#942/#943: stamp one adjudicated enumeration cell into data/empty_regions.jsonl.

The stamp carries the exact method scope (the settings.json of the run, the pre-registration
in docs/notes/2026-10-05-942-943-two-working-body-generator.md sec. 6) and the gauntlet
verdict per gate-passing cycler. "Empty" (or "no new cycler") is conditional on that scope.

Usage: uv run python scripts/stamp_942_cell.py --gauntlet FILE --settings RUN/settings.json
        --region-id ID --family TEXT --verdict TEXT --interpretation TEXT --git-sha SHA
        [--dry-run]
"""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from cyclerfinder.data.empty_regions import EmptyRegionReport, append_empty_region
from cyclerfinder.data.method_capability import MethodCapability

REPO = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gauntlet", type=Path, required=True)
    ap.add_argument("--settings", type=Path, required=True)
    ap.add_argument("--region-id", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--verdict", required=True)
    ap.add_argument("--interpretation", required=True)
    ap.add_argument("--git-sha", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    g = json.loads(args.gauntlet.read_text())
    st = json.loads(args.settings.read_text())
    centre = "Jupiter" if g["cell"] in ("gc", "ge", "gc1", "ge1") else "Sun"
    report = EmptyRegionReport(
        region_id=args.region_id,
        family=args.family,
        centre=centre,
        topologies=(
            {
                "template": "[A-block, A->B, B-block, B->A], one visit per cycle",
                "k_synodic": st["k"],
                "max_returns_per_block": st["max_returns"],
                "transfer_revs": st["transfer_revs"],
                "generic_revs": st["generic_revs"],
                "resonant_only": st.get("resonant_only", ""),
                "catalogue": "resonant 1:1, 2:1, 1:2, 3:2, 2:3; half-rev (1,0,p/a), (3,1,p/a); "
                "generic same-body Lambert",
            },
        ),
        method_capability=MethodCapability(
            genome="two-working-body cycle templates (scripts/run_942_enumerate.py)",
            corrector="two_working_body.correct_dates (date residual) + minimax turn gate",
            capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
            git_sha=args.git_sha,
        ),
        search_extent={
            "points_total": g["n_structures"],
            "points_total_meaning": "cycle templates (structures), each seeded on a grid",
            "seeds": {
                "n_phase": st["n_phase"],
                "n_split": st["n_split"],
                "n_refine": st["n_refine"],
                "min_sep_frac": 0.03,
            },
            "ephem_model": "circular coplanar ideal model",
            "center": centre,
            "n_zeros": g["n_zeros"],
            "n_physical_cyclers": g["n_physical"],
        },
        prune_gates=(
            "exact zero: max |V_inf magnitude residual| < 1e-8 km/s",
            "#888/#937 demanded-turn gate pass at every massive flyby (registry floors)",
            "near-180 demanded turn >= 175 deg rejected (provisional)",
            "independent re-propagation miss < 1 km at every encounter",
        ),
        result={
            "n_gate_passing": g["n_gate_passing"],
            "gate_passing": [
                {
                    "k": c["k"],
                    "key": c["key"],
                    "vinf_kms": c["vinf_kms"],
                    "collisions": c["collisions"],
                    "r_min_km": c["r_min_km"],
                    "r_max_km": c["r_max_km"],
                    "dop853_miss_km": c["cross_check"]["max_arrival_miss_km"],
                }
                for c in g["candidates"]
            ],
        },
        verdict=args.verdict,
        interpretation=args.interpretation,
        source_anchors="docs/notes/2026-10-05-942-943-two-working-body-generator.md sec. 6; "
        f"data/942_cell_{g['cell']}_gauntlet.json",
        run={"date": date.today().isoformat(), "task": 942, "git_sha": args.git_sha},
    )
    if args.dry_run:
        print(report)
        return
    append_empty_region(REPO / "data" / "empty_regions.jsonl", report)
    print(f"stamped {args.region_id}")


if __name__ == "__main__":
    main()
