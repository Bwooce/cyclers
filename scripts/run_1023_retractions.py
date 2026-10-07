"""#1023 (c): append retraction re-stamps for the n-body part of the seven void Jovian stamps.

Lead ruling 2026-10-08: the #318 / #501 real-ephemeris n-body stages are NOT RE-RUNNABLE until
#1039 (open real-ephemeris chain formulation). Their "0 closed" came from ``jovian_shoot`` with a
translation-only periodicity wrap (#968 note sec. 2.1) on a real configuration that does not repeat.
Each stamp gets an APPEND-ONLY companion record (the registry is never edited): the same region,
the patched-conic prefilter result kept, the n-body capability tags removed, and a
retraction verdict.
Idempotent: a region already re-stamped is skipped.

    uv run python scripts/run_1023_retractions.py
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]
VOID = (
    "jovian-cgcec-sobol-smoke-318-2026-06-30",
    "jovian-ege-sobol-broadened-501-2026-06-30",
    "jovian-gcg-sobol-broadened-501-2026-06-30",
    "jovian-egce-sobol-broadened-501-2026-06-30",
    "jovian-iei-sobol-broadened-501-2026-06-30",
    "jovian-iegi-sobol-broadened-501-2026-06-30",
    "jovian-egcge-sobol-broadened-501-2026-06-30",
)
NBODY_TAGS = {"n-body-shoot", "analytic-stm"}
REASON = (
    "N-BODY STAGE RETRACTED (#1023, lead ruling 2026-10-08): its 0 closures came from jovian_shoot "
    "with a translation-only periodicity wrap (fixed in b9c27f77; #968 note sec. 2.1) on a real "
    "ephemeris whose configuration does not repeat, so the n-body stage could not have closed a "
    "cycler; NOT RE-RUNNABLE until #1039 (open real-ephemeris chain formulation). The "
    "patched-conic prefilter result of the original stamp stands; this record carries "
    "only that capability."
)


def main() -> None:
    sha = subprocess.run(
        ["git", "rev-parse", "--short=8", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    preflight_search(
        task_no=1023,
        region_id="jovian-318-501-nbody-retractions",
        method=MethodCapability(
            genome="retraction re-stamps (no search)",
            corrector="none",
            capability_tags=frozenset({"patched-conic"}),
            git_sha=sha,
        ),
        script_path=Path(__file__),
        n_points=len(VOID),
    )
    text = DEFAULT_EMPTY_REGIONS_PATH.read_text()
    lines = [json.loads(x) for x in text.splitlines() if x.strip()]
    present = {d["region_id"] for d in lines}
    by_id = {d["region_id"]: d for d in lines}
    for rid in VOID:
        new_id = f"{rid}-nbody-retracted-1023"
        if new_id in present:
            print(f"skip (already re-stamped): {new_id}")
            continue
        old = by_id[rid]
        mc = old["method_capability"]
        tags = frozenset(t for t in mc["capability_tags"] if t not in NBODY_TAGS)
        res = old["result"]
        extent = dict(old["search_extent"])
        extent["points_total"] = int(res.get("n_cells") or extent.get("n_sobol_samples") or 0)
        report = EmptyRegionReport(
            region_id=new_id,
            family=old["family"],
            centre=old["centre"],
            topologies=tuple(old["topologies"]),
            method_capability=MethodCapability(
                genome=mc["genome"].replace(" + jovian_shoot(jacobian=stm)", ""),
                corrector="evaluate_joint_cell prefilter only (n-body stage retracted, #1023)",
                capability_tags=tags,
                git_sha=sha,
            ),
            search_extent=extent,
            prune_gates=tuple(g for g in old["prune_gates"] if "jovian_shoot" not in g),
            result={
                "n_cells": res.get("n_cells"),
                "n_feasible": res.get("n_feasible"),
                "n_shot_retracted": res.get("n_shot"),
                "n_close_retracted": res.get("n_close"),
                "original_git_sha": res.get("git_sha"),
            },
            verdict=REASON,
            interpretation=(
                f"Supersedes the n-body part of {rid}. The prefilter part (n_feasible "
                f"{res.get('n_feasible')} of {res.get('n_cells')}) is unchanged; no statement "
                "is made about n-body closure in this region until #1039."
            ),
            source_anchors=(
                f"supersedes {rid} (n-body stage only); "
                "docs/notes/2026-10-08-1023-jovian-void-rerun.md sec. 4; "
                "docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 2.1"
            ),
            run={"date": "2026-10-08", "git_sha": sha, "task": 1023},
        )
        append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
        print(f"appended: {new_id}")


if __name__ == "__main__":
    main()
