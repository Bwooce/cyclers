"""#829 — re-classify EVERY historical ``russell12-*`` runlog under the bend gate.

Why this exists (`[[feedback_bugfix_invalidates_past_searches]]`)
----------------------------------------------------------------
``campaign_russell12.py`` filtered its epoch-grid results on
``BallisticClosureResult.converged`` ALONE until #829. Spec §14 **V0** requires
the hard constraints (V_inf cap, r_p >= r_p_min, **bend <= max**) on top of the
closure residual, so a ``converged``-only pick is not a V0-admissible closure —
and because magnitude mode drives the powered/degenerate basin to ~1e-14 while a
bend-feasible solution sits at a comparable residual, ranking by residual over
the converged-only set SYSTEMATICALLY returned a bend-infeasible trajectory.
Every CLOSE-AND-MATCH / CLOSE-OFF-ANCHOR verdict this campaign has ever emitted
was therefore bend-blind, not just #820's — so the historical runlogs need
re-classifying, not only the forward code path.

What this script can and cannot conclude
----------------------------------------
A runlog record stores only the BEST result the campaign picked, not the whole
epoch grid. So for each historical record this script can say:

* the recorded verdict is **RETRACTED** — the stored best violates a V0 hard
  constraint, so that CLOSE-* label was never admissible; or
* the recorded verdict **STANDS** — the stored best satisfies them.

It can NOT conclude that no admissible closure exists at some other epoch of that
run: that information was never persisted. Answering THAT requires re-running,
which #829 did for the current (#820-reposed) genome —
``russell12-circular-20260812T-829-bendgate.jsonl``. The pre-#820 runs cannot be
meaningfully re-run at all: their genome was mis-posed (#820), so a re-run would
answer a different question. Record-level retraction is the honest instrument
there.

A retraction is NOT a validation-level downgrade. It withdraws a campaign label;
any ``validation_level`` change in either direction is an adjudication task.

Usage::

    uv run python scripts/reclassify_russell12_runlogs.py [--out FILE]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = REPO_ROOT / "data" / "runs"
GLOB = "russell12-*.jsonl"

# Verdicts that assert a closure was found (the ones a bend-blind filter could
# wrongly award). NO-CLOSE asserts the opposite and cannot be over-claimed.
CLOSURE_OUTCOMES = {"CLOSE-AND-MATCH", "CLOSE-OFF-ANCHOR", "CLOSE-MATCH-SYMMETRIC-ONLY"}


def reclassify_record(rec: dict[str, Any]) -> dict[str, Any]:
    """Apply the #829 hard-constraint gate to one persisted runlog record."""
    audit = rec.get("solver_audit") or {}
    outcome = str(rec.get("outcome"))
    bend = audit.get("bend_feasible")
    cap = audit.get("vinf_cap_ok")
    violations = [
        name
        for name, ok in (("bend", bend), ("vinf_cap", cap))
        if ok is False  # None = unrecorded
    ]
    if outcome not in CLOSURE_OUTCOMES:
        status, revised = "N/A", outcome
    elif bend is None and cap is None:
        status, revised = "UNKNOWN", outcome
    elif violations:
        status, revised = "RETRACTED", "CLOSE-INADMISSIBLE"
    else:
        status, revised = "STANDS", outcome
    return {
        "row_id": rec.get("row_id"),
        "genome": rec.get("genome"),
        "model": rec.get("model"),
        "code_version": rec.get("code_version"),
        "recorded_outcome": outcome,
        "bend_feasible": bend,
        "vinf_cap_ok": cap,
        "violations": violations,
        "status": status,
        "revised_outcome": revised,
        "measured_turn_ratio": audit.get("measured_turn_ratio"),
        "sourced_turn_ratio": (rec.get("sourced_anchors") or {}).get("sourced_turn_ratio"),
    }


def reclassify_runlogs(runs_dir: Path = RUNS_DIR) -> dict[str, Any]:
    """Re-classify every ``russell12-*.jsonl`` runlog in ``runs_dir``."""
    files: dict[str, list[dict[str, Any]]] = {}
    for path in sorted(runs_dir.glob(GLOB)):
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        files[path.name] = [reclassify_record(r) for r in rows]
    counts: dict[str, int] = {}
    for entries in files.values():
        for e in entries:
            counts[e["status"]] = counts.get(e["status"], 0) + 1
    return {"files": files, "counts": counts}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", type=str, default=str(RUNS_DIR))
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    result = reclassify_runlogs(Path(args.runs_dir))
    for name, entries in result["files"].items():
        print(f"=== {name}")
        for e in entries:
            note = f" violated={','.join(e['violations'])}" if e["violations"] else ""
            print(
                f"  {e['row_id']!s:24s} {e['recorded_outcome']:26s} -> "
                f"{e['status']:9s} {e['revised_outcome']}{note}"
            )
    print("\ncounts:", result["counts"])
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
