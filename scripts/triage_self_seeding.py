#!/usr/bin/env python3
"""Self-seeding REACHABLE-vs-OFF-FAMILY triage of the unsourced catalogue rows (#177).

The cheap per-row pass of task #177 build 2: for every UNSOURCED catalogue row (the
russell-ch4 rows without a confirmed App-C / V3 closure, plus the ~200 russell-ocampo
members) decide REACHABLE vs OFF-FAMILY by the #173 gating condition — does any
coplanar G-arc *branch* transit (the #177 multi-rev Stage-A enumeration) land within
``tol_days`` of the row's tabulated (Russell simple-model) E->M transit signature?

This is the bulk, kept FAST (no n-body): a row is OFF-FAMILY here when NO branch of
its coplanar descriptor shape reproduces its transit signature — exactly the 6.44Gg3
mechanism #173 found (coplanar short-way 131 d vs the row's 262 d), now tested against
ALL branches (short / long / k-rev) before declaring OFF-FAMILY.

Honesty (binding, brief): the row's tabulated transit is the EXPECTED side; the
emerged branch transits are EVIDENCE. ``tol_days`` is the GATE and is held FIXED
across all rows (NOT loosened per row to inflate REACHABLE). A large OFF-FAMILY count
is a valid, important result. The OFF-FAMILY negatives written here are method-versioned
("self-seed multirev v1, coplanar-branch transit gate") and #172's empty-region
registry should ingest this runlog.

ocampo rows carry NO ``free_return_arcs[]`` g/G ToF descriptor (only aphelion + the
simple-model transit + v_inf anchors), and the aphelion-ratio<1 members are sub-Mars
(reach Mars only near Mars perihelion) — the single coplanar G-arc shape cannot be
DERIVED for them without the free-return arc ToFs. They are recorded honestly as
``OFF-FAMILY-NO-DESCRIPTOR`` (the coplanar-branch gate is inapplicable), NOT silently
dropped and NOT counted as REACHABLE. Reaching them needs a descriptor-recovery step
out of this task's scope; the negative is recorded for the registry.

Resumable: re-running skips rows already present in the target runlog.
"""

from __future__ import annotations

import argparse
import subprocess
import warnings
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.runlog import RunLog, RunRecord, default_runlog_path
from cyclerfinder.search.self_seeding import g_arc_branches, triage_transit_match

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
DEFAULT_RUNLOG_DIR = REPO_ROOT / "data" / "runs"

# The transit-match gate (days). FIXED across all rows. Calibrated, NOT tuned per row:
# the validated S1L1 PASS has coplanar short-way 169 d vs tabulated 150 d (delta +19);
# 6.44Gg3's long branch 292 vs 262 (delta +30). Held at 30 d so S1L1 (the known PASS)
# is admitted; NOT widened to manufacture REACHABLE for any other row (brief honesty).
TOL_DAYS = 30.0
MAX_G_REVS = 2
METHOD = "self-seed-multirev-v1/coplanar-branch-transit-gate"


def _code_version() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def _unsourced_rows() -> list[dict[str, Any]]:
    """Unsourced ch4 (non-V3) + all ocampo rows, sorted by id (resumable order)."""
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    out = [
        r
        for r in rows
        if (str(r["id"]).startswith("russell-ch4") and r.get("validation_level") != "V3")
        or str(r["id"]).startswith("russell-ocampo")
    ]
    out.sort(key=lambda r: str(r["id"]))
    return out


def _triage_row(row: dict[str, Any], code_version: str) -> RunRecord:
    """Cheap coplanar-branch transit-match triage of one row -> a RunRecord."""
    rid = str(row["id"])
    aph = (row.get("orbit_elements") or {}).get("aphelion_au")
    vinf = {e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])}
    transit = (row.get("invariants") or {}).get("transit_times_days")
    fra = row.get("free_return_arcs") or []
    g_tofs = [a.get("tof_years") for a in fra if a.get("tof_years") is not None]

    base = dict(
        row_id=rid,
        genome="self-seed-multirev",
        model="circular-coplanar",
        code_version=code_version,
        sourced_vinf_kms=vinf,
        sourced_anchors={
            "aphelion_au": aph,
            "transit_times_days": transit,
            "turn_ratio": (row.get("invariants") or {}).get("turn_ratio"),
            "method": METHOD,
            "tol_days": TOL_DAYS,
        },
    )

    # Rows without the g/G free-return arc descriptor cannot have a coplanar G-arc
    # shape DERIVED (ocampo members; the few ch4 rows with <2 arc ToFs). Honest
    # inapplicable-gate negative — recorded, never counted REACHABLE.
    if len(g_tofs) < 2 or aph is None or "E" not in vinf or "M" not in vinf or not transit:
        return RunRecord(
            outcome="OFF-FAMILY-NO-DESCRIPTOR",
            seed={"g_tofs_years": g_tofs, "reason": "no 2-arc g/G free-return descriptor"},
            solver_audit={"gate": "transit-match", "applicable": False},
            **base,
        )

    real_transit = float(transit[0])
    try:
        branches = g_arc_branches(
            aph, g_tofs[0], g_tofs[1], vinf["E"], vinf["M"], max_g_revs=MAX_G_REVS
        )
    except (ValueError, RuntimeError) as exc:
        return RunRecord(
            outcome="OFF-FAMILY-NO-CLOSE",
            seed={"g_tofs_years": g_tofs},
            solver_audit={"gate": "transit-match", "applicable": True, "error": str(exc)},
            **base,
        )

    tr = triage_transit_match(branches, real_transit, tol_days=TOL_DAYS)
    outcome = "REACHABLE" if tr.reachable else "OFF-FAMILY"
    return RunRecord(
        outcome=outcome,
        residual_kms=None,
        seed={
            "g_tofs_years": g_tofs,
            "g_arc_a_au": branches[0].a_au,
            "g_arc_e": branches[0].e,
        },
        solver_audit={
            "gate": "transit-match",
            "applicable": True,
            "real_eph_transit_days": tr.real_eph_transit_days,
            "best_branch": tr.best_branch,
            "best_g_revs": tr.best_g_revs,
            "best_tof_days": tr.best_tof_days,
            "delta_days": tr.delta_days,
            "branch_tofs_days": tr.branch_tofs,
        },
        **base,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runlog-dir", default=str(DEFAULT_RUNLOG_DIR))
    ap.add_argument("--timestamp", default=None, help="filesystem-safe stamp for the runlog name")
    ap.add_argument("--limit", type=int, default=None, help="cap rows processed (resumable)")
    args = ap.parse_args()

    warnings.filterwarnings("ignore")
    code_version = _code_version()
    rows = _unsourced_rows()

    path = default_runlog_path(args.runlog_dir, "self-seeding-triage", args.timestamp)
    runlog = RunLog(path)
    done = {rec.row_id for rec in runlog.read()}
    pending = [r for r in rows if str(r["id"]) not in done]
    if args.limit is not None:
        pending = pending[: args.limit]

    counts: dict[str, int] = {}
    for row in pending:
        rec = _triage_row(row, code_version)
        runlog.append(rec)
        counts[rec.outcome] = counts.get(rec.outcome, 0) + 1
        if rec.outcome in ("REACHABLE",):
            print(
                f"  {rec.row_id:36s} {rec.outcome}  {rec.solver_audit.get('best_branch')}@"
                f"{rec.solver_audit.get('best_tof_days'):.0f}d "
                f"d={rec.solver_audit.get('delta_days'):+.0f}"
            )

    total_done = len(done) + len(pending)
    print(f"\n[{datetime.now(UTC).isoformat()}] triage runlog: {path}")
    print(f"processed {len(pending)} new rows (now {total_done}/{len(rows)} total)")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
