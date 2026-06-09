#!/usr/bin/env python3
"""Validate the REACHABLE subset of the #177 self-seeding triage (build 3).

For every row the cheap transit-match triage flagged REACHABLE, run the FULL tail:
self-seed (multi-rev best branch) -> corrected-topology longitude rendezvous ->
INDEPENDENT n-body confirm (REBOUND/IAS15, Sun-only — Russell's patched-conic cruise
model, the #167 arbiter) -> a V3-CANDIDATE recommendation (the family-correct window
+ the n-body in-band verdict) or an honest negative.

The decisive scientific term is the EMERGED v_inf vs the row's SOURCED v_inf anchor
(the term that failed for 6.44Gg3 in #173). The transit gate only proves a branch
reaches Mars at the right transit; on-family requires the emerged v_inf to also match
the anchors. PER-ROW, NO BATCH-TRUST: each row is confirmed (or not) on its own; one
row's CONFIRM never transfers to another (brief rule).

NO catalogue writeback — V3-CANDIDATE is a RECOMMENDATION held for main-session review.
Bands / tolerances are the #165/#167 constants, NEVER loosened to inflate a PASS.

Honesty: the row's SOURCED v_inf anchors + tabulated transit are EXPECTED; the emerged
epoch / v_inf / miss are EVIDENCE. A clean OFF-FAMILY-AT-ANCHOR-VINF is a first-class
result, recorded for the registry.

Resumable: re-running skips rows already present in the target runlog.
"""

from __future__ import annotations

import argparse
import subprocess
import warnings
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.runlog import RunLog, RunRecord, default_runlog_path
from cyclerfinder.nbody.propagator import RestrictedNBody
from cyclerfinder.search.self_seeding import (
    FamilyAnchors,
    g_arc_branches,
    on_family,
    self_seed_g_leg,
    triage_transit_match,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
DEFAULT_RUNLOG_DIR = REPO_ROOT / "data" / "runs"
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

TOL_DAYS = 30.0
MAX_G_REVS = 2
VINF_BAND_KMS = 1.5  # the #173 breathing band; NOT loosened
METHOD = "self-seed-multirev-v1/full-tail"

# The 3-Mars-SOI encounter band — the #165/#167 constant, NEVER loosened.
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
BAND_AU = 3.0 * _MARS_SOI_AU
# Real-eph v_inf breathing ceiling at Mars (S1L1 App-C 3.2-8.0); EVIDENCE band only.
REAL_EPH_VINF_CEILING = 8.2


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


def _reachable_ids(triage_runlog: Path) -> list[str]:
    """Row ids the triage runlog flagged REACHABLE (sorted)."""
    rl = RunLog(triage_runlog)
    ids = sorted({rec.row_id for rec in rl.read() if rec.outcome == "REACHABLE"})
    return ids


def _row_by_id(rows: list[dict[str, Any]], rid: str) -> dict[str, Any]:
    return next(r for r in rows if str(r["id"]) == rid)


def _validate_row(row: dict[str, Any], ephem: Ephemeris, code_version: str) -> RunRecord:
    """Full self-seed + independent n-body confirm of one REACHABLE row."""
    rid = str(row["id"])
    aph = float((row.get("orbit_elements") or {}).get("aphelion_au"))
    vinf = {e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])}
    transit = (row.get("invariants") or {}).get("transit_times_days")
    fra = row.get("free_return_arcs") or []
    g_tofs = [a.get("tof_years") for a in fra if a.get("tof_years") is not None]

    branches = g_arc_branches(
        aph, g_tofs[0], g_tofs[1], vinf["E"], vinf["M"], max_g_revs=MAX_G_REVS
    )
    tr = triage_transit_match(branches, float(transit[0]), tol_days=TOL_DAYS)
    best = next(b for b in branches if b.branch == tr.best_branch)

    # Generic 2027 launch window (no published seed centres it).
    t_center = (datetime(2027, 1, 1, tzinfo=UTC) - J2000).total_seconds()
    results = self_seed_g_leg(best, ephem, t_center, refine=True)
    anchors = FamilyAnchors(vinf_e=vinf["E"], vinf_m=vinf["M"], vinf_band_kms=VINF_BAND_KMS)

    prop = RestrictedNBody("rebound")
    best_rec: dict[str, Any] | None = None
    for r in results:
        verdict = on_family(r, anchors)
        # Independent integrator confirms the geometric arrival (the arbiter).
        r_e, v_e = ephem.state("E", r.t_depart_sec)
        r0 = np.asarray(r_e, dtype=np.float64)
        v0 = np.asarray(v_e, dtype=np.float64) + r.vinf_vec
        out = prop.propagate(r0, v0, r.t_depart_sec, r.t_arrive_sec, accuracy=1e-11)
        r_m, v_m = ephem.state("M", r.t_arrive_sec)
        nbody_miss = float(np.linalg.norm(out.r_km - np.asarray(r_m)) / AU_KM)
        nbody_vinf_m = float(np.linalg.norm(out.v_km_s - np.asarray(v_m)))
        cand = {
            "branch": best.branch,
            "tof_g_days": best.tof_g_days,
            "depart_offset_days": (r.t_depart_sec - t_center) / 86400.0,
            "vinf_e_kms": r.vinf_e_kms,
            "vinf_m_kms": r.vinf_m_kms,
            "residual_lon_deg": r.residual_lon_deg,
            "lambert_miss_au": r.mars_miss_au,
            "nbody_converged": bool(out.converged),
            "nbody_energy_drift": float(out.energy_rel_drift),
            "nbody_miss_au": nbody_miss,
            "nbody_vinf_m_kms": nbody_vinf_m,
            "on_family": bool(verdict.on_family),
            "vinf_e_ok": bool(verdict.vinf_e_ok),
            "vinf_m_ok": bool(verdict.vinf_m_ok),
            "miss_ok": bool(verdict.miss_ok),
            "in_real_eph_vinf_band": bool(nbody_vinf_m < REAL_EPH_VINF_CEILING),
        }
        # Keep the candidate closest to the v_inf_M anchor (the decisive term).
        if best_rec is None or abs(cand["vinf_m_kms"] - vinf["M"]) < abs(
            best_rec["vinf_m_kms"] - vinf["M"]
        ):
            best_rec = cand

    assert best_rec is not None
    # Verdict: V3-CANDIDATE only if on-family AND the independent n-body confirms the
    # encounter in-band at a real-eph v_inf. Otherwise the honest three-way negative.
    nbody_in_band = best_rec["nbody_converged"] and best_rec["nbody_miss_au"] < BAND_AU
    if best_rec["on_family"] and nbody_in_band:
        outcome = "V3-CANDIDATE"
    elif nbody_in_band:
        # Longitude rendezvous + in-band geometric arrival, but v_inf off the anchor.
        outcome = "OFF-FAMILY-AT-ANCHOR-VINF"
    else:
        outcome = "OFF-FAMILY-NO-CONFIRM"

    achieved = {"E": best_rec["vinf_e_kms"], "M": best_rec["vinf_m_kms"]}
    return RunRecord(
        row_id=rid,
        genome="self-seed-multirev",
        outcome=outcome,
        model="astropy",
        code_version=code_version,
        achieved_vinf_kms=achieved,
        sourced_vinf_kms=vinf,
        sourced_anchors={
            "aphelion_au": aph,
            "transit_times_days": transit,
            "method": METHOD,
            "vinf_band_kms": VINF_BAND_KMS,
            "band_au": BAND_AU,
        },
        seed={"g_tofs_years": g_tofs, "best_branch": best.branch, "tof_g_days": best.tof_g_days},
        residual_kms=None,
        solver_audit=best_rec,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--triage-runlog", required=True, help="the REACHABLE source runlog")
    ap.add_argument("--runlog-dir", default=str(DEFAULT_RUNLOG_DIR))
    ap.add_argument("--timestamp", default=None)
    args = ap.parse_args()

    warnings.filterwarnings("ignore")
    code_version = _code_version()
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    reachable = _reachable_ids(Path(args.triage_runlog))

    ephem = Ephemeris("astropy")
    path = default_runlog_path(args.runlog_dir, "self-seeding-reachable", args.timestamp)
    runlog = RunLog(path)
    done = {rec.row_id for rec in runlog.read()}

    counts: dict[str, int] = {}
    for rid in reachable:
        if rid in done:
            continue
        print(f"[{datetime.now(UTC).isoformat()}] validating {rid} ...")
        rec = _validate_row(_row_by_id(rows, rid), ephem, code_version)
        runlog.append(rec)
        counts[rec.outcome] = counts.get(rec.outcome, 0) + 1
        sa = rec.solver_audit
        print(
            f"   {rec.outcome}: branch={sa['branch']} vinfE={sa['vinf_e_kms']:.2f}"
            f"(anc {rec.sourced_vinf_kms['E']}) vinfM={sa['vinf_m_kms']:.2f}"
            f"(anc {rec.sourced_vinf_kms['M']}) nbody_miss={sa['nbody_miss_au']:.2e}AU"
            f" nbody_vinfM={sa['nbody_vinf_m_kms']:.2f}"
        )

    print(f"\nreachable-validation runlog: {path}")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
