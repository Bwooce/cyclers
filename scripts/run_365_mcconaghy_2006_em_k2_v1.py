"""#365 Phase D — McConaghy 2006 em-k2 (S1L1) V0→V1 attempt + honest-negative.

The McConaghy 2006 JSR digest (``docs/notes/2026-06-17-digest-mcconaghy-2006.md``)
identified ``mcconaghy-2006-em-k2`` as the 8-vehicle V0→V1 candidate sink for
Tables 2-9 of the McConaghy/Landau/Yam 2006 paper (the four cycler-Δv-only
itineraries + four cycler+taxi itineraries). The McConaghy 2004 JSR S1L1 22-
encounter DE405 itinerary (Table 6) also lands at this row (same physical
cycler, different launch epoch).

The same-model §14 V1 like-for-like ask: does the closer_sweep_v1 substrate
(single circular-coplanar free-return ellipse) reproduce the row's abstract
V_inf anchor (E=4.7, M=5.0 km/s, transit 153 d) within 0.5 km/s AND clear the
V_inf-continuity gate?

VERDICT (2026-06-17): the corrector closes with EMERGED V_inf E=4.771,
M=5.036 km/s (|ΔE|=0.07, |ΔM|=0.04, BOTH < 0.5 km/s V1 floor), BUT the
§14 V1 mechanics gate REJECTS the row: V_inf-continuity break of 24 km/s at
the Mars flyby, the reconstructed return leg does NOT close to Earth — the
S1L1 cycler is MULTI-ARC (two generic Earth-Earth arcs joined at the Mars
flyby, McConaghy-Russell-Longuski 2005 label
``2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)``) and a single radial-
crossing ellipse fundamentally cannot represent it.

HONEST NEGATIVE per ``feedback_orbit_closure_discipline``: V0 stays V0; the
multi-arc topology is the verdict, not a corrector deficiency. The 8
McConaghy 2006 Tables 2-9 vehicles + the McConaghy 2004 Table 6 itinerary
all condense into this single honest-negative gate at the V1 tier. The
em-k2 row's REAL elevated-tier evidence lives at the sibling
``russell-ch4-4.991gG2`` row (V3, S1L1 corrected closure on real DE440
ephemeris with REBOUND/IAS15 cross-check, see #167/#94) — em-k2's V0 is by
PUBLICATION GAP at THIS row's metric (a circular-coplanar single-ellipse
anchor), not by validation infrastructure.

Usage::

    uv run python scripts/run_365_mcconaghy_2006_em_k2_v1.py
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.free_return import _residuals as _fr_residuals
from cyclerfinder.search.free_return import (
    free_return_correct,
    seed_ae_from_aphelion_transit,
)
from cyclerfinder.search.free_return_v1 import free_return_v1_mechanics

REPO_ROOT = Path(__file__).resolve().parent.parent
DAY_S = 86400.0
MU = MU_SUN_KM3_S2

TOL_VINF_KMS = 0.5
CORRECTOR_TOL_KMS = 0.1
PHASE_EPOCHS = 4096

ROW_ID = "mcconaghy-2006-em-k2"


def _git_sha() -> str:
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT)
            .decode()
            .strip()
        )
    except Exception:
        return "UNKNOWN"


def main() -> None:
    catalogue_path = REPO_ROOT / "data" / "catalogue.yaml"
    rows = yaml.safe_load(catalogue_path.read_text())
    row = next(r for r in rows if r["id"] == ROW_ID)
    aphelion = float(row["orbit_elements"]["aphelion_au"])
    transit = row["invariants"]["transit_times_days"]
    sourced = {e["body"]: float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"]}
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S

    ephem = Ephemeris("circular")
    a_seed, e_seed = seed_ae_from_aphelion_transit(aphelion, float(transit[0]), mu=MU)

    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, PHASE_EPOCHS, endpoint=False):
        t0 = float(frac) * period_sec
        try:
            res = _fr_residuals(
                np.array([a_seed, e_seed, t0]),
                period_days=period_sec / DAY_S,
                ephem=ephem,
                bodies=("E", "M"),
                mu=MU,
            )
        except Exception:
            continue
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0

    sol = free_return_correct(
        t0_seed_sec=best_t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=ephem,
        mu=MU,
        tol_kms=CORRECTOR_TOL_KMS,
    )
    derived = {k: float(v) for k, v in (sol.vinf_kms or {}).items()}
    converged = bool(sol.converged)
    d_e = abs(derived.get("E", float("nan")) - sourced["E"]) if "E" in derived else float("nan")
    d_m = abs(derived.get("M", float("nan")) - sourced["M"]) if "M" in derived else float("nan")
    vinf_gate_e = d_e <= TOL_VINF_KMS if "E" in derived else False
    vinf_gate_m = d_m <= TOL_VINF_KMS if "M" in derived else False

    v1_passed = False
    vinf_continuous = False
    v1_detail = ""
    if converged:
        try:
            v1 = free_return_v1_mechanics(sol, ephem, period_sec, mu=MU)
            v1_passed = bool(v1.v1_passed)
            vinf_continuous = bool(v1.vinf_continuous)
            v1_detail = str(v1.detail)
        except Exception as exc:
            v1_detail = f"v1 check raised {type(exc).__name__}: {exc}"

    passes_v1 = bool(converged and vinf_gate_e and vinf_gate_m and v1_passed)
    if passes_v1:
        verdict = "PASS"
        fail_mode = "—"
    elif converged and vinf_gate_e and vinf_gate_m:
        verdict = "FAIL"
        fail_mode = "CLOSE-NOT-V1 (vinf-continuity break — multi-arc, single ellipse insufficient)"
    elif not converged:
        verdict = "FAIL"
        fail_mode = "NO-CLOSE (single ellipse cannot close at the sourced geometry)"
    else:
        verdict = "FAIL"
        fail_mode = f"VINF-MISS (|ΔE|={d_e:.3f}, |ΔM|={d_m:.3f}, floor 0.5 km/s)"

    out = {
        "kind": "v1_verdict_mcconaghy_2006_em_k2",
        "task": "#365 Phase D — McConaghy 2006 em-k2 (S1L1) V0→V1 attempt",
        "candidate_id": ROW_ID,
        "paper_cycler": "S1L1 (ballistic two-synodic Earth-Mars)",
        "paper_table": (
            "McConaghy 2006 JSR 43(2) Tables 2-9 (8 DE405 itineraries); "
            "McConaghy 2004 JSR 41(4) Table 6 (22-encounter DE405 itinerary)"
        ),
        "passes_v1": passes_v1,
        "verdict": verdict,
        "fail_mode": fail_mode,
        "v1_floor_kms": TOL_VINF_KMS,
        "corrector_tol_kms": CORRECTOR_TOL_KMS,
        "sourced_vinf_kms": {"E": sourced["E"], "M": sourced["M"]},
        "derived_vinf_kms": derived,
        "delta_vinf_kms": {
            "E": d_e if "E" in derived else None,
            "M": d_m if "M" in derived else None,
        },
        "vinf_gate_E_passed": vinf_gate_e,
        "vinf_gate_M_passed": vinf_gate_m,
        "v1_mechanics_passed": v1_passed,
        "vinf_continuous": vinf_continuous,
        "v1_detail": v1_detail,
        "seed": {
            "a_seed_au": a_seed,
            "e_seed": e_seed,
            "best_t0_sec": best_t0,
            "best_phase_residual_kms": best_res,
        },
        "corrector": {
            "converged": converged,
            "max_residual_kms": float(sol.max_residual_kms),
        },
        "source_quotes": {
            "vinf": (
                "McConaghy 2006 JSR 43(2) p.458 Table 1 (circular-coplanar S1L1): "
                "'Δv required 0 km/s; Earth V_inf 4.7 km/s; Earth flyby altitude "
                "31,809 km; Mars V_inf 5.0 km/s; Earth-Mars transfer time 153.15 days; "
                "Repeat time 4-2/7 yr.'"
            ),
            "topology": (
                "McConaghy-Russell-Longuski 2005 JSR 42(4) Table 2: 'Ballistic S1L1 "
                "cycler = 2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)' — TWO "
                "generic Earth-Earth arcs, total t_f = 4.2785 yr; multi-arc by "
                "construction. McConaghy 2006 p.457: 'the first leg is the short-"
                "period solution making 1-2 revs, and the second leg is the "
                "long-period solution making 1-2 revs.'"
            ),
            "honest_negative_disposition": (
                "feedback_orbit_closure_discipline: the math gives the verdict, "
                "not a tuned pass. em-k2's REAL elevated-tier evidence lives at the "
                "sibling row russell-ch4-4.991gG2 (V3, real-eph DE440 closure with "
                "REBOUND/IAS15 cross-check per #167/#94). em-k2 STAYS V0 at the "
                "single-ellipse V1 gate by topology — multi-arc construction is the "
                "anchor's identity, not a corrector deficiency."
            ),
        },
        "git_sha": _git_sha(),
        "generated_at_utc": datetime.now(UTC).isoformat(),
    }

    out_path = REPO_ROOT / "data" / f"{ROW_ID}_v1_verdict.jsonl"
    with out_path.open("w") as fh:
        fh.write(json.dumps(out) + "\n")

    print(
        f"{ROW_ID}: {verdict} ({fail_mode}) — "
        f"derived E={derived.get('E')}, M={derived.get('M')}; "
        f"sourced E={sourced['E']}, M={sourced['M']}"
    )


if __name__ == "__main__":
    main()
