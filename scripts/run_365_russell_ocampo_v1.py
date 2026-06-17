"""#365 Phase D — Russell-Ocampo 2003 Tables 5-8 V0→V1 promotion driver.

Runs the §14 V1 like-for-like same-model corrector against the 4 cyclers for
which Russell-Ocampo 2003 (AAS-03-145) publishes per-encounter Δv vectors
in Tables 5-8:

* ``Cycler-2-5-1-3`` → catalogue ``russell-ocampo-2.5.1+0``
* ``Cycler-3-1-2-11`` → catalogue ``russell-ocampo-3.1.2+1``
* ``Cycler-4-3-1-20`` → catalogue ``russell-ocampo-4.3.1-5``
* ``Cycler-4-5-2-12`` → catalogue ``russell-ocampo-4.5.2-2``

§14 V1 like-for-like gate (per ``tests/search/test_closer_sweep_v1.py``):

1. The single circular-coplanar free-return ellipse seeded from the row's
   sourced ``(aphelion_au, transit_em_days)`` closes the corrector residual
   below ``CORRECTOR_TOL_KMS``;
2. The EMERGED per-body V_inf magnitudes at Earth and Mars match the
   independently sourced Russell 2004 Table 3.4 anchor within
   ``TOL_VINF_KMS = 0.5 km/s`` (the spec §14 V1 floor);
3. The §14 V1 mechanics gate passes: closed, V_inf-continuous single
   free-return ellipse.

HONEST NEGATIVES discipline (project memory ``feedback_orbit_closure_discipline``):
nothing here is tuned to force a pass. The single-ellipse-cannot-close
NEGATIVE for the multi-arc rows IS the verdict.

Emits ``data/<row-id>_v1_verdict.jsonl`` per candidate with the per-gate
breakdown and source quotes (verbatim from Russell-Ocampo 2003 Tables 5-8 and
Russell 2004 Table 3.4 — sourced-only EXPECTED values, never our own
computation).

Usage::

    uv run python scripts/run_365_russell_ocampo_v1.py
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

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

# Spec §14 V1 floor — mirrors tests/search/test_closer_sweep_v1.py.
TOL_VINF_KMS = 0.5
CORRECTOR_TOL_KMS = 0.1
PHASE_EPOCHS = 4096

# Russell-Ocampo 2003 Tables 5-8 — the 4 cyclers for which per-encounter Δv
# vectors are published, enabling V1 like-for-like reproduction. Catalogue
# mapping cross-referenced against the digest verdict note
# ``docs/notes/2026-06-17-digest-russell-ocampo-2003.md`` §3 (catalogue's
# signed-i form vs the paper's unsigned-i; verified by V_inf + ToF match).
CANDIDATES: list[dict[str, Any]] = [
    {
        "row_id": "russell-ocampo-2.5.1+0",
        "paper_cycler": "Cycler-2-5-1-3",
        "paper_table": "Table 5 (p.15)",
        "expected_vinf_e": 7.8,  # Russell-Ocampo 2003 Table 4 + Russell 2004 Table 3.4
        "expected_vinf_m": 9.9,
        "expected_transit_em_days": 94,
        "source_quote_vinf": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 4 p.14 row '2-5-1-3': "
            "'AR=1.44 TR=1.12 94 d 7.8 km/s @E 9.9 km/s @M'; "
            "corroborated by Russell 2004 dissertation Table 3.4 row 2.5.1.+0."
        ),
        "source_quote_state": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 5 p.15: per-encounter "
            "Δv vectors (km/s) and time (days) for one complete cycle plus "
            "first leg of second cycle, initial r_mars at t_0 = [1.41, 0.57, 0] AU."
        ),
    },
    {
        "row_id": "russell-ocampo-3.1.2+1",
        "paper_cycler": "Cycler-3-1-2-11",
        "paper_table": "Table 6 (p.15)",
        # Russell 2004 Table 3.4 catalogue row; consistent with paper Table 4
        "expected_vinf_e": 3.4,
        "expected_vinf_m": 4.6,
        "expected_transit_em_days": 181,
        "source_quote_vinf": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 6 p.15 per-encounter "
            "Δv vectors; corroborated by Russell 2004 dissertation Table 3.4 "
            "row 3.1.2.+1 (V_inf at Earth 3.4 km/s, at Mars 4.6 km/s)."
        ),
        "source_quote_state": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 6 p.15: per-encounter "
            "Δv vectors, time (days) 0/181/1083/1265/2348/2529, initial "
            "r_mars at t_0 = [1.15, 0.99, 0] AU."
        ),
    },
    {
        "row_id": "russell-ocampo-4.3.1-5",
        "paper_cycler": "Cycler-4-3-1-20",
        "paper_table": "Table 7 (p.15)",
        "expected_vinf_e": 3.1,  # Russell 2004 Table 3.4 catalogue row; paper Table 4 (3.10/2.53)
        "expected_vinf_m": 2.5,
        "expected_transit_em_days": 268,
        "source_quote_vinf": (
            "Russell-Ocampo 2003 (AAS-03-145) p.15: 'Cycler-4-3-1-20 has "
            "remarkably low energy requirements at Earth and Mars ... "
            "V_inf of 3.10 km/s at Earth ... V_inf of 2.53 km/s at Mars'; "
            "Table 4 p.14: 'AR=0.99 TR=1.29 268 d 3.1 km/s @E 2.5 km/s @M'; "
            "corroborated by Russell 2004 dissertation Table 3.4 row 4.3.1.-5."
        ),
        "source_quote_state": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 7 p.15: per-encounter "
            "Δv vectors, time (days) 0/268/2583/3131/3399, initial "
            "r_mars at t_0 = [0.93, 1.20, 0] AU."
        ),
    },
    {
        "row_id": "russell-ocampo-4.5.2-2",
        "paper_cycler": "Cycler-4-5-2-12",
        "paper_table": "Table 8 (p.15)",
        "expected_vinf_e": 3.4,  # Russell 2004 Table 3.4 catalogue row
        "expected_vinf_m": 4.6,
        "expected_transit_em_days": 191,
        "source_quote_vinf": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 8 p.15 per-encounter "
            "Δv vectors; corroborated by Russell 2004 dissertation Table 3.4 "
            "row 4.5.2.-2."
        ),
        "source_quote_state": (
            "Russell-Ocampo 2003 (AAS-03-145) Table 8 p.15: per-encounter "
            "Δv vectors, time (days) 0/191/1109/1474/1657/2022/3131/3322, "
            "initial r_mars at t_0 = [1.03, 1.12, 0] AU."
        ),
    },
]


def _row(rid: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    for r in rows:
        if r["id"] == rid:
            return r
    raise KeyError(rid)


def _best_phase_t0(
    a_seed: float, e_seed: float, period_sec: float, ephem: Ephemeris
) -> tuple[float, float]:
    """Scan t0 over one period; return the (best_t0, best_residual) at the
    SOURCED ``(a, e)`` (mirrors closer_sweep_v1)."""
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
    return best_t0, best_res


def _git_sha() -> str:
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT)
            .decode()
            .strip()
        )
    except Exception:
        return "UNKNOWN"


def _run_candidate(spec: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    row = _row(spec["row_id"], rows)
    aphelion = float(row["orbit_elements"]["aphelion_au"])
    transit = row["invariants"]["transit_times_days"]
    sourced = {e["body"]: float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"]}
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S

    # Independent sourced anchor (from the digest's verbatim PDF read) vs
    # the catalogue's stored anchor (also from Russell 2004 Table 3.4).
    # These should match - any discrepancy is itself a row defect.
    sourced_e_expected = spec["expected_vinf_e"]
    sourced_m_expected = spec["expected_vinf_m"]
    sourced_anchor_consistent = (
        abs(sourced.get("E", -999.0) - sourced_e_expected) < 0.01
        and abs(sourced.get("M", -999.0) - sourced_m_expected) < 0.01
    )

    ephem = Ephemeris("circular")
    seed_ok = True
    seed_failure: str | None = None
    try:
        a_seed, e_seed = seed_ae_from_aphelion_transit(aphelion, float(transit[0]), mu=MU)
    except Exception as exc:
        seed_ok = False
        seed_failure = f"{type(exc).__name__}: {exc}"
        a_seed = float("nan")
        e_seed = float("nan")

    derived: dict[str, float] = {}
    converged = False
    max_residual_kms = float("inf")
    v1_passed = False
    vinf_continuous = False
    v1_detail = ""
    best_t0 = float("nan")
    best_phase_residual = float("inf")

    if seed_ok:
        best_t0, best_phase_residual = _best_phase_t0(a_seed, e_seed, period_sec, ephem)
        try:
            sol = free_return_correct(
                t0_seed_sec=best_t0,
                a_seed_au=a_seed,
                e_seed=e_seed,
                period_sec=period_sec,
                ephem=ephem,
                mu=MU,
                tol_kms=CORRECTOR_TOL_KMS,
            )
            converged = bool(sol.converged)
            max_residual_kms = float(sol.max_residual_kms)
            if sol.vinf_kms:
                derived = {k: float(v) for k, v in sol.vinf_kms.items()}
            if converged:
                try:
                    v1 = free_return_v1_mechanics(sol, ephem, period_sec, mu=MU)
                    v1_passed = bool(v1.v1_passed)
                    vinf_continuous = bool(v1.vinf_continuous)
                    v1_detail = str(v1.detail)
                except Exception as exc:
                    v1_detail = f"v1 check raised {type(exc).__name__}: {exc}"
        except Exception as exc:
            v1_detail = f"corrector raised {type(exc).__name__}: {exc}"

    # Per-encounter V_inf comparison against sourced anchor.
    d_e = abs(derived["E"] - sourced.get("E", float("nan"))) if "E" in derived else float("nan")
    d_m = abs(derived["M"] - sourced.get("M", float("nan"))) if "M" in derived else float("nan")
    vinf_gate_e = d_e <= TOL_VINF_KMS if "E" in derived else False
    vinf_gate_m = d_m <= TOL_VINF_KMS if "M" in derived else False

    passes_v1 = bool(converged and vinf_gate_e and vinf_gate_m and v1_passed)
    verdict = "PASS" if passes_v1 else "FAIL"

    if not seed_ok:
        fail_mode = "NOT-REACHABLE (seed failed — aphelion below Mars sma or geometry)"
    elif not converged:
        fail_mode = "NO-CLOSE (single ellipse cannot close at the sourced geometry)"
    elif not (vinf_gate_e and vinf_gate_m):
        fail_mode = f"VINF-MISS (|ΔE|={d_e:.3f}, |ΔM|={d_m:.3f}, floor 0.5 km/s)"
    elif not v1_passed:
        fail_mode = "CLOSE-NOT-V1 (vinf-continuity break — multi-arc, single ellipse insufficient)"
    else:
        fail_mode = "—"

    return {
        "kind": "v1_verdict_russell_ocampo_2003_tables_5_8",
        "task": "#365 Phase D — Russell-Ocampo 2003 Tables 5-8 V0→V1 promotion",
        "candidate_id": spec["row_id"],
        "paper_cycler": spec["paper_cycler"],
        "paper_table": spec["paper_table"],
        "passes_v1": passes_v1,
        "verdict": verdict,
        "fail_mode": fail_mode,
        "v1_floor_kms": TOL_VINF_KMS,
        "corrector_tol_kms": CORRECTOR_TOL_KMS,
        "sourced_vinf_kms": {"E": sourced.get("E"), "M": sourced.get("M")},
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
        "sourced_anchor_consistent_with_paper": sourced_anchor_consistent,
        "expected_vinf_kms_from_paper": {
            "E": sourced_e_expected,
            "M": sourced_m_expected,
        },
        "seed": {
            "a_seed_au": a_seed,
            "e_seed": e_seed,
            "best_t0_sec": best_t0,
            "best_phase_residual_kms": best_phase_residual,
            "seed_failure": seed_failure,
        },
        "corrector": {
            "converged": converged,
            "max_residual_kms": max_residual_kms,
        },
        "source_quotes": {
            "vinf": spec["source_quote_vinf"],
            "state": spec["source_quote_state"],
        },
        "git_sha": _git_sha(),
        "generated_at_utc": datetime.now(UTC).isoformat(),
    }


def main() -> None:
    catalogue_path = REPO_ROOT / "data" / "catalogue.yaml"
    rows = yaml.safe_load(catalogue_path.read_text())

    for spec in CANDIDATES:
        verdict = _run_candidate(spec, rows)
        out_path = REPO_ROOT / "data" / f"{spec['row_id']}_v1_verdict.jsonl"
        with out_path.open("w") as fh:
            fh.write(json.dumps(verdict) + "\n")
        print(
            f"{spec['row_id']}: {verdict['verdict']} ({verdict['fail_mode']}) — "
            f"derived E={verdict['derived_vinf_kms'].get('E')}, "
            f"M={verdict['derived_vinf_kms'].get('M')}; "
            f"sourced E={verdict['sourced_vinf_kms']['E']}, "
            f"M={verdict['sourced_vinf_kms']['M']}"
        )


if __name__ == "__main__":
    main()
