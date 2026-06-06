"""#125 Part 3 — gauntlet tier sweep over all 237 catalogue rows -> the ledger.

Runs :func:`cyclerfinder.verify.gauntlet.run_gauntlet` over every catalogue row
TO THE EXTENT each row's data supports, and writes one
:class:`cyclerfinder.data.ledger.LedgerEntry` per row (carrying ``verdict_tier``
and the full ``verdict_audit``) to ``data/gauntlet_ledger.jsonl``.

Axes supplied per row
---------------------
* **Axis C (provenance)** — :func:`classify_validation` over the row's
  back-filled provenance tags (``orbit_source`` / ``vinf_source`` /
  ``*_fidelity``), plus a corroboration label. Every row has these.
* **Axis A (agreement)** — supplied ONLY for the rows with RECORDED real-closure
  code-path agreement evidence (the Aldrin pair: the §14 V1 lamberthub + Kepler
  + coplanar paths exercised by ``tests/verify/test_agreement_lamberthub.py``).
  For all other rows Axis A is unavailable — no test builds/cross-checks them on
  real ephemeris, so the gauntlet must not pretend they were machine-confirmed.
* **Axis B / D** — not run / not falsified for any row (no adversarial pass).

Consequence of the gauntlet rules (``verify/gauntlet.py``): GOLD/SILVER need
machine-confirmation (Axis A available AND agreed). Only the Aldrin INBOUND row
has that (agreed=True) -> cross_validated -> GOLD. The Aldrin OUTBOUND row has
Axis A available-but-FAILING (path b, the coplanar-vs-real construction, vetoes:
agreed=False) -> REJECTED — the honest multi-axis verdict, stricter than the
single-axis §14 V1 gate the row still legitimately holds. Every other row has
Axis A unavailable and no failing axis -> BRONZE.

NO catalogue field changes from this part — verdicts live in the ledger only.

The Axis-A agreement classifications for the Aldrin pair are pinned constants
here, traceable to the recorded test evidence (not recomputed at sweep time, so
the sweep stays fast and offline). They mirror exactly what
``crosscheck_code_paths`` returns on the real-DE440 build (verified 2026-06-06).

Usage::

    uv run python scripts/sweep_gauntlet_ledger.py
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.ledger import Ledger, LedgerEntry
from cyclerfinder.data.provenance import Corroboration, Tier, classify_validation
from cyclerfinder.verify.agreement import (
    AgreementReport,
    ConstructionOptimiserPathResult,
    KeplerRepropPathResult,
    LamberthubPathResult,
)
from cyclerfinder.verify.gauntlet import run_gauntlet, validate_verdict

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
LEDGER_PATH = REPO_ROOT / "data" / "gauntlet_ledger.jsonl"


def _empty_paths() -> tuple:
    a = LamberthubPathResult(available=False, per_leg=(), max_diff_mps=0.0, passed=False)
    b = ConstructionOptimiserPathResult(
        available=False,
        resonant_vinf_kms={},
        cycler_vinf_kms={},
        construction_max_diff_kms=float("inf"),
        optimiser_available=False,
        optimiser_vinf_kms={},
        optimiser_max_diff_kms=None,
        max_diff_kms=float("inf"),
        passed=False,
    )
    c = KeplerRepropPathResult(
        available=False, per_leg_residual_km=(), max_residual_km=float("inf"), passed=False
    )
    return a, b, c


def _agreement(*, agreed: bool, n_available: int, n_passed: int) -> AgreementReport:
    a, b, c = _empty_paths()
    return AgreementReport(
        lamberthub=a,
        construction_optimiser=b,
        kepler_reprop=c,
        n_paths_available=n_available,
        n_paths_passed=n_passed,
        agreed=agreed,
    )


# Recorded Axis-A agreement evidence for the Aldrin pair (real-DE440), traceable
# to tests/verify/test_agreement_lamberthub.py (verified 2026-06-06):
#   outbound — 3 paths available, 2 passed, agreed False (path b vetoes).
#   inbound  — 2 paths available, 2 passed, agreed True.
_AXIS_A_AGREEMENT: dict[str, AgreementReport] = {
    "aldrin-classic-em-k1-outbound": _agreement(agreed=False, n_available=3, n_passed=2),
    "aldrin-classic-em-k1-inbound": _agreement(agreed=True, n_available=2, n_passed=2),
}


def _row_tier(row: dict) -> Tier:  # type: ignore[type-arg]
    orbit_fid = row.get("orbit_fidelity")
    vinf_fid = row.get("vinf_fidelity")
    same_fid = orbit_fid is not None and vinf_fid is not None and orbit_fid == vinf_fid
    return classify_validation(
        row.get("orbit_source"),
        row.get("vinf_source"),
        same_fidelity=same_fid,
    )


def _corroboration_for(tier: Tier) -> Corroboration:
    # Axis C corroboration label: a cross_validated row pairs two independent
    # sources -> STRONGLY_SOURCED; otherwise a single source -> SINGLE_SOURCED.
    # (Only consequential when Axis A is also machine-confirmed; BRONZE rows are
    # unaffected by this label.)
    return (
        Corroboration.STRONGLY_SOURCED
        if tier is Tier.CROSS_VALIDATED
        else Corroboration.SINGLE_SOURCED
    )


def main() -> None:
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    # Fresh ledger each run (deterministic artifact; idempotent rewrite).
    if LEDGER_PATH.exists():
        LEDGER_PATH.unlink()
    ledger = Ledger(LEDGER_PATH)

    now = datetime.now(UTC).isoformat()
    counts: dict[str, int] = {}
    for row in rows:
        rid = str(row["id"])
        tier = _row_tier(row)
        verdict = run_gauntlet(
            rid,
            agreement=_AXIS_A_AGREEMENT.get(rid),
            provenance_tier=tier,
            corroboration=_corroboration_for(tier),
            known_id=rid,
            notes="catalogue-row gauntlet sweep (#125 Part 3); Axis C + (Aldrin) Axis A only",
        )
        validate_verdict(verdict)  # teeth: refuse any over-claiming verdict
        counts[verdict.tier.value] = counts.get(verdict.tier.value, 0) + 1
        ledger.record(
            LedgerEntry(
                cell_id=rid,
                status="searched",
                n_solutions=0,
                best_dv_kms=None,
                signature_hashes=(),
                validation_level=row.get("validation_level"),
                t_done=now,
                host="gauntlet-sweep",
                verdict_tier=verdict.tier.value,
                verdict_audit={
                    "axis_results": verdict.axis_results,
                    "provenance": verdict.provenance,
                },
            )
        )

    print(f"rows: {len(rows)}")
    print(f"verdict census: {counts}")
    print(f"ledger -> {LEDGER_PATH}")


if __name__ == "__main__":
    main()
