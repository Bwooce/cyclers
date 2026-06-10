"""CR3BP backfill script: Earth-Moon catalogue rows.

Runs the STM periodic-orbit corrector on the three catalogue rows that carry
model_assumption: cr3bp (arenstorf-em-figure8-1963, genova-aldrin-2015-em-3petal-cycler,
wittal-2022-em-cycler-family).

Usage:
    uv run python scripts/cr3bp_backfill.py
    uv run python scripts/cr3bp_backfill.py --report /tmp/cr3bp_backfill.txt

Outputs:
    docs/notes/2026-06-10-cr3bp-backfill-results.md  (always written)
    --report file if supplied (plain text runlog)

IMPORTANT: NO catalogue writeback.  This script is read-only w.r.t. data/catalogue.yaml
and validate.py.  The results note contains PROPOSED fields for human review only.
"""

from __future__ import annotations

import argparse
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# ---------------------------------------------------------------------------
# Earth-Moon physical CR3BP system (from the satellite registry)
# ---------------------------------------------------------------------------
EM_SYSTEM = cr3bp.cr3bp_system("Earth", "Moon")
EM_MU_PHYSICAL = EM_SYSTEM.mu  # ~0.01200 (JPL GMs)

# ---------------------------------------------------------------------------
# Arenstorf sourced golden IC
# Source: Arenstorf 1963 / Hairer, Nørsett & Wanner "Solving ODEs I",
#         Springer (1993 revised ed.), p. 129, test problem B5.
# mu = 0.012277471 (the value used in the Hairer reference, slightly different
# from the physical Earth-Moon mu; this is the standard CR3BP test-problem value).
# x0 = 0.994, vy0 = -2.0015851063790825, T = 17.0652165601579625.
# ---------------------------------------------------------------------------
ARENSTORF_MU = 0.012277471
ARENSTORF_X0 = 0.994
ARENSTORF_VY0 = -2.0015851063790825
ARENSTORF_PERIOD = 17.0652165601579625


def run_arenstorf() -> dict[str, object]:
    """Backfill arenstorf-em-figure8-1963 using the published sourced IC."""
    sysm = cr3bp.CR3BPSystem(
        mu=ARENSTORF_MU,
        primary="Earth",
        secondary="Moon",
        l_km=EM_SYSTEM.l_km,
        t_s=EM_SYSTEM.t_s,
    )
    s0 = np.array([ARENSTORF_X0, 0.0, 0.0, 0.0, ARENSTORF_VY0, 0.0])
    res = cp.correct_periodic(sysm, s0, ARENSTORF_PERIOD)
    return {
        "row_id": "arenstorf-em-figure8-1963",
        "sourced_ic": True,
        "ic_source": "Hairer, Nørsett & Wanner 'Solving ODEs I' p.129 (B5); Arenstorf 1963",
        "mu_used": ARENSTORF_MU,
        "initial_state_nd": s0.tolist(),
        "period_guess_nd": ARENSTORF_PERIOD,
        "converged": res.converged,
        "closure_residual": res.closure_residual,
        "period_nd": res.period,
        "jacobi": res.jacobi,
        "state0_nd": res.state0.tolist(),
        "lunit_km": EM_SYSTEM.l_km,
        "tunit_s": EM_SYSTEM.t_s,
    }


def run_genova() -> dict[str, object]:
    """Attempt backfill of genova-aldrin-2015-em-3petal-cycler.

    The catalogue row carries state_nd: null — no sourced IC is available from
    the accessible NTRS abstract (full AAS-15 PDF was inaccessible at ingest).
    Per the plan's honesty rules we CANNOT fabricate an IC; this row is marked
    NO_SOURCED_IC and skipped.
    """
    return {
        "row_id": "genova-aldrin-2015-em-3petal-cycler",
        "sourced_ic": False,
        "reason": (
            "No published initial conditions available from the accessible NTRS abstract "
            "(NTRS 20150018049). Full AAS-15 PDF was inaccessible at ingest; state_nd is "
            "null in the catalogue row. Per honesty rules, NO IC is fabricated. "
            "Backfill requires the full Genova & Aldrin 2015 paper PDF."
        ),
        "converged": None,
        "closure_residual": None,
        "period_nd": None,
        "jacobi": None,
        "state0_nd": None,
    }


def run_wittal() -> dict[str, object]:
    """Attempt backfill of wittal-2022-em-cycler-family.

    The catalogue row carries state_nd: null — no sourced IC is available from
    the accessible NTRS abstract (full IAC-22-C1.6.6 PDF was inaccessible at ingest).
    This is a family seed entry; per honesty rules we CANNOT fabricate an IC.
    """
    return {
        "row_id": "wittal-2022-em-cycler-family",
        "sourced_ic": False,
        "reason": (
            "No published initial conditions available from the accessible NTRS abstract "
            "(NTRS 20220013595 / IAC-22-C1.6.6). The catalogue row is a family seed "
            "(state_nd: null). Per honesty rules, NO IC is fabricated. "
            "Backfill requires the full Wittal, Miaule & Asher 2022 paper PDF."
        ),
        "converged": None,
        "closure_residual": None,
        "period_nd": None,
        "jacobi": None,
        "state0_nd": None,
    }


def build_report(results: list[dict[str, object]], elapsed_s: float) -> str:
    """Build the markdown results note."""
    ts = datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines: list[str] = []
    lines.append("# CR3BP Earth-Moon Backfill Results")
    lines.append("")
    lines.append(f"**Run timestamp:** {ts}")
    lines.append("**Script:** `scripts/cr3bp_backfill.py`")
    lines.append(f"**Elapsed:** {elapsed_s:.1f} s")
    lines.append("")
    lines.append("**Status:** NO catalogue writeback. Results are PROPOSED for human review only.")
    lines.append("")
    lines.append("---")
    lines.append("")

    for r in results:
        row_id = str(r["row_id"])
        lines.append(f"## Row: `{row_id}`")
        lines.append("")

        if not r["sourced_ic"]:
            lines.append("**Backfill outcome: NO_SOURCED_IC — skipped**")
            lines.append("")
            lines.append(f"> {r['reason']}")
            lines.append("")
            lines.append("No CR3BP fields can be proposed for this row.")
            lines.append("")
            lines.append("---")
            lines.append("")
            continue

        converged = r["converged"]
        closure = r["closure_residual"]
        lines.append(f"**Backfill outcome: {'CONVERGED' if converged else 'DID NOT CONVERGE'}**")
        lines.append("")
        lines.append(f"- IC source: {r['ic_source']}")
        lines.append(f"- μ used: `{r['mu_used']}`")
        lines.append(f"- Initial state (nd): `{r['initial_state_nd']}`")
        lines.append(f"- Period guess (nd): `{r['period_guess_nd']}`")
        lines.append("")
        lines.append("**Corrector output:**")
        lines.append("")
        lines.append("| Field | Value |")
        lines.append("|---|---|")
        lines.append(f"| converged | `{converged}` |")
        lines.append(f"| closure_residual | `{closure:.3e}` |")
        lines.append(f"| period_nd (corrected) | `{r['period_nd']:.15g}` |")
        lines.append(f"| jacobi_constant | `{r['jacobi']:.10g}` |")
        lines.append(f"| state0_nd (corrected) | `{r['state0_nd']}` |")
        lines.append(f"| lunit_km | `{r['lunit_km']}` |")
        lines.append(f"| tunit_s | `{r['tunit_s']:.6g}` |")
        lines.append("")

        if converged:
            lines.append("**PROPOSED `orbit_elements.cr3bp` fields (review-gated, NO writeback):**")
            lines.append("")
            lines.append("```yaml")
            lines.append("    cr3bp:")
            lines.append(f"      mass_ratio: {r['mu_used']}")
            lines.append(f"      jacobi_constant: {r['jacobi']:.10g}")
            lines.append(f"      period_nd: {r['period_nd']:.15g}")
            lines.append(f"      state_nd: {r['state0_nd']}")
            lines.append(f"      lunit_km: {r['lunit_km']}")
            lines.append(f"      tunit_s: {r['tunit_s']:.6g}")
            lines.append("```")
            lines.append("")
            lines.append("**PROPOSED `_LEVEL_EVIDENCE` line (Arenstorf row only):**")
            lines.append("")
            lines.append(
                "The Arenstorf IC is sourced from Hairer et al. (1993), a citable "
                "published reference, and the corrector converges to closure < 1e-10. "
                "This meets the criteria for promotion from V0 (citation-only) to a "
                "higher validation level once the Jacobi constant is cross-checked against "
                "an independent source (e.g. the JPL three-body periodic-orbit catalog)."
            )
            lines.append("")
            lines.append(
                "_Proposed level: **V1** (computed from sourced IC; pending independent "
                "Jacobi cross-check for V2/V3)._"
            )
            lines.append("")

        lines.append("---")
        lines.append("")

    lines.append("## Summary")
    lines.append("")
    converged_rows = [r["row_id"] for r in results if r.get("converged")]
    no_ic_rows = [r["row_id"] for r in results if not r["sourced_ic"]]
    lines.append(f"- **Converged (sourced IC + periodic orbit found):** {len(converged_rows)}")
    if converged_rows:
        for rid in converged_rows:
            lines.append(f"  - `{rid}`")
    lines.append(f"- **No sourced IC (skipped, not fabricated):** {len(no_ic_rows)}")
    if no_ic_rows:
        for rid in no_ic_rows:
            lines.append(f"  - `{rid}`")
    lines.append("")
    lines.append(
        "**NO writeback to `data/catalogue.yaml` or `validate.py`.**  "
        "Promotion is review-gated (separate step)."
    )
    lines.append("")

    return "\n".join(lines)


def build_runlog(results: list[dict[str, object]], elapsed_s: float) -> str:
    """Build a plain-text runlog."""
    lines: list[str] = []
    ts = datetime.now(tz=UTC).isoformat()
    lines.append(f"cr3bp_backfill runlog  {ts}  elapsed={elapsed_s:.1f}s")
    lines.append("=" * 70)
    for r in results:
        lines.append(f"row: {r['row_id']}")
        if not r["sourced_ic"]:
            lines.append("  status: NO_SOURCED_IC")
            lines.append(f"  reason: {str(r['reason'])[:120]}...")
        else:
            status = "CONVERGED" if r["converged"] else "DID_NOT_CONVERGE"
            lines.append(f"  status: {status}")
            lines.append(f"  closure_residual: {r['closure_residual']:.3e}")
            lines.append(f"  period_nd: {r['period_nd']:.10g}")
            lines.append(f"  jacobi: {r['jacobi']:.10g}")
        lines.append("")
    lines.append("NO catalogue writeback performed.")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CR3BP Earth-Moon catalogue backfill (no writeback)"
    )
    parser.add_argument("--report", type=Path, default=None, help="Write plain-text runlog here")
    args = parser.parse_args()

    t0 = datetime.now(tz=UTC)

    print("cr3bp_backfill: running Earth-Moon CR3BP corrector on 3 catalogue rows...")
    print(f"  Earth-Moon physical mu (from registry): {EM_MU_PHYSICAL:.8f}")
    print()

    results: list[dict[str, object]] = []

    print("  [1/3] arenstorf-em-figure8-1963  (sourced IC: YES)")
    r1 = run_arenstorf()
    results.append(r1)
    if r1["converged"]:
        print(
            f"    -> CONVERGED  closure={r1['closure_residual']:.3e}  "
            f"jacobi={r1['jacobi']:.8g}  period={r1['period_nd']:.10g}"
        )
    else:
        print(f"    -> DID NOT CONVERGE  closure={r1['closure_residual']:.3e}")

    print("  [2/3] genova-aldrin-2015-em-3petal-cycler  (sourced IC: NO)")
    r2 = run_genova()
    results.append(r2)
    print("    -> NO_SOURCED_IC — skipped (not fabricated)")

    print("  [3/3] wittal-2022-em-cycler-family  (sourced IC: NO)")
    r3 = run_wittal()
    results.append(r3)
    print("    -> NO_SOURCED_IC — skipped (not fabricated)")

    t1 = datetime.now(tz=UTC)
    elapsed = (t1 - t0).total_seconds()

    md = build_report(results, elapsed)
    note_path = Path("docs/notes/2026-06-10-cr3bp-backfill-results.md")
    note_path.write_text(md, encoding="utf-8")
    print(f"\nResults note written: {note_path}")

    if args.report is not None:
        runlog = build_runlog(results, elapsed)
        args.report.write_text(runlog, encoding="utf-8")
        print(f"Runlog written:       {args.report}")

    print("\nNO catalogue writeback performed.")


if __name__ == "__main__":
    main()
