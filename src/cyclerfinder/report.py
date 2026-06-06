# AP_FLAKE8_CLEAN
"""Campaign-report assembly — markdown + JSON, provenance-separated.

GOLDEN DISCIPLINE (spec §11.3/§17, M8-Core golden discipline). The report
artifact keeps **sourced** facts (catalogue-traceable: period / sequence /
sourced V∞) strictly separate from **computed** results (our optimiser's /
ledger's V∞ / ΔV / closure residual / validation level). The split is
*structural*: ``_SOURCED_KEYS`` and ``_COMPUTED_KEYS`` are disjoint by
construction (asserted in tests), so a computed value can never be laundered
into a sourced anchor.

Each candidate also carries its gauntlet :class:`VerdictTier` verbatim
(``--with-verdicts``) so the trust ladder (spec §14, V0-V5) is visible and never
upgraded by the report.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder import __version__
from cyclerfinder.data.ledger import LedgerEntry, LedgerLoader

# Disjoint key partitions — the structural guarantee behind golden discipline.
_COMPUTED_KEYS: frozenset[str] = frozenset(
    {"best_dv_kms", "n_solutions", "validation_level", "signature_hashes"}
)
_SOURCED_KEYS: frozenset[str] = frozenset({"period_yr", "sequence", "vinf_multiset_kms"})


def _now_iso() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def _candidate_from_entry(entry: LedgerEntry) -> dict[str, Any]:
    """Partition one ledger entry into disjoint sourced / computed blocks.

    The ledger holds only computed facts; the ``sourced`` block is left empty
    here (a future catalogue join populates period / sequence / sourced V∞ when
    a cell matched a published row). The two sub-dicts never share a key.
    """
    computed: dict[str, Any] = {
        "best_dv_kms": entry.best_dv_kms,
        "n_solutions": entry.n_solutions,
        "validation_level": entry.validation_level,
        "signature_hashes": list(entry.signature_hashes),
    }
    sourced: dict[str, Any] = {}
    return {
        "cell_id": entry.cell_id,
        "status": entry.status,
        "sourced": sourced,
        "computed": computed,
    }


def build_campaign_report(ledger_path: Path | str) -> dict[str, Any]:
    """Read the ledger and assemble the provenance-separated report dict."""
    entries = list(LedgerLoader(ledger_path))
    candidates = [_candidate_from_entry(e) for e in entries]

    level_counts: dict[str, int] = {}
    n_solved = 0
    for e in entries:
        if e.status == "solved":
            n_solved += 1
        if e.validation_level is not None:
            level_counts[e.validation_level] = level_counts.get(e.validation_level, 0) + 1

    return {
        "generated_at": _now_iso(),
        "ledger": str(ledger_path),
        "provenance": {"tool": "cyclerfinder", "version": __version__},
        "summary": {
            "n_cells": len(entries),
            "n_solved": n_solved,
            "level_counts": level_counts,
        },
        # Top-level markers so a reader (and the tests) can see the discipline
        # is present even when an individual candidate's sourced block is empty.
        "sourced": sorted(_SOURCED_KEYS),
        "computed": sorted(_COMPUTED_KEYS),
        "candidates": candidates,
    }


def attach_verdicts(report: dict[str, Any]) -> dict[str, Any]:
    """Attach a verbatim gauntlet :class:`VerdictTier` to each candidate.

    The verdict is carried **verbatim** — the report never upgrades a tier. The
    only axes derivable from a bare ledger entry are provenance notes; without
    Axis-A machine-confirmation the combiner returns BRONZE (a weak signal),
    which is the honest floor for an unsourced ledger row.
    """
    from cyclerfinder.verify.gauntlet import run_gauntlet

    for cand in report["candidates"]:
        verdict = run_gauntlet(
            candidate_id=cand["cell_id"],
            notes=f"validation_level={cand['computed'].get('validation_level')}",
        )
        cand["verdict_tier"] = verdict.tier.name
        cand["confidence"] = verdict.confidence
    return report


def render_json(report: dict[str, Any]) -> str:
    """Serialise the report to deterministic JSON."""
    return json.dumps(report, indent=2, sort_keys=True)


def render_markdown(report: dict[str, Any]) -> str:
    """Render the report as a markdown document with sourced/computed groups."""
    prov = report["provenance"]
    summary = report["summary"]
    lines: list[str] = []
    lines.append("# Cycler campaign report")
    lines.append("")
    lines.append(f"- Generated: `{report['generated_at']}`")
    lines.append(f"- Ledger: `{report['ledger']}`")
    lines.append(f"- Tool: `{prov['tool']} {prov['version']}`")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Cells: {summary['n_cells']}")
    lines.append(f"- Solved: {summary['n_solved']}")
    if summary["level_counts"]:
        levels = ", ".join(f"{k}: {v}" for k, v in sorted(summary["level_counts"].items()))
        lines.append(f"- Validation levels: {levels}")
    lines.append("")
    lines.append("## Candidates")
    lines.append("")
    lines.append(
        "Each candidate separates **Sourced** (catalogue-traceable) from "
        "**Computed** (our optimiser / ledger) facts; a computed value is never "
        "presented as a sourced anchor."
    )
    lines.append("")
    has_verdicts = any("verdict_tier" in c for c in report["candidates"])
    header = "| cell_id | validation_level | Computed best_dv_kms | Sourced period_yr |"
    sep = "|---|---|---|---|"
    if has_verdicts:
        header += " verdict_tier |"
        sep += "---|"
    lines.append(header)
    lines.append(sep)
    for cand in report["candidates"]:
        computed = cand["computed"]
        sourced = cand.get("sourced", {})
        row = (
            f"| `{cand['cell_id']}` "
            f"| {computed.get('validation_level')} "
            f"| {computed.get('best_dv_kms')} "
            f"| {sourced.get('period_yr', '—')} |"
        )
        if has_verdicts:
            row += f" {cand.get('verdict_tier')} |"
        lines.append(row)
    lines.append("")
    return "\n".join(lines)
