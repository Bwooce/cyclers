"""M8-UX Phase 4: campaign report assembly. GOLDEN DISCIPLINE — the artifact
keeps sourced facts (catalogue period/sequence/sourced Vinf) separate from
computed results (our optimiser's Vinf/dV/closure). A report never presents a
computed value as sourced."""

from __future__ import annotations

from pathlib import Path

from cyclerfinder.data.ledger import Ledger, LedgerEntry
from cyclerfinder.report import (
    build_campaign_report,
    render_markdown,
)


def _seed_ledger(path: Path) -> None:
    led = Ledger(path)
    led.record(
        LedgerEntry(  # minimal solved entry; fields per ledger.py:73
            cell_id="EM|E-M-E|k2|r00|bss",
            status="solved",
            n_solutions=1,
            best_dv_kms=0.012,
            signature_hashes=("sha1:deadbeef",),
            validation_level="V2",
            t_done="2026-06-06T00:00:00Z",
            host="test",
        )
    )


def test_report_has_sourced_and_computed_sections(tmp_path: Path) -> None:
    ledger = tmp_path / "c.jsonl"
    _seed_ledger(ledger)
    rep = build_campaign_report(ledger)
    assert "sourced" in rep and "computed" in rep
    # the computed dV is under `computed`, never under `sourced`
    cand = rep["candidates"][0]
    assert "best_dv_kms" in cand["computed"]
    assert "best_dv_kms" not in cand.get("sourced", {})


def test_sourced_and_computed_keys_are_disjoint(tmp_path: Path) -> None:
    """Structural golden discipline: the two key sets never overlap."""
    from cyclerfinder.report import _COMPUTED_KEYS, _SOURCED_KEYS

    assert _COMPUTED_KEYS.isdisjoint(_SOURCED_KEYS)


def test_render_markdown_carries_provenance_header(tmp_path: Path) -> None:
    ledger = tmp_path / "c.jsonl"
    _seed_ledger(ledger)
    md = render_markdown(build_campaign_report(ledger))
    assert "# Cycler campaign report" in md
    assert "Sourced" in md and "Computed" in md
    assert "validation_level" in md.lower() or "V2" in md


def test_attach_verdicts_carries_tier(tmp_path: Path) -> None:
    from cyclerfinder.report import attach_verdicts

    ledger = tmp_path / "c.jsonl"
    _seed_ledger(ledger)
    rep = attach_verdicts(build_campaign_report(ledger))
    cand = rep["candidates"][0]
    assert cand["verdict_tier"] in {"GOLD", "SILVER", "BRONZE", "REJECTED"}
