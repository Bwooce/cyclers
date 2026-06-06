"""M8-UX Phase 3: `cyclerfinder discover` wires data.discover.discover into a
ledger and surfaces V-level gates. VEM/ephemeris/V3 paths are slow or M-ED."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cyclerfinder.cli import main


def test_discover_em_writes_ledger(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    ledger = tmp_path / "em.jsonl"
    code = main(
        [
            "discover",
            "--bodies",
            "E,M",
            "--k",
            "2",
            "--vinf-cap",
            "7.0",
            "--ledger",
            str(ledger),
            "--l-max",
            "3",
            "--max-cells",
            "3",
            "--n-starts",
            "2",
            "--no-de",
            "--format",
            "json",
        ]
    )
    assert code in (0, 5)  # 0 = solved cells; 5 = none solved (still valid)
    assert ledger.exists()  # the ledger is always written
    summary = json.loads(capsys.readouterr().out)
    assert {"cells_searched", "cells_solved", "ledger"} <= set(summary)


def test_discover_enable_v3_flag_accepted(tmp_path: Path) -> None:
    """3.3 [M-ED]: --enable-v3 plumbs through; flag is accepted + passed."""
    ledger = tmp_path / "em_v3.jsonl"
    code = main(
        [
            "discover",
            "--bodies",
            "E,M",
            "--k",
            "2",
            "--vinf-cap",
            "7.0",
            "--ledger",
            str(ledger),
            "--l-max",
            "3",
            "--max-cells",
            "2",
            "--n-starts",
            "1",
            "--no-de",
            "--enable-v3",
            "--format",
            "json",
        ]
    )
    assert code in (0, 5)
    assert ledger.exists()


@pytest.mark.skip(
    reason="OPEN RESEARCH (was M-ED-BLOCKED): M-ED landed but VEM ballistic "
    "closure to the sourced Jones multisets does NOT converge — the 2026-06-06 "
    "dense 16-core scan (2816 points/row, task #110) floors at max-V_inf "
    "~17.9-18.5 km/s vs sourced 2.4-7.0, zero bend-feasible. V3 promotion for "
    "VEM therefore cannot occur yet. The 3D-inclination hypothesis was tested "
    "and REFUTED (task #120): the same grid on the inclined-circular backend "
    "moved the floor only ~0.1-0.4 km/s (still zero bend-feasible) and a DE440 "
    "control floors at ~18.2. Un-skip when a corrector variant (multi-arc-per-leg "
    "seeding or eccentric intermediate flybys — NOT 3D-inclination-only) reaches "
    "the family; see data/OUTSTANDING.md M-ED open-research entry."
)
@pytest.mark.slow
def test_discover_vem_reaches_v3(tmp_path: Path) -> None:  # pragma: no cover - skipped
    ledger = tmp_path / "vem.jsonl"
    code = main(
        [
            "discover",
            "--bodies",
            "V,E,M",
            "--k",
            "3",
            "--ledger",
            str(ledger),
            "--fidelity",
            "ephemeris",
            "--enable-v3",
            "--priority-date",
            "2032-01-01",
            "--vinf-targets",
            "E=5.65,M=3.05",
            "--format",
            "json",
        ]
    )
    assert code == 0
