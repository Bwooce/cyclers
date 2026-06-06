"""M8-UX Phase 1: `cyclerfinder enumerate` wires enumerate_cells +
tisserand_feasible. Default body set is the M8 VEM anchor (spec §8 line 152)."""

from __future__ import annotations

import json

import pytest

from cyclerfinder.cli import main


def test_enumerate_vem_defaults_emit_cells(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(
        ["enumerate", "--bodies", "V,E,M", "--l-max", "3", "--k-max", "1", "--format", "json"]
    )
    assert code == 0
    rows = json.loads(capsys.readouterr().out)
    assert rows, "expected at least one VEM cell"
    # every row carries the cell id and a feasibility verdict
    assert {"cell_id", "feasible"} <= set(rows[0])
    assert all(r["cell_id"].startswith("VEM|") for r in rows)


def test_enumerate_feasible_only_filters(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(
        [
            "enumerate",
            "--bodies",
            "V,E,M",
            "--l-max",
            "3",
            "--k-max",
            "1",
            "--vinf-cap",
            "7.0",
            "--feasible-only",
            "--format",
            "json",
        ]
    )
    assert code == 0
    rows = json.loads(capsys.readouterr().out)
    assert all(r["feasible"] is True for r in rows)


def test_enumerate_unknown_body_is_usage_error(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["enumerate", "--bodies", "V,E,Z"])
    assert exc.value.code == 2
    assert "Z" in capsys.readouterr().err
