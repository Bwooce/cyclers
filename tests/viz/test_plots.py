"""M8-UX Phase 5: optional-extra viz. Headless (Agg), file-output, no display.
Skipped wholesale when matplotlib (the [viz] extra) is not installed."""

from __future__ import annotations

from pathlib import Path

import pytest

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")  # before pyplot import anywhere

from cyclerfinder.viz import plots  # noqa: E402

_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def test_beat_diagram_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "beat.png"
    plots.beat_diagram(["V", "E", "M"], out)
    assert out.exists() and out.stat().st_size > 0
    assert out.read_bytes()[:8] == _PNG_MAGIC


def test_porkchop_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "pork.png"
    plots.porkchop(
        "E",
        "M",
        epoch_range=("2032-01-01", "2034-01-01"),
        tof_range=(100.0, 400.0),
        out_path=out,
        n_epoch=6,
        n_tof=6,
    )
    assert out.exists() and out.stat().st_size > 0
    assert out.read_bytes()[:8] == _PNG_MAGIC


def test_trajectory_writes_png(tmp_path: Path) -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_idealized
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    result = optimise_cell_idealized(
        cell, Ephemeris(model="circular"), vinf_cap=7.0, n_starts=1, use_de=False
    )
    out = tmp_path / "traj.png"
    plots.trajectory(result, out)
    assert out.exists() and out.stat().st_size > 0
    assert out.read_bytes()[:8] == _PNG_MAGIC
