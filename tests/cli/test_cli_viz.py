"""M8-UX Phase 5: `cyclerfinder viz <kind>` dispatch.

The beat smoke needs matplotlib (skipped without the [viz] extra). The
missing-extra exit-3 path is covered with a monkeypatch so it runs even on a
viz-installed runner."""

from __future__ import annotations

from pathlib import Path

import pytest

from cyclerfinder.cli import main


def test_viz_beat_writes_file(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    out = tmp_path / "b.png"
    code = main(["viz", "beat", "--bodies", "V,E,M", "--out", str(out)])
    assert code == 0
    assert out.exists() and out.stat().st_size > 0


def test_viz_missing_extra_exits_three(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    pytest.importorskip("matplotlib")
    from cyclerfinder.viz import MissingVizExtra, plots

    def _boom(*_a: object, **_k: object) -> None:
        raise MissingVizExtra("matplotlib is required; install '.[viz]'")

    monkeypatch.setattr(plots, "beat_diagram", _boom)
    out = tmp_path / "b.png"
    code = main(["viz", "beat", "--bodies", "V,E,M", "--out", str(out)])
    assert code == 3
    assert "viz" in capsys.readouterr().err.lower()


def test_viz_no_kind_is_usage_error() -> None:
    with pytest.raises(SystemExit) as exc:
        main(["viz"])
    assert exc.value.code == 2
