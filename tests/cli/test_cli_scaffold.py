"""M8-UX Phase 0: the CLI entry point exists, is typed, and dispatches.

These exercise the parser/dispatch shell only — no physics. Subprocess smoke
runs use the installed `cyclerfinder` console script (pyproject [project.scripts]).
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from cyclerfinder.cli import build_parser, main


def test_version_flag_exits_zero(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "cyclerfinder" in capsys.readouterr().out.lower()


def test_no_subcommand_prints_usage_and_exits_nonzero(
    capsys: pytest.CaptureFixture[str],
) -> None:
    code = main([])
    assert code == 2  # argparse-style "no command" usage error
    assert "usage" in capsys.readouterr().err.lower()


def test_parser_lists_all_five_subcommands() -> None:
    parser = build_parser()
    # argparse stores subparser choices on the subparsers action.
    sub = next(a for a in parser._actions if a.dest == "command")
    assert sub.choices is not None
    assert set(sub.choices) == {"enumerate", "solve", "discover", "report", "viz"}


def test_console_script_version_subprocess() -> None:
    """The installed console script runs end-to-end (uv-resolved env)."""
    out = subprocess.run(
        [sys.executable, "-m", "cyclerfinder.cli", "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert out.returncode == 0
    assert "cyclerfinder" in out.stdout.lower()
