"""#1063: the catalogue-window lock refuses commits while a window is open.

Each test builds a throwaway git repository, copies the three window scripts into it, and runs
them as the hook and the lead would. Expected behaviour (from the #1063 specification, not from
the scripts): no lock -> commit allowed; lock -> refused with the lead, time and reason named;
lock plus CYCLERS_WINDOW_LEAD=1 -> allowed; the githooks wrapper refuses before any installed
pre-commit hook runs.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = ("scripts/window_lock_check.sh", "scripts/window.sh", "scripts/githooks/pre-commit")

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")


def _run(
    args: list[str], cwd: Path, env_extra: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.pop("CYCLERS_WINDOW_LEAD", None)
    env.update(
        {
            "GIT_AUTHOR_NAME": "t",
            "GIT_AUTHOR_EMAIL": "t@example.invalid",
            "GIT_COMMITTER_NAME": "t",
            "GIT_COMMITTER_EMAIL": "t@example.invalid",
        }
    )
    if env_extra:
        env.update(env_extra)
    return subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True, check=False)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    _run(["git", "init", "-q", str(tmp_path)], cwd=tmp_path)
    for rel in SCRIPTS:
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / rel, dst)
        dst.chmod(0o755)
    (tmp_path / "a.txt").write_text("a\n")
    return tmp_path


def _lock(repo: Path) -> Path:
    return repo / ".git" / "CYCLERS_WINDOW_OPEN"


def test_no_lock_allows(repo: Path) -> None:
    proc = _run(["scripts/window_lock_check.sh"], cwd=repo)
    assert proc.returncode == 0, proc.stderr


def test_lock_refuses_and_names_lead_time_reason(repo: Path) -> None:
    opened = _run(["scripts/window.sh", "open", "window 4: ev rows"], cwd=repo)
    assert opened.returncode == 0, opened.stderr
    assert _lock(repo).read_text().startswith("lead ")
    proc = _run(["scripts/window_lock_check.sh"], cwd=repo)
    assert proc.returncode == 1
    assert "COMMIT REFUSED" in proc.stderr
    assert "opened by: lead" in proc.stderr
    assert "window 4: ev rows" in proc.stderr
    stamp = _lock(repo).read_text().split()[1]
    assert stamp in proc.stderr


def test_lead_variable_allows_inside_window(repo: Path) -> None:
    _run(["scripts/window.sh", "open", "w"], cwd=repo)
    proc = _run(["scripts/window_lock_check.sh"], cwd=repo, env_extra={"CYCLERS_WINDOW_LEAD": "1"})
    assert proc.returncode == 0, proc.stderr
    other = _run(
        ["scripts/window_lock_check.sh"], cwd=repo, env_extra={"CYCLERS_WINDOW_LEAD": "yes"}
    )
    assert other.returncode == 1


def test_close_reopens_commits_and_double_open_refused(repo: Path) -> None:
    assert _run(["scripts/window.sh", "open", "w"], cwd=repo).returncode == 0
    assert _run(["scripts/window.sh", "open", "w2"], cwd=repo).returncode == 1
    assert _run(["scripts/window.sh", "close"], cwd=repo).returncode == 0
    assert not _lock(repo).exists()
    assert _run(["scripts/window_lock_check.sh"], cwd=repo).returncode == 0
    assert _run(["scripts/window.sh", "close"], cwd=repo).returncode == 1


def test_githooks_wrapper_blocks_a_real_commit_before_the_installed_hook(repo: Path) -> None:
    # An installed hook that would leave a marker if it ran: a refused commit must not reach it.
    installed = repo / ".git" / "hooks" / "pre-commit"
    installed.parent.mkdir(parents=True, exist_ok=True)
    installed.write_text("#!/usr/bin/env bash\ntouch .git/installed_hook_ran\nexit 0\n")
    installed.chmod(0o755)
    _run(["git", "config", "core.hooksPath", "scripts/githooks"], cwd=repo)
    _run(["git", "add", "a.txt"], cwd=repo)

    _run(["scripts/window.sh", "open", "w"], cwd=repo)
    refused = _run(["git", "commit", "-q", "-m", "x", "--", "a.txt"], cwd=repo)
    assert refused.returncode != 0
    assert "COMMIT REFUSED" in refused.stderr
    assert not (repo / ".git" / "installed_hook_ran").exists()
    assert _run(["git", "rev-parse", "--verify", "-q", "HEAD"], cwd=repo).returncode != 0

    lead = _run(
        ["git", "commit", "-q", "-m", "x", "--", "a.txt"],
        cwd=repo,
        env_extra={"CYCLERS_WINDOW_LEAD": "1"},
    )
    assert lead.returncode == 0, lead.stderr
    assert (repo / ".git" / "installed_hook_ran").exists()

    _run(["scripts/window.sh", "close"], cwd=repo)
    (repo / "a.txt").write_text("b\n")
    after = _run(["git", "commit", "-q", "-m", "y", "--", "a.txt"], cwd=repo)
    assert after.returncode == 0, after.stderr
