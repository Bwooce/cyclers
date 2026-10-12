"""#1070: scripts/commit.sh commits by pathspec and reports its own hash.

Each test builds a throwaway git repository with a copy of the wrapper. Expected behaviour (from
the #1070 specification, not from the script): success prints exactly one line
`COMMITTED <short hash> <subject>` naming the commit the call made, and only the named paths are
committed; a refused commit prints `NOT COMMITTED` and the hook's last lines and leaves HEAD and
the index as they were.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WRAPPER = "scripts/commit.sh"

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")


def _run(args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
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
    return subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True, check=False)


def _git(repo: Path, *args: str) -> str:
    proc = _run(["git", *args], cwd=repo)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "config", "core.hooksPath", ".git/hooks")
    dst = tmp_path / WRAPPER
    dst.parent.mkdir(parents=True)
    shutil.copy2(REPO / WRAPPER, dst)
    dst.chmod(0o755)
    (tmp_path / "a.txt").write_text("a\n")
    (tmp_path / "b.txt").write_text("b\n")
    _git(tmp_path, "add", "a.txt", "b.txt", WRAPPER)
    _git(tmp_path, "commit", "-q", "-m", "base")
    return tmp_path


def _hook(repo: Path, name: str, body: str) -> None:
    h = repo / ".git" / "hooks" / name
    h.write_text("#!/bin/sh\n" + body)
    h.chmod(0o755)


def _committed(proc: subprocess.CompletedProcess[str]) -> tuple[str, str]:
    lines = proc.stdout.splitlines()
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert len(lines) == 1, lines
    m = re.fullmatch(r"COMMITTED ([0-9a-f]{7,}) (.*)", lines[0])
    assert m, lines[0]
    return m.group(1), m.group(2)


def test_success_one_line_and_only_named_paths(repo: Path) -> None:
    (repo / "a.txt").write_text("a2\n")
    (repo / "b.txt").write_text("b2\n")
    proc = _run([WRAPPER, "#1: change a\n\nbody line", "a.txt"], cwd=repo)
    short, subject = _committed(proc)
    assert subject == "#1: change a"
    assert _git(repo, "rev-parse", "HEAD").startswith(short)
    assert _git(repo, "show", "--name-only", "--format=", "HEAD") == "a.txt"
    assert _git(repo, "status", "--porcelain") == "M b.txt"


def test_untracked_path_is_committed(repo: Path) -> None:
    (repo / "new.txt").write_text("n\n")
    (repo / "stray.txt").write_text("s\n")
    short, _ = _committed(_run([WRAPPER, "#2: new file", "new.txt"], cwd=repo))
    assert _git(repo, "show", "--name-only", "--format=", short) == "new.txt"
    assert _git(repo, "status", "--porcelain") == "?? stray.txt"


def test_hook_refusal_prints_not_committed_and_hook_tail(repo: Path) -> None:
    head = _git(repo, "rev-parse", "HEAD")
    _hook(repo, "pre-commit", "echo early line\necho 'COMMIT REFUSED: window open'\nexit 1\n")
    (repo / "a.txt").write_text("a3\n")
    (repo / "new.txt").write_text("n\n")
    proc = _run([WRAPPER, "#3: refused", "a.txt", "new.txt"], cwd=repo)
    assert proc.returncode != 0
    lines = proc.stdout.splitlines()
    assert lines[0].startswith("NOT COMMITTED")
    assert "COMMIT REFUSED: window open" in lines
    assert not any(line.startswith("COMMITTED") for line in lines)
    assert _git(repo, "rev-parse", "HEAD") == head
    # the untracked file is not left staged; both edits stay in the working tree
    assert _git(repo, "diff", "--cached", "--name-only") == ""
    assert _git(repo, "diff", "--name-only") == "a.txt"
    assert _git(repo, "ls-files", "--others") == "new.txt"


def test_reports_own_commit_when_head_moves_on(repo: Path) -> None:
    # a post-commit hook commits again before the wrapper's commit returns, as a sibling agent
    # committing in between would: the wrapper must name its own commit, not HEAD
    _hook(
        repo,
        "post-commit",
        'if [ -z "${SIBLING:-}" ]; then echo s >> b.txt; '
        "SIBLING=1 git commit -q -m sibling -- b.txt; fi\n",
    )
    (repo / "a.txt").write_text("a4\n")
    short, subject = _committed(_run([WRAPPER, "#4: mine", "a.txt"], cwd=repo))
    assert subject == "#4: mine"
    assert _git(repo, "log", "-1", "--format=%s") == "sibling"
    assert _git(repo, "log", "-1", "--format=%s", short) == "#4: mine"


@pytest.mark.parametrize(
    "args", [[], ["#5: no paths"], ["", "a.txt"], ["#5: dot", "."], ["#5: flag", "-a"]]
)
def test_usage_and_whole_tree_pathspecs_refused(repo: Path, args: list[str]) -> None:
    (repo / "a.txt").write_text("a5\n")
    head = _git(repo, "rev-parse", "HEAD")
    proc = _run([WRAPPER, *args], cwd=repo)
    assert proc.returncode == 2
    assert "COMMITTED " not in proc.stdout
    assert _git(repo, "rev-parse", "HEAD") == head
