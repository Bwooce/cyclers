"""#964 -- scripts/repo_sync_status.py reports ahead/behind and refuses a missing upstream.

Offline: builds a bare "remote" and clones in a temp directory; no network.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "repo_sync_status.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("repo_sync_status", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses resolve their module via sys.modules
    spec.loader.exec_module(mod)
    return mod


rs = _load()


def _run(cwd: Path, *args: str) -> None:
    subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True)


def _commit(repo: Path, name: str) -> None:
    (repo / name).write_text(name)
    _run(repo, "git", "add", name)
    _run(repo, "git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", name)


@pytest.fixture
def remote_and_clones(tmp_path: Path) -> tuple[Path, Path, Path]:
    remote = tmp_path / "remote.git"
    _run(tmp_path, "git", "init", "-q", "--bare", "-b", "main", str(remote))
    a = tmp_path / "a"
    _run(tmp_path, "git", "clone", "-q", str(remote), str(a))
    _run(a, "git", "checkout", "-q", "-b", "main")
    _commit(a, "one")
    _run(a, "git", "push", "-q", "-u", "origin", "main")
    b = tmp_path / "b"
    _run(tmp_path, "git", "clone", "-q", str(remote), str(b))
    return remote, a, b


def test_up_to_date_clone_is_ok(remote_and_clones: tuple[Path, Path, Path]) -> None:
    _, _, b = remote_and_clones
    st = rs.repo_status(b, fetch=True)
    assert st.ok and (st.ahead, st.behind) == (0, 0)


def test_behind_is_detected_only_after_fetch(remote_and_clones: tuple[Path, Path, Path]) -> None:
    _, a, b = remote_and_clones
    _commit(a, "two")
    _commit(a, "three")
    _run(a, "git", "push", "-q")
    assert rs.repo_status(b, fetch=False).behind == 0  # stale refs: looks fine
    st = rs.repo_status(b, fetch=True)
    assert not st.ok and st.behind == 2
    assert "behind" in (st.problem or "")
    assert (b / "two").exists() is False  # fetch never touches the working tree


def test_ahead_and_dirty_are_reported(remote_and_clones: tuple[Path, Path, Path]) -> None:
    _, _, b = remote_and_clones
    _commit(b, "local")
    (b / "scratch.txt").write_text("x")
    st = rs.repo_status(b, fetch=False)
    assert st.ok and st.ahead == 1 and st.dirty == 1


def test_missing_upstream_is_a_problem(tmp_path: Path) -> None:
    repo = tmp_path / "solo"
    _run(tmp_path, "git", "init", "-q", "-b", "main", str(repo))
    _commit(repo, "one")
    st = rs.repo_status(repo, fetch=False)
    assert not st.ok and st.upstream is None
    assert "no upstream" in (st.problem or "")


def test_missing_repo_and_exit_code(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert rs.main([str(tmp_path / "nope"), "--no-fetch"]) == 1
    assert "not a git checkout" in capsys.readouterr().out
