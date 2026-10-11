"""#1069: the ratchet hooks trust a fresh window-suite record and otherwise run pytest.

Expected results are those pre-registered in sec. 4 (T1-T9) of
docs/notes/2026-10-12-1069-hook-trusts-fresh-suite-log.md, written before the code; T10 adds the
shared-tree case of sec. 3. Each test builds a throwaway git repository, copies
scripts/suite_record.py and scripts/ratchet_hook.sh into it, and runs them with the system python3
as the hook does.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = ("scripts/suite_record.py", "scripts/ratchet_hook.sh")
GOOD_LOG = "....................                                    [ 50%]\n" + (
    "....................                                    [100%]\n"
)

pytestmark = pytest.mark.skipif(
    shutil.which("git") is None or shutil.which("python3") is None, reason="git/python3 missing"
)


def _run(
    args: list[str], cwd: Path, env_extra: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
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


def _rec(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return _run(["python3", "scripts/suite_record.py", *args], cwd=repo)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = tmp_path / "r"
    root.mkdir()
    _run(["git", "init", "-q", str(root)], cwd=root)
    for rel in SCRIPTS:
        dst = root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / rel, dst)
        dst.chmod(0o755)
    (root / "data").mkdir()
    (root / "data/catalogue.yaml").write_text("- id: a\n")
    (root / "data/OUTSTANDING.md").write_text("notes\n")
    (root / "tests/data").mkdir(parents=True)
    (root / "tests/data/test_x.py").write_text("def test_x():\n    pass\n")
    _run(["git", "add", "-A"], cwd=root)
    assert _run(["git", "commit", "-q", "-m", "init"], cwd=root).returncode == 0
    return root


def _logs(tmp: Path, text: str = GOOD_LOG) -> list[str]:
    log = tmp / "c1.log"
    log.write_text(text)
    return ["--log", f"{log}:0"]


def _green_record(repo: Path, covers: tuple[str, ...] = ("tests/data",)) -> None:
    assert _rec(repo, "start").returncode == 0
    fin = _rec(repo, "finish", *_logs(repo.parent), "--covers", *covers)
    assert fin.returncode == 0, fin.stdout + fin.stderr
    assert "OK" in fin.stdout


def _record(repo: Path) -> dict[str, object]:
    data: dict[str, object] = json.loads((repo / ".git/CYCLERS_SUITE_OK").read_text())
    return data


def test_t1_no_record_runs(repo: Path) -> None:
    proc = _rec(repo, "check", "tests/data", "-q")
    assert proc.returncode == 1
    assert "no window-suite record" in proc.stdout


def test_t2_green_record_on_same_tree_is_trusted(repo: Path) -> None:
    _green_record(repo)
    rec = _record(repo)
    assert rec["ok"] is True and rec["failed"] == []
    proc = _rec(repo, "check", "tests/data/test_x.py", "-q")
    assert proc.returncode == 0, proc.stdout
    assert "ratchets: trusted window-suite log" in proc.stdout


def test_t3_covered_edit_after_finish_runs(repo: Path) -> None:
    _green_record(repo)
    (repo / "data/catalogue.yaml").write_text("- id: a\n- id: b\n")
    _run(["git", "add", "data/catalogue.yaml"], cwd=repo)
    proc = _rec(repo, "check", "tests/data")
    assert proc.returncode == 1
    assert "differs" in proc.stdout


def test_t4_stale_record_runs(repo: Path) -> None:
    _green_record(repo)
    rec = _record(repo)
    old = dt.datetime.now(dt.timezone.utc) - dt.timedelta(hours=1, minutes=1)  # noqa: UP017
    rec["started"] = old.isoformat()
    (repo / ".git/CYCLERS_SUITE_OK").write_text(json.dumps(rec))
    proc = _rec(repo, "check", "tests/data")
    assert proc.returncode == 1
    assert "stale" in proc.stdout


@pytest.mark.parametrize(
    ("text", "code", "reason"),
    [
        (GOOD_LOG + "FAILED tests/data/test_x.py::test_x - assert 0\n", 0, "FAILED/ERROR"),
        (GOOD_LOG + "ERROR tests/data/test_x.py - ImportError\n", 0, "FAILED/ERROR"),
        (GOOD_LOG, 1, "exit status 1"),
        ("....                                    [ 40%]\n", 0, "[100%]"),
    ],
)
def test_t5_bad_suite_is_not_ok(repo: Path, text: str, code: int, reason: str) -> None:
    assert _rec(repo, "start").returncode == 0
    log = repo.parent / "c1.log"
    log.write_text(text)
    fin = _rec(repo, "finish", "--log", f"{log}:{code}", "--covers", "tests/data")
    assert fin.returncode == 1
    assert reason in fin.stdout
    assert _record(repo)["ok"] is False
    proc = _rec(repo, "check", "tests/data")
    assert proc.returncode == 1
    assert "record not ok" in proc.stdout


def test_t6_covered_edit_during_suite_is_not_ok(repo: Path) -> None:
    assert _rec(repo, "start").returncode == 0
    (repo / "tests/data/test_x.py").write_text("def test_x():\n    assert 1\n")
    fin = _rec(repo, "finish", *_logs(repo.parent), "--covers", "tests/data")
    assert fin.returncode == 1
    assert "changed during the suite" in fin.stdout


def test_t7_target_outside_covers_runs(repo: Path) -> None:
    _green_record(repo, covers=("tests/data",))
    proc = _rec(repo, "check", "tests/data/test_x.py", "tests/verify/test_y.py::test_z")
    assert proc.returncode == 1
    assert "did not cover tests/verify/test_y.py::test_z" in proc.stdout


def test_t8_excluded_md_edit_still_trusted(repo: Path) -> None:
    _green_record(repo)
    (repo / "data/OUTSTANDING.md").write_text("notes, edited\n")
    _run(["git", "add", "data/OUTSTANDING.md"], cwd=repo)
    proc = _rec(repo, "check", "tests/data")
    assert proc.returncode == 0, proc.stdout


def _install_hook_with_fake_uv(repo: Path) -> dict[str, str]:
    hook = repo / ".git/hooks/pre-commit"
    hook.parent.mkdir(parents=True, exist_ok=True)
    hook.write_text("#!/usr/bin/env bash\nexec scripts/ratchet_hook.sh tests/data -q\n")
    hook.chmod(0o755)
    fake = repo.parent / "bin"
    fake.mkdir()
    uv = fake / "uv"
    uv.write_text(f'#!/usr/bin/env bash\ntouch "{repo.parent}/pytest_ran"\nexit 0\n')
    uv.chmod(0o755)
    return {"PATH": f"{fake}:{os.environ['PATH']}"}


def test_t9_real_commit_skips_or_runs(repo: Path) -> None:
    env = _install_hook_with_fake_uv(repo)
    marker = repo.parent / "pytest_ran"

    # The window edit is in the working tree when the suite runs; the commit is by pathspec.
    (repo / "data/catalogue.yaml").write_text("- id: a\n- id: b\n")
    _green_record(repo)
    ok = _run(
        ["git", "commit", "-q", "-m", "w", "--", "data/catalogue.yaml"], cwd=repo, env_extra=env
    )
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert "trusted window-suite log" in ok.stdout + ok.stderr
    assert not marker.exists()

    # Edited after the suite: the hook runs pytest (the fake uv leaves the marker).
    (repo / "data/catalogue.yaml").write_text("- id: a\n- id: b\n- id: c\n")
    again = _run(
        ["git", "commit", "-q", "-m", "w2", "--", "data/catalogue.yaml"], cwd=repo, env_extra=env
    )
    assert again.returncode == 0, again.stdout + again.stderr
    assert marker.exists()


def test_t10_other_agents_unstaged_covered_edit_runs(repo: Path) -> None:
    # An unstaged edit by someone else was in the tree the suite read, but is not committed.
    env = _install_hook_with_fake_uv(repo)
    (repo / "data/catalogue.yaml").write_text("- id: a\n- id: b\n")
    (repo / "tests/data/test_x.py").write_text("def test_x():\n    assert 1\n")
    _green_record(repo)
    proc = _run(
        ["git", "commit", "-q", "-m", "w", "--", "data/catalogue.yaml"], cwd=repo, env_extra=env
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (repo.parent / "pytest_ran").exists()
