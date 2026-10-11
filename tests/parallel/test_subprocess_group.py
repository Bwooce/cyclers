"""#1068: run_group kills the whole process group on a timeout (no orphaned grandchild).

The child stands in for ``uv run``: it starts a grandchild (the "python" that would write a result),
records the grandchild's pid, and waits. The control shows the defect with plain
``subprocess.run(timeout=...)``: the grandchild survives.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from cyclerfinder.parallel.subprocess_group import run_group

_CHILD = (
    "import subprocess, sys, time\n"
    "g = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
    "open(sys.argv[1], 'w').write(str(g.pid))\n"
    "time.sleep(60)\n"
)


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # a zombie still answers kill(0); treat "not running" via waitpid when it is ours
    return True


def _grandchild_pid(path: Path) -> int:
    for _ in range(100):
        if path.exists() and path.read_text().strip():
            return int(path.read_text())
        time.sleep(0.05)
    raise AssertionError("grandchild pid not written")


def _gone(pid: int, wait_s: float = 3.0) -> bool:
    t0 = time.time()
    while time.time() - t0 < wait_s:
        if not _alive(pid):
            return True
        time.sleep(0.05)
    return False


def test_run_group_kills_the_grandchild(tmp_path: Path) -> None:
    pidfile = tmp_path / "g.pid"
    with pytest.raises(subprocess.TimeoutExpired):
        run_group([sys.executable, "-c", _CHILD, str(pidfile)], timeout=3)
    assert _gone(_grandchild_pid(pidfile)), "grandchild survived run_group's timeout"


def test_control_plain_run_leaves_the_grandchild(tmp_path: Path) -> None:
    """The defect run_group fixes (positive control for the test above)."""
    pidfile = tmp_path / "g.pid"
    with pytest.raises(subprocess.TimeoutExpired):
        subprocess.run([sys.executable, "-c", _CHILD, str(pidfile)], timeout=3)
    pid = _grandchild_pid(pidfile)
    try:
        assert _alive(pid), "control: plain subprocess.run should leave the grandchild running"
    finally:
        os.kill(pid, 9)


def test_run_group_normal_exit_and_capture() -> None:
    r = run_group([sys.executable, "-c", "print('ok')"], timeout=30, capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == "ok"
    with pytest.raises(subprocess.CalledProcessError):
        run_group([sys.executable, "-c", "raise SystemExit(3)"], timeout=30, check=True)
