"""Process-lifecycle test: orphaned workers are reaped on parent SIGKILL (#479).

The PR_SET_PDEATHSIG initializer installed in every worker process ensures the
kernel delivers SIGTERM to the worker the instant its parent exits — even when
the parent is SIGKILLed and no Python handler can run.

Test strategy
-------------
1. Spawn a subprocess that runs a long parallel_sweep (cells that sleep).
2. Collect the subprocess's child worker PIDs via /proc (Linux-only, no psutil).
3. SIGKILL the subprocess.
4. Poll up to ~3 s; assert that every worker PID has exited.

The test is skipped on non-Linux because:
  (a) PR_SET_PDEATHSIG is a Linux kernel feature, and
  (b) the /proc-based PID enumeration is Linux-only.
"""

from __future__ import annotations

import contextlib
import os
import pathlib
import signal
import subprocess
import sys
import textwrap
import time

import pytest

# ---------------------------------------------------------------------------
# Platform guard
# ---------------------------------------------------------------------------

pytestmark = pytest.mark.skipif(
    sys.platform != "linux",
    reason="PR_SET_PDEATHSIG and /proc PID enumeration are Linux-only",
)

# ---------------------------------------------------------------------------
# The subprocess script (executed as a child process by the test)
# ---------------------------------------------------------------------------

# This script is run as a separate Python process.  It starts a parallel_sweep
# over sleep cells so it will still be running when the test SIGKILLs it.
# It writes its own PID to a file so the test can find it, then writes each
# worker PID (one per line) to a second file once the pool is started.
_WORKER_SCRIPT = textwrap.dedent(
    """\
    import os, sys, time, multiprocessing

    pid_file = sys.argv[1]
    worker_pid_file = sys.argv[2]
    src = sys.argv[3]

    # Write our own PID so the test knows we started.
    with open(pid_file, "w") as f:
        f.write(str(os.getpid()))

    sys.path.insert(0, src)

    from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep

    def _long_sleep(x):
        # Write this worker's PID into the shared file once it starts.
        with open(worker_pid_file, "a") as f:
            f.write(str(os.getpid()) + "\\n")
        time.sleep(30)
        return x

    cfg = ParallelSweepConfig(n_workers=2, backend="loky")
    parallel_sweep(list(range(4)), _long_sleep, config=cfg)
    """
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _find_src_path() -> str:
    """Return the absolute path to the project's src/ directory."""
    here = os.path.dirname(__file__)
    src = os.path.normpath(os.path.join(here, "..", "..", "src"))
    assert os.path.isdir(src), f"src dir not found: {src}"
    return src


def _pid_alive(pid: int) -> bool:
    """Return True if *pid* is a live (non-zombie) process on Linux."""
    status_path = f"/proc/{pid}/status"
    try:
        with open(status_path) as f:
            for line in f:
                if line.startswith("State:"):
                    # e.g. "State:\tZ (zombie)"
                    state = line.split()[1]
                    return state != "Z"
        return True  # file exists but no State line — treat as alive
    except FileNotFoundError:
        return False  # process gone


def _wait_all_dead(pids: list[int], timeout: float = 4.0) -> list[int]:
    """Poll until all *pids* have exited or *timeout* elapses.

    Returns the list of PIDs that are still alive after *timeout*.
    """
    deadline = time.monotonic() + timeout
    remaining = list(pids)
    while remaining and time.monotonic() < deadline:
        remaining = [pid for pid in remaining if _pid_alive(pid)]
        if remaining:
            time.sleep(0.1)
    return remaining


# ---------------------------------------------------------------------------
# The actual test
# ---------------------------------------------------------------------------


def test_workers_exit_after_parent_sigkill(tmp_path: pathlib.Path) -> None:
    """Workers die within ~4 s of their parent being SIGKILLed.

    The PR_SET_PDEATHSIG initializer installed by #479 wires the kernel's
    "deliver SIGTERM on parent exit" mechanism.  We verify it fires even
    under SIGKILL (which bypasses all Python signal handlers).
    """
    pid_file = str(tmp_path / "parent.pid")
    worker_pid_file = str(tmp_path / "workers.pid")
    src_path = _find_src_path()

    # Create the worker PID file up front so the subprocess can append to it.
    open(worker_pid_file, "w").close()

    # Launch the subprocess.
    proc = subprocess.Popen(
        [sys.executable, "-c", _WORKER_SCRIPT, pid_file, worker_pid_file, src_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Wait until at least 2 worker PIDs have been written (pool is running).
    deadline = time.monotonic() + 20.0
    worker_pids: list[int] = []
    while time.monotonic() < deadline:
        time.sleep(0.2)
        if not os.path.exists(pid_file):
            continue
        with open(worker_pid_file) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if len(lines) >= 2:
            worker_pids = [int(p) for p in lines]
            break

    with open(worker_pid_file) as _f:
        _wpf_contents = _f.read()
    assert worker_pids, (
        "Could not enumerate worker PIDs within 20 s — "
        "did the subprocess start correctly? "
        f"pid_file exists={os.path.exists(pid_file)}, "
        f"worker_pid_file contents={_wpf_contents!r}"
    )

    # Verify workers are actually alive before killing the parent.
    assert any(_pid_alive(pid) for pid in worker_pids), (
        f"All worker PIDs {worker_pids} already dead before SIGKILL — test setup issue."
    )

    # SIGKILL the parent (no Python handler can run). The parent may already
    # have exited on its own by the time we get here (CI scheduling/resource
    # pressure can widen this race well beyond what's seen locally) -- that's
    # fine for the property under test (workers get cleaned up once the
    # parent is dead), so a ProcessLookupError here is not a test failure.
    with contextlib.suppress(ProcessLookupError):
        os.kill(proc.pid, signal.SIGKILL)
    proc.wait()

    # Assert: all workers should exit within the polling window.
    surviving = _wait_all_dead(worker_pids, timeout=4.0)
    assert not surviving, (
        f"Orphaned worker PIDs still alive after parent SIGKILL: {surviving}. "
        "PR_SET_PDEATHSIG initializer may not have been installed correctly."
    )
