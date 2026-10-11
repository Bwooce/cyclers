"""Run a child process with a timeout that kills its WHOLE process group (#1068).

``subprocess.run(cmd, timeout=...)`` kills only the direct child on a timeout. When that child is a
launcher such as ``uv run python ...``, the python grandchild survives as an orphan and can write
its result after the caller has recorded the run as NOT RUN
(``docs/papercuts/2026-10-11-subprocess-timeout-leaves-uv-orphan.md``). :func:`run_group` starts the
child in its own session and, on a timeout, sends SIGKILL to the whole group before re-raising
:class:`subprocess.TimeoutExpired`, the same contract as ``subprocess.run``.
"""

from __future__ import annotations

import os
import signal
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import IO, Any


def run_group(
    cmd: Sequence[str],
    *,
    timeout: float,
    cwd: str | Path | None = None,
    capture_output: bool = False,
    text: bool = False,
    check: bool = False,
    stdout: int | IO[Any] | None = None,
    stderr: int | IO[Any] | None = None,
) -> subprocess.CompletedProcess[Any]:
    """``subprocess.run`` with ``timeout`` that kills the child's whole process group.

    Raises :class:`subprocess.TimeoutExpired` after the group is dead, and
    :class:`subprocess.CalledProcessError` on a non-zero exit when ``check`` is set.
    """
    if capture_output:
        stdout = stderr = subprocess.PIPE
    p = subprocess.Popen(
        list(cmd), cwd=cwd, stdout=stdout, stderr=stderr, text=text, start_new_session=True
    )
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        raise
    done = subprocess.CompletedProcess(list(cmd), p.returncode, out, err)
    if check:
        done.check_returncode()
    return done
