"""#964 -- fetch the three cycler repos and report ahead/behind before work starts.

A clean ``git status`` can hide a checkout that is far behind its remote. On
2026-10-05 the private corpus clone was 87 commits behind and had no upstream
set, so ``git status`` showed nothing behind and 71 indexed PDFs looked missing.
This script runs ``git fetch`` (which only updates remote-tracking refs; it never
touches the working tree, the index or local branches) and then reports, per
repo: branch, upstream, ahead/behind counts and the number of uncommitted
paths. A missing upstream is reported as a problem, not as "up to date".

Usage::

    python scripts/repo_sync_status.py              # the default three repos
    python scripts/repo_sync_status.py --no-fetch   # report from existing refs
    python scripts/repo_sync_status.py PATH [PATH ...]

The defaults are the siblings of this checkout: ``cyclers`` (this repo),
``cyclers.space`` and ``cyclers_pdf``. Exit status 0 when every repo exists,
has an upstream and is not behind; 1 otherwise.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPOS = (REPO_ROOT, REPO_ROOT.parent / "cyclers.space", REPO_ROOT.parent / "cyclers_pdf")


@dataclass
class RepoStatus:
    path: Path
    ok: bool
    branch: str = "-"
    upstream: str | None = None
    ahead: int = 0
    behind: int = 0
    dirty: int = 0
    problem: str | None = None

    def line(self) -> str:
        if self.problem and self.upstream is None and self.branch == "-":
            return f"PROBLEM  {self.path}: {self.problem}"
        tag = "OK     " if self.ok else "PROBLEM"
        up = self.upstream or "(no upstream)"
        extra = f"; {self.problem}" if self.problem else ""
        return (
            f"{tag}  {self.path.name}: {self.branch} -> {up}, ahead {self.ahead}, "
            f"behind {self.behind}, uncommitted {self.dirty}{extra}"
        )


def _git(path: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-C", str(path), *args], capture_output=True, text=True, check=False
    )


def repo_status(path: Path, fetch: bool = True) -> RepoStatus:
    if not (path / ".git").exists():
        return RepoStatus(path, ok=False, problem="not a git checkout (missing or not cloned)")
    if fetch:
        f = _git(path, "fetch", "--quiet")
        fetch_err = f.stderr.strip() if f.returncode != 0 else None
    else:
        fetch_err = None
    branch = _git(path, "rev-parse", "--abbrev-ref", "HEAD").stdout.strip() or "-"
    up = _git(path, "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
    upstream = up.stdout.strip() if up.returncode == 0 else None
    dirty = len([ln for ln in _git(path, "status", "--porcelain").stdout.splitlines() if ln])
    st = RepoStatus(path, ok=True, branch=branch, upstream=upstream, dirty=dirty)
    if upstream is None:
        st.ok = False
        st.problem = "no upstream set: ahead/behind unknown (git branch -u origin/<branch>)"
        return st
    counts = _git(path, "rev-list", "--left-right", "--count", "HEAD...@{u}").stdout.split()
    if len(counts) == 2:
        st.ahead, st.behind = int(counts[0]), int(counts[1])
    if st.behind > 0:
        st.ok = False
        st.problem = f"{st.behind} commit(s) behind {upstream}: pull before building on it"
    if fetch_err:
        st.ok = False
        st.problem = (st.problem + "; " if st.problem else "") + f"fetch failed: {fetch_err}"
    return st


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("paths", nargs="*", type=Path, help="repos to check (default: the three)")
    ap.add_argument("--no-fetch", action="store_true", help="do not run git fetch")
    args = ap.parse_args(argv)
    paths = args.paths or list(DEFAULT_REPOS)
    statuses = [repo_status(p, fetch=not args.no_fetch) for p in paths]
    for s in statuses:
        print(s.line())
    return 0 if all(s.ok for s in statuses) else 1


if __name__ == "__main__":
    sys.exit(main())
