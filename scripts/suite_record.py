#!/usr/bin/env python3
"""#1069: record a green window suite so the ratchet hooks can trust it, not re-run pytest.

Owner decision 2026-10-12: "Hook trusts a fresh suite log". Design and expected test results:
docs/notes/2026-10-12-1069-hook-trusts-fresh-suite-log.md.

    scripts/suite_record.py start
    scripts/suite_record.py finish --log c1.log:0 --log c2.log:0 --covers tests/data ...
    scripts/suite_record.py check <pytest target> ...   (run by scripts/ratchet_hook.sh)

The record lives in <git common dir>/CYCLERS_SUITE_OK (inside .git, never tracked). `check`
exits 0 only when the record is ok, under an hour old, its content hash equals that of the tree
being committed, and the hook's targets are inside the record's covers. Stdlib only and
Python 3.9 compatible, because the hook runs it with the system python3.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

MAX_AGE_S = 3600.0
COVERED_DIRS = ("src/", "tests/", "data/", "scripts/")
COVERED_FILES = ("pyproject.toml", "uv.lock", "conftest.py", ".pre-commit-config.yaml")
EXCLUDED_SUFFIXES = (".md", ".txt", ".log")
FAILED_RE = re.compile(r"^(FAILED|ERROR) ")
COMPLETE_RE = re.compile(r"\[\s*100%\]")


def _git(*args: str, env: dict[str, str] | None = None) -> str:
    return subprocess.run(
        ["git", *args], capture_output=True, text=True, check=True, env=env
    ).stdout.strip()


def _common_dir() -> Path:
    return Path(_git("rev-parse", "--git-common-dir")).resolve()


def record_path() -> Path:
    return _common_dir() / "CYCLERS_SUITE_OK"


def pending_path() -> Path:
    return _common_dir() / "CYCLERS_SUITE_PENDING"


def covered(path: str) -> bool:
    if path.endswith(EXCLUDED_SUFFIXES):
        return False
    return path.startswith(COVERED_DIRS) or path in COVERED_FILES


def tree_content_hash(tree: str) -> str:
    """SHA-256 over the sorted (mode, blob id, path) of the covered files in a git tree."""
    out = subprocess.run(
        ["git", "ls-tree", "-r", "-z", tree], capture_output=True, text=True, check=True
    ).stdout
    lines = []
    for entry in out.split("\0"):
        if not entry:
            continue
        meta, path = entry.split("\t", 1)
        mode, _kind, blob = meta.split()
        if covered(path):
            lines.append(f"{mode} {blob} {path}")
    return hashlib.sha256("\n".join(sorted(lines)).encode()).hexdigest()


def worktree_hash() -> str:
    """Content hash of the working tree's tracked files (what pytest read), via a temp index."""
    real_index = Path(_git("rev-parse", "--git-path", "index"))
    toplevel = Path(_git("rev-parse", "--show-toplevel"))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_index = Path(tmp) / "index"
        if real_index.exists():
            shutil.copy2(real_index, tmp_index)
        env = dict(os.environ, GIT_INDEX_FILE=str(tmp_index))
        paths = [p.rstrip("/") for p in COVERED_DIRS + COVERED_FILES]
        present = [p for p in paths if (toplevel / p).exists()]
        if present:
            subprocess.run(
                ["git", "-C", str(toplevel), "add", "-u", "--", *present],
                capture_output=True,
                text=True,
                check=True,
                env=env,
            )
        tree = _git("-C", str(toplevel), "write-tree", env=env)
    return tree_content_hash(tree)


def index_hash() -> str:
    """Content hash of the tree git is about to commit (the index the hook sees)."""
    return tree_content_hash(_git("write-tree"))


def _now() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc)  # noqa: UP017 (the hook runs system Python 3.9)


def _local(ts: str) -> str:
    when = dt.datetime.fromisoformat(ts)
    try:
        from zoneinfo import ZoneInfo

        return when.astimezone(ZoneInfo("Australia/Sydney")).strftime("%Y-%m-%d %H:%M %Z")
    except Exception:
        return ts


def cmd_start(_args: argparse.Namespace) -> int:
    data = {"started": _now().isoformat(), "hash": worktree_hash()}
    pending_path().write_text(json.dumps(data, indent=1) + "\n")
    print(f"suite record: started {_local(data['started'])}, hash {data['hash'][:12]}")
    return 0


def _parse_log_arg(spec: str) -> tuple[Path, int]:
    path, sep, code = spec.rpartition(":")
    if not sep or not path:
        raise SystemExit(f"--log wants <file>:<exit status>, got {spec!r}")
    return Path(path), int(code)


def cmd_finish(args: argparse.Namespace) -> int:
    pending = pending_path()
    if not pending.exists():
        print(
            "suite record: no 'start' for this suite; run 'scripts/suite_record.py start' first",
            file=sys.stderr,
        )
        return 2
    start = json.loads(pending.read_text())
    end_hash = worktree_hash()
    reasons: list[str] = []
    logs = []
    failed: list[str] = []
    for spec in args.log:
        path, code = _parse_log_arg(spec)
        text = path.read_text(errors="replace") if path.exists() else ""
        complete = bool(COMPLETE_RE.search(text))
        bad = [ln.rstrip() for ln in text.splitlines() if FAILED_RE.match(ln)]
        failed.extend(bad)
        logs.append({"path": str(path), "exit": code, "complete": complete, "failed": len(bad)})
        if not path.exists():
            reasons.append(f"{path}: log missing")
        if code != 0:
            reasons.append(f"{path}: exit status {code}")
        if not complete:
            reasons.append(f"{path}: no [100%] progress line (suite did not finish)")
    if failed:
        reasons.append(f"{len(failed)} FAILED/ERROR lines")
    if end_hash != start["hash"]:
        reasons.append("covered files changed during the suite")
    if not args.log:
        reasons.append("no logs given")
    if not args.covers:
        reasons.append("no covers given")
    rec = {
        "started": start["started"],
        "finished": _now().isoformat(),
        "hash": start["hash"],
        "hash_at_finish": end_hash,
        "logs": logs,
        "failed": failed,
        "covers": [c.rstrip("/") for c in args.covers],
        "ok": not reasons,
        "reasons": reasons,
    }
    record_path().write_text(json.dumps(rec, indent=1) + "\n")
    pending.unlink()
    if rec["ok"]:
        print(f"suite record: OK, started {_local(rec['started'])}, hash {rec['hash'][:12]}")
        return 0
    print("suite record: NOT OK (the hooks will run pytest): " + "; ".join(reasons))
    return 1


def _target_path(target: str) -> str:
    return target.split("::", 1)[0].rstrip("/")


def _inside(path: str, cover: str) -> bool:
    return path == cover or path.startswith(cover + "/")


def check_reason(targets: list[str]) -> str | None:
    """None when the record can be trusted for these hook targets, else the reason it cannot."""
    path = record_path()
    if not path.exists():
        return "no window-suite record"
    try:
        rec = json.loads(path.read_text())
    except (OSError, ValueError):
        return "record unreadable"
    if not rec.get("ok"):
        return "record not ok: " + "; ".join(rec.get("reasons", []))
    try:
        age = (_now() - dt.datetime.fromisoformat(rec["started"])).total_seconds()
    except (KeyError, ValueError):
        return "record has no valid start time"
    if not 0 <= age < MAX_AGE_S:
        return f"record is stale ({age / 60:.0f} min since the suite started)"
    if rec.get("hash") != index_hash():
        return "tree being committed differs from the tree the suite ran on"
    covers = rec.get("covers", [])
    outside = [t for t in targets if not any(_inside(_target_path(t), c) for c in covers)]
    if outside:
        return "suite did not cover " + ", ".join(outside)
    return None


def cmd_check(args: argparse.Namespace) -> int:
    targets = [t for t in args.targets if not t.startswith("-")]
    reason = check_reason(targets)
    if reason is None:
        started = json.loads(record_path().read_text())["started"]
        print(f"ratchets: trusted window-suite log {_local(started)}")
        return 0
    print(f"ratchets: running pytest ({reason})")
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("start", help="before the first suite chunk").set_defaults(func=cmd_start)
    fin = sub.add_parser("finish", help="after the last suite chunk")
    fin.add_argument("--log", action="append", default=[], help="<tee'd log>:<exit status>")
    fin.add_argument("--covers", nargs="+", default=[], help="pytest targets the suite ran")
    fin.set_defaults(func=cmd_finish)
    chk = sub.add_parser("check", help="used by scripts/ratchet_hook.sh")
    chk.add_argument("targets", nargs=argparse.REMAINDER)
    chk.set_defaults(func=cmd_check)
    args = parser.parse_args(argv)
    result: int = args.func(args)
    return result


if __name__ == "__main__":
    sys.exit(main())
