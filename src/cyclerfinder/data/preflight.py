"""Mandatory pre-flight gate for one-off search scripts (#521 phase 2).

Phase 1 (2026-07-02) proved the negative registry CAN mechanically catch a
redundant re-run (see ``tests/data/test_method_capability.py``). This module
makes that check ACTUALLY BLOCK one: every ``scripts/run_*.py`` is required
(by ``tests/scripts/test_scripts_call_preflight.py``, an AST-based ratchet)
to call :func:`preflight_search` near the top of ``main()``, before building
its grid or spawning any workers.

Three checks, each grounded in a real incident from the 2026-07-01/02
session this module exists to stop from recurring:

1. **Task-number hygiene** (the #513/#516 double-claim incident): ``task_no``
   must already be recorded in ``data/OUTSTANDING.md``'s TASK ALLOCATIONS
   section, and must match the number embedded in the calling script's own
   filename (``run_NNN_*.py``) if present — the exact self-inconsistency
   #515-517 shipped with this session after a rename.
2. **Registry subsumption** (the "#411 says NOT more grids, #515-517 did it
   anyway" incident): :func:`~cyclerfinder.data.method_capability.should_sweep`
   against the real on-disk ``empty_regions.jsonl`` refuses to re-run a
   search an equal-or-stronger prior method already emptied.
3. **Timing-pilot requirement** (the #520 12-hour, zero-output abort): a grid
   above :data:`LARGE_GRID_THRESHOLD` points must supply a measured
   ``timing_pilot_seconds_per_point`` before it is allowed to run at full
   size.

``override_reason`` downgrades every failing check to a printed + logged
WARNING instead of a raise — the escape hatch is explicit and audited, never
silent. Every invocation (blocked, overridden, or clean) is appended to the
pre-flight runlog (default ``data/runlogs/preflight_runlog.jsonl``,
gitignored) for the audit trail.
"""

from __future__ import annotations

import datetime
import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    load_empty_regions_list,
)
from cyclerfinder.data.method_capability import MethodCapability, should_sweep

_REPO_ROOT: Path = Path(__file__).resolve().parent.parent.parent.parent
DEFAULT_OUTSTANDING_PATH: Path = _REPO_ROOT / "data" / "OUTSTANDING.md"
DEFAULT_PREFLIGHT_RUNLOG_PATH: Path = _REPO_ROOT / "data" / "runlogs" / "preflight_runlog.jsonl"

# A grid at or below this many points does not need a measured timing pilot.
# The #520 incident (8,640 points, 12+ hours, zero persisted output) is the
# reason this exists — pick a threshold small enough that even a *wrong*
# estimate still fails cheap. 500 points at a generous 60s/point worst case
# is ~8 hours: already large enough to want a pilot first.
LARGE_GRID_THRESHOLD = 500

_TASK_ALLOCATION_RE = re.compile(r"^-\s+\*\*#(\d+)\*\*", re.MULTILINE)
_SCRIPT_FILENAME_TASK_RE = re.compile(r"run_(\d+)_")


class PreflightBlockedError(RuntimeError):
    """Raised by :func:`preflight_search` when a check fails with no override."""


@dataclass(frozen=True)
class PreflightResult:
    """The gate's verdict — always returned, even when overridden."""

    proceed: bool
    task_no: int
    region_id: str
    warnings: tuple[str, ...]
    override_reason: str | None


def _task_numbers_in_outstanding(outstanding_path: Path) -> set[int]:
    if not outstanding_path.exists():
        return set()
    text = outstanding_path.read_text(encoding="utf-8")
    return {int(m) for m in _TASK_ALLOCATION_RE.findall(text)}


def _filename_declared_task_no(script_path: Path) -> int | None:
    m = _SCRIPT_FILENAME_TASK_RE.search(script_path.name)
    return int(m.group(1)) if m is not None else None


def _run_checks(
    *,
    task_no: int,
    region_id: str,
    method: MethodCapability,
    script_path: Path,
    n_points: int,
    timing_pilot_seconds_per_point: float | None,
    outstanding_path: Path,
    registry: Sequence[EmptyRegionReport],
) -> list[str]:
    """Return a list of failure messages (empty = every check passed)."""
    failures: list[str] = []

    allocated = _task_numbers_in_outstanding(outstanding_path)
    if task_no not in allocated:
        failures.append(
            f"task #{task_no} is not recorded in {outstanding_path}'s TASK "
            f"ALLOCATIONS section (per [[project_task_numbering_convention]]). "
            f"Register it there before running a search that claims it — an "
            f"unregistered number risks colliding with work a concurrent agent "
            f"has already claimed (this happened twice in 24h on 2026-07-01/02)."
        )

    declared = _filename_declared_task_no(script_path)
    if declared is not None and declared != task_no:
        failures.append(
            f"script {script_path.name} is named for #{declared} but calls "
            f"preflight_search(task_no={task_no}) — filename and declared task "
            f"number must match (exactly the #515-517 rename bug from this "
            f"session, where a filename got renumbered but the internal #NNN "
            f"references did not)."
        )

    if not should_sweep(region_id=region_id, method=method, registry=registry):
        failures.append(
            f"region {region_id!r} is already covered by an equal-or-stronger "
            f"prior negative in the registry (should_sweep() returned False). "
            f"Re-running an equal-or-weaker method here learns nothing new — "
            f"this is the exact #411/#515-517 re-run pattern this gate exists "
            f"to stop. If your method is genuinely more capable, its "
            f"capability_tags should reflect that (see "
            f"src/cyclerfinder/data/method_capability.py's partial order)."
        )

    if n_points > LARGE_GRID_THRESHOLD and timing_pilot_seconds_per_point is None:
        failures.append(
            f"n_points={n_points} exceeds LARGE_GRID_THRESHOLD="
            f"{LARGE_GRID_THRESHOLD} with no timing_pilot_seconds_per_point "
            f"supplied. Time a small pilot (~32-64 points spanning the grid's "
            f"parameter ranges) and pass the measured per-point cost — this is "
            f"the exact #520 incident (an unbudgeted 8,640-point grid ran 16 "
            f"cores for 12+ hours and produced zero persisted output)."
        )

    return failures


def _append_runlog(
    path: Path,
    *,
    task_no: int,
    region_id: str,
    proceed: bool,
    failures: Sequence[str],
    override_reason: str | None,
    n_points: int,
    timing_pilot_seconds_per_point: float | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "task_no": task_no,
        "region_id": region_id,
        "proceed": proceed,
        "failures": list(failures),
        "override_reason": override_reason,
        "n_points": n_points,
        "timing_pilot_seconds_per_point": timing_pilot_seconds_per_point,
    }
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")


def preflight_search(
    *,
    task_no: int,
    region_id: str,
    method: MethodCapability,
    script_path: Path | str,
    n_points: int,
    timing_pilot_seconds_per_point: float | None = None,
    override_reason: str | None = None,
    outstanding_path: Path | None = None,
    empty_regions_path: Path | None = None,
    registry: Sequence[EmptyRegionReport] | None = None,
    runlog_path: Path | None = None,
) -> PreflightResult:
    """The mandatory gate every ``scripts/run_*.py`` must call before sweeping.

    Runs three checks (task-number hygiene, registry subsumption, timing-pilot
    requirement — see module docstring) and either returns a clean
    ``PreflightResult(proceed=True, warnings=())`` or raises
    :class:`PreflightBlockedError` naming every failing check. Pass
    ``override_reason`` to downgrade every failure to a logged warning instead
    of a raise (``PreflightResult.proceed`` stays ``True``,
    ``PreflightResult.warnings`` carries the messages) — the explicit, audited
    escape hatch.

    Every call is appended to the pre-flight runlog (default
    ``data/runlogs/preflight_runlog.jsonl``, gitignored) regardless of outcome.
    """
    resolved_outstanding = outstanding_path or DEFAULT_OUTSTANDING_PATH
    resolved_empty_regions = empty_regions_path or DEFAULT_EMPTY_REGIONS_PATH
    resolved_registry = (
        registry if registry is not None else load_empty_regions_list(resolved_empty_regions)
    )
    resolved_runlog = runlog_path or DEFAULT_PREFLIGHT_RUNLOG_PATH
    resolved_script_path = Path(script_path)

    failures = _run_checks(
        task_no=task_no,
        region_id=region_id,
        method=method,
        script_path=resolved_script_path,
        n_points=n_points,
        timing_pilot_seconds_per_point=timing_pilot_seconds_per_point,
        outstanding_path=resolved_outstanding,
        registry=resolved_registry,
    )

    proceed = not failures or override_reason is not None
    _append_runlog(
        resolved_runlog,
        task_no=task_no,
        region_id=region_id,
        proceed=proceed,
        failures=failures,
        override_reason=override_reason,
        n_points=n_points,
        timing_pilot_seconds_per_point=timing_pilot_seconds_per_point,
    )

    if failures and override_reason is None:
        raise PreflightBlockedError(
            f"preflight_search blocked task #{task_no} region {region_id!r} "
            f"({len(failures)} check(s) failed):\n"
            + "\n".join(f"  - {f}" for f in failures)
            + "\n\nIf you are confident this is a genuine exception, re-call "
            "with override_reason=<why> to proceed anyway (logged)."
        )

    if failures:
        for f in failures:
            print(f"[preflight WARNING, overridden: {override_reason}] {f}")

    return PreflightResult(
        proceed=proceed,
        task_no=task_no,
        region_id=region_id,
        warnings=tuple(failures),
        override_reason=override_reason,
    )


__all__ = [
    "DEFAULT_OUTSTANDING_PATH",
    "DEFAULT_PREFLIGHT_RUNLOG_PATH",
    "LARGE_GRID_THRESHOLD",
    "MethodCapability",
    "PreflightBlockedError",
    "PreflightResult",
    "preflight_search",
]
