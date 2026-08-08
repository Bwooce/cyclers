"""Generic checkpointed cell-grid campaign runner (#788).

Combines two pieces of infrastructure that existed separately before this task:

* :mod:`cyclerfinder.search.discovery_campaign` — the checkpointed-restart
  pattern (append-only checkpoint, ``index % n_workers`` sharding across
  independently-launched OS processes, method-versioned empty-region
  writeback). That module hosts the ``SearchTarget`` protocol tied to one
  concrete genome (repeated-moon multi-rev cyclers over one primary body).
* :mod:`cyclerfinder.parallel.parallel_sweep` — the joblib parallel substrate
  (measured 5x speedup on 8 workers), previously wired into one-shot sweep
  scripts with no supervising resumable loop around it.

Nothing before this task combined the two: every long-running discovery
script in this project was "launch N processes by hand, re-run until the
checkpoint file says done." This module is the missing supervising loop:
**one process**, parallel across cores via ``parallel_sweep``, checkpointed
durably so it can be **killed at any point (including SIGKILL) and resumed**
without corrupting state or re-doing completed work.

Design choices worth being explicit about
------------------------------------------

**Generic over "cell".** :func:`run_grid_campaign` takes an arbitrary
``Sequence[Any]`` of picklable cells (dicts, dataclasses, namedtuples — same
pickle-safety contract as :func:`~cyclerfinder.parallel.parallel_sweep.
parallel_sweep`) and a top-level ``worker`` callable that maps one cell to a
:class:`CellOutcome`. The future campaigns (#789-#792) each supply their own
grid + worker; this module has zero knowledge of resonance ratios, moon
sequences, or Tisserand pruning.

**Durability granularity is per-BATCH, not literally per-cell.**
``parallel_sweep`` only surfaces results after its whole ``Parallel()`` call
returns (its own docstring: ``on_cell_complete`` "invoked AFTER joblib
returns (not mid-flight)"). True per-cell streaming would require a custom
joblib backend, which is out of scope. Instead this runner submits
``checkpoint_batch_size`` cells per ``parallel_sweep`` call and writes every
cell in that batch durably (flushed + fsynced) the instant the batch returns,
before starting the next one. A kill loses at most one in-flight batch, never
completed batches — and ``checkpoint_batch_size`` is the caller's knob to
bound that loss window (small for expensive cells, larger for cheap ones).

**``results.jsonl`` is its own checkpoint (self-describing, single file).**
``discovery_campaign.py`` uses two files (a plain checkpoint of indices, plus
a separate review-queue / empty-region sink) and accepts a narrow race:
if the process dies between writing the sink entry and appending the
checkpoint index, that candidate re-runs on resume (harmless there only
because the domain-level dedup hash makes a re-run idempotent). This module
does not have a domain-level dedup key to lean on, so it collapses "the
result" and "the checkpoint" into ONE file: an index is "done" iff a
complete, parseable JSON line for it exists in ``results.jsonl``. This
removes the two-file race entirely. ``discovery_campaign.py`` is left
UNCHANGED — its checkpoint format genuinely differs (plain index-per-line vs.
this module's self-describing result-per-line) and it is tested, stable, and
may have a sibling task committing to it concurrently; there is no shared
logic here worth extracting beyond what already lives in
``data.empty_regions`` / ``data.method_capability``, which this module reuses
directly rather than reimplementing.

**A durable line means "this cell was actually evaluated."** ``parallel_sweep``
has two distinct failure modes that must NOT be confused:

1. A single cell's ``worker`` call raised — ``parallel_sweep`` reports that
   cell's slot as ``None`` and the rest of the batch's results are real.
2. The whole batch submission failed at the joblib/executor level (a
   ``PicklingError`` at submission time, or a worker process dying opaquely —
   OOM-killed, segfault — which loky surfaces as a whole-``Parallel()``
   exception). ``parallel_sweep`` then reports **every** cell in that call as
   ``None`` with ``notes`` prefixed ``"joblib-level failure"`` — those cells
   were NEVER evaluated.

Treating case 2 the same as case 1 would durably (and permanently, since
resume skips anything with a results line) record cells that never actually
ran as ``"error"`` — silent data loss on a weeks-long campaign. This module
detects the ``"joblib-level failure"`` prefix and refuses to write ANY result
line for that batch; it raises instead, so the batch is retried in full on
the next invocation (see :func:`run_grid_campaign`).

**Worker contract.** A worker should catch its own exceptions and return
``CellOutcome(status="error", error=<message>)`` rather than let them
propagate — ``parallel_sweep`` only preserves the first three per-cell error
messages project-wide (its ``notes`` field), so exceptions escaping the
worker lose their detail past the third failure in a batch. The ``None``
path (worker call raised and was swallowed by ``parallel_sweep``) is treated
as a durable ``"error"`` outcome by this module, but callers should not rely
on that path for diagnostics.

**Empty-region emission is an aggregate, not per-cell.** A "miss" cell is
recorded as such in ``results.jsonl``; the method-versioned
``EmptyRegionReport`` (``data/empty_regions.jsonl`` schema, reused verbatim
from :mod:`cyclerfinder.data.empty_regions`) is written ONCE, only when the
entire cell grid has been evaluated (across however many resumed invocations
that took), every cell converged to "miss", and zero cells errored. This
matches ``discovery_campaign.py``'s own convention (one empty-region record
per swept region, not one per candidate).

**Widening a completed grid needs a new ``region_id`` (or new capability
tags).** If a campaign's ``EmptyRegionSpec.method_capability`` does not
change, :func:`~cyclerfinder.data.method_capability.should_sweep` will no-op
a re-invocation over the same ``region_id`` once an empty-region record for
it exists — this is the intended re-sweep gate (a weaker-or-equal method
re-running a region a capability-comparable method already emptied learns
nothing new), not a bug. Widen the grid under a new ``region_id`` if you mean
to search adjacent territory.

**Set ``timeout_seconds_per_cell`` for real campaigns.** With no timeout, one
hung cell stalls its entire batch indefinitely and nothing in that batch
checkpoints until it returns (or is killed externally).
"""

from __future__ import annotations

import json
import os
import subprocess
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from cyclerfinder.data.empty_regions import (
    EmptyRegionReport,
    append_empty_region,
    load_empty_regions_list,
)
from cyclerfinder.data.method_capability import MethodCapability, should_sweep
from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep

CellStatus = Literal["hit", "miss", "error"]
"""A cell's outcome, in this project's own vocabulary:

* ``"hit"``  — the worker found something the campaign cares about (a
  candidate to route onward — the generic runner does not know the domain's
  routing schema, so it is left in ``CellOutcome.payload`` for the caller to
  act on after the run, e.g. append to ``data/review_queue.jsonl``).
* ``"miss"`` — the cell was evaluated cleanly and found nothing. If EVERY
  cell in the grid misses, the aggregate is a genuine negative-results-registry
  entry (``empty_regions.jsonl``).
* ``"error"`` — the worker raised (or a joblib-level failure was detected and
  the batch retried — see the module docstring). An errored cell blocks
  empty-region emission until it is re-evaluated cleanly (a kill-signal is
  not evidence of "empty").
"""


@dataclass(frozen=True)
class CellOutcome:
    """One cell's durable result. Must be pickle-safe (crosses a process boundary)."""

    status: CellStatus
    payload: dict[str, Any] = field(default_factory=dict)
    error: str = ""


@dataclass(frozen=True)
class CampaignRunnerConfig:
    """Execution knobs for :func:`run_grid_campaign`.

    Attributes
    ----------
    n_workers:
        Forwarded to :class:`~cyclerfinder.parallel.parallel_sweep.
        ParallelSweepConfig`. ``-1`` (default) = all cores.
    backend:
        Forwarded to ``ParallelSweepConfig`` (``"loky"`` default).
    checkpoint_batch_size:
        Number of cells submitted per ``parallel_sweep`` call. This is the
        durability granularity: a kill loses at most this many in-flight
        cells, never a completed batch. Deliberately distinct from
        ``ParallelSweepConfig.chunk_size`` (joblib's own dispatch-batching
        knob) to avoid confusing the two.
    timeout_seconds_per_cell:
        Forwarded to ``ParallelSweepConfig``. ``None`` = no timeout (see the
        module docstring's warning about stalled batches).
    max_batches:
        Cap the number of batches run by THIS invocation (``None`` = exhaust
        the pending cells). Lets a caller deliberately run in bounded
        segments even when a single invocation could finish the whole grid.
    pause_seconds_per_batch:
        Sleep this long between batches (0.0 = no pause, the default). A
        sensor-independent duty-cycle knob: pinning ``n_workers`` below
        ``os.cpu_count()`` already leaves headroom, but a long unattended
        campaign on shared hardware (e.g. a laptop also running other work)
        benefits from a deliberate breather too. Skipped after the final
        batch (no point delaying an already-finished run).
    thermal_backoff_seconds:
        If > 0, poll :func:`_os_thermal_throttled` between batches and sleep
        this long (in addition to ``pause_seconds_per_batch``) whenever
        macOS itself reports it is already thermally limiting CPU speed.
        0.0 (default) disables the check entirely — safe on Linux/other
        platforms where ``pmset`` doesn't exist, since the check is opt-in.
    """

    n_workers: int = -1
    backend: Literal["loky", "threading", "multiprocessing"] = "loky"
    checkpoint_batch_size: int = 25
    timeout_seconds_per_cell: float | None = None
    max_batches: int | None = None
    pause_seconds_per_batch: float = 0.0
    thermal_backoff_seconds: float = 0.0


@dataclass(frozen=True)
class CampaignRunnerRouting:
    """Where a run's artefacts live.

    ``results_path`` is BOTH the durable per-cell result log and the
    checkpoint (see the module docstring); it should live under a gitignored
    ``out/`` tree for a real campaign, matching ``discovery_campaign_daemon.py``'s
    own convention, or ``tmp_path`` in tests. ``empty_regions_path`` defaults
    to the real registry but tests must override it to a temp path — this
    module never hardcodes ``data/empty_regions.jsonl``.
    """

    results_path: Path
    empty_regions_path: Path | None = None

    def __post_init__(self) -> None:
        # Defensive coercion: CLI-driven callers (argparse) hand these in as
        # plain strings; normalise to Path once here so every downstream use
        # site (file I/O, ``.parent.mkdir``) can assume a real Path.
        object.__setattr__(self, "results_path", Path(self.results_path))
        if self.empty_regions_path is not None:
            object.__setattr__(self, "empty_regions_path", Path(self.empty_regions_path))


@dataclass(frozen=True)
class EmptyRegionSpec:
    """Everything :func:`run_grid_campaign` needs to build an ``EmptyRegionReport``.

    Mirrors ``discovery_campaign._empty_region_report`` but generic: the
    caller supplies the domain description, this module supplies the
    ``search_extent`` / ``result`` counts from the actual run.
    """

    region_id: str
    family: str
    centre: str
    method_capability: MethodCapability
    prune_gates: tuple[str, ...]
    topologies: tuple[dict[str, Any], ...] = ()
    source_anchors: str = ""
    interpretation: str = ""
    run_extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class GridCampaignStats:
    """Real, reportable counters — always derived from ``results.jsonl`` itself.

    Recomputed at the end of every invocation by re-scanning the results
    file, so a stats object is correct regardless of how many prior
    invocations contributed to it (resumed campaigns accumulate truthfully).
    """

    total_cells: int = 0
    evaluated_total: int = 0
    hits: int = 0
    misses: int = 0
    errors: int = 0
    hit_indices: tuple[int, ...] = ()
    batches_run_this_invocation: int = 0
    skipped_capability_subsumed: bool = False
    empty_region_emitted: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "total_cells": self.total_cells,
            "evaluated_total": self.evaluated_total,
            "hits": self.hits,
            "misses": self.misses,
            "errors": self.errors,
            "hit_indices": list(self.hit_indices),
            "batches_run_this_invocation": self.batches_run_this_invocation,
            "skipped_capability_subsumed": self.skipped_capability_subsumed,
            "empty_region_emitted": self.empty_region_emitted,
        }


# ---------------------------------------------------------------------------
# results.jsonl: the self-describing checkpoint
# ---------------------------------------------------------------------------


def _parse_result_line(line: str) -> dict[str, Any] | None:
    """Parse one ``results.jsonl`` line; ``None`` on any malformed line.

    A malformed line can only arise from a write interrupted mid-flight
    (SIGKILL between the OS accepting part of a ``write()`` and the rest —
    see the kill/resume test for the empirical check that this does not
    actually happen with the fsync'd single-``write()`` append below, but
    tolerating it costs nothing and turns a hypothetical corruption into "this
    cell is just not done yet, re-run it" rather than a crash).
    """
    try:
        payload = json.loads(line)
    except json.JSONDecodeError:
        return None
    if not isinstance(payload, dict) or "index" not in payload or "status" not in payload:
        return None
    return payload


def _load_results(path: Path) -> dict[int, dict[str, Any]]:
    """Every completed cell's payload, keyed by index (last line wins on dup)."""
    if not path.exists():
        return {}
    out: dict[int, dict[str, Any]] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        parsed = _parse_result_line(line)
        if parsed is None:
            continue
        out[int(parsed["index"])] = parsed
    return out


def _append_cell_result(path: Path, index: int, outcome: CellOutcome) -> None:
    """Durably append one cell's result: single write + flush + fsync.

    Writing before the next batch starts (and never batching multiple cells'
    lines into one buffered flush) is what bounds a kill's damage to the
    current in-flight batch rather than the whole run.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "index": index,
        "status": outcome.status,
        "payload": outcome.payload,
        "error": outcome.error,
    }
    line = json.dumps(record, sort_keys=True, ensure_ascii=True) + "\n"
    with path.open("a", encoding="utf-8") as fh:
        fh.write(line)
        fh.flush()
        os.fsync(fh.fileno())


# ---------------------------------------------------------------------------
# The runner
# ---------------------------------------------------------------------------


def _sleep(seconds: float) -> None:
    """Indirection so tests can stub the campaign's own pauses without
    patching the process-wide ``time.sleep`` (which loky's internal polling
    also calls -- patching that globally makes the executor spin)."""
    time.sleep(seconds)


def _os_thermal_throttled() -> bool:
    """Best-effort check: is macOS already limiting CPU speed for heat?

    Parses ``pmset -g therm``'s ``CPU_Speed_Limit``/``CPU_Scheduler_Limit``
    fields (percentages; <100 means the OS itself is actively throttling).
    Fails open (returns ``False``) on any error — missing binary, non-macOS,
    unparseable output, or no warning ever recorded (the common idle case,
    printed as prose rather than a field) — so this can never hang a
    campaign or misbehave on a platform without ``pmset``.
    """
    try:
        out = subprocess.run(
            ["pmset", "-g", "therm"],
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout
    except (OSError, subprocess.TimeoutExpired):
        return False
    for line in out.splitlines():
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        if key.strip() in ("CPU_Speed_Limit", "CPU_Scheduler_Limit"):
            try:
                if float(value.strip()) < 100.0:
                    return True
            except ValueError:
                continue
    return False


def run_grid_campaign(
    cells: Sequence[Any],
    worker: Callable[[Any], CellOutcome],
    *,
    routing: CampaignRunnerRouting,
    config: CampaignRunnerConfig | None = None,
    empty_region_spec: EmptyRegionSpec | None = None,
) -> GridCampaignStats:
    """Run (or resume) a generic checkpointed, parallel cell-grid campaign.

    ``cells[i]`` is candidate index ``i`` (the resumability key — deterministic
    ordering is the caller's responsibility, same contract as
    ``discovery_campaign.SearchTarget.enumerate_candidates``). Already-done
    indices (a complete line in ``routing.results_path``) are skipped. Pending
    indices are submitted to :func:`~cyclerfinder.parallel.parallel_sweep.
    parallel_sweep` in batches of ``config.checkpoint_batch_size``; every
    cell in a batch is durably appended the instant that batch returns.

    If ``empty_region_spec`` is given and ``routing.empty_regions_path`` is
    set, the capability-subsumption gate
    (:func:`~cyclerfinder.data.method_capability.should_sweep`) is checked
    FIRST: if a prior empty-region record for the same ``region_id`` was
    produced by a method that subsumes this one, the whole invocation is a
    no-op (``skipped_capability_subsumed=True``). Otherwise, once every cell
    in ``cells`` has a durable "miss" result (across however many invocations
    that took) with zero errors and zero hits, one ``EmptyRegionReport`` is
    appended (idempotent: skipped if a report with this ``region_id`` already
    exists in the file).
    """
    cfg = config if config is not None else CampaignRunnerConfig()
    n_total = len(cells)

    if empty_region_spec is not None and routing.empty_regions_path is not None:
        registry = load_empty_regions_list(routing.empty_regions_path)
        if not should_sweep(
            region_id=empty_region_spec.region_id,
            method=empty_region_spec.method_capability,
            registry=registry,
        ):
            done_now = _load_results(routing.results_path)
            return GridCampaignStats(
                total_cells=n_total,
                evaluated_total=len(done_now),
                skipped_capability_subsumed=True,
            )

    done = _load_results(routing.results_path)
    pending = [i for i in range(n_total) if i not in done]

    psweep_cfg = ParallelSweepConfig(
        n_workers=cfg.n_workers,
        backend=cfg.backend,
        chunk_size=1,
        timeout_seconds_per_cell=cfg.timeout_seconds_per_cell,
    )

    batches_run = 0
    for start in range(0, len(pending), cfg.checkpoint_batch_size):
        if cfg.max_batches is not None and batches_run >= cfg.max_batches:
            break
        batch_indices = pending[start : start + cfg.checkpoint_batch_size]
        batch_cells = [cells[i] for i in batch_indices]

        result = parallel_sweep(batch_cells, worker, config=psweep_cfg)

        if result.notes.startswith("joblib-level failure"):
            # Whole-batch submission/executor failure (PicklingError, a worker
            # OOM-killed/segfaulted, ...): NONE of these cells were actually
            # evaluated. Writing durable "error" lines here would permanently
            # (resume treats any durable line as done) hide never-run cells.
            # Refuse to write anything for this batch and surface the failure
            # loudly; the batch is retried in full on the next invocation.
            raise RuntimeError(
                f"run_grid_campaign: batch starting at pending-index {start} "
                f"({len(batch_indices)} cells) failed at the joblib/executor "
                f"level, not per-cell: {result.notes}. No result lines were "
                f"written for this batch; re-invoke to retry it."
            )

        for local_i, idx in enumerate(batch_indices):
            outcome = result.results[local_i]
            if not isinstance(outcome, CellOutcome):
                # The worker call raised and parallel_sweep swallowed it (or a
                # contract-violating worker returned something else). Either
                # way the cell WAS submitted and DID return control to the
                # parent, so recording it as a durable error (rather than
                # silently retrying forever) is correct here — unlike the
                # joblib-level case above, this is isolated to one cell.
                outcome = CellOutcome(
                    status="error",
                    error=f"worker did not return a CellOutcome (batch notes: {result.notes})",
                )
            _append_cell_result(routing.results_path, idx, outcome)
        batches_run += 1

        more_pending = start + cfg.checkpoint_batch_size < len(pending)
        more_allowed = cfg.max_batches is None or batches_run < cfg.max_batches
        if more_pending and more_allowed:
            if cfg.pause_seconds_per_batch > 0:
                _sleep(cfg.pause_seconds_per_batch)
            if cfg.thermal_backoff_seconds > 0 and _os_thermal_throttled():
                _sleep(cfg.thermal_backoff_seconds)

    final = _load_results(routing.results_path)
    hits = sum(1 for r in final.values() if r["status"] == "hit")
    misses = sum(1 for r in final.values() if r["status"] == "miss")
    errors = sum(1 for r in final.values() if r["status"] == "error")
    hit_indices = tuple(sorted(i for i, r in final.items() if r["status"] == "hit"))

    stats = GridCampaignStats(
        total_cells=n_total,
        evaluated_total=len(final),
        hits=hits,
        misses=misses,
        errors=errors,
        hit_indices=hit_indices,
        batches_run_this_invocation=batches_run,
    )

    if (
        empty_region_spec is not None
        and routing.empty_regions_path is not None
        and stats.evaluated_total == n_total
        and stats.hits == 0
        and stats.errors == 0
        and n_total > 0
    ):
        stats.empty_region_emitted = _maybe_emit_empty_region(
            routing.empty_regions_path, empty_region_spec, stats
        )

    return stats


def _maybe_emit_empty_region(
    path: Path,
    spec: EmptyRegionSpec,
    stats: GridCampaignStats,
) -> bool:
    """Append the aggregate empty-region report, unless already present.

    Idempotent on ``region_id`` so re-invoking a fully-complete, all-miss
    campaign never double-writes the registry.
    """
    existing_ids = {r.region_id for r in load_empty_regions_list(path)}
    if spec.region_id in existing_ids:
        return False
    report = EmptyRegionReport(
        region_id=spec.region_id,
        family=spec.family,
        centre=spec.centre,
        topologies=spec.topologies,
        method_capability=spec.method_capability,
        search_extent={
            "points_total": stats.total_cells,
            "candidates_evaluated": stats.evaluated_total,
        },
        prune_gates=spec.prune_gates,
        result={
            "hits": stats.hits,
            "misses": stats.misses,
            "errors": stats.errors,
        },
        verdict="empty",
        interpretation=spec.interpretation
        or (
            "No hit cell found across the full swept grid; empty as far as this method could reach."
        ),
        source_anchors=spec.source_anchors,
        run={"region_id": spec.region_id, **spec.run_extra},
    )
    append_empty_region(path, report)
    return True


__all__ = [
    "CampaignRunnerConfig",
    "CampaignRunnerRouting",
    "CellOutcome",
    "CellStatus",
    "EmptyRegionSpec",
    "GridCampaignStats",
    "run_grid_campaign",
]
