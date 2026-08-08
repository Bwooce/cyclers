"""Tests for the generic checkpointed cell-grid campaign runner (#788).

Covers the mandatory-per-task properties:

* durable per-cell writes (results land on disk before the next batch starts),
* resume-from-checkpoint (in-process: a capped invocation + a follow-up
  invocation skips completed indices and finishes the grid exactly once),
* resume-after-a-REAL-kill (a subprocess is SIGKILLed mid-run and a second
  subprocess invocation resumes it correctly -- not mocked),
* parallel execution (distinct worker PIDs across cells in one batch),
* empty-region emission (schema round-trips through the real
  ``data.empty_regions`` validator, and is capability-subsumption-gated),
* the joblib-level (whole-batch) failure path does NOT durably mark
  never-evaluated cells as done.

All top-level worker closures live at module scope so they pickle under
loky (the project's own established contract, see
``tests/parallel/test_parallel_sweep.py``).
"""

from __future__ import annotations

import contextlib
import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from cyclerfinder.data.empty_regions import load_empty_regions_list, validate_empty_region
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.search.campaign_runner import (
    CampaignRunnerConfig,
    CampaignRunnerRouting,
    CellOutcome,
    EmptyRegionSpec,
    run_grid_campaign,
)

# ---------------------------------------------------------------------------
# Top-level worker closures (pickle-safe)
# ---------------------------------------------------------------------------


def _worker_hit_on_multiple_of_7(cell: dict[str, int]) -> CellOutcome:
    n = cell["n"]
    if n % 7 == 0 and n != 0:
        return CellOutcome(status="hit", payload={"n": n, "pid": os.getpid()})
    return CellOutcome(status="miss", payload={"n": n, "pid": os.getpid()})


def _worker_always_miss(cell: dict[str, int]) -> CellOutcome:
    return CellOutcome(status="miss", payload={"n": cell["n"], "pid": os.getpid()})


def _worker_pid_only(cell: dict[str, int]) -> CellOutcome:
    """Records the worker PID and does a small amount of real work, no sleep."""
    total = sum(i * i for i in range(2000))
    return CellOutcome(status="miss", payload={"n": cell["n"], "pid": os.getpid(), "chk": total})


def _worker_raises_on_13(cell: dict[str, int]) -> CellOutcome:
    if cell["n"] == 13:
        raise ValueError("boom at 13")
    return CellOutcome(status="miss", payload={"n": cell["n"]})


def _worker_crashes_the_process_on_2(cell: dict[str, int]) -> CellOutcome:
    """Kills its OWN worker process abruptly for cell n==2 (loky's own
    ``TerminatedWorkerError`` trigger -- a segfault/OOM-kill analog), so the
    whole batch's ``Parallel()`` call fails at the joblib/executor level, not
    just this one cell."""
    if cell["n"] == 2:
        os._exit(1)  # pragma: no cover -- terminates this worker process
    return CellOutcome(status="miss", payload={"n": cell["n"]})


# ---------------------------------------------------------------------------
# Basic durability + parallel-execution
# ---------------------------------------------------------------------------


def _method_capability(tag: str = "toy") -> MethodCapability:
    return MethodCapability(
        genome="toy-grid (#788 test)",
        corrector="none",
        capability_tags=frozenset({tag}),
        git_sha="test",
    )


def test_durable_per_cell_writes_and_hit_detection(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(20)]
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")
    config = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=5)

    stats = run_grid_campaign(cells, _worker_hit_on_multiple_of_7, routing=routing, config=config)

    assert stats.total_cells == 20
    assert stats.evaluated_total == 20
    assert stats.hits == 2  # 7 and 14
    assert stats.hit_indices == (7, 14)
    assert stats.misses == 18
    assert stats.errors == 0

    # Every cell landed on disk, one JSON line each, in a well-formed shape.
    lines = routing.results_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 20
    seen_indices = set()
    for line in lines:
        rec = json.loads(line)
        assert set(rec) == {"index", "status", "payload", "error"}
        seen_indices.add(rec["index"])
    assert seen_indices == set(range(20))


def test_parallel_execution_uses_multiple_worker_pids(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(40)]
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")
    config = CampaignRunnerConfig(n_workers=4, checkpoint_batch_size=40, backend="loky")

    run_grid_campaign(cells, _worker_pid_only, routing=routing, config=config)

    lines = routing.results_path.read_text(encoding="utf-8").splitlines()
    pids = {json.loads(line)["payload"]["pid"] for line in lines}
    n_cpus = os.cpu_count() or 1
    if n_cpus >= 2:
        assert len(pids) >= 2, f"expected work spread across >1 worker process, got pids={pids}"


# ---------------------------------------------------------------------------
# Resume-from-checkpoint (in-process, capped invocation then a follow-up)
# ---------------------------------------------------------------------------


def test_resume_skips_completed_cells_and_finishes_exactly_once(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(30)]
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")

    # First invocation: cap at 2 batches of 5 -> 10 cells done, 20 pending.
    cfg1 = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=5, max_batches=2)
    stats1 = run_grid_campaign(cells, _worker_always_miss, routing=routing, config=cfg1)
    assert stats1.evaluated_total == 10
    assert stats1.batches_run_this_invocation == 2

    lines_after_first = routing.results_path.read_text(encoding="utf-8").splitlines()
    assert len(lines_after_first) == 10
    done_after_first = {json.loads(line)["index"] for line in lines_after_first}
    assert done_after_first == set(range(10))

    # Second invocation: no cap -> finishes the remaining 20, does NOT redo the first 10.
    cfg2 = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=5)
    stats2 = run_grid_campaign(cells, _worker_always_miss, routing=routing, config=cfg2)
    assert stats2.evaluated_total == 30
    assert stats2.batches_run_this_invocation == 4  # 20 remaining / 5 per batch

    lines_final = routing.results_path.read_text(encoding="utf-8").splitlines()
    indices_final = [json.loads(line)["index"] for line in lines_final]
    assert len(indices_final) == 30, "no cell should be re-run (would create a duplicate line)"
    assert sorted(indices_final) == list(range(30))


# ---------------------------------------------------------------------------
# Empty-region emission: format + capability-subsumption gate
# ---------------------------------------------------------------------------


def test_empty_region_emitted_only_on_full_all_miss_completion(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(6)]
    empty_regions_path = tmp_path / "empty_regions.jsonl"
    routing = CampaignRunnerRouting(
        results_path=tmp_path / "results.jsonl",
        empty_regions_path=empty_regions_path,
    )
    spec = EmptyRegionSpec(
        region_id="toy-region-788",
        family="toy grid (#788 test)",
        centre="nowhere",
        method_capability=_method_capability(),
        prune_gates=("n % 7 == 0",),
    )

    # Partial run: no empty-region yet (grid not fully evaluated).
    cfg_partial = CampaignRunnerConfig(n_workers=1, checkpoint_batch_size=3, max_batches=1)
    stats_partial = run_grid_campaign(
        cells, _worker_always_miss, routing=routing, config=cfg_partial, empty_region_spec=spec
    )
    assert stats_partial.empty_region_emitted is False
    assert not empty_regions_path.exists()

    # Finish the grid: now the aggregate empty-region record is written once.
    stats_full = run_grid_campaign(
        cells, _worker_always_miss, routing=routing, empty_region_spec=spec
    )
    assert stats_full.empty_region_emitted is True
    reports = load_empty_regions_list(empty_regions_path)
    assert len(reports) == 1
    report = reports[0]
    assert report.region_id == "toy-region-788"
    assert report.search_extent["points_total"] == 6
    assert report.search_extent["candidates_evaluated"] == 6
    # The real project-wide validator must accept the emitted record.
    validate_empty_region(report)

    # Re-invoking after full completion must NOT double-write the registry
    # (capability-subsumption gate: same method, region already recorded empty).
    stats_again = run_grid_campaign(
        cells, _worker_always_miss, routing=routing, empty_region_spec=spec
    )
    assert stats_again.skipped_capability_subsumed is True
    assert len(load_empty_regions_list(empty_regions_path)) == 1


def test_a_single_hit_blocks_empty_region_emission(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(10)]
    empty_regions_path = tmp_path / "empty_regions.jsonl"
    routing = CampaignRunnerRouting(
        results_path=tmp_path / "results.jsonl",
        empty_regions_path=empty_regions_path,
    )
    spec = EmptyRegionSpec(
        region_id="toy-region-788-hit",
        family="toy grid (#788 test)",
        centre="nowhere",
        method_capability=_method_capability(),
        prune_gates=("n % 7 == 0",),
    )
    stats = run_grid_campaign(
        cells, _worker_hit_on_multiple_of_7, routing=routing, empty_region_spec=spec
    )
    assert stats.hits == 1  # n=7
    assert stats.empty_region_emitted is False
    assert not empty_regions_path.exists()


def test_capability_subsumption_gate_skips_a_subsumed_resweep(tmp_path: Path) -> None:
    """A stronger prior empty-region blocks a weaker-or-equal method's re-sweep."""
    empty_regions_path = tmp_path / "empty_regions.jsonl"
    from cyclerfinder.data.empty_regions import EmptyRegionReport, append_empty_region

    strong = MethodCapability(
        genome="toy strong",
        corrector="none",
        capability_tags=frozenset({"multi-arc"}),
        git_sha="prior",
    )
    weak = MethodCapability(
        genome="toy weak",
        corrector="none",
        capability_tags=frozenset({"single-arc"}),
        git_sha="test",
    )
    append_empty_region(
        empty_regions_path,
        EmptyRegionReport(
            region_id="toy-subsumed-region",
            family="toy",
            centre="nowhere",
            topologies=(),
            method_capability=strong,
            search_extent={"points_total": 5},
            prune_gates=("gate",),
            result={},
            verdict="empty",
            interpretation="prior sweep",
            source_anchors="",
            run={},
        ),
    )

    cells = [{"n": i} for i in range(5)]
    routing = CampaignRunnerRouting(
        results_path=tmp_path / "results.jsonl", empty_regions_path=empty_regions_path
    )
    spec = EmptyRegionSpec(
        region_id="toy-subsumed-region",
        family="toy",
        centre="nowhere",
        method_capability=weak,
        prune_gates=("gate",),
    )
    stats = run_grid_campaign(cells, _worker_always_miss, routing=routing, empty_region_spec=spec)
    assert stats.skipped_capability_subsumed is True
    assert stats.evaluated_total == 0
    assert not routing.results_path.exists()


# ---------------------------------------------------------------------------
# Per-cell error vs. whole-batch (joblib-level) failure
# ---------------------------------------------------------------------------


def test_per_cell_exception_is_recorded_as_a_durable_error(tmp_path: Path) -> None:
    cells = [{"n": i} for i in range(20)]
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")
    config = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=20)

    stats = run_grid_campaign(cells, _worker_raises_on_13, routing=routing, config=config)
    assert stats.errors == 1
    assert stats.evaluated_total == 20
    lines = {
        json.loads(line)["index"]: json.loads(line)
        for line in routing.results_path.read_text(encoding="utf-8").splitlines()
    }
    assert lines[13]["status"] == "error"
    assert all(lines[i]["status"] == "miss" for i in range(20) if i != 13)


def test_whole_batch_joblib_level_failure_writes_nothing_and_raises(tmp_path: Path) -> None:
    """A worker process dying abruptly (segfault/OOM-kill analog) must NOT
    durably mark ANY cell in that batch done.

    ``_worker_crashes_the_process_on_2`` calls ``os._exit(1)`` inside its own
    worker process for cell ``n==2``. loky detects the abrupt worker death as
    a ``TerminatedWorkerError`` and fails the WHOLE ``Parallel()`` call for
    that batch (every cell submitted with it, not just cell 2).
    ``parallel_sweep`` reports this as a joblib-level failure; the runner
    must refuse to write result lines for that batch (they were never
    actually evaluated) and raise loudly instead of silently marking them
    done.
    """
    cells = [{"n": i} for i in range(4)]
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")
    config = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=4, backend="loky")

    with pytest.raises(RuntimeError, match="joblib/executor level"):
        run_grid_campaign(cells, _worker_crashes_the_process_on_2, routing=routing, config=config)

    assert not routing.results_path.exists() or routing.results_path.read_text() == ""


# ---------------------------------------------------------------------------
# Resume-after-a-REAL-kill: a subprocess is SIGKILLed mid-run, restarted, and
# must resume correctly. NOT mocked -- a real process, a real signal.
# ---------------------------------------------------------------------------

_KILL_RESUME_SCRIPT = textwrap.dedent(
    """\
    import json, os, sys, time

    src = sys.argv[1]
    results_path = sys.argv[2]
    pid_log_path = sys.argv[3]
    n_cells = int(sys.argv[4])

    sys.path.insert(0, src)

    from cyclerfinder.search.campaign_runner import (
        CampaignRunnerConfig,
        CampaignRunnerRouting,
        CellOutcome,
        run_grid_campaign,
    )

    def worker(cell):
        # Real (non-trivial) CPU work per cell so a batch takes long enough
        # to reliably straddle the test's SIGKILL, without a fixed sleep.
        total = sum(i * i for i in range(3_000_000))
        with open(pid_log_path, "a") as f:
            f.write(json.dumps({"pid": os.getpid(), "index": cell["n"]}) + "\\n")
        return CellOutcome(status="miss", payload={"n": cell["n"], "chk": total})

    cells = [{"n": i} for i in range(n_cells)]
    routing = CampaignRunnerRouting(results_path=results_path)
    config = CampaignRunnerConfig(n_workers=2, checkpoint_batch_size=3, backend="loky")
    stats = run_grid_campaign(cells, worker, routing=routing, config=config)
    print("DONE", stats.evaluated_total, flush=True)
    """
)


def _find_src_path() -> str:
    here = Path(__file__).resolve()
    src = here.parents[2] / "src"
    assert src.is_dir(), f"src dir not found: {src}"
    return str(src)


def test_resume_after_real_sigkill(tmp_path: Path) -> None:
    """Kill the runner mid-campaign with SIGKILL; a fresh invocation must
    finish the grid with no duplicate work and no corruption.

    This is the single most important property per #788's own task spec:
    genuinely restartable across process kills, not just clean shutdowns.
    """
    results_path = tmp_path / "results.jsonl"
    pid_log_path = tmp_path / "pids.jsonl"
    n_cells = 30
    src_path = _find_src_path()

    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            _KILL_RESUME_SCRIPT,
            src_path,
            str(results_path),
            str(pid_log_path),
            str(n_cells),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Poll until SOME (but not all) cells have durably completed, then kill.
    deadline = time.monotonic() + 60.0
    partial_indices_before_kill: set[int] = set()
    while time.monotonic() < deadline:
        if results_path.exists():
            lines = [ln for ln in results_path.read_text(encoding="utf-8").splitlines() if ln]
            n_done = len(lines)
            if 0 < n_done < n_cells:
                partial_indices_before_kill = {json.loads(ln)["index"] for ln in lines}
                break
        if proc.poll() is not None:
            break  # finished (or died) before we ever caught it mid-run
        time.sleep(0.05)

    assert partial_indices_before_kill, (
        "never observed a partial (0 < n < n_cells) results.jsonl before the "
        "process finished -- cannot exercise the kill-mid-run path; widen "
        "n_cells or the per-cell work if this proves flaky"
    )
    assert proc.poll() is None, "process finished before we could kill it mid-run"

    with contextlib.suppress(ProcessLookupError):
        os.kill(proc.pid, signal.SIGKILL)
    proc.wait(timeout=10)

    # Sanity: fewer lines than the full grid at kill time (we actually caught
    # it mid-run, not after natural completion).
    lines_after_kill = [ln for ln in results_path.read_text(encoding="utf-8").splitlines() if ln]
    assert len(lines_after_kill) < n_cells

    # Evidence of no-redo (requirement 3): the PIDs of the SECOND run must
    # never process an index that was ALREADY durable before the kill.
    pid_log_before_resume = pid_log_path.read_text(encoding="utf-8").splitlines()
    pids_before_resume = {json.loads(ln)["pid"] for ln in pid_log_before_resume if ln}

    # Resume: a fresh invocation of the identical command.
    proc2 = subprocess.Popen(
        [
            sys.executable,
            "-c",
            _KILL_RESUME_SCRIPT,
            src_path,
            str(results_path),
            str(pid_log_path),
            str(n_cells),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = proc2.communicate(timeout=120)
    assert proc2.returncode == 0, f"resume run failed: stdout={out!r} stderr={err!r}"
    assert f"DONE {n_cells}" in out

    # Final state: every index 0..39 present EXACTLY ONCE (no duplicate
    # lines, i.e. no cell was re-run after already being durable).
    final_lines = [ln for ln in results_path.read_text(encoding="utf-8").splitlines() if ln]
    final_indices = [json.loads(ln)["index"] for ln in final_lines]
    assert len(final_indices) == n_cells, (
        f"expected exactly {n_cells} result lines, got {len(final_indices)} "
        "(duplicates would mean a completed cell was re-run after resume)"
    )
    assert sorted(final_indices) == list(range(n_cells))

    # Cross-check via the worker's own side-log (independent of results.jsonl):
    # every index that was ALREADY durable before the kill must NOT appear
    # again in the pid log entries written after the resume started (i.e. the
    # worker function itself was never invoked a second time for it).
    all_pid_log_lines = [json.loads(ln) for ln in pid_log_path.read_text().splitlines() if ln]
    seen_after_first_batch_of_resume = [
        rec for rec in all_pid_log_lines if rec["pid"] not in pids_before_resume
    ]
    resumed_indices = {rec["index"] for rec in seen_after_first_batch_of_resume}
    redone = partial_indices_before_kill & resumed_indices
    assert not redone, f"resume re-ran already-durable indices: {redone}"
