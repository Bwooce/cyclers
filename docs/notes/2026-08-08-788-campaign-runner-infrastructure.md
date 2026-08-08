# #788 — generic checkpointed cell-grid campaign runner

2026-08-08. Foundational infrastructure task, unblocking `#789`-`#792` (the
combinatorial-search campaign menu). No science in this task — a toy worker
function only, per the task's own explicit scope boundary.

## What existed before

Two pieces of infrastructure, never combined:

* `src/cyclerfinder/search/discovery_campaign.py` — a checkpointed-restart
  pattern for ONE concrete genome (repeated-moon multi-rev cyclers): an
  append-only checkpoint file of completed candidate indices, `index %
  n_workers` sharding across N independently-launched OS processes (see
  `scripts/discovery_campaign_daemon.py`'s own docstring: `for w in 0 1 2 3;
  do ... --worker-id "$w" & done`), and method-versioned empty-region
  writeback.
* `src/cyclerfinder/parallel/parallel_sweep.py` — a joblib wrapper
  (measured 5x speedup on 8 workers) used by one-shot sweep scripts
  (`scripts/run_567_epoch_robustness_scan.py`, `scripts/run_338_parallel_demo.py`,
  etc.) with no supervising loop around it — every long sweep in this
  project has been "launch it, watch it, re-launch by hand if it dies."

Nothing combined them, and nothing in the project could be killed and
resumed as a single supervised process. That is what this task builds.

## What was built

`src/cyclerfinder/search/campaign_runner.py` — one new function,
`run_grid_campaign(cells, worker, *, routing, config, empty_region_spec)`,
plus its supporting dataclasses (`CellOutcome`, `CampaignRunnerConfig`,
`CampaignRunnerRouting`, `EmptyRegionSpec`, `GridCampaignStats`).

**Generic over "cell".** `cells` is any `Sequence[Any]` of picklable objects
(same pickle-safety contract as `parallel_sweep` itself); `worker` is a
top-level `Callable[[cell], CellOutcome]`. The module has zero knowledge of
resonance ratios, moon sequences, or any domain content — `#789`-`#792` each
supply their own grid + worker, unmodified runner.

**One process, parallel across cores.** Internally batches the pending cells
(`config.checkpoint_batch_size` per round-trip, default 25) through
`parallel_sweep`, so a single `run_grid_campaign()` invocation IS the
parallel pool — no more hand-launched sibling processes.

**`results.jsonl` is its own checkpoint.** A JSONL line
`{"index", "status", "payload", "error"}` per cell, written with
`flush()` + `os.fsync()` the instant each batch returns, BEFORE the next
batch starts. An index is "done" iff a complete, parseable line for it
exists in the file — checkpoint and result are the same durable write, so
there is no two-file race between "wrote the result" and "marked it done".
This is a deliberate departure from `discovery_campaign.py`'s two-file
scheme (plain-index checkpoint + separate sink file), which is left
UNCHANGED: its format genuinely differs, it is tested and stable, a sibling
task (`#793`) may be committing to the repo concurrently, and there was no
generic logic left to extract beyond what already lives in
`data.empty_regions` / `data.method_capability` — both reused directly
(not reimplemented) by the new module.

**Durability granularity is per-BATCH, not literally per-cell** — documented
plainly rather than oversold. `parallel_sweep` only surfaces results after
its whole `Parallel()` call returns (mid-flight streaming would need a
custom joblib backend, out of scope). A kill loses at most
`checkpoint_batch_size` in-flight cells, never a completed batch; that knob
is the caller's dial for the loss window.

**The joblib-level (whole-batch) failure trap, and how it's avoided.** Two
distinct `parallel_sweep` failure modes exist: (1) one cell's worker raised
— that cell's slot is `None`, the rest of the batch is real; (2) the WHOLE
`Parallel()` call failed at the executor level (a worker process
OOM-killed/segfaulted, a `PicklingError` at submission) — `parallel_sweep`
then reports **every** cell in that call as `None`, prefixed
`"joblib-level failure"` in `notes`. Treating case 2 like case 1 would
durably record cells that were NEVER evaluated as `"error"` — permanent
silent data loss on a weeks-long campaign, since resume skips anything with
a results line. `run_grid_campaign` detects the prefix, refuses to write
ANY line for that batch, and raises loudly so the batch retries in full next
invocation. Covered by
`test_whole_batch_joblib_level_failure_writes_nothing_and_raises`, which
triggers a REAL `TerminatedWorkerError` by having one cell call
`os._exit(1)` inside its own worker process.

**Empty-region emission is an aggregate**, matching `discovery_campaign.py`'s
own convention: one `EmptyRegionReport` (unchanged schema from
`data.empty_regions`, validated by the project's real
`validate_empty_region`) written once the ENTIRE grid has a durable "miss"
result (across however many resumed invocations that took) with zero hits
and zero errors — not one record per cell. Idempotent on `region_id` (a
completed, all-miss campaign re-invoked a second time does not double-write
the registry). Gated up front by the existing capability-subsumption check
(`data.method_capability.should_sweep`): a weaker-or-equal method
re-sweeping a region a stronger method already emptied is a no-op.

## Worked example (from the test suite)

```python
from cyclerfinder.search.campaign_runner import (
    CampaignRunnerConfig, CampaignRunnerRouting, CellOutcome,
    EmptyRegionSpec, run_grid_campaign,
)
from cyclerfinder.data.method_capability import MethodCapability

def worker(cell: dict) -> CellOutcome:
    n = cell["n"]
    if n % 7 == 0 and n != 0:
        return CellOutcome(status="hit", payload={"n": n})
    return CellOutcome(status="miss", payload={"n": n})

cells = [{"n": i} for i in range(1000)]
routing = CampaignRunnerRouting(
    results_path=Path("out/campaign_788_demo/results.jsonl"),
    empty_regions_path=Path("data/empty_regions.jsonl"),  # real registry, real campaigns only
)
spec = EmptyRegionSpec(
    region_id="demo-region-788",
    family="toy demo grid",
    centre="nowhere",
    method_capability=MethodCapability(
        genome="toy", corrector="none", capability_tags=frozenset({"toy"}), git_sha="...",
    ),
    prune_gates=("n % 7 == 0",),
)
config = CampaignRunnerConfig(n_workers=-1, checkpoint_batch_size=50,
                               timeout_seconds_per_cell=300.0)
stats = run_grid_campaign(cells, worker, routing=routing, config=config,
                           empty_region_spec=spec)
print(stats.as_dict())
```

Re-invoking the identical call (same `results_path`) after a kill at any
point — clean exit, SIGKILL, `kill -9` from another terminal, machine reboot
— resumes from exactly where it left off. A future campaign script wraps
this in an `argparse` CLI + its own grid/worker, the same shape as
`scripts/discovery_campaign_daemon.py` today but WITHOUT the `for w in
0 1 2 3; do ... & done` hand-launch loop.

## Kill/resume test: real evidence

`tests/search/test_campaign_runner.py::test_resume_after_real_sigkill`
spawns a real Python subprocess running `run_grid_campaign` over 30 cells
(loky backend, 2 workers, `checkpoint_batch_size=3`), polls
`results.jsonl` until it observes a genuine partial state (more than zero,
fewer than 30 durable lines), sends a REAL `SIGKILL` (`os.kill(pid,
signal.SIGKILL)`, not a mock), then launches a second, independent
subprocess invocation of the identical command and asserts: (a) it exits 0
and reports all 30 cells done, (b) `results.jsonl` has EXACTLY 30 lines
covering indices 0..29 with no duplicates, and (c) an independent
worker-side PID log proves no index that was already durable before the
kill was ever re-processed by the resumed run's worker PIDs. This test
passed on 4 consecutive runs during development (3 repeats of the isolated
test plus the full-file run), including under `pytest-xdist` parallel
collection (8 workers on this machine).

A hand-run, narrated version of the same scenario (outside pytest, so the
`print()` timeline is visible) — this is the actual transcript, not a
paraphrase:

```
[t=539797.52] caught mid-run: 3/30 cells done, indices=[0, 1, 2]
SIGKILL pid=90247 now
after kill: 3/30 lines in results.jsonl
resuming with a fresh process invocation...
resume stdout: DONE 30
resume returncode: 0
final: 30 lines, indices == range(n_cells)? True
duplicate indices present: False
```

3 cells were durable at the moment of the SIGKILL (out of a
`checkpoint_batch_size=3` first batch — exactly the "at most one in-flight
batch is at risk" bound this design targets); the process was actually
killed mid-flight (`3/30`, not `30/30` — the test itself asserts this, so a
future accidental widening of the timing window that makes the process
finish before the kill lands would FAIL the test rather than silently pass);
the resumed process picked up cleanly and finished with no duplicate or
missing indices.

`tests/search/test_campaign_runner.py` also covers (all passing):

* `test_durable_per_cell_writes_and_hit_detection` — every cell lands as its
  own well-formed JSON line; hits are correctly identified.
* `test_parallel_execution_uses_multiple_worker_pids` — work is actually
  spread across >1 OS process when `os.cpu_count() >= 2`.
* `test_resume_skips_completed_cells_and_finishes_exactly_once` — an
  in-process capped-then-uncapped resume completes the grid with exactly
  one line per index (no duplicate work), without needing a real kill.
* `test_empty_region_emitted_only_on_full_all_miss_completion` +
  `test_a_single_hit_blocks_empty_region_emission` — the aggregate
  empty-region record is emitted exactly once, only on full all-miss
  completion, and its payload round-trips through the real
  `validate_empty_region`.
* `test_capability_subsumption_gate_skips_a_subsumed_resweep` — a prior
  stronger-method empty-region record correctly no-ops a weaker resweep.
* `test_per_cell_exception_is_recorded_as_a_durable_error` — an isolated
  per-cell exception is durable and does not affect sibling cells.
* `test_whole_batch_joblib_level_failure_writes_nothing_and_raises` — the
  joblib-level-failure trap described above, verified with a real
  `os._exit(1)` inside a worker process, not a mock.

## Verification

* `uv run pytest tests/search/test_campaign_runner.py -q` — 9/9 passed
  (repeated 4x during development, no flakes observed).
* `uv run pytest tests/data tests/search -q` — the mandated ratchet. 2
  failures observed (`test_eggie_ballistic.py::
  test_gate_b_table4_vinf_reached_but_subsurface`,
  `test_504_pluto_charon_kk_sweep.py::test_504_sweep_33`), both traced to
  the concurrently-running `#793` session's own uncommitted mid-edit of
  `lambert.py`/`catalogue.yaml` at the moment the suite ran, not this task's
  diff — neither failing test's import graph reaches `campaign_runner.py`
  (a brand-new, unreferenced module), and this task touched no other
  existing source file. Confirmed clean via `git log` that `#793` committed
  its `lambert.py` fix (`efbe7239`) in the same window; re-running
  `tests/search/test_campaign_runner.py` alone post-commit still passed
  9/9. This task makes no `catalogue.yaml` change.
* `uv run ruff check .` / `uv run ruff format --check .` — clean on the
  new/changed files.
* `uv run mypy src tests` — clean (0 errors across all 839 checked source
  files, the full canonical invocation, not a single-file check).

## What this is NOT

Per the task's own explicit scope boundary: no real campaign was built or
run. `#789` (Resonant Atlas), `#790` (chain-itinerary enumeration), `#791`
(moon-tour sequences), `#792` (Uranus census) each still need their own
grid + worker function written against this runner — none of that domain
work happened here. Multi-MACHINE sharding (beyond multi-core, within one
process) is not built either; it was deliberately left out as unneeded
complexity for now — a future need can partition `cells` per machine with
distinct `results_path`/`empty_regions_path` values, which the current
design already supports without modification (each machine's `results.jsonl`
is independent; a merge step would just union the index sets before
computing final stats, or `#789`-style campaigns can simply set `region_id`
per shard).

## Follow-on task registered

`#795` — filed for a small `argparse` CLI wrapper / worked example script
demonstrating `run_grid_campaign` end-to-end against a real (non-toy) toy
grid under `scripts/`, once the first real campaign (`#789`+) needs one; not
built here since the task explicitly scoped "demonstrated with a
toy/example worker function" as the test suite, not a new script (a new
`scripts/run_*.py` file would also trip the `tests/scripts` preflight AST
ratchet, out of this task's mandated verification scope).
