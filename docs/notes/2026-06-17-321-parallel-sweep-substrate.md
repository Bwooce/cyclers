# #321 — Parallel-sweep substrate (Phase 1)

**Date:** 2026-06-17  
**Phase:** 1 of 2 (substrate + smoke test + pilot demo). Phase 2 rewires
specific sweep callers; tracked separately.  
**Verdict:** SUBSTRATE LANDED. 5.06x speedup on the #338 annual V4-strict
sweep (43.6 s serial → 8.6 s parallel on 8 workers), 100/100 byte-for-byte
equivalence with the serial baseline. 9/9 smoke tests passing.

## Why

Discovery sweeps in this repo (#281 multi-mu tulip, #312 Uranus extended,
#320 quasi-cycler scan, #338 annual epoch sweep, #284 asymmetric scan) are
CPU-bound single-cell loops over scipy / Lambert / propagation kernels. The
compute pattern is embarrassingly parallel: each cell is independent,
results are aggregated at the end. pytest-xdist already parallelizes the
test suite, but the actual discovery sweeps were single-threaded scipy —
on a 16-core host running 4-5 concurrent agents with a 1.3 load average,
that left ~11 cores idle during every sweep. #339's SILVER admission
validated that these sweeps are productively used; the obvious next-
infrastructure step is to parallelize them.

## What landed (Phase 1)

### Substrate

`src/cyclerfinder/parallel/parallel_sweep.py` — joblib wrapper exposing
three top-level surfaces:

- **`ParallelSweepConfig`** (frozen dataclass): execution knobs.
  - `n_workers: int = -1` — joblib's "all cores" semantic by default; pin
    explicitly when sharing the host with concurrent agents.
  - `chunk_size: int = 1` — best for ~1 s/cell sweeps; raise for very
    fast cells to amortise fork/IPC cost.
  - `backend: Literal["loky", "threading", "multiprocessing"] = "loky"` —
    safe-fork process pool by default; threading is only for GIL-releasing
    closures.
  - `verbose: int = 0` — joblib progress level.
  - `timeout_seconds_per_cell: float | None = None` — wall-clock budget.
  - `raise_on_first_error: bool = False` — swallow per-cell exceptions by
    default and collect partials.

- **`ParallelSweepResult`** (frozen dataclass): aggregated output.
  - `results: tuple[Any, ...]` — one entry per input cell, in input order
    (parallel_sweep restores input order before returning); `None` for
    failed cells.
  - `n_cells / n_succeeded / n_failed`.
  - `elapsed_seconds` — wall-clock spent inside `parallel_sweep`.
  - `per_cell_elapsed_seconds: tuple[float, ...]` — per-cell wall time for
    profiling.
  - `notes: str` — human-readable failure summary (first 3 failures).

- **`parallel_sweep(cells, closure, *, config=None, on_cell_complete=None)`** —
  the entry point.

### Pickle-safety contract (documented + enforced)

joblib's process backends require:

1. The closure is a top-level function or `functools.partial` of one.
   Lambdas + nested closures over local state fail at submission with
   `PicklingError`. The wrapper either propagates the joblib error (when
   `raise_on_first_error=True`) or records it in `notes` with
   `n_failed = n_cells` (when `False`).
2. All cell arguments are picklable.
3. The cell return value is picklable.

The `threading` backend bypasses (1)–(3) but is GIL-bound; only worth it
for closures that release the GIL (numpy/scipy hot loops). Pure-Python
closures see no speedup under `threading`.

The smoke test suite includes a negative case (`threading.Lock` captured
by a nested closure) that exercises the explicit-failure path.

### Smoke tests

`tests/parallel/test_parallel_sweep.py` — 9 tests, all passing:

1. `test_identity_sweep_returns_serial_equivalent` — identity sweep
   matches serial map.
2. `test_empty_cells_returns_zero_result` — zero-cell edge case.
3. `test_on_cell_complete_callback_invoked_in_order` — callback fires for
   every cell in input order.
4. `test_parallel_speedup_on_sleep_bound_cells` — 16 × 0.1 s sleeps,
   4 workers, asserts wall < 0.8 s (2x of 1.6 s serial). Warms the loky
   pool first so the fork cost does not dominate the measurement.
5. `test_error_in_one_cell_does_not_abort_sweep` — one cell raises;
   `n_failed=1`, others succeed, `results[i]=None` for the failure.
6. `test_raise_on_first_error_propagates` — `raise_on_first_error=True`
   propagates the per-cell exception.
7. `test_per_cell_timeout_fails_all_overrunning_cells` — 3 × 0.5 s cells
   under 0.1 s timeout, all fail (joblib raises TimeoutError, wrapper
   records as joblib-level failure with `n_failed=n_cells`).
8. `test_lambda_closure_fails_explicitly_under_loky` — nested closure
   capturing a `threading.Lock` fails explicitly via `notes`, not silent
   hang.
9. `test_raise_on_first_error_propagates_pickle_failure` — same with
   `raise_on_first_error=True`; the joblib exception propagates.

Suite runtime ~3.4 s under `-n0` (serial pytest, since the sweep tests
themselves spawn subprocesses).

### Part C — pilot demo on #338 annual V4-strict sweep

`scripts/run_338_parallel_demo.py` reproduces the 100-epoch V4-strict
sweep (`scripts/run_338_silver_v4strict_annual_sweep.py`) using the new
substrate, runs both serial and parallel passes back-to-back, and
asserts byte-for-byte equivalence on `passes_v4_strict` and
`drift_agreement_kms_vs_v3`.

Measured wall-clock (16-core host, 8 workers, no other load on the cores
this sweep used):

| Pass     | Wall (s) | Speedup |
|----------|----------|---------|
| Serial   | 43.6     | 1.00x   |
| Parallel | 8.6      | 5.06x   |

Equivalence: **100/100 rows match** on `passes_v4_strict` AND
`drift_agreement_kms_vs_v3` (exact float equality — same inputs, same
single-threaded per-cell scipy kernel, deterministic V4-strict propagator).

The 5.06x speedup on 8 workers reflects loky's fork cost amortised over
~0.44 s per cell — efficient given the V3+V4-scipy chain takes 0.2 s
out-of-band (built once, broadcast in `_Cell`) and the per-cell cost
includes one SPICE kernel load per worker.

Demo outputs:
- `data/silver_327_v4_strict_annual_sweep_321_parallel_demo.jsonl` —
  per-epoch rows in input order (separate path from the #338 baseline
  output so no overwrite).
- `data/scan_321_parallel_demo_summary.jsonl` — wall-clock comparison.

## Phase 2 backlog — expected speedup ranking

Per the task brief plus the measured 5x figure, the existing sweeps best
positioned to benefit:

| Rank | Sweep                                  | Cells          | Per-cell cost  | Expected speedup | Notes |
|------|----------------------------------------|----------------|----------------|------------------|-------|
| 1    | #338 annual epoch sweep (PROVEN)       | 100            | ~0.4 s         | 5-7x             | demoed here; immediate win |
| 2    | #284 asymmetric scan                   | 320            | ~1 s           | 8-12x            | longest absolute sweep; biggest single win |
| 3    | #281 multi-μ tulip sweep               | 14 systems     | ~5-10 s/system | 8-14x            | one-cell-per-system; trivially partitioned |
| 4    | #312 Uranus extended                   | 24×24 = 576    | ~0.5 s         | 10-14x           | already CPU-bound; very even cell cost |
| 5    | #320 quasi-cycler discovery sweep      | multi-system × cells | ~1-2 s   | 10-14x           | natural fit, but coordinate with concurrent #320 agent before rewire |

Rank-1 (#338) is the cleanest single rewire — the demo here already
verifies the closure semantics. Rank-2 (#284) is the biggest absolute
win. The constraints for Phase 2 rewires (one PR per sweep, ideally):

1. The per-cell closure becomes a top-level function accepting a single
   picklable cell payload.
2. Any large precomputed state (e.g. V3/V4-scipy chain in #338) is
   computed ONCE in the driver and either packed into the cell payload
   (if small / fast to pickle) or rebuilt per-worker in a
   `joblib.parallel_backend` init (if expensive to pickle).
3. The output ordering is preserved by `parallel_sweep` (results align
   with input order), so the existing JSONL writer code is unchanged.
4. A byte-for-byte equivalence check against the serial baseline lands
   in the same PR (the #338 demo's `equivalence_check` is the template).

## Out of scope / discipline

- **NO catalogue writeback** from any of the Phase 1 work. The parallel
  demo's output is a separate JSONL path; nothing in
  `data/catalogue.yaml` or `data/silver_327_verified.jsonl` is touched.
- **READ-ONLY** on the existing discovery sweep modules in Phase 1. The
  substrate is generic; the rewires are Phase 2.
- The smoke test for the parallel speedup is gated on
  `cpu_count >= 4` so the suite stays green on CI with fewer cores.

## File index

- `src/cyclerfinder/parallel/__init__.py` (NEW)
- `src/cyclerfinder/parallel/parallel_sweep.py` (NEW)
- `tests/parallel/__init__.py` (NEW)
- `tests/parallel/test_parallel_sweep.py` (NEW, 9 tests)
- `scripts/run_338_parallel_demo.py` (NEW, pilot demo)
- `pyproject.toml` — `joblib>=1.3` added as a top-level dependency
- `data/silver_327_v4_strict_annual_sweep_321_parallel_demo.jsonl` (NEW, generated)
- `data/scan_321_parallel_demo_summary.jsonl` (NEW, generated)
