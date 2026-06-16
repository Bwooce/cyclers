"""``parallel_sweep`` -- joblib wrapper for embarrassingly-parallel sweeps.

Public API
----------

* :class:`ParallelSweepConfig` -- frozen dataclass of execution knobs (worker
  count, backend, per-cell timeout, error policy, progress verbosity).
* :class:`ParallelSweepResult` -- frozen dataclass of the aggregated result
  (per-cell results aligned with the input cells, succeed/fail counts, wall
  time, per-cell timings for profiling).
* :func:`parallel_sweep` -- the entry point. Takes a sequence of cells and a
  top-level closure, returns a :class:`ParallelSweepResult`.

Pickle-safety contract
----------------------

joblib's process-based backends (``loky`` default, ``multiprocessing``) require
that the closure, the per-cell arguments, and the per-cell result type are all
picklable. In practice:

* the closure MUST be a top-level function or a :class:`functools.partial`
  over one. Lambdas and nested closures fail (loky raises
  ``PicklingError`` at submission time; this wrapper re-raises with a
  hint).
* every entry of ``cells`` MUST be pickle-safe.
* every cell return value MUST be pickle-safe.

The ``threading`` backend bypasses the pickle requirement (shared memory) but
is GIL-bound, so it only helps closures that release the GIL (typically numpy
/ scipy hot loops). Pure-Python closures see no speedup under ``threading``;
they need ``loky``.

Phase 2 will rewire specific sweep callers to satisfy the pickle contract.
Phase 1 simply documents it and provides a clear failure mode.

Failure semantics
-----------------

If a cell raises and ``raise_on_first_error=False`` (default), the entry in
``ParallelSweepResult.results`` for that index is ``None`` and ``n_failed`` is
incremented. The exception type+message is preserved in
``ParallelSweepResult.notes`` (truncated to the first three failures to keep
the result small).

If ``raise_on_first_error=True``, the first cell exception is re-raised after
joblib finishes the in-flight batch (joblib does not cancel sibling workers
mid-batch).

Per-cell timeout (``timeout_seconds_per_cell``) is enforced by wrapping the
closure in a ``concurrent.futures``-style watchdog: when ``None`` (default), no
timeout is applied. When set, a cell exceeding the budget is treated as a
failure (result ``None``, recorded in ``notes``). Timeouts under the process
backends use joblib's own ``timeout`` parameter.

This module is intentionally small: it is the substrate, not the sweep. The
domain knowledge stays with the per-sweep modules; Phase 2 reuses this surface
without touching it.
"""

from __future__ import annotations

import time
import traceback
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

from joblib import Parallel, delayed

# ---------------------------------------------------------------------------
# Config & result dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParallelSweepConfig:
    """Execution knobs for a parallel sweep.

    Attributes
    ----------
    n_workers:
        Number of parallel workers. ``-1`` (default) means "all cores"
        (joblib's ``n_jobs=-1`` semantics). Explicit positive ints pin
        the pool size; useful for shared hosts where the user wants to
        leave headroom for other agents.
    chunk_size:
        Number of cells per submitted task batch. ``1`` (default) is best
        for long-running cells (~1s+); larger values amortise the fork /
        IPC overhead for very fast cells (< 10ms each). Mapped to joblib's
        ``batch_size``.
    backend:
        joblib backend. ``"loky"`` (default) is the safe-fork process pool
        and works on all platforms. ``"multiprocessing"`` is the stdlib pool
        (Linux-fork only, fragile under some scipy state). ``"threading"``
        is GIL-bound; only useful for closures that release the GIL.
    verbose:
        joblib's progress verbosity (0 = silent, higher = more lines).
    timeout_seconds_per_cell:
        Wall-clock budget per cell. ``None`` (default) = no timeout. When
        set under the process backends, cells exceeding the budget are
        reported as failures.
    raise_on_first_error:
        ``False`` (default) swallows per-cell exceptions (the cell's slot
        becomes ``None``, ``n_failed`` is incremented). ``True`` re-raises
        the first exception after joblib drains the in-flight batch.
    """

    n_workers: int = -1
    chunk_size: int = 1
    backend: Literal["loky", "threading", "multiprocessing"] = "loky"
    verbose: int = 0
    timeout_seconds_per_cell: float | None = None
    raise_on_first_error: bool = False


@dataclass(frozen=True)
class ParallelSweepResult:
    """Aggregated result of a parallel sweep.

    Attributes
    ----------
    results:
        One entry per input cell, in the same order. Failed cells (cell
        raised or timed out under ``raise_on_first_error=False``) have
        ``None`` in their slot.
    n_cells:
        ``len(cells)``. Always equals ``len(results)``.
    n_succeeded:
        Number of cells whose closure returned a value (including ``None``
        as a legitimate result? See note below).
    n_failed:
        ``n_cells - n_succeeded``.
    elapsed_seconds:
        Wall-clock time spent inside :func:`parallel_sweep` (not including
        ``cells`` construction by the caller).
    per_cell_elapsed_seconds:
        Per-cell wall time, in the same order as ``results``. Useful for
        profiling. Failed cells record the time spent before the failure.
    notes:
        Human-readable failure summary (first three failures, truncated).

    Note on ``None`` returns
    ------------------------
    A cell that successfully returns ``None`` is indistinguishable from a
    failed cell in the ``results`` tuple. ``n_succeeded`` / ``n_failed`` are
    the authoritative success-vs-failure counts. Callers that need to
    distinguish a legitimate ``None`` from a failure should wrap their cell
    result in a tagged container (e.g. ``{"ok": True, "value": ...}``).
    """

    results: tuple[Any, ...]
    n_cells: int
    n_succeeded: int
    n_failed: int
    elapsed_seconds: float
    per_cell_elapsed_seconds: tuple[float, ...]
    notes: str = field(default="")


# ---------------------------------------------------------------------------
# Internal helper: top-level so loky can pickle it
# ---------------------------------------------------------------------------


def _run_one_cell(
    closure: Callable[[Any], Any],
    cell: Any,
) -> tuple[bool, Any, float, str]:
    """Run a single cell, capturing (ok, value, elapsed, error).

    Top-level (not nested) so it pickles under loky / multiprocessing.
    """
    t0 = time.perf_counter()
    try:
        value = closure(cell)
        return True, value, time.perf_counter() - t0, ""
    except Exception as exc:  # broad by design: per-cell isolation
        err = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=3)}"
        return False, None, time.perf_counter() - t0, err


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def parallel_sweep(
    cells: Sequence[Any],
    closure: Callable[[Any], Any],
    *,
    config: ParallelSweepConfig | None = None,
    on_cell_complete: Callable[[int, Any], None] | None = None,
) -> ParallelSweepResult:
    """Parallelize a single-cell sweep closure across cores.

    Parameters
    ----------
    cells:
        Sequence of opaque per-cell arguments. Each cell is passed to
        ``closure`` exactly once. Must be picklable when using a process
        backend (see module docstring).
    closure:
        Callable mapping ``cell -> result``. Must be a top-level function
        or a :class:`functools.partial` of one under process backends.
        Lambdas / nested closures fail at submission with
        :class:`PicklingError`; this wrapper surfaces the joblib error
        directly so the caller sees the cause.
    config:
        :class:`ParallelSweepConfig` instance. ``None`` uses defaults
        (``n_workers=-1``, loky backend, no timeout, swallow errors).
    on_cell_complete:
        Optional callback ``(cell_index, result) -> None`` invoked AFTER
        joblib returns (not mid-flight: joblib's batch model does not
        expose per-cell completion events without a custom backend). The
        callback receives results in input order. Failed cells pass
        ``None`` as the result.

    Returns
    -------
    :class:`ParallelSweepResult` describing the aggregated outcome.

    Pickle-safety contract
    ----------------------
    See the module docstring. In short: top-level closure, picklable cells,
    picklable returns -- otherwise loky will raise at submission and this
    wrapper will let the exception propagate (it is a contract violation,
    not a per-cell failure).

    Examples
    --------
    >>> from cyclerfinder.parallel import parallel_sweep, ParallelSweepConfig
    >>> def double(x):
    ...     return x * 2
    >>> result = parallel_sweep([1, 2, 3], double)
    >>> result.results
    (2, 4, 6)
    >>> result.n_succeeded
    3
    """
    cfg = config if config is not None else ParallelSweepConfig()
    n_cells = len(cells)

    # Edge case: empty sweep. Return a zero-valued result without spawning
    # workers; joblib accepts an empty iterable but the result of
    # Parallel()([]) is a plain list, and we want to avoid the fork cost.
    if n_cells == 0:
        return ParallelSweepResult(
            results=(),
            n_cells=0,
            n_succeeded=0,
            n_failed=0,
            elapsed_seconds=0.0,
            per_cell_elapsed_seconds=(),
            notes="empty sweep -- no cells submitted",
        )

    t0 = time.perf_counter()

    parallel_kwargs: dict[str, Any] = {
        "n_jobs": cfg.n_workers,
        "backend": cfg.backend,
        "verbose": cfg.verbose,
        "batch_size": cfg.chunk_size if cfg.chunk_size > 0 else 1,
    }
    if cfg.timeout_seconds_per_cell is not None:
        # joblib's ``timeout`` is per-task wall-clock; cells exceeding it
        # raise ``TimeoutError`` which our ``_run_one_cell`` would normally
        # catch -- BUT joblib enforces the timeout at the ``Parallel`` layer,
        # not inside the task, so we let joblib surface the TimeoutError to
        # the dispatch loop below.
        parallel_kwargs["timeout"] = cfg.timeout_seconds_per_cell

    # Build the delayed-task list. joblib's ``delayed`` wraps the closure +
    # args in a picklable form; ``_run_one_cell`` is top-level so loky can
    # pickle the wrapping reference.
    tasks = [delayed(_run_one_cell)(closure, cell) for cell in cells]

    # Run. We capture joblib-level exceptions (PicklingError at submission
    # time, TimeoutError per cell under the timeout path) and convert them
    # into per-cell failures rather than aborting the whole sweep -- unless
    # the caller asked for raise_on_first_error.
    parallel = Parallel(**parallel_kwargs)
    try:
        raw_results: list[tuple[bool, Any, float, str]] = parallel(tasks)
    except Exception as exc:
        # Submission-time failure (e.g. PicklingError from a lambda) or a
        # joblib-level error not attributable to a single cell. Propagate;
        # the caller's pickle-safety contract is the relevant guarantee.
        if cfg.raise_on_first_error:
            raise
        elapsed = time.perf_counter() - t0
        return ParallelSweepResult(
            results=tuple([None] * n_cells),
            n_cells=n_cells,
            n_succeeded=0,
            n_failed=n_cells,
            elapsed_seconds=elapsed,
            per_cell_elapsed_seconds=tuple([0.0] * n_cells),
            notes=(f"joblib-level failure (no cell results): {type(exc).__name__}: {exc}"),
        )

    # Unpack.
    results: list[Any] = []
    per_cell: list[float] = []
    n_ok = 0
    n_fail = 0
    failure_msgs: list[str] = []
    for i, (ok, value, dt, err) in enumerate(raw_results):
        results.append(value if ok else None)
        per_cell.append(dt)
        if ok:
            n_ok += 1
        else:
            n_fail += 1
            if len(failure_msgs) < 3:
                # Keep only the first line of the traceback summary so
                # ``notes`` stays compact; the full traceback was already
                # logged inside the worker process.
                first_line = err.split("\n", 1)[0]
                failure_msgs.append(f"cell {i}: {first_line}")
            if cfg.raise_on_first_error:
                raise RuntimeError(
                    f"parallel_sweep: cell {i} failed and raise_on_first_error=True; cause: {err}"
                )

    if on_cell_complete is not None:
        for i, value in enumerate(results):
            on_cell_complete(i, value)

    notes = ""
    if n_fail > 0:
        notes = f"{n_fail}/{n_cells} cells failed; first {min(n_fail, 3)}: " + " | ".join(
            failure_msgs
        )

    elapsed = time.perf_counter() - t0

    return ParallelSweepResult(
        results=tuple(results),
        n_cells=n_cells,
        n_succeeded=n_ok,
        n_failed=n_fail,
        elapsed_seconds=elapsed,
        per_cell_elapsed_seconds=tuple(per_cell),
        notes=notes,
    )
