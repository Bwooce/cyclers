"""Smoke tests for ``cyclerfinder.parallel.parallel_sweep`` (#321 Part B).

These tests verify the substrate behaviour:

* identity sweep semantics (same output as serial map),
* parallel speedup on a sleep-bound closure (skipped if cpu_count < 4),
* per-cell error swallowing under ``raise_on_first_error=False``,
* per-cell timeout policy,
* pickle-safety negative test (lambda must fail explicitly under loky).

They do NOT exercise any sweep callers; that is Phase 2.

All helper closures are top-level functions so loky can pickle them. This
is the same contract the Phase 2 rewires will satisfy.
"""

from __future__ import annotations

import os
import time

import pytest

from cyclerfinder.parallel import (
    ParallelSweepConfig,
    ParallelSweepResult,
    parallel_sweep,
)

# ---------------------------------------------------------------------------
# Top-level helper closures (pickle-safe)
# ---------------------------------------------------------------------------


def _double(x: int) -> int:
    return x * 2


def _sleep_then_return(x: float) -> float:
    """Sleep ``x`` seconds, return ``x``. Used for speedup + timeout tests."""
    time.sleep(x)
    return x


def _raise_on_index_2(x: int) -> int:
    """Raise ValueError for x == 2; otherwise return x*10."""
    if x == 2:
        raise ValueError("cell index 2 is poisoned")
    return x * 10


# ---------------------------------------------------------------------------
# Test 1 — identity sweep
# ---------------------------------------------------------------------------


def test_identity_sweep_returns_serial_equivalent() -> None:
    cells = [1, 2, 3, 4, 5]
    result = parallel_sweep(cells, _double)
    assert isinstance(result, ParallelSweepResult)
    assert result.results == (2, 4, 6, 8, 10)
    assert result.n_cells == 5
    assert result.n_succeeded == 5
    assert result.n_failed == 0
    assert result.notes == ""
    # Each per-cell time recorded; non-negative.
    assert len(result.per_cell_elapsed_seconds) == 5
    assert all(t >= 0.0 for t in result.per_cell_elapsed_seconds)


def test_empty_cells_returns_zero_result() -> None:
    result = parallel_sweep([], _double)
    assert result.results == ()
    assert result.n_cells == 0
    assert result.n_succeeded == 0
    assert result.n_failed == 0
    assert result.elapsed_seconds == 0.0


def test_on_cell_complete_callback_invoked_in_order() -> None:
    seen: list[tuple[int, int]] = []

    def cb(idx: int, value: object) -> None:
        # Cast for the test; the API uses ``Any``.
        seen.append((idx, value))  # type: ignore[arg-type]

    result = parallel_sweep([10, 20, 30], _double, on_cell_complete=cb)
    assert result.results == (20, 40, 60)
    assert seen == [(0, 20), (1, 40), (2, 60)]


# ---------------------------------------------------------------------------
# Test 2 — parallel speedup
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    (os.cpu_count() or 1) < 4,
    reason="needs >= 4 logical cores for the 2x speedup assertion",
)
def test_parallel_speedup_on_sleep_bound_cells() -> None:
    """16 cells of 0.1s each: serial 1.6s, 4-worker target < 0.8s.

    We warm the loky pool with a small no-op sweep before timing so the
    fork cost does not dominate. The actual cost we measure is the steady
    state per-cell cost across workers.
    """
    cfg = ParallelSweepConfig(n_workers=4, backend="loky")
    # Warm up the loky pool (first call pays the worker-spawn cost).
    parallel_sweep([0.01] * 4, _sleep_then_return, config=cfg)

    cells = [0.1] * 16
    serial_budget = 0.1 * 16  # 1.6s
    speedup_target = serial_budget / 2.0  # 2x = 0.8s

    t0 = time.perf_counter()
    result = parallel_sweep(cells, _sleep_then_return, config=cfg)
    wall = time.perf_counter() - t0

    assert result.n_succeeded == 16
    assert result.n_failed == 0
    assert result.results == tuple([0.1] * 16)
    assert wall < speedup_target, (
        f"parallel wall {wall:.3f}s exceeded {speedup_target:.3f}s budget "
        f"(2x speedup on {serial_budget:.3f}s serial). "
        f"per_cell sum = {sum(result.per_cell_elapsed_seconds):.3f}s"
    )


# ---------------------------------------------------------------------------
# Test 3 — error handling (swallow + collect partials)
# ---------------------------------------------------------------------------


def test_error_in_one_cell_does_not_abort_sweep() -> None:
    cells = [0, 1, 2, 3, 4]
    result = parallel_sweep(cells, _raise_on_index_2)
    assert result.n_cells == 5
    assert result.n_succeeded == 4
    assert result.n_failed == 1
    # The failing cell's slot is None, others are x*10.
    assert result.results == (0, 10, None, 30, 40)
    assert "ValueError" in result.notes
    assert "cell 2" in result.notes


def test_raise_on_first_error_propagates() -> None:
    cells = [0, 1, 2, 3, 4]
    cfg = ParallelSweepConfig(raise_on_first_error=True)
    with pytest.raises(RuntimeError, match="cell 2 failed"):
        parallel_sweep(cells, _raise_on_index_2, config=cfg)


# ---------------------------------------------------------------------------
# Test 4 — per-cell timeout
# ---------------------------------------------------------------------------


def test_per_cell_timeout_fails_all_overrunning_cells() -> None:
    """Each cell sleeps 0.5s; timeout budget is 0.1s. Expect all failures.

    Under loky, the timeout terminates the worker; joblib surfaces the
    failure to our submission loop and we convert it to a per-cell None.
    The exact failure path lands in ``notes`` (joblib-level message or
    per-cell traceback) -- we assert only the counts here.
    """
    cells = [0.5] * 3
    cfg = ParallelSweepConfig(
        n_workers=2,
        backend="loky",
        timeout_seconds_per_cell=0.1,
        raise_on_first_error=False,
    )
    result = parallel_sweep(cells, _sleep_then_return, config=cfg)
    # All cells should be failures: each cell sleeps 0.5s but the worker
    # is killed at 0.1s. Under joblib's process-pool semantics this is a
    # whole-sweep failure (joblib raises TimeoutError once any task overruns
    # and ``n_jobs > 1``), so the wrapper records it as a joblib-level
    # failure rather than per-cell. Either way: n_succeeded == 0.
    assert result.n_succeeded == 0
    assert result.n_failed == 3
    assert all(r is None for r in result.results)
    assert "Timeout" in result.notes or "timeout" in result.notes.lower()


# ---------------------------------------------------------------------------
# Test 5 — pickle-safety negative (lambda must fail explicitly)
# ---------------------------------------------------------------------------


def test_lambda_closure_fails_explicitly_under_loky() -> None:
    """Lambdas are not picklable by stdlib pickle; loky has cloudpickle and
    actually accepts them. But nested closures over local state DO fail.
    The contract we surface is that the failure is explicit, not a silent
    hang.

    For the negative test we use a closure that references a non-picklable
    local object (``threading.Lock``) -- loky's cloudpickle cannot serialise
    this and raises at submission time.
    """
    import threading

    lock = threading.Lock()

    def closure_with_unpicklable_capture(x: int) -> int:
        # Reference ``lock`` in the closure so cloudpickle has to capture it.
        with lock:
            return x * 2

    # raise_on_first_error=False is the default; we expect the wrapper to
    # surface the joblib-level error in ``notes``, not raise. n_failed
    # should equal n_cells (no per-cell results were produced).
    result = parallel_sweep([1, 2, 3], closure_with_unpicklable_capture)
    assert result.n_failed == 3
    assert result.n_succeeded == 0
    # notes should mention the picklability failure or contain a joblib
    # error string -- we keep the assertion loose because the exact error
    # text varies across joblib versions.
    assert len(result.notes) > 0


def test_raise_on_first_error_propagates_pickle_failure() -> None:
    """Same as above but with ``raise_on_first_error=True``: the joblib
    PicklingError (or equivalent) should propagate."""
    import threading

    lock = threading.Lock()

    def closure_with_unpicklable_capture(x: int) -> int:
        with lock:
            return x * 2

    cfg = ParallelSweepConfig(raise_on_first_error=True)
    with pytest.raises(Exception):  # noqa: B017 -- exact type varies by joblib
        parallel_sweep([1, 2, 3], closure_with_unpicklable_capture, config=cfg)
