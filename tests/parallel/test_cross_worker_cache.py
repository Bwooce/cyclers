"""Regression-safety tests for cross-worker cache pre-warming (task #474).

The #472 commit memoized the deterministic per-leg cost functions with
:func:`functools.cache` for a ~28x single-process speedup. That cache is
PER-PROCESS, so under joblib's default ``"loky"`` backend (fresh-interpreter
workers) every worker rebuilds it cold and the speedup evaporates. Task #474's
fix: run the campaign under the ``"multiprocessing"`` backend — which on Linux
forks the pool — with the parent's caches pre-warmed *before* the fork, so the
workers inherit them copy-on-write.

These tests guard the two correctness properties (mirroring #472's discipline):

1. PARITY — the warmed/forked-worker path returns BIT-IDENTICAL results to the
   direct single-process path over representative + edge inputs. Pre-warming
   reads only frozen module constants, so it can change no value.
2. SHARING — the cache is actually shared/effective across the forked workers:
   a worker reports a NON-EMPTY cache at its very first call (it inherited the
   parent's warmed entries via COW), whereas under ``"loky"`` the same worker
   starts COLD. This is the analogue of #472's verified-cache-hit assertion.

The tests are Linux-fork specific (``multiprocessing`` forks on Linux). On a
platform whose default start method is not ``fork`` the SHARING test is skipped;
the PARITY test is platform-independent.
"""

from __future__ import annotations

import functools
import multiprocessing as mp
from collections.abc import Generator

import pytest

from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
from cyclerfinder.search import vilm
from cyclerfinder.search.cache_warm import DEFAULT_LEGS, warm_moon_leg_caches
from cyclerfinder.search.moon_prune import moon_leg_admissible

# ---------------------------------------------------------------------------
# Test isolation: reset loky's reusable executor before each test
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_loky_executor_cold() -> Generator[None, None, None]:
    """Shut down loky's global reusable executor before each test.

    The loky backend maintains a process pool with a reusable executor that
    persists across test runs. Worker processes retain their per-process cache
    state (e.g. functools.cache entries), so running a loky-based test after
    another one with warm caches causes the worker pool to be reused with warm
    caches already present. This fixture ensures a cold slate for each test by
    forcing executor shutdown before the test body runs, so subsequent loky
    calls spawn fresh workers.

    See task #481 and test_loky_workers_start_cold_control for the motivation.
    """
    try:
        from joblib.externals.loky import get_reusable_executor

        executor = get_reusable_executor()
        executor.shutdown(wait=True)
    except (ImportError, AttributeError):
        # loky not available; test will skip anyway if loky is needed
        pass
    yield


# Representative + edge domain (mirrors the bench skeleton; small for CI speed).
_LEGS = DEFAULT_LEGS
_VINF_GRID = [round(2.0 + 0.25 * k, 4) for k in range(6)]  # 2.0 .. 3.25
_BUDGET_GRID = [1.0, 3.0]

# Flatten to (primary, a, b, vinf, budget) cells for the sweep.
_CELLS = [
    (primary, a, b, vinf, budget)
    for (primary, a, b) in _LEGS
    for vinf in _VINF_GRID
    for budget in _BUDGET_GRID
]

_FORK_DEFAULT = mp.get_start_method(allow_none=False) == "fork"


# ---------------------------------------------------------------------------
# Top-level closure (pickle-safe under process backends)
# ---------------------------------------------------------------------------


def _admit_cell(
    cell: tuple[str, str, str, float, float],
) -> tuple[str, str, str, float, float, bool, str]:
    primary, a, b, vinf, budget = cell
    ok, reason = moon_leg_admissible(a, b, vinf_kms=vinf, budget_kms=budget, primary=primary)
    return (primary, a, b, vinf, budget, ok, reason)


def _worker_cache_currsize(cell: tuple[str, str, str, float, float]) -> int:
    """Report the worker's VILM floor-cache size at first call, THEN run the cell.

    Read the cache size BEFORE doing any work so a forked (warm) worker reports
    > 0 (inherited entries) while a spawned (cold) worker reports 0.
    """
    return vilm.vilm_dv_floor.cache_info().currsize


def _serial_reference() -> list[tuple[str, str, str, float, float, bool, str]]:
    """Direct single-process map — the bit-identical reference."""
    return [_admit_cell(c) for c in _CELLS]


# ---------------------------------------------------------------------------
# Test 1 — PARITY: warmed/forked path is bit-identical to the direct path
# ---------------------------------------------------------------------------


def test_prewarm_multiprocessing_bit_identical() -> None:
    """The pre-warmed multiprocessing sweep returns BIT-IDENTICAL results."""
    reference = _serial_reference()

    cfg = ParallelSweepConfig(
        n_workers=2,
        backend="multiprocessing",
        prewarm=functools.partial(
            warm_moon_leg_caches,
            legs=_LEGS,
            vinf_grid=_VINF_GRID,
            budget_grid=_BUDGET_GRID,
        ),
    )
    result = parallel_sweep(_CELLS, _admit_cell, config=cfg)

    assert result.n_failed == 0, result.notes
    assert list(result.results) == reference  # bit-identical, in order


def test_prewarm_loky_also_bit_identical() -> None:
    """Same parity holds under loky (cold workers) — value-correctness is
    backend-independent; only the SHARING benefit differs."""
    reference = _serial_reference()

    cfg = ParallelSweepConfig(
        n_workers=2,
        backend="loky",
        prewarm=functools.partial(
            warm_moon_leg_caches,
            legs=_LEGS,
            vinf_grid=_VINF_GRID,
            budget_grid=_BUDGET_GRID,
        ),
    )
    result = parallel_sweep(_CELLS, _admit_cell, config=cfg)

    assert result.n_failed == 0, result.notes
    assert list(result.results) == reference


def test_prewarm_none_is_noop_bit_identical() -> None:
    """``prewarm=None`` (legacy default) is unchanged and bit-identical."""
    reference = _serial_reference()
    cfg = ParallelSweepConfig(n_workers=2, backend="multiprocessing", prewarm=None)
    result = parallel_sweep(_CELLS, _admit_cell, config=cfg)
    assert result.n_failed == 0, result.notes
    assert list(result.results) == reference


# ---------------------------------------------------------------------------
# Test 2 — SHARING: forked workers inherit the warmed cache (COW)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _FORK_DEFAULT,
    reason="cross-worker COW cache sharing requires a fork start method (Linux)",
)
def test_forked_workers_inherit_warm_cache() -> None:
    """Under multiprocessing+prewarm a worker starts with a NON-EMPTY cache.

    The parent warms ``vilm.vilm_dv_floor`` before the fork; the children share
    those pages copy-on-write, so each worker's FIRST observation of
    ``cache_info().currsize`` is already > 0.
    """
    cfg = ParallelSweepConfig(
        n_workers=2,
        backend="multiprocessing",
        prewarm=functools.partial(
            warm_moon_leg_caches,
            legs=_LEGS,
            vinf_grid=_VINF_GRID,
            budget_grid=_BUDGET_GRID,
        ),
    )
    result = parallel_sweep(_CELLS, _worker_cache_currsize, config=cfg)

    assert result.n_failed == 0, result.notes
    # At least one leg-pair was warmed, so every worker saw a warm cache.
    assert min(result.results) > 0, (
        f"expected warm caches in every forked worker, got currsizes={result.results}"
    )


@pytest.mark.skipif(
    not _FORK_DEFAULT,
    reason="control assertion compares against the forked-warm case",
)
def test_loky_workers_start_cold_control() -> None:
    """CONTROL: loky workers (fresh interpreters) start COLD even WITH prewarm.

    This is the positive control proving the SHARING test above measures a real
    effect: the same pre-warm hook does NOT reach loky workers (they spawn fresh
    interpreters), so their first cache observation is 0. (The parent runs the
    prewarm, but loky's workers do not inherit the parent address space.)
    """
    cfg = ParallelSweepConfig(
        n_workers=2,
        backend="loky",
        prewarm=functools.partial(
            warm_moon_leg_caches,
            legs=_LEGS,
            vinf_grid=_VINF_GRID,
            budget_grid=_BUDGET_GRID,
        ),
    )
    result = parallel_sweep(_CELLS, _worker_cache_currsize, config=cfg)

    assert result.n_failed == 0, result.notes
    # Loky workers do not inherit the parent's warmed cache: cold start.
    assert min(result.results) == 0, f"expected cold loky workers, got currsizes={result.results}"
