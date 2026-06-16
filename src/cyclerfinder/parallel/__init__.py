"""Parallel-sweep substrate (#321 Phase 1).

Discovery sweeps in this repo (#281 multi-mu tulip, #312 Uranus extended, #320
quasi-cycler scan, #338 annual epoch sweep, #284 asymmetric scan) are CPU-bound
single-cell loops over scipy / Lambert / propagation kernels. The compute pattern
is embarrassingly parallel: each cell is independent, results are aggregated at
the end. Until #321 these sweeps were serial because no shared substrate
existed.

This package provides the ``parallel_sweep`` wrapper -- a thin layer over
``joblib.Parallel`` that:

* fans a single-cell closure across cores (configurable backend / worker count),
* surfaces per-cell timings for profiling,
* swallows per-cell exceptions by default (so a bad cell does not abort the
  whole sweep) while reporting an explicit failure count,
* documents the pickle-safety contract that all process-pool backends impose
  (top-level closure, picklable args, picklable return).

Phase 2 (separate task) rewires specific sweep callers to use this substrate.
Phase 1 only lays the substrate plus a smoke-test suite and a side-by-side
speedup demo on one existing sweep.
"""

from cyclerfinder.parallel.parallel_sweep import (
    ParallelSweepConfig,
    ParallelSweepResult,
    parallel_sweep,
)

__all__ = [
    "ParallelSweepConfig",
    "ParallelSweepResult",
    "parallel_sweep",
]
