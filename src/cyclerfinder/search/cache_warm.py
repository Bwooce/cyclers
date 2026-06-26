"""Cross-worker cache pre-warming for the #472 per-leg cost memoizations (#474).

The #472 commit (54f8252) memoized the deterministic per-leg cost functions with
:func:`functools.cache` for a ~28x single-process speedup. But ``functools.cache``
is PER-PROCESS, and the discovery campaigns fan out under joblib. The default
``"loky"`` backend spawns FRESH interpreters (cold caches), so every worker
rebuilds the cache from scratch and the speedup evaporates.

The fix (task #474): run a campaign under the ``"multiprocessing"`` backend,
which on Linux forks the worker pool from the parent. Forked children share the
parent's pages copy-on-write — so if the parent's ``functools.cache`` tables are
already populated BEFORE the fork, every worker inherits them for free (the
read-only cache pages are shared, never duplicated; ~free with the host's RAM).

This module is the pre-warm hook. :func:`warm_moon_leg_caches` drives the public
composite :func:`cyclerfinder.search.moon_prune.moon_leg_admissible` over a
campaign's discrete ``(moon-pair, V∞, budget)`` domain, which transitively
populates EVERY #472 cache (VILM, Tisserand, bend, correct) without this module
having to hand-replicate each memoized signature — so it cannot drift out of sync
when a new cached function is added downstream.

Purity contract
---------------
The warmed functions read only the *frozen* module-constant body tables
(``SATELLITES`` / ``PRIMARIES`` / ``PLANETS``), so a warmed cache is
BIT-IDENTICAL to a cold-built one. Pre-warming changes no result; it only moves
the cost build out of every worker and into the shared parent once. This is the
same purity contract the #472 memoizations already satisfy and that the parity
tests in ``tests/parallel/test_cross_worker_cache.py`` assert.

Usage
-----
Pass the bound warmer as the ``prewarm`` hook of a multiprocessing-backed
:class:`~cyclerfinder.parallel.ParallelSweepConfig`::

    import functools
    from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
    from cyclerfinder.search.cache_warm import warm_moon_leg_caches

    cfg = ParallelSweepConfig(
        backend="multiprocessing",  # forks on Linux -> COW-shared caches
        prewarm=functools.partial(
            warm_moon_leg_caches,
            legs=legs, vinf_grid=vinf_grid, budget_grid=budget_grid,
        ),
    )
    result = parallel_sweep(cells, closure, config=cfg)
"""

from __future__ import annotations

from collections.abc import Sequence

from cyclerfinder.search.moon_prune import moon_leg_admissible

# A representative default leg set spanning the Jovian + Saturnian adjacent-moon
# pairs the moon-tour campaigns sweep (mirrors scripts/bench_472_memoization.py).
DEFAULT_LEGS: tuple[tuple[str, str, str], ...] = (
    ("Jupiter", "Io", "Europa"),
    ("Jupiter", "Europa", "Ganymede"),
    ("Jupiter", "Ganymede", "Callisto"),
    ("Saturn", "Titan", "Rhea"),
    ("Saturn", "Rhea", "Dione"),
    ("Saturn", "Dione", "Tethys"),
)


def warm_moon_leg_caches(
    *,
    legs: Sequence[tuple[str, str, str]] = DEFAULT_LEGS,
    vinf_grid: Sequence[float],
    budget_grid: Sequence[float],
) -> None:
    """Populate the #472 per-leg caches over the given discrete arg domain.

    Drives :func:`moon_leg_admissible` over the cartesian product of
    ``legs`` x ``vinf_grid`` x ``budget_grid``. Each call transitively touches
    every #472-memoized cost function (VILM floor, Tisserand linkability, bend
    feasibility), so on return the parent process holds a fully-warmed cache for
    that domain. Forked children then inherit it copy-on-write.

    Parameters
    ----------
    legs:
        Sequence of ``(primary, moon_a, moon_b)`` triples. Defaults to
        :data:`DEFAULT_LEGS`.
    vinf_grid:
        Discrete V∞ values (km/s) the campaign evaluates per leg.
    budget_grid:
        Discrete ΔV budgets (km/s) the campaign evaluates per leg.

    Notes
    -----
    Pure / idempotent: re-running it only produces cache HITS after the first
    pass. Exceptions from any single ``(leg, V∞, budget)`` cell are swallowed —
    a body absent from the registry must not abort warming of the rest; an
    un-warmed cell simply pays its (deterministic, identical) cost lazily in the
    worker.
    """
    for primary, moon_a, moon_b in legs:
        for vinf in vinf_grid:
            for budget in budget_grid:
                try:
                    moon_leg_admissible(
                        moon_a,
                        moon_b,
                        vinf_kms=vinf,
                        budget_kms=budget,
                        primary=primary,
                    )
                except Exception:  # warming is best-effort
                    # A missing body / out-of-domain cell must not abort the
                    # rest of the warm-up. The value is recomputed (identically)
                    # in the worker if ever needed.
                    continue
