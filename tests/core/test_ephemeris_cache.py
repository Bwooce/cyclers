"""State-cache + vectorised batch tests for the Ephemeris layer (Task #128-S2).

These guard the two performance levers added in
``docs/notes/2026-06-07-ephemeris-perf-s2.md``:

1. Per-instance ``state()`` memoisation (default ON) — must be byte-identical to
   the uncached path and must NEVER share entries across differently-
   parameterised instances (circular / astropy / ramped-continuation backends).
2. The vectorised ``states(bodies, epochs)`` batch API — per-element
   byte-identical to repeated ``state()`` calls.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris

# A replay fixture: a fixed list of (body, epoch) queries with deliberate exact
# duplicates, so the cache hit AND miss paths are both exercised on replay.
_REPLAY: list[tuple[str, float]] = [
    ("E", 0.0),
    ("M", 1.0e7),
    ("E", 2.5e7),
    ("M", 1.0e7),  # dup of #2
    ("E", 0.0),  # dup of #1
    ("V", -3.0e7),
    ("E", 2.5e7),  # dup of #3
]


@pytest.mark.parametrize("model", ["circular", "inclined-circular", "astropy"])
def test_cache_is_byte_identical_to_uncached(model: str) -> None:
    """Cache ON must return byte-identical floats to cache OFF, on miss AND hit."""
    cached = Ephemeris(model, cache=True)
    plain = Ephemeris(model, cache=False)
    # Two passes so the second pass hits the cache for every entry.
    for _ in range(2):
        for body, t in _REPLAY:
            r_c, v_c = cached.state(body, t)
            r_p, v_p = plain.state(body, t)
            assert np.array_equal(r_c, r_p), (model, body, t)
            assert np.array_equal(v_c, v_p), (model, body, t)


def test_cache_returns_independent_writeable_copies() -> None:
    """Mutating a returned array must not corrupt the cached value."""
    e = Ephemeris("astropy", cache=True)
    r, v = e.state("E", 0.0)
    assert r.flags.writeable and v.flags.writeable
    r[0] += 1234.5
    v[1] -= 9.9
    r2, v2 = e.state("E", 0.0)
    assert r2[0] != r[0]
    assert v2[1] != v[1]


def test_cache_disabled_path_matches() -> None:
    e_on = Ephemeris("astropy", cache=True)
    e_off = Ephemeris("astropy", cache=False)
    r1, v1 = e_on.state("M", 7.0e7)
    r2, v2 = e_off.state("M", 7.0e7)
    assert np.array_equal(r1, r2)
    assert np.array_equal(v1, v2)


def test_cache_lru_eviction_bounds_memory() -> None:
    """LRU caps at cache_size; eviction does not change returned values."""
    e = Ephemeris("circular", cache=True, cache_size=3)
    plain = Ephemeris("circular", cache=False)
    epochs = [float(i) * 1.0e6 for i in range(10)]
    for t in epochs:
        e.state("E", t)
    # Cache holds at most cache_size entries after the sweep.
    assert len(e._state_cache) <= 3
    # Re-query an evicted early epoch: still byte-identical (recomputed).
    r, v = e.state("E", epochs[0])
    rp, vp = plain.state("E", epochs[0])
    assert np.array_equal(r, rp) and np.array_equal(v, vp)


def test_clear_cache_empties_it() -> None:
    e = Ephemeris("astropy", cache=True)
    e.state("E", 0.0)
    assert len(e._state_cache) == 1
    e.clear_cache()
    assert len(e._state_cache) == 0


# ---------------------------------------------------------------------------
# Cache isolation: two differently-parameterised instances must NOT share
# entries (PROOF the cache key is per-instance, never global).
# ---------------------------------------------------------------------------


def test_circular_and_astropy_caches_do_not_cross_contaminate() -> None:
    circ = Ephemeris("circular", cache=True)
    astro = Ephemeris("astropy", cache=True)
    body, t = "M", 5.0e7
    rc, _ = circ.state(body, t)
    ra, _ = astro.state(body, t)
    # Each instance serves its OWN cached value on replay.
    rc2, _ = circ.state(body, t)
    ra2, _ = astro.state(body, t)
    assert np.array_equal(rc, rc2)
    assert np.array_equal(ra, ra2)
    # The two backends genuinely disagree at this (body, epoch); if the cache
    # leaked one into the other this assertion would fail.
    assert not np.allclose(rc, ra)


def test_ramped_continuation_backends_do_not_cross_contaminate() -> None:
    """continuation.ramped_ephemeris swaps _backend post-construction; two
    differently-ramped instances must each cache their own states."""
    from cyclerfinder.search.continuation import (
        _RampedElementsBackend,
        ramped_ephemeris,
    )

    body, t = "M", 5.0e7
    r0 = ramped_ephemeris(0.0, 0.0, 0.0)  # bit-identical to circular
    r1 = ramped_ephemeris(1.0, 1.0, 1.0)  # full e+i+phase ellipse
    a0, _ = r0.state(body, t)
    a1, _ = r1.state(body, t)
    # Replay: each serves its own cached value.
    assert np.array_equal(a0, r0.state(body, t)[0])
    assert np.array_equal(a1, r1.state(body, t)[0])
    assert not np.allclose(a0, a1)
    # lam=0 ramped is byte-identical to a fresh circular (and the circular cache
    # never leaked into the ramped instance).
    assert np.array_equal(a0, Ephemeris("circular", cache=True).state(body, t)[0])

    # Post-construction backend swap on the SAME instance invalidates the cache.
    e = Ephemeris("circular", cache=True)
    pre, _ = e.state(body, t)
    e._backend = _RampedElementsBackend(1.0, 1.0, 1.0)
    post, _ = e.state(body, t)
    assert not np.allclose(pre, post)


# ---------------------------------------------------------------------------
# Vectorised batch API
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", ["circular", "inclined-circular", "astropy"])
def test_states_batch_matches_scalar_state(model: str) -> None:
    batch = Ephemeris(model, cache=False)
    scalar = Ephemeris(model, cache=False)
    bodies = [b for b, _ in _REPLAY]
    epochs = [t for _, t in _REPLAY]
    out = batch.states(bodies, epochs)
    assert len(out) == len(_REPLAY)
    for (rb, vb), (body, t) in zip(out, _REPLAY, strict=True):
        rs, vs = scalar.state(body, t)
        assert np.array_equal(rb, rs), (model, body, t)
        assert np.array_equal(vb, vs), (model, body, t)


def test_states_single_element_shape() -> None:
    e = Ephemeris("astropy", cache=False)
    out = e.states(["E"], [3.3e7])
    assert len(out) == 1
    r, v = out[0]
    assert r.shape == (3,) and v.shape == (3,)
    rs, vs = e.state("E", 3.3e7)
    assert np.array_equal(r, rs) and np.array_equal(v, vs)


def test_states_uses_cache_for_hits() -> None:
    e = Ephemeris("astropy", cache=True)
    # Prime the cache with one entry.
    r_primed, _ = e.state("E", 0.0)
    # Batch including the primed entry + a fresh one.
    out = e.states(["E", "M"], [0.0, 1.0e7])
    assert np.array_equal(out[0][0], r_primed)
    # Cache now holds both.
    assert ("E", 0.0) in e._state_cache
    assert ("M", 1.0e7) in e._state_cache


def test_states_length_mismatch_raises() -> None:
    e = Ephemeris("circular")
    with pytest.raises(ValueError, match="length mismatch"):
        e.states(["E", "M"], [0.0])
