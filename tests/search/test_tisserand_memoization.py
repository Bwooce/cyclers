"""Regression-safety tests for the #472 Tisserand memoization.

:func:`cyclerfinder.search.tisserand.linkable`, :func:`~.vinf_to_tisserand`,
:func:`~.tisserand_to_vinf` and :func:`~._a_p_km` are memoized with
``lru_cache``. All arguments are hashable (str / float / tuple) and the
functions read only the *frozen* ``PLANETS`` / ``SATELLITES`` body tables plus a
scalar ``mu``, so the cache is exact. Guards parity, key-correctness, purity.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest

import cyclerfinder.search.tisserand as tisserand
from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.satellites import PRIMARIES


@pytest.fixture(autouse=True)
def _clear() -> Iterator[None]:
    for fn in (
        tisserand._a_p_km,
        tisserand.vinf_to_tisserand,
        tisserand.tisserand_to_vinf,
        tisserand.linkable,
    ):
        fn.cache_clear()
    yield


def test_parity_conversions() -> None:
    for body in ("E", "M", "V", "J"):
        for v in (0.0, 2.0, 5.0, 9.0):
            assert tisserand.vinf_to_tisserand(body, v) == (
                tisserand.vinf_to_tisserand.__wrapped__(body, v)
            )
        for t in (2.0, 2.5, 3.5):
            assert tisserand.tisserand_to_vinf(body, t) == (
                tisserand.tisserand_to_vinf.__wrapped__(body, t)
            )


def test_parity_linkable_heliocentric_and_jovicentric() -> None:
    # heliocentric planet pairs
    for a, b, v in [("E", "M", 3.0), ("E", "V", 5.0), ("V", "M", 7.0)]:
        assert tisserand.linkable(a, b, v) == tisserand.linkable.__wrapped__(a, b, v)
    # Jovicentric moon pair with the same (str, float, tuple, mu) keys the
    # moon-prune gate feeds.
    mu = PRIMARIES["Jupiter"]
    arange = (5.0e-3, 2.0e-2)
    cached = tisserand.linkable("Io", "Europa", 4.0, a_range_au=arange, mu=mu)
    uncached = tisserand.linkable.__wrapped__("Io", "Europa", 4.0, a_range_au=arange, mu=mu)
    assert cached == uncached


def test_key_correctness_distinct_vinf_and_mu() -> None:
    # Different V∞ must not alias.
    t1 = tisserand.vinf_to_tisserand("E", 3.0)
    t2 = tisserand.vinf_to_tisserand("E", 7.0)
    assert t1 != t2
    assert tisserand.vinf_to_tisserand.cache_info().currsize == 2
    # mu is part of the key: heliocentric vs a primary GM give distinct entries.
    tisserand.vinf_to_tisserand.cache_clear()
    helio = tisserand.vinf_to_tisserand("E", 5.0)
    jovi = tisserand.vinf_to_tisserand("E", 5.0, mu=PRIMARIES["Jupiter"])
    assert helio != jovi
    assert tisserand.vinf_to_tisserand.cache_info().currsize == 2


def test_key_correctness_linkable_repeat_is_a_hit() -> None:
    tisserand.linkable.cache_clear()
    first = tisserand.linkable("E", "M", 3.0)
    assert tisserand.linkable.cache_info().hits == 0
    second = tisserand.linkable("E", "M", 3.0)
    assert second == first
    assert tisserand.linkable.cache_info().hits == 1


def test_key_correctness_a_range_tuple_keys_distinctly() -> None:
    tisserand.linkable.cache_clear()
    tisserand.linkable("E", "M", 4.0, a_range_au=(0.3, 5.0))
    tisserand.linkable("E", "M", 4.0, a_range_au=(0.5, 3.0))
    # two different a_range tuples -> two distinct cache entries (no tuple
    # aliasing).
    assert tisserand.linkable.cache_info().currsize == 2


def test_purity_body_tables_immutable() -> None:
    from dataclasses import FrozenInstanceError

    with pytest.raises(FrozenInstanceError):
        PLANETS["E"].sma_au = 0.0  # type: ignore[misc]
