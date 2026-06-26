"""Regression-safety tests for the #472 flyby scalar-cost memoization.

:func:`cyclerfinder.core.flyby.max_bend`, :func:`~.dv_from_turn_deficit` and
:func:`~.dv_powered_flyby_periapsis` are memoized with ``lru_cache``. They take
only float arguments (no ndarray keys, no module-constant reads), so the cache
is exact. These tests guard parity, key-correctness, and that the vector
functions were deliberately left un-cached.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest

import cyclerfinder.core.flyby as flyby

# Earth-like μ and a safe periapsis; representative + edge V∞.
_MU = 3.986004418e5
_RP = 6578.0
_VINFS = [0.0, 1.0, 3.5, 6.9, 12.0]


@pytest.fixture(autouse=True)
def _clear() -> Iterator[None]:
    flyby.max_bend.cache_clear()
    flyby.dv_from_turn_deficit.cache_clear()
    flyby.dv_powered_flyby_periapsis.cache_clear()
    yield


def test_parity_max_bend() -> None:
    for v in _VINFS:
        assert flyby.max_bend(_MU, _RP, v) == flyby.max_bend.__wrapped__(_MU, _RP, v)


def test_parity_dv_from_turn_deficit() -> None:
    for v in _VINFS:
        for dreq, dmax in [(0.5, 1.0), (1.2, 0.8), (3.0, 0.4)]:
            assert flyby.dv_from_turn_deficit(v, dreq, dmax) == (
                flyby.dv_from_turn_deficit.__wrapped__(v, dreq, dmax)
            )


def test_parity_dv_powered_flyby_periapsis() -> None:
    for v in _VINFS:
        for dreq, dmax in [(0.5, 1.0), (1.2, 0.8), (3.2, 0.4)]:
            assert flyby.dv_powered_flyby_periapsis(v, dreq, dmax, _MU, _RP) == (
                flyby.dv_powered_flyby_periapsis.__wrapped__(v, dreq, dmax, _MU, _RP)
            )


def test_key_correctness_distinct_and_hit() -> None:
    a = flyby.max_bend(_MU, _RP, 3.5)
    b = flyby.max_bend(_MU, _RP, 6.9)
    assert a != b  # distinct V∞ -> distinct results, no aliasing
    assert flyby.max_bend.cache_info().misses == 2
    # repeat is a verified hit
    again = flyby.max_bend(_MU, _RP, 3.5)
    assert again == a
    assert flyby.max_bend.cache_info().hits == 1


def test_error_inputs_not_cached_as_values() -> None:
    """A ValueError input re-raises every call (functools.cache does not cache it)."""
    with pytest.raises(ValueError):
        flyby.max_bend(_MU, _RP, -1.0)
    with pytest.raises(ValueError):
        flyby.max_bend(_MU, _RP, -1.0)
    # no value cached for the bad key
    assert flyby.max_bend.cache_info().currsize == 0


def test_vector_flyby_dv_is_not_cached() -> None:
    """Vector entrypoints keep their ndarray contract — no lru_cache wrapper."""
    assert not hasattr(flyby.flyby_dv, "cache_info")
    assert not hasattr(flyby.bend_angle, "cache_info")
    assert not hasattr(flyby.is_ballistic_feasible, "cache_info")
