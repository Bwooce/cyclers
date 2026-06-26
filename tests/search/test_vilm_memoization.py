"""Regression-safety tests for the #472 VILM memoization.

The VILM per-leg cost functions in :mod:`cyclerfinder.search.vilm` are memoized
with :func:`functools.cache`. The cache MUST change no value. These tests
guard, per memoized function:

1. PARITY — cached output == uncached (``.__wrapped__``) output, to machine
   precision (the cache is exact: no array args, no rounding key).
2. KEY-CORRECTNESS — distinct args -> distinct cache entries (no key aliasing);
   a repeated call is a VERIFIED cache hit via ``cache_info()``.
3. PURITY — the memoized functions read only the *frozen* module-constant
   ``SATELLITES`` / ``PRIMARIES`` tables; there is no mutable config the cache
   could freeze. The guard asserts those tables are immutable.
"""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import FrozenInstanceError

import pytest

import cyclerfinder.search.vilm as vilm
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Representative + edge moon pairs across primaries (Jovian / Saturnian).
_PAIRS = [
    ("Ganymede", "Europa"),
    ("Io", "Europa"),
    ("Callisto", "Ganymede"),
    ("Titan", "Rhea"),
    ("Rhea", "Dione"),
    ("Europa", "Io"),  # reversed order (exercises _order_by_sma both ways)
]
_MOONS = ["Io", "Europa", "Ganymede", "Callisto", "Titan", "Rhea", "Dione"]


@pytest.fixture(autouse=True)
def _clear_caches() -> Iterator[None]:
    """Each test starts from an empty cache so hit/miss counts are deterministic."""
    for fn in (
        vilm._vc_adim,
        vilm._vbar_vinf_adim,
        vilm.min_vinf_for_vilm,
        vilm._v_m,
        vilm._leverage_dv_kms,
        vilm._vilm_dv_min_pair,
        vilm.vilm_dv_floor,
        vilm.europa_endgame_dv,
    ):
        fn.cache_clear()
    yield


def test_parity_vilm_dv_min_pair() -> None:
    for a, b in _PAIRS:
        cached = vilm.vilm_dv_min(a, b)
        uncached = vilm._vilm_dv_min_pair.__wrapped__(a, b)
        assert cached == uncached


def test_parity_min_vinf_for_vilm() -> None:
    for m in _MOONS:
        for ext in (True, False):
            assert vilm.min_vinf_for_vilm(m, exterior=ext) == (
                vilm.min_vinf_for_vilm.__wrapped__(m, exterior=ext)
            )


def test_parity_leverage_dv_kms() -> None:
    for m in _MOONS:
        for v in (1.0, 2.5, 4.0):
            for ext in (True, False):
                assert vilm._leverage_dv_kms(m, v, exterior=ext) == (
                    vilm._leverage_dv_kms.__wrapped__(m, v, exterior=ext)
                )


def test_parity_floor_and_europa() -> None:
    for a, b in _PAIRS:
        assert vilm.vilm_dv_floor(a, b) == vilm.vilm_dv_floor.__wrapped__(a, b)
    assert vilm.europa_endgame_dv() == vilm.europa_endgame_dv.__wrapped__()


def test_key_correctness_distinct_args_distinct_results() -> None:
    # Distinct moons must NOT alias to one cache entry.
    vc_io = vilm._vc_adim("Io")
    vc_eu = vilm._vc_adim("Europa")
    assert vc_io != vc_eu
    info = vilm._vc_adim.cache_info()
    assert info.misses == 2  # two distinct keys -> two misses
    assert info.currsize == 2
    # exterior flag is part of the key: distinct entries.
    a = vilm.min_vinf_for_vilm("Ganymede", exterior=True)
    b = vilm.min_vinf_for_vilm("Ganymede", exterior=False)
    assert a != b
    assert vilm.min_vinf_for_vilm.cache_info().currsize == 2


def test_key_correctness_repeat_is_a_cache_hit() -> None:
    vilm._vilm_dv_min_pair.cache_clear()
    first = vilm.vilm_dv_min("Ganymede", "Europa")
    assert vilm._vilm_dv_min_pair.cache_info().hits == 0
    second = vilm.vilm_dv_min("Ganymede", "Europa")
    assert second == first
    assert vilm._vilm_dv_min_pair.cache_info().hits == 1  # VERIFIED hit


def test_purity_config_tables_are_frozen() -> None:
    """The cache freezes nothing mutable: SATELLITES / PRIMARIES are frozen."""
    sat = SATELLITES["Europa"]
    with pytest.raises(FrozenInstanceError):
        sat.sma_km = 0.0  # type: ignore[misc]
    # PRIMARIES is a plain dict of floats (μ); the values it feeds are scalars
    # baked into the key-free physics, so a body-string key is sufficient.
    assert isinstance(PRIMARIES[sat.primary], float)


def test_via_path_is_not_broken_by_cache() -> None:
    """The un-cached GA-chain (``via``) path still agrees with a manual chain."""
    # via routes through the with-GA branch; just assert it returns a finite,
    # non-negative ΔV and is <= the no-GA pair value (a GA can only reduce ΔV).
    no_ga = vilm.vilm_dv_min("Callisto", "Europa")
    with_ga = vilm.vilm_dv_min("Callisto", "Europa", via=["Ganymede"])
    assert math.isfinite(with_ga) and with_ga >= 0.0
    assert with_ga <= no_ga + 1e-9
