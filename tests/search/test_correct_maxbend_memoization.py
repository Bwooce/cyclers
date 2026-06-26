"""Regression-safety tests for the #472 ``_max_bend_deg`` memoization.

:func:`cyclerfinder.search.correct._max_bend_deg` now delegates its nominal
(``rp_factors=None``) path to the memoized
:func:`~cyclerfinder.search.correct._max_bend_deg_nominal`. The ``rp_factors``
(mutable dict) path is INTENTIONALLY left un-cached — its dict argument is
unhashable and, more importantly, caching it would freeze a per-call config
override the cache must never alias. These tests pin that behaviour.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import pytest

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.search.correct import _max_bend_deg, _max_bend_deg_nominal


@pytest.fixture(autouse=True)
def _clear() -> Iterator[None]:
    _max_bend_deg_nominal.cache_clear()
    yield


def _manual(vinf_kms: float, body: str, factor: float = 1.0) -> float:
    pl = PLANETS[body]
    r_p = pl.radius_eq_km + pl.safe_alt_km * factor
    e = 1.0 + r_p * vinf_kms * vinf_kms / pl.mu_km3_s2
    return float(np.degrees(2.0 * np.arcsin(1.0 / e)))


def test_parity_nominal_path() -> None:
    for body in ("E", "M", "V"):
        for v in (1.0, 3.5, 6.9):
            cached = _max_bend_deg(v, body)
            assert cached == _max_bend_deg_nominal.__wrapped__(v, body)
            assert cached == pytest.approx(_manual(v, body), rel=0, abs=1e-12)


def test_rp_factors_path_is_correct_and_uncached() -> None:
    """A scaled rp_factors must change the result AND not poison the cache."""
    nominal = _max_bend_deg(5.0, "E")
    scaled = _max_bend_deg(5.0, "E", {"E": 0.5})
    assert scaled != nominal
    assert scaled == pytest.approx(_manual(5.0, "E", 0.5), rel=0, abs=1e-12)
    # The scaled (config-override) call must NOT have written a cache entry that
    # could later be returned for the nominal key.
    _max_bend_deg_nominal.cache_clear()
    again_scaled = _max_bend_deg(5.0, "E", {"E": 0.5})
    assert again_scaled == scaled
    assert _max_bend_deg_nominal.cache_info().currsize == 0  # scaled path never cached
    # nominal still correct after the scaled calls
    assert _max_bend_deg(5.0, "E") == nominal


def test_rp_factors_absent_body_uses_nominal() -> None:
    """rp_factors that does not name ``body`` must equal the nominal value."""
    nominal = _max_bend_deg(4.0, "E")
    assert _max_bend_deg(4.0, "E", {"M": 0.5}) == nominal


def test_key_correctness_distinct_and_hit() -> None:
    a = _max_bend_deg(3.5, "E")
    b = _max_bend_deg(6.9, "E")
    assert a != b
    assert _max_bend_deg_nominal.cache_info().misses == 2
    assert _max_bend_deg(3.5, "E") == a
    assert _max_bend_deg_nominal.cache_info().hits == 1


def test_config_change_would_not_be_stale() -> None:
    """Purity guard: the nominal cache is keyed on (vinf, body) only because the
    body record is frozen config. If a caller needs a DIFFERENT periapsis it
    must use rp_factors (the un-cached path) — proven by the scaled value
    differing from the cached nominal one for the same (vinf, body) key."""
    nominal = _max_bend_deg(5.0, "E")  # cached under ("E", 5.0)
    scaled = _max_bend_deg(5.0, "E", {"E": 2.0})  # same key, different config
    assert scaled != nominal  # the override is honoured, not served from cache
