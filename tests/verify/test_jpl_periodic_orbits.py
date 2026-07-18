"""JPL periodic_orbits client — parser tests (no network) + mu reconciliation.

The parser is exercised against a saved 4-row Earth-Moon L1-halo fixture
captured from the live API (tests/verify/fixtures/). The live ``query`` is NOT
tested here for the actual network fetch (no network in CI); its CACHING
layer (#647) and filter-parameter construction ARE tested by stubbing
``urllib.request.urlopen`` so a cache-hit path proves it never touches the
network.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.verify.jpl_periodic_orbits import (
    FAMILIES_REQUIRING_BRANCH,
    FAMILIES_REQUIRING_LIBR,
    OUR_EM_MU,
    SUPPORTED_FAMILIES,
    SUPPORTED_SYSTEMS,
    _cache_path,
    parse_payload,
    query,
    reconcile_mu,
)

_FIXTURE = Path(__file__).parent / "fixtures" / "jpl_earth_moon_l1_halo_sample.json"


@pytest.fixture
def payload() -> dict[str, Any]:
    data: dict[str, Any] = json.loads(_FIXTURE.read_text())
    return data


def test_parse_system_constants(payload: dict[str, Any]) -> None:
    constants, _ = parse_payload(payload)
    assert constants.name == "Earth-Moon"
    # JPL's reported Earth-Moon mass ratio (the value to reconcile against ours).
    assert constants.mu == pytest.approx(1.215058560962404e-2, rel=0, abs=1e-15)
    assert constants.lunit_km > 0
    assert constants.tunit_s > 0
    # L1/L2/L3 lie on the x-axis (y=z=0); L1 is between the primaries (0<x<1).
    assert "L1" in constants.libration_points
    l1 = constants.libration_points["L1"]
    assert 0.0 < l1[0] < 1.0
    assert l1[1] == 0.0 and l1[2] == 0.0


def test_parse_orbits_shape_and_fields(payload: dict[str, Any]) -> None:
    _, orbits = parse_payload(payload)
    assert len(orbits) == 4
    for o in orbits:
        assert o.state0.shape == (6,)
        assert o.state0.dtype == np.float64
        assert np.isfinite(o.state0).all()
        assert o.period > 0
        assert np.isfinite(o.jacobi)
        assert np.isfinite(o.stability)


def test_parse_is_field_name_keyed_not_positional(payload: dict[str, Any]) -> None:
    """A column re-order in the payload must not silently mis-map values.

    Shuffle ``fields`` (and each data row in lockstep) and confirm the parser
    recovers the SAME orbits — proving it keys by name, not position.
    """
    _, orbits_ref = parse_payload(payload)

    fields = list(payload["fields"])
    perm = [fields.index(c) for c in reversed(fields)]
    shuffled = {
        **payload,
        "fields": [fields[i] for i in perm],
        "data": [[row[i] for i in perm] for row in payload["data"]],
    }
    _, orbits_shuf = parse_payload(shuffled)

    assert len(orbits_shuf) == len(orbits_ref)
    for a, b in zip(orbits_ref, orbits_shuf, strict=True):
        np.testing.assert_array_equal(a.state0, b.state0)
        assert a.jacobi == b.jacobi
        assert a.period == b.period
        assert a.stability == b.stability


def test_parse_raises_on_missing_field(payload: dict[str, Any]) -> None:
    broken = {**payload, "fields": [f for f in payload["fields"] if f != "jacobi"]}
    with pytest.raises(ValueError, match="missing expected fields"):
        parse_payload(broken)


def test_reconcile_mu_reports_the_gap() -> None:
    # The convention gap that gates any JPL-vs-ours cross-check: ~1e-7 relative.
    rec = reconcile_mu(1.215058560962404e-2)
    assert rec["our_mu"] == OUR_EM_MU
    assert rec["abs_diff"] == pytest.approx(1.215e-9, abs=5e-10)
    assert 1e-8 < rec["rel_diff"] < 1e-6


# ---------------------------------------------------------------------------
# Supported systems/families (#647) -- confirmed live 2026-07-18 via a direct
# 400-error parameter listing (curl https://ssd-api.jpl.nasa.gov/periodic_
# orbits.api?sys=bogus) and the API's own documentation page.
# ---------------------------------------------------------------------------


def test_supported_systems_is_the_confirmed_seven() -> None:
    assert {
        "sun-earth",
        "earth-moon",
        "sun-mars",
        "jupiter-europa",
        "saturn-enceladus",
        "saturn-titan",
        "mars-phobos",
    } == SUPPORTED_SYSTEMS
    assert "sun-jupiter" not in SUPPORTED_SYSTEMS  # the #641 correction


def test_supported_families_is_the_documented_twelve() -> None:
    assert {
        "halo",
        "vertical",
        "axial",
        "lyapunov",
        "longp",
        "short",
        "butterfly",
        "dragonfly",
        "resonant",
        "dro",
        "dpo",
        "lpo",
    } == SUPPORTED_FAMILIES


def test_families_requiring_libr_and_branch() -> None:
    assert FAMILIES_REQUIRING_LIBR <= SUPPORTED_FAMILIES
    assert FAMILIES_REQUIRING_BRANCH <= FAMILIES_REQUIRING_LIBR
    assert "halo" in FAMILIES_REQUIRING_BRANCH


# ---------------------------------------------------------------------------
# Caching layer (#647): a cache-hit must never touch the network. Stub
# urllib.request.urlopen to raise if called, prove the cache path alone
# satisfies the request.
# ---------------------------------------------------------------------------


def test_cache_path_is_deterministic_and_param_sensitive(tmp_path: Path) -> None:
    p1 = _cache_path(tmp_path, {"sys": "earth-moon", "family": "halo", "jacobimin": "3.0"})
    p2 = _cache_path(tmp_path, {"sys": "earth-moon", "family": "halo", "jacobimin": "3.0"})
    p3 = _cache_path(tmp_path, {"sys": "earth-moon", "family": "halo", "jacobimin": "3.1"})
    assert p1 == p2
    assert p1 != p3
    assert p1.parent == tmp_path


def test_query_cache_hit_never_touches_the_network(
    tmp_path: Path, payload: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    params = {"sys": "earth-moon", "family": "halo", "libr": "1", "branch": "N"}
    cache_file = _cache_path(tmp_path, params)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(payload))

    def _boom(*args: object, **kwargs: object) -> object:
        raise AssertionError("cache hit must not call urlopen")

    monkeypatch.setattr("urllib.request.urlopen", _boom)

    constants, orbits = query("earth-moon", "halo", libr=1, branch="N", cache_dir=tmp_path)
    assert constants.name == "Earth-Moon"
    assert len(orbits) == len(payload["data"])


def test_query_cache_miss_writes_the_cache_file(
    tmp_path: Path, payload: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:

    class _FakeResponse:
        def __enter__(self) -> _FakeResponse:
            return self

        def __exit__(self, *exc: object) -> None:
            return None

        def read(self) -> bytes:
            return json.dumps(payload).encode("utf-8")

    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: _FakeResponse())

    constants, orbits = query("earth-moon", "halo", libr=1, branch="N", cache_dir=tmp_path)
    assert constants.name == "Earth-Moon"
    assert len(orbits) == len(payload["data"])

    cached_files = list(tmp_path.glob("*.json"))
    assert len(cached_files) == 1
    on_disk = json.loads(cached_files[0].read_text())
    assert on_disk["system"]["name"] == "Earth-Moon"

    # A second call with an identical cache_dir must now hit the cache, not
    # urlopen (proves the write-then-read round trip actually works).
    monkeypatch.setattr(
        "urllib.request.urlopen",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("must use the cache")),
    )
    constants2, orbits2 = query("earth-moon", "halo", libr=1, branch="N", cache_dir=tmp_path)
    assert constants2.name == constants.name
    assert len(orbits2) == len(orbits)


def test_query_without_cache_dir_preserves_prior_always_live_behaviour(
    payload: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    """cache_dir=None (the default) must still hit urlopen every time -- the
    historical behaviour every pre-existing call site relies on."""
    calls = {"n": 0}

    class _FakeResponse:
        def __enter__(self) -> _FakeResponse:
            return self

        def __exit__(self, *exc: object) -> None:
            return None

        def read(self) -> bytes:
            calls["n"] += 1
            return json.dumps(payload).encode("utf-8")

    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: _FakeResponse())
    query("earth-moon", "halo", libr=1, branch="N")
    query("earth-moon", "halo", libr=1, branch="N")
    assert calls["n"] == 2  # no caching => the (fake) network was hit twice


def test_query_forwards_filter_params_into_the_url(
    payload: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, str] = {}

    class _FakeResponse:
        def __enter__(self) -> _FakeResponse:
            return self

        def __exit__(self, *exc: object) -> None:
            return None

        def read(self) -> bytes:
            return json.dumps(payload).encode("utf-8")

    def _fake_urlopen(req: Any, timeout: float = 60.0) -> _FakeResponse:
        import urllib.parse as up

        parsed = up.urlparse(req.full_url)
        captured.update(dict(up.parse_qsl(parsed.query)))
        return _FakeResponse()

    monkeypatch.setattr("urllib.request.urlopen", _fake_urlopen)
    query(
        "saturn-titan",
        "dro",
        jacobimin=2.9,
        jacobimax=3.1,
        periodmin=1.0,
        periodmax=5.0,
        periodunits="TU",
        stabmin=0.5,
        stabmax=100.0,
    )
    assert captured["sys"] == "saturn-titan"
    assert captured["family"] == "dro"
    assert captured["jacobimin"] == repr(2.9)
    assert captured["jacobimax"] == repr(3.1)
    assert captured["periodmin"] == repr(1.0)
    assert captured["periodmax"] == repr(5.0)
    assert captured["periodunits"] == "TU"
    assert captured["stabmin"] == repr(0.5)
    assert captured["stabmax"] == repr(100.0)
