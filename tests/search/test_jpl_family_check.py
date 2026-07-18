"""Tests for the JPL SSD Three-Body Periodic Orbit catalog GATE (#647).

Two kinds of coverage:

1. **Golden tests** (sourced, no network): fixture payloads captured VERBATIM
   from live ``https://ssd-api.jpl.nasa.gov/periodic_orbits.api`` queries on
   2026-07-18 (see each fixture's own ``filter``/``system`` block for the
   exact request), fed through :func:`check_jpl_family` via an injected
   ``query_fn`` that loads the local fixture instead of the network. The
   EXPECTED jacobi/period/stability values below are copied verbatim from
   those captured payloads (JPL's own published catalog data), never from
   this codebase's own computed output, per
   ``[[feedback_golden_tests_sourced_only]]``.
2. **#641 reproduction**: the 5 real Sun-Jupiter periodic-orbit family
   clusters #641 found (see ``data/OUTSTANDING.md`` `#641`'s RESULT block)
   all fed the same spurious "published" verdict citing one cycler paper
   through ``literature_check.check_literature``'s cycler-keyword matcher.
   This module must instead return an HONEST ``"not-covered"`` for every one
   of them (Sun-Jupiter is not one of JPL's 7 indexed systems) — different,
   more honest output than the old spurious shared match.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from cyclerfinder.search.jpl_family_check import (
    DEFAULT_JACOBI_TOL,
    DEFAULT_PERIOD_REL_TOL,
    QueryFn,
    check_jpl_family,
    normalize_family,
    normalize_system,
)
from cyclerfinder.verify.jpl_periodic_orbits import JplOrbit, JplSystemConstants, parse_payload

_FIXTURES = Path(__file__).parent / "fixtures"


def _fixture_query_fn(fixture_name: str) -> QueryFn:
    """Build a fake ``jpo.query``-shaped callable that ignores its params and
    always returns the same locally-parsed fixture payload (no network)."""
    payload = json.loads((_FIXTURES / fixture_name).read_text())

    def _fn(system: str, family: str, **kwargs: Any) -> tuple[JplSystemConstants, list[JplOrbit]]:
        return parse_payload(payload)

    return _fn


def _empty_query_fn(
    system: str, family: str, **kwargs: Any
) -> tuple[JplSystemConstants, list[JplOrbit]]:
    """A fake query that ran cleanly but found nothing (a real 'no-match')."""
    payload = {
        "system": {
            "name": "Earth-Moon",
            "mass_ratio": "1.215058560962404e-02",
            "lunit": 389703.264829278,
            "tunit": 382981.289129055,
        },
        "fields": ["x", "y", "z", "vx", "vy", "vz", "jacobi", "period", "stability"],
        "data": [],
    }
    return parse_payload(payload)


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------


def test_normalize_system_and_family() -> None:
    assert normalize_system("Earth-Moon") == "earth-moon"
    assert normalize_system("Earth_Moon") == "earth-moon"
    assert normalize_system(" Sun Jupiter ") == "sun-jupiter"
    assert normalize_family("Halo") == "halo"
    assert normalize_family("DRO") == "dro"


# ---------------------------------------------------------------------------
# Golden tests: sourced JPL-catalogued members at 3 different supported
# systems. Each EXPECTED value below is copied verbatim from the fixture's
# own captured JPL response (see the fixture file's "filter"/"system" block).
# ---------------------------------------------------------------------------


def test_golden_earth_moon_halo_l1_north_matches() -> None:
    """Earth-Moon L1 halo (N branch), JPL jacobi=3.00000362096528,
    period=1.8077229598885376 TU, stability=2.59397548287644 -- captured live
    2026-07-18 via ``sys=earth-moon&family=halo&libr=1&branch=N&jacobimin=3.0
    &jacobimax=3.02`` (see ``fixtures/jpl_earth_moon_halo_l1n_c3.json``)."""
    result = check_jpl_family(
        "Earth-Moon",
        "halo",
        jacobi=3.00000362096528,
        period=1.8077229598885376,
        mu=0.01215058439469525,
        libr=1,
        branch="N",
        query_fn=_fixture_query_fn("jpl_earth_moon_halo_l1n_c3.json"),
    )
    assert result.status == "matched"
    assert result.system == "earth-moon"
    assert result.family == "halo"
    assert result.matched_jacobi == pytest.approx(3.00000362096528)
    assert result.matched_period == pytest.approx(1.8077229598885376)
    assert result.matched_stability == pytest.approx(2.59397548287644)
    assert result.jacobi_diff == pytest.approx(0.0, abs=1e-9)
    assert result.period_diff == pytest.approx(0.0, abs=1e-9)
    assert result.confidence == pytest.approx(1.0, abs=1e-6)
    assert result.jpl_mu == pytest.approx(1.215058560962404e-2)
    assert result.mu_reconciliation is not None
    assert result.citation  # a real citation string, not empty


def test_golden_saturn_titan_dro_matches() -> None:
    """Saturn-Titan DRO, JPL jacobi=2.99991417019945,
    period=2.2365274428311821 TU, stability=1.00000000079466 -- captured live
    2026-07-18 via ``sys=saturn-titan&family=dro&jacobimin=2.9999&
    jacobimax=3.0006`` (see ``fixtures/jpl_saturn_titan_dro_c3.json``)."""
    result = check_jpl_family(
        "saturn-titan",
        "dro",
        jacobi=2.99991417019945,
        period=2.2365274428311821,
        query_fn=_fixture_query_fn("jpl_saturn_titan_dro_c3.json"),
    )
    assert result.status == "matched"
    assert result.matched_jacobi == pytest.approx(2.99991417019945)
    assert result.matched_period == pytest.approx(2.2365274428311821)
    assert result.matched_stability == pytest.approx(1.00000000079466)
    assert result.confidence == pytest.approx(1.0, abs=1e-6)
    # Saturn-Titan mass ratio as JPL reports it (system constants block).
    assert result.jpl_mu == pytest.approx(2.366393158331484e-4)


def test_golden_mars_phobos_dro_matches() -> None:
    """Mars-Phobos DRO, JPL jacobi=2.99999004459787,
    period=4.7406333626886079 TU, stability=1.0 -- captured live 2026-07-18
    via ``sys=mars-phobos&family=dro&jacobimin=2.99999&jacobimax=3.00001``
    (see ``fixtures/jpl_mars_phobos_dro_c3.json``)."""
    result = check_jpl_family(
        "mars-phobos",
        "dro",
        jacobi=2.99999004459787,
        period=4.7406333626886079,
        query_fn=_fixture_query_fn("jpl_mars_phobos_dro_c3.json"),
    )
    assert result.status == "matched"
    assert result.matched_jacobi == pytest.approx(2.99999004459787)
    assert result.matched_period == pytest.approx(4.7406333626886079)
    assert result.matched_stability == pytest.approx(1.0)
    assert result.confidence == pytest.approx(1.0, abs=1e-6)
    assert result.jpl_mu == pytest.approx(1.611081404409632e-8)


def test_golden_picks_closest_of_multiple_candidates() -> None:
    """The fixture carries 2 rows; a jacobi/period aimed at the SECOND row
    must pick that one, not just always return data[0] -- proves the matcher
    actually compares invariants rather than returning an arbitrary member."""
    result = check_jpl_family(
        "earth-moon",
        "halo",
        jacobi=3.00002077318162,
        period=2.0285868576890547,
        libr=1,
        branch="N",
        query_fn=_fixture_query_fn("jpl_earth_moon_halo_l1n_c3.json"),
    )
    assert result.status == "matched"
    assert result.matched_jacobi == pytest.approx(3.00002077318162)
    assert result.matched_period == pytest.approx(2.0285868576890547)
    assert result.matched_stability == pytest.approx(1.52299191695338)


# ---------------------------------------------------------------------------
# Honest scope handling: unsupported system / family / no-match / errors.
# ---------------------------------------------------------------------------


def test_unsupported_system_returns_not_covered_without_network() -> None:
    def _boom(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("must not call the network/query_fn for an unsupported system")

    result = check_jpl_family("Sun-Jupiter", "halo", jacobi=3.0, period=2.0, query_fn=_boom)
    assert result.status == "not-covered"
    assert result.system == "sun-jupiter"
    assert "sun-jupiter" in result.notes.lower()
    assert result.n_candidates == 0
    assert result.confidence == 0.0


def test_unsupported_family_returns_not_covered_without_network() -> None:
    def _boom(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("must not call the network/query_fn for an unsupported family")

    result = check_jpl_family("earth-moon", "cycler", jacobi=3.0, period=2.0, query_fn=_boom)
    assert result.status == "not-covered"
    assert result.family == "cycler"


def test_family_requiring_libr_without_one_is_error_not_no_match() -> None:
    result = check_jpl_family("earth-moon", "halo", jacobi=3.0, period=2.0)
    assert result.status == "error"
    assert "libr" in result.notes.lower()


def test_family_requiring_branch_without_one_is_error_not_no_match() -> None:
    result = check_jpl_family("earth-moon", "halo", jacobi=3.0, period=2.0, libr=1)
    assert result.status == "error"
    assert "branch" in result.notes.lower()


def test_clean_query_with_no_candidates_is_no_match_not_error() -> None:
    result = check_jpl_family(
        "earth-moon",
        "dro",
        jacobi=1.5,
        period=99.0,
        query_fn=_empty_query_fn,
    )
    assert result.status == "no-match"
    assert result.n_candidates == 0
    assert result.jpl_mu is not None  # system constants still came back


def test_query_exception_is_error_not_crash() -> None:
    def _raises(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError("simulated network failure")

    result = check_jpl_family("earth-moon", "dro", jacobi=3.0, period=2.0, query_fn=_raises)
    assert result.status == "error"
    assert "simulated network failure" in result.notes


def test_nonphysical_inputs_are_error() -> None:
    result = check_jpl_family("earth-moon", "dro", jacobi=float("nan"), period=2.0)
    assert result.status == "error"
    result = check_jpl_family("earth-moon", "dro", jacobi=3.0, period=0.0)
    assert result.status == "error"


def test_defaults_are_sane_tolerances() -> None:
    assert 0 < DEFAULT_JACOBI_TOL < 1
    assert 0 < DEFAULT_PERIOD_REL_TOL < 1


def test_to_review_block_shape() -> None:
    result = check_jpl_family(
        "earth-moon",
        "halo",
        jacobi=3.00000362096528,
        period=1.8077229598885376,
        libr=1,
        branch="N",
        query_fn=_fixture_query_fn("jpl_earth_moon_halo_l1n_c3.json"),
    )
    block = result.to_review_block()
    assert block["status"] == "matched"
    assert block["system"] == "earth-moon"
    assert isinstance(block["query_params"], dict)


# ---------------------------------------------------------------------------
# #641 reproduction: the 5 real Sun-Jupiter clusters that ALL got the same
# spurious "published" verdict from literature_check.check_literature. This
# module must return an honest, DISTINCT-in-kind "not-covered" for each --
# not a fake shared match, not a crash, not a false "no-match".
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Sj641Cluster:
    jacobi: float
    period: float
    tag: str


# The 5 genuine periodic-orbit clusters #641 found at Sun-Jupiter
# mu=9.538811521016751e-4 (data/OUTSTANDING.md #641's RESULT block, verbatim
# jacobi/period/tag).
_SUN_JUPITER_641_CLUSTERS: tuple[_Sj641Cluster, ...] = (
    _Sj641Cluster(jacobi=3.0088, period=1.884, tag="halo"),
    _Sj641Cluster(jacobi=2.9896, period=3.594, tag="dro-like"),
    _Sj641Cluster(jacobi=3.0372, period=2.069, tag="lyapunov-like"),
    _Sj641Cluster(jacobi=3.0379, period=1.615, tag="lyapunov-like"),
    _Sj641Cluster(jacobi=3.4548, period=8.876, tag="far-family"),
)


def _boom(*args: Any, **kwargs: Any) -> Any:
    raise AssertionError("Sun-Jupiter is not JPL-covered; must never reach the network")


def test_641_sun_jupiter_clusters_all_honestly_not_covered() -> None:
    results = [
        check_jpl_family(
            "Sun-Jupiter",
            "halo" if c.tag == "halo" else "resonant",
            jacobi=c.jacobi,
            period=c.period,
            libr=1 if c.tag == "halo" else None,
            branch="N" if c.tag == "halo" else None,
            query_fn=_boom,
        )
        for c in _SUN_JUPITER_641_CLUSTERS
    ]
    assert len(results) == 5
    for r, c in zip(results, _SUN_JUPITER_641_CLUSTERS, strict=True):
        assert r.status == "not-covered", (c, r)
        assert r.system == "sun-jupiter"
        assert r.confidence == 0.0
        # Honest: no citation is fabricated for an out-of-scope system.
        assert r.citation == ""

    # The OLD literature_check.py path returned the SAME single citation
    # ("Strange-Russell 2007 AAS 07-277") for all 5 physically distinct
    # clusters -- a spurious identical "published" match. This path instead
    # gives every cluster the SAME status ("not-covered", correctly, since
    # they share one out-of-scope system) but crucially never asserts a
    # false positive citation/match for any of them -- the core defect this
    # task fixes (a real numeric verdict engine, not a keyword collision)
    # is that none of the 5 gets a fabricated "matched" family identity.
    assert all(r.status == "not-covered" for r in results)
    assert all(r.matched_jacobi is None for r in results)
    assert all(r.matched_period is None for r in results)
    assert all(not r.citation for r in results)
