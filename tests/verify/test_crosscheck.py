"""M7 V1 Lambert cross-check tests — spec §14 V1, plan §3.4.

Covers :mod:`cyclerfinder.verify.crosscheck`:

* Gate — :func:`test_v1_lambert_crosscheck_aldrin` (M7 binding gate):
  every single-rev leg of a real-ephemeris Aldrin cycler agrees with
  both lamberthub solvers to < :data:`V1_TOLERANCE_MPS`.
* Helper-level tests for :func:`crosscheck_leg` /
  :func:`crosscheck_cycler` behaviour (multi-rev routing, endpoint
  extraction, aggregation).
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler, Leg
from cyclerfinder.search.phase_match import phase_signature_from_catalogue_entry
from cyclerfinder.verify.crosscheck import (
    V1_TOLERANCE_MPS,
    LambertCrosscheckResult,
    crosscheck_cycler,
    crosscheck_leg,
)
from cyclerfinder.verify.real_closure import (
    _resolve_real_t_start,
    construct_real_ephemeris_cycler,
)
from tests.data._catalogue_loader_m6b import load_m6b_entries

ALDRIN_PRIORITY = datetime(1985, 10, 28, tzinfo=UTC)


@pytest.fixture(scope="module")
def aldrin_entry() -> dict[str, object]:
    """Return the loaded Aldrin outbound catalogue entry."""
    entries = load_m6b_entries()
    return next(e for e in entries if e["id"] == "aldrin-classic-em-k1-outbound")


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


@pytest.fixture(scope="module")
def aldrin_cycler(
    aldrin_entry: dict[str, object],
    astropy_ephem: Ephemeris,
) -> Cycler:
    """Real-ephemeris Aldrin cycler (single-rev E->M->E Lambert chain)."""
    signature = phase_signature_from_catalogue_entry(aldrin_entry)
    t_start = _resolve_real_t_start(signature, astropy_ephem, ALDRIN_PRIORITY)
    assert t_start is not None
    return construct_real_ephemeris_cycler(aldrin_entry, astropy_ephem, t_start)


# ---------------------------------------------------------------------------
# M7 BINDING GATE — V1 lamberthub agreement on a real Aldrin cycler
# ---------------------------------------------------------------------------


def test_v1_lambert_crosscheck_aldrin(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """M7 GATE — spec §14 V1, plan §3.4.

    Every single-rev leg of the real-ephemeris Aldrin cycler re-solves
    to within :data:`V1_TOLERANCE_MPS` between the in-house Lambert and
    both lamberthub solvers.
    """
    results = crosscheck_cycler(aldrin_cycler, astropy_ephem)

    assert results, "Aldrin cycler has no single-rev legs to cross-check"
    assert len(results) == sum(1 for leg in aldrin_cycler.legs if leg.n_revs == 0)
    for r in results:
        assert isinstance(r, LambertCrosscheckResult)
        assert r.passed is True, (
            f"leg {r.leg_index} disagreement {r.max_diff_mps} m/s >= {V1_TOLERANCE_MPS} m/s"
        )
        assert r.max_diff_mps < V1_TOLERANCE_MPS


# ---------------------------------------------------------------------------
# crosscheck_leg helper behaviour
# ---------------------------------------------------------------------------


def test_crosscheck_leg_indexes_and_passes(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """``crosscheck_leg`` carries the supplied index and agrees on leg 0."""
    leg = aldrin_cycler.legs[0]
    assert leg.n_revs == 0
    result = crosscheck_leg(leg, aldrin_cycler, astropy_ephem, leg_index=7)
    assert result.leg_index == 7
    assert result.passed is True
    assert result.max_diff_mps < V1_TOLERANCE_MPS
    # The in-house departure velocity is reported, not zeroed.
    assert any(abs(c) > 0.0 for c in result.mine_v1_kms)


def test_crosscheck_leg_multirev_raises(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """A multi-rev leg routes to ``NotImplementedError`` (single-rev only)."""
    base = aldrin_cycler.legs[0]
    multirev = Leg(
        from_body=base.from_body,
        to_body=base.to_body,
        t_depart=base.t_depart,
        t_arrive=base.t_arrive,
        v_depart=base.v_depart,
        v_arrive=base.v_arrive,
        n_revs=1,
        branch=base.branch,
    )
    with pytest.raises(NotImplementedError):
        crosscheck_leg(multirev, aldrin_cycler, astropy_ephem, leg_index=0)


def test_crosscheck_leg_missing_endpoint_raises(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """A leg whose epochs match no encounter raises ``ValueError``."""
    base = aldrin_cycler.legs[0]
    orphan = Leg(
        from_body=base.from_body,
        to_body=base.to_body,
        t_depart=base.t_depart + 1.0e9,
        t_arrive=base.t_arrive + 1.0e9,
        v_depart=base.v_depart,
        v_arrive=base.v_arrive,
        n_revs=0,
        branch=base.branch,
    )
    with pytest.raises(ValueError, match="no matching encounter"):
        crosscheck_leg(orphan, aldrin_cycler, astropy_ephem)


# ---------------------------------------------------------------------------
# crosscheck_cycler aggregation
# ---------------------------------------------------------------------------


def test_crosscheck_cycler_skips_multirev_legs(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """Multi-rev legs are skipped; only single-rev legs are represented."""
    legs = list(aldrin_cycler.legs)
    multirev = Leg(
        from_body=legs[0].from_body,
        to_body=legs[0].to_body,
        t_depart=legs[0].t_depart,
        t_arrive=legs[0].t_arrive,
        v_depart=legs[0].v_depart,
        v_arrive=legs[0].v_arrive,
        n_revs=2,
        branch=legs[0].branch,
    )
    mixed = Cycler(
        bodies=aldrin_cycler.bodies,
        period=aldrin_cycler.period,
        encounters=aldrin_cycler.encounters,
        legs=[*legs, multirev],
        sense=aldrin_cycler.sense,
    )
    results = crosscheck_cycler(mixed, astropy_ephem)
    assert len(results) == sum(1 for leg in mixed.legs if leg.n_revs == 0)
    assert all(r.passed for r in results)


def test_crosscheck_cycler_all_multirev_returns_empty(
    aldrin_cycler: Cycler,
    astropy_ephem: Ephemeris,
) -> None:
    """A cycler of only multi-rev legs yields an empty tuple ("N/A")."""
    legs = tuple(
        Leg(
            from_body=leg.from_body,
            to_body=leg.to_body,
            t_depart=leg.t_depart,
            t_arrive=leg.t_arrive,
            v_depart=leg.v_depart,
            v_arrive=leg.v_arrive,
            n_revs=1,
            branch=leg.branch,
        )
        for leg in aldrin_cycler.legs
    )
    all_multirev = Cycler(
        bodies=aldrin_cycler.bodies,
        period=aldrin_cycler.period,
        encounters=aldrin_cycler.encounters,
        legs=list(legs),
        sense=aldrin_cycler.sense,
    )
    assert crosscheck_cycler(all_multirev, astropy_ephem) == ()
