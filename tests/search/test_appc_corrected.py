"""Construction mechanics for the generic App-C reconstruction (#170).

Fast (no REBOUND) checks that the reusable App-C primitive reconstructs the reachable
batch parents from their sourced per-leg blocks: correct leg count, 7 Mars-transit
legs, two-body-Sun intercept of the TRUE DE440 Mars at the published per-leg v_inf
with a true-longitude rendezvous. The independent n-body gate + the PARTIAL verdict
live in ``tests/nbody/test_appc_batch_nbody.py``.

EXPECTED = Russell 2004 App-C printed per-leg transit time / v_inf (sourced); the
achieved miss / v_inf / longitude are EVIDENCE.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.appc_corrected import (
    APPC_8049GGF2,
    APPC_8165GFHF2,
    REACHABLE_BLOCKS,
    build_seeded_arcs,
    reconstruct_mars_encounters,
)


@pytest.fixture(scope="module")
def ephem() -> Ephemeris:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Ephemeris("astropy")


def test_reachable_blocks_are_the_two_sourced_parents() -> None:
    assert set(REACHABLE_BLOCKS) == {
        "russell-ch4-8.049gGf2",
        "russell-ch4-8.165Gfh-f2",
    }
    assert APPC_8049GGF2.parent_number == 188
    assert APPC_8165GFHF2.parent_number == 192
    # These are POWERED parents (App-C total Δv far above the 120 m/s V3 budget).
    assert APPC_8049GGF2.total_dv_kms == pytest.approx(0.436091, abs=1e-6)
    assert APPC_8165GFHF2.total_dv_kms == pytest.approx(1.677496, abs=1e-6)


@pytest.mark.parametrize("catalogue_id", sorted(REACHABLE_BLOCKS))
def test_seven_mars_transit_legs(catalogue_id: str, ephem: Ephemeris) -> None:
    block = REACHABLE_BLOCKS[catalogue_id]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        arcs = build_seeded_arcs(block, ephem)
    mars = [a for a in arcs if a.is_mars_transit]
    assert len(mars) == 7
    assert {a.leg_no + 1 for a in mars} == set(block.mars_transit)


@pytest.mark.parametrize("catalogue_id", sorted(REACHABLE_BLOCKS))
def test_mars_encounters_intercept_true_mars_at_published_vinf(
    catalogue_id: str, ephem: Ephemeris
) -> None:
    block = REACHABLE_BLOCKS[catalogue_id]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        encs = reconstruct_mars_encounters(block, ephem)
    assert len(encs) == 7
    for e in encs:
        # Two-body intercept of the TRUE DE440 Mars (well inside an SOI).
        assert e.miss_au < 5e-5, f"{catalogue_id} M{e.arrival_leg_no}: miss {e.miss_au:.2e} AU"
        # v_inf matches the published per-leg value (real-eph, breathes).
        assert e.vinf_kms == pytest.approx(e.pub_vinf_kms, abs=2e-3)
        # True-longitude rendezvous with real Mars.
        dlon = abs(((e.sc_lon_deg - e.mars_lon_deg + 180.0) % 360.0) - 180.0)
        assert dlon < 0.05, f"{catalogue_id} M{e.arrival_leg_no}: Δlon {dlon:.3f} deg"


@pytest.mark.parametrize("catalogue_id", sorted(REACHABLE_BLOCKS))
def test_mars_vinf_breathes_not_a_single_anchor(catalogue_id: str, ephem: Ephemeris) -> None:
    """The real-eph Mars v_inf is epoch-dependent (breathes), not one coplanar value."""
    block = REACHABLE_BLOCKS[catalogue_id]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        encs = reconstruct_mars_encounters(block, ephem)
    vinfs = [e.vinf_kms for e in encs]
    assert max(vinfs) - min(vinfs) > 1.0, "Mars v_inf should breathe across cycles"
    assert np.mean(vinfs) == pytest.approx(block.mars_vinf_avg, abs=0.6)
