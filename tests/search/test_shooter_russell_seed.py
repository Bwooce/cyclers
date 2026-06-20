from __future__ import annotations

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import russell_parent_to_ballistic_seed
from cyclerfinder.search.shooter_russell_seed import russell_shooting_seed


def test_adapter_builds_shooting_seed() -> None:
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    ephem = Ephemeris("astropy")
    sseed = russell_shooting_seed(seed, t0_sec=0.0, ephem=ephem)
    n = len(seed.sequence)
    assert len(sseed.node_states) == n
    assert all(tuple(s.shape) == (6,) for s in sseed.node_states)
    assert len(sseed.tofs) == n - 1
    assert sseed.sequence == seed.sequence
