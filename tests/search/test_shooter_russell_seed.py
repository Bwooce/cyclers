from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.nbody import shooter
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


@pytest.mark.slow
def test_shoot_from_russell_seed_returns_result() -> None:
    """#388 smoke: the constructed Russell parent feeds the full n-body shooter.

    Asserts only that ``shooter.shoot`` returns a non-None :class:`ShootResult`
    (converged or not) — the batch (``scripts/shooter_russell_batch.py``) is the
    real validation; this guards the seed bridge -> shoot() plumbing.
    """
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    ephem = Ephemeris("astropy")
    sseed = russell_shooting_seed(seed, t0_sec=0.0, ephem=ephem)
    bodies = tuple(dict.fromkeys(seed.sequence))
    res = shooter.shoot(
        sseed,
        ephem=ephem,
        bodies=bodies,
        accuracy=1e-9,
        max_nfev=2,
        max_wall_sec=5.0,
        n_jobs=16,
    )
    assert res is not None
    assert isinstance(res.converged, bool)
    assert len(res.vinf_per_encounter_kms) == len(seed.sequence)
