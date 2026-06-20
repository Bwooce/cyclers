from __future__ import annotations

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import NarcSeed, russell_parent_to_ballistic_seed


def test_bridge_builds_real_seed() -> None:
    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    assert isinstance(seed, NarcSeed)
    assert seed.sequence == ("E", "E", "M", "M")
    nlegs = len(seed.sequence) - 1
    assert len(seed.per_leg_revs) == nlegs
    assert len(seed.per_leg_branch) == nlegs
    assert len(seed.tof_seed_days) == nlegs
    assert abs(seed.tof_seed_days[0] - 533.70) < 1.0
    assert abs(seed.tof_seed_days[1] - 150.0) < 1.0
    assert abs(seed.tof_seed_days[2] - 1026.21) < 1.0
    real_syn_yr = 1.0 / (1.0 / 1.0 - 1.0 / 1.8808)
    assert abs(seed.period_sec - cyc.p * real_syn_yr * 365.25 * 86400.0) < 1.0e7
    assert seed.vinf_anchor_e_kms > 0 and seed.vinf_anchor_m_kms > 0
    assert all(b in ("single", "low", "high") for b in seed.per_leg_branch)
