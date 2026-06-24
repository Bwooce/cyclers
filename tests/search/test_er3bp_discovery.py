from __future__ import annotations

from cyclerfinder.search.er3bp_discovery import Er3bpSeed, standard_family_seeds


def test_standard_family_seeds_returns_usable_floor() -> None:
    seeds = standard_family_seeds(target_e=0.0549)
    assert len(seeds) >= 1
    for s in seeds:
        assert isinstance(s, Er3bpSeed)
        assert s.state0.shape == (6,)
        assert s.period_f > 0.0
        assert 0.0 < s.target_e < 1.0
        assert s.system.primary_name and s.system.secondary_name
        assert s.source
