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


def test_parent_phase_angle_range() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.narc_continuation import parent_phase_angle

    ph = parent_phase_angle(Ephemeris("astropy"), 0.0)
    assert -3.14159 - 1e-6 <= ph <= 3.14159 + 1e-6


def test_candidate_epochs_match_phase() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.narc_continuation import candidate_epochs, parent_phase_angle

    ephem = Ephemeris("astropy")
    target_phase = 1.0  # rad
    epochs = candidate_epochs(ephem, target_phase, launch_window_synodics=range(1, 6), grid=50)
    assert epochs  # at least one epoch found
    for t0 in epochs:
        # wrapped phase error within ~9 deg of target
        delta = parent_phase_angle(ephem, t0) - target_phase
        err = abs(((delta + 3.14159265) % 6.2831853) - 3.14159265)
        assert err < 0.16


def test_narc_driver_circular_rung_smoke() -> None:
    import cyclerfinder.search.continuation as cont
    from cyclerfinder.data.catalog import load_catalog
    from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
    from cyclerfinder.search.generic_return import RussellModel
    from cyclerfinder.search.narc_continuation import (
        NarcContinuationResult,
        narc_continuation_correct,
        russell_parent_to_ballistic_seed,
    )

    row = load_catalog().by_id["russell-ch4-4.991gG2"].raw
    m = RussellModel()
    phsi = descriptor_to_phsi(row)
    assert phsi is not None
    cyc = assemble_cycler(m, phsi)
    assert cyc is not None
    seed = russell_parent_to_ballistic_seed(m, cyc, row)
    circ = cont.ramped_ephemeris(0.0, 0.0)  # circular-coplanar backend
    res = narc_continuation_correct(
        seed, ladder=(1,), final_ephemeris=circ, epochs=[0.0], vinf_cap=99.0
    )
    assert isinstance(res, NarcContinuationResult)
    assert isinstance(res.converged, bool)
    assert res.max_residual_kms >= 0.0
    assert isinstance(res.emerged_vinf_kms, tuple)
    assert isinstance(res.bend_feasible, bool)
