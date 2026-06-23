from __future__ import annotations

import pytest

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.global_precursor_engine import (
    decision_cost,
    eccentric_tp_linkable_radius_au,
    eccentric_tp_seeds,
    rank_band,
)


def test_eccentric_radius_reduces_to_mean_a_for_circular_backend() -> None:
    """The 'circular' ephemeris backend places bodies on circles at their mean a,
    so the eccentric-aware encounter radius must equal the body's sma_au."""
    eph = Ephemeris("circular")
    # t_sec arbitrary; circular backend radius is epoch-invariant in magnitude.
    r_au = eccentric_tp_linkable_radius_au("E", t_sec=0.0, ephemeris=eph)
    assert abs(r_au - PLANETS["E"].sma_au) < 1e-6


def test_eccentric_tp_seeds_returns_candidates_terminating_at_target() -> None:
    """The seeder returns MGAChainCandidates whose final body is the cycler's
    first encounter body and whose terminal V∞ bin is within tol of the seed."""
    eph = Ephemeris("astropy")
    seeds = eccentric_tp_seeds(
        first_body="E",
        seed_vinf_kms=6.5,
        launch_window=("2030-01-01T00:00:00", "2032-12-31T00:00:00"),
        ephemeris=eph,
        intermediate_bodies=("V", "E"),
        max_legs=3,
        vinf_grid_kms=(6.0, 7.0),
        tof_box_days_per_leg=(80.0, 500.0),
        epoch_step_days=120.0,
        vinf_terminal_tol_kms=0.8,
    )
    assert len(seeds) > 0
    for s in seeds:
        assert s.sequence[-1] == "E"
        assert abs(s.vinf_tuple_kms[-1] - 6.5) <= 0.8 + 1e-9


def test_zero_dsm_vector_matches_plain_ballistic_closure() -> None:
    """A decision vector with all DSM magnitudes 0 closes identically to a
    plain EpochLockedTrajectory with no dsm_specs (the DSM layer is a no-op)."""
    from cyclerfinder.genome.epoch_aware_genome import (
        EpochLockedTrajectory,
        close_epoch_locked,
    )
    from cyclerfinder.search.global_precursor_engine import evaluate_decision_vector

    eph = Ephemeris("astropy")
    sequence = ("E", "M")
    leg_tofs = (250.0,)
    vinf_expected = (6.5, 9.7)
    launch = "2031-03-01T00:00:00"

    plain = close_epoch_locked(
        EpochLockedTrajectory(
            sequence=sequence,
            leg_tofs_days=leg_tofs,
            vinf_kms_at_encounters=vinf_expected,
            launch_epoch_utc=launch,
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc=launch,
            validity_window_end_utc="2032-03-01T00:00:00",
            inserts_into="aldrin-classic-em-k1-outbound",
        ),
        eph,
        closure_tol_kms=1.0e6,
        flyby_continuity_tol_kms=1.0e6,
        independent_cross_check=False,
        independent_tol_kms=1.0e6,
    )
    x = [0.0, 250.0, 0.5, 0.0, 0.0, 0.0]
    result = evaluate_decision_vector(
        x,
        sequence=sequence,
        seed_launch_epoch_utc=launch,
        vinf_expected_kms=vinf_expected,
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
    )
    assert abs(result.closure.closure_residual_kms - plain.closure_residual_kms) < 1e-6
    assert result.total_dsm_dv_kms == 0.0


def test_negative_tof_vector_is_infeasible_with_cost_floor() -> None:
    from cyclerfinder.search.global_precursor_engine import evaluate_decision_vector

    eph = Ephemeris("astropy")
    x = [0.0, -10.0, 0.5, 0.0, 0.0, 0.0]
    result = evaluate_decision_vector(
        x,
        sequence=("E", "M"),
        seed_launch_epoch_utc="2031-03-01T00:00:00",
        vinf_expected_kms=(6.5, 9.7),
        ephemeris=eph,
        inserts_into="aldrin-classic-em-k1-outbound",
    )
    assert result.feasible is False
    assert result.closure.closure_residual_kms >= 1.0e6
    assert result.closure.flyby_continuity_max_dv_kms >= 1.0e6
    assert result.closure.converged is False
    assert result.total_dsm_dv_kms >= 1.0e6
    assert result.per_leg_dsm_kms == ()


def test_ballistic_scores_below_equal_powered() -> None:
    """Two evals with identical closure+continuity but different DSM cost:
    the ballistic (0 DSM) one must have the strictly lower cost."""

    class _FakeClosure:
        closure_residual_kms = 2.0
        flyby_continuity_max_dv_kms = 1.0

    class _FakeEval:
        def __init__(self, dsm: float) -> None:
            self.closure = _FakeClosure()
            self.total_dsm_dv_kms = dsm

    ballistic = decision_cost(_FakeEval(0.0))  # type: ignore[arg-type]
    powered = decision_cost(_FakeEval(0.4))  # type: ignore[arg-type]
    assert ballistic < powered


@pytest.mark.parametrize(
    "dv_kms,expected",
    [
        (0.0005, "strictly_ballistic"),  # 0.5 m/s < 1
        (0.005, "essentially_ballistic"),  # 5 m/s < 10
        (0.250, "low_maintenance"),  # 250 m/s < 300
        (0.500, "powered_dsm"),  # 500 m/s >= 300
    ],
)
def test_rank_band_boundaries(dv_kms: float, expected: str) -> None:
    assert rank_band(dv_kms) == expected
