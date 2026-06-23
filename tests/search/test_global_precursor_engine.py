from __future__ import annotations

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.global_precursor_engine import (
    eccentric_tp_linkable_radius_au,
    eccentric_tp_seeds,
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
