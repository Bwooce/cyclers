"""#307 Task 1: multi-rev Lambert branch selection in close_epoch_locked.

The closure driver historically used single-rev Lambert only (max_revs=0). Phase 5
adds a ``max_revs`` knob: each non-DSM leg enumerates its multi-rev Lambert branches
and the driver returns the branch combination with the lowest closure residual.

Invariant under test: the single-rev solution is always in the candidate set, so
``max_revs > 0`` can NEVER yield a worse closure residual than ``max_revs = 0``
(monotone-improvement). And ``max_revs = 0`` must be byte-identical to the historical
path (no regression).
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import (
    EpochLockedTrajectory,
    close_epoch_locked,
)


@pytest.fixture(scope="module")
def ephem() -> Ephemeris:
    return Ephemeris("astropy")


def _eme_trajectory() -> EpochLockedTrajectory:
    """A plain Earth-Mars-Earth epoch-locked trajectory (closure need not converge;
    the multi-rev invariant is structural and holds for any valid trajectory)."""
    return EpochLockedTrajectory(
        sequence=("E", "M", "E"),
        leg_tofs_days=(227.0, 274.0),
        vinf_kms_at_encounters=(3.0, 3.0, 3.0),
        launch_epoch_utc="2027-01-01T00:00:00",
        orbit_class="mga_tour",
        n_returns=1,
        validity_window_start_utc="2027-01-01T00:00:00",
        validity_window_end_utc="2027-12-31T00:00:00",
        ephemeris="DE440",
        notes="#307 Task 1 multi-rev invariant test",
    )


def test_max_revs_kwarg_accepted(ephem: Ephemeris) -> None:
    traj = _eme_trajectory()
    c = close_epoch_locked(traj, ephem, independent_cross_check=False, max_revs=1)
    assert c.closure_residual_kms >= 0.0


def test_multirev_never_worse_than_single_rev(ephem: Ephemeris) -> None:
    traj = _eme_trajectory()
    base = close_epoch_locked(traj, ephem, independent_cross_check=False, max_revs=0)
    multi = close_epoch_locked(traj, ephem, independent_cross_check=False, max_revs=2)
    # The single-rev combination is always enumerated at max_revs>0, so the
    # selected closure residual is never worse than the single-rev-only residual.
    assert multi.closure_residual_kms <= base.closure_residual_kms + 1e-9


def test_max_revs_zero_is_unchanged(ephem: Ephemeris) -> None:
    """max_revs=0 (default) must equal an explicit max_revs=0 — the historical path."""
    traj = _eme_trajectory()
    default = close_epoch_locked(traj, ephem, independent_cross_check=False)
    explicit = close_epoch_locked(traj, ephem, independent_cross_check=False, max_revs=0)
    assert default.closure_residual_kms == explicit.closure_residual_kms
    assert default.per_encounter_vinf_kms == explicit.per_encounter_vinf_kms
