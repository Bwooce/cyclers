"""#430 reproducibility guard: the global-engine rebuild reuses close_epoch_locked
unchanged, so the closure PHYSICS must stay byte-stable. Pins a fixed candidate's
output (values captured from HEAD before/during the rebuild)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import EpochLockedTrajectory, close_epoch_locked

# Captured from close_epoch_locked at the time the rebuild guard was written.
CAPTURED_CLOSURE = 3.972747771731968
CAPTURED_CONTINUITY = 0.0


def test_close_epoch_locked_pinned_em_candidate() -> None:
    eph = Ephemeris("astropy")
    c = close_epoch_locked(
        EpochLockedTrajectory(
            sequence=("E", "M"),
            leg_tofs_days=(250.0,),
            vinf_kms_at_encounters=(6.5, 9.7),
            launch_epoch_utc="2031-03-01T00:00:00",
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc="2031-03-01T00:00:00",
            validity_window_end_utc="2032-03-01T00:00:00",
            inserts_into="aldrin-classic-em-k1-outbound",
        ),
        eph,
        closure_tol_kms=1.0e6,
        flyby_continuity_tol_kms=1.0e6,
        independent_cross_check=False,
        independent_tol_kms=1.0e6,
    )
    assert c.closure_residual_kms == pytest.approx(CAPTURED_CLOSURE, rel=1e-9)
    assert c.flyby_continuity_max_dv_kms == pytest.approx(CAPTURED_CONTINUITY, rel=1e-9)
