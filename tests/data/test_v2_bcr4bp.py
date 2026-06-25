"""V2 long-span BCR4BP multi-lap drift gauntlet tests (#305 Part C, Option A).

Gates:
  1. The sourced Andreu POL1 substitute (Sun-commensurate n=1) propagates >=3
     laps with bounded position drift -> PASSES V2.
  2. A non-commensurate #303 L1 member (drift ~ 3.55 rad) FAILS V2 with the
     explicit reason ``non_commensurate_no_strict_period`` (Option A honesty).

Validation against a SOURCED family member; runs in the DEFAULT suite (not
slow): 3 laps of POL1 at T ~ 6.79 TU complete in a few seconds.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v2_bcr4bp import (
    V2_BCR4BP_DRIFT_FLOOR_KMS,
    V2_BCR4BP_N_CYCLES_MIN,
    V2VerdictBCR4BP,
    run_v2_bcr4bp,
)
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    BCR4BPPeriodicOrbit,
    correct_bcr4bp_periodic,
)

_POL1_X = -0.8369141677649317
_POL1_PY = -0.8391311559808445
_POL1_VY = _POL1_PY - _POL1_X
_POL1_SEED = np.array([_POL1_X, 0.0, 0.0, 0.0, _POL1_VY, 0.0], dtype=np.float64)


def _close_pol1() -> BCR4BPPeriodicOrbit:
    sys_bcr = bcr4bp.andreu_default()
    period_fixed = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)
    return correct_bcr4bp_periodic(
        sys_bcr,
        _POL1_SEED,
        period_fixed,
        sun_commensurate_n=1,
        free_vars=(0, IDX_YDOT),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-10,
        independent_tol=1e-6,
        state_step_cap=0.2,
        require_monotone_decrease=False,
        max_iter=80,
    )


def test_v2_bcr4bp_pol1_commensurate_bounded_passes() -> None:
    """POL1 (n=1, lap on T) stays bounded over 3 laps and PASSES V2."""
    orbit = _close_pol1()
    assert orbit.converged
    verdict = run_v2_bcr4bp("andreu-pol1-bcr4bp", orbit, n_cycles=3)

    assert isinstance(verdict, V2VerdictBCR4BP)
    assert verdict.is_commensurate, (
        f"POL1 (T fixed at commensurate value) should be commensurate; reason={verdict.reason}"
    )
    assert verdict.converged_each_cycle
    assert len(verdict.per_cycle_drift_km) == 3
    assert verdict.max_drift_km <= V2_BCR4BP_DRIFT_FLOOR_KMS, (
        f"POL1 drift {verdict.max_drift_km:.1f} km > floor {verdict.drift_floor_km:.1f} km"
    )
    assert verdict.passes_v2_bcr4bp
    # V3 chains on per_cycle length.
    assert len(verdict.per_cycle_drift_km) >= V2_BCR4BP_N_CYCLES_MIN


def test_v2_bcr4bp_non_commensurate_member_fails_with_reason() -> None:
    """A non-commensurate #303 L1 member FAILS V2 with the explicit reason."""
    sys_bcr = bcr4bp.andreu_default()
    member_ic = np.array(
        [0.811525646809614, 0.0, 0.0, 0.0, 0.2561842600515311, 0.0], dtype=np.float64
    )
    orbit = BCR4BPPeriodicOrbit(
        state_initial=member_ic,
        period_nondim=2.946253088497441,
        sun_commensurate_n=1,
        sun_phase_drift=3.55732377737483,  # NOT commensurate
        converged=True,
        corrector_residual=2.868800169455521e-15,
        independent_closure_residual=9.680441313880588e-09,
        n_iter=1,
        system=sys_bcr,
        free_vars=(0, IDX_YDOT, 6),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        notes="#303 L1 member (non-commensurate)",
    )
    verdict = run_v2_bcr4bp("l1-303-member", orbit, n_cycles=3)

    assert not verdict.passes_v2_bcr4bp
    assert not verdict.is_commensurate
    assert verdict.reason == "non_commensurate_no_strict_period"
    # Honest: no laps were propagated on a meaningless period.
    assert verdict.n_cycles_propagated == 0


def test_v2_bcr4bp_rejects_too_few_cycles() -> None:
    """V2 enforces the spec >=3 lap minimum at the API boundary."""
    orbit = _close_pol1()
    import pytest

    with pytest.raises(ValueError, match="n_cycles >= 3"):
        run_v2_bcr4bp("too-few", orbit, n_cycles=2)
