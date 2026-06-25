"""V0 internal-consistency BCR4BP gauntlet tests (#305 Part B).

Gates:
  1. The sourced Andreu POL1 substitute (Sun-commensurate n=1) PASSES V0:
     converged, finite, periapsis above the sourced Earth/Moon floors, and is
     NOT tagged quasi_periodic (drift ~ 0).
  2. An IC placed inside the Earth flyby floor FAILS the periapsis gate.
  3. A non-commensurate orbit (a #303 L1 family member, drift ~ 3.55 rad) is
     TAGGED quasi_periodic and does NOT silently pass as Sun-commensurate.

Validation against a SOURCED family member; the published POL1 IC is the SEED,
not the EXPECTED side of any equality (feedback_golden_tests_sourced_only).
Runs in the DEFAULT suite (not slow).
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v0_bcr4bp import (
    V0_BCR4BP_PHASE_DRIFT_CONVENTION,
    V0VerdictBCR4BP,
    run_v0_bcr4bp,
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


def test_v0_bcr4bp_pol1_passes_and_is_commensurate() -> None:
    """POL1 (n=1) passes V0 and is NOT tagged quasi_periodic."""
    orbit = _close_pol1()
    assert orbit.converged
    verdict = run_v0_bcr4bp("andreu-pol1-bcr4bp", orbit)

    assert isinstance(verdict, V0VerdictBCR4BP)
    assert verdict.converged_corrector
    assert verdict.state_finite
    assert verdict.periapsis_ok, (
        f"POL1 periapsis below floor: earth={verdict.min_periapsis_earth_km:.1f} "
        f"(floor {verdict.earth_floor_km:.1f}), moon={verdict.min_periapsis_moon_km:.1f} "
        f"(floor {verdict.moon_floor_km:.1f})"
    )
    assert verdict.sun_commensurate_n_is_positive_int
    assert not verdict.quasi_periodic, (
        f"POL1 (T fixed at commensurate value) should not be quasi_periodic; "
        f"drift={verdict.sun_phase_drift:.3e}"
    )
    assert verdict.passes_v0_bcr4bp


def test_v0_bcr4bp_sub_floor_earth_periapsis_fails() -> None:
    """An IC inside the Earth flyby floor FAILS the periapsis gate.

    Place a near-circular low orbit a few hundred km above the Earth surface in
    nondim units; its periapsis is below the 6578 km Earth flyby floor, so V0
    must fail on periapsis even though the state is finite.
    """
    sys_bcr = bcr4bp.andreu_default()
    mu = sys_bcr.mu
    # Position ~200 km above Earth (well below the 6578 km floor at periapsis if
    # we make the orbit dip; place it AT ~5000 km from Earth center -> below floor).
    r_earth_nondim = 5000.0 / 384400.0  # ~0.013 LU, < 6578 km floor
    bad_ic = np.array([-mu + r_earth_nondim, 0.0, 0.0, 0.0, 0.5, 0.0], dtype=np.float64)
    bad_orbit = BCR4BPPeriodicOrbit(
        state_initial=bad_ic,
        period_nondim=0.5,
        sun_commensurate_n=1,
        sun_phase_drift=0.0,
        converged=True,
        corrector_residual=0.0,
        independent_closure_residual=0.0,
        n_iter=1,
        system=sys_bcr,
        free_vars=(0, IDX_YDOT),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=False,
        notes="sub-floor Earth periapsis",
    )
    verdict = run_v0_bcr4bp("sub-floor-earth", bad_orbit)
    assert not verdict.periapsis_ok, (
        f"min Earth periapsis {verdict.min_periapsis_earth_km:.1f} km should be "
        f"below floor {verdict.earth_floor_km:.1f} km"
    )
    assert not verdict.passes_v0_bcr4bp


def test_v0_bcr4bp_non_commensurate_member_tagged_quasi_periodic() -> None:
    """A #303-style non-commensurate L1 member is TAGGED quasi_periodic.

    Build a BCR4BPPeriodicOrbit with a large sun_phase_drift (the #303 L1
    family members carry drift ~ 3.55 rad because T is free and not
    Sun-commensurate). V0 must surface this, not hide it.
    """
    sys_bcr = bcr4bp.andreu_default()
    # A #303 L1 member (step 0) IC + period (data/bcr4bp_l1_family_303.jsonl).
    member_ic = np.array(
        [0.811525646809614, 0.0, 0.0, 0.0, 0.2561842600515311, 0.0], dtype=np.float64
    )
    orbit = BCR4BPPeriodicOrbit(
        state_initial=member_ic,
        period_nondim=2.946253088497441,
        sun_commensurate_n=1,
        sun_phase_drift=3.55732377737483,  # echoed from the family jsonl
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
    verdict = run_v0_bcr4bp("l1-303-member", orbit)
    assert verdict.quasi_periodic, (
        f"a member with drift {verdict.sun_phase_drift:.3e} rad >> "
        f"{V0_BCR4BP_PHASE_DRIFT_CONVENTION:.0e} convention must be tagged quasi_periodic"
    )
    # quasi_periodic does NOT veto V0 internal consistency (if periapsis is OK
    # and corrector converged, the orbit is still internally consistent over
    # one nominal period) — but the flag must be set for V2+ to refuse it.
    assert verdict.sun_commensurate_n_is_positive_int


def test_v0_bcr4bp_floors_are_sourced_not_hardcoded() -> None:
    """The Earth/Moon floors trace to the live sourced constants tables."""
    from cyclerfinder.core.constants import PLANETS
    from cyclerfinder.core.satellites import SATELLITES

    orbit = _close_pol1()
    verdict = run_v0_bcr4bp("floor-check", orbit)
    assert verdict.earth_floor_km == PLANETS["E"].radius_eq_km + PLANETS["E"].safe_alt_km
    assert verdict.moon_floor_km == SATELLITES["Moon"].radius_eq_km + SATELLITES["Moon"].safe_alt_km
