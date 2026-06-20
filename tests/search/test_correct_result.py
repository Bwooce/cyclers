"""M-ED Phase 1: BallisticClosureResult shape (plan Phase 1)."""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import BallisticClosureResult, ballistic_correct


def test_result_fields_present() -> None:
    r = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0, 379.0, 1027.0),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4, 5.6, 5.6),
        converged=True,
        bend_feasible=True,
    )
    assert r.converged is True
    assert r.bend_feasible is True
    assert r.max_residual_kms == 0.04
    assert len(r.tof_days) == 3
    assert len(r.vinf_per_encounter_kms) == 4


def test_constraints_satisfied_is_converged_and_feasible() -> None:
    """constraints_satisfied = converged AND bend_feasible AND vinf-cap met
    (spec §2.2)."""
    r = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0,),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4),
        converged=True,
        bend_feasible=True,
        vinf_cap_ok=True,
    )
    assert r.constraints_satisfied is True
    r2 = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0,),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4),
        converged=True,
        bend_feasible=False,
        vinf_cap_ok=True,
    )
    assert r2.constraints_satisfied is False


def test_lm_raises_on_underdetermined() -> None:
    """Default ``method="lm"`` reproduces the #388 m<n ValueError unchanged.

    A 2-encounter (E-M) vector-mode call has n_res = 2*(2-2)+1 = 1 residual but
    x0 = [t0, tof] = 2 free vars, i.e. m=1 < n=2 (under-determined); ``lm`` raises.
    """
    with pytest.raises(ValueError, match="lm"):
        ballistic_correct(
            sequence=("E", "M"),
            per_leg_revs=(0,),
            per_leg_branch=("single",),
            t0_seed_sec=0.0,
            tof_seed_days=(250.0,),
            period_sec=500.0 * 86400.0,
            ephem=Ephemeris("circular"),
            vinf_cap=99.0,
            residual_mode="vector",
        )


def test_trf_handles_underdetermined() -> None:
    """``method="trf"`` returns a result on the same m<n call (no ValueError)."""
    r = ballistic_correct(
        sequence=("E", "M"),
        per_leg_revs=(0,),
        per_leg_branch=("single",),
        t0_seed_sec=0.0,
        tof_seed_days=(250.0,),
        period_sec=500.0 * 86400.0,
        ephem=Ephemeris("circular"),
        vinf_cap=99.0,
        residual_mode="vector",
        method="trf",
    )
    assert isinstance(r, BallisticClosureResult)
    assert isinstance(r.converged, bool)
    assert r.max_residual_kms >= 0.0
