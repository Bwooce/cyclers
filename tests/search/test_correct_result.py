"""M-ED Phase 1: BallisticClosureResult shape (plan Phase 1)."""

from __future__ import annotations

from cyclerfinder.search.correct import BallisticClosureResult


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
