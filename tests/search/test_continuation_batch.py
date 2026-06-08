"""Fast mechanics tests for the continuation BATCH driver (no DE440).

The slow per-row V3 evidence gate lives in
``tests/verify/test_continuation_v3_batch.py`` (DE440, marked slow). This module
covers only the cheap mechanics — seed derivation from the catalogue, the
result-packaging, and the basin-selection bookkeeping — against a ramped-element
final target (fast, deterministic, no astropy), so the batch plumbing is gated
without paying the DE440 cost.
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.search import continuation as cont
from cyclerfinder.search import continuation_batch as cb


def test_lift_set_is_the_v1_rows_plus_aldrin() -> None:
    """The lift set is exactly the four Russell V1 rows + the Aldrin outbound."""
    assert cb.RUSSELL_V1_RIDS == (
        "russell-ch4-5.30gGf3",
        "russell-ch4-9.94Gg3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-9.353Gg2",
    )
    assert (*cb.RUSSELL_V1_RIDS, "aldrin-classic-em-k1-outbound") == cb.LIFT_SET_RIDS


def test_seeds_derive_from_catalogue() -> None:
    """Every lift-set row yields a finite (a, e) seed, a sourced launch, and a
    two-body V_inf anchor (the CONSTRAINED inputs)."""
    seeds = {s.rid: s for s in cb.seeds_for_ids()}
    assert set(seeds) == set(cb.LIFT_SET_RIDS)
    for s in seeds.values():
        assert 1.0 < s.a_seed_au < 2.0, (s.rid, s.a_seed_au)
        assert 0.0 < s.e_seed < 0.95, (s.rid, s.e_seed)
        assert s.period_sec > 0.0
        assert s.sourced_launch.tzinfo is not None
        assert set(s.sourced_vinf) == {"E", "M"}
    # Aldrin's seed is the sourced (a, e) taken directly (not aphelion+transit).
    aldrin = seeds["aldrin-classic-em-k1-outbound"]
    assert aldrin.a_seed_au == cb._ALDRIN_A_AU
    assert aldrin.e_seed == cb._ALDRIN_E
    assert aldrin.sourced_launch == datetime(2003, 8, 6, tzinfo=UTC)


def test_russell_row_seed_uses_aphelion_transit() -> None:
    """A Russell row's seed comes from the #106 aphelion+transit derivation."""
    from cyclerfinder.search.free_return import seed_ae_from_aphelion_transit

    seeds = {s.rid: s for s in cb.seeds_for_ids(("russell-ch4-5.30gGf3",))}
    s = seeds["russell-ch4-5.30gGf3"]
    # aphelion 2.17, transit 118 d (catalogue) -> the shared derivation.
    a_exp, e_exp = seed_ae_from_aphelion_transit(2.17, 118.0)
    assert s.a_seed_au == pytest.approx(a_exp)
    assert s.e_seed == pytest.approx(e_exp)


def test_non_lift_row_has_no_seed() -> None:
    """A row without a circular-coplanar seed (no aphelion/transit) raises — it is
    not in the lift set (the honest refusal: continuation cannot seed it)."""
    row = {"id": "no-seed-row", "period": {"years": 4.0}, "orbit_elements": {}, "invariants": {}}
    with pytest.raises(ValueError, match="not in the continuation lift set"):
        cb.seed_from_catalogue_row(row)


def test_run_at_seed_against_ramped_target_is_finite() -> None:
    """A single continuation run at an explicit seed against a ramped (1,1,1)
    target (no DE440) produces a finite result and the full audit trail."""
    seed = cb.BatchRowSeed(
        rid="probe",
        a_seed_au=1.30,
        e_seed=0.257,
        period_sec=4.27 * 365.25 * 86400.0,
        sourced_launch=datetime(2025, 6, 15, tzinfo=UTC),
        sourced_vinf={"E": 4.99, "M": 5.10},
    )
    result = cb.run_continuation_at_seed(
        seed,
        cb.epoch_sec(seed.sourced_launch),
        ladder=(1,),
        final_ephemeris=cont.ramped_ephemeris(1.0, 1.0, 1.0),
    )
    import numpy as np

    assert np.isfinite(result.best_final.max_residual_kms)
    assert 243 in result.skipped
    # The audit trail is the p/e/i/ephemeris steps.
    assert [s.phase for s in result.rungs[0].steps] == ["p-ramp", "e-ramp", "i-ramp", "ephemeris"]


def test_basin_selection_packages_evidence_fields() -> None:
    """The basin sweep against a ramped target packages every V3-evidence field
    (residual, winning rung, emerged V_inf vs anchor, phase-match) — the
    bookkeeping the slow gate consumes. Uses a tiny sweep + ramped target (fast)."""
    seed = cb.BatchRowSeed(
        rid="probe",
        a_seed_au=1.30,
        e_seed=0.257,
        period_sec=4.27 * 365.25 * 86400.0,
        sourced_launch=datetime(2025, 6, 15, tzinfo=UTC),
        sourced_vinf={"E": 4.99, "M": 5.10},
    )
    res = cb.run_continuation_for_seed(
        seed,
        ladder=(1,),
        final_ephemeris=cont.ramped_ephemeris(1.0, 1.0, 1.0),
        sweep_days=10,
        step_days=10,
    )
    # All evidence fields present and self-consistent.
    assert res.seed is seed
    assert set(res.emerged_vinf_kms) == {"E", "M"}
    assert set(res.vinf_offset_kms) == {"E", "M"}
    assert res.window_offset_days >= 0.0
    assert isinstance(res.phase_matched, bool)
    assert isinstance(res.vinf_matched, bool)
    assert isinstance(res.selected, bool)
    # reaches_v3_closure is the conjunction of the four halves.
    assert res.reaches_v3_closure == (
        res.selected and res.ballistic_within_tol and res.vinf_matched and res.phase_matched
    )
