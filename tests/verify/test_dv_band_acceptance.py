"""Band-aware maintenance-ΔV acceptance — task #420.

Asserts the coupling of the v4.8 Axis-B ``dv_band`` taxonomy
(``docs/notes/2026-06-22-dv-band-definitions.md``) to the §14 V2 closure
acceptance: a strictly-ballistic row is judged against < 1 m/s / 7 cycles, an
essentially-ballistic row against < 10, a powered-DSM row against its
STRICTLY-POSITIVE budget (not zero), and a ``dv_band=None`` row falls back to
the existing generic criterion.

These are fast, pure-Python tests (no DE440 / BVP solve) — they exercise the
threshold map, the scaled acceptance function, and the V2-powered gate's
per-cycle ΔV decision helper directly, using the four banded catalogue rows as
fixtures where helpful.
"""

from __future__ import annotations

import pytest

from cyclerfinder.verify.dv_band_acceptance import (
    RUSSELL_BASIS_CYCLES,
    accept_maintenance_dv,
    assign_dv_band_from_measurement,
    classify_dv_band,
    dv_band_threshold,
)
from cyclerfinder.verify.v2_powered import _MAINTENANCE_DV_SANITY_MAX_KMS, _maintenance_dv_ok

# The four banded catalogue rows (id -> dv_band) used as fixtures.
BANDED_ROWS: dict[str, str] = {
    "aldrin-classic-em-k1-outbound": "powered_dsm",
    "aldrin-classic-em-k1-inbound": "powered_dsm",
    "mcconaghy-2006-em-k2": "essentially_ballistic",
    "liang-2024-cgcec-ephemeris-2033": "strictly_ballistic",
}


# ---------------------------------------------------------------------------
# Threshold map
# ---------------------------------------------------------------------------


def test_strictly_ballistic_window_is_sub_1mps() -> None:
    w = dv_band_threshold("strictly_ballistic")
    assert w is not None
    assert w.lower_mps == 0.0
    assert w.upper_mps == pytest.approx(1.0)


def test_essentially_ballistic_window_is_sub_10mps() -> None:
    w = dv_band_threshold("essentially_ballistic")
    assert w is not None
    assert w.lower_mps == 0.0
    assert w.upper_mps == pytest.approx(10.0)


def test_low_maintenance_window_is_sub_300mps() -> None:
    w = dv_band_threshold("low_maintenance")
    assert w is not None
    assert w.lower_mps == 0.0
    assert w.upper_mps == pytest.approx(300.0)


def test_powered_dsm_window_is_positive_floor_not_zero() -> None:
    """powered_dsm must close to a STRICTLY-POSITIVE budget, not zero. Its lower
    bound is the top of Russell's net (300 m/s / 7 cycles); its upper bound is
    the project-convention V2-powered sanity ceiling (3.5 km/s/cycle in the
    7-cycle basis), NOT a sourced tight tier."""
    w = dv_band_threshold("powered_dsm")
    assert w is not None
    assert w.lower_mps == pytest.approx(300.0)  # top of Russell's net, not 0
    assert w.upper_mps == pytest.approx(3_500.0 * RUSSELL_BASIS_CYCLES)
    assert w.upper_mps > w.lower_mps


def test_null_band_has_no_window() -> None:
    assert dv_band_threshold(None) is None


def test_low_thrust_sep_has_no_impulsive_window() -> None:
    """low_thrust_sep is a regime, not an impulsive m/s ceiling -> no window."""
    assert dv_band_threshold("low_thrust_sep") is None


def test_unknown_band_has_no_window() -> None:
    assert dv_band_threshold("not-a-band") is None


# ---------------------------------------------------------------------------
# accept_maintenance_dv — scaled to the propagated cycle count
# ---------------------------------------------------------------------------


def test_strictly_ballistic_accepts_below_1mps_per_7cycles() -> None:
    # liang-2024 reports ~1e-6 m/s/cycle -> ~7e-6 m/s over 7 cycles: accepted.
    assert accept_maintenance_dv(
        7e-6, dv_band="strictly_ballistic", n_cycles=7, generic_max_mps=120.0
    )
    # 5 m/s over 7 cycles is essentially-ballistic territory: NOT strictly.
    assert not accept_maintenance_dv(
        5.0, dv_band="strictly_ballistic", n_cycles=7, generic_max_mps=120.0
    )


def test_essentially_ballistic_accepts_below_10_rejects_above() -> None:
    # mcconaghy-2006 ~10 m/s tier: 8 m/s / 7 cycles accepted, 50 m/s rejected.
    assert accept_maintenance_dv(
        8.0, dv_band="essentially_ballistic", n_cycles=7, generic_max_mps=120.0
    )
    assert not accept_maintenance_dv(
        50.0, dv_band="essentially_ballistic", n_cycles=7, generic_max_mps=120.0
    )


def test_powered_dsm_rejects_zero_accepts_budget() -> None:
    """A ΔV≈0 'ballistic neighbour' must NOT satisfy powered_dsm; a stated
    powered budget (above Russell's 300 m/s / 7-cycle net) must — and crucially
    is NOT failed for being non-ballistic."""
    # Zero close: NOT a powered cycle.
    assert not accept_maintenance_dv(0.0, dv_band="powered_dsm", n_cycles=7, generic_max_mps=120.0)
    # Below the 300 m/s/7-cycle floor: still NOT powered.
    assert not accept_maintenance_dv(
        100.0, dv_band="powered_dsm", n_cycles=7, generic_max_mps=120.0
    )
    # A sourced-scale powered budget (~1.73-2.04 km/s / 7 cycles,
    # Byrnes-Longuski-Aldrin 1993) is accepted: above the floor, below the
    # sanity ceiling. NOT failed for being non-ballistic.
    assert accept_maintenance_dv(1900.0, dv_band="powered_dsm", n_cycles=7, generic_max_mps=120.0)
    # The V2-powered gate's over-estimating surrogate (~2.9 km/s/cycle ->
    # ~20 km/s / 7 cycles) is still within the project-convention sanity window
    # (it is a degenerate rejector at 3.5 km/s/cycle), so the recorded-V2 Aldrin
    # is NOT regressed by the band coupling.
    assert accept_maintenance_dv(
        7 * 2900.0, dv_band="powered_dsm", n_cycles=7, generic_max_mps=120.0
    )


def test_null_band_uses_generic_path() -> None:
    """dv_band=None falls back to the caller's generic criterion verbatim."""
    assert accept_maintenance_dv(119.0, dv_band=None, n_cycles=7, generic_max_mps=120.0)
    assert not accept_maintenance_dv(121.0, dv_band=None, n_cycles=7, generic_max_mps=120.0)
    # A row is NOT promoted on a band it does not carry: a tiny ΔV that would
    # pass strictly_ballistic still rides the generic ceiling when band is None.
    assert accept_maintenance_dv(0.0, dv_band=None, n_cycles=7, generic_max_mps=120.0)


def test_window_scales_with_cycle_count() -> None:
    """The 7-cycle window pro-ratas to the propagated cycle count."""
    # essentially_ballistic = 10 m/s / 7 cycles -> 10/7 m/s per single cycle.
    per_cycle_ceiling = 10.0 / RUSSELL_BASIS_CYCLES
    assert accept_maintenance_dv(
        per_cycle_ceiling - 1e-6,
        dv_band="essentially_ballistic",
        n_cycles=1,
        generic_max_mps=120.0,
    )
    assert not accept_maintenance_dv(
        per_cycle_ceiling + 1e-6,
        dv_band="essentially_ballistic",
        n_cycles=1,
        generic_max_mps=120.0,
    )


def test_negative_dv_raises() -> None:
    with pytest.raises(ValueError):
        accept_maintenance_dv(-1.0, dv_band=None, n_cycles=7, generic_max_mps=120.0)


def test_nonpositive_cycle_count_raises() -> None:
    with pytest.raises(ValueError):
        accept_maintenance_dv(1.0, dv_band=None, n_cycles=0, generic_max_mps=120.0)


# ---------------------------------------------------------------------------
# V2-powered gate per-cycle decision helper (_maintenance_dv_ok)
# ---------------------------------------------------------------------------


def test_v2powered_generic_path_unchanged_when_band_none() -> None:
    """dv_band=None keeps the exact pre-#420 criterion: 0 < dv < 3.5 km/s."""
    assert _maintenance_dv_ok(2.9, dv_band=None)  # the Aldrin in-family value
    assert not _maintenance_dv_ok(0.0, dv_band=None)  # ballistic neighbour rejected
    assert not _maintenance_dv_ok(_MAINTENANCE_DV_SANITY_MAX_KMS + 0.1, dv_band=None)


def test_v2powered_powered_dsm_rejects_zero_accepts_aldrin_budget() -> None:
    """Under powered_dsm, a ΔV≈0 close fails (not powered) but the Aldrin
    in-family surrogate ~2.9 km/s/cycle passes — within the per-cycle powered
    window [300/7, 24500/7] = [42.9, 3500] m/s, so the recorded-V2 Aldrin is
    NOT regressed."""
    assert not _maintenance_dv_ok(0.0, dv_band="powered_dsm")
    # 2.9 km/s = 2900 m/s/cycle: above the 42.9 m/s floor, below the 3500 m/s
    # sanity ceiling -> accepted.
    assert _maintenance_dv_ok(2.9, dv_band="powered_dsm")
    # A ΔV just below the per-cycle floor (42.9 m/s) is rejected as not powered.
    assert not _maintenance_dv_ok(0.030, dv_band="powered_dsm")
    # Equivalence to the scaled helper at n_cycles=1.
    assert _maintenance_dv_ok(2.9, dv_band="powered_dsm") == accept_maintenance_dv(
        2900.0, dv_band="powered_dsm", n_cycles=1, generic_max_mps=3500.0
    )


def test_v2powered_strictly_ballistic_accepts_near_zero() -> None:
    """A ballistic band correctly accepts a ΔV≈0 close (the generic > 0 rule is
    NOT re-applied under a band)."""
    assert _maintenance_dv_ok(0.0, dv_band="strictly_ballistic")
    # but rejects an Aldrin-scale powered value.
    assert not _maintenance_dv_ok(2.9, dv_band="strictly_ballistic")


# ---------------------------------------------------------------------------
# classify_dv_band — dv_band as a reproduction OUTPUT (task #422)
# ---------------------------------------------------------------------------


def test_classify_boundary_values_default_7cycle_basis() -> None:
    """The four specified boundary measurements bin to the four bands (7-cycle)."""
    assert classify_dv_band(0.5) == "strictly_ballistic"
    assert classify_dv_band(5.0) == "essentially_ballistic"
    assert classify_dv_band(150.0) == "low_maintenance"
    assert classify_dv_band(1000.0) == "powered_dsm"


def test_classify_is_half_open_at_the_top() -> None:
    """Bins are ``< ceiling`` — a value exactly ON a ballistic ceiling falls into
    the next band (matches the literature 'less than 1 m/s' wording and the
    inclusive-upper dv_band_threshold windows)."""
    assert classify_dv_band(0.999) == "strictly_ballistic"
    assert classify_dv_band(1.0) == "essentially_ballistic"
    assert classify_dv_band(9.999) == "essentially_ballistic"
    assert classify_dv_band(10.0) == "low_maintenance"
    assert classify_dv_band(299.999) == "low_maintenance"
    assert classify_dv_band(300.0) == "powered_dsm"


def test_classify_zero_is_strictly_ballistic() -> None:
    assert classify_dv_band(0.0) == "strictly_ballistic"


def test_classify_scales_to_7cycle_basis() -> None:
    """A measurement over n_cycles is scaled to Russell's 7-cycle basis before
    binning. 1 m/s over 1 cycle == 7 m/s over 7 cycles -> essentially_ballistic;
    over 7 cycles the same 1 m/s total is exactly the strictly/essentially
    boundary -> essentially_ballistic."""
    # 1 m/s over a SINGLE cycle = 7 m/s / 7 cycles -> essentially_ballistic.
    assert classify_dv_band(1.0, n_cycles=1) == "essentially_ballistic"
    # 0.1 m/s over a single cycle = 0.7 m/s / 7 cycles -> strictly_ballistic.
    assert classify_dv_band(0.1, n_cycles=1) == "strictly_ballistic"
    # 700 m/s over 14 cycles = 350 m/s / 7 cycles -> powered_dsm.
    assert classify_dv_band(700.0, n_cycles=14) == "powered_dsm"
    # 200 m/s over 14 cycles = 100 m/s / 7 cycles -> low_maintenance.
    assert classify_dv_band(200.0, n_cycles=14) == "low_maintenance"


def test_classify_never_returns_low_thrust_sep() -> None:
    """low_thrust_sep is a regime, not an inferable magnitude — classify_dv_band
    only ever returns the four impulsive bands."""
    bands = {classify_dv_band(v, n_cycles=7) for v in (0.0, 0.5, 5.0, 150.0, 1000.0, 1e6)}
    assert "low_thrust_sep" not in bands
    assert bands <= {
        "strictly_ballistic",
        "essentially_ballistic",
        "low_maintenance",
        "powered_dsm",
    }


def test_classify_negative_dv_raises() -> None:
    with pytest.raises(ValueError):
        classify_dv_band(-1.0)


def test_classify_nonpositive_cycle_count_raises() -> None:
    with pytest.raises(ValueError):
        classify_dv_band(1.0, n_cycles=0)


# ---------------------------------------------------------------------------
# assign_dv_band_from_measurement — provenance + mismatch flagging (task #422)
# ---------------------------------------------------------------------------


def test_assign_null_band_acquires_computed_band() -> None:
    """A null-band row (the ~215 case) acquires the measured band, marked
    computed-v3."""
    c = assign_dv_band_from_measurement(5.0, n_cycles=7, sourced_dv_band=None)
    assert c.dv_band == "essentially_ballistic"
    assert c.dv_band_source == "computed-v3"
    assert c.measured_band == "essentially_ballistic"
    assert c.sourced_band is None
    assert c.mismatch is False
    assert c.measured_total_dv_mps == pytest.approx(5.0)
    assert c.n_cycles == 7


def test_assign_sourced_band_is_never_overwritten() -> None:
    """A sourced band is kept verbatim even when the measurement classifies
    differently; provenance stays 'sourced'."""
    # Measured is WORSE than sourced -> kept + flagged.
    c = assign_dv_band_from_measurement(1000.0, n_cycles=7, sourced_dv_band="strictly_ballistic")
    assert c.dv_band == "strictly_ballistic"  # NOT overwritten
    assert c.dv_band_source == "sourced"
    assert c.measured_band == "powered_dsm"


def test_assign_mismatch_flagged_when_measured_more_powered() -> None:
    """orbit-closure discipline: a measured band STRICTLY worse (more powered)
    than the sourced band is flagged for human review."""
    c = assign_dv_band_from_measurement(50.0, n_cycles=7, sourced_dv_band="strictly_ballistic")
    assert c.measured_band == "low_maintenance"
    assert c.mismatch is True
    assert "MISMATCH" in c.detail


def test_assign_no_mismatch_when_measured_cheaper_or_equal() -> None:
    """A measured band cheaper than (or equal to) the sourced band is consistent
    — bounded by the sourced ceiling — and is NOT flagged."""
    # Cheaper: sourced powered_dsm, measured strictly_ballistic.
    cheaper = assign_dv_band_from_measurement(0.5, n_cycles=7, sourced_dv_band="powered_dsm")
    assert cheaper.measured_band == "strictly_ballistic"
    assert cheaper.mismatch is False
    # Equal: sourced essentially_ballistic, measured essentially_ballistic.
    equal = assign_dv_band_from_measurement(
        5.0, n_cycles=7, sourced_dv_band="essentially_ballistic"
    )
    assert equal.measured_band == "essentially_ballistic"
    assert equal.mismatch is False


def test_assign_sourced_low_thrust_sep_never_mismatches() -> None:
    """A sourced low_thrust_sep (a regime off the impulsive cost scale) is not
    comparable to an impulsive measurement -> never a mismatch, kept verbatim."""
    c = assign_dv_band_from_measurement(2000.0, n_cycles=7, sourced_dv_band="low_thrust_sep")
    assert c.dv_band == "low_thrust_sep"
    assert c.dv_band_source == "sourced"
    assert c.mismatch is False


def test_assign_uses_7cycle_scaling_for_mismatch() -> None:
    """The mismatch decision uses the SAME 7-cycle-scaled band as classify, so a
    single-cycle measurement is judged on the right basis."""
    # 1 m/s over 1 cycle = 7 m/s / 7 cycles -> essentially_ballistic; worse than
    # a sourced strictly_ballistic -> mismatch.
    c = assign_dv_band_from_measurement(1.0, n_cycles=1, sourced_dv_band="strictly_ballistic")
    assert c.measured_band == "essentially_ballistic"
    assert c.mismatch is True
