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
