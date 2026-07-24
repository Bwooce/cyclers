"""Tests for `#705`'s epoch-robustness scan (`scripts/run_705_epoch_robustness_scan.py`).

`#705` reuses `#704`'s own module
(:mod:`cyclerfinder.search.ccr4bp_real_ephemeris_consistency`) UNMODIFIED as
a library -- the physics/propagation pipeline itself is already covered by
`#704`'s own test suite
(``tests/search/test_ccr4bp_real_ephemeris_consistency.py``). This test file
covers only the NEW orchestration/aggregation logic `#705` adds: running
`#704`'s own dense-synodic-scan methodology at several epochs and aggregating
duty-cycle / local-minimum results across them.

SPICE-gated (skipped if the URA111 kernel is not installed), matching
`#704`'s own test module's skip convention.
"""

from __future__ import annotations

import datetime as _dt
import math

import numpy as np
import pytest

import scripts.run_705_epoch_robustness_scan as run
from cyclerfinder.search.ccr4bp_real_ephemeris_consistency import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
)

_KERNELS_PRESENT = (
    DEFAULT_LSK_PATH.exists() and DEFAULT_PCK_PATH.exists() and DEFAULT_URA_PATH.exists()
)
_SKIP_REASON = (
    f"URA111 SPICE kernel not installed (looked at {DEFAULT_URA_PATH}); "
    "see docs/notes or scripts/install_uranian_spice.sh"
)


def test_duty_cycle_window_counting_matches_704_own_logic() -> None:
    """`#705`'s own :func:`_duty_cycle` must count distinct sub-threshold
    windows the same way `#704`'s own driver script did: contiguous runs of
    ``True`` in the boolean threshold mask, including WRAP-adjacent runs
    treated as separate (matching `#704`'s own non-circular counting)."""
    # Two candidate windows (positions 1-2, and position 5), tuned so the two
    # available thresholds (500/1000/2000/5000/10000/20000/50000 km) resolve
    # them differently: below_1000km only catches the first window (800 km
    # values), below_2000km catches both (1500 km also qualifies).
    pos = np.array([60000.0, 800.0, 800.0, 60000.0, 60000.0, 1500.0, 60000.0])
    duty = run._duty_cycle(pos)
    below_1000 = duty["below_1000km"]
    assert below_1000["n_distinct_windows"] == 1
    assert below_1000["fraction_of_period"] == pytest.approx(2.0 / 7.0)
    below_2000 = duty["below_2000km"]
    assert below_2000["n_distinct_windows"] == 2
    assert below_2000["fraction_of_period"] == pytest.approx(3.0 / 7.0)


def test_comparable_threshold_is_ten_times_704_own_best_point() -> None:
    """The 'comparable tightness' gate is EXACTLY 10x `#704`'s own committed
    best point (``pos_gap_km=84.460...``, ``vel_gap_km_s=0.005905...``) --
    sourced, not re-derived, and explicit about the order-of-magnitude
    factor used."""
    assert pytest.approx(10.0) == run._ORDER_OF_MAGNITUDE_FACTOR
    assert (
        pytest.approx(10.0 * run._REF_704_BEST_POS_GAP_KM) == run._COMPARABLE_POS_GAP_THRESHOLD_KM
    )
    assert (
        pytest.approx(10.0 * run._REF_704_BEST_VEL_GAP_KM_S)
        == run._COMPARABLE_VEL_GAP_THRESHOLD_KM_S
    )
    # These are #704's own committed data/found/704.../result.json:summary
    # values, reproduced here as a literal anti-drift check.
    assert pytest.approx(84.46019822482435) == run._REF_704_BEST_POS_GAP_KM
    assert pytest.approx(0.005904706092575937) == run._REF_704_BEST_VEL_GAP_KM_S


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_epoch_robustness_scan_runs_end_to_end_and_is_finite() -> None:
    """End-to-end smoke test at reduced scale (2 epochs, 12 synodic points
    each -- not `#705`'s own full-resolution production run) against `#701`'s
    OWN actual saved connection, deterministically reconstructed the same
    way `#704`'s own driver script does. Does not assert genuineness or
    recurrence -- only that the NEW multi-epoch orchestration runs end to
    end and produces finite, physically-scaled numbers for every epoch (the
    regression floor for the aggregation logic actually working), matching
    `#704`'s own e2e smoke test's spirit for the single-epoch case.
    """
    t0 = 0.0
    _torus, departure_u, target_s, t_u_tu, refined = run._reconstruct_701_connection(t0)
    assert refined.residual_norm < 1e-10

    import cyclerfinder.core.ccr4bp_umbriel_titania as ut

    system = ut.uranus_umbriel_titania_default()
    l_km = ut.L_KM
    v_unit_km_s = ut.v_unit_km_s()
    synodic_period_days = 7.909  # #704's own committed value, ballpark check only

    base_dts = [
        _dt.datetime(2000, 1, 1, tzinfo=_dt.UTC),
        _dt.datetime(2083, 1, 1, tzinfo=_dt.UTC),
    ]
    n_synodic_saved = run.N_SYNODIC
    try:
        run.N_SYNODIC = 12
        scans = []
        for i, base_dt in enumerate(base_dts):
            scan = run._scan_one_epoch(
                base_dt,
                synodic_period_days,
                departure_u,
                target_s,
                t_u_tu,
                l_km,
                v_unit_km_s,
                system.mu,
                t0,
                i,
            )
            scans.append(scan)
    finally:
        run.N_SYNODIC = n_synodic_saved

    assert len(scans) == 2
    for scan in scans:
        assert math.isfinite(scan["synodic_scan_pos_gap_km_min"])
        assert math.isfinite(scan["synodic_scan_pos_gap_km_max"])
        assert math.isfinite(scan["synodic_scan_pos_gap_km_median"])
        lm = scan["local_minimum"]
        assert lm["propagation_success"]
        assert math.isfinite(lm["pos_gap_km"])
        assert math.isfinite(lm["vel_gap_km_s"])
        # Sanity scale: the local-minimum gap should be well within the
        # system's own physical extent (a total blowup would indicate a
        # wiring bug in the new orchestration, not a genuine physics
        # finding -- matching #704's own e2e smoke test's scale check).
        assert lm["pos_gap_km"] < 5.0e6
        assert len(scan["synodic_phase_scan"]) == 12
        assert isinstance(scan["comparable_to_704_2030_window"], bool)
