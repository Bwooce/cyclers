"""STAGE 3 — robust multi-seed window finder tests.

GOLDEN-TEST DISCIPLINE: leg-duration arithmetic is deterministic geometry
(`# INVARIANT`). Window-ranking / dedup / count assertions follow the
function contract (`# COMPUTED`); they are NOT re-asserted sourced V_inf
values. Sourced V_inf anchors (Aldrin 6.5/9.7) are used as inputs only.
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.phase_match import (
    LaunchWindow,
    PhaseSignature,
    find_candidate_windows,
    leg_duration_seeds,
)

SPD = SECONDS_PER_DAY


# ---------------------------------------------------------------------------
# leg_duration_seeds
# ---------------------------------------------------------------------------


def test_leg_duration_seeds_returns_symmetric_plus_perturbations() -> None:
    """2-leg (146 d, 634 d) + 5 perturb_fracs → 5 signatures, legs > 0,
    total period conserved. `# INVARIANT` (period conservation is geometry)."""
    d1, d2 = 146.0 * SPD, 634.0 * SPD
    period = d1 + d2
    seeds = leg_duration_seeds(
        bodies=("E", "M", "E"),
        primary_leg_durations_s=(d1, d2),
        vinf_target_kms=(6.5, 9.7, 6.5),
        period_s=period,
        perturb_fracs=(0.0, 0.10, -0.10, 0.20, -0.20),
    )
    assert len(seeds) == 5
    for sig in seeds:
        assert all(leg > 0.0 for leg in sig.leg_durations_s)
        # INVARIANT: period conservation under the redistribution arithmetic.
        assert sum(sig.leg_durations_s) == pytest.approx(period, rel=1e-12)
        assert sig.bodies == ("E", "M", "E")
        assert sig.vinf_target_kms == (6.5, 9.7, 6.5)


def test_leg_duration_seeds_clips_to_min_leg() -> None:
    """A -50% perturbation on a 60-d leg clips to min_leg_days=30 d rather
    than going negative. `# COMPUTED` (clip boundary)."""
    d1, d2 = 60.0 * SPD, 540.0 * SPD
    period = d1 + d2
    seeds = leg_duration_seeds(
        bodies=("E", "M", "E"),
        primary_leg_durations_s=(d1, d2),
        vinf_target_kms=(6.5, 9.7, 6.5),
        period_s=period,
        perturb_fracs=(0.0, -0.5),
        min_leg_days=30.0,
    )
    for sig in seeds:
        # COMPUTED: every leg respects the 30-d floor.
        assert all(leg >= 30.0 * SPD - 1.0 for leg in sig.leg_durations_s)


def test_leg_duration_seeds_three_leg_conserves_period() -> None:
    """N=3 family: perturbation spreads across the remaining legs and total
    period is conserved when no clip triggers. `# INVARIANT`."""
    d = (200.0 * SPD, 300.0 * SPD, 400.0 * SPD)
    period = sum(d)
    seeds = leg_duration_seeds(
        bodies=("E", "M", "E", "E"),
        primary_leg_durations_s=d,
        vinf_target_kms=(5.0, 3.0, 5.0, 5.0),
        period_s=period,
        perturb_fracs=(0.0, 0.10, -0.10),
    )
    assert len(seeds) == 3
    for sig in seeds:
        assert all(leg > 0.0 for leg in sig.leg_durations_s)
        assert sum(sig.leg_durations_s) == pytest.approx(period, rel=1e-12)


# ---------------------------------------------------------------------------
# find_candidate_windows
# ---------------------------------------------------------------------------


def test_find_candidate_windows_deduplicates_close_windows() -> None:
    """Two synthetic seeds producing windows within the dedup band collapse to
    one survivor with the lower mismatch_kms. `# COMPUTED` (dedup contract).

    Uses a monkeypatch-free synthetic by driving find_real_windows through a
    tiny ephemeris range where both seeds land near-identical windows; instead
    we exercise the dedup directly with a stub via the public function on a
    stubbed find_real_windows.
    """
    import cyclerfinder.search.phase_match as pm

    base = datetime(2030, 1, 1, tzinfo=UTC)
    sig_a = PhaseSignature(("E", "M"), (146.0 * SPD,), (6.5, 9.7))
    sig_b = PhaseSignature(("E", "M"), (150.0 * SPD,), (6.5, 9.7))

    def fake_find_real_windows(sig, ephem, date_range, n, step_days, mismatch_cap_kms):  # type: ignore[no-untyped-def]
        if sig is sig_a:
            return [LaunchWindow(base, 2.0, (6.5, 9.7))]
        return [LaunchWindow(base.replace(day=11), 0.7, (6.5, 9.7))]  # +10 d, lower mismatch

    orig = pm.find_real_windows
    pm.find_real_windows = fake_find_real_windows  # type: ignore[assignment]
    try:
        windows = find_candidate_windows(
            [sig_a, sig_b],
            Ephemeris(model="circular"),
            (base, base.replace(month=2)),
            dedup_window_days=30.0,
        )
    finally:
        pm.find_real_windows = orig

    # COMPUTED: within 30 d → single survivor, the lower-mismatch one.
    assert len(windows) == 1
    assert windows[0].mismatch_kms == 0.7


def test_find_candidate_windows_ranks_by_mismatch_not_calendar() -> None:
    """Primary regression test for the proximity-bias bug. seed_A window near
    priority with mismatch 3.5; seed_B window 4 yr away with mismatch 0.8 →
    windows[0] is seed_B (0.8). `# COMPUTED` (ranking contract; 3.5/0.8 are
    synthetic inputs, not sourced values)."""
    import cyclerfinder.search.phase_match as pm

    near = datetime(2030, 1, 1, tzinfo=UTC)
    far = datetime(2034, 1, 1, tzinfo=UTC)
    sig_a = PhaseSignature(("E", "M"), (146.0 * SPD,), (6.5, 9.7))
    sig_b = PhaseSignature(("E", "M"), (160.0 * SPD,), (6.5, 9.7))

    def fake_find_real_windows(sig, ephem, date_range, n, step_days, mismatch_cap_kms):  # type: ignore[no-untyped-def]
        if sig is sig_a:
            return [LaunchWindow(near, 3.5, (6.5, 9.7))]
        return [LaunchWindow(far, 0.8, (6.5, 9.7))]

    orig = pm.find_real_windows
    pm.find_real_windows = fake_find_real_windows  # type: ignore[assignment]
    try:
        windows = find_candidate_windows(
            [sig_a, sig_b],
            Ephemeris(model="circular"),
            (near, far),
        )
    finally:
        pm.find_real_windows = orig

    # COMPUTED: lowest mismatch wins regardless of calendar proximity.
    assert windows[0].mismatch_kms == 0.8


def test_find_candidate_windows_aldrin_merges_seeds() -> None:
    """Aldrin sig (E-M, 146 d, vinf 6.5/9.7) + 2 perturbed variants over
    2026-2036 → len(windows) >= 3, all mismatch_kms ascending. `# COMPUTED`
    (count and ordering, not magnitudes; 6.5/9.7 are sourced inputs)."""
    ephem = Ephemeris(model="astropy")
    seeds = leg_duration_seeds(
        bodies=("E", "M", "E"),
        primary_leg_durations_s=(146.0 * SPD, (780.0 - 146.0) * SPD),
        vinf_target_kms=(6.5, 9.7, 6.5),
        period_s=780.0 * SPD,
        perturb_fracs=(0.0, 0.10, -0.10),
    )
    windows = find_candidate_windows(
        seeds,
        ephem,
        (datetime(2026, 1, 1, tzinfo=UTC), datetime(2036, 1, 1, tzinfo=UTC)),
        n=10,
        step_days=10.0,
        mismatch_cap_kms=5.0,
    )
    assert len(windows) >= 3, f"expected >= 3 merged windows, got {len(windows)}"
    mismatches = [w.mismatch_kms for w in windows]
    assert mismatches == sorted(mismatches)
