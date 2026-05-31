"""Tests for phase_match: catalogue-entry signatures + real-window finder."""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.phase_match import (
    LaunchWindow,
    PhaseSignature,
    find_real_windows,
    phase_signature_from_catalogue_entry,
)

# A minimal Aldrin-shaped catalogue entry; matches the real catalogue row
# shape (subset of fields) for entry `aldrin-classic-em-k1-outbound`.
ALDRIN_ENTRY = {
    "id": "aldrin-classic-em-k1-outbound",
    "bodies": ["E", "M"],
    "vinf_kms_at_encounters": [
        {"body": "E", "vinf_kms": 6.5},
        {"body": "M", "vinf_kms": 9.7},
    ],
    "legs": [{"from": "E", "to": "M", "tof_days": 146}],
}


def test_phase_signature_from_aldrin_entry() -> None:
    sig = phase_signature_from_catalogue_entry(ALDRIN_ENTRY)
    assert sig.bodies == ("E", "M")
    assert len(sig.leg_durations_s) == 1
    # 146 d = 12_614_400 s.
    assert sig.leg_durations_s[0] == pytest.approx(146 * 86400.0)
    assert sig.vinf_target_kms == (6.5, 9.7)
    assert sig.primary == "Sun"


def test_phase_signature_rejects_null_tof() -> None:
    bad = dict(ALDRIN_ENTRY)
    bad["legs"] = [{"from": "E", "to": "M", "tof_days": None}]
    with pytest.raises(ValueError, match="null tof_days"):
        phase_signature_from_catalogue_entry(bad)


def test_phase_signature_rejects_null_vinf() -> None:
    bad = dict(ALDRIN_ENTRY)
    bad["vinf_kms_at_encounters"] = [
        {"body": "E", "vinf_kms": None},
        {"body": "M", "vinf_kms": 9.7},
    ]
    with pytest.raises(ValueError, match="null vinf_kms"):
        phase_signature_from_catalogue_entry(bad)


def test_phase_signature_rejects_empty_bodies() -> None:
    bad = dict(ALDRIN_ENTRY)
    bad["bodies"] = []
    with pytest.raises(ValueError, match="empty bodies"):
        phase_signature_from_catalogue_entry(bad)


def test_phase_signature_constructor_validates_lengths() -> None:
    with pytest.raises(ValueError, match="leg_durations_s length"):
        PhaseSignature(
            bodies=("E", "M"),
            leg_durations_s=(),  # should have len 1
            vinf_target_kms=(6.5, 9.7),
        )
    with pytest.raises(ValueError, match="vinf_target_kms length"):
        PhaseSignature(
            bodies=("E", "M"),
            leg_durations_s=(146 * 86400.0,),
            vinf_target_kms=(6.5,),  # should have len 2
        )


def test_find_real_windows_aldrin_2026_2036() -> None:
    """The Aldrin signature should match real Earth-Mars opposition geometry
    multiple times in 2026-2036 (synodic ≈ 2.135 yr → ~4-5 windows fit)."""
    sig = phase_signature_from_catalogue_entry(ALDRIN_ENTRY)
    ephem = Ephemeris(model="astropy")
    windows = find_real_windows(
        sig,
        ephem,
        (datetime(2026, 1, 1, tzinfo=UTC), datetime(2036, 1, 1, tzinfo=UTC)),
        n=5,
        step_days=10.0,
        mismatch_cap_kms=3.0,
    )
    assert len(windows) >= 3, f"expected at least 3 windows in 10-yr range, got {len(windows)}"
    assert len(windows) <= 5
    for w in windows:
        assert isinstance(w, LaunchWindow)
        # All windows must fall inside the requested date range.
        assert (
            datetime(2026, 1, 1, tzinfo=UTC) <= w.departure_date <= datetime(2036, 1, 1, tzinfo=UTC)
        )
        # V∞ at Earth in the actual Lambert should be in a plausible range
        # for a real Earth-Mars 146-d transfer (the Aldrin target is 6.5 km/s;
        # real-ephemeris geometry can shift this by 1-2 km/s).
        assert 4.0 < w.vinf_actual_kms[0] < 10.0
    # Returned windows must be sorted by ascending mismatch.
    mismatches = [w.mismatch_kms for w in windows]
    assert mismatches == sorted(mismatches)


def test_find_real_windows_rejects_non_heliocentric_primary() -> None:
    sig = PhaseSignature(
        bodies=("E", "Moon"),
        leg_durations_s=(86400.0 * 3,),
        vinf_target_kms=(1.0, 1.0),
        primary="Earth",
    )
    ephem = Ephemeris(model="astropy")
    with pytest.raises(NotImplementedError, match="primary='Sun'"):
        find_real_windows(
            sig, ephem, (datetime(2026, 1, 1, tzinfo=UTC), datetime(2026, 2, 1, tzinfo=UTC))
        )
