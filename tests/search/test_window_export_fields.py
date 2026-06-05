"""Tests for the launch-window export-derived columns (C_3, first-leg ToF).

The windows.json exporter (cyclers.space ``scripts/compute_windows.py``)
derives two cheap columns from data this repo already produces:

* per-window ``c3_km2_s2 = |V_inf,depart|^2`` — the square of the first
  encounter's actual V_inf on a :class:`LaunchWindow`.
* per-cycler ``tof_first_leg_days`` — the first-leg ToF of the cycler's
  :class:`PhaseSignature`, in days.

These tests pin the relationships the exporter relies on, so a change to
:class:`LaunchWindow` / :class:`PhaseSignature` that would silently break the
exported columns fails here instead. They are fast (no ephemeris).
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from cyclerfinder.search.phase_match import (
    LaunchWindow,
    phase_signature_from_catalogue_entry,
)

SECONDS_PER_DAY = 86400.0

_ENTRY = {
    "id": "aldrin-classic-em-k1-outbound",
    "bodies": ["E", "M"],
    "vinf_kms_at_encounters": [
        {"body": "E", "vinf_kms": 6.5},
        {"body": "M", "vinf_kms": 9.7},
    ],
    "legs": [{"from": "E", "to": "M", "tof_days": 146}],
}


def _c3_from_window(w: LaunchWindow) -> float:
    """Reproduce the exporter's C_3 derivation for one window."""
    return w.vinf_actual_kms[0] ** 2


def test_c3_equals_vinf_depart_squared() -> None:
    # Synthetic window: departure V_inf 6.4 km/s -> C_3 = 40.96 km^2/s^2.
    w = LaunchWindow(
        departure_date=datetime(2030, 1, 1, tzinfo=UTC),
        mismatch_kms=0.1,
        vinf_actual_kms=(6.4, 9.7),
    )
    assert _c3_from_window(w) == pytest.approx(6.4**2)
    assert _c3_from_window(w) == pytest.approx(40.96)


def test_c3_zero_for_zero_vinf() -> None:
    w = LaunchWindow(
        departure_date=datetime(2030, 1, 1, tzinfo=UTC),
        mismatch_kms=0.0,
        vinf_actual_kms=(0.0, 5.0),
    )
    assert _c3_from_window(w) == pytest.approx(0.0)


def test_tof_first_leg_days_from_signature() -> None:
    sig = phase_signature_from_catalogue_entry(_ENTRY)
    tof_first_days = sig.leg_durations_s[0] / SECONDS_PER_DAY
    assert tof_first_days == pytest.approx(146.0)


def test_export_field_shape_alignment() -> None:
    # The exporter emits one C_3 per window, aligned with vinf_actual_kms.
    windows = [
        LaunchWindow(datetime(2030, 1, 1, tzinfo=UTC), 0.1, (6.4, 9.7)),
        LaunchWindow(datetime(2032, 1, 1, tzinfo=UTC), 0.2, (5.1, 8.3)),
    ]
    c3 = [round(_c3_from_window(w), 3) for w in windows]
    vinf = [[round(v, 3) for v in w.vinf_actual_kms] for w in windows]
    assert len(c3) == len(vinf) == len(windows)
    for c, v in zip(c3, vinf, strict=True):
        assert c == pytest.approx(v[0] ** 2, abs=1e-3)
