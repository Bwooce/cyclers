"""Tier-1 Phase 2: planet-centred circular moon ephemeris (plan Phase 2 Task 2.0)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import SATELLITES


def test_europa_state_is_about_jupiter_at_its_sma() -> None:
    ephem = Ephemeris(model="circular", center="Jupiter")
    r, v = ephem.state("Europa", 0.0)
    # |r| ~ Europa's SMA about Jupiter (km), NOT a heliocentric AU-scaled radius.
    assert np.linalg.norm(r) == pytest.approx(SATELLITES["Europa"].sma_km, rel=1e-9)
    # circular speed |v| = sqrt(mu_jup / a); ~13.7 km/s.
    assert 13.0 < np.linalg.norm(v) < 14.5


def test_heliocentric_default_unchanged() -> None:
    # No center -> the existing heliocentric circular backend, byte-identical.
    helio = Ephemeris(model="circular")
    r, _ = helio.state("E", 0.0)
    assert 1.40e8 < np.linalg.norm(r) < 1.55e8  # ~1 AU in km
