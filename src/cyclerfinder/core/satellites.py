"""Planet-centric satellite (moon) registry — Tier-1 patched-conic moon tours.

A moon cannot live in :data:`cyclerfinder.core.constants.PLANETS`:
``PlanetData.sma_au`` is intrinsically heliocentric (constants.py:139) and its
mean motion derives from ``MU_SUN``. This is the about-the-primary sibling.

Body-code scheme: FULL MOON NAMES (data/README.md:65-79) —
``"Io"``/``"Europa"``/``"Ganymede"``/``"Callisto"``/``"Titan"``/… — to avoid
collision with the heliocentric V/E/M planet codes. ``mean_motion_deg_day`` is
DERIVED at import from ``sma_km`` + the primary's mu (Kepler III), never
hand-copied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SatelliteData:
    name: str
    code: str  # full moon name (scheme: data/README.md:65-79)
    primary: str  # "Jupiter" | "Saturn" | "Earth" | "Mars"
    mu_km3_s2: float  # moon GM
    radius_eq_km: float
    sma_km: float  # SMA ABOUT THE PRIMARY (km, not AU, not Sun-relative)
    mean_motion_deg_day: float  # about the primary (use mean_motion_deg_day_about)
    safe_alt_km: float


def mean_motion_deg_day_about(sma_km: float, *, mu_primary: float) -> float:
    """Mean motion (deg/day) about a primary, Kepler III (cf. constants.py:149-159)."""
    period_s = 2.0 * math.pi * math.sqrt(sma_km**3 / mu_primary)
    return 360.0 / (period_s / 86400.0)
