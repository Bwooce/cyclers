"""Russell (2004) circular-coplanar generic-return model constants.

Russell's simplified Earth-Mars model uses an idealised Mars period of
1.875 yr (chosen so the geometry repeats every 15 yr -- NOT the real
1.881 yr) and a canonical unit system with mu_sun = 1 AU^3/TU^2. These
constants are load-bearing: later golden tests against Russell's tables
only match in this exact model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt

import numpy as np

YEAR_DAYS = 365.25


@dataclass(frozen=True)
class RussellModel:
    """Circular-coplanar two-body Sun model in Russell's canonical units."""

    tu_days: float = 58.1324409
    au_km: float = 149597871.0
    mu_sun: float = 1.0
    mu_earth: float = 3.00348960e-6
    periods_yr: dict[str, float] = field(default_factory=lambda: {"E": 1.0, "M": 1.875})

    def period_yr(self, body: str) -> float:
        """Orbital period of ``body`` in years."""
        return self.periods_yr[body]

    def synodic_yr(self, a: str, b: str) -> float:
        """Synodic period between bodies ``a`` and ``b`` in years."""
        return 1.0 / abs(1.0 / self.period_yr(a) - 1.0 / self.period_yr(b))

    def sma_au(self, body: str) -> float:
        """Semi-major axis (AU) via Kepler III in heliocentric units.

        With P in years and a in AU (Sun-centred), Kepler III is P^2 = a^3, so
        ``a = P_yr^(2/3)`` directly -- Earth (P=1 yr) is a=1.0 AU EXACTLY by
        definition, Mars (P=1.875 yr) is 1.875^(2/3). Periods enter only as a
        ratio, so this is self-consistent independent of the Julian/sidereal year
        and the canonical TU (which is retained for absolute-time conversions).
        """
        return float(self.period_yr(body) ** (2.0 / 3.0))

    def body_circular_speed(self, body: str) -> float:
        """Circular orbital speed (AU/TU) = sqrt(mu_sun / a)."""
        return sqrt(self.mu_sun / self.sma_au(body))

    def body_state(self, body: str, theta: float) -> tuple[np.ndarray, np.ndarray]:
        """Coplanar circular state at in-plane angle ``theta`` (radians).

        Returns prograde position and velocity as float64 ``(3,)`` arrays.
        """
        a = self.sma_au(body)
        v_circ = self.body_circular_speed(body)
        r = np.array([a * np.cos(theta), a * np.sin(theta), 0.0], dtype=np.float64)
        v = np.array([-v_circ * np.sin(theta), v_circ * np.cos(theta), 0.0], dtype=np.float64)
        return r, v
