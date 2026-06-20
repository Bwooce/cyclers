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


@dataclass(frozen=True)
class GenericReturn:
    """A single generic-return solution in Russell's referencing-angle space.

    Fields
    ------
    psi_deg:
        Referencing angle psi (degrees) of the v_inf vector -- see
        :func:`psi_of_vinf_vec`.
    tof_body_periods:
        Time of flight in body orbital periods.
    a_au:
        Heliocentric transfer semi-major axis (AU).
    n_revs:
        Number of complete heliocentric revolutions on the transfer.
    branch:
        Lambert/transfer branch label, e.g. ``"slow"`` or ``"fast"``.
    vinf:
        Hyperbolic excess speed magnitude (canonical AU/TU).
    """

    psi_deg: float
    tof_body_periods: float
    a_au: float
    n_revs: int
    branch: str
    vinf: float


def psi_of_vinf_vec(
    vinf_vec: np.typing.ArrayLike,
    r_B: np.typing.ArrayLike,  # noqa: N803 -- Russell's body-vector notation r_B/v_B
    v_B: np.typing.ArrayLike,  # noqa: N803 -- Russell's body-vector notation r_B/v_B
) -> float:
    """Referencing angle psi of a v_inf vector (radians).

    Russell (2004) Sec. 2.7.3, verbatim: "the angular coordinate ... [of the
    v_inf solution] location in the ecliptic plane on the v_inf sphere
    referenced to v_B, with positive being aligned with r_B".

    That is, psi is the in-plane angle of the v_inf vector measured FROM the
    body velocity direction ``v_hat_B``, POSITIVE rotating toward the body
    position direction ``r_hat_B``.

    The in-plane orthonormal basis is built as::

        e1 = v_hat_B
        e2 = (r_hat_B - (r_hat_B . e1) e1) normalised

    so ``e2`` is the component of ``r_hat_B`` orthogonal to ``e1`` (pointing
    toward ``r_B``), giving positive psi toward ``r_B`` as Russell defines.
    Then ``psi = atan2(vinf . e2, vinf . e1)``.
    """
    vinf = np.asarray(vinf_vec, dtype=np.float64)
    r = np.asarray(r_B, dtype=np.float64)
    v = np.asarray(v_B, dtype=np.float64)

    e1 = v / np.linalg.norm(v)
    rhat = r / np.linalg.norm(r)
    e2 = rhat - np.dot(rhat, e1) * e1
    e2 = e2 / np.linalg.norm(e2)

    return float(np.arctan2(np.dot(vinf, e2), np.dot(vinf, e1)))
