"""Russell (2004) circular-coplanar generic-return model constants.

Russell's simplified Earth-Mars model uses an idealised Mars period of
1.875 yr (chosen so the geometry repeats every 15 yr -- NOT the real
1.881 yr) and a canonical unit system with mu_sun = 1 AU^3/TU^2. These
constants are load-bearing: later golden tests against Russell's tables
only match in this exact model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import degrees, pi, sqrt

import numpy as np

from cyclerfinder.core.lambert import LambertError, lambert

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


# Provisional Lambert-branch -> Russell-branch label mapping (Russell §2.7.4).
# Maps the multi-rev energy branch ("low"->"slow", "high"->"fast") and the
# single-rev transfer ("single"->"slow"). This is a provisional assignment; the
# later golden test A5 (against Russell's tables) will confirm or correct it.
_BRANCH_MAP: dict[str, str] = {"single": "slow", "low": "slow", "high": "fast"}


def generate_generic_returns(
    model: RussellModel,
    body: str,
    *,
    max_tof_body_periods: float = 6.0,
    dtheta_deg: float = 0.5,
    refine_dtheta_deg: float = 1.0 / 24.0,
    max_revs_cap: int = 15,
) -> list[GenericReturn]:
    """Multi-rev Lambert grid of same-body generic returns (Russell §2.7.3-2.7.5).

    A generic return is a free return to the same body: depart ``body`` at
    in-plane angle 0 and return to the body's later position after some time of
    flight. The grid steps the time of flight (measured in body orbital
    periods) and, at each step, solves the multi-revolution Lambert problem
    Sun-to-Sun between the body's departure and arrival positions, recording one
    :class:`GenericReturn` per Lambert solution.

    Parameters
    ----------
    model:
        The :class:`RussellModel` providing the circular-coplanar body states
        and canonical ``mu_sun``.
    body:
        Body key (e.g. ``"E"``) understood by :meth:`RussellModel.body_state`.
    max_tof_body_periods:
        Upper bound on time of flight, in body orbital periods.
    dtheta_deg:
        Grid step in body angle (degrees); the TOF step in body periods is
        ``dtheta_deg / 360``.
    refine_dtheta_deg:
        Reserved for a later refinement pass; accepted but unused here.
    max_revs_cap:
        Maximum heliocentric revolution count passed to :func:`lambert`.

    Returns
    -------
    list[GenericReturn]
        One entry per Lambert solution across the grid (single-rev plus any
        feasible multi-rev branches).
    """
    p_tu = 2.0 * pi * model.sma_au(body) ** 1.5
    dtof = dtheta_deg / 360.0
    step = dtheta_deg / 360.0

    r1, v_b0 = model.body_state(body, 0.0)
    r = float(np.linalg.norm(r1))

    returns: list[GenericReturn] = []

    n_steps = int((max_tof_body_periods - dtof) / step) + 1
    for i in range(n_steps):
        tof_bp = dtof + i * step
        if tof_bp > max_tof_body_periods:
            break

        # Body position is periodic in its angle, so the wrapped angle is fine.
        r2, _ = model.body_state(body, 2.0 * pi * tof_bp)
        tof_canonical = tof_bp * p_tu

        try:
            sols = lambert(
                r1, r2, tof_canonical, mu=model.mu_sun, prograde=True, max_revs=max_revs_cap
            )
        except LambertError:
            continue

        for sol in sols:
            vinf_vec = sol.v1 - v_b0
            vinf = float(np.linalg.norm(vinf_vec))
            if vinf <= 0.0:
                # Degenerate: the transfer coincides with the body's own
                # circular orbit (zero hyperbolic excess) -- not a return.
                continue
            psi_deg = degrees(psi_of_vinf_vec(vinf_vec, r1, v_b0))

            speed = float(np.linalg.norm(sol.v1))
            denom = 2.0 / r - speed * speed / model.mu_sun
            if abs(denom) < 1.0e-12:
                # Parabolic transfer -- semi-major axis undefined; skip.
                continue
            a_au = 1.0 / denom

            branch = _BRANCH_MAP.get(sol.branch, sol.branch)

            returns.append(
                GenericReturn(
                    psi_deg=psi_deg,
                    tof_body_periods=tof_bp,
                    a_au=a_au,
                    n_revs=int(sol.n_revs),
                    branch=branch,
                    vinf=vinf,
                )
            )

    return returns
