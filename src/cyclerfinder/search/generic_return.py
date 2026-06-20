"""Russell (2004) circular-coplanar generic-return model constants.

Russell's simplified Earth-Mars model uses an idealised Mars period of
1.875 yr (chosen so the geometry repeats every 15 yr -- NOT the real
1.881 yr) and a canonical unit system with mu_sun = 1 AU^3/TU^2. These
constants are load-bearing: later golden tests against Russell's tables
only match in this exact model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise
from math import degrees, pi, sqrt

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.kepler import propagate
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
        Russell's tabulated transfer "a" value (his Tables 2.2/2.3 "a (AU)"
        column). Reproduction of those tables (golden test A5) shows this
        column is the transfer *orbital period* in canonical Earth-period
        units, i.e. ``a_sma ** 1.5`` (Kepler III with ``mu_sun = 1``), NOT the
        semi-major axis itself. We therefore report that same quantity here so
        the field mirrors Russell's column exactly; the heliocentric semi-major
        axis is ``a_au ** (2/3)``.
    n_revs:
        Number of complete heliocentric revolutions on the transfer.
    branch:
        Russell slow/fast transfer label (``"slow"`` for his signed ``N > 0``,
        ``"fast"`` for ``N < 0``). Determined by the phase of the final partial
        revolution: the spacecraft past apoapsis on arrival (mean-anomaly
        fraction ``> 1/2``) is the slow branch, before apoapsis is the fast
        branch -- see :func:`_russell_branch`.
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


def _russell_branch(tof_canonical: float, a_sma: float, mu: float) -> str:
    """Russell slow/fast label of a transfer from its final-rev phase.

    Russell (2004) signs his revolution count ``N`` as ``+`` for slow and ``-``
    for fast returns. Reproducing his Tables 2.2 and 2.3 (golden test A5) shows
    the Lambert energy branch ("low"/"high") does NOT track this sign -- the
    energy ordering swaps as the min-energy boundary is crossed. The invariant
    that does track it is the phase of the transfer's final (partial)
    revolution: with ``M`` the mean anomaly swept over the flight time,
    ``frac = (M / 2pi) mod 1`` is the fraction of the last loop completed.
    ``frac > 1/2`` (spacecraft past apoapsis, descending toward periapsis on
    arrival) is Russell's slow branch; ``frac < 1/2`` (still outbound, before
    apoapsis) is the fast branch. This matches all seven anchor rows across both
    tables.
    """
    period = 2.0 * pi * a_sma**1.5
    frac = (tof_canonical / period) % 1.0
    return "slow" if frac > 0.5 else "fast"


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
            a_sma = 1.0 / denom
            if a_sma <= 0.0:
                # Hyperbolic transfer -- no orbital period; not a Russell return.
                continue
            # Russell's "a (AU)" column is the transfer period in canonical
            # Earth-period units (a_sma ** 1.5); see GenericReturn.a_au.
            a_au = a_sma**1.5

            branch = _russell_branch(tof_canonical, a_sma, model.mu_sun)

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


def bin_sub_families(
    returns: list[GenericReturn],
) -> dict[tuple[int, str], list[GenericReturn]]:
    """Group generic returns into sub-families keyed by ``(n_revs, branch)``.

    Each value list is sorted ascending by hyperbolic-excess speed ``.vinf``,
    so that consecutive entries bracket intermediate target speeds for the
    interpolation step in :func:`returns_at_vinf` (Russell §A.2).
    """
    bins: dict[tuple[int, str], list[GenericReturn]] = {}
    for g in returns:
        bins.setdefault((g.n_revs, g.branch), []).append(g)
    for lst in bins.values():
        lst.sort(key=lambda g: g.vinf)
    return bins


def returns_at_vinf(
    model: RussellModel,
    body: str,
    vinf: float,
    *,
    max_tof_body_periods: float = 6.0,
    dtheta_deg: float = 0.5,
    max_revs_cap: int = 15,
) -> list[GenericReturn]:
    """Generic returns at an exact target ``|v_inf|`` (Russell §A.2).

    Generates the coarse grid, bins it into ``(n_revs, branch)`` sub-families,
    and for each consecutive pair bracketing the target ``vinf`` linearly
    interpolates a seed ``(psi, tof)``. The seed is then refined with a 1-D
    Kepler corrector that fixes ``|v_inf| = vinf`` and drives the two in-plane
    position residuals (re-encounter of the body) to zero. One refined
    :class:`GenericReturn` is returned per converged bracket.

    Parameters
    ----------
    model:
        The :class:`RussellModel`.
    body:
        Body key understood by :meth:`RussellModel.body_state`.
    vinf:
        Target hyperbolic-excess speed magnitude (canonical AU/TU).
    max_tof_body_periods, dtheta_deg, max_revs_cap:
        Passed through to :func:`generate_generic_returns`.

    Returns
    -------
    list[GenericReturn]
        Refined returns (``vinf`` exact by construction) for every bracket that
        converged to a sub-micro position miss; non-converged brackets skipped.
    """
    rs = generate_generic_returns(
        model,
        body,
        max_tof_body_periods=max_tof_body_periods,
        dtheta_deg=dtheta_deg,
        max_revs_cap=max_revs_cap,
    )
    bins = bin_sub_families(rs)

    p_tu = 2.0 * pi * model.sma_au(body) ** 1.5
    r1, v_b0 = model.body_state(body, 0.0)
    r1_n = float(np.linalg.norm(r1))

    e1 = v_b0 / np.linalg.norm(v_b0)
    rhat = r1 / r1_n
    e2 = rhat - np.dot(rhat, e1) * e1
    e2 = e2 / np.linalg.norm(e2)

    def residual(x: np.ndarray) -> np.ndarray:
        tof_canonical, psi = float(x[0]), float(x[1])
        vinf_vec = vinf * (np.cos(psi) * e1 + np.sin(psi) * e2)
        v0 = v_b0 + vinf_vec
        r_prop, _ = propagate(r1, v0, tof_canonical, model.mu_sun)
        theta_arr = 2.0 * pi * tof_canonical / p_tu
        r_body, _ = model.body_state(body, theta_arr)
        miss: np.ndarray = (r_prop - r_body)[:2]
        return miss

    refined: list[GenericReturn] = []
    seen: list[tuple[float, float]] = []  # (tof_bp, psi_deg) of accepted returns
    for (n_revs, _branch), lst in bins.items():
        # Bracket along the time-of-flight grid: within a single (n_revs, branch)
        # sub-family the v_inf-vs-ToF curve is non-monotonic, so consecutive
        # ToF-ordered grid points -- not the v_inf-sorted order used for the
        # public binning -- bracket EVERY crossing of the target |v_inf|.
        ordered = sorted(lst, key=lambda g: g.tof_body_periods)
        for g_lo, g_hi in pairwise(ordered):
            lo, hi = sorted((g_lo.vinf, g_hi.vinf))
            if not (lo <= vinf <= hi):
                continue
            span = g_hi.vinf - g_lo.vinf
            t = 0.0 if span == 0.0 else (vinf - g_lo.vinf) / span
            psi_seed_rad = np.radians(g_lo.psi_deg + t * (g_hi.psi_deg - g_lo.psi_deg))
            tof_bp_seed = g_lo.tof_body_periods + t * (
                g_hi.tof_body_periods - g_lo.tof_body_periods
            )

            sol = least_squares(
                residual,
                [tof_bp_seed * p_tu, psi_seed_rad],
                xtol=1e-12,
                ftol=1e-12,
            )
            if float(np.max(np.abs(sol.fun))) >= 1e-6:
                continue

            tof_solved, psi_solved = float(sol.x[0]), float(sol.x[1])
            vinf_vec = vinf * (np.cos(psi_solved) * e1 + np.sin(psi_solved) * e2)
            v0 = v_b0 + vinf_vec
            speed = float(np.linalg.norm(v0))
            denom = 2.0 / r1_n - speed * speed / model.mu_sun
            if abs(denom) < 1.0e-12:
                continue
            a_sma = 1.0 / denom
            if a_sma <= 0.0:
                continue
            a_au = a_sma**1.5

            tof_bp_solved = tof_solved / p_tu
            psi_deg_solved = degrees(psi_solved)
            # Adjacent ToF brackets can converge to the same crossing; dedupe.
            if any(
                abs(tof_bp_solved - tt) < 1e-4 and abs(psi_deg_solved - pp) < 1e-3
                for tt, pp in seen
            ):
                continue
            seen.append((tof_bp_solved, psi_deg_solved))

            refined.append(
                GenericReturn(
                    psi_deg=psi_deg_solved,
                    tof_body_periods=tof_bp_solved,
                    a_au=a_au,
                    n_revs=n_revs,
                    branch=_russell_branch(tof_solved, a_sma, model.mu_sun),
                    vinf=vinf,
                )
            )

    return refined
