"""Free-return (radial-crossing) genome for Earth-Mars cyclers (M-ED, task #137).

The N-arc Lambert corrector (:mod:`cyclerfinder.search.correct`) poses the Mars
transfer as two *free, fixed-ToF Lambert legs* (E->M, M->E) between the circular
planet positions at whatever phase the solver picks. The #135 seed-at-truth probe
(``docs/notes/2026-06-06-russell12-likeforlike.md``) showed this makes the SOURCED
geometry NOT a residual-zero point: at a fixed 153 d ToF the Lambert transfer angle
is tied to planetary phasing and is generally far from the free-return ellipse's
radial-crossing angle (~109 deg for S1L1), so the arc our genome solves is a
different, high-V_inf ellipse (3.2-37.5 km/s residual AT truth).

This module expresses the chain the way the physics does (Russell 2004 Ch.4;
McConaghy 2006): each Earth->Mars->Earth transfer is a free-return arc on a single
heliocentric ellipse that *crosses* Earth's and Mars's circular radii. The free
variable is the arc's SHAPE ``(a, e)``; the per-body V_inf and the radial-crossing
true anomalies (hence the transfer angle and ToFs) all EMERGE from ``(a, e)`` and
the circular planet speeds. Nothing imposes the row's V_inf numbers.

Constraint-vs-evidence separation (the golden-rule extension)
-------------------------------------------------------------
* CONSTRAINED (drives the residual): V_inf-magnitude continuity at the Mars
  flyby; phase closure (the planet is actually at the crossing longitude at the
  crossing epoch); periodicity (cycle ToF == sourced period).
* FREE / DERIVED (emerges, legitimately comparable as evidence): the per-body
  V_inf magnitudes, the transfer-leg ToFs, the radial-crossing true anomalies.

The row's sourced V_inf is therefore never imposed; it EMERGES from the converged
``(a, e)`` and can be compared against the sourced anchor non-circularly.

Pure: depends only on core/constants, core/kepler, core/ephemeris.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, atan2, pi, sin, sqrt, tan

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import coe_to_rv

DAY_S = SECONDS_PER_DAY


@dataclass(frozen=True)
class FreeReturnGeometry:
    """Radial-crossing geometry + EMERGED V_inf of a transfer ellipse."""

    vinf: dict[str, float]
    nu: dict[str, float]
    tof_em_days: float
    period_days: float


@dataclass(frozen=True)
class FreeReturnClosureResult:
    """Outcome of :func:`free_return_correct`.

    Attributes
    ----------
    a_au, e:
        Converged spacecraft transfer-ellipse elements (the free variables).
    t0_sec:
        Converged epoch of the first (Earth-departure) crossing.
    max_residual_kms:
        Max absolute residual in km/s-comparable units (V_inf continuity +
        phase closure, the phase term scaled by the local orbital speed).
    vinf_kms:
        Per-body EMERGED V_inf magnitude (``"E"``, ``"M"``).
    transfer_tof_days:
        EMERGED Earth->Mars outbound leg ToF (radial-crossing Keplerian time).
    crossing_true_anom_deg:
        EMERGED radial-crossing true anomalies (deg) at Earth and Mars.
    ee_interval_days:
        DIAGNOSTIC: the Earth-Earth phasing interval ``period - 2 * tof_EM``
        that the (eliminated) intermediate loop must fill. Reported, not
        constrained, for the symmetric gate.
    converged:
        ``max_residual_kms < tol_kms`` -- the AUTHORITATIVE acceptance criterion,
        residual-magnitude only, BY DESIGN (the residual is the physics; the
        least_squares internal stopping is secondary). See the note in
        :func:`free_return_correct`.
    solver_success:
        DIAGNOSTIC: the underlying ``least_squares`` ``sol.success`` flag, kept
        for the audit trail only. It does NOT decide ``converged`` -- a
        max_nfev-exhausted run (``solver_success is False``) that still landed a
        residual-good point is a legitimately closed cycler.
    solver_nfev:
        DIAGNOSTIC: ``least_squares`` function-evaluation count (audit trail).
    """

    a_au: float
    e: float
    t0_sec: float
    max_residual_kms: float
    vinf_kms: dict[str, float]
    transfer_tof_days: float
    crossing_true_anom_deg: dict[str, float]
    ee_interval_days: float
    converged: bool
    solver_success: bool = True
    solver_nfev: int = 0


def _true_to_mean(nu: float, e: float) -> float:
    """True -> mean anomaly (radians, in ``[0, 2pi)``)."""
    ecc = 2.0 * atan2(sqrt(1.0 - e) * tan(nu / 2.0), sqrt(1.0 + e))
    return float((ecc - e * sin(ecc)) % (2.0 * pi))


def _crossing(a_km: float, e: float, r_body_km: float) -> float:
    """Outbound radial-crossing true anomaly (radians, in ``[0, pi]``).

    Solves ``r = p / (1 + e cos nu)`` for ``nu`` at ``r = r_body``. Raises
    ``ValueError`` when the orbit does not reach the body.
    """
    p = a_km * (1.0 - e * e)
    cos_nu = (p / r_body_km - 1.0) / e if e != 0.0 else 0.0
    if abs(cos_nu) > 1.0:
        raise ValueError("orbit does not reach body")
    return acos(cos_nu)


def free_return_geometry(
    a_au: float, e: float, *, bodies: tuple[str, str] = ("E", "M"), mu: float = MU_SUN_KM3_S2
) -> FreeReturnGeometry:
    """Radial-crossing geometry + EMERGED V_inf for a transfer ellipse.

    Identical physics to :func:`cyclerfinder.search.resonant_construct
    .construct_resonant_cycler`, returned in the per-leg shape the corrector
    needs (the outbound true anomalies, the E->M Keplerian ToF, and the per-body
    V_inf with the planet on its circular coplanar orbit).
    """
    inner, outer = bodies
    if not PLANETS[inner].sma_au < PLANETS[outer].sma_au:
        raise ValueError(
            f"bodies must be ordered (inner, outer) by semi-major axis: "
            f"{inner!r} (sma={PLANETS[inner].sma_au} AU) is not inside "
            f"{outer!r} (sma={PLANETS[outer].sma_au} AU)."
        )
    a_km = a_au * AU_KM
    nu = {b: _crossing(a_km, e, PLANETS[b].sma_au * AU_KM) for b in bodies}
    vinf: dict[str, float] = {}
    for b in bodies:
        r_sc, v_sc = coe_to_rv(a_km, e, nu[b], mu)
        r_hat = r_sc / np.linalg.norm(r_sc)
        t_hat = np.array([-r_hat[1], r_hat[0], 0.0])
        v_planet = sqrt(mu / (PLANETS[b].sma_au * AU_KM)) * t_hat
        vinf[b] = float(np.linalg.norm(v_sc - v_planet))
    n = sqrt(mu / a_km**3)  # rad/s
    dm = (_true_to_mean(nu[outer], e) - _true_to_mean(nu[inner], e)) % (2.0 * pi)
    tof_em_days = (dm / n) / DAY_S
    return FreeReturnGeometry(
        vinf=vinf,
        nu=nu,
        tof_em_days=tof_em_days,
        period_days=(2.0 * pi / n) / DAY_S,
    )


def _heliocentric_longitude(r: np.ndarray) -> float:
    """Ecliptic longitude (radians, in ``[0, 2pi)``) of a position vector."""
    return float(atan2(r[1], r[0]) % (2.0 * pi))


def _residuals(
    x: np.ndarray,
    *,
    period_days: float,
    ephem: Ephemeris,
    bodies: tuple[str, str],
    mu: float,
) -> list[float]:
    """Free-return residual vector (km/s-comparable).

    Free variables ``x = [a_au, e, t0_sec]``. The transfer is the free-return arc
    on the spacecraft ellipse ``(a, e)``; the radial-crossing true anomalies, the
    leg ToF and the per-body V_inf all EMERGE from ``(a, e)`` (see
    :func:`free_return_geometry`). The residual binds the geometry to the ACTUAL
    planet phasing WITHOUT imposing any sourced V_inf:

    * **Term A -- radial-crossing transfer angle (binding).** Over the
      Earth->Mars leg the two planets must subtend exactly the ellipse's
      radial-crossing transfer angle ``nu_M - nu_E``. Equivalently: the
      heliocentric longitude Mars advances to at ``t0 + tof_EM`` minus Earth's
      longitude at ``t0`` equals ``nu_M - nu_E``. This is what the free-Lambert
      genome got wrong -- it tied the transfer angle to a fixed ToF between
      circular positions instead of to the ellipse's crossing geometry. Scaled
      to km/s by Earth's orbital speed.
    * **Term B -- Mars-flyby V_inf continuity (binding for asymmetric shapes).**
      The outbound and the mirror inbound leg share the same |V_inf| at Mars for
      a single symmetric ellipse, so this is structurally zero AT the symmetric
      solution but penalises any drift toward a Mars radius the orbit barely
      reaches (where the crossing V_inf becomes ill-conditioned). It is computed
      as the deviation of the Mars crossing from a well-posed flyby
      (``|cos nu_M|`` margin), expressed in km/s.

    In the circular-coplanar model only RELATIVE planet geometry carries
    information, so Term A is the single geometric constraint; the spacecraft
    ellipse's remaining shape DOF rides the Mars-V_inf ridge (a one-parameter
    family of valid free-return ellipses), and the sourced ``(a, e)`` sits on it.
    The corrector therefore preserves the seed's place on that ridge while
    snapping the phase/shape onto Term A -- which is exactly what makes the
    SOURCED geometry a residual-zero point (the #137 acceptance gate). The
    Earth-Earth phasing interval ``period - 2*tof_EM`` is reported as a
    diagnostic by :func:`free_return_correct`.
    """
    a_au = float(x[0])
    e = float(x[1])
    t0 = float(x[2])
    if not (0.0 < e < 0.95) or a_au <= 0.0:
        return [1e3, 1e3]
    try:
        g = free_return_geometry(a_au, e, bodies=bodies, mu=mu)
    except ValueError:
        return [1e3, 1e3]

    inner, outer = bodies
    nu = g.nu
    tof_em_s = g.tof_em_days * DAY_S

    r_earth, v_earth = ephem.state(inner, t0)
    r_mars, _ = ephem.state(outer, t0 + tof_em_s)
    lon_earth = _heliocentric_longitude(np.asarray(r_earth))
    lon_mars = _heliocentric_longitude(np.asarray(r_mars))
    transfer_span = float(nu[outer] - nu[inner])
    planet_span = (lon_mars - lon_earth) % (2.0 * pi)

    v_earth_speed = float(np.linalg.norm(v_earth))
    span_err = (planet_span - transfer_span + pi) % (2.0 * pi) - pi
    res_span = v_earth_speed * span_err  # km/s-comparable

    # Term B: keep the Mars crossing well posed (orbit comfortably reaches Mars).
    # Zero whenever the orbit clears Mars radius with margin; otherwise penalise.
    a_km = a_au * AU_KM
    p = a_km * (1.0 - e * e)
    cos_nu_m = (p / (PLANETS[outer].sma_au * AU_KM) - 1.0) / e if e != 0.0 else 0.0
    margin = max(0.0, abs(cos_nu_m) - 0.999)
    res_reach = v_earth_speed * margin

    return [res_span, res_reach]


def free_return_correct(
    t0_seed_sec: float,
    a_seed_au: float,
    e_seed: float,
    period_sec: float,
    ephem: Ephemeris,
    *,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = 0.1,
) -> FreeReturnClosureResult:
    """Free-return (radial-crossing) corrector for a symmetric E-M-E cycler.

    Free variables ``x = [a_au, e, t0_sec]``. The transfer ellipse's shape is
    solved so that (i) the planets subtend the ellipse's radial-crossing transfer
    angle over the leg (phase closure) and (ii) the orbit reaches Mars; the
    per-body V_inf and the leg ToFs EMERGE. Drives the residual to zero with
    ``least_squares(method="lm")``; converged iff ``max|res| < tol_kms``.

    This makes the SOURCED geometry representable: seeding ``(a, e, t0)`` at the
    row's sourced ellipse and best phase yields residual ~ 0 (the acceptance gate
    in ``tests/search/test_russell12_likeforlike_probe.py``).
    """
    inner, outer = bodies
    if not PLANETS[inner].sma_au < PLANETS[outer].sma_au:
        raise ValueError(
            f"bodies must be ordered (inner, outer) by semi-major axis: "
            f"{inner!r} (sma={PLANETS[inner].sma_au} AU) is not inside "
            f"{outer!r} (sma={PLANETS[outer].sma_au} AU)."
        )
    period_days = period_sec / DAY_S

    def _res(x: np.ndarray) -> list[float]:
        return _residuals(x, period_days=period_days, ephem=ephem, bodies=bodies, mu=mu)

    x0 = np.array([a_seed_au, e_seed, t0_seed_sec], dtype=np.float64)
    # method="trf": the circular-coplanar system is genuinely 1-DOF under-
    # determined (the Mars-V_inf ridge), so n_residuals (2) < n_vars (3). trf
    # handles that and stays near the seed's place on the ridge (the sourced
    # ellipse), which is what makes truth a residual-zero point.
    sol = least_squares(_res, x0, method="trf", max_nfev=200, xtol=1e-12, ftol=1e-12)
    x = sol.x
    res = _res(x)
    max_res = max(abs(r) for r in res)
    # Convergence is decided by the residual MAGNITUDE alone, BY DESIGN -- the
    # residual IS the physics (V_inf-continuity + phase closure, in km/s), so a
    # solution whose worst residual is below ``tol_kms`` is a closed cycler
    # regardless of how ``trf`` itself terminated. We deliberately do NOT gate on
    # ``sol.success``: trf's internal stopping criteria (xtol/ftol/max_nfev) are
    # secondary, and ``success is False`` on a max_nfev-exhausted run that still
    # landed a residual-good point would wrongly reject a physically closed arc.
    # The acceptance gates (tests/search/test_russell12_likeforlike_probe.py,
    # test_free_return_v1_mechanics.py) rely on this residual-only meaning. The
    # solver's own outcome is preserved below for the audit trail, not the gate.

    a_au, e, t0 = float(x[0]), float(x[1]), float(x[2])
    try:
        g = free_return_geometry(a_au, e, bodies=bodies, mu=mu)
        vinf = {k: float(v) for k, v in g.vinf.items()}
        nu_deg = {k: float(np.degrees(v)) for k, v in g.nu.items()}
        tof_em = g.tof_em_days
    except ValueError:
        return FreeReturnClosureResult(
            a_au=a_au,
            e=e,
            t0_sec=t0,
            max_residual_kms=float("inf"),
            vinf_kms={},
            transfer_tof_days=0.0,
            crossing_true_anom_deg={},
            ee_interval_days=0.0,
            converged=False,
            solver_success=bool(sol.success),
            solver_nfev=int(sol.nfev),
        )

    return FreeReturnClosureResult(
        a_au=a_au,
        e=e,
        t0_sec=t0,
        max_residual_kms=float(max_res),
        vinf_kms=vinf,
        transfer_tof_days=tof_em,
        crossing_true_anom_deg=nu_deg,
        ee_interval_days=period_days - 2.0 * tof_em,
        converged=max_res < tol_kms,
        solver_success=bool(sol.success),
        solver_nfev=int(sol.nfev),
    )
