"""Sun-Mars weak-stability-boundary (WSB) / ballistic-capture machinery (#681).

The interplanetary analog of the cislunar WSB substrate (`core/wsb.py`,
`genome/bct_transfer.py`, `search/cislunar_bct_search.py`, task #378), built for
the Sun-Mars system to reproduce Topputo & Belbruno, "Earth-Mars transfers with
ballistic capture," Celestial Mechanics and Dynamical Astronomy 121:329-346
(2015), doi:10.1007/s10569-015-9605-8, and to search for a *repeating*
capture<->escape chain (a `quasi_cycler`-class object whose return leg
re-acquires the Mars WSB set each cycle).

Model
-----
Planar restricted three-body problem, Sun + Mars + massless spacecraft, in a
**heliocentric inertial** frame (Sun fixed at the origin). Mars moves on a
Keplerian ellipse about the Sun (semi-major axis ``a_M``, eccentricity ``e_M``,
argument of periapsis 0 -- perihelion on +x), which reproduces the paper's
*elliptic* restricted problem to O(mu) = O(3e-7) (the barycentre/Sun-centre
offset is negligible against the paper's own ~2 km/s cost figures). Working
units are physical: km, s, km/s -- so Mars-relative Kepler energies and speeds
are directly interpretable and directly comparable to the paper's tables.

Sourcing / honesty
------------------
* Physical constants are the paper's own Table 4 values (module constants below),
  used verbatim so the reproduction is checkable against the paper's numbers
  rather than against a slightly different in-repo constant set.
* The Mars-relative Kepler energy and the periapsis predicate are closed forms
  (Belbruno 2004 Def 3.10 / eq 3.9, transcribed to the Sun-Mars system).
* This module asserts NO catalogue object. It makes the repeating-capture
  hypothesis *searchable*; the honest prior (task #378 cislunar clean negative;
  Belbruno 2004 Thm 3.58 chaotic capture; the paper's own Sect. 4.2 50-revolution
  non-recurrence) is that no repeating object exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

# --- Physical constants: Topputo & Belbruno 2015, Table 4 (verbatim) ---------
MU_SUN_KM3_S2: float = 1.32712e11
"""Solar gravitational parameter, km^3/s^2 (Table 4)."""

MU_MARS_KM3_S2: float = 4.28280e4
"""Mars gravitational parameter, km^3/s^2 (Table 4)."""

AU_KM: float = 149_597_870.66
"""Astronomical unit, km (Table 4)."""

MARS_A_KM: float = 1.523688399 * AU_KM
"""Mars orbital semi-major axis, km (Table 4: a = 1.523688399 AU)."""

MARS_E: float = 0.093418671
"""Mars orbital eccentricity (Table 4)."""

MARS_RADIUS_KM: float = 3389.5
"""Mars mean radius, km (IAU) -- collision floor for the integrator."""

EARTH_A_KM: float = 1.000000230 * AU_KM
"""Earth orbital semi-major axis, km (Table 4)."""

EARTH_E: float = 0.016751040
"""Earth orbital eccentricity (Table 4)."""

# Derived Mars mean motion (Sun-centred two-body, mu = MU_SUN; the O(mu_M)
# correction to Mars's own orbit is ~3e-7 and irrelevant here).
MARS_MEAN_MOTION: float = math.sqrt(MU_SUN_KM3_S2 / MARS_A_KM**3)
"""Mars mean motion, rad/s."""

MARS_PERIOD_S: float = 2.0 * math.pi / MARS_MEAN_MOTION
"""Mars sidereal period, s (~1.881 yr)."""

Branch = Literal["prograde", "retrograde"]


def _solve_kepler_e(mean_anom: float, ecc: float, *, tol: float = 1e-14) -> float:
    """Solve Kepler's equation ``M = E - e sin E`` for the eccentric anomaly."""
    m = math.remainder(mean_anom, 2.0 * math.pi)
    e_anom = m if ecc < 0.8 else math.pi
    for _ in range(80):
        f = e_anom - ecc * math.sin(e_anom) - m
        fp = 1.0 - ecc * math.cos(e_anom)
        step = f / fp
        e_anom -= step
        if abs(step) < tol:
            break
    return e_anom


def _true_anom_to_mean(f: float, ecc: float) -> float:
    """True anomaly -> mean anomaly for eccentricity ``ecc``."""
    e_anom = 2.0 * math.atan2(
        math.sqrt(1.0 - ecc) * math.sin(f / 2.0),
        math.sqrt(1.0 + ecc) * math.cos(f / 2.0),
    )
    return e_anom - ecc * math.sin(e_anom)


def body_state(
    t: float, *, a_km: float, ecc: float, f0: float, mu: float = MU_SUN_KM3_S2
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Heliocentric position/velocity (km, km/s) of a body on a Kepler ellipse.

    The body is at true anomaly ``f0`` at ``t = 0`` (perihelion on +x, prograde
    about +z) and advances Keplerian in time ``t`` (s).
    """
    n = math.sqrt(mu / a_km**3)
    m0 = _true_anom_to_mean(f0, ecc)
    m = m0 + n * t
    e_anom = _solve_kepler_e(m, ecc)
    cos_e = math.cos(e_anom)
    sin_e = math.sin(e_anom)
    r = a_km * (1.0 - ecc * cos_e)
    x = a_km * (cos_e - ecc)
    y = a_km * math.sqrt(1.0 - ecc**2) * sin_e
    vfac = math.sqrt(mu * a_km) / r
    vx = -vfac * sin_e
    vy = vfac * math.sqrt(1.0 - ecc**2) * cos_e
    return (
        np.array([x, y], dtype=np.float64),
        np.array([vx, vy], dtype=np.float64),
    )


def mars_state(t: float, f0: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Heliocentric Mars position/velocity (km, km/s) at time ``t`` (s)."""
    return body_state(t, a_km=MARS_A_KM, ecc=MARS_E, f0=f0)


def sunmars_eom(t: float, state4: NDArray[np.float64], f0: float) -> NDArray[np.float64]:
    """Planar restricted 3-body EOM (Sun + Mars), heliocentric inertial frame.

    ``state4 = [x, y, vx, vy]`` in km, km/s. ``f0`` is Mars's true anomaly at
    ``t = 0``.
    """
    x, y, vx, vy = (float(v) for v in state4)
    r_s = np.array([x, y], dtype=np.float64)
    r_m, _ = mars_state(t, f0)
    d_s = float(np.hypot(x, y))
    rel = r_s - r_m
    d_m = float(np.hypot(rel[0], rel[1]))
    acc = -MU_SUN_KM3_S2 * r_s / d_s**3 - MU_MARS_KM3_S2 * rel / d_m**3
    return np.array([vx, vy, acc[0], acc[1]], dtype=np.float64)


def mars_kepler_energy(state4: NDArray[np.float64], t: float, f0: float) -> float:
    """Two-body Kepler energy E_2 of the spacecraft w.r.t. Mars (Belbruno Def 3.10).

    ``E_2 = 1/2 |v - v_M|^2 - mu_M / |r - r_M|``. ``E_2 <= 0`` is ballistic
    capture (Def 3.11).
    """
    r_m, v_m = mars_state(t, f0)
    rel = np.asarray(state4[:2], dtype=np.float64) - r_m
    vrel = np.asarray(state4[2:], dtype=np.float64) - v_m
    d = float(np.hypot(rel[0], rel[1]))
    if d == 0.0:
        return math.inf
    return 0.5 * float(vrel @ vrel) - MU_MARS_KM3_S2 / d


def mars_distance(state4: NDArray[np.float64], t: float, f0: float) -> float:
    """Spacecraft-Mars distance (km)."""
    r_m, _ = mars_state(t, f0)
    rel = np.asarray(state4[:2], dtype=np.float64) - r_m
    return float(np.hypot(rel[0], rel[1]))


def mars_radial_rate(state4: NDArray[np.float64], t: float, f0: float) -> float:
    """Mars-relative radial rate d/dt|r - r_M| (km/s). Zero at peri/apoapsis."""
    r_m, v_m = mars_state(t, f0)
    rel = np.asarray(state4[:2], dtype=np.float64) - r_m
    vrel = np.asarray(state4[2:], dtype=np.float64) - v_m
    d = float(np.hypot(rel[0], rel[1]))
    if d == 0.0:
        return 0.0
    return float(rel @ vrel) / d


def capture_periapsis_state(
    *, r_p_km: float, ecc: float, theta: float, f0: float, branch: Branch = "prograde"
) -> NDArray[np.float64]:
    """Build a spacecraft state at a Mars periapsis on the osculating ellipse.

    Position: at angle ``theta`` (rad) around Mars at radius ``r_p_km``.
    Velocity: Mars-relative periapsis speed ``sqrt(mu_M (1+e)/r_p)`` tangential
    to the radial direction, plus Mars's heliocentric velocity. By construction
    the Mars-relative radial rate is 0 and ``E_2 = mu_M (e-1) / (2 r_p) < 0``
    (a bound periapsis) -- the paper's grid seed (Sect. 4).
    """
    r_m, v_m = mars_state(0.0, f0)
    radial = np.array([math.cos(theta), math.sin(theta)], dtype=np.float64)
    r_s = r_m + r_p_km * radial
    v_peri = math.sqrt(MU_MARS_KM3_S2 * (1.0 + ecc) / r_p_km)
    sign = 1.0 if branch == "prograde" else -1.0
    tangent = np.array([-math.sin(theta), math.cos(theta)], dtype=np.float64) * sign
    v_s = v_m + v_peri * tangent
    return np.array([r_s[0], r_s[1], v_s[0], v_s[1]], dtype=np.float64)


@dataclass(frozen=True)
class PeriapsisEvent:
    """A Mars-relative periapsis crossing on an integrated arc."""

    t_s: float
    dist_km: float
    e2: float


@dataclass(frozen=True)
class StabilityResult:
    """Outcome of a forward/backward n-stability integration (Belbruno Sect. 4).

    Attributes
    ----------
    n_captured_revs :
        Number of Mars periapsis passages with ``E_2 <= 0`` reached before the
        first escape (the paper's forward n-stability count).
    escaped :
        True iff the spacecraft left the Mars vicinity unbound (a periapsis with
        ``E_2 > 0``, i.e. unbound even at its closest Mars approach).
    recaptured_after_escape :
        True iff, AFTER an escape, ANY periapsis is bound (``E_2 <= 0``) within
        ``hill_km`` -- the LOOSE re-acquisition flag (kept as a diagnostic; a
        single such dip during a co-orbital conjunction is NOT a genuine
        recapture -- see ``n_recapture_episodes``).
    n_recapture_episodes :
        Number of distinct post-escape *episodes* of SUSTAINED capture -- maximal
        runs of >= 2 consecutive bound periapses (``E_2 <= 0``) within
        ``hill_km`` (i.e. at least one genuine bound Mars revolution). A
        *repeating* capture chain (the quasi-cycler object class) requires this
        to be >= 2; a single episode is an isolated temporary capture, not a
        cycler.
    max_sustained_bound_revs :
        Longest run of consecutive bound periapses within ``hill_km`` after an
        escape (the depth of the deepest temporary recapture).
    min_recapture_dist_km :
        Closest Mars distance of any post-escape bound periapsis within
        ``hill_km`` (``None`` if there is none).
    collided :
        True iff the spacecraft fell below the Mars radius.
    periapses :
        All Mars-relative periapsis events on the arc, in integration order.
    t_final_s :
        Integration end time (s, signed by direction).
    """

    n_captured_revs: int
    escaped: bool
    recaptured_after_escape: bool
    n_recapture_episodes: int
    max_sustained_bound_revs: int
    min_recapture_dist_km: float | None
    collided: bool
    periapses: tuple[PeriapsisEvent, ...]
    t_final_s: float


def integrate_stability(
    state0: NDArray[np.float64],
    f0: float,
    *,
    direction: Literal["forward", "backward"] = "forward",
    horizon_revs: float = 6.0,
    hill_km: float = 1.0e6,
    rtol: float = 1e-11,
    atol: float = 1e-9,
    max_step_frac: float = 0.02,
) -> StabilityResult:
    """Integrate a capture state and classify its Mars WSB stability.

    ``horizon_revs`` is measured in *Mars* sidereal periods. Mars-relative
    periapses are found EXACTLY via a radial-rate zero-crossing event (robust to
    step size, unlike sampling). ``hill_km`` (~Mars Hill radius) bounds the
    "captured near Mars" region. Pre-escape consecutive bound periapses are the
    forward n-stability count; the first unbound periapsis (``E_2 > 0``) marks
    escape; post-escape, maximal runs of >= 2 consecutive bound periapses within
    ``hill_km`` are counted as sustained-recapture *episodes* (a repeating chain
    needs >= 2 episodes).
    """
    sign = 1.0 if direction == "forward" else -1.0
    horizon_tau = horizon_revs * MARS_PERIOD_S

    # Integrate in a monotone-increasing pseudo-time tau = |t|; map back via sign.
    def rhs_tau(tau: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return sign * sunmars_eom(sign * tau, y, f0)

    def collision_event(tau: float, y: NDArray[np.float64]) -> float:
        return mars_distance(y, sign * tau, f0) - MARS_RADIUS_KM

    collision_event.terminal = True  # type: ignore[attr-defined]
    collision_event.direction = -1.0  # type: ignore[attr-defined]

    def periapsis_event(tau: float, y: NDArray[np.float64]) -> float:
        return mars_radial_rate(y, sign * tau, f0)

    periapsis_event.direction = 1.0  # type: ignore[attr-defined]  # - -> + is periapsis

    max_step = max_step_frac * MARS_PERIOD_S
    sol = solve_ivp(
        rhs_tau,
        (0.0, horizon_tau),
        state0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=False,
        max_step=max_step,
        events=(collision_event, periapsis_event),
    )
    collided = bool(sol.t_events is not None and len(sol.t_events[0]) > 0)

    periapses: list[PeriapsisEvent] = []
    if sol.t_events is not None and sol.y_events is not None:
        for tau_ev, y_ev in zip(sol.t_events[1], sol.y_events[1], strict=True):
            t = sign * float(tau_ev)
            periapses.append(
                PeriapsisEvent(
                    t_s=t,
                    dist_km=mars_distance(y_ev, t, f0),
                    e2=mars_kepler_energy(y_ev, t, f0),
                )
            )

    n_captured = 0
    escaped = False
    recaptured = False
    n_episodes = 0
    max_run = 0
    run = 0
    min_recap: float | None = None
    for p in periapses:
        if not escaped:
            if p.e2 > 0.0:
                escaped = True
            elif p.dist_km < hill_km:
                n_captured += 1
            continue
        # Post-escape: track sustained-recapture episodes within the Hill sphere.
        if p.e2 <= 0.0 and p.dist_km < hill_km:
            recaptured = True
            run += 1
            if run == 2:
                # A run reaching length 2 = >= 1 genuine bound revolution = a
                # SUSTAINED recapture episode (counted once per maximal run).
                n_episodes += 1
            max_run = max(max_run, run)
            min_recap = p.dist_km if min_recap is None else min(min_recap, p.dist_km)
        else:
            run = 0

    return StabilityResult(
        n_captured_revs=n_captured,
        escaped=escaped,
        recaptured_after_escape=recaptured,
        n_recapture_episodes=n_episodes,
        max_sustained_bound_revs=max_run,
        min_recapture_dist_km=min_recap,
        collided=collided,
        periapses=tuple(periapses),
        t_final_s=sign * float(sol.t[-1]),
    )


@dataclass(frozen=True)
class HohmannBaseline:
    """A bitangential Hohmann Earth->Mars reference (paper Table 5)."""

    dv1_kms: float
    dv2_inf_kms: float
    dv_total_kms: float
    tof_days: float


def hohmann_baseline(
    *, earth_apsis: Literal["peri", "apo"], mars_apsis: Literal["peri", "apo"]
) -> HohmannBaseline:
    """Bitangential Hohmann transfer between apsidal points (paper Table 5).

    Departs an Earth apsidal radius, arrives at a Mars apsidal radius on an
    apse-aligned transfer ellipse. ``dv1`` is the heliocentric departure burn;
    ``dv2_inf`` the arrival relative (v-infinity) speed at Mars.
    """
    r1 = EARTH_A_KM * (1.0 - EARTH_E if earth_apsis == "peri" else 1.0 + EARTH_E)
    r2 = MARS_A_KM * (1.0 - MARS_E if mars_apsis == "peri" else 1.0 + MARS_E)
    a_t = 0.5 * (r1 + r2)

    def vis(r: float, a: float) -> float:
        return math.sqrt(MU_SUN_KM3_S2 * (2.0 / r - 1.0 / a))

    dv1 = abs(vis(r1, a_t) - vis(r1, EARTH_A_KM))
    dv2 = abs(vis(r2, MARS_A_KM) - vis(r2, a_t))
    tof = math.pi * math.sqrt(a_t**3 / MU_SUN_KM3_S2) / 86400.0
    return HohmannBaseline(dv1_kms=dv1, dv2_inf_kms=dv2, dv_total_kms=dv1 + dv2, tof_days=tof)
